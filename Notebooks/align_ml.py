"""Alignment-aware feature extraction + classical ML for per-position
error detection. This is the part of the original notebook (cell 19)
that was actually well-designed but never got wired up correctly.
"""
import math, random
from collections import Counter
import numpy as np

BASES = list('ACGT')
LABEL_NAMES = {0: 'clean', 1: 'substitution', 2: 'insertion', 3: 'deletion'}


def needleman_wunsch(ref, query, match=2, mismatch=-1, gap=-2):
    n, m = len(ref), len(query)
    dp = [[0]*(m+1) for _ in range(n+1)]
    for i in range(1, n+1): dp[i][0] = i*gap
    for j in range(1, m+1): dp[0][j] = j*gap
    for i in range(1, n+1):
        rb = ref[i-1]
        for j in range(1, m+1):
            s = match if rb == query[j-1] else mismatch
            dp[i][j] = max(dp[i-1][j-1]+s, dp[i-1][j]+gap, dp[i][j-1]+gap)
    ra, qa = [], []
    i, j = n, m
    while i > 0 or j > 0:
        if i > 0 and j > 0:
            s = match if ref[i-1] == query[j-1] else mismatch
            if dp[i][j] == dp[i-1][j-1] + s:
                ra.append(ref[i-1]); qa.append(query[j-1]); i -= 1; j -= 1; continue
        if i > 0 and dp[i][j] == dp[i-1][j] + gap:
            ra.append(ref[i-1]); qa.append('-'); i -= 1; continue
        if j > 0 and dp[i][j] == dp[i][j-1] + gap:
            ra.append('-'); qa.append(query[j-1]); j -= 1; continue
        break
    ra.reverse(); qa.reverse()
    return ''.join(ra), ''.join(qa)


def classify_alignment(ref_a, q_a):
    labels = []
    for r, q in zip(ref_a, q_a):
        if r == '-' and q != '-': labels.append(2)      # insertion
        elif r != '-' and q == '-': labels.append(3)    # deletion
        elif r != q: labels.append(1)                    # substitution
        else: labels.append(0)                            # clean
    return labels


def one_hot(b):
    return [1 if b == x else 0 for x in BASES]


def gc_win(seq, pos, w=7):
    if not seq: return 0.0
    pos = max(0, min(pos, len(seq)-1))
    s = seq[max(0, pos-w//2):min(len(seq), pos+w//2+1)]
    return (s.count('G')+s.count('C'))/len(s) if s else 0.0


def hp_run(seq, pos):
    if not seq or pos < 0 or pos >= len(seq): return 0
    b = seq[pos]
    run = 1; i = pos-1
    while i >= 0 and seq[i] == b: run += 1; i -= 1
    i = pos+1
    while i < len(seq) and seq[i] == b: run += 1; i += 1
    return run


def extract_features(aligned_ref, aligned_noisy, pos):
    """One feature vector per alignment column. Uses only information
    that is genuinely available in the noisy read (local base context,
    homopolymer run, GC window) plus the alignment gap structure -
    this is what makes it learnable at all for substitutions/insertions.
    Deletion *detection* is learnable from neighbouring noisy context;
    deletion *base recovery* fundamentally cannot come from the noisy
    read alone (the base was never observed) - see README note in the
    notebook about oracle vs blind correction.
    """
    ref_b = aligned_ref[pos]
    noisy_b = aligned_noisy[pos]
    ungapped = aligned_noisy.replace('-', '')
    q_before = aligned_noisy[:pos].replace('-', '')
    qpos = min(len(q_before), max(len(ungapped)-1, 0))

    def at(seq, i, default='A'):
        return seq[i] if 0 <= i < len(seq) else default

    b = noisy_b if noisy_b != '-' else '-'
    prev, nxt = at(ungapped, qpos-1), at(ungapped, qpos)
    p2, n2 = at(ungapped, qpos-2), at(ungapped, qpos+1)

    hp = hp_run(ungapped, min(qpos, max(len(ungapped)-1, 0))) if ungapped else 0
    gc7 = gc_win(ungapped, qpos, 7)
    gc15 = gc_win(ungapped, qpos, 15)

    feats = (
        one_hot(b if b in BASES else 'A') +
        one_hot(prev) + one_hot(nxt) + one_hot(p2) + one_hot(n2) +
        [gc7, gc15, hp/10.0, 1.0 if hp >= 3 else 0.0,
         pos/max(len(aligned_ref), 1),
         1.0 if b == prev else 0.0, 1.0 if b == nxt else 0.0,
         1.0 if ref_b == '-' else 0.0, 1.0 if noisy_b == '-' else 0.0]
    )
    return feats


def build_dataset(strand_read_pairs, samples_per_class=4000, seed=0):
    """strand_read_pairs: list of (clean_strand, noisy_read_seq)."""
    rng = random.Random(seed)
    by_class = {0: [], 1: [], 2: [], 3: []}
    for clean, noisy in strand_read_pairs:
        ra, qa = needleman_wunsch(clean, noisy)
        labels = classify_alignment(ra, qa)
        for pos, lab in enumerate(labels):
            by_class[lab].append((ra, qa, pos))
    X, y = [], []
    for lab, items in by_class.items():
        take = items if len(items) <= samples_per_class else rng.sample(items, samples_per_class)
        for ra, qa, pos in take:
            X.append(extract_features(ra, qa, pos))
            y.append(lab)
    return np.array(X, dtype=np.float32), np.array(y, dtype=np.int64), {k: len(v) for k, v in by_class.items()}


def correct_with_labels(aligned_ref, aligned_noisy, labels, oracle=False):
    """Rebuild a corrected sequence from predicted per-column labels.
    oracle=True fills deletions with the TRUE reference base (upper
    bound on what detection-only correction can achieve - the base
    itself was never observed in the noisy read, so this is not a
    blind result, just a diagnostic ceiling).
    oracle=False drops deletion columns instead of guessing (honest
    blind behaviour without a base-calling model)."""
    out = []
    for i, lab in enumerate(labels):
        r, q = aligned_ref[i], aligned_noisy[i]
        if lab == 0 or lab == 1:
            if q != '-':
                out.append(q)
        elif lab == 2:
            continue  # insertion: drop the extra base
        elif lab == 3:
            if oracle and r != '-':
                out.append(r)
            # blind: nothing we can honestly add
    return ''.join(out)
