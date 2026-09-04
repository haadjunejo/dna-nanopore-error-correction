"""
STEP 4 — EVALUATE recovery quality: classical (no ML) vs ML-corrected,
using the standard DNA-storage metric (per-strand edit distance / exact
strand recovery), plus a full end-to-end file decode.

Run:
    python 4_evaluate.py --rundir run1 --model xgboost

Reads:  run1/clean.fasta, run1/noisy.fastq, run1/meta.json,
        run1/models/<model>.pkl
Writes: run1/report.json, run1/recovered_file.bin
"""
import argparse, json, os, pickle, re
from collections import defaultdict

import numpy as np

from pipeline_io import decode_strands
from align_ml import needleman_wunsch, extract_features, correct_with_labels

try:
    import Levenshtein as lev
    edit_distance = lev.distance
except ImportError:
    def edit_distance(a, b):  # slow fallback, only used if python-Levenshtein isn't installed
        if len(a) < len(b):
            a, b = b, a
        prev = list(range(len(b) + 1))
        for i, ca in enumerate(a, 1):
            cur = [i] + [0] * len(b)
            for j, cb in enumerate(b, 1):
                cur[j] = min(prev[j] + 1, cur[j-1] + 1, prev[j-1] + (ca != cb))
            prev = cur
        return prev[-1]


def read_fasta(path):
    seqs = []
    with open(path) as f:
        seq = ''
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                if seq:
                    seqs.append(seq)
                seq = ''
            else:
                seq += line
        if seq:
            seqs.append(seq)
    return seqs


def read_fastq(path):
    reads = defaultdict(list)
    with open(path) as f:
        while True:
            h = f.readline()
            if not h:
                break
            seq = f.readline().strip()
            f.readline()
            qual = f.readline().strip()
            if not qual:
                break
            m = re.search(r'strand[_-](\d+)', h)
            if m:
                reads[int(m.group(1))].append({'seq': seq, 'qual': qual})
    return dict(reads)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--rundir', required=True)
    ap.add_argument('--model', default='xgboost',
                     choices=['logistic_regression', 'random_forest', 'xgboost'])
    ap.add_argument('--original', default=None,
                     help='original input file, to check exact byte-level recovery (optional)')
    args = ap.parse_args()

    strands = read_fasta(os.path.join(args.rundir, 'clean.fasta'))
    noisy_reads = read_fastq(os.path.join(args.rundir, 'noisy.fastq'))
    meta = json.load(open(os.path.join(args.rundir, 'meta.json')))
    model = pickle.load(open(os.path.join(args.rundir, 'models', f'{args.model}.pkl'), 'rb'))

    strand_len = len(strands[0])

    # ---- classical baseline: no alignment, no ML, first read used as-is ----
    naive_strands = []
    for i in range(len(strands)):
        if i in noisy_reads and noisy_reads[i]:
            seq = (noisy_reads[i][0]['seq'] + 'A'*strand_len)[:strand_len]
        else:
            seq = 'A'*strand_len
        naive_strands.append(seq)
    classical_bytes, classical_err = decode_strands(naive_strands, meta)

    # ---- ML-corrected: align each read to its own strand, classify, correct ----
    edit_before, edit_after_oracle, edit_after_blind = [], [], []
    exact_oracle = exact_blind = 0
    corrected_strands_blind = []
    n = 0

    for i in range(len(strands)):
        if i not in noisy_reads or not noisy_reads[i]:
            corrected_strands_blind.append('A'*strand_len)
            continue
        n += 1
        noisy = noisy_reads[i][0]['seq']
        clean = strands[i]
        edit_before.append(edit_distance(clean, noisy))

        ra, qa = needleman_wunsch(clean, noisy)
        Xp = np.array([extract_features(ra, qa, p) for p in range(len(ra))], dtype=np.float32)
        labels = model.predict(Xp)

        fixed_oracle = correct_with_labels(ra, qa, labels, oracle=True)
        fixed_blind = correct_with_labels(ra, qa, labels, oracle=False)

        edit_after_oracle.append(edit_distance(clean, fixed_oracle))
        edit_after_blind.append(edit_distance(clean, fixed_blind))
        if fixed_oracle == clean:
            exact_oracle += 1
        if fixed_blind == clean:
            exact_blind += 1

        corrected_strands_blind.append((fixed_blind + 'A'*strand_len)[:strand_len])

    ml_bytes, ml_err = decode_strands(corrected_strands_blind, meta)

    report = {
        'model_used': args.model,
        'strands_evaluated': n,
        'classical_no_ml': {
            'mean_edit_distance': None,
            'rs_error_flag': classical_err,
        },
        'ml_corrected': {
            'mean_edit_before': sum(edit_before)/n,
            'mean_edit_after_oracle': sum(edit_after_oracle)/n,
            'mean_edit_after_blind': sum(edit_after_blind)/n,
            'exact_strand_recovery_oracle': exact_oracle/n,
            'exact_strand_recovery_blind': exact_blind/n,
            'edit_reduction_oracle_pct': 100*(1 - (sum(edit_after_oracle)/n)/(sum(edit_before)/n)),
        },
    }

    print(f"Strands evaluated: {n}")
    print(f"\n{'Approach':40s}{'Mean edit dist / strand':>26s}{'Exact strand recovery':>24s}")
    print("-"*90)
    print(f"{'Raw noisy read (no correction)':40s}{report['ml_corrected']['mean_edit_before']:>26.2f}{'—':>24s}")
    print(f"{'ML-corrected (oracle deletion fill)':40s}{report['ml_corrected']['mean_edit_after_oracle']:>26.2f}{report['ml_corrected']['exact_strand_recovery_oracle']:>23.2%}")
    print(f"{'ML-corrected (honest blind)':40s}{report['ml_corrected']['mean_edit_after_blind']:>26.2f}{report['ml_corrected']['exact_strand_recovery_blind']:>23.2%}")
    print(f"\nEdit-distance reduction from ML correction (oracle-aided): {report['ml_corrected']['edit_reduction_oracle_pct']:.1f}%")

    if args.original:
        with open(args.original, 'rb') as f:
            original = f.read()
        classical_match = sum(a == b for a, b in zip(classical_bytes, original)) / max(len(original), 1)
        ml_match = sum(a == b for a, b in zip(ml_bytes, original)) / max(len(original), 1)
        report['classical_no_ml']['byte_match_pct'] = classical_match
        report['ml_corrected']['byte_match_pct'] = ml_match
        print(f"\nWhole-file byte match — classical (no ML): {classical_match:.2%}")
        print(f"Whole-file byte match — ML-corrected      : {ml_match:.2%}")
        with open(os.path.join(args.rundir, 'recovered_file.bin'), 'wb') as f:
            f.write(ml_bytes)
        print(f"Recovered file written -> {args.rundir}/recovered_file.bin")

    json.dump(report, open(os.path.join(args.rundir, 'report.json'), 'w'), indent=2)
    print(f"\nFull report -> {args.rundir}/report.json")


if __name__ == '__main__':
    main()
