"""
Core encode/decode logic for the DNA storage pipeline.
Fixes vs. the original notebook:
  - strand index width is derived from the actual number of strands
    (was hard-coded to 4 nt = 8 bits = 256 strands max)
  - index block carries its own parity nt so a corrupted index is
    DETECTED instead of silently misordering strands
  - marker stripping is search-based (tolerant of small shifts) instead
    of assuming markers sit at a fixed offset
"""
import struct, heapq, random, math
from collections import Counter

B2D = {'00':'A','01':'C','10':'G','11':'T'}
D2B = {v:k for k,v in B2D.items()}
BASES = list('ACGT')
MARKER = "ACGTACGT"

GC_MIN, GC_MAX = 0.40, 0.60
MAX_HP = 3

# ---------------- Huffman ----------------
class _HNode:
    def __init__(self, c, f): self.c=c; self.f=f; self.l=self.r=None
    def __lt__(self, o): return self.f < o.f

def _huff_tree(data):
    heap = [_HNode(c,f) for c,f in Counter(data).items()]
    heapq.heapify(heap)
    if len(heap) == 1:
        # single-symbol edge case
        only = heap[0]
        return _HNode(None, only.f, ) if False else only
    while len(heap) > 1:
        l=heapq.heappop(heap); r=heapq.heappop(heap)
        m=_HNode(None,l.f+r.f); m.l=l; m.r=r; heapq.heappush(heap,m)
    return heap[0]

def _huff_codes(node, prefix=''):
    if node is None: return {}
    if node.c is not None: return {node.c: prefix or '0'}
    return {**_huff_codes(node.l, prefix+'0'), **_huff_codes(node.r, prefix+'1')}

def huffman_compress(data: bytes):
    if not data: return '', {}
    codes = _huff_codes(_huff_tree(data))
    return ''.join(codes[b] for b in data), codes

def huffman_decompress(bits: str, codes: dict) -> bytes:
    rev = {v:k for k,v in codes.items()}
    out, cur = [], ''
    for b in bits:
        cur += b
        if cur in rev:
            out.append(rev[cur]); cur=''
    return bytes(out)

# ---------------- LZ77 ----------------
def lz77_compress(data: bytes, win=128, look=12):
    tokens, i = [], 0
    n = len(data)
    while i < n:
        best_off, best_len = 0, 0
        for j in range(max(0,i-win), i):
            l = 0
            d = i - j
            while l < look and i+l < n and data[j+(l % d)] == data[i+l]:
                l += 1
            if l > best_len: best_len=l; best_off=d
        if best_len > 2:
            nc = data[i+best_len] if i+best_len < n else 0
            tokens.append((best_off, best_len, nc)); i += best_len+1
        else:
            tokens.append((0, 0, data[i])); i += 1
    return tokens

def lz77_decompress(tokens):
    """Defensive: if upstream corruption produced an invalid (off, len)
    pair, stop cleanly instead of raising, so one bad strand doesn't
    crash the whole file recovery. Callers should check whether the
    output length looks right and treat short output as partial loss."""
    out = bytearray()
    for off, ln, ch in tokens:
        try:
            if off == 0 and ln == 0:
                out.append(ch)
            else:
                s = len(out) - off
                if s < 0:
                    break
                for k in range(ln):
                    out.append(out[s+k])
                out.append(ch)
        except IndexError:
            break
    return bytes(out)

def compress_file_bytes(raw: bytes):
    tok = lz77_compress(raw)
    tb  = b''.join(struct.pack('BBB',o,l,c) for o,l,c in tok)
    bits, codes = huffman_compress(tb)
    if len(bits) % 2:
        bits += '0'
    meta = {'orig_size': len(raw), 'bits': len(bits),
            'codes': {str(k): v for k, v in codes.items()}, 'n_tokens': len(tok)}
    return bits, meta

def decompress_bits(bits, meta):
    codes = {int(k): v for k, v in meta['codes'].items()}
    if not codes:
        return b''
    tb = huffman_decompress(bits[:meta['bits']], codes)
    tok = [struct.unpack('BBB', tb[i:i+3]) for i in range(0, len(tb)-2, 3)]
    out = lz77_decompress(tok)
    # lz77_compress appends a placeholder literal byte even when a match
    # runs exactly to the end of the input (there's no real "next byte"
    # to encode); that adds one spurious trailing byte. orig_size is
    # known, so trim to it rather than trusting the token stream's length.
    return out[:meta['orig_size']]

# ---------------- Binary <-> DNA ----------------
def bits_to_dna(bits):
    if len(bits) % 2:
        bits += '0'
    return ''.join(B2D[bits[i:i+2]] for i in range(0, len(bits), 2))

def dna_to_bits(dna):
    return ''.join(D2B.get(b, '00') for b in dna)

# ---------------- Biological constraints ----------------
def gc(seq):
    return (seq.count('G') + seq.count('C')) / max(len(seq), 1)

def fix_hp(seq, max_hp=MAX_HP):
    s = list(seq); run = 1
    for i in range(1, len(s)):
        if s[i] == s[i-1]:
            run += 1
        else:
            run = 1
        if run > max_hp:
            s[i] = random.choice([b for b in BASES if b != s[i]])
            run = 1
    return ''.join(s)

def fix_gc(seq, gc_min=GC_MIN, gc_max=GC_MAX):
    s = list(seq)
    for _ in range(300):
        g = gc(''.join(s))
        if gc_min <= g <= gc_max:
            break
        if g < gc_min:
            at_idx = [i for i, b in enumerate(s) if b in 'AT']
            if not at_idx: break
            s[random.choice(at_idx)] = random.choice('GC')
        else:
            gc_idx = [i for i, b in enumerate(s) if b in 'GC']
            if not gc_idx: break
            s[random.choice(gc_idx)] = random.choice('AT')
    return ''.join(s)

def constrain(seq):
    seq = fix_hp(seq)
    seq = fix_gc(seq)
    seq = fix_hp(seq)
    return seq

# ---------------- Strand indexing (FIXED) ----------------
def index_width_nt(n_strands):
    """Nucleotides needed to address n_strands (2 bits/nt), min 4nt."""
    if n_strands <= 1:
        return 4
    bits_needed = max(8, math.ceil(math.log2(n_strands)))
    if bits_needed % 2:
        bits_needed += 1
    nt = bits_needed // 2
    return max(nt, 4)

def _int_to_dna(val, nt):
    bits = format(val, f'0{nt*2}b')
    return ''.join(B2D[bits[i:i+2]] for i in range(0, len(bits), 2))

def _dna_to_int(seq):
    bits = ''.join(D2B.get(b, '00') for b in seq)
    return int(bits, 2) if bits else 0

def _parity_nt(index_seq):
    """1-nt parity over the index block: lets the decoder detect
    (not silently ignore) a corrupted address."""
    v = 0
    for b in index_seq:
        v ^= D2B.get(b, '00') and int(D2B[b], 2)
    return B2D[format(v & 3, '02b')]

def to_strands(dna, payload_len):
    """Split dna into strands of a fixed payload length, each carrying
    an index block (width chosen to fit the actual strand count) plus
    a 1nt parity guard on that index."""
    n_strands = math.ceil(len(dna) / payload_len) if dna else 1
    idx_nt = index_width_nt(n_strands)
    strands = []
    for i in range(0, len(dna), payload_len) if dna else [0]:
        chunk = dna[i:i+payload_len].ljust(payload_len, 'A')
        idx = i // payload_len
        idx_seq = _int_to_dna(idx, idx_nt)
        parity = _parity_nt(idx_seq)
        strands.append(idx_seq + parity + chunk)
    return strands, idx_nt

def parse_strand_header(strand, idx_nt):
    idx_seq = strand[:idx_nt]
    parity = strand[idx_nt:idx_nt+1]
    payload = strand[idx_nt+1:]
    ok = (parity == _parity_nt(idx_seq))
    return _dna_to_int(idx_seq), ok, payload

# ---------------- Markers (search-based, tolerant) ----------------
def add_markers(seq, every=40, marker=MARKER):
    out = ''
    for i in range(0, len(seq), every):
        out += seq[i:i+every]
        if i + every < len(seq):
            out += marker
    return out

def strip_markers_by_search(seq, marker=MARKER, every=40, search_radius=3, max_mismatch=1):
    """Remove markers by checking near their EXPECTED stride position
    (every `every` data-nt), not by scanning the whole sequence for the
    marker string. A global scan is unsafe: an 8-mer as short as
    'ACGTACGT' collides with real (compressed/ECC) data content purely
    by chance every ~40-60kb, and one false hit desyncs everything after
    it. Anchoring the search to the expected offset (with a small
    +/- radius to tolerate an upstream indel) avoids that failure mode
    and still gives some indel tolerance for the noisy/blind-decode case.
    """
    mlen = len(marker)
    out = []
    i = 0
    while i < len(seq):
        out.append(seq[i:i+every])
        i += every
        if i >= len(seq):
            break
        best_off, best_mm = 0, mlen + 1
        for off in range(-search_radius, search_radius + 1):
            start = i + off
            window = seq[start:start+mlen]
            if len(window) != mlen:
                continue
            mm = sum(1 for a, b in zip(window, marker) if a != b)
            if mm < best_mm:
                best_mm, best_off = mm, off
        if best_mm <= max_mismatch:
            i += best_off + mlen
        else:
            # marker not found near expected slot (heavy corruption) -
            # assume no shift rather than derailing the whole strand
            i += mlen
    return ''.join(out)

# ---------------- Error-Correcting Code (Reed-Solomon) ----------------
import reedsolo

def rs_protect(data: bytes, ecc_bytes=32):
    """Add Reed-Solomon parity bytes so the payload survives both the
    constraint-fixer's own edits and downstream sequencing errors.
    RSCodec works in blocks of <=255 bytes total, so we chunk."""
    rsc = reedsolo.RSCodec(ecc_bytes)
    block = 255 - ecc_bytes
    out = bytearray()
    for i in range(0, len(data), block):
        chunk = data[i:i+block]
        out += rsc.encode(chunk)
    return bytes(out), block, ecc_bytes

def rs_recover(data: bytes, block, ecc_bytes, n_original_blocks):
    rsc = reedsolo.RSCodec(ecc_bytes)
    full_block = block + ecc_bytes
    out = bytearray()
    errors_corrected = 0
    for i in range(n_original_blocks):
        chunk = data[i*full_block:(i+1)*full_block]
        try:
            decoded, _, _ = rsc.decode(bytes(chunk))
            out += decoded
        except reedsolo.ReedSolomonError:
            # Uncorrectable block: keep raw (still gives partial recovery)
            out += bytes(chunk[:block])
            errors_corrected = -1
    return bytes(out), errors_corrected
