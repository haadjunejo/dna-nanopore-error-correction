"""Full encode/decode using core.py, with RS-ECC protecting the payload
before it ever touches the constraint-fixer or DNA mapping."""
import json, math
import core

def encode_bytes(raw: bytes, payload_nt=80, ecc_bytes=48, marker_every=40):
    bits_raw, cmeta = core.compress_file_bytes(raw)
    # bits -> bytes for RS (pad to byte boundary)
    pad = (-len(bits_raw)) % 8
    bits_padded = bits_raw + '0'*pad
    byte_data = int(bits_padded, 2).to_bytes(len(bits_padded)//8, 'big') if bits_padded else b''
    protected, block, ecc = core.rs_protect(byte_data, ecc_bytes=ecc_bytes)
    n_blocks = math.ceil(len(byte_data)/block) if byte_data else 0
    protected_bits = ''.join(format(b, '08b') for b in protected)
    dna = core.bits_to_dna(protected_bits)
    dna_c = core.constrain(dna)
    marked = core.add_markers(dna_c, every=marker_every)
    strands, idx_nt = core.to_strands(marked, payload_len=payload_nt)
    meta = {
        'cmeta': cmeta, 'bits_len': len(bits_raw), 'pad': pad,
        'byte_len': len(byte_data), 'block': block, 'ecc_bytes': ecc,
        'n_blocks': n_blocks, 'dna_len': len(dna_c), 'marked_len': len(marked),
        'idx_nt': idx_nt, 'marker_every': marker_every, 'n_strands': len(strands),
    }
    return strands, meta

def decode_strands(strands, meta):
    idx_nt = meta['idx_nt']
    parsed = [core.parse_strand_header(s, idx_nt) for s in strands]
    good = {idx: payload for idx, ok, payload in parsed if ok}
    n_expected = meta['n_strands']
    full = ''.join(good.get(i, 'A'*(len(strands[0])-idx_nt-1)) for i in range(n_expected))
    full = full[:meta['marked_len']]
    raw_dna = core.strip_markers_by_search(full, every=meta['marker_every'])
    raw_dna = raw_dna[:meta['dna_len']]
    protected_bits = core.dna_to_bits(raw_dna)
    protected_bytes = int(protected_bits, 2).to_bytes(len(protected_bits)//8, 'big') if protected_bits else b''
    byte_data, err = core.rs_recover(protected_bytes, meta['block'], meta['ecc_bytes'], meta['n_blocks'])
    byte_data = byte_data[:meta['byte_len']]
    bits_padded = ''.join(format(b, '08b') for b in byte_data)
    bits_raw = bits_padded[:meta['bits_len']]
    original = core.decompress_bits(bits_raw, meta['cmeta'])
    return original, err
