"""
STEP 1 — ENCODE a file into DNA strands.

Run:
    python 1_encode.py --input myfile.pdf --outdir run1

Produces (inside --outdir):
    clean.fasta   -> the DNA strands (what you'd "synthesize")
    meta.json     -> everything decode.py needs to reverse the process

This is the only file that needs core.py and pipeline_io.py to encode.
"""
import argparse, json, os
from pipeline_io import encode_bytes


def write_fasta(strands, path):
    with open(path, 'w') as f:
        for i, s in enumerate(strands):
            f.write(f'>strand_{i:06d}\n{s}\n')


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--input', required=True, help='file to encode')
    ap.add_argument('--outdir', required=True, help='folder to write clean.fasta + meta.json into')
    ap.add_argument('--payload-nt', type=int, default=100, help='data nucleotides per strand')
    ap.add_argument('--ecc-bytes', type=int, default=48, help='Reed-Solomon parity bytes per 255-byte block')
    ap.add_argument('--marker-every', type=int, default=40, help='insert a marker every N nt')
    args = ap.parse_args()

    os.makedirs(args.outdir, exist_ok=True)

    with open(args.input, 'rb') as f:
        raw = f.read()

    print(f"Input file   : {args.input} ({len(raw):,} bytes)")
    strands, meta = encode_bytes(raw, payload_nt=args.payload_nt,
                                  ecc_bytes=args.ecc_bytes, marker_every=args.marker_every)

    fasta_path = os.path.join(args.outdir, 'clean.fasta')
    meta_path = os.path.join(args.outdir, 'meta.json')
    write_fasta(strands, fasta_path)
    json.dump(meta, open(meta_path, 'w'), indent=2)

    print(f"Strands      : {len(strands):,}")
    print(f"Strand length: {len(strands[0])} nt")
    print(f"Index width  : {meta['idx_nt']} nt (supports up to {4**meta['idx_nt']:,} strands)")
    print(f"Written -> {fasta_path}")
    print(f"Written -> {meta_path}")

    # Sanity check right here: does this round-trip with ZERO sequencing
    # errors? If this fails, nothing downstream can possibly work, so we
    # check it immediately instead of finding out three steps later.
    from pipeline_io import decode_strands
    recovered, rs_err = decode_strands(strands, meta)
    ok = recovered == raw
    print(f"Zero-error round-trip check: {'PASSED' if ok else 'FAILED'}")
    if not ok:
        print("WARNING: even with no sequencing errors, decode did not reproduce "
              "the input exactly. Try increasing --ecc-bytes.")


if __name__ == '__main__':
    main()
