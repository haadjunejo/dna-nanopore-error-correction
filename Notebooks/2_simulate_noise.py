"""
STEP 2 — SIMULATE nanopore sequencing noise on the encoded strands.

Requires Badread:
    pip install git+https://github.com/rrwick/Badread.git

Run:
    python 2_simulate_noise.py --rundir run1

Reads:  run1/clean.fasta
Writes: run1/noisy.fastq
"""
import argparse, json, os, subprocess


def run_badread(fasta_in, fastq_out, strand_len, identity='90,98,4', depth='3x'):
    cmd = [
        'badread', 'simulate',
        '--reference', fasta_in,
        '--quantity', depth,
        '--identity', identity,
        '--length', f'{strand_len + 10},5',
        '--junk_reads', '0',
        '--random_reads', '0',
        '--chimeras', '0',
    ]
    print('Running:', ' '.join(cmd))
    with open(fastq_out, 'w') as out_f:
        r = subprocess.run(cmd, stdout=out_f, stderr=subprocess.PIPE, text=True)
    if r.returncode != 0:
        raise RuntimeError(r.stderr[-2000:])
    print(r.stderr.strip().splitlines()[-1] if r.stderr.strip() else '')


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--rundir', required=True, help='folder created by 1_encode.py')
    ap.add_argument('--identity', default='90,98,4', help='badread "mean,max,stdev" read identity')
    ap.add_argument('--depth', default='3x', help='coverage, e.g. 3x = ~3 reads per strand')
    args = ap.parse_args()

    fasta_in = os.path.join(args.rundir, 'clean.fasta')

    # read actual strand length straight from the fasta (simplest source of truth)
    with open(fasta_in) as f:
        f.readline()
        strand_len = len(f.readline().strip())

    fastq_out = os.path.join(args.rundir, 'noisy.fastq')
    run_badread(fasta_in, fastq_out, strand_len, identity=args.identity, depth=args.depth)
    print(f"Written -> {fastq_out}")


if __name__ == '__main__':
    main()
