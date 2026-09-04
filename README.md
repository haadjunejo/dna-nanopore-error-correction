# DNA Storage Pipeline — Phase 1 (fixed, working, split into separate files)

Each step is its own file. Run them **in order** — each one reads files
that the previous step wrote. Everything for one experiment lives in a
single `--rundir` folder (e.g. `run1/`), so you can run several
experiments side by side just by using different folder names.

## Files

| File | Job | Reads | Writes |
|---|---|---|---|
| `core.py` | Shared building blocks (compression, Reed-Solomon ECC, DNA mapping, GC/homopolymer constraints, strand indexing, markers). Not run directly. | — | — |
| `pipeline_io.py` | Full encode/decode logic, built on `core.py`. Not run directly. | — | — |
| `align_ml.py` | Alignment (Needleman-Wunsch) + feature extraction for ML. Not run directly. | — | — |
| `1_encode.py` | Turns your file into DNA strands | your input file | `clean.fasta`, `meta.json` |
| `2_simulate_noise.py` | Simulates real nanopore sequencing errors | `clean.fasta` | `noisy.fastq` |
| `3_train_ml.py` | Trains the error-detector models | `clean.fasta`, `noisy.fastq` | `models/*.pkl` |
| `4_evaluate.py` | Corrects noisy reads with the trained model, measures recovery | everything above | `report.json`, `recovered_file.bin` |

You never run `core.py`, `pipeline_io.py`, or `align_ml.py` yourself —
the numbered scripts import them.

## 0. One-time setup

```bash
pip install -r requirements.txt
pip install git+https://github.com/rrwick/Badread.git
```

Badread is the tool that simulates real Oxford Nanopore-style sequencing
errors (substitutions, insertions, deletions). It's a separate install
because it comes from GitHub, not PyPI.

## 1. Encode your file into DNA

```bash
python 1_encode.py --input myfile.pdf --outdir run1
```

This compresses your file, protects it with Reed-Solomon error-correcting
code, converts it to a DNA sequence, enforces GC-content and homopolymer
limits, adds markers, and splits it into strands with a protected index.

It also immediately checks that decoding gives back your exact original
file with **zero** sequencing errors introduced — if that check fails,
something is wrong before any noise is even involved, and you'll see
`FAILED` printed instead of `PASSED`. Don't move on to step 2 until this
says `PASSED`.

Useful options:
- `--payload-nt 100` — data nucleotides per strand (bigger = fewer strands, more risk per strand)
- `--ecc-bytes 48` — how much Reed-Solomon protection per 255-byte block (bigger = more robust, more overhead)

## 2. Simulate nanopore sequencing

```bash
python 2_simulate_noise.py --rundir run1
```

Runs Badread against `run1/clean.fasta` and writes `run1/noisy.fastq` —
this is what your data would look like after actually being sequenced.

Useful options:
- `--identity "90,98,4"` — mean/max/stdev read identity (lower = noisier reads)
- `--depth 3x` — how many reads per strand (like sequencing coverage)

## 3. Train the ML error detectors

```bash
python 3_train_ml.py --rundir run1
```

Aligns every noisy read back to its own clean strand (Needleman-Wunsch),
labels every position as clean/substitution/insertion/deletion, extracts
features, and trains three classifiers: Logistic Regression, Random
Forest, and XGBoost. Saves all three plus a report of which did best.

This step can take a minute or two depending on how many strands/reads
you generated in step 2.

## 4. Evaluate recovery

```bash
python 4_evaluate.py --rundir run1 --model xgboost --original myfile.pdf
```

`--model` picks which trained classifier to use (`logistic_regression`,
`random_forest`, or `xgboost` — check `run1/models/training_report.json`
to see which scored best).

`--original` is optional but recommended — pass your original input file
again here and it will also report whole-file byte match and write out
`run1/recovered_file.bin` so you can inspect what actually came back.

This prints (and saves to `run1/report.json`):
- **Classical baseline**: what you get with no alignment/ML at all (Reed-Solomon only). Expect this to be very poor — that's the real, measured reason DNA storage research needs more than classical codes once deletions are involved.
- **ML-corrected (oracle)**: uses the true reference to fill in detected deletions. This measures how good *detection* is — it's an upper bound, not a blind result.
- **ML-corrected (honest blind)**: never looks at the reference. This is the fair number for "what can this system currently do on its own."

## Running everything for a new file

```bash
python 1_encode.py --input newfile.txt --outdir run2
python 2_simulate_noise.py --rundir run2
python 3_train_ml.py --rundir run2
python 4_evaluate.py --rundir run2 --model xgboost --original newfile.txt
```

## What's next (Phase 2)

Steps 1–2 (encoding, noise simulation) don't need to change. To bring in
CNN/LSTM/Transformer models, replace step 3's feature extraction +
classical models with your deep learning models, trained on the exact
same `(clean_strand, noisy_read)` pairs that `3_train_ml.py` builds —
that keeps the comparison against this baseline fair.

## Known limitations (stated honestly, not hidden)

- Badread generates multiple reads per strand (`--depth`), but this
  pipeline currently only uses the *first* read per strand for
  correction — there's no majority-vote consensus across reads yet.
  That's a clear, relatively easy next improvement.
- Detecting *that* a deletion happened is learnable from a single noisy
  read; recovering *which* base was deleted is not — the base was never
  observed. Only "oracle" mode can do that, and it's not a blind result.
  A Transformer that reasons across the whole strand is the real fix.
- The interior marker scheme is a research prototype: a short repeated
  marker can coincidentally match real data content in large files.
