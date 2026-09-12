[README.md](https://github.com/user-attachments/files/32148372/README.md)
# DNA Storage Pipeline — Phase 1 (fixed, working, split into separa# DNA Nanopore Error Correction Pipeline

## Quick start - train once, test anytime

**Train the model (run once, on your big dataset):**
```bash
python test_pipeline.py --input clean_dataset_200kb.txt --outdir run_final --train \
    --depth 15x --identity "97,99.5,1" --chunk-bytes 24 --ecc-bytes 16
```

**Test any file afterward (reuses the trained model, no retraining):**
```bash
python test_pipeline.py --input test_text.txt --outdir run_test_text \
    --model-dir run_final/models --model random_forest \
    --identity "97,99.5,1" --chunk-bytes 24 --ecc-bytes 16
```

---

A DNA-storage encode/decode pipeline that simulates Oxford Nanopore
sequencing noise (via Badread), detects and corrects errors using a
trained ML classifier, and reconstructs the original file using
erasure-aware Reed-Solomon decoding.

---

## 1. How the pipeline works

```
 original file
      |
      v
 [1_encode.py]  ->  splits file into small chunks, compresses + Reed-Solomon
      |              protects EACH chunk independently, maps to DNA, one
      |              chunk = one strand  ->  clean.fasta + meta.json
      v
 [2_simulate_noise.py]  ->  Badread simulates nanopore sequencing noise on
      |                      each strand  ->  noisy.fastq (multiple noisy
      |                      reads per strand)
      v
 [3_train_ml.py]  ->  aligns noisy reads back to their clean strand,
      |                labels every position (clean / substitution /
      |                insertion / deletion), trains Logistic Regression /
      |                Random Forest / XGBoost classifiers on ~33 features
      |                per position (sequence context, GC%, homopolymer
      |                run, quality score, alignment gap structure)
      v
 [4_evaluate.py]  ->  for each strand: builds a multi-read CONSENSUS,
      |                runs it through the trained classifier, corrects
      |                predicted errors, flags low-confidence positions as
      |                Reed-Solomon ERASURES (correctable at double the
      |                rate of blind errors), decodes each chunk
      |                independently, reassembles the file
      v
 recovered_file.bin (+ full report.json)
```

### Why "chunked" instead of one continuous compressed stream

The original design compressed the *entire file* as one Huffman+LZ77
stream split across strands. That has a fatal property: if even ONE
strand fails to recover perfectly, every symbol boundary after it
shifts, corrupting the rest of the file regardless of how many other
strands were fine. **Whole-file recovery was literally 0% or 100%, never
in between.**

The current design compresses and Reed-Solomon-protects each small chunk
**independently**, with each chunk mapped to exactly one strand. A bad
strand now only costs its own ~24 bytes, not the rest of the file — so
recovery scales **proportionally** with how many strands you actually
get right.

### Why erasure-aware decoding matters

Reed-Solomon has two correction modes:
- **Blind error correction** (position unknown): fixes up to `ecc_bytes / 2` bad bytes
- **Erasure correction** (position known): fixes up to `ecc_bytes` bad bytes — **double**

The ML model's per-position predictions (substitution / deletion) are
used as erasure hints, capped and ranked by confidence so the erasure
budget is never overrun by false positives.

---

## 2. Files in this repo

| File | Purpose |
|---|---|
| `core_chunked.py` | Compression (LZ77 + shared Huffman codebook), Reed-Solomon protect/recover, DNA<->binary mapping, biological constraints (GC%, homopolymer limits) |
| `pipeline_io_chunked.py` | Ties chunking + compression + RS + DNA mapping together: `encode_bytes_chunked()` / `decode_strands_chunked()` |
| `align_ml.py` | Needleman-Wunsch alignment, per-position feature extraction (33 features), multi-read consensus voting, ML-based correction, erasure position scoring |
| `1_encode.py` | STEP 1 - encode a file into `clean.fasta` + `meta.json` |
| `2_simulate_noise.py` | STEP 2 - run Badread to produce `noisy.fastq` |
| `3_train_ml.py` | STEP 3 - train Logistic Regression / Random Forest / XGBoost on labeled alignment data |
| `4_evaluate.py` | STEP 4 - consensus + ML-correct + erasure-decode + report |
| `test_pipeline.py` | **All-in-one**: encode -> noise -> (train OR reuse existing model) -> decode -> report, in a single command |
| `check_feature_importance.py` | Diagnostic: shows which of the 33 features the trained Random Forest actually relies on |

---

## 3. Setup

```bash
pip install -r requirements.txt
pip install git+https://github.com/rrwick/Badread.git
```

On Windows, if `edlib` (a Badread dependency) fails to build, install the
"Desktop development with C++" workload via the Visual Studio Installer
first.

---

## 4. Commands

### 4a. Recommended - single command via `test_pipeline.py`

**Train once, on your main dataset:**
```bash
python test_pipeline.py --input clean_dataset_200kb.txt --outdir run_final --train \
    --depth 15x --identity "97,99.5,1" --chunk-bytes 24 --ecc-bytes 16
```

**Test any other file afterward, reusing that trained model (no retraining):**
```bash
python test_pipeline.py --input test_text.txt --outdir run_test_text \
    --model-dir run_final/models --model random_forest \
    --identity "97,99.5,1" --chunk-bytes 24 --ecc-bytes 16
```

Key flags:
| Flag | Meaning |
|---|---|
| `--chunk-bytes N` | original file bytes per strand (default 24) |
| `--ecc-bytes N` | Reed-Solomon parity bytes per chunk (default 16) |
| `--depth 15x` | simulated reads per strand |
| `--identity "97,99.5,1"` | Badread mean,max,stdev identity - see note below, this matters more than any other setting |
| `--max-consensus-reads N` | reads combined per strand via majority vote before ML correction (default 8) |
| `--train` | train a fresh model on this run's own noise |
| `--model-dir DIR --model NAME` | reuse an existing trained model instead of retraining |

### 4b. Manual step-by-step (equivalent, if you want to inspect intermediate files)

```bash
python 1_encode.py --input myfile.txt --outdir run1 --chunk-bytes 24 --ecc-bytes 16
python 2_simulate_noise.py --rundir run1 --depth 15x --identity "97,99.5,1"
python 3_train_ml.py --rundir run1
python 4_evaluate.py --rundir run1 --model random_forest --original myfile.txt
```

---

## 5. IMPORTANT: the `--identity` setting matters more than anything else

Every experiment in this project converged on one finding: **Reed-Solomon
recovery has a hard, information-theoretic ceiling that no amount of
chunk-size or ECC tuning can move past.**

Because 4 nucleotides pack into 1 byte, a single wrong base corrupts the
*entire* byte - so a "modest" 10% per-base error rate becomes roughly
**1-(1-0.10)^4 = ~34% byte-level corruption**, and a 27% per-base error
rate becomes **~72% byte corruption**. Since noise hits the ECC bytes at
the same rate as the data bytes, adding more ECC does not help once
corruption exceeds a fixed ratio - it's not a tunable bug, it's basic
coding theory.

Measured relationship (chunk_bytes=24, ecc_bytes matched to error rate,
perfect erasure knowledge / theoretical upper bound):

| Raw per-base error rate | Whole-file byte match achieved |
|---|---|
| 3% | **52.2%** |
| 5% | 30.8% |
| 8% | 14.8% |
| 10% (Badread's noisy default `90,98,4`) | 12.4% |
| 15% | 6.1% |

`--identity "90,98,4"` (Badread's default) simulates a **worst-case,
very noisy** nanopore run. `--identity "97,99.5,1"` represents a
realistic **high-quality nanopore / duplex sequencing** scenario, which
modern basecallers routinely achieve. Use the latter to let the
pipeline's actual design (chunking + consensus + ML + erasure decoding)
demonstrate its value, rather than fighting noise that is mathematically
unrecoverable at reasonable redundancy overhead.

---

## 6. Known limitations / open problems for future work

1. **Frameshift risk**: if the ML model misclassifies an insertion or
   deletion (false positive/negative), the corrected sequence's length
   can drift out of alignment with the true byte grid, corrupting
   everything downstream *within that chunk*. Indel detection is
   currently ~1.00 precision/recall, so this is rare but not zero.
2. **Substitution-vs-clean detection ceiling**: unlike indels (which
   leave a structural alignment gap), substitutions have no direct
   structural signal. Even with quality-score features, precision/recall
   plateaus around 0.65-0.70. This is the primary target for the CNN /
   LSTM / Transformer models planned in the main project.
3. **Consensus limited by correlated nanopore errors**: real nanopore
   errors are context-dependent (specific k-mers / homopolymer runs are
   *systematically* hard to call), so multiple reads of the same strand
   often share the same error at the same position - majority voting
   helps less than it would against purely random, independent noise.
4. **1 XOR-parity bit protects the strand index** - sufficient to detect
   (not correct) a corrupted index and prevent misrouting a strand's
   data into the wrong chunk, but weak compared to a real ECC on the
   index itself.

---

## 7. Roadmap toward the main project (CNN / LSTM / Transformer)

The classical ML baseline (Logistic Regression / Random Forest /
XGBoost) in this repo exists specifically to give the deep-learning
models something concrete to beat:

- **CNN**: 1D convolutions over local sequence windows should capture
  the k-mer/homopolymer error patterns this project's hand-crafted
  +/-2-base one-hot features only approximate.
- **LSTM**: sequence-aware modeling of insertion/deletion patterns,
  targeting the frameshift-risk limitation above.
- **Transformer**: full-strand attention-based reconstruction, ideally
  replacing the consensus + classify + erasure-decode chain with a
  single learned correction step.

See `EXPERIMENT_REPORT.pdf` for the full history of what was tried,
what worked, what didn't, and why - intended as source material for the
main project's research paper.
te files)

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
