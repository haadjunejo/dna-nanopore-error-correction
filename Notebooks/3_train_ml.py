"""
STEP 3 — TRAIN the ML error-type classifiers (Logistic Regression,
Random Forest, XGBoost) on real Badread noise.

Run:
    python 3_train_ml.py --rundir run1

Reads:  run1/clean.fasta, run1/noisy.fastq
Writes: run1/models/*.pkl, run1/models/training_report.json
"""
import argparse, json, os, pickle, re
from collections import defaultdict

import numpy as np
from sklearn.model_selection import train_test_split
from sklearn.linear_model import LogisticRegression
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import accuracy_score, f1_score, classification_report
import xgboost as xgb

from align_ml import build_dataset, LABEL_NAMES


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
    ap.add_argument('--samples-per-class', type=int, default=6000,
                     help='alignment columns sampled per error class for training')
    args = ap.parse_args()

    strands = read_fasta(os.path.join(args.rundir, 'clean.fasta'))
    noisy_reads = read_fastq(os.path.join(args.rundir, 'noisy.fastq'))

    pairs = [(strands[sid], r['seq']) for sid, rl in noisy_reads.items()
             if sid < len(strands) for r in rl]
    print(f"(clean, noisy) read pairs available: {len(pairs):,}")
    if not pairs:
        raise RuntimeError("No reads found. Did 2_simulate_noise.py run successfully?")

    X, y, class_counts = build_dataset(pairs, samples_per_class=args.samples_per_class, seed=42)
    print("Alignment columns available per class:",
          {LABEL_NAMES[k]: v for k, v in class_counts.items()})
    print("Feature matrix:", X.shape)

    Xtr, Xte, ytr, yte = train_test_split(X, y, test_size=0.2, random_state=42, stratify=y)

    models = {
        'logistic_regression': LogisticRegression(max_iter=2000, class_weight='balanced'),
        'random_forest': RandomForestClassifier(n_estimators=200, class_weight='balanced',
                                                 n_jobs=-1, random_state=42),
        'xgboost': xgb.XGBClassifier(n_estimators=200, max_depth=6,
                                      eval_metric='mlogloss', random_state=42),
    }

    os.makedirs(os.path.join(args.rundir, 'models'), exist_ok=True)
    report = {}
    for name, m in models.items():
        m.fit(Xtr, ytr)
        pred = m.predict(Xte)
        acc = accuracy_score(yte, pred)
        f1m = f1_score(yte, pred, average='macro')
        report[name] = {'accuracy': float(acc), 'macro_f1': float(f1m)}
        print(f"{name:22s} accuracy={acc:.4f}  macro-F1={f1m:.4f}")
        with open(os.path.join(args.rundir, 'models', f'{name}.pkl'), 'wb') as f:
            pickle.dump(m, f)

    best_name = max(report, key=lambda k: report[k]['macro_f1'])
    report['best_model'] = best_name
    print(f"\nBest model: {best_name}")
    print("\n" + classification_report(
        yte, models[best_name].predict(Xte),
        target_names=[LABEL_NAMES[i] for i in range(4)], zero_division=0))

    json.dump(report, open(os.path.join(args.rundir, 'models', 'training_report.json'), 'w'), indent=2)
    print(f"\nModels saved -> {args.rundir}/models/*.pkl")
    print(f"Report saved -> {args.rundir}/models/training_report.json")


if __name__ == '__main__':
    main()
