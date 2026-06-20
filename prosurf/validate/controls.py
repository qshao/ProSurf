import numpy as np
from prosurf.metric.engine_a import score_charges_a

def protein_zwit_score(charges, cfg):
    scores = score_charges_a(charges, cfg)
    return max((s.z for s in scores), default=0.0)

def auroc(pos_scores, neg_scores):
    pos = np.asarray(pos_scores); neg = np.asarray(neg_scores)
    if len(pos) == 0 or len(neg) == 0:
        return float("nan")
    wins = sum((p > n) + 0.5 * (p == n) for p in pos for n in neg)
    return float(wins / (len(pos) * len(neg)))
