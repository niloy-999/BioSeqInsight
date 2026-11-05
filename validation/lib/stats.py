"""Statistics used by the validation study (standard library only).

Everything here is deliberately simple and independently testable:
Wilson intervals for proportions, an exact McNemar test for paired
decisions, cluster-bootstrap intervals for derived metrics (clustered by
query accession, because several candidates can share one query), and a
rank-based AUC.
"""

from __future__ import annotations

import math
import random
from typing import Callable, Iterable, Sequence


def wilson(k: int, n: int, z: float = 1.959963984540054) -> tuple[float, float, float]:
    """Point estimate and Wilson score interval for k successes out of n."""
    if n <= 0:
        return (float("nan"), float("nan"), float("nan"))
    p = k / n
    denom = 1 + z * z / n
    centre = (p + z * z / (2 * n)) / denom
    half = z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / denom
    return (p, max(0.0, centre - half), min(1.0, centre + half))


def mcnemar_exact(b: int, c: int) -> float:
    """Two-sided exact McNemar p-value from the discordant counts b and c."""
    n = b + c
    if n == 0:
        return 1.0
    k = min(b, c)
    tail = sum(math.comb(n, i) for i in range(0, k + 1)) / (2 ** n)
    return min(1.0, 2 * tail)


def confusion(truth: Sequence[bool], decision: Sequence[bool]) -> dict[str, int]:
    """2x2 counts where *positive* means "the candidate is correct"."""
    tp = fp = tn = fn = 0
    for t, d in zip(truth, decision):
        if d and t:
            tp += 1
        elif d and not t:
            fp += 1
        elif not d and not t:
            tn += 1
        else:
            fn += 1
    return {"tp": tp, "fp": fp, "tn": tn, "fn": fn}


def metrics(c: dict[str, int]) -> dict[str, float]:
    """Precision, recall, specificity, F1, FAR, FNR, MCC, accuracy."""
    tp, fp, tn, fn = c["tp"], c["fp"], c["tn"], c["fn"]

    def div(a: float, b: float) -> float:
        return a / b if b else float("nan")

    precision = div(tp, tp + fp)
    recall = div(tp, tp + fn)
    specificity = div(tn, tn + fp)
    f1 = div(2 * precision * recall, precision + recall) if tp else (0.0 if (fp or fn) else float("nan"))
    denom = math.sqrt((tp + fp) * (tp + fn) * (tn + fp) * (tn + fn))
    mcc = div(tp * tn - fp * fn, denom)
    return {
        "precision": precision,
        "recall": recall,
        "specificity": specificity,
        "f1": f1,
        "false_association_rate": div(fp, tp + fp),  # = 1 - precision
        "false_negative_rate": div(fn, tp + fn),
        "false_positive_rate": div(fp, fp + tn),
        "mcc": mcc,
        "accuracy": div(tp + tn, tp + fp + tn + fn),
    }


def auc(scores_pos: Sequence[float], scores_neg: Sequence[float]) -> float:
    """Area under the ROC curve via the Mann-Whitney U statistic (ties = 0.5)."""
    if not scores_pos or not scores_neg:
        return float("nan")
    ranked = sorted([(s, 1) for s in scores_pos] + [(s, 0) for s in scores_neg])
    rank_sum_pos = 0.0
    i = 0
    n = len(ranked)
    while i < n:
        j = i
        while j + 1 < n and ranked[j + 1][0] == ranked[i][0]:
            j += 1
        avg_rank = (i + j) / 2 + 1
        rank_sum_pos += avg_rank * sum(1 for k in range(i, j + 1) if ranked[k][1] == 1)
        i = j + 1
    n1, n0 = len(scores_pos), len(scores_neg)
    return (rank_sum_pos - n1 * (n1 + 1) / 2) / (n1 * n0)


def cluster_bootstrap(
    rows: Sequence[dict],
    cluster_key: str,
    statistic: Callable[[list[dict]], float],
    n_boot: int = 2000,
    seed: int = 20260927,
    alpha: float = 0.05,
) -> tuple[float, float, float]:
    """Percentile bootstrap of ``statistic`` resampling whole clusters."""
    clusters: dict[str, list[dict]] = {}
    for row in rows:
        clusters.setdefault(str(row[cluster_key]), []).append(row)
    keys = sorted(clusters)
    point = statistic(list(rows))
    if len(keys) < 2:
        return (point, float("nan"), float("nan"))
    rng = random.Random(seed)
    values: list[float] = []
    for _ in range(n_boot):
        sample: list[dict] = []
        for _k in keys:
            sample.extend(clusters[rng.choice(keys)])
        v = statistic(sample)
        if not (isinstance(v, float) and math.isnan(v)):
            values.append(v)
    if not values:
        return (point, float("nan"), float("nan"))
    values.sort()
    lo = values[int((alpha / 2) * len(values))]
    hi = values[min(len(values) - 1, int((1 - alpha / 2) * len(values)))]
    return (point, lo, hi)


def pearson(x: Sequence[float], y: Sequence[float]) -> float:
    n = len(x)
    if n < 2:
        return float("nan")
    mx, my = sum(x) / n, sum(y) / n
    sxx = sum((a - mx) ** 2 for a in x)
    syy = sum((b - my) ** 2 for b in y)
    sxy = sum((a - mx) * (b - my) for a, b in zip(x, y))
    return sxy / math.sqrt(sxx * syy) if sxx and syy else float("nan")


def median(values: Iterable[float]) -> float:
    v = sorted(values)
    if not v:
        return float("nan")
    m = len(v) // 2
    return v[m] if len(v) % 2 else (v[m - 1] + v[m]) / 2


__all__ = ["wilson", "mcnemar_exact", "confusion", "metrics", "auc", "cluster_bootstrap", "pearson", "median"]
