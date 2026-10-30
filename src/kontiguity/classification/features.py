"""The 12 'simple' multi-scale features of a contact profile (identical to the training code).

  log10_depth      log10 of the total trans contacts with the host
  exvar_<s>        log mean squared Pearson residual after aggregation in windows of s
                   (5, 25, 100, 500 kb): how far the profile departs from the expectation at that scale
  maxz_<s>         log(1 + highest positive residual) at that scale: strongest peak
  top1pct_oe_25kb  log mean obs/exp of the 1 % best 25-kb windows
  max_window_frac  share of all contacts falling in the single best 25-kb window (cis anchoring)
  chrom_end_oe     log obs/exp in the 25 kb at both ends of each host chromosome
"""
import numpy as np

EPS = 1e-6
NAMES = ["log10_depth", "exvar_5kb", "maxz_5kb", "exvar_25kb", "maxz_25kb", "top1pct_oe_25kb", "max_window_frac",
         "exvar_100kb", "maxz_100kb", "exvar_500kb", "maxz_500kb", "chrom_end_oe"]


def hand_feats(R, w, mask, chrom):
    R = np.atleast_2d(np.asarray(R, float)) * mask
    tot = R.sum(1, keepdims=True)
    E = tot * w[None]
    f = [np.log10(tot[:, 0] + 1)]
    n = R.shape[1]
    for k in (1, 5, 20, 100):
        m = n // k
        Rk = R[:, :m * k].reshape(len(R), m, k).sum(2)
        Ek = E[:, :m * k].reshape(len(R), m, k).sum(2)
        ok = Ek > 0
        z = np.where(ok, (Rk - Ek) / np.sqrt(Ek + EPS), 0)
        f += [np.log((z ** 2).sum(1) / np.maximum(ok.sum(1), 1) + EPS), np.log(np.maximum(z.max(1), 0) + 1)]
        if k == 5:
            q = np.sort(np.where(ok, Rk / (Ek + EPS), 0), 1)[:, -max(1, m // 100):]
            f.append(np.log(q.mean(1) + EPS))
            f.append(Rk.max(1) / (tot[:, 0] + 1))
    ends = np.zeros(n, bool)
    b = np.r_[0, np.where(np.diff(chrom))[0] + 1, n]
    for s, e in zip(b[:-1], b[1:]):
        ends[s:s + 5] = True
        ends[max(s, e - 5):e] = True
    ends &= mask > 0
    f.append(np.log((R[:, ends].sum(1) + 1) / (E[:, ends].sum(1) + 1)))
    return np.stack(f, 1)


def residual(r, w, mask):
    """Freeman-Tukey residual vs expectation (~N(0,1) under the null at any depth)."""
    r = np.where(mask, r, 0.0)
    e = r.sum() * w
    z = np.sqrt(r) + np.sqrt(r + 1) - np.sqrt(4 * e + 1)
    return np.where(mask, z, 0.0)


def smooth_oe(r, w, mask, k=10):
    r = np.where(mask, r, 0.0)
    e = r.sum() * w
    ks = np.ones(k)
    rs, es = np.convolve(r, ks, "same"), np.convolve(e, ks, "same")
    return np.where(es > 0, np.log2((rs + .5) / (es + .5)), 0)


def effect_size(r, w, mask, n_null=8, seed=0):
    """Amplitude of the signal beyond Poisson noise: sqrt(var_obs - var_null) of the 50-kb smoothed log2 O/E.
    Known hitchhiking plasmids: 0.35-0.71; organelles / rDNA episome: 0.17-0.23."""
    rng = np.random.default_rng(seed)
    m = mask > 0
    v = smooth_oe(r, w, mask)[m].var()
    vn = np.mean([smooth_oe(rng.multinomial(int(r.sum()), w).astype(float), w, mask)[m].var() for _ in range(n_null)])
    return float(np.sqrt(max(v - vn, 0)))


def autocorr1(r, w, mask):
    """Lag-1 autocorrelation of the residual: true hotspots span several bins, artefactual spikes do not."""
    z = residual(r, w, mask)[mask > 0]
    if z.std() == 0:
        return float("nan")
    return float(np.corrcoef(z[:-1], z[1:])[0, 1])
