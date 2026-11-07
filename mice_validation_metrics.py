"""
MICE-based temporal-validation metrics for the frozen M1 model, C1-unflagged headline stratum.

Inputs (from the user's pipeline, unchanged):
  data/processed/val_headline_mi_preds.rds   -- output of R/11_mice_imputation.R: headline-stratum rows
                                                 with the predict-then-pool MI prediction m1_mi_pred (n=3,064)
  data/processed/analytic_pooled.rds         -- full pooled analytic file (for the validation design: strata/PSU)

Everything here is a design-based, survey-weighted quantity (WTPH2YR primary weight, SDMVSTRA/SDMVPSU
design), with 95% confidence intervals from a Rao-Wu rescaled bootstrap (R = 500 replicates, the
protocol's pre-specified number), resampling n_h - 1 PSUs with replacement within each stratum.
With NHANES' 2 PSUs per stratum the replicate weight is 2*w for the selected PSU and 0 for the other
(Rao, Wu & Yue 1992), which is also what survey::as.svrepdesign(type="mrbbootstrap") produces for a
single-stage PSU design.

Mirrors the logic of R/08_decision_curve_calibration.R and R/13_subgroup_equity_mi.R but on the MI
predictions, and extends it with: calibration slope/intercept, weighted C-statistic, CIs on O/E and on
every subgroup O/E and NB, the full decision curve with CIs for all strategies, the calibration slope
under the pre-specified >=50-event rule + collapsing hierarchy, and sensitivity analyses S1/S3/S9.
"""
import json, sys
import numpy as np, pandas as pd
sys.path.insert(0, '.')
from rds_reader import read_rds

UP = '../../data/processed/'   # run from results/mice_metrics_python/
OUT = './'
R_BOOT = 500
SEED = 20260921
THRESHOLDS = np.round(np.arange(0.02, 0.2001, 0.01), 2)
PT_PRIMARY = 0.05
MIN_EVENTS_SLOPE = 50
RACE_LABELS = {1: "Mexican American", 2: "Other Hispanic", 3: "Non-Hispanic White",
               4: "Non-Hispanic Black", 6: "Non-Hispanic Asian", 7: "Other/Multiracial"}
INCOME_BREAKS = [-np.inf, 1.0, 2.0, 4.0, np.inf]
INCOME_LABELS = ["<1.0 (below FPL)", "1.0-1.99", "2.0-3.99", ">=4.0"]

# ---------------------------------------------------------------- data
mi = read_rds(UP + 'val_headline_mi_preds.rds')
pooled = read_rds(UP + 'analytic_pooled.rds')
val = pooled[(pooled.role == 'validation') & pooled.WTPH2YR.notna()].copy()
val['combined_stratum'] = val.cycle + '_' + val.sdmvstra.astype(int).astype(str)
val['combined_psu'] = val.combined_stratum + '_' + val.sdmvpsu.astype(int).astype(str)
assert mi.SEQN.isin(val.SEQN).all()
d = mi.copy()
assert (d.combined_stratum == val.set_index('SEQN').loc[d.SEQN, 'combined_stratum'].values).all()

y = d.ckd_case.astype(bool).to_numpy()
p = d.m1_mi_pred.to_numpy().astype(float)
w = d.WTPH2YR.to_numpy().astype(float)
w_mec = d.wt_mec.to_numpy().astype(float)
c2 = d.c2_flag.astype(bool).to_numpy()
y_s1 = ((d.egfr_race_free.to_numpy() < 60) | (d.acr_uncorrected.to_numpy() >= 30))
n = len(d)
print(f"headline MI rows n={n}, events={y.sum()}")

# ---------------------------------------------------------------- bootstrap replicate weights
rng = np.random.default_rng(SEED)
strata = val.combined_stratum.unique()
psu_of_stratum = {s: sorted(val.loc[val.combined_stratum == s, 'combined_psu'].unique()) for s in strata}
assert all(len(v) == 2 for v in psu_of_stratum.values()), "expected 2 PSUs per stratum"
psu_idx = {ps: i for i, ps in enumerate(sorted(val.combined_psu.unique()))}
row_psu = d.combined_psu.map(psu_idx).to_numpy()
mult = np.zeros((R_BOOT, len(psu_idx)))          # replicate multiplier per PSU
for r in range(R_BOOT):
    for s, psus in psu_of_stratum.items():
        n_h = len(psus); m_h = n_h - 1
        draws = rng.choice(n_h, size=m_h, replace=True)
        counts = np.bincount(draws, minlength=n_h)
        lam = np.sqrt(m_h / (n_h - 1))               # Rao-Wu-Yue rescaling (fpc ~ 0)
        for i, ps in enumerate(psus):
            mult[r, psu_idx[ps]] = 1 - lam + lam * (n_h / m_h) * counts[i]
W = w[None, :] * mult[:, row_psu]                  # R x n replicate weights (WTPH2YR)
W_mec = w_mec[None, :] * mult[:, row_psu]

# ---------------------------------------------------------------- statistics
def expit(x): return 1 / (1 + np.exp(-x))
lp = np.log(np.clip(p, 1e-6, 1 - 1e-6) / (1 - np.clip(p, 1e-6, 1 - 1e-6)))

def oe(y, p, w):
    return (w * y).sum() / (w * p).sum()

def cal_intercept(y, lp, w, iters=50):
    a = 0.0
    for _ in range(iters):
        mu = expit(a + lp)
        g = (w * (y - mu)).sum(); h = (w * mu * (1 - mu)).sum()
        step = g / h; a += step
        if abs(step) < 1e-10: break
    return a

def cal_slope(y, lp, w, iters=50):
    X = np.column_stack([np.ones_like(lp), lp]); b = np.array([0.0, 1.0])
    for _ in range(iters):
        mu = expit(X @ b); g = X.T @ (w * (y - mu)); H = (X * (w * mu * (1 - mu))[:, None]).T @ X
        step = np.linalg.solve(H, g); b += step
        if np.abs(step).max() < 1e-10: break
    return b  # (intercept, slope)

def weighted_c(y, p, w):
    """Design-weighted concordance: P(p_case > p_control) with weights w_i*w_j, ties count 1/2."""
    order = np.argsort(p, kind='mergesort'); ps, ys, ws = p[order], y[order], w[order]
    # group ties
    uniq, inv = np.unique(ps, return_inverse=True)
    wc = np.bincount(inv, weights=ws * ys, minlength=len(uniq))      # case weight per unique value
    wn = np.bincount(inv, weights=ws * (~ys), minlength=len(uniq))   # control weight per unique value
    cum_wn = np.cumsum(wn) - wn                                      # controls strictly below
    conc = (wc * (cum_wn + 0.5 * wn)).sum()
    return conc / (wc.sum() * wn.sum())

def nb(treat, y, w, pt):
    return ((w * (treat & y)).sum() - (w * (treat & ~y)).sum() * pt / (1 - pt)) / w.sum()

def all_stats(y, p, lp, w, c2, thresholds=THRESHOLDS):
    out = {'obs': (w * y).sum() / w.sum(), 'exp': (w * p).sum() / w.sum()}
    out['oe'] = out['obs'] / out['exp']
    out['cal_int'] = cal_intercept(y, lp, w)
    out['cal_slope'] = cal_slope(y, lp, w)[1]
    out['c'] = weighted_c(y, p, w)
    for pt in thresholds:
        k = f"{pt:.2f}"
        out[f'nb_model_{k}'] = nb(p >= pt, y, w, pt)
        out[f'nb_all_{k}'] = nb(np.ones(len(y), bool), y, w, pt)
        out[f'nb_c2_{k}'] = nb(c2, y, w, pt)
    return out

def boot(fn, W):
    reps = [fn(W[r]) for r in range(W.shape[0])]
    keys = reps[0].keys()
    return {k: np.array([r[k] for r in reps]) for k in keys}

def ci(arr, est=None):
    lo, hi = np.nanpercentile(arr, [2.5, 97.5]); se = np.nanstd(arr, ddof=1)
    return {'lo': float(lo), 'hi': float(hi), 'se': float(se)}

results = {}

# ---- primary: headline stratum, WTPH2YR, bridging-corrected outcome, MI predictions
est = all_stats(y, p, lp, w, c2)
rep = boot(lambda ww: all_stats(y, p, lp, ww, c2), W)
results['primary'] = {k: {'est': float(v), **ci(rep[k])} for k, v in est.items()}
print("\n=== PRIMARY (MICE, WTPH2YR, headline stratum) ===")
for k in ['obs', 'exp', 'oe', 'cal_int', 'cal_slope', 'c', 'nb_model_0.05', 'nb_all_0.05', 'nb_c2_0.05', 'nb_model_0.10']:
    r = results['primary'][k]; print(f"{k:14s} {r['est']:.4f}  95% CI {r['lo']:.4f} to {r['hi']:.4f}  (SE {r['se']:.4f})")

# crossing points of comparators (threshold where NB turns negative)
def first_negative(prefix):
    for pt in THRESHOLDS:
        if results['primary'][f'{prefix}_{pt:.2f}']['est'] < 0: return float(pt)
    return None
results['crossing'] = {'test_all_negative_from': first_negative('nb_all'), 'c2_negative_from': first_negative('nb_c2'),
                       'model_min_nb': min(results['primary'][f'nb_model_{pt:.2f}']['est'] for pt in THRESHOLDS)}
# model minus comparator difference at 5% with CI (paired within replicate)
for comp in ['all', 'c2']:
    diff = rep['nb_model_0.05'] - rep[f'nb_{comp}_0.05']
    results['primary'][f'nb_diff_model_minus_{comp}_0.05'] = {'est': float(est['nb_model_0.05'] - est[f'nb_{comp}_0.05']), **ci(diff)}
    print(f"model - {comp} NB@5%: {results['primary'][f'nb_diff_model_minus_{comp}_0.05']}")
print("crossing:", results['crossing'])

# ---- subgroups: O/E, NB@5%, slope where >=50 events (full granularity + collapsing hierarchy)
race = d.race_eth.to_numpy().astype(int)
income_cat = pd.cut(d.income_pir, bins=INCOME_BREAKS, labels=INCOME_LABELS)
ins = d.hiq011.to_numpy()
groups = {}
for code, lab in RACE_LABELS.items(): groups[('race', lab)] = race == code
groups[('race_collapsed', 'Hispanic (combined)')] = (race == 1) | (race == 2)
for lab in INCOME_LABELS: groups[('income', lab)] = (income_cat == lab).to_numpy()
groups[('income_collapsed', '<2.0x FPL')] = (income_cat.isin(INCOME_LABELS[:2])).to_numpy()
groups[('income_collapsed', '>=2.0x FPL')] = (income_cat.isin(INCOME_LABELS[2:])).to_numpy()
groups[('insurance', 'Insured')] = ins == 1
groups[('insurance', 'Uninsured')] = ins == 2
results['subgroups'] = []
for (dim, lab), m in groups.items():
    ev = int(y[m].sum()); nn = int(m.sum())
    e_oe = oe(y[m], p[m], w[m]); e_nb = nb(p[m] >= PT_PRIMARY, y[m], w[m], PT_PRIMARY)
    r_oe = np.array([oe(y[m], p[m], W[r][m]) for r in range(R_BOOT)])
    r_nb = np.array([nb(p[m] >= PT_PRIMARY, y[m], W[r][m], PT_PRIMARY) for r in range(R_BOOT)])
    row = {'dimension': dim, 'level': lab, 'n': nn, 'events': ev,
           'oe': {'est': float(e_oe), **ci(r_oe)}, 'nb_5pct': {'est': float(e_nb), **ci(r_nb)}}
    if ev >= MIN_EVENTS_SLOPE:
        e_sl = cal_slope(y[m], lp[m], w[m])[1]; e_in = cal_intercept(y[m], lp[m], w[m])
        r_sl = np.array([cal_slope(y[m], lp[m], W[r][m])[1] for r in range(R_BOOT)])
        r_in = np.array([cal_intercept(y[m], lp[m], W[r][m]) for r in range(R_BOOT)])
        row['cal_slope'] = {'est': float(e_sl), **ci(r_sl)}; row['cal_int'] = {'est': float(e_in), **ci(r_in)}
    results['subgroups'].append(row)
    s = f"{dim:17s} {lab:22s} n={nn:5d} ev={ev:3d}  O/E {e_oe:.2f} ({row['oe']['lo']:.2f}-{row['oe']['hi']:.2f})  NB5 {e_nb*1000:.1f} ({row['nb_5pct']['lo']*1000:.1f}-{row['nb_5pct']['hi']*1000:.1f})"
    if 'cal_slope' in row: s += f"  slope {row['cal_slope']['est']:.2f} ({row['cal_slope']['lo']:.2f}-{row['cal_slope']['hi']:.2f})"
    print(s)

# weighted SD of log(O/E) across race levels (weights = unweighted n, as in R/13) + max abs difference, with percentile intervals
race_masks = [race == c for c in RACE_LABELS]
def sd_log_oe(ww):
    logs = np.array([np.log(oe(y[m], p[m], ww[m])) for m in race_masks]); ns = np.array([m.sum() for m in race_masks])
    mu = (ns * logs).sum() / ns.sum(); return np.sqrt((ns * (logs - mu) ** 2).sum() / ns.sum())
def max_abs_diff(ww):
    oes = np.array([oe(y[m], p[m], ww[m]) for m in race_masks]); return oes.max() - oes.min()
results['race_summary'] = {'sd_log_oe': {'est': float(sd_log_oe(w)), **ci(np.array([sd_log_oe(W[r]) for r in range(R_BOOT)]))},
                           'max_abs_diff_oe': {'est': float(max_abs_diff(w)), **ci(np.array([max_abs_diff(W[r]) for r in range(R_BOOT)]))}}
print("race summary:", results['race_summary'])
# NH Black vs NH White O/E ratio with CI (paired)
mb, mw_ = race == 4, race == 3
ratio = np.array([oe(y[mb], p[mb], W[r][mb]) / oe(y[mw_], p[mw_], W[r][mw_]) for r in range(R_BOOT)])
results['race_summary']['oe_ratio_black_vs_white'] = {'est': float(oe(y[mb], p[mb], w[mb]) / oe(y[mw_], p[mw_], w[mw_])), **ci(ratio)}
print("NHB/NHW O/E ratio:", results['race_summary']['oe_ratio_black_vs_white'])

# ---- sensitivity analyses under MICE: S1 (no bridging correction), S3 (WTMEC2YR), S9 (C2-unflagged stratum)
def sens(label, yy, pp, ww, WW, cc2, mask=None):
    if mask is not None:
        yy, pp, ww, WW, cc2 = yy[mask], pp[mask], ww[mask], WW[:, mask], cc2[mask]
    lpp = np.log(np.clip(pp, 1e-6, 1 - 1e-6) / (1 - np.clip(pp, 1e-6, 1 - 1e-6)))
    e = all_stats(yy, pp, lpp, ww, cc2, thresholds=[PT_PRIMARY])
    rr = boot(lambda x: all_stats(yy, pp, lpp, x, cc2, thresholds=[PT_PRIMARY]), WW)
    res = {k: {'est': float(v), **ci(rr[k])} for k, v in e.items()}
    res['n'] = int(len(yy)); res['events'] = int(yy.sum())
    print(f"{label}: n={len(yy)} ev={yy.sum()} O/E {res['oe']['est']:.3f} ({res['oe']['lo']:.2f}-{res['oe']['hi']:.2f}) "
          f"slope {res['cal_slope']['est']:.2f} C {res['c']['est']:.3f} NB5 {res['nb_model_0.05']['est']*1000:.1f} ({res['nb_model_0.05']['lo']*1000:.1f}-{res['nb_model_0.05']['hi']*1000:.1f}) test-all {res['nb_all_0.05']['est']*1000:.1f} C2 {res['nb_c2_0.05']['est']*1000:.1f}")
    return res
print("\n=== SENSITIVITY (MICE) ===")
results['sensitivity'] = {
    'S1_no_bridging': sens('S1 no bridging', y_s1, p, w, W, c2),
    'S3_wtmec2yr': sens('S3 WTMEC2YR', y, p, w_mec, W_mec, c2),
    'S9_c2_unflagged': sens('S9 C2-unflagged', y, p, w, W, c2, mask=~c2),
}
# bridging flip counts within the MI analytic rows
results['bridging_flips'] = {'newly_case': int((y & ~y_s1).sum()), 'newly_noncase': int((~y & y_s1).sum()), 'n': int(n)}
print("bridging flips (MI rows):", results['bridging_flips'])

# ---- decision-curve table
dca = pd.DataFrame({'threshold': THRESHOLDS})
for strat in ['model', 'all', 'c2']:
    dca[f'nb_{strat}'] = [results['primary'][f'nb_{strat}_{pt:.2f}']['est'] for pt in THRESHOLDS]
    dca[f'nb_{strat}_lo'] = [results['primary'][f'nb_{strat}_{pt:.2f}']['lo'] for pt in THRESHOLDS]
    dca[f'nb_{strat}_hi'] = [results['primary'][f'nb_{strat}_{pt:.2f}']['hi'] for pt in THRESHOLDS]
dca.to_csv(OUT + 'decision_curve_mice.csv', index=False)
print("\n", dca[['threshold', 'nb_model', 'nb_model_lo', 'nb_model_hi', 'nb_all', 'nb_c2']].round(4).to_string(index=False))

# weighted proportion flagged by the model at 5% (tests per 1,000) -- useful for "tests performed per 1,000"
for pt in [0.05, 0.10]:
    flag = p >= pt
    results['primary'][f'tested_per_1000_{pt:.2f}'] = {'est': float((w * flag).sum() / w.sum() * 1000)}
    results['primary'][f'cases_found_per_1000_{pt:.2f}'] = {'est': float((w * (flag & y)).sum() / w.sum() * 1000)}
results['primary']['events_per_1000_total'] = {'est': float((w * y).sum() / w.sum() * 1000)}
results['primary']['c2_tested_per_1000'] = {'est': float((w * c2).sum() / w.sum() * 1000)}
results['primary']['c2_cases_found_per_1000'] = {'est': float((w * (c2 & y)).sum() / w.sum() * 1000)}
print({k: v for k, v in results['primary'].items() if 'per_1000' in k})

results['meta'] = {'n': int(n), 'events': int(y.sum()), 'R_boot': R_BOOT, 'seed': SEED, 'strata': len(strata), 'psus': len(psu_idx)}
json.dump(results, open(OUT + 'mice_validation_results.json', 'w'), indent=1)
pd.DataFrame([{'dimension': r['dimension'], 'level': r['level'], 'n': r['n'], 'events': r['events'],
               'oe': r['oe']['est'], 'oe_lo': r['oe']['lo'], 'oe_hi': r['oe']['hi'],
               'nb_5pct': r['nb_5pct']['est'], 'nb_lo': r['nb_5pct']['lo'], 'nb_hi': r['nb_5pct']['hi'],
               'cal_slope': r.get('cal_slope', {}).get('est'), 'slope_lo': r.get('cal_slope', {}).get('lo'), 'slope_hi': r.get('cal_slope', {}).get('hi')}
              for r in results['subgroups']]).to_csv(OUT + 'subgroup_oe_mice_ci.csv', index=False)
print("saved.")
