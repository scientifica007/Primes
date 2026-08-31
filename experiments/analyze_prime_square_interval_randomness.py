#!/usr/bin/env python3
"""Randomness experiment for prime counts between consecutive prime squares.

Tests a specific baseline, not an absolute notion of randomness:
    X_n ~ independent Poisson(alpha * mu_n)
    mu_n = integral_{p_n^2}^{p_(n+1)^2} dt/log(t)

The single factor alpha calibrates the total expected count to the observed total,
so the main tests focus on local variation and dependence rather than a tiny
finite-range density bias in the prime number theorem approximation.
"""

import math
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.special import expi
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parent
INPUT = ROOT / "first_1000_prime_square_interval_index.csv"
REPS = 50000
SEED = 20260831


def acf(v, lag):
    return np.corrcoef(v[:-lag], v[lag:])[0, 1]


def p_upper(obs, sims):
    return (1 + np.sum(sims >= obs)) / (len(sims) + 1)


def p_lower(obs, sims):
    return (1 + np.sum(sims <= obs)) / (len(sims) + 1)


def p_two_abs(obs, sims):
    return (1 + np.sum(np.abs(sims) >= abs(obs))) / (len(sims) + 1)


df = pd.read_csv(INPUT)
n = df["الترتيب_n"].to_numpy(dtype=int)
p = df["p_n"].to_numpy(dtype=int)
p_next = df["p_التالي"].to_numpy(dtype=int)
a = df["p_n^2"].to_numpy(dtype=np.float64)
b = df["p_التالي^2"].to_numpy(dtype=np.float64)
x = df["عدد_الأوليات_بين_p2_ومربع_التالي"].to_numpy(dtype=np.float64)

gap = p_next - p
interval_length = b - a

# li(b)-li(a) = integral_a^b dt/log(t), for a,b > 1.
mu = expi(np.log(b)) - expi(np.log(a))
alpha = x.sum() / mu.sum()
mu_cal = alpha * mu
z = (x - mu_cal) / np.sqrt(mu_cal)

acf_obs = np.array([acf(z, k) for k in range(1, 21)])
D_obs = np.sum((x - mu_cal) ** 2 / mu_cal)
dispersion_ratio = D_obs / (len(x) - 1)
lag1_obs = acf_obs[0]
maxacf_obs = np.max(np.abs(acf_obs))

signs = z >= 0
runs_obs = 1 + np.sum(signs[1:] != signs[:-1])
npos = int(signs.sum())
nneg = len(signs) - npos
runs_mean_classical = 1 + 2 * npos * nneg / len(signs)
runs_var_classical = (
    2 * npos * nneg * (2 * npos * nneg - npos - nneg)
    / ((len(signs) ** 2) * (len(signs) - 1))
)
runs_z_classical = (runs_obs - runs_mean_classical) / math.sqrt(runs_var_classical)

corr_raw_length = np.corrcoef(x, interval_length)[0, 1]
corr_raw_gap = np.corrcoef(x, gap)[0, 1]
corr_z_n = np.corrcoef(z, n)[0, 1]
corr_z_gap = np.corrcoef(z, gap)[0, 1]
corr_z_length = np.corrcoef(z, interval_length)[0, 1]
global_bias_sigma = (x.sum() - mu.sum()) / math.sqrt(mu.sum())

# Monte Carlo under the calibrated independent-Poisson baseline.
rng = np.random.default_rng(SEED)
D_sim = np.empty(REPS)
lag1_sim = np.empty(REPS)
maxacf_sim = np.empty(REPS)
runs_sim = np.empty(REPS, dtype=np.int32)
corr_n_sim = np.empty(REPS)
corr_gap_sim = np.empty(REPS)

n_center = n - n.mean()
gap_center = gap - gap.mean()
n_norm = np.sqrt(np.sum(n_center ** 2))
gap_norm = np.sqrt(np.sum(gap_center ** 2))

batch = 500
for start in range(0, REPS, batch):
    m = min(batch, REPS - start)
    Y = rng.poisson(mu_cal, size=(m, len(mu_cal))).astype(float)
    Z = (Y - mu_cal) / np.sqrt(mu_cal)

    D_sim[start:start+m] = np.sum((Y - mu_cal) ** 2 / mu_cal, axis=1)

    A, B = Z[:, :-1], Z[:, 1:]
    Ac = A - A.mean(axis=1, keepdims=True)
    Bc = B - B.mean(axis=1, keepdims=True)
    lag1_sim[start:start+m] = np.sum(Ac * Bc, axis=1) / np.sqrt(
        np.sum(Ac ** 2, axis=1) * np.sum(Bc ** 2, axis=1)
    )

    mx = np.zeros(m)
    for k in range(1, 21):
        A, B = Z[:, :-k], Z[:, k:]
        Ac = A - A.mean(axis=1, keepdims=True)
        Bc = B - B.mean(axis=1, keepdims=True)
        r = np.sum(Ac * Bc, axis=1) / np.sqrt(
            np.sum(Ac ** 2, axis=1) * np.sum(Bc ** 2, axis=1)
        )
        mx = np.maximum(mx, np.abs(r))
    maxacf_sim[start:start+m] = mx

    S = Z >= 0
    runs_sim[start:start+m] = 1 + np.sum(S[:, 1:] != S[:, :-1], axis=1)

    Zc = Z - Z.mean(axis=1, keepdims=True)
    denZ = np.sqrt(np.sum(Zc ** 2, axis=1))
    corr_n_sim[start:start+m] = (Zc @ n_center) / (denZ * n_norm)
    corr_gap_sim[start:start+m] = (Zc @ gap_center) / (denZ * gap_norm)

pD_lo = p_lower(D_obs, D_sim)
pD_hi = p_upper(D_obs, D_sim)
pD = min(1.0, 2 * min(pD_lo, pD_hi))
p_lag1 = p_two_abs(lag1_obs, lag1_sim)
p_maxacf = p_upper(maxacf_obs, maxacf_sim)
p_runs = (
    1
    + np.sum(np.abs(runs_sim - runs_sim.mean()) >= abs(runs_obs - runs_sim.mean()))
) / (REPS + 1)
p_corr_n = p_two_abs(corr_z_n, corr_n_sim)
p_corr_gap = p_two_abs(corr_z_gap, corr_gap_sim)

rows = pd.DataFrame({
    "n": n,
    "p_n": p,
    "p_next": p_next,
    "prime_gap": gap,
    "p_n_squared": a.astype(np.int64),
    "p_next_squared": b.astype(np.int64),
    "interval_length": interval_length.astype(np.int64),
    "observed_prime_count": x.astype(int),
    "cramer_expected_count": mu,
    "global_density_scale_alpha": alpha,
    "calibrated_expected_count": mu_cal,
    "residual_observed_minus_expected": x - mu_cal,
    "standardized_residual_z": z,
})
rows.to_csv(ROOT / "prime_square_interval_randomness_results.csv", index=False, encoding="utf-8-sig")

top = rows.iloc[np.argsort(np.abs(z))[::-1][:15]].copy()
top.to_csv(ROOT / "prime_square_interval_randomness_top_deviations.csv", index=False, encoding="utf-8-sig")

summary = pd.DataFrame([
    ["sample_size_intervals", len(x), "", ""],
    ["observed_total_primes", int(x.sum()), "", ""],
    ["cramer_expected_total", float(mu.sum()), "", ""],
    ["global_scale_alpha_observed_over_expected", float(alpha), "", ""],
    ["uncalibrated_total_bias_in_poisson_sigma", float(global_bias_sigma), "", ""],
    ["raw_corr_count_vs_interval_length", float(corr_raw_length), "", ""],
    ["raw_corr_count_vs_prime_gap", float(corr_raw_gap), "", ""],
    ["pearson_dispersion_D", float(D_obs), float(pD), "two-sided Monte Carlo; calibrated Poisson"],
    ["dispersion_ratio_D_over_df", float(dispersion_ratio), "", "df=999"],
    ["lag1_autocorrelation_z", float(lag1_obs), float(p_lag1), "two-sided Monte Carlo"],
    ["max_abs_autocorrelation_lags_1_20", float(maxacf_obs), float(p_maxacf), "upper-tail Monte Carlo with multiple-lag correction"],
    ["runs_count_residual_signs", int(runs_obs), float(p_runs), "Monte Carlo"],
    ["runs_test_z_classical", float(runs_z_classical), "", ""],
    ["corr_z_vs_index_n", float(corr_z_n), float(p_corr_n), "two-sided Monte Carlo"],
    ["corr_z_vs_prime_gap", float(corr_z_gap), float(p_corr_gap), "two-sided Monte Carlo"],
    ["corr_z_vs_interval_length", float(corr_z_length), "", ""],
], columns=["metric", "observed_value", "monte_carlo_p_value", "notes"])
summary.to_csv(ROOT / "prime_square_interval_randomness_summary.csv", index=False, encoding="utf-8-sig")

# Reproducible SVG plots.
fig, ax = plt.subplots(figsize=(10, 6))
ax.scatter(mu_cal, x, s=10, alpha=0.7)
lo, hi = min(mu_cal.min(), x.min()), max(mu_cal.max(), x.max())
ax.plot([lo, hi], [lo, hi], linewidth=1)
ax.set_xlabel("Calibrated expected prime count")
ax.set_ylabel("Observed prime count")
ax.set_title("Prime counts between consecutive prime squares: observed vs expected")
ax.grid(True, alpha=0.25)
fig.tight_layout()
fig.savefig(ROOT / "prime_square_counts_observed_vs_expected.svg")
plt.close(fig)

fig, ax = plt.subplots(figsize=(10, 6))
ax.plot(n, z, linewidth=0.8)
ax.axhline(0, linewidth=1)
ax.set_xlabel("n")
ax.set_ylabel("Standardized residual z")
ax.set_title("Standardized residuals under calibrated Cramer-Poisson baseline")
ax.grid(True, alpha=0.25)
fig.tight_layout()
fig.savefig(ROOT / "prime_square_interval_standardized_residuals.svg")
plt.close(fig)

# Human-readable report.
min_mc = 1 / (REPS + 1)
report = f"""# تجربة العشوائية: عدد الأوليات بين مربعين أوليين متتاليين

## السؤال

ندرس السلسلة:

$$
X_n=\#\{{q\text{{ أولي}}:p_n^2<q<p_{{n+1}}^2\}},\qquad n=1,\ldots,1000.
$$

السؤال ليس «هل الأوليات عشوائية مطلقًا؟»، بل:

> هل تتصرف أعداد الأوليات في هذه الفترات كما لو كانت ناتجة عن نموذج عشوائي مستقل بسيط بعد أخذ تغير طول الفترات وكثافة الأوليات في الحسبان؟

## النموذج المرجعي

استعملنا نموذج كرامر/بواسون مبسطًا. لكل فترة حسبنا:

$$
\mu_n=\int_{{p_n^2}}^{{p_{{n+1}}^2}}\frac{{dt}}{{\log t}}.
$$

ثم استعملنا معامل معايرة عالمي واحد:

$$
\alpha=\frac{{\sum X_n}}{{\sum\mu_n}}={alpha:.12f}.
$$

وأصبح التوقع المحلي `alpha * mu_n`. هذه المعايرة تمنع فرقًا إجماليًا صغيرًا في تقريب كثافة الأوليات من السيطرة على اختبارات البنية المحلية.

بعدها عرّفنا الباقي المعياري:

$$
Z_n=\frac{{X_n-\alpha\mu_n}}{{\sqrt{{\alpha\mu_n}}}}.
$$

وأجرينا {REPS:,} محاكاة Monte Carlo مستقلة من توزيع بواسون بهذه المتوسطات، ببذرة ثابتة `{SEED}`.

## تحقق أولي من الكثافة

- مجموع الأوليات المرصود في الفترات: **{int(x.sum()):,}**.
- مجموع توقع كرامر غير المعاير: **{mu.sum():,.3f}**.
- الفرق الكلي يساوي فقط **{global_bias_sigma:.3f}** انحراف معياري بواسوني تقريبي.
- الارتباط بين العدد الخام للأوليات وطول الفترة: **{corr_raw_length:.6f}**.
- الارتباط بين العدد الخام وفجوة الأوليين `p_(n+1)-p_n`: **{corr_raw_gap:.6f}**.

إذن معظم تغير العدد الخام مفهوم مباشرة من اختلاف أطوال الفترات؛ لذلك لا يجوز اختبار السلسلة الخام كما لو كانت ثابتة المتوسط.

## النتائج الأساسية بعد التطبيع

| الاختبار | القيمة المرصودة | قيمة p من Monte Carlo | القراءة الأولية |
|---|---:|---:|---|
| تشتت Pearson | {D_obs:.6f} | {pD:.8f} | أقل بكثير من بواسون |
| نسبة التشتت `D/999` | {dispersion_ratio:.6f} | — | نحو 24% من التشتت البواسوني |
| الارتباط الذاتي عند lag=1 | {lag1_obs:.6f} | {p_lag1:.8f} | ارتباط سلبي واضح |
| أكبر `|ACF|` بين lag 1..20 | {maxacf_obs:.6f} | {p_maxacf:.8f} | غير معتاد تحت النموذج، مع تصحيح البحث عبر 20 lag |
| عدد التتابعات لإشارة الباقي | {runs_obs} | {p_runs:.8f} | تبدلات أكثر من المتوقع نسبيًا |
| ارتباط الباقي بالترتيب n | {corr_z_n:.6f} | {p_corr_n:.8f} | لا دليل واضح على اتجاه زمني بسيط |
| ارتباط الباقي بفجوة الأوليين | {corr_z_gap:.6f} | {p_corr_gap:.8f} | صغير، لكنه غير معتاد تحت النموذج |

مع {REPS:,} محاكاة، أصغر قيمة p تجريبية يمكن تمييزها هي تقريبًا `{min_mc:.8f}`؛ لذلك القيم عند هذا الحد تعني «أصغر من قدرة هذه المحاكاة على القياس بدقة أكبر»، لا قيمة تحليلية دقيقة.

## أهم ما تعلمناه

### 1. النموذج العشوائي البسيط مرفوض بقوة

أوضح إشارة هي التشتت. لو كانت الأعداد داخل الفترات تتصرف كمتغيرات بواسون مستقلة، نتوقع تقريبًا أن تكون قيمة Pearson قرب عدد درجات الحرية، أي قرب 999. المرصود هو:

$$
D={D_obs:.3f},\qquad D/999={dispersion_ratio:.3f}.
$$

هذا يعني أن العد الحقيقي **أكثر انتظامًا بكثير** من نموذج «رميات مستقلة» بسيط.

هذا لا يعني أن لدينا «نسبة عشوائية 24%». الرقم 24% هنا هو نسبة تشتت محددة بالنسبة إلى نموذج بواسون محدد فقط.

### 2. يوجد أثر اعتماد بين الفترات المتجاورة

وجدنا:

$$
r_1={lag1_obs:.3f}.
$$

أي أن الفترة التي تأتي فوق توقع النموذج قليلًا تميل، إحصائيًا، إلى أن يتبعها انحراف في الاتجاه المعاكس أكثر مما يسمح به نموذج الاستقلال البسيط. ويؤيد ذلك اختبار التتابع الذي وجد تبدلات أكثر من المتوقع.

### 3. هذا ليس اكتشافًا أن «الأوليات غير عشوائية»

الأعداد الأولية أصلًا حتمية. ما أثبتناه تجريبيًا أضيق بكثير:

> نموذج بواسون المستقل بكثافة `1/log(x)` أبسط من أن يصف التقلبات المحلية لهذه البيانات.

وهذا متوقع نوعيًا لأن الأوليات ليست نقاطًا مستقلة: القيود الناتجة عن القسمة على 2 و3 و5 و7... تمنع كثيرًا من المواقع بصورة دورية. أي أن نموذج العجلات الذي ندرسه نفسه يعطينا سببًا فوريًا للشك في استقلال بواسون.

## السؤال الذي ولدته التجربة

الخطوة التالية الطبيعية ليست إعلان «وجود نمط جديد»، بل بناء نموذج عشوائي أكثر عدلًا يحترم جزءًا من بنية الغربال.

مثلًا: نزيل مسبقًا المواقع التي تضربها عجلات 2 و3 و5 و7... ثم نختار نقاطًا عشوائية فقط من المواضع الباقية مع معايرة الكثافة، وبعد ذلك نعيد الاختبارات نفسها.

السؤال الجديد هو:

> كم من الانخفاض الشديد في التشتت والارتباط السلبي الذي رأيناه تفسره قيود الغربال/العجلات المعروفة، وكم يبقى بعد إدخالها في النموذج العشوائي؟

هذا هو الانتقال الصحيح من «رفض نموذج بسيط» إلى محاولة تحديد مصدر البنية.

## الملفات

- `prime_square_interval_randomness_results.csv`: النتائج التفصيلية لكل فترة.
- `prime_square_interval_randomness_summary.csv`: ملخص الاختبارات.
- `prime_square_interval_randomness_top_deviations.csv`: أكبر 15 انحرافًا معياريًا.
- `prime_square_counts_observed_vs_expected.svg`: المرصود مقابل المتوقع.
- `prime_square_interval_standardized_residuals.svg`: البواقي المعيارية بالتتابع.
- `analyze_prime_square_interval_randomness.py`: سكربت إعادة التجربة.
"""
(ROOT / "PRIME_SQUARE_INTERVAL_RANDOMNESS_EXPERIMENT.md").write_text(report, encoding="utf-8")

print(summary.to_string(index=False))
