Flexible dose-response models
================

- [Clinical outcome, safety measure, time-varying confounders, and
  dataset](#clinical-outcome-safety-measure-time-varying-confounders-and-dataset)
  - [Empirical support and positivity
    diagnostics](#empirical-support-and-positivity-diagnostics)
- [Notation](#notation)
- [Step 1 - Inverse probability of dose and censoring weighting
  (IPDCW)](#step-1---inverse-probability-of-dose-and-censoring-weighting-ipdcw)
  - [Inverse probability of dose weights
    (IPTW)](#inverse-probability-of-dose-weights-iptw)
  - [Inverse probability of censoring weights
    (IPCW)](#inverse-probability-of-censoring-weights-ipcw)
  - [Total stabilized weights and
    truncation](#total-stabilized-weights-and-truncation)
- [Step 2 - Weighted repeated-measures marginal structural model
  MSM](#step-2---weighted-repeated-measures-marginal-structural-model-msm)
- [Step 3 - Predictions](#step-3---predictions)
- [References](#references)
- [Appendix](#appendix)
  - [Predicted probabilities under the ordinal dose
    model](#predicted-probabilities-under-the-ordinal-dose-model)
  - [Alternative multinomial dose model for stabilized dose
    weights](#alternative-multinomial-dose-model-for-stabilized-dose-weights)
  - [Predicted trajectories using multinomial
    IPTW](#predicted-trajectories-using-multinomial-iptw)

## Clinical outcome, safety measure, time-varying confounders, and dataset

The clinical outcome of the included study is the Hamilton Depression
Rating Scale (HAMD) score, denoted by $Y_{it}$, measured for patient $i$
at visit $t$. The baseline HAMD score was denoted by $Y_{i0}$. Treatment
response was expressed as improvement from baseline:

``` math
\Delta Y_{it} = Y_{i0} - Y_{it}
```

Thus, positive values of $ΔY_{it}$ indicate improvement. Safety was
summarized at each visit using a side-effect score $S_{it}$, measured on
a 0–10 scale, where higher values indicate a greater number of side
effects reported by the patient.

In this setting, we distinguished between previous dose history and
time-varying confounder history. Previous dose $d_{it-1}$ was treated as
part of the treatment process, because it defines the dose trajectory
and is a strong predictor of the next assigned dose. <u>Time-varying
confounders</u> were defined as post-baseline clinical variables, other
than dose itself, that may be affected by previous dose and may also
influence subsequent dose assignment, future outcomes, or
discontinuation. Specifically, efficacy (HAMD improvement) $ΔY_{it}$,
and side-effect severity $S_{it}$ were treated as time-varying
confounders.

The dataset includes 5 active trial-by-treatment arms. Arm-level patient
counts, patient-visit counts, and observed dose levels are shown below.
Summary statistics for efficacy and side-effects are shown below.

| arm_name | n_rows | n_patients | n_treatment_weight_rows | n_censoring_weight_rows | n_msm_rows | dose_levels |
|:---|---:|---:|---:|---:|---:|:---|
| 29060_002_PAROXETINE | 1056 | 170 | 564 | 525 | 564 | 10, 20, 30, 40, 50 |
| 29060_003_IMIPRAMINE | 1105 | 241 | 684 | 682 | 691 | 20, 65, 80, 145, 210, 275 |
| 29060_003_PAROXETINE | 1192 | 241 | 758 | 756 | 758 | 10, 20, 30, 40, 50 |
| 29060_007_AMITRIPTYLINE | 95 | 13 | 52 | 52 | 52 | 20, 150, 200 |
| 29060_007_PAROXETINE | 91 | 13 | 43 | 43 | 43 | 20, 30, 40, 50 |

<small><em>Arm-level analysis summary</em></small>

Side-effect values were simulated and intentionally generated to have a
weak negative association with dose.
![](README_files/figure-gfm/side-effect-grid-1.png)<!-- -->

### Empirical support and positivity diagnostics

Before interpreting the weighted MSM, we assessed empirical support for
the observed dose process. Positivity requires that, within the
histories represented in the analysis, the dose strategies being
compared have adequate observed support. In practice, very sparse dose
levels, empty dose-transition cells, or near-deterministic dose
assignment can lead to unstable treatment weights and unreliable
extrapolation. The following grid plots summarize the observed dose
distribution by visit, previous-to-current dose transitions, and
sparse-dose grouping for each active arm.

![](README_files/figure-gfm/dose-by-visit-grid-1.png)<!-- -->

![](README_files/figure-gfm/dose-transition-grid-1.png)<!-- -->

![](README_files/figure-gfm/dose-grouping-grid-1.png)<!-- -->

| arm_name | n_patients | n_treatment_rows | n_dose_levels | min_n_per_dose | dose_levels | modelable |
|:---|---:|---:|---:|---:|:---|:---|
| 29060_002_PAROXETINE | 165 | 564 | 5 | 25 | 10, 20, 30, 40, 50 | TRUE |
| 29060_003_IMIPRAMINE | 229 | 684 | 6 | 25 | 20, 65, 80, 145, 210, 275 | TRUE |
| 29060_003_PAROXETINE | 237 | 758 | 5 | 27 | 10, 20, 30, 40, 50 | TRUE |
| 29060_007_AMITRIPTYLINE | 13 | 52 | 3 | 15 | 20, 150, 200 | FALSE |
| 29060_007_PAROXETINE | 12 | 43 | 4 | 7 | 20, 30, 40, 50 | FALSE |

<small><em>Arm-level support check used to identify modelable active
arms</em></small>

## Notation

<table class="notation-table">
<thead>
<tr>
<th>
Notation
</th>
<th>
Explanation
</th>
</tr>
</thead>
<tbody>
<tr>
<td>
$i = 1,\ldots,n$
</td>
<td>
Patient index.
</td>
</tr>
<tr>
<td>
$t = 0,1,\ldots,T_i$
</td>
<td>
Ordered visit index for patient $i$; $t = 0$ denotes baseline.
</td>
</tr>
<tr>
<td>
$\mathrm{week}_{it}$
</td>
<td>
Actual study week of visit $t$ for patient $i$.
</td>
</tr>
<tr>
<td>
$Y_{it}$
</td>
<td>
HAMD score for patient $i$ at visit $t$.
</td>
</tr>
<tr>
<td>
$\Delta Y_{it} = Y_{i0} - Y_{it}$
</td>
<td>
HAMD improvement from baseline at visit $t$. Positive values indicate
improvement.
</td>
</tr>
<tr>
<td>
$S_{it}$
</td>
<td>
Side-effect score for patient $i$ at visit $t$, measured on a 0–10
scale.
</td>
</tr>
<tr>
<td>
$d_{it}$
</td>
<td>
Observed dose actually assigned to patient $i$ at visit $t$.
</td>
</tr>
<tr>
<td>
$H_{it}$
</td>
<td>
Full observed history used in the denominator dose model.
</td>
</tr>
<tr>
<td>
$H_{it}^{*}$
</td>
<td>
Reduced observed history used in the numerator dose model.
</td>
</tr>
<tr>
<td>
$R_{it}$
</td>
<td>
Indicator for remaining observed at the next visit;<br> $R_{it} = 1$ if
patient $i$ remains observed at visit $t$, and 0 otherwise.
</td>
</tr>
<tr>
<td>
$H_{it}^{C}$
</td>
<td>
Full observed history used in the denominator censoring model.
</td>
</tr>
<tr>
<td>
$H_{it}^{C*}$
</td>
<td>
Reduced observed history used in the numerator censoring model.
</td>
</tr>
<tr>
<td>
$\hat{p}_{it}^{D}$
</td>
<td>
Fitted denominator probability of the observed dose.
</td>
</tr>
<tr>
<td>
$\hat{q}_{it}^{D}$
</td>
<td>
Fitted numerator probability of the observed dose.
</td>
</tr>
<tr>
<td>
$\mathrm{SW}_{it}^{D}$
</td>
<td>
Visit-specific stabilized dose weight.
</td>
</tr>
<tr>
<td>
$\mathrm{cSW}_{it}^{D}$
</td>
<td>
Cumulative stabilized dose weight through visit $t$.
</td>
</tr>
<tr>
<td>
$\hat{p}_{it}^{C}$
</td>
<td>
Fitted denominator probability of remaining observed.
</td>
</tr>
<tr>
<td>
$\hat{q}_{it}^{C}$
</td>
<td>
Fitted numerator probability of remaining observed.
</td>
</tr>
<tr>
<td>
$\mathrm{SW}_{it}^{C}$
</td>
<td>
Visit-specific stabilized censoring weight for patient $i$.
</td>
</tr>
<tr>
<td>
$\mathrm{cSW}_{it}^{C}$
</td>
<td>
Cumulative stabilized censoring weight through visit $t$ for patient
$i$.
</td>
</tr>
<tr>
<td>
$\mathrm{SW}_{it}^{\mathrm{total}}$
</td>
<td>
Final stabilized weight for patient $i$ at visit $t$.
</td>
</tr>
</tbody>
</table>

## Step 1 - Inverse probability of dose and censoring weighting (IPDCW)

### Inverse probability of dose weights (IPTW)

The first step for estimating marginal structural models (MSMs) was to
construct stabilized inverse probability of dose weights that adjust for
time-varying confounding affected by prior dose, following Robins,
Hernán and Brumback (2000).

The **denominator model** represents the full dose assignment mechanism
and estimates the probability of receiving the observed dose conditional
on past dose and the full history of time-varying covariates. The
**numerator model** is a simplified version of the dose assignment
mechanism and includes only time and past dose history. Its purpose is
to stabilize the weights and reduce variance without reintroducing
confounding. The stabilized weight at visit $t$ is defined as the ratio
of the probability of receiving the observed dose estimated from the
numerator model to the probability from the denominator model.

#### Dose weight denominator

The **denominator model** the dose weight was fitted at the
patient-visit level using ordinal logistic regression, treating dose as
an ordered categorical variable (20\<30\<40\<50). For each threshold
d∈{20,30,40}, we modelled the cumulative probability of receiving a dose
less than or equal to d.

``` math
\begin{aligned}
\mathrm{logit}
\left\{
\Pr(D_{it} \le c \mid H_{it})
\right\}
&=
\log
\left[
\frac{
\Pr(D_{it} \le c \mid H_{it})
}{
1 - \Pr(D_{it} \le c \mid H_{it})
}
\right] \\
&=
\alpha_{0c}
+ \alpha_1 f(\mathrm{week}_{it})
+ \alpha_2 \Delta Y_{it}
+ \alpha_3 S_{it}
+ \alpha_4 d_{it-1} \\
&\quad
+ \alpha_5 d_{it-1}\Delta Y_{it}
+ \alpha_6 d_{it-1}S_{it}
+ \alpha_7 Y_{i0}.
\end{aligned}
```

where

``` math
H_{it}
=
\left\{
\mathrm{week}_{it},
\Delta Y_{it},
S_{it},
d_{it-1},
d_{it-1}S_{it},
d_{it-1}\Delta Y_{it},
Y_{i0}
\right\}.
```

is the observed history available for patient $i$ before the dose
assignment at visit $t$. The term $`f(\mathrm{week}_{it})`$ denotes a
restricted cubic spline function of actual study week, specified with
three knots at the 10th, 50th, and 90th percentiles of the observed
visit distribution. In the ordinal dose-assignment model, the current
dose $D_{it}$ was treated as an ordered categorical variable, and
previous dose-history was represented categorically through $d_{it-1}$.
We included interactions between previous dose, HAMD improvement, and
side-effects in the denominator dose model. This allowed the
dose-assignment model to reflect that the same clinical information may
lead to different dose decisions depending on the dose the patient was
already taking. For example, poor improvement may lead to dose
escalation for a patient on a low dose, whereas side effects may lead to
dose reduction for a patient already on a high dose.

|                             |  Value | Std. Error | t value | p_value |
|:----------------------------|-------:|-----------:|--------:|--------:|
| rms::rcs(visit, 3)visit     | -0.002 |      0.107 |  -0.016 |   0.987 |
| rms::rcs(visit, 3)visit’    | -0.057 |      0.110 |  -0.515 |   0.607 |
| delta_outcome               | -0.048 |      0.041 |  -1.159 |   0.246 |
| side.effects                | -0.677 |      0.129 |  -5.239 |   0.000 |
| dose_lag1_f20               |  0.606 |      1.005 |   0.603 |   0.547 |
| dose_lag1_f30               |  0.714 |      1.014 |   0.704 |   0.481 |
| dose_lag1_f40               |  0.596 |      0.995 |   0.599 |   0.549 |
| dose_lag1_f50               |  1.730 |      1.003 |   1.726 |   0.084 |
| outcome_0                   |  0.005 |      0.021 |   0.240 |   0.810 |
| delta_outcome:dose_lag1_f20 |  0.053 |      0.044 |   1.218 |   0.223 |
| delta_outcome:dose_lag1_f30 |  0.059 |      0.044 |   1.347 |   0.178 |
| delta_outcome:dose_lag1_f40 |  0.041 |      0.043 |   0.954 |   0.340 |
| delta_outcome:dose_lag1_f50 |  0.044 |      0.045 |   0.967 |   0.333 |
| side.effects:dose_lag1_f20  |  0.204 |      0.134 |   1.527 |   0.127 |
| side.effects:dose_lag1_f30  |  0.281 |      0.138 |   2.038 |   0.042 |
| side.effects:dose_lag1_f40  |  0.385 |      0.139 |   2.775 |   0.006 |
| side.effects:dose_lag1_f50  |  0.193 |      0.141 |   1.366 |   0.172 |
| 10\|20                      | -5.684 |      1.168 |  -4.865 |   0.000 |
| 20\|30                      | -1.916 |      1.122 |  -1.707 |   0.088 |
| 30\|40                      | -0.469 |      1.121 |  -0.418 |   0.676 |
| 40\|50                      |  0.995 |      1.125 |   0.884 |   0.377 |

<small><em>Coefficient table for the dose weight denominator:
29060_002_PAROXETINE</em></small>

|                              |  Value | Std. Error | t value | p_value |
|:-----------------------------|-------:|-----------:|--------:|--------:|
| rms::rcs(visit, 3)visit      | -0.076 |      0.125 |  -0.609 |   0.543 |
| rms::rcs(visit, 3)visit’     |  0.171 |      0.162 |   1.060 |   0.289 |
| delta_outcome                |  0.026 |      0.012 |   2.207 |   0.027 |
| side.effects                 | -0.549 |      0.048 | -11.500 |   0.000 |
| dose_lag1_f65                | -0.362 |      0.912 |  -0.397 |   0.692 |
| dose_lag1_f80                |  1.112 |      0.796 |   1.397 |   0.163 |
| dose_lag1_f145               | -0.505 |      0.406 |  -1.245 |   0.213 |
| dose_lag1_f210               |  0.992 |      0.464 |   2.138 |   0.032 |
| dose_lag1_f275               |  3.049 |      0.714 |   4.268 |   0.000 |
| outcome_0                    | -0.018 |      0.020 |  -0.935 |   0.350 |
| delta_outcome:dose_lag1_f65  | -0.044 |      0.029 |  -1.515 |   0.130 |
| delta_outcome:dose_lag1_f80  | -0.058 |      0.038 |  -1.531 |   0.126 |
| delta_outcome:dose_lag1_f145 | -0.014 |      0.016 |  -0.838 |   0.402 |
| delta_outcome:dose_lag1_f210 | -0.022 |      0.023 |  -0.945 |   0.345 |
| delta_outcome:dose_lag1_f275 | -0.035 |      0.042 |  -0.836 |   0.403 |
| side.effects:dose_lag1_f65   |  0.218 |      0.109 |   2.004 |   0.045 |
| side.effects:dose_lag1_f80   |  0.090 |      0.100 |   0.901 |   0.368 |
| side.effects:dose_lag1_f145  |  0.210 |      0.059 |   3.560 |   0.000 |
| side.effects:dose_lag1_f210  |  0.015 |      0.069 |   0.222 |   0.824 |
| side.effects:dose_lag1_f275  | -0.218 |      0.102 |  -2.140 |   0.032 |
| 20\|65                       | -4.001 |      0.694 |  -5.763 |   0.000 |
| 65\|80                       | -3.687 |      0.692 |  -5.327 |   0.000 |
| 80\|145                      | -3.418 |      0.690 |  -4.954 |   0.000 |
| 145\|210                     | -1.202 |      0.677 |  -1.776 |   0.076 |
| 210\|275                     |  0.833 |      0.681 |   1.223 |   0.221 |

<small><em>Coefficient table for the dose weight denominator:
29060_003_IMIPRAMINE</em></small>

|                             |  Value | Std. Error | t value | p_value |
|:----------------------------|-------:|-----------:|--------:|--------:|
| rms::rcs(visit, 3)visit     |  0.028 |      0.110 |   0.254 |   0.799 |
| rms::rcs(visit, 3)visit’    | -0.119 |      0.185 |  -0.645 |   0.519 |
| delta_outcome               |  0.054 |      0.047 |   1.156 |   0.248 |
| side.effects                | -0.787 |      0.109 |  -7.190 |   0.000 |
| dose_lag1_f20               | -1.792 |      0.820 |  -2.184 |   0.029 |
| dose_lag1_f30               | -1.698 |      0.819 |  -2.073 |   0.038 |
| dose_lag1_f40               | -0.813 |      0.822 |  -0.989 |   0.323 |
| dose_lag1_f50               |  1.105 |      0.897 |   1.231 |   0.218 |
| outcome_0                   |  0.013 |      0.018 |   0.730 |   0.465 |
| delta_outcome:dose_lag1_f20 | -0.042 |      0.049 |  -0.870 |   0.385 |
| delta_outcome:dose_lag1_f30 | -0.045 |      0.049 |  -0.916 |   0.360 |
| delta_outcome:dose_lag1_f40 | -0.041 |      0.050 |  -0.832 |   0.405 |
| delta_outcome:dose_lag1_f50 | -0.070 |      0.052 |  -1.338 |   0.181 |
| side.effects:dose_lag1_f20  |  0.417 |      0.113 |   3.697 |   0.000 |
| side.effects:dose_lag1_f30  |  0.542 |      0.115 |   4.708 |   0.000 |
| side.effects:dose_lag1_f40  |  0.445 |      0.117 |   3.793 |   0.000 |
| side.effects:dose_lag1_f50  |  0.141 |      0.130 |   1.080 |   0.280 |
| 10\|20                      | -6.921 |      0.968 |  -7.149 |   0.000 |
| 20\|30                      | -3.313 |      0.932 |  -3.554 |   0.000 |
| 30\|40                      | -1.500 |      0.927 |  -1.618 |   0.106 |
| 40\|50                      |  0.154 |      0.925 |   0.166 |   0.868 |

<small><em>Coefficient table for the dose weight denominator:
29060_003_PAROXETINE</em></small>

A **positive coefficient** means that higher values of that variable are
associated with a higher probability of receiving a higher dose, whereas
a **negative coefficient** means that higher values of that variable are
associated with a higher probability of receiving a lower dose. After
fitting the ordinal model, predicted probabilities were obtained for
each dose category (see Appendix). The denominator probability used in
the dose weight was the fitted probability of the dose actually
received:

``` math
\hat{p}_{it}^{D}
=
\widehat{\Pr}\!\left(D_{it} = d_{it} \mid H_{it}\right).
```

![](README_files/figure-gfm/iptw-dose-model-profile-grid-1.png)<!-- -->

![](README_files/figure-gfm/iptw-dose-model-observed-grid-1.png)<!-- -->

#### Dose weight numerator

The **numerator model** is like the denominator model, except that all
potential time-dependent confounders were excluded, namely improvement
from baseline and side-effect severity. Let

``` math
H_{it}^{*}
=
\left\{
\mathrm{week}_{it},
d_{it-1},
Y_{i0}
\right\}
```

denotes the reduced observed history.

``` math
\begin{aligned}
\mathrm{logit}
\left\{
\Pr\left(D_{it} \le c \mid H_{it}^{*}\right)
\right\}
&=
\log
\left[
\frac{
\Pr\left(D_{it} \le c \mid H_{it}^{*}\right)
}{
1 - \Pr\left(D_{it} \le c \mid H_{it}^{*}\right)
}
\right] \\
&=
\beta_{0c}
+ \beta_{1} f(\mathrm{week}_{it})
+ \beta_{2} d_{it-1}
+ \beta_{3} Y_{i0}.
\end{aligned}
```

|                          |  Value | Std. Error | t value | p_value |
|:-------------------------|-------:|-----------:|--------:|--------:|
| rms::rcs(visit, 3)visit  |  0.017 |      0.099 |   0.177 |   0.860 |
| rms::rcs(visit, 3)visit’ | -0.095 |      0.105 |  -0.909 |   0.363 |
| dose_lag1_f20            |  2.266 |      0.454 |   4.995 |   0.000 |
| dose_lag1_f30            |  2.821 |      0.464 |   6.075 |   0.000 |
| dose_lag1_f40            |  3.244 |      0.468 |   6.931 |   0.000 |
| dose_lag1_f50            |  3.830 |      0.477 |   8.027 |   0.000 |
| outcome_0                |  0.013 |      0.017 |   0.738 |   0.461 |
| 10\|20                   | -0.412 |      0.668 |  -0.617 |   0.537 |
| 20\|30                   |  2.510 |      0.689 |   3.642 |   0.000 |
| 30\|40                   |  3.530 |      0.694 |   5.088 |   0.000 |
| 40\|50                   |  4.625 |      0.700 |   6.609 |   0.000 |

<small><em>Coefficient table for the dose weight numerator:
29060_002_PAROXETINE</em></small>

|                          |  Value | Std. Error | t value | p_value |
|:-------------------------|-------:|-----------:|--------:|--------:|
| rms::rcs(visit, 3)visit  | -0.025 |      0.113 |  -0.225 |   0.822 |
| rms::rcs(visit, 3)visit’ |  0.114 |      0.147 |   0.779 |   0.436 |
| dose_lag1_f65            |  0.482 |      0.301 |   1.600 |   0.110 |
| dose_lag1_f80            |  1.034 |      0.345 |   2.997 |   0.003 |
| dose_lag1_f145           |  0.855 |      0.180 |   4.763 |   0.000 |
| dose_lag1_f210           |  1.569 |      0.216 |   7.271 |   0.000 |
| dose_lag1_f275           |  2.028 |      0.319 |   6.358 |   0.000 |
| outcome_0                | -0.006 |      0.016 |  -0.354 |   0.724 |
| 20\|65                   | -0.141 |      0.522 |  -0.270 |   0.787 |
| 65\|80                   |  0.058 |      0.522 |   0.111 |   0.911 |
| 80\|145                  |  0.228 |      0.522 |   0.437 |   0.662 |
| 145\|210                 |  1.679 |      0.525 |   3.196 |   0.001 |
| 210\|275                 |  3.174 |      0.537 |   5.914 |   0.000 |

<small><em>Coefficient table for the dose weight numerator:
29060_003_IMIPRAMINE</em></small>

|                          |  Value | Std. Error | t value | p_value |
|:-------------------------|-------:|-----------:|--------:|--------:|
| rms::rcs(visit, 3)visit  |  0.073 |      0.104 |   0.701 |   0.483 |
| rms::rcs(visit, 3)visit’ | -0.166 |      0.179 |  -0.928 |   0.353 |
| dose_lag1_f20            |  0.068 |      0.447 |   0.152 |   0.879 |
| dose_lag1_f30            |  0.884 |      0.447 |   1.977 |   0.048 |
| dose_lag1_f40            |  1.435 |      0.458 |   3.136 |   0.002 |
| dose_lag1_f50            |  1.713 |      0.486 |   3.528 |   0.000 |
| outcome_0                |  0.014 |      0.015 |   0.922 |   0.357 |
| 10\|20                   | -2.241 |      0.639 |  -3.508 |   0.000 |
| 20\|30                   |  0.753 |      0.620 |   1.214 |   0.225 |
| 30\|40                   |  2.145 |      0.622 |   3.448 |   0.001 |
| 40\|50                   |  3.407 |      0.629 |   5.412 |   0.000 |

<small><em>Coefficient table for the dose weight numerator:
29060_003_PAROXETINE</em></small>

The fitted probability corresponding to the observed dose $d_{it}$ was
extracted as

``` math
\hat{q}_{it}^{D}
=
\widehat{\Pr}\left(D_{it} = d_{it} \mid H_{it}^{*}\right).
```

#### Stabilized dose weights

The visit-specific stabilized dose weight was then calculated as

``` math
\mathrm{SW}_{it}^{D}
=
\frac{\hat{q}_{it}^{D}}{\hat{p}_{it}^{D}}.
```

The cumulative stabilized dose weight through visit $t$ was obtained by
multiplying the visit-specific weights from the first post-baseline
visit up to the current visit:

``` math
\mathrm{cSW}_{it}^{D}
=
\prod_{s=1}^{t} \mathrm{SW}_{is}^{D}.
```

| arm_name | n | n_patients | mean_SW_treatment | sd_SW_treatment | min_SW_treatment | p1_SW_treatment | p50_SW_treatment | p99_SW_treatment | max_SW_treatment | mean_cSW_treatment | sd_cSW_treatment | min_cSW_treatment | p1_cSW_treatment | p50_cSW_treatment | p99_cSW_treatment | max_cSW_treatment | ESS_cSW_treatment |
|:---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 29060_002_PAROXETINE | 564 | 165 | 1.007 | 1.251 | 0.050 | 0.217 | 0.737 | 5.275 | 22.825 | 0.893 | 1.032 | 0.020 | 0.047 | 0.600 | 5.458 | 8.406 | 241.575 |
| 29060_003_IMIPRAMINE | 684 | 229 | 1.153 | 3.099 | 0.119 | 0.140 | 0.690 | 11.093 | 68.989 | 1.482 | 5.394 | 0.005 | 0.029 | 0.587 | 23.512 | 74.955 | 48.058 |
| 29060_003_PAROXETINE | 758 | 237 | 1.120 | 2.587 | 0.060 | 0.111 | 0.791 | 5.139 | 61.007 | 2.064 | 19.838 | 0.003 | 0.016 | 0.684 | 10.899 | 486.519 | 8.126 |

<small><em>Summary of visit-specific and cumulative stabilized dose
weights</em></small>

Note: The ordinal dose model assumes proportional odds, meaning that the
covariates have a common effect across the cumulative dose-category
thresholds. In other words, covariates shift the probability toward
higher or lower dose categories similarly across cut-points. This
assumption will be examined in sensitivity analysis by comparing the
results with those from a multinomial dose model; see Appendix.

### Inverse probability of censoring weights (IPCW)

#### Censoring weight denominator

The **denominator model** described the dropout process conditional on
the observed history. It included current dose, current and previous
HAMD improvement, current and previous side-effect scores, and baseline
covariates:

``` math
\begin{aligned}
\mathrm{logit}
\left\{
\Pr\left(R_{it+1} = 1 \mid H_{it}^{C}\right)
\right\}
&=
\gamma_{0}
+ \gamma_{1} f(\mathrm{week}_{it})
+ \gamma_{2} d_{it}
+ \gamma_{3} \Delta Y_{it}
+ \gamma_{4} \Delta Y_{it-1} \\
&\quad
+ \gamma_{5} S_{it}
+ \gamma_{6} S_{it-1}
+ \gamma_{7} Y_{i0}
+ \gamma_{8} \mathrm{age}_{i}
+ \gamma_{9} \mathrm{sex}_{i}.
\end{aligned}
```

where

``` math
H_{it}^{C}
=
\left\{
\mathrm{week}_{it},
d_{it},
\Delta Y_{it},
\Delta Y_{it-1},
S_{it},
S_{it-1},
Y_{i0},
\mathrm{age}_{i},
\mathrm{sex}_{i}
\right\}.
```

|                    | Estimate | Std. Error | z value | Pr(\>\|z\|) | p_value |
|:-------------------|---------:|-----------:|--------:|------------:|--------:|
| (Intercept)        |    5.015 |      1.635 |   3.067 |       0.002 |   0.002 |
| visit              |   -0.898 |      0.108 |  -8.318 |       0.000 |   0.000 |
| dose_current_f20   |    0.002 |      0.718 |   0.002 |       0.998 |   0.998 |
| dose_current_f30   |    0.002 |      0.773 |   0.002 |       0.998 |   0.998 |
| dose_current_f40   |    0.282 |      0.833 |   0.339 |       0.735 |   0.735 |
| dose_current_f50   |    0.304 |      0.872 |   0.349 |       0.727 |   0.727 |
| delta_outcome      |   -0.022 |      0.019 |  -1.128 |       0.259 |   0.259 |
| delta_outcome_lag1 |    0.022 |      0.019 |   1.135 |       0.256 |   0.256 |
| side.effects       |   -0.053 |      0.057 |  -0.928 |       0.354 |   0.354 |
| side.effects_lag1  |    0.097 |      0.044 |   2.183 |       0.029 |   0.029 |
| outcome_0          |   -0.013 |      0.043 |  -0.301 |       0.763 |   0.763 |
| age                |    0.033 |      0.015 |   2.137 |       0.033 |   0.033 |
| sexM               |   -0.061 |      0.308 |  -0.197 |       0.844 |   0.844 |

<small><em>Coefficient table for the censoring weight denominator:
29060_002_PAROXETINE</em></small>

|                    | Estimate | Std. Error | z value | Pr(\>\|z\|) | p_value |
|:-------------------|---------:|-----------:|--------:|------------:|--------:|
| (Intercept)        |    5.744 |      1.140 |   5.037 |       0.000 |   0.000 |
| visit              |   -1.189 |      0.101 | -11.722 |       0.000 |   0.000 |
| dose_current_f65   |    1.153 |      0.690 |   1.671 |       0.095 |   0.095 |
| dose_current_f80   |    0.903 |      0.766 |   1.179 |       0.238 |   0.238 |
| dose_current_f145  |    0.774 |      0.362 |   2.142 |       0.032 |   0.032 |
| dose_current_f210  |    0.349 |      0.438 |   0.798 |       0.425 |   0.425 |
| dose_current_f275  |    0.610 |      0.549 |   1.111 |       0.267 |   0.267 |
| delta_outcome      |    0.006 |      0.018 |   0.341 |       0.733 |   0.733 |
| delta_outcome_lag1 |    0.023 |      0.017 |   1.344 |       0.179 |   0.179 |
| side.effects       |    0.071 |      0.049 |   1.442 |       0.149 |   0.149 |
| side.effects_lag1  |   -0.035 |      0.032 |  -1.092 |       0.275 |   0.275 |
| outcome_0          |   -0.048 |      0.035 |  -1.344 |       0.179 |   0.179 |
| age                |    0.006 |      0.010 |   0.622 |       0.534 |   0.534 |
| sexM               |    0.132 |      0.238 |   0.553 |       0.580 |   0.580 |

<small><em>Coefficient table for the censoring weight denominator:
29060_003_IMIPRAMINE</em></small>

|                    | Estimate | Std. Error | z value | Pr(\>\|z\|) | p_value |
|:-------------------|---------:|-----------:|--------:|------------:|--------:|
| (Intercept)        |    6.462 |      1.265 |   5.108 |       0.000 |   0.000 |
| visit              |   -1.107 |      0.089 | -12.393 |       0.000 |   0.000 |
| dose_current_f20   |   -0.549 |      0.717 |  -0.766 |       0.444 |   0.444 |
| dose_current_f30   |   -0.043 |      0.725 |  -0.059 |       0.953 |   0.953 |
| dose_current_f40   |   -0.440 |      0.763 |  -0.576 |       0.564 |   0.564 |
| dose_current_f50   |    0.385 |      0.809 |   0.476 |       0.634 |   0.634 |
| delta_outcome      |    0.022 |      0.017 |   1.300 |       0.194 |   0.194 |
| delta_outcome_lag1 |    0.009 |      0.015 |   0.597 |       0.550 |   0.550 |
| side.effects       |    0.029 |      0.040 |   0.715 |       0.475 |   0.475 |
| side.effects_lag1  |    0.032 |      0.034 |   0.962 |       0.336 |   0.336 |
| outcome_0          |   -0.043 |      0.033 |  -1.321 |       0.187 |   0.187 |
| age                |    0.004 |      0.011 |   0.395 |       0.693 |   0.693 |
| sexM               |   -0.123 |      0.233 |  -0.526 |       0.599 |   0.599 |

<small><em>Coefficient table for the censoring weight denominator:
29060_003_PAROXETINE</em></small>

The denominator predicted probability is

``` math
\hat{p}_{it}^{C}
=
\widehat{\Pr}
\left(
R_{it+1} = 1 \mid H_{it}^{C}
\right).
```

#### Censoring weight numerator

The **numerator model** was used only to stabilize the censoring
weights. It excluded the time-varying confounders, while retaining
visit, current dose, and baseline covariates:

``` math
\begin{aligned}
\mathrm{logit}
\left\{
\Pr\left(R_{it+1} = 1 \mid H_{it}^{C*}\right)
\right\}
&=
\delta_{0}
+ \delta_{1} \mathrm{week}_{it}
+ \delta_{2} d_{it}
+ \delta_{3} Y_{i0}
+ \delta_{4} \mathrm{age}_{i}
+ \delta_{5} \mathrm{sex}_{i}.
\end{aligned}
```

where

``` math
H_{it}^{C*}
=
\left\{
\mathrm{week}_{it},
d_{it},
Y_{i0},
\mathrm{age}_{i},
\mathrm{sex}_{i}
\right\}.
```

|                  | Estimate | Std. Error | z value | Pr(\>\|z\|) | p_value |
|:-----------------|---------:|-----------:|--------:|------------:|--------:|
| (Intercept)      |    5.076 |      1.392 |   3.648 |       0.000 |   0.000 |
| visit            |   -0.856 |      0.097 |  -8.798 |       0.000 |   0.000 |
| dose_current_f20 |   -0.343 |      0.701 |  -0.489 |       0.625 |   0.625 |
| dose_current_f30 |   -0.283 |      0.727 |  -0.389 |       0.697 |   0.697 |
| dose_current_f40 |    0.022 |      0.751 |   0.030 |       0.976 |   0.976 |
| dose_current_f50 |    0.232 |      0.762 |   0.305 |       0.761 |   0.761 |
| outcome_0        |   -0.009 |      0.032 |  -0.282 |       0.778 |   0.778 |
| age              |    0.030 |      0.015 |   2.041 |       0.041 |   0.041 |
| sexM             |   -0.054 |      0.301 |  -0.179 |       0.858 |   0.858 |

<small><em>Coefficient table for the censoring weight numerator:
29060_002_PAROXETINE</em></small>

|                   | Estimate | Std. Error | z value | Pr(\>\|z\|) | p_value |
|:------------------|---------:|-----------:|--------:|------------:|--------:|
| (Intercept)       |    5.417 |      0.919 |   5.894 |       0.000 |   0.000 |
| visit             |   -1.122 |      0.085 | -13.157 |       0.000 |   0.000 |
| dose_current_f65  |    1.045 |      0.695 |   1.503 |       0.133 |   0.133 |
| dose_current_f80  |    0.691 |      0.729 |   0.948 |       0.343 |   0.343 |
| dose_current_f145 |    0.502 |      0.299 |   1.680 |       0.093 |   0.093 |
| dose_current_f210 |    0.013 |      0.332 |   0.038 |       0.969 |   0.969 |
| dose_current_f275 |    0.139 |      0.417 |   0.334 |       0.738 |   0.738 |
| outcome_0         |   -0.024 |      0.028 |  -0.843 |       0.399 |   0.399 |
| age               |    0.007 |      0.010 |   0.748 |       0.454 |   0.454 |
| sexM              |    0.143 |      0.236 |   0.604 |       0.546 |   0.546 |

<small><em>Coefficient table for the censoring weight numerator:
29060_003_IMIPRAMINE</em></small>

|                  | Estimate | Std. Error | z value | Pr(\>\|z\|) | p_value |
|:-----------------|---------:|-----------:|--------:|------------:|--------:|
| (Intercept)      |    6.044 |      1.105 |   5.471 |       0.000 |   0.000 |
| visit            |   -1.035 |      0.075 | -13.837 |       0.000 |   0.000 |
| dose_current_f20 |   -0.556 |      0.705 |  -0.788 |       0.431 |   0.431 |
| dose_current_f30 |   -0.117 |      0.707 |  -0.166 |       0.868 |   0.868 |
| dose_current_f40 |   -0.613 |      0.733 |  -0.836 |       0.403 |   0.403 |
| dose_current_f50 |    0.143 |      0.759 |   0.188 |       0.851 |   0.851 |
| outcome_0        |   -0.015 |      0.026 |  -0.597 |       0.551 |   0.551 |
| age              |    0.005 |      0.010 |   0.452 |       0.651 |   0.651 |
| sexM             |   -0.086 |      0.231 |  -0.375 |       0.708 |   0.708 |

<small><em>Coefficient table for the censoring weight numerator:
29060_003_PAROXETINE</em></small>

The numerator predicted probability was

``` math
\hat{q}_{it}^{C}
=
\widehat{\Pr}
\left(
R_{it+1} = 1 \mid H_{it}^{C*}
\right).
```

#### Stabilized censoring weights

The visit-specific stabilized censoring weight was then calculated as

``` math
\mathrm{SW}_{it}^{C}
=
\frac{\hat{q}_{it}^{C}}{\hat{p}_{it}^{C}}.
```

The cumulative stabilized censoring weight through visit $t$ was
obtained by multiplying the visit-specific censoring weights from the
first post-baseline interval up to the current visit:

``` math
\mathrm{cSW}_{it}^{C}
=
\prod_{s=1}^{t} \mathrm{SW}_{is}^{C}.
```

The effective sample size after weighting was calculated using the Kish
formula:

``` math
\mathrm{ESS}
=
\frac{
\left(
\sum_{i,t} \mathrm{SW}_{it}^{\mathrm{total}}
\right)^2
}{
\sum_{i,t}
\left(
\mathrm{SW}_{it}^{\mathrm{total}}
\right)^2
}.
```

| arm_name | n | n_patients | mean_SW_censoring | sd_SW_censoring | min_SW_censoring | p1_SW_censoring | p50_SW_censoring | p99_SW_censoring | max_SW_censoring | mean_cSW_censoring | sd_cSW_censoring | min_cSW_censoring | p1_cSW_censoring | p50_cSW_censoring | p99_cSW_censoring | max_cSW_censoring | ESS_cSW_censoring |
|:---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 29060_002_PAROXETINE | 525 | 163 | 1.012 | 0.124 | 0.616 | 0.731 | 1.000 | 1.490 | 2.162 | 1.012 | 0.123 | 0.631 | 0.738 | 1.000 | 1.508 | 2.132 | 517.408 |
| 29060_003_IMIPRAMINE | 682 | 229 | 1.011 | 0.115 | 0.604 | 0.738 | 0.998 | 1.476 | 2.007 | 1.012 | 0.123 | 0.559 | 0.715 | 0.998 | 1.589 | 2.032 | 672.000 |
| 29060_003_PAROXETINE | 756 | 237 | 1.009 | 0.087 | 0.660 | 0.804 | 1.001 | 1.349 | 1.620 | 1.009 | 0.102 | 0.596 | 0.727 | 1.001 | 1.400 | 1.626 | 748.430 |

<small><em>Summary of visit-specific and cumulative stabilized censoring
weights</em></small>

### Total stabilized weights and truncation

The final stabilized weight used in the marginal structural model was
obtained by multiplying the cumulative stabilized dose weight and the
cumulative stabilized censoring weight for each patient-visit row. This
accounts jointly for non-random dose titration and informative censoring
up to visit $t$:

``` math
\mathrm{SW}_{it}^{\mathrm{total}}
=
\mathrm{cSW}_{it}^{D}
\times
\mathrm{cSW}_{it}^{C}.
```

Weight stability was assessed by examining the distribution of the dose,
censoring, and total weights, including their mean, standard deviation,
selected percentiles, maximum value, and the effective sample size after
truncation. The effective sample size after weighting was calculated
using the Kish formula:

``` math
\mathrm{ESS}
=
\frac{
\left(
\sum_{i,t} \mathrm{SW}_{it}^{\mathrm{total}}
\right)^2
}{
\sum_{i,t}
\left(
\mathrm{SW}_{it}^{\mathrm{total}}
\right)^2
}.
```

The total-weight summary describes the final stabilized weights,
$\mathrm{SW}_{it}^{\mathrm{total}}$, obtained by multiplying the
cumulative dose and censoring weights.

| arm_name | n | n_patients | mean_SW_total | sd_SW_total | min_SW_total | p1_SW_total | p50_SW_total | p99_SW_total | max_SW_total | ESS_SW_total |
|:---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 29060_002_PAROXETINE | 525 | 163 | 0.911 | 1.088 | 0.020 | 0.051 | 0.616 | 5.757 | 8.781 | 216.645 |
| 29060_003_IMIPRAMINE | 682 | 229 | 1.494 | 5.386 | 0.003 | 0.026 | 0.592 | 21.814 | 74.782 | 48.773 |
| 29060_003_PAROXETINE | 756 | 237 | 1.846 | 15.156 | 0.003 | 0.015 | 0.687 | 11.971 | 349.965 | 11.071 |

<small><em>Summary of total stabilized weights</em></small>

To reduce the influence of extreme weights and improve precision, total
stabilized weights were truncated at the 1st and 99th percentiles. The
distribution of the truncated weights is summarized below.

| arm_name | n | n_patients | mean_SW_total_trunc | sd_SW_total_trunc | min_SW_total_trunc | p1_SW_total_trunc | p50_SW_total_trunc | p99_SW_total_trunc | max_SW_total_trunc | ESS_SW_total_trunc |
|:---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 29060_002_PAROXETINE | 525 | 163 | 0.891 | 0.970 | 0.051 | 0.052 | 0.616 | 5.751 | 5.757 | 240.551 |
| 29060_003_IMIPRAMINE | 682 | 229 | 1.239 | 2.879 | 0.026 | 0.026 | 0.592 | 21.787 | 21.814 | 106.649 |
| 29060_003_PAROXETINE | 756 | 237 | 0.989 | 1.461 | 0.015 | 0.015 | 0.687 | 11.850 | 11.971 | 237.723 |

<small><em>Summary of truncated total stabilized weights</em></small>

## Step 2 - Weighted repeated-measures marginal structural model MSM

After estimating the total weights for each patient-visit row, we fitted
a weighted repeated-measures marginal structural model to estimate the
effect of dose history on HAMD improvement over time. The model was
estimated using generalized estimating equations with an identity link
and Gaussian working variance. We used an independence working
correlation structure and robust sandwich standard errors clustered by
patient. The independence structure was used as a working correlation
assumption only. Inference was based on the robust variance estimator,
which accounts for within-patient correlation in repeated HAMD
measurements and provides valid standard errors even if the working
correlation structure is misspecified. The final stabilized weights,
$\mathrm{SW}_{it}^{\mathrm{total}}$, were used as observation-level
weights.

Following the dose-history structure used by Lipkovich et al., the
outcome model included recent dose history. The most recent dose
$d_{it-1}$, the dose one visit earlier $d_{it-2}$, and the dose two
visits earlier $d_{it-3}$ were represented categorically. A separate 0
category was retained where a lagged dose was not yet available. The
average earlier dose, $\bar{d}_{i<t-3}$, was kept continuous because it
can take values that are not actual dose levels, for example 25 or 33.3
mg.

The weighted MSM was specified as

``` math
\begin{aligned}
\Delta Y_{it}
&=
\eta_{0}
+ \eta_{1} f(\mathrm{week}_{it})
+ \eta_{2} Y_{i0}
+ \eta_{3} d_{it-1}
+ \eta_{4} d_{it-2} \\
&\quad
+ \eta_{5} d_{it-3}
+ \eta_{6} \overline{d}_{i\lt t-3}
+ \eta_{7} Y_{i0} f(\mathrm{week}_{it}) \\ 
&\quad 
+ \eta_{8}\!\left(d_{it-1}, f(\mathrm{week}_{it})\right)
+ \varepsilon_{it}.
\end{aligned}
```

The term $`f(\mathrm{week}_{it})`$ denotes a restricted cubic spline
function of actual study week, with three knots at the 10th, 50th, and
90th percentiles of the observed visit distribution. $Y_{i0}$ is the
baseline HAMD score. Following the dose-history structure used by
Lipkovich et al., the outcome model included recent dose history. The
most recent dose $d_{it-1}$, the dose one visit earlier $d_{it-2}$, and
the dose two visits earlier $d_{it-3}$ were represented categorically.
The average earlier dose, $\bar d_{i<t-3}$, was kept continuous because
it can take values that are not actual dose levels, for example 25 or
33.3 mg. A separate 0 category was retained where a lagged dose was not
yet available. For example, at the first post-baseline visit,
$d_{it-2}$, $d_{it-3}$, and $\bar d_{i<t-3}$ were set to 0 mg. We used
interactions between study week and the most recent dose category,
allowing the effect of the most recent dose to vary over follow-up, and
between baseline HAMD and week to allow patients with different baseline
severity to have different improvement trajectories over time.

This model estimates the mean HAMD improvement trajectory under
alternative dose histories in the weighted pseudo-population, where dose
titration and censoring are no longer driven by the measured clinical
history included in the weight models.

| Term                                   | Estimate | Std.err |   Wald | Pr(\>\|W\|) |
|:---------------------------------------|---------:|--------:|-------:|------------:|
| (Intercept)                            |  -27.417 |   8.117 | 11.409 |       0.001 |
| rms::rcs(visit, 3)visit                |    5.324 |   4.152 |  1.644 |       0.200 |
| rms::rcs(visit, 3)visit’               |   -1.778 |   7.090 |  0.063 |       0.802 |
| outcome_0                              |    0.770 |   0.246 |  9.828 |       0.002 |
| dose_lag1_f20                          |    1.787 |   5.083 |  0.124 |       0.725 |
| dose_lag1_f30                          |    9.644 |   5.131 |  3.533 |       0.060 |
| dose_lag1_f40                          |    8.210 |   4.855 |  2.860 |       0.091 |
| dose_lag1_f50                          |    2.106 |   5.204 |  0.164 |       0.686 |
| dose_lag2_f10                          |   -5.370 |   2.534 |  4.491 |       0.034 |
| dose_lag2_f20                          |   -2.294 |   1.706 |  1.809 |       0.179 |
| dose_lag2_f30                          |   -0.307 |   1.840 |  0.028 |       0.868 |
| dose_lag2_f40                          |   -0.890 |   2.158 |  0.170 |       0.680 |
| dose_lag2_f50                          |   -0.025 |   1.998 |  0.000 |       0.990 |
| dose_lag3_f10                          |   -4.607 |   2.713 |  2.883 |       0.089 |
| dose_lag3_f20                          |    0.534 |   1.727 |  0.096 |       0.757 |
| dose_lag3_f30                          |    1.137 |   2.069 |  0.302 |       0.583 |
| dose_lag3_f40                          |   -0.262 |   2.074 |  0.016 |       0.900 |
| dose_lag3_f50                          |    0.284 |   1.807 |  0.025 |       0.875 |
| avg_dose_before_lag3                   |   -0.063 |   0.031 |  4.091 |       0.043 |
| rms::rcs(visit, 3)visit:outcome_0      |    0.128 |   0.124 |  1.074 |       0.300 |
| rms::rcs(visit, 3)visit’:outcome_0     |   -0.295 |   0.213 |  1.911 |       0.167 |
| rms::rcs(visit, 3)visit:dose_lag1_f20  |   -2.514 |   2.366 |  1.129 |       0.288 |
| rms::rcs(visit, 3)visit’:dose_lag1_f20 |    3.054 |   3.805 |  0.644 |       0.422 |
| rms::rcs(visit, 3)visit:dose_lag1_f30  |   -6.149 |   2.421 |  6.449 |       0.011 |
| rms::rcs(visit, 3)visit’:dose_lag1_f30 |    8.391 |   4.089 |  4.212 |       0.040 |
| rms::rcs(visit, 3)visit:dose_lag1_f40  |   -6.434 |   2.355 |  7.464 |       0.006 |
| rms::rcs(visit, 3)visit’:dose_lag1_f40 |    9.143 |   3.915 |  5.453 |       0.020 |
| rms::rcs(visit, 3)visit:dose_lag1_f50  |   -3.318 |   2.338 |  2.014 |       0.156 |
| rms::rcs(visit, 3)visit’:dose_lag1_f50 |    3.329 |   3.777 |  0.777 |       0.378 |

<small><em>Weighted marginal structural model for HAMD improvement:
29060_002_PAROXETINE</em></small>

| Term                                    | Estimate | Std.err |  Wald | Pr(\>\|W\|) |
|:----------------------------------------|---------:|--------:|------:|------------:|
| (Intercept)                             |  -14.837 |   9.969 | 2.215 |       0.137 |
| rms::rcs(visit, 3)visit                 |   -1.396 |   4.729 | 0.087 |       0.768 |
| rms::rcs(visit, 3)visit’                |    4.116 |   5.944 | 0.479 |       0.489 |
| outcome_0                               |    0.723 |   0.354 | 4.167 |       0.041 |
| dose_lag1_f65                           |   -2.676 |   4.929 | 0.295 |       0.587 |
| dose_lag1_f80                           |    0.860 |   9.339 | 0.008 |       0.927 |
| dose_lag1_f145                          |   -5.963 |   3.426 | 3.029 |       0.082 |
| dose_lag1_f210                          |   -2.774 |   3.528 | 0.618 |       0.432 |
| dose_lag1_f275                          |    0.634 |   4.711 | 0.018 |       0.893 |
| dose_lag2_f20                           |    2.031 |   1.228 | 2.732 |       0.098 |
| dose_lag2_f65                           |   -1.919 |   2.128 | 0.813 |       0.367 |
| dose_lag2_f80                           |    0.041 |   2.097 | 0.000 |       0.984 |
| dose_lag2_f145                          |   -0.320 |   1.220 | 0.069 |       0.793 |
| dose_lag2_f210                          |   -1.182 |   1.733 | 0.465 |       0.495 |
| dose_lag2_f275                          |   -0.885 |   1.730 | 0.262 |       0.609 |
| dose_lag3_f20                           |    1.327 |   2.639 | 0.253 |       0.615 |
| dose_lag3_f65                           |   -2.075 |   2.030 | 1.044 |       0.307 |
| dose_lag3_f80                           |    1.411 |   2.447 | 0.332 |       0.564 |
| dose_lag3_f145                          |    2.507 |   1.662 | 2.276 |       0.131 |
| dose_lag3_f210                          |    2.034 |   1.561 | 1.699 |       0.192 |
| dose_lag3_f275                          |    3.822 |   2.576 | 2.202 |       0.138 |
| avg_dose_before_lag3                    |   -0.002 |   0.007 | 0.086 |       0.770 |
| rms::rcs(visit, 3)visit:outcome_0       |    0.078 |   0.162 | 0.233 |       0.629 |
| rms::rcs(visit, 3)visit’:outcome_0      |   -0.091 |   0.200 | 0.205 |       0.651 |
| rms::rcs(visit, 3)visit:dose_lag1_f65   |    3.875 |   2.624 | 2.180 |       0.140 |
| rms::rcs(visit, 3)visit’:dose_lag1_f65  |   -6.119 |   3.862 | 2.510 |       0.113 |
| rms::rcs(visit, 3)visit:dose_lag1_f80   |    1.363 |   4.134 | 0.109 |       0.742 |
| rms::rcs(visit, 3)visit’:dose_lag1_f80  |   -2.472 |   5.077 | 0.237 |       0.626 |
| rms::rcs(visit, 3)visit:dose_lag1_f145  |    3.590 |   1.622 | 4.898 |       0.027 |
| rms::rcs(visit, 3)visit’:dose_lag1_f145 |   -4.571 |   2.374 | 3.707 |       0.054 |
| rms::rcs(visit, 3)visit:dose_lag1_f210  |    2.159 |   1.627 | 1.761 |       0.184 |
| rms::rcs(visit, 3)visit’:dose_lag1_f210 |   -3.448 |   1.953 | 3.118 |       0.077 |
| rms::rcs(visit, 3)visit:dose_lag1_f275  |    1.077 |   2.232 | 0.233 |       0.630 |
| rms::rcs(visit, 3)visit’:dose_lag1_f275 |   -2.195 |   3.389 | 0.419 |       0.517 |

<small><em>Weighted marginal structural model for HAMD improvement:
29060_003_IMIPRAMINE</em></small>

| Term                                   | Estimate | Std.err |   Wald | Pr(\>\|W\|) |
|:---------------------------------------|---------:|--------:|-------:|------------:|
| (Intercept)                            |  -19.701 |   8.493 |  5.381 |       0.020 |
| rms::rcs(visit, 3)visit                |    0.353 |   3.292 |  0.011 |       0.915 |
| rms::rcs(visit, 3)visit’               |    1.753 |   5.813 |  0.091 |       0.763 |
| outcome_0                              |    1.194 |   0.249 | 23.073 |       0.000 |
| dose_lag1_f20                          |  -10.952 |   6.762 |  2.624 |       0.105 |
| dose_lag1_f30                          |  -10.514 |   6.600 |  2.538 |       0.111 |
| dose_lag1_f40                          |   -8.432 |   7.101 |  1.410 |       0.235 |
| dose_lag1_f50                          |  -14.975 |   7.337 |  4.165 |       0.041 |
| dose_lag2_f10                          |    0.017 |   2.541 |  0.000 |       0.995 |
| dose_lag2_f20                          |   -0.151 |   1.403 |  0.012 |       0.914 |
| dose_lag2_f30                          |    0.329 |   1.323 |  0.062 |       0.804 |
| dose_lag2_f40                          |    1.919 |   1.477 |  1.687 |       0.194 |
| dose_lag2_f50                          |    1.859 |   1.421 |  1.712 |       0.191 |
| dose_lag3_f10                          |    2.703 |   2.387 |  1.283 |       0.257 |
| dose_lag3_f20                          |   -0.908 |   1.092 |  0.691 |       0.406 |
| dose_lag3_f30                          |   -0.719 |   1.009 |  0.508 |       0.476 |
| dose_lag3_f40                          |   -0.422 |   1.413 |  0.089 |       0.765 |
| dose_lag3_f50                          |   -2.218 |   1.755 |  1.597 |       0.206 |
| avg_dose_before_lag3                   |    0.031 |   0.028 |  1.227 |       0.268 |
| rms::rcs(visit, 3)visit:outcome_0      |   -0.101 |   0.102 |  0.982 |       0.322 |
| rms::rcs(visit, 3)visit’:outcome_0     |    0.148 |   0.167 |  0.792 |       0.373 |
| rms::rcs(visit, 3)visit:dose_lag1_f20  |    5.026 |   2.901 |  3.003 |       0.083 |
| rms::rcs(visit, 3)visit’:dose_lag1_f20 |   -6.278 |   5.553 |  1.278 |       0.258 |
| rms::rcs(visit, 3)visit:dose_lag1_f30  |    4.734 |   2.758 |  2.946 |       0.086 |
| rms::rcs(visit, 3)visit’:dose_lag1_f30 |   -7.189 |   5.306 |  1.836 |       0.175 |
| rms::rcs(visit, 3)visit:dose_lag1_f40  |    3.886 |   2.957 |  1.728 |       0.189 |
| rms::rcs(visit, 3)visit’:dose_lag1_f40 |   -6.135 |   5.577 |  1.210 |       0.271 |
| rms::rcs(visit, 3)visit:dose_lag1_f50  |    6.758 |   3.058 |  4.882 |       0.027 |
| rms::rcs(visit, 3)visit’:dose_lag1_f50 |  -10.299 |   5.700 |  3.265 |       0.071 |

<small><em>Weighted marginal structural model for HAMD improvement:
29060_003_PAROXETINE</em></small>

Interpretation: If a dose coefficient is positive: then higher previous
dose is associated with greater HAMD improvement. If a dose coefficient
is negative: then higher previous dose is associated with lower HAMD
improvement.

## Step 3 - Predictions

After fitting the weighted MSM, we predicted HAMD improvement
trajectories under sustained dose strategies, setting the dose-history
variables to the values corresponding to each strategy. Predictions were
obtained at the mean baseline HAMD score. Pointwise 95% confidence
intervals were calculated using the robust sandwich covariance matrix
from the fitted GEE. These intervals reflect uncertainty in the fitted
MSM coefficients; they do not incorporate additional uncertainty from
estimating the treatment and censoring weights unless a full bootstrap
is used.

![](README_files/figure-gfm/strategy-prediction-grid-1.png)<!-- -->

| weight_model | arm_name | trial_name | week | dose_strategy | strategy_dose | active_n_patients | active_n_visit_rows | placebo_n_patients | placebo_n_visit_rows | predicted_active_improvement | SE_active_prediction | observed_placebo_improvement | SE_observed_placebo | dose_vs_placebo_difference | SE | lower_95 | upper_95 |
|:---|:---|:---|---:|:---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| Ordinal IPTW | 29060_002_PAROXETINE | 29060_002 | 1 | Target 10 mg | 10 | 86 | 86 | 73 | 73 | 2.091 | 2.791 | 3.753 | 0.839 | -1.662 | 2.915 | -7.375 | 4.051 |
| Ordinal IPTW | 29060_002_PAROXETINE | 29060_002 | 2 | Target 10 mg | 10 | 91 | 91 | 99 | 99 | 5.231 | 2.126 | 4.322 | 0.827 | 0.909 | 2.281 | -3.561 | 5.380 |
| Ordinal IPTW | 29060_002_PAROXETINE | 29060_002 | 3 | Target 10 mg | 10 | 89 | 89 | 89 | 89 | 6.881 | 3.192 | 7.746 | 0.933 | -0.865 | 3.326 | -7.384 | 5.653 |
| Ordinal IPTW | 29060_002_PAROXETINE | 29060_002 | 4 | Target 10 mg | 10 | 81 | 81 | 86 | 86 | 10.936 | 3.068 | 9.505 | 1.155 | 1.431 | 3.278 | -4.993 | 7.856 |
| Ordinal IPTW | 29060_002_PAROXETINE | 29060_002 | 5 | Target 10 mg | 10 | 50 | 50 | 37 | 37 | 12.562 | 2.622 | 9.757 | 1.674 | 2.805 | 3.111 | -3.292 | 8.902 |
| Ordinal IPTW | 29060_002_PAROXETINE | 29060_002 | 6 | Target 10 mg | 10 | 55 | 55 | 49 | 49 | 12.569 | 2.248 | 10.814 | 1.163 | 1.754 | 2.531 | -3.206 | 6.715 |
| Ordinal IPTW | 29060_002_PAROXETINE | 29060_002 | 7 | Target 10 mg | 10 | 62 | 62 | 45 | 45 | 11.765 | 2.638 | 12.941 | 1.436 | -1.176 | 3.003 | -7.063 | 4.710 |
| Ordinal IPTW | 29060_002_PAROXETINE | 29060_002 | 8 | Target 10 mg | 10 | 11 | 11 | 6 | 6 | 10.827 | 3.697 | 19.861 | 4.748 | -9.035 | 6.017 | -20.828 | 2.758 |
| Ordinal IPTW | 29060_002_PAROXETINE | 29060_002 | 1 | Target 20 mg | 20 | 86 | 86 | 73 | 73 | 1.365 | 1.561 | 3.753 | 0.839 | -2.388 | 1.772 | -5.862 | 1.085 |
| Ordinal IPTW | 29060_002_PAROXETINE | 29060_002 | 2 | Target 20 mg | 20 | 91 | 91 | 99 | 99 | 5.151 | 1.409 | 4.322 | 0.827 | 0.830 | 1.633 | -2.372 | 4.031 |
| Ordinal IPTW | 29060_002_PAROXETINE | 29060_002 | 3 | Target 20 mg | 20 | 89 | 89 | 89 | 89 | 9.388 | 1.649 | 7.746 | 0.933 | 1.642 | 1.895 | -2.072 | 5.356 |
| Ordinal IPTW | 29060_002_PAROXETINE | 29060_002 | 4 | Target 20 mg | 20 | 81 | 81 | 86 | 86 | 12.415 | 1.506 | 9.505 | 1.155 | 2.910 | 1.898 | -0.810 | 6.630 |
| Ordinal IPTW | 29060_002_PAROXETINE | 29060_002 | 5 | Target 20 mg | 20 | 50 | 50 | 37 | 37 | 13.775 | 1.318 | 9.757 | 1.674 | 4.018 | 2.130 | -0.158 | 8.193 |
| Ordinal IPTW | 29060_002_PAROXETINE | 29060_002 | 6 | Target 20 mg | 20 | 55 | 55 | 49 | 49 | 14.025 | 1.264 | 10.814 | 1.163 | 3.211 | 1.718 | -0.156 | 6.577 |
| Ordinal IPTW | 29060_002_PAROXETINE | 29060_002 | 7 | Target 20 mg | 20 | 62 | 62 | 45 | 45 | 13.719 | 1.604 | 12.941 | 1.436 | 0.778 | 2.153 | -3.441 | 4.997 |
| Ordinal IPTW | 29060_002_PAROXETINE | 29060_002 | 8 | Target 20 mg | 20 | 11 | 11 | 6 | 6 | 13.321 | 2.222 | 19.861 | 4.748 | -6.540 | 5.242 | -16.814 | 3.733 |
| Ordinal IPTW | 29060_002_PAROXETINE | 29060_002 | 1 | Target 30 mg | 30 | 86 | 86 | 73 | 73 | 5.586 | 1.403 | 3.753 | 0.839 | 1.833 | 1.635 | -1.371 | 5.037 |
| Ordinal IPTW | 29060_002_PAROXETINE | 29060_002 | 2 | Target 30 mg | 30 | 91 | 91 | 99 | 99 | 7.872 | 1.428 | 4.322 | 0.827 | 3.551 | 1.650 | 0.317 | 6.785 |
| Ordinal IPTW | 29060_002_PAROXETINE | 29060_002 | 3 | Target 30 mg | 30 | 89 | 89 | 89 | 89 | 9.481 | 1.790 | 7.746 | 0.933 | 1.735 | 2.019 | -2.222 | 5.692 |
| Ordinal IPTW | 29060_002_PAROXETINE | 29060_002 | 4 | Target 30 mg | 30 | 81 | 81 | 86 | 86 | 11.466 | 1.726 | 9.505 | 1.155 | 1.962 | 2.076 | -2.108 | 6.031 |
| Ordinal IPTW | 29060_002_PAROXETINE | 29060_002 | 5 | Target 30 mg | 30 | 50 | 50 | 37 | 37 | 13.120 | 1.839 | 9.757 | 1.674 | 3.363 | 2.487 | -1.511 | 8.236 |
| Ordinal IPTW | 29060_002_PAROXETINE | 29060_002 | 6 | Target 30 mg | 30 | 55 | 55 | 49 | 49 | 14.552 | 2.433 | 10.814 | 1.163 | 3.738 | 2.696 | -1.547 | 9.023 |
| Ordinal IPTW | 29060_002_PAROXETINE | 29060_002 | 7 | Target 30 mg | 30 | 62 | 62 | 45 | 45 | 15.874 | 3.456 | 12.941 | 1.436 | 2.933 | 3.742 | -4.401 | 10.267 |
| Ordinal IPTW | 29060_002_PAROXETINE | 29060_002 | 8 | Target 30 mg | 30 | 11 | 11 | 6 | 6 | 17.178 | 4.657 | 19.861 | 4.748 | -2.684 | 6.650 | -15.718 | 10.351 |
| Ordinal IPTW | 29060_002_PAROXETINE | 29060_002 | 1 | Target 40 mg | 40 | 86 | 86 | 73 | 73 | 3.868 | 1.189 | 3.753 | 0.839 | 0.115 | 1.455 | -2.738 | 2.967 |
| Ordinal IPTW | 29060_002_PAROXETINE | 29060_002 | 2 | Target 40 mg | 40 | 91 | 91 | 99 | 99 | 5.308 | 1.720 | 4.322 | 0.827 | 0.986 | 1.908 | -2.753 | 4.726 |
| Ordinal IPTW | 29060_002_PAROXETINE | 29060_002 | 3 | Target 40 mg | 40 | 89 | 89 | 89 | 89 | 4.746 | 1.938 | 7.746 | 0.933 | -3.000 | 2.151 | -7.215 | 1.216 |
| Ordinal IPTW | 29060_002_PAROXETINE | 29060_002 | 4 | Target 40 mg | 40 | 81 | 81 | 86 | 86 | 6.813 | 1.798 | 9.505 | 1.155 | -2.692 | 2.137 | -6.880 | 1.496 |
| Ordinal IPTW | 29060_002_PAROXETINE | 29060_002 | 5 | Target 40 mg | 40 | 50 | 50 | 37 | 37 | 8.736 | 1.465 | 9.757 | 1.674 | -1.021 | 2.224 | -5.380 | 3.338 |
| Ordinal IPTW | 29060_002_PAROXETINE | 29060_002 | 6 | Target 40 mg | 40 | 55 | 55 | 49 | 49 | 10.563 | 1.200 | 10.814 | 1.163 | -0.251 | 1.671 | -3.526 | 3.025 |
| Ordinal IPTW | 29060_002_PAROXETINE | 29060_002 | 7 | Target 40 mg | 40 | 62 | 62 | 45 | 45 | 12.343 | 1.499 | 12.941 | 1.436 | -0.599 | 2.075 | -4.666 | 3.469 |
| Ordinal IPTW | 29060_002_PAROXETINE | 29060_002 | 8 | Target 40 mg | 40 | 11 | 11 | 6 | 6 | 14.114 | 2.233 | 19.861 | 4.748 | -5.747 | 5.246 | -16.030 | 4.535 |
| Ordinal IPTW | 29060_002_PAROXETINE | 29060_002 | 1 | Target 50 mg | 50 | 86 | 86 | 73 | 73 | 0.879 | 1.427 | 3.753 | 0.839 | -2.874 | 1.655 | -6.118 | 0.370 |
| Ordinal IPTW | 29060_002_PAROXETINE | 29060_002 | 2 | Target 50 mg | 50 | 91 | 91 | 99 | 99 | 6.137 | 1.991 | 4.322 | 0.827 | 1.816 | 2.156 | -2.410 | 6.042 |
| Ordinal IPTW | 29060_002_PAROXETINE | 29060_002 | 3 | Target 50 mg | 50 | 89 | 89 | 89 | 89 | 7.472 | 2.005 | 7.746 | 0.933 | -0.274 | 2.211 | -4.608 | 4.061 |
| Ordinal IPTW | 29060_002_PAROXETINE | 29060_002 | 4 | Target 50 mg | 50 | 81 | 81 | 86 | 86 | 9.828 | 1.933 | 9.505 | 1.155 | 0.323 | 2.252 | -4.090 | 4.737 |
| Ordinal IPTW | 29060_002_PAROXETINE | 29060_002 | 5 | Target 50 mg | 50 | 50 | 50 | 37 | 37 | 10.586 | 1.778 | 9.757 | 1.674 | 0.829 | 2.442 | -3.957 | 5.615 |
| Ordinal IPTW | 29060_002_PAROXETINE | 29060_002 | 6 | Target 50 mg | 50 | 55 | 55 | 49 | 49 | 10.279 | 1.670 | 10.814 | 1.163 | -0.535 | 2.035 | -4.524 | 3.454 |
| Ordinal IPTW | 29060_002_PAROXETINE | 29060_002 | 7 | Target 50 mg | 50 | 62 | 62 | 45 | 45 | 9.440 | 1.845 | 12.941 | 1.436 | -3.501 | 2.338 | -8.084 | 1.081 |
| Ordinal IPTW | 29060_002_PAROXETINE | 29060_002 | 8 | Target 50 mg | 50 | 11 | 11 | 6 | 6 | 8.512 | 2.306 | 19.861 | 4.748 | -11.349 | 5.278 | -21.694 | -1.005 |
| Ordinal IPTW | 29060_003_IMIPRAMINE | 29060_003 | 1 | Target 20 mg | 20 | 120 | 120 | 125 | 125 | 4.922 | 1.969 | 3.113 | 0.645 | 1.808 | 2.072 | -2.252 | 5.869 |
| Ordinal IPTW | 29060_003_IMIPRAMINE | 29060_003 | 2 | Target 20 mg | 20 | 144 | 144 | 144 | 144 | 7.686 | 2.107 | 4.042 | 0.681 | 3.644 | 2.214 | -0.696 | 7.984 |
| Ordinal IPTW | 29060_003_IMIPRAMINE | 29060_003 | 3 | Target 20 mg | 20 | 132 | 132 | 144 | 144 | 10.118 | 1.506 | 5.359 | 0.671 | 4.759 | 1.648 | 1.528 | 7.990 |
| Ordinal IPTW | 29060_003_IMIPRAMINE | 29060_003 | 4 | Target 20 mg | 20 | 108 | 108 | 134 | 134 | 11.980 | 1.443 | 7.574 | 0.824 | 4.405 | 1.662 | 1.149 | 7.662 |
| Ordinal IPTW | 29060_003_IMIPRAMINE | 29060_003 | 5 | Target 20 mg | 20 | 53 | 53 | 58 | 58 | 14.393 | 1.379 | 8.501 | 1.149 | 5.892 | 1.795 | 2.374 | 9.410 |
| Ordinal IPTW | 29060_003_IMIPRAMINE | 29060_003 | 6 | Target 20 mg | 20 | 73 | 73 | 68 | 68 | 17.083 | 1.855 | 10.208 | 1.163 | 6.875 | 2.190 | 2.583 | 11.167 |
| Ordinal IPTW | 29060_003_IMIPRAMINE | 29060_003 | 7 | Target 20 mg | 20 | 50 | 50 | 82 | 82 | 19.819 | 2.730 | 13.486 | 1.137 | 6.333 | 2.958 | 0.536 | 12.130 |
| Ordinal IPTW | 29060_003_IMIPRAMINE | 29060_003 | 8 | Target 20 mg | 20 | 2 | 2 | 1 | 1 | 22.555 | 3.733 | 12.693 | NA | 9.862 | NA | NA | NA |
| Ordinal IPTW | 29060_003_IMIPRAMINE | 29060_003 | 1 | Target 65 mg | 65 | 120 | 120 | 125 | 125 | 6.121 | 2.102 | 3.113 | 0.645 | 3.008 | 2.199 | -1.302 | 7.317 |
| Ordinal IPTW | 29060_003_IMIPRAMINE | 29060_003 | 2 | Target 65 mg | 65 | 144 | 144 | 144 | 144 | 8.566 | 1.803 | 4.042 | 0.681 | 4.524 | 1.927 | 0.747 | 8.301 |
| Ordinal IPTW | 29060_003_IMIPRAMINE | 29060_003 | 3 | Target 65 mg | 65 | 132 | 132 | 144 | 144 | 9.660 | 2.635 | 5.359 | 0.671 | 4.301 | 2.719 | -1.029 | 9.630 |
| Ordinal IPTW | 29060_003_IMIPRAMINE | 29060_003 | 4 | Target 65 mg | 65 | 108 | 108 | 134 | 134 | 11.154 | 2.676 | 7.574 | 0.824 | 3.580 | 2.800 | -1.908 | 9.068 |
| Ordinal IPTW | 29060_003_IMIPRAMINE | 29060_003 | 5 | Target 65 mg | 65 | 53 | 53 | 58 | 58 | 11.242 | 2.451 | 8.501 | 1.149 | 2.741 | 2.707 | -2.565 | 8.046 |
| Ordinal IPTW | 29060_003_IMIPRAMINE | 29060_003 | 6 | Target 65 mg | 65 | 73 | 73 | 68 | 68 | 10.628 | 3.111 | 10.208 | 1.163 | 0.420 | 3.321 | -6.090 | 6.929 |
| Ordinal IPTW | 29060_003_IMIPRAMINE | 29060_003 | 7 | Target 65 mg | 65 | 50 | 50 | 82 | 82 | 9.896 | 4.574 | 13.486 | 1.137 | -3.590 | 4.713 | -12.828 | 5.648 |
| Ordinal IPTW | 29060_003_IMIPRAMINE | 29060_003 | 8 | Target 65 mg | 65 | 2 | 2 | 1 | 1 | 9.165 | 6.314 | 12.693 | NA | -3.528 | NA | NA | NA |
| Ordinal IPTW | 29060_003_IMIPRAMINE | 29060_003 | 1 | Target 80 mg | 80 | 120 | 120 | 125 | 125 | 7.145 | 5.084 | 3.113 | 0.645 | 4.032 | 5.125 | -6.013 | 14.077 |
| Ordinal IPTW | 29060_003_IMIPRAMINE | 29060_003 | 2 | Target 80 mg | 80 | 144 | 144 | 144 | 144 | 9.185 | 2.445 | 4.042 | 0.681 | 5.143 | 2.538 | 0.169 | 10.118 |
| Ordinal IPTW | 29060_003_IMIPRAMINE | 29060_003 | 3 | Target 80 mg | 80 | 132 | 132 | 144 | 144 | 12.242 | 3.107 | 5.359 | 0.671 | 6.883 | 3.179 | 0.653 | 13.114 |
| Ordinal IPTW | 29060_003_IMIPRAMINE | 29060_003 | 4 | Target 80 mg | 80 | 108 | 108 | 134 | 134 | 13.753 | 3.461 | 7.574 | 0.824 | 6.179 | 3.558 | -0.794 | 13.152 |
| Ordinal IPTW | 29060_003_IMIPRAMINE | 29060_003 | 5 | Target 80 mg | 80 | 53 | 53 | 58 | 58 | 15.025 | 3.287 | 8.501 | 1.149 | 6.524 | 3.482 | -0.301 | 13.349 |
| Ordinal IPTW | 29060_003_IMIPRAMINE | 29060_003 | 6 | Target 80 mg | 80 | 73 | 73 | 68 | 68 | 16.178 | 3.753 | 10.208 | 1.163 | 5.970 | 3.929 | -1.732 | 13.671 |
| Ordinal IPTW | 29060_003_IMIPRAMINE | 29060_003 | 7 | Target 80 mg | 80 | 50 | 50 | 82 | 82 | 17.310 | 5.156 | 13.486 | 1.137 | 3.824 | 5.280 | -6.524 | 14.173 |
| Ordinal IPTW | 29060_003_IMIPRAMINE | 29060_003 | 8 | Target 80 mg | 80 | 2 | 2 | 1 | 1 | 18.443 | 6.975 | 12.693 | NA | 5.750 | NA | NA | NA |
| Ordinal IPTW | 29060_003_IMIPRAMINE | 29060_003 | 1 | Target 145 mg | 145 | 120 | 120 | 125 | 125 | 2.548 | 0.943 | 3.113 | 0.645 | -0.565 | 1.142 | -2.804 | 1.674 |
| Ordinal IPTW | 29060_003_IMIPRAMINE | 29060_003 | 2 | Target 145 mg | 145 | 144 | 144 | 144 | 144 | 6.369 | 0.858 | 4.042 | 0.681 | 2.327 | 1.095 | 0.181 | 4.473 |
| Ordinal IPTW | 29060_003_IMIPRAMINE | 29060_003 | 3 | Target 145 mg | 145 | 132 | 132 | 144 | 144 | 12.020 | 1.536 | 5.359 | 0.671 | 6.661 | 1.677 | 3.375 | 9.947 |
| Ordinal IPTW | 29060_003_IMIPRAMINE | 29060_003 | 4 | Target 145 mg | 145 | 108 | 108 | 134 | 134 | 14.301 | 1.387 | 7.574 | 0.824 | 6.727 | 1.613 | 3.566 | 9.889 |
| Ordinal IPTW | 29060_003_IMIPRAMINE | 29060_003 | 5 | Target 145 mg | 145 | 53 | 53 | 58 | 58 | 15.673 | 1.365 | 8.501 | 1.149 | 7.171 | 1.784 | 3.674 | 10.669 |
| Ordinal IPTW | 29060_003_IMIPRAMINE | 29060_003 | 6 | Target 145 mg | 145 | 73 | 73 | 68 | 68 | 16.589 | 1.767 | 10.208 | 1.163 | 6.381 | 2.116 | 2.234 | 10.527 |
| Ordinal IPTW | 29060_003_IMIPRAMINE | 29060_003 | 7 | Target 145 mg | 145 | 50 | 50 | 82 | 82 | 17.429 | 2.465 | 13.486 | 1.137 | 3.943 | 2.715 | -1.378 | 9.263 |
| Ordinal IPTW | 29060_003_IMIPRAMINE | 29060_003 | 8 | Target 145 mg | 145 | 2 | 2 | 1 | 1 | 18.269 | 3.274 | 12.693 | NA | 5.576 | NA | NA | NA |
| Ordinal IPTW | 29060_003_IMIPRAMINE | 29060_003 | 1 | Target 210 mg | 210 | 120 | 120 | 125 | 125 | 4.307 | 0.974 | 3.113 | 0.645 | 1.193 | 1.168 | -1.096 | 3.483 |
| Ordinal IPTW | 29060_003_IMIPRAMINE | 29060_003 | 2 | Target 210 mg | 210 | 144 | 144 | 144 | 144 | 5.879 | 1.287 | 4.042 | 0.681 | 1.837 | 1.456 | -1.016 | 4.690 |
| Ordinal IPTW | 29060_003_IMIPRAMINE | 29060_003 | 3 | Target 210 mg | 210 | 132 | 132 | 144 | 144 | 9.800 | 1.668 | 5.359 | 0.671 | 4.441 | 1.798 | 0.918 | 7.964 |
| Ordinal IPTW | 29060_003_IMIPRAMINE | 29060_003 | 4 | Target 210 mg | 210 | 108 | 108 | 134 | 134 | 11.430 | 1.503 | 7.574 | 0.824 | 3.856 | 1.714 | 0.496 | 7.216 |
| Ordinal IPTW | 29060_003_IMIPRAMINE | 29060_003 | 5 | Target 210 mg | 210 | 53 | 53 | 58 | 58 | 12.509 | 1.328 | 8.501 | 1.149 | 4.008 | 1.756 | 0.565 | 7.450 |
| Ordinal IPTW | 29060_003_IMIPRAMINE | 29060_003 | 6 | Target 210 mg | 210 | 73 | 73 | 68 | 68 | 13.313 | 1.243 | 10.208 | 1.163 | 3.105 | 1.702 | -0.232 | 6.441 |
| Ordinal IPTW | 29060_003_IMIPRAMINE | 29060_003 | 7 | Target 210 mg | 210 | 50 | 50 | 82 | 82 | 14.070 | 1.394 | 13.486 | 1.137 | 0.584 | 1.799 | -2.941 | 4.110 |
| Ordinal IPTW | 29060_003_IMIPRAMINE | 29060_003 | 8 | Target 210 mg | 210 | 2 | 2 | 1 | 1 | 14.828 | 1.733 | 12.693 | NA | 2.135 | NA | NA | NA |
| Ordinal IPTW | 29060_003_IMIPRAMINE | 29060_003 | 1 | Target 275 mg | 275 | 120 | 120 | 125 | 125 | 6.633 | 2.165 | 3.113 | 0.645 | 3.520 | 2.259 | -0.907 | 7.946 |
| Ordinal IPTW | 29060_003_IMIPRAMINE | 29060_003 | 2 | Target 275 mg | 275 | 144 | 144 | 144 | 144 | 7.471 | 2.628 | 4.042 | 0.681 | 3.429 | 2.715 | -1.892 | 8.750 |
| Ordinal IPTW | 29060_003_IMIPRAMINE | 29060_003 | 3 | Target 275 mg | 275 | 132 | 132 | 144 | 144 | 12.307 | 3.726 | 5.359 | 0.671 | 6.948 | 3.786 | -0.472 | 14.368 |
| Ordinal IPTW | 29060_003_IMIPRAMINE | 29060_003 | 4 | Target 275 mg | 275 | 108 | 108 | 134 | 134 | 13.723 | 3.603 | 7.574 | 0.824 | 6.149 | 3.696 | -1.094 | 13.392 |
| Ordinal IPTW | 29060_003_IMIPRAMINE | 29060_003 | 5 | Target 275 mg | 275 | 53 | 53 | 58 | 58 | 14.990 | 3.231 | 8.501 | 1.149 | 6.488 | 3.429 | -0.232 | 13.209 |
| Ordinal IPTW | 29060_003_IMIPRAMINE | 29060_003 | 6 | Target 275 mg | 275 | 73 | 73 | 68 | 68 | 16.181 | 3.200 | 10.208 | 1.163 | 5.973 | 3.404 | -0.700 | 12.646 |
| Ordinal IPTW | 29060_003_IMIPRAMINE | 29060_003 | 7 | Target 275 mg | 275 | 50 | 50 | 82 | 82 | 17.360 | 3.793 | 13.486 | 1.137 | 3.874 | 3.960 | -3.887 | 11.635 |
| Ordinal IPTW | 29060_003_IMIPRAMINE | 29060_003 | 8 | Target 275 mg | 275 | 2 | 2 | 1 | 1 | 18.539 | 4.800 | 12.693 | NA | 5.846 | NA | NA | NA |
| Ordinal IPTW | 29060_003_PAROXETINE | 29060_003 | 1 | Target 10 mg | 10 | 126 | 126 | 125 | 125 | 9.634 | 4.035 | 3.113 | 0.645 | 6.521 | 4.086 | -1.488 | 14.529 |
| Ordinal IPTW | 29060_003_PAROXETINE | 29060_003 | 2 | Target 10 mg | 10 | 153 | 153 | 144 | 144 | 7.487 | 2.897 | 4.042 | 0.681 | 3.445 | 2.976 | -2.387 | 9.277 |
| Ordinal IPTW | 29060_003_PAROXETINE | 29060_003 | 3 | Target 10 mg | 10 | 145 | 145 | 144 | 144 | 9.286 | 1.549 | 5.359 | 0.671 | 3.927 | 1.689 | 0.617 | 7.236 |
| Ordinal IPTW | 29060_003_PAROXETINE | 29060_003 | 4 | Target 10 mg | 10 | 115 | 115 | 134 | 134 | 9.727 | 1.221 | 7.574 | 0.824 | 2.153 | 1.473 | -0.734 | 5.040 |
| Ordinal IPTW | 29060_003_PAROXETINE | 29060_003 | 5 | Target 10 mg | 10 | 58 | 58 | 58 | 58 | 11.589 | 2.488 | 8.501 | 1.149 | 3.088 | 2.740 | -2.283 | 8.459 |
| Ordinal IPTW | 29060_003_PAROXETINE | 29060_003 | 6 | Target 10 mg | 10 | 79 | 79 | 68 | 68 | 14.399 | 4.902 | 10.208 | 1.163 | 4.191 | 5.038 | -5.683 | 14.066 |
| Ordinal IPTW | 29060_003_PAROXETINE | 29060_003 | 7 | Target 10 mg | 10 | 74 | 74 | 82 | 82 | 17.683 | 7.819 | 13.486 | 1.137 | 4.197 | 7.902 | -11.290 | 19.683 |
| Ordinal IPTW | 29060_003_PAROXETINE | 29060_003 | 8 | Target 10 mg | 10 | 6 | 6 | 1 | 1 | 21.045 | 10.849 | 12.693 | NA | 8.353 | NA | NA | NA |
| Ordinal IPTW | 29060_003_PAROXETINE | 29060_003 | 1 | Target 20 mg | 20 | 126 | 126 | 125 | 125 | 3.708 | 0.969 | 3.113 | 0.645 | 0.595 | 1.164 | -1.686 | 2.876 |
| Ordinal IPTW | 29060_003_PAROXETINE | 29060_003 | 2 | Target 20 mg | 20 | 153 | 153 | 144 | 144 | 6.245 | 1.052 | 4.042 | 0.681 | 2.203 | 1.253 | -0.252 | 4.658 |
| Ordinal IPTW | 29060_003_PAROXETINE | 29060_003 | 3 | Target 20 mg | 20 | 145 | 145 | 144 | 144 | 8.550 | 1.218 | 5.359 | 0.671 | 3.191 | 1.391 | 0.465 | 5.917 |
| Ordinal IPTW | 29060_003_PAROXETINE | 29060_003 | 4 | Target 20 mg | 20 | 115 | 115 | 134 | 134 | 10.966 | 1.151 | 7.574 | 0.824 | 3.392 | 1.416 | 0.617 | 6.167 |
| Ordinal IPTW | 29060_003_PAROXETINE | 29060_003 | 5 | Target 20 mg | 20 | 58 | 58 | 58 | 58 | 13.234 | 0.985 | 8.501 | 1.149 | 4.732 | 1.513 | 1.766 | 7.699 |
| Ordinal IPTW | 29060_003_PAROXETINE | 29060_003 | 6 | Target 20 mg | 20 | 79 | 79 | 68 | 68 | 15.402 | 1.002 | 10.208 | 1.163 | 5.194 | 1.535 | 2.185 | 8.204 |
| Ordinal IPTW | 29060_003_PAROXETINE | 29060_003 | 7 | Target 20 mg | 20 | 74 | 74 | 82 | 82 | 17.522 | 1.430 | 13.486 | 1.137 | 4.035 | 1.827 | 0.455 | 7.616 |
| Ordinal IPTW | 29060_003_PAROXETINE | 29060_003 | 8 | Target 20 mg | 20 | 6 | 6 | 1 | 1 | 19.633 | 2.067 | 12.693 | NA | 6.940 | NA | NA | NA |
| Ordinal IPTW | 29060_003_PAROXETINE | 29060_003 | 1 | Target 30 mg | 30 | 126 | 126 | 125 | 125 | 3.853 | 1.006 | 3.113 | 0.645 | 0.740 | 1.195 | -1.602 | 3.082 |
| Ordinal IPTW | 29060_003_PAROXETINE | 29060_003 | 2 | Target 30 mg | 30 | 153 | 153 | 144 | 144 | 6.552 | 0.833 | 4.042 | 0.681 | 2.510 | 1.076 | 0.402 | 4.618 |
| Ordinal IPTW | 29060_003_PAROXETINE | 29060_003 | 3 | Target 30 mg | 30 | 145 | 145 | 144 | 144 | 8.888 | 1.100 | 5.359 | 0.671 | 3.529 | 1.289 | 1.004 | 6.055 |
| Ordinal IPTW | 29060_003_PAROXETINE | 29060_003 | 4 | Target 30 mg | 30 | 115 | 115 | 134 | 134 | 10.568 | 1.063 | 7.574 | 0.824 | 2.994 | 1.345 | 0.359 | 5.630 |
| Ordinal IPTW | 29060_003_PAROXETINE | 29060_003 | 5 | Target 30 mg | 30 | 58 | 58 | 58 | 58 | 11.872 | 0.951 | 8.501 | 1.149 | 3.371 | 1.491 | 0.448 | 6.294 |
| Ordinal IPTW | 29060_003_PAROXETINE | 29060_003 | 6 | Target 30 mg | 30 | 79 | 79 | 68 | 68 | 12.926 | 0.911 | 10.208 | 1.163 | 2.717 | 1.478 | -0.179 | 5.614 |
| Ordinal IPTW | 29060_003_PAROXETINE | 29060_003 | 7 | Target 30 mg | 30 | 74 | 74 | 82 | 82 | 13.853 | 1.167 | 13.486 | 1.137 | 0.367 | 1.630 | -2.827 | 3.561 |
| Ordinal IPTW | 29060_003_PAROXETINE | 29060_003 | 8 | Target 30 mg | 30 | 6 | 6 | 1 | 1 | 14.760 | 1.630 | 12.693 | NA | 2.068 | NA | NA | NA |
| Ordinal IPTW | 29060_003_PAROXETINE | 29060_003 | 1 | Target 40 mg | 40 | 126 | 126 | 125 | 125 | 5.088 | 1.805 | 3.113 | 0.645 | 1.975 | 1.916 | -1.781 | 5.731 |
| Ordinal IPTW | 29060_003_PAROXETINE | 29060_003 | 2 | Target 40 mg | 40 | 153 | 153 | 144 | 144 | 8.559 | 1.311 | 4.042 | 0.681 | 4.517 | 1.477 | 1.622 | 7.411 |
| Ordinal IPTW | 29060_003_PAROXETINE | 29060_003 | 3 | Target 40 mg | 40 | 145 | 145 | 144 | 144 | 10.862 | 1.643 | 5.359 | 0.671 | 5.503 | 1.775 | 2.024 | 8.982 |
| Ordinal IPTW | 29060_003_PAROXETINE | 29060_003 | 4 | Target 40 mg | 40 | 115 | 115 | 134 | 134 | 12.207 | 1.564 | 7.574 | 0.824 | 4.633 | 1.768 | 1.168 | 8.099 |
| Ordinal IPTW | 29060_003_PAROXETINE | 29060_003 | 5 | Target 40 mg | 40 | 58 | 58 | 58 | 58 | 13.440 | 1.317 | 8.501 | 1.149 | 4.939 | 1.748 | 1.514 | 8.364 |
| Ordinal IPTW | 29060_003_PAROXETINE | 29060_003 | 6 | Target 40 mg | 40 | 79 | 79 | 68 | 68 | 14.598 | 1.137 | 10.208 | 1.163 | 4.390 | 1.627 | 1.202 | 7.579 |
| Ordinal IPTW | 29060_003_PAROXETINE | 29060_003 | 7 | Target 40 mg | 40 | 74 | 74 | 82 | 82 | 15.719 | 1.457 | 13.486 | 1.137 | 2.233 | 1.848 | -1.389 | 5.855 |
| Ordinal IPTW | 29060_003_PAROXETINE | 29060_003 | 8 | Target 40 mg | 40 | 6 | 6 | 1 | 1 | 16.833 | 2.139 | 12.693 | NA | 4.140 | NA | NA | NA |
| Ordinal IPTW | 29060_003_PAROXETINE | 29060_003 | 1 | Target 50 mg | 50 | 126 | 126 | 125 | 125 | 1.417 | 2.253 | 3.113 | 0.645 | -1.696 | 2.344 | -6.290 | 2.898 |
| Ordinal IPTW | 29060_003_PAROXETINE | 29060_003 | 2 | Target 50 mg | 50 | 153 | 153 | 144 | 144 | 7.583 | 1.249 | 4.042 | 0.681 | 3.541 | 1.422 | 0.753 | 6.329 |
| Ordinal IPTW | 29060_003_PAROXETINE | 29060_003 | 3 | Target 50 mg | 50 | 145 | 145 | 144 | 144 | 10.464 | 1.979 | 5.359 | 0.671 | 5.105 | 2.090 | 1.009 | 9.201 |
| Ordinal IPTW | 29060_003_PAROXETINE | 29060_003 | 4 | Target 50 mg | 50 | 115 | 115 | 134 | 134 | 12.656 | 1.806 | 7.574 | 0.824 | 5.082 | 1.985 | 1.191 | 8.973 |
| Ordinal IPTW | 29060_003_PAROXETINE | 29060_003 | 5 | Target 50 mg | 50 | 58 | 58 | 58 | 58 | 13.695 | 1.620 | 8.501 | 1.149 | 5.194 | 1.986 | 1.302 | 9.086 |
| Ordinal IPTW | 29060_003_PAROXETINE | 29060_003 | 6 | Target 50 mg | 50 | 79 | 79 | 68 | 68 | 13.965 | 1.584 | 10.208 | 1.163 | 3.757 | 1.965 | -0.096 | 7.609 |
| Ordinal IPTW | 29060_003_PAROXETINE | 29060_003 | 7 | Target 50 mg | 50 | 74 | 74 | 82 | 82 | 13.850 | 1.950 | 13.486 | 1.137 | 0.363 | 2.257 | -4.061 | 4.788 |
| Ordinal IPTW | 29060_003_PAROXETINE | 29060_003 | 8 | Target 50 mg | 50 | 6 | 6 | 1 | 1 | 13.671 | 2.603 | 12.693 | NA | 0.978 | NA | NA | NA |

<small><em>Dose-strategy predictions versus observed placebo by
week</em></small>

| weight_model | trial_name | arm_name | week | source | dose_strategy | strategy_dose | mean_improvement | SE | lower_95 | upper_95 | n_patients | n_visit_rows |
|:---|:---|:---|---:|:---|:---|---:|---:|---:|---:|---:|---:|---:|
| Ordinal IPTW | 29060_002 | 29060_002_PAROXETINE | 1 | Observed placebo | Observed placebo | NA | 3.753 | 0.839 | 2.108 | 5.398 | 73 | 73 |
| Ordinal IPTW | 29060_002 | 29060_002_PAROXETINE | 1 | Active GEE/MSM prediction | Target 10 mg | 10 | 2.091 | 2.791 | -3.380 | 7.562 | 86 | 86 |
| Ordinal IPTW | 29060_002 | 29060_002_PAROXETINE | 1 | Active GEE/MSM prediction | Target 20 mg | 20 | 1.365 | 1.561 | -1.695 | 4.424 | 86 | 86 |
| Ordinal IPTW | 29060_002 | 29060_002_PAROXETINE | 1 | Active GEE/MSM prediction | Target 30 mg | 30 | 5.586 | 1.403 | 2.837 | 8.335 | 86 | 86 |
| Ordinal IPTW | 29060_002 | 29060_002_PAROXETINE | 1 | Active GEE/MSM prediction | Target 40 mg | 40 | 3.868 | 1.189 | 1.537 | 6.198 | 86 | 86 |
| Ordinal IPTW | 29060_002 | 29060_002_PAROXETINE | 1 | Active GEE/MSM prediction | Target 50 mg | 50 | 0.879 | 1.427 | -1.918 | 3.675 | 86 | 86 |
| Ordinal IPTW | 29060_002 | 29060_002_PAROXETINE | 2 | Observed placebo | Observed placebo | NA | 4.322 | 0.827 | 2.701 | 5.942 | 99 | 99 |
| Ordinal IPTW | 29060_002 | 29060_002_PAROXETINE | 2 | Active GEE/MSM prediction | Target 10 mg | 10 | 5.231 | 2.126 | 1.065 | 9.397 | 91 | 91 |
| Ordinal IPTW | 29060_002 | 29060_002_PAROXETINE | 2 | Active GEE/MSM prediction | Target 20 mg | 20 | 5.151 | 1.409 | 2.390 | 7.912 | 91 | 91 |
| Ordinal IPTW | 29060_002 | 29060_002_PAROXETINE | 2 | Active GEE/MSM prediction | Target 30 mg | 30 | 7.872 | 1.428 | 5.073 | 10.671 | 91 | 91 |
| Ordinal IPTW | 29060_002 | 29060_002_PAROXETINE | 2 | Active GEE/MSM prediction | Target 40 mg | 40 | 5.308 | 1.720 | 1.938 | 8.678 | 91 | 91 |
| Ordinal IPTW | 29060_002 | 29060_002_PAROXETINE | 2 | Active GEE/MSM prediction | Target 50 mg | 50 | 6.137 | 1.991 | 2.234 | 10.041 | 91 | 91 |
| Ordinal IPTW | 29060_002 | 29060_002_PAROXETINE | 3 | Observed placebo | Observed placebo | NA | 7.746 | 0.933 | 5.918 | 9.574 | 89 | 89 |
| Ordinal IPTW | 29060_002 | 29060_002_PAROXETINE | 3 | Active GEE/MSM prediction | Target 10 mg | 10 | 6.881 | 3.192 | 0.624 | 13.138 | 89 | 89 |
| Ordinal IPTW | 29060_002 | 29060_002_PAROXETINE | 3 | Active GEE/MSM prediction | Target 20 mg | 20 | 9.388 | 1.649 | 6.156 | 12.621 | 89 | 89 |
| Ordinal IPTW | 29060_002 | 29060_002_PAROXETINE | 3 | Active GEE/MSM prediction | Target 30 mg | 30 | 9.481 | 1.790 | 5.972 | 12.990 | 89 | 89 |
| Ordinal IPTW | 29060_002 | 29060_002_PAROXETINE | 3 | Active GEE/MSM prediction | Target 40 mg | 40 | 4.746 | 1.938 | 0.948 | 8.545 | 89 | 89 |
| Ordinal IPTW | 29060_002 | 29060_002_PAROXETINE | 3 | Active GEE/MSM prediction | Target 50 mg | 50 | 7.472 | 2.005 | 3.542 | 11.402 | 89 | 89 |
| Ordinal IPTW | 29060_002 | 29060_002_PAROXETINE | 4 | Observed placebo | Observed placebo | NA | 9.505 | 1.155 | 7.242 | 11.768 | 86 | 86 |
| Ordinal IPTW | 29060_002 | 29060_002_PAROXETINE | 4 | Active GEE/MSM prediction | Target 10 mg | 10 | 10.936 | 3.068 | 4.923 | 16.949 | 81 | 81 |
| Ordinal IPTW | 29060_002 | 29060_002_PAROXETINE | 4 | Active GEE/MSM prediction | Target 20 mg | 20 | 12.415 | 1.506 | 9.462 | 15.367 | 81 | 81 |
| Ordinal IPTW | 29060_002 | 29060_002_PAROXETINE | 4 | Active GEE/MSM prediction | Target 30 mg | 30 | 11.466 | 1.726 | 8.084 | 14.849 | 81 | 81 |
| Ordinal IPTW | 29060_002 | 29060_002_PAROXETINE | 4 | Active GEE/MSM prediction | Target 40 mg | 40 | 6.813 | 1.798 | 3.289 | 10.337 | 81 | 81 |
| Ordinal IPTW | 29060_002 | 29060_002_PAROXETINE | 4 | Active GEE/MSM prediction | Target 50 mg | 50 | 9.828 | 1.933 | 6.039 | 13.617 | 81 | 81 |
| Ordinal IPTW | 29060_002 | 29060_002_PAROXETINE | 5 | Observed placebo | Observed placebo | NA | 9.757 | 1.674 | 6.477 | 13.037 | 37 | 37 |
| Ordinal IPTW | 29060_002 | 29060_002_PAROXETINE | 5 | Active GEE/MSM prediction | Target 10 mg | 10 | 12.562 | 2.622 | 7.423 | 17.702 | 50 | 50 |
| Ordinal IPTW | 29060_002 | 29060_002_PAROXETINE | 5 | Active GEE/MSM prediction | Target 20 mg | 20 | 13.775 | 1.318 | 11.191 | 16.359 | 50 | 50 |
| Ordinal IPTW | 29060_002 | 29060_002_PAROXETINE | 5 | Active GEE/MSM prediction | Target 30 mg | 30 | 13.120 | 1.839 | 9.515 | 16.725 | 50 | 50 |
| Ordinal IPTW | 29060_002 | 29060_002_PAROXETINE | 5 | Active GEE/MSM prediction | Target 40 mg | 40 | 8.736 | 1.465 | 5.865 | 11.607 | 50 | 50 |
| Ordinal IPTW | 29060_002 | 29060_002_PAROXETINE | 5 | Active GEE/MSM prediction | Target 50 mg | 50 | 10.586 | 1.778 | 7.101 | 14.071 | 50 | 50 |
| Ordinal IPTW | 29060_002 | 29060_002_PAROXETINE | 6 | Observed placebo | Observed placebo | NA | 10.814 | 1.163 | 8.534 | 13.094 | 49 | 49 |
| Ordinal IPTW | 29060_002 | 29060_002_PAROXETINE | 6 | Active GEE/MSM prediction | Target 10 mg | 10 | 12.569 | 2.248 | 8.163 | 16.974 | 55 | 55 |
| Ordinal IPTW | 29060_002 | 29060_002_PAROXETINE | 6 | Active GEE/MSM prediction | Target 20 mg | 20 | 14.025 | 1.264 | 11.548 | 16.502 | 55 | 55 |
| Ordinal IPTW | 29060_002 | 29060_002_PAROXETINE | 6 | Active GEE/MSM prediction | Target 30 mg | 30 | 14.552 | 2.433 | 9.785 | 19.320 | 55 | 55 |
| Ordinal IPTW | 29060_002 | 29060_002_PAROXETINE | 6 | Active GEE/MSM prediction | Target 40 mg | 40 | 10.563 | 1.200 | 8.211 | 12.916 | 55 | 55 |
| Ordinal IPTW | 29060_002 | 29060_002_PAROXETINE | 6 | Active GEE/MSM prediction | Target 50 mg | 50 | 10.279 | 1.670 | 7.005 | 13.553 | 55 | 55 |
| Ordinal IPTW | 29060_002 | 29060_002_PAROXETINE | 7 | Observed placebo | Observed placebo | NA | 12.941 | 1.436 | 10.127 | 15.755 | 45 | 45 |
| Ordinal IPTW | 29060_002 | 29060_002_PAROXETINE | 7 | Active GEE/MSM prediction | Target 10 mg | 10 | 11.765 | 2.638 | 6.595 | 16.936 | 62 | 62 |
| Ordinal IPTW | 29060_002 | 29060_002_PAROXETINE | 7 | Active GEE/MSM prediction | Target 20 mg | 20 | 13.719 | 1.604 | 10.576 | 16.863 | 62 | 62 |
| Ordinal IPTW | 29060_002 | 29060_002_PAROXETINE | 7 | Active GEE/MSM prediction | Target 30 mg | 30 | 15.874 | 3.456 | 9.101 | 22.647 | 62 | 62 |
| Ordinal IPTW | 29060_002 | 29060_002_PAROXETINE | 7 | Active GEE/MSM prediction | Target 40 mg | 40 | 12.343 | 1.499 | 9.406 | 15.280 | 62 | 62 |
| Ordinal IPTW | 29060_002 | 29060_002_PAROXETINE | 7 | Active GEE/MSM prediction | Target 50 mg | 50 | 9.440 | 1.845 | 5.824 | 13.056 | 62 | 62 |
| Ordinal IPTW | 29060_002 | 29060_002_PAROXETINE | 8 | Observed placebo | Observed placebo | NA | 19.861 | 4.748 | 10.556 | 29.166 | 6 | 6 |
| Ordinal IPTW | 29060_002 | 29060_002_PAROXETINE | 8 | Active GEE/MSM prediction | Target 10 mg | 10 | 10.827 | 3.697 | 3.582 | 18.072 | 11 | 11 |
| Ordinal IPTW | 29060_002 | 29060_002_PAROXETINE | 8 | Active GEE/MSM prediction | Target 20 mg | 20 | 13.321 | 2.222 | 8.967 | 17.675 | 11 | 11 |
| Ordinal IPTW | 29060_002 | 29060_002_PAROXETINE | 8 | Active GEE/MSM prediction | Target 30 mg | 30 | 17.178 | 4.657 | 8.050 | 26.306 | 11 | 11 |
| Ordinal IPTW | 29060_002 | 29060_002_PAROXETINE | 8 | Active GEE/MSM prediction | Target 40 mg | 40 | 14.114 | 2.233 | 9.738 | 18.490 | 11 | 11 |
| Ordinal IPTW | 29060_002 | 29060_002_PAROXETINE | 8 | Active GEE/MSM prediction | Target 50 mg | 50 | 8.512 | 2.306 | 3.993 | 13.032 | 11 | 11 |
| Ordinal IPTW | 29060_003 | 29060_003_IMIPRAMINE | 1 | Observed placebo | Observed placebo | NA | 3.113 | 0.645 | 1.849 | 4.377 | 125 | 125 |
| Ordinal IPTW | 29060_003 | 29060_003_IMIPRAMINE | 1 | Active GEE/MSM prediction | Target 20 mg | 20 | 4.922 | 1.969 | 1.063 | 8.780 | 120 | 120 |
| Ordinal IPTW | 29060_003 | 29060_003_IMIPRAMINE | 1 | Active GEE/MSM prediction | Target 65 mg | 65 | 6.121 | 2.102 | 2.001 | 10.240 | 120 | 120 |
| Ordinal IPTW | 29060_003 | 29060_003_IMIPRAMINE | 1 | Active GEE/MSM prediction | Target 80 mg | 80 | 7.145 | 5.084 | -2.820 | 17.111 | 120 | 120 |
| Ordinal IPTW | 29060_003 | 29060_003_IMIPRAMINE | 1 | Active GEE/MSM prediction | Target 145 mg | 145 | 2.548 | 0.943 | 0.700 | 4.396 | 120 | 120 |
| Ordinal IPTW | 29060_003 | 29060_003_IMIPRAMINE | 1 | Active GEE/MSM prediction | Target 210 mg | 210 | 4.307 | 0.974 | 2.398 | 6.215 | 120 | 120 |
| Ordinal IPTW | 29060_003 | 29060_003_IMIPRAMINE | 1 | Active GEE/MSM prediction | Target 275 mg | 275 | 6.633 | 2.165 | 2.390 | 10.875 | 120 | 120 |
| Ordinal IPTW | 29060_003 | 29060_003_IMIPRAMINE | 2 | Observed placebo | Observed placebo | NA | 4.042 | 0.681 | 2.708 | 5.376 | 144 | 144 |
| Ordinal IPTW | 29060_003 | 29060_003_IMIPRAMINE | 2 | Active GEE/MSM prediction | Target 20 mg | 20 | 7.686 | 2.107 | 3.556 | 11.816 | 144 | 144 |
| Ordinal IPTW | 29060_003 | 29060_003_IMIPRAMINE | 2 | Active GEE/MSM prediction | Target 65 mg | 65 | 8.566 | 1.803 | 5.032 | 12.099 | 144 | 144 |
| Ordinal IPTW | 29060_003 | 29060_003_IMIPRAMINE | 2 | Active GEE/MSM prediction | Target 80 mg | 80 | 9.185 | 2.445 | 4.393 | 13.978 | 144 | 144 |
| Ordinal IPTW | 29060_003 | 29060_003_IMIPRAMINE | 2 | Active GEE/MSM prediction | Target 145 mg | 145 | 6.369 | 0.858 | 4.688 | 8.051 | 144 | 144 |
| Ordinal IPTW | 29060_003 | 29060_003_IMIPRAMINE | 2 | Active GEE/MSM prediction | Target 210 mg | 210 | 5.879 | 1.287 | 3.358 | 8.401 | 144 | 144 |
| Ordinal IPTW | 29060_003 | 29060_003_IMIPRAMINE | 2 | Active GEE/MSM prediction | Target 275 mg | 275 | 7.471 | 2.628 | 2.320 | 12.622 | 144 | 144 |
| Ordinal IPTW | 29060_003 | 29060_003_IMIPRAMINE | 3 | Observed placebo | Observed placebo | NA | 5.359 | 0.671 | 4.043 | 6.675 | 144 | 144 |
| Ordinal IPTW | 29060_003 | 29060_003_IMIPRAMINE | 3 | Active GEE/MSM prediction | Target 20 mg | 20 | 10.118 | 1.506 | 7.167 | 13.069 | 132 | 132 |
| Ordinal IPTW | 29060_003 | 29060_003_IMIPRAMINE | 3 | Active GEE/MSM prediction | Target 65 mg | 65 | 9.660 | 2.635 | 4.495 | 14.824 | 132 | 132 |
| Ordinal IPTW | 29060_003 | 29060_003_IMIPRAMINE | 3 | Active GEE/MSM prediction | Target 80 mg | 80 | 12.242 | 3.107 | 6.152 | 18.332 | 132 | 132 |
| Ordinal IPTW | 29060_003 | 29060_003_IMIPRAMINE | 3 | Active GEE/MSM prediction | Target 145 mg | 145 | 12.020 | 1.536 | 9.009 | 15.031 | 132 | 132 |
| Ordinal IPTW | 29060_003 | 29060_003_IMIPRAMINE | 3 | Active GEE/MSM prediction | Target 210 mg | 210 | 9.800 | 1.668 | 6.532 | 13.068 | 132 | 132 |
| Ordinal IPTW | 29060_003 | 29060_003_IMIPRAMINE | 3 | Active GEE/MSM prediction | Target 275 mg | 275 | 12.307 | 3.726 | 5.005 | 19.609 | 132 | 132 |
| Ordinal IPTW | 29060_003 | 29060_003_IMIPRAMINE | 4 | Observed placebo | Observed placebo | NA | 7.574 | 0.824 | 5.960 | 9.189 | 134 | 134 |
| Ordinal IPTW | 29060_003 | 29060_003_IMIPRAMINE | 4 | Active GEE/MSM prediction | Target 20 mg | 20 | 11.980 | 1.443 | 9.151 | 14.808 | 108 | 108 |
| Ordinal IPTW | 29060_003 | 29060_003_IMIPRAMINE | 4 | Active GEE/MSM prediction | Target 65 mg | 65 | 11.154 | 2.676 | 5.909 | 16.399 | 108 | 108 |
| Ordinal IPTW | 29060_003 | 29060_003_IMIPRAMINE | 4 | Active GEE/MSM prediction | Target 80 mg | 80 | 13.753 | 3.461 | 6.969 | 20.537 | 108 | 108 |
| Ordinal IPTW | 29060_003 | 29060_003_IMIPRAMINE | 4 | Active GEE/MSM prediction | Target 145 mg | 145 | 14.301 | 1.387 | 11.583 | 17.020 | 108 | 108 |
| Ordinal IPTW | 29060_003 | 29060_003_IMIPRAMINE | 4 | Active GEE/MSM prediction | Target 210 mg | 210 | 11.430 | 1.503 | 8.484 | 14.376 | 108 | 108 |
| Ordinal IPTW | 29060_003 | 29060_003_IMIPRAMINE | 4 | Active GEE/MSM prediction | Target 275 mg | 275 | 13.723 | 3.603 | 6.662 | 20.784 | 108 | 108 |
| Ordinal IPTW | 29060_003 | 29060_003_IMIPRAMINE | 5 | Observed placebo | Observed placebo | NA | 8.501 | 1.149 | 6.250 | 10.753 | 58 | 58 |
| Ordinal IPTW | 29060_003 | 29060_003_IMIPRAMINE | 5 | Active GEE/MSM prediction | Target 20 mg | 20 | 14.393 | 1.379 | 11.691 | 17.096 | 53 | 53 |
| Ordinal IPTW | 29060_003 | 29060_003_IMIPRAMINE | 5 | Active GEE/MSM prediction | Target 65 mg | 65 | 11.242 | 2.451 | 6.438 | 16.046 | 53 | 53 |
| Ordinal IPTW | 29060_003 | 29060_003_IMIPRAMINE | 5 | Active GEE/MSM prediction | Target 80 mg | 80 | 15.025 | 3.287 | 8.582 | 21.468 | 53 | 53 |
| Ordinal IPTW | 29060_003 | 29060_003_IMIPRAMINE | 5 | Active GEE/MSM prediction | Target 145 mg | 145 | 15.673 | 1.365 | 12.997 | 18.349 | 53 | 53 |
| Ordinal IPTW | 29060_003 | 29060_003_IMIPRAMINE | 5 | Active GEE/MSM prediction | Target 210 mg | 210 | 12.509 | 1.328 | 9.905 | 15.113 | 53 | 53 |
| Ordinal IPTW | 29060_003 | 29060_003_IMIPRAMINE | 5 | Active GEE/MSM prediction | Target 275 mg | 275 | 14.990 | 3.231 | 8.657 | 21.322 | 53 | 53 |
| Ordinal IPTW | 29060_003 | 29060_003_IMIPRAMINE | 6 | Observed placebo | Observed placebo | NA | 10.208 | 1.163 | 7.928 | 12.488 | 68 | 68 |
| Ordinal IPTW | 29060_003 | 29060_003_IMIPRAMINE | 6 | Active GEE/MSM prediction | Target 20 mg | 20 | 17.083 | 1.855 | 13.447 | 20.720 | 73 | 73 |
| Ordinal IPTW | 29060_003 | 29060_003_IMIPRAMINE | 6 | Active GEE/MSM prediction | Target 65 mg | 65 | 10.628 | 3.111 | 4.530 | 16.725 | 73 | 73 |
| Ordinal IPTW | 29060_003 | 29060_003_IMIPRAMINE | 6 | Active GEE/MSM prediction | Target 80 mg | 80 | 16.178 | 3.753 | 8.821 | 23.534 | 73 | 73 |
| Ordinal IPTW | 29060_003 | 29060_003_IMIPRAMINE | 6 | Active GEE/MSM prediction | Target 145 mg | 145 | 16.589 | 1.767 | 13.125 | 20.052 | 73 | 73 |
| Ordinal IPTW | 29060_003 | 29060_003_IMIPRAMINE | 6 | Active GEE/MSM prediction | Target 210 mg | 210 | 13.313 | 1.243 | 10.876 | 15.749 | 73 | 73 |
| Ordinal IPTW | 29060_003 | 29060_003_IMIPRAMINE | 6 | Active GEE/MSM prediction | Target 275 mg | 275 | 16.181 | 3.200 | 9.910 | 22.452 | 73 | 73 |
| Ordinal IPTW | 29060_003 | 29060_003_IMIPRAMINE | 7 | Observed placebo | Observed placebo | NA | 13.486 | 1.137 | 11.257 | 15.715 | 82 | 82 |
| Ordinal IPTW | 29060_003 | 29060_003_IMIPRAMINE | 7 | Active GEE/MSM prediction | Target 20 mg | 20 | 19.819 | 2.730 | 14.468 | 25.171 | 50 | 50 |
| Ordinal IPTW | 29060_003 | 29060_003_IMIPRAMINE | 7 | Active GEE/MSM prediction | Target 65 mg | 65 | 9.896 | 4.574 | 0.931 | 18.861 | 50 | 50 |
| Ordinal IPTW | 29060_003 | 29060_003_IMIPRAMINE | 7 | Active GEE/MSM prediction | Target 80 mg | 80 | 17.310 | 5.156 | 7.205 | 27.416 | 50 | 50 |
| Ordinal IPTW | 29060_003 | 29060_003_IMIPRAMINE | 7 | Active GEE/MSM prediction | Target 145 mg | 145 | 17.429 | 2.465 | 12.598 | 22.260 | 50 | 50 |
| Ordinal IPTW | 29060_003 | 29060_003_IMIPRAMINE | 7 | Active GEE/MSM prediction | Target 210 mg | 210 | 14.070 | 1.394 | 11.339 | 16.802 | 50 | 50 |
| Ordinal IPTW | 29060_003 | 29060_003_IMIPRAMINE | 7 | Active GEE/MSM prediction | Target 275 mg | 275 | 17.360 | 3.793 | 9.926 | 24.794 | 50 | 50 |
| Ordinal IPTW | 29060_003 | 29060_003_IMIPRAMINE | 8 | Observed placebo | Observed placebo | NA | 12.693 | NA | NA | NA | 1 | 1 |
| Ordinal IPTW | 29060_003 | 29060_003_IMIPRAMINE | 8 | Active GEE/MSM prediction | Target 20 mg | 20 | 22.555 | 3.733 | 15.239 | 29.872 | 2 | 2 |
| Ordinal IPTW | 29060_003 | 29060_003_IMIPRAMINE | 8 | Active GEE/MSM prediction | Target 65 mg | 65 | 9.165 | 6.314 | -3.210 | 21.539 | 2 | 2 |
| Ordinal IPTW | 29060_003 | 29060_003_IMIPRAMINE | 8 | Active GEE/MSM prediction | Target 80 mg | 80 | 18.443 | 6.975 | 4.772 | 32.114 | 2 | 2 |
| Ordinal IPTW | 29060_003 | 29060_003_IMIPRAMINE | 8 | Active GEE/MSM prediction | Target 145 mg | 145 | 18.269 | 3.274 | 11.852 | 24.686 | 2 | 2 |
| Ordinal IPTW | 29060_003 | 29060_003_IMIPRAMINE | 8 | Active GEE/MSM prediction | Target 210 mg | 210 | 14.828 | 1.733 | 11.430 | 18.225 | 2 | 2 |
| Ordinal IPTW | 29060_003 | 29060_003_IMIPRAMINE | 8 | Active GEE/MSM prediction | Target 275 mg | 275 | 18.539 | 4.800 | 9.130 | 27.948 | 2 | 2 |
| Ordinal IPTW | 29060_003 | 29060_003_PAROXETINE | 1 | Observed placebo | Observed placebo | NA | 3.113 | 0.645 | 1.849 | 4.377 | 125 | 125 |
| Ordinal IPTW | 29060_003 | 29060_003_PAROXETINE | 1 | Active GEE/MSM prediction | Target 10 mg | 10 | 9.634 | 4.035 | 1.726 | 17.542 | 126 | 126 |
| Ordinal IPTW | 29060_003 | 29060_003_PAROXETINE | 1 | Active GEE/MSM prediction | Target 20 mg | 20 | 3.708 | 0.969 | 1.810 | 5.607 | 126 | 126 |
| Ordinal IPTW | 29060_003 | 29060_003_PAROXETINE | 1 | Active GEE/MSM prediction | Target 30 mg | 30 | 3.853 | 1.006 | 1.882 | 5.824 | 126 | 126 |
| Ordinal IPTW | 29060_003 | 29060_003_PAROXETINE | 1 | Active GEE/MSM prediction | Target 40 mg | 40 | 5.088 | 1.805 | 1.551 | 8.625 | 126 | 126 |
| Ordinal IPTW | 29060_003 | 29060_003_PAROXETINE | 1 | Active GEE/MSM prediction | Target 50 mg | 50 | 1.417 | 2.253 | -3.000 | 5.834 | 126 | 126 |
| Ordinal IPTW | 29060_003 | 29060_003_PAROXETINE | 2 | Observed placebo | Observed placebo | NA | 4.042 | 0.681 | 2.708 | 5.376 | 144 | 144 |
| Ordinal IPTW | 29060_003 | 29060_003_PAROXETINE | 2 | Active GEE/MSM prediction | Target 10 mg | 10 | 7.487 | 2.897 | 1.810 | 13.165 | 153 | 153 |
| Ordinal IPTW | 29060_003 | 29060_003_PAROXETINE | 2 | Active GEE/MSM prediction | Target 20 mg | 20 | 6.245 | 1.052 | 4.184 | 8.306 | 153 | 153 |
| Ordinal IPTW | 29060_003 | 29060_003_PAROXETINE | 2 | Active GEE/MSM prediction | Target 30 mg | 30 | 6.552 | 0.833 | 4.920 | 8.185 | 153 | 153 |
| Ordinal IPTW | 29060_003 | 29060_003_PAROXETINE | 2 | Active GEE/MSM prediction | Target 40 mg | 40 | 8.559 | 1.311 | 5.990 | 11.128 | 153 | 153 |
| Ordinal IPTW | 29060_003 | 29060_003_PAROXETINE | 2 | Active GEE/MSM prediction | Target 50 mg | 50 | 7.583 | 1.249 | 5.135 | 10.032 | 153 | 153 |
| Ordinal IPTW | 29060_003 | 29060_003_PAROXETINE | 3 | Observed placebo | Observed placebo | NA | 5.359 | 0.671 | 4.043 | 6.675 | 144 | 144 |
| Ordinal IPTW | 29060_003 | 29060_003_PAROXETINE | 3 | Active GEE/MSM prediction | Target 10 mg | 10 | 9.286 | 1.549 | 6.249 | 12.322 | 145 | 145 |
| Ordinal IPTW | 29060_003 | 29060_003_PAROXETINE | 3 | Active GEE/MSM prediction | Target 20 mg | 20 | 8.550 | 1.218 | 6.163 | 10.938 | 145 | 145 |
| Ordinal IPTW | 29060_003 | 29060_003_PAROXETINE | 3 | Active GEE/MSM prediction | Target 30 mg | 30 | 8.888 | 1.100 | 6.732 | 11.044 | 145 | 145 |
| Ordinal IPTW | 29060_003 | 29060_003_PAROXETINE | 3 | Active GEE/MSM prediction | Target 40 mg | 40 | 10.862 | 1.643 | 7.641 | 14.082 | 145 | 145 |
| Ordinal IPTW | 29060_003 | 29060_003_PAROXETINE | 3 | Active GEE/MSM prediction | Target 50 mg | 50 | 10.464 | 1.979 | 6.585 | 14.343 | 145 | 145 |
| Ordinal IPTW | 29060_003 | 29060_003_PAROXETINE | 4 | Observed placebo | Observed placebo | NA | 7.574 | 0.824 | 5.960 | 9.189 | 134 | 134 |
| Ordinal IPTW | 29060_003 | 29060_003_PAROXETINE | 4 | Active GEE/MSM prediction | Target 10 mg | 10 | 9.727 | 1.221 | 7.334 | 12.120 | 115 | 115 |
| Ordinal IPTW | 29060_003 | 29060_003_PAROXETINE | 4 | Active GEE/MSM prediction | Target 20 mg | 20 | 10.966 | 1.151 | 8.709 | 13.223 | 115 | 115 |
| Ordinal IPTW | 29060_003 | 29060_003_PAROXETINE | 4 | Active GEE/MSM prediction | Target 30 mg | 30 | 10.568 | 1.063 | 8.486 | 12.651 | 115 | 115 |
| Ordinal IPTW | 29060_003 | 29060_003_PAROXETINE | 4 | Active GEE/MSM prediction | Target 40 mg | 40 | 12.207 | 1.564 | 9.141 | 15.274 | 115 | 115 |
| Ordinal IPTW | 29060_003 | 29060_003_PAROXETINE | 4 | Active GEE/MSM prediction | Target 50 mg | 50 | 12.656 | 1.806 | 9.116 | 16.197 | 115 | 115 |
| Ordinal IPTW | 29060_003 | 29060_003_PAROXETINE | 5 | Observed placebo | Observed placebo | NA | 8.501 | 1.149 | 6.250 | 10.753 | 58 | 58 |
| Ordinal IPTW | 29060_003 | 29060_003_PAROXETINE | 5 | Active GEE/MSM prediction | Target 10 mg | 10 | 11.589 | 2.488 | 6.713 | 16.466 | 58 | 58 |
| Ordinal IPTW | 29060_003 | 29060_003_PAROXETINE | 5 | Active GEE/MSM prediction | Target 20 mg | 20 | 13.234 | 0.985 | 11.303 | 15.165 | 58 | 58 |
| Ordinal IPTW | 29060_003 | 29060_003_PAROXETINE | 5 | Active GEE/MSM prediction | Target 30 mg | 30 | 11.872 | 0.951 | 10.009 | 13.736 | 58 | 58 |
| Ordinal IPTW | 29060_003 | 29060_003_PAROXETINE | 5 | Active GEE/MSM prediction | Target 40 mg | 40 | 13.440 | 1.317 | 10.859 | 16.022 | 58 | 58 |
| Ordinal IPTW | 29060_003 | 29060_003_PAROXETINE | 5 | Active GEE/MSM prediction | Target 50 mg | 50 | 13.695 | 1.620 | 10.521 | 16.870 | 58 | 58 |
| Ordinal IPTW | 29060_003 | 29060_003_PAROXETINE | 6 | Observed placebo | Observed placebo | NA | 10.208 | 1.163 | 7.928 | 12.488 | 68 | 68 |
| Ordinal IPTW | 29060_003 | 29060_003_PAROXETINE | 6 | Active GEE/MSM prediction | Target 10 mg | 10 | 14.399 | 4.902 | 4.791 | 24.007 | 79 | 79 |
| Ordinal IPTW | 29060_003 | 29060_003_PAROXETINE | 6 | Active GEE/MSM prediction | Target 20 mg | 20 | 15.402 | 1.002 | 13.438 | 17.367 | 79 | 79 |
| Ordinal IPTW | 29060_003 | 29060_003_PAROXETINE | 6 | Active GEE/MSM prediction | Target 30 mg | 30 | 12.926 | 0.911 | 11.140 | 14.712 | 79 | 79 |
| Ordinal IPTW | 29060_003 | 29060_003_PAROXETINE | 6 | Active GEE/MSM prediction | Target 40 mg | 40 | 14.598 | 1.137 | 12.369 | 16.828 | 79 | 79 |
| Ordinal IPTW | 29060_003 | 29060_003_PAROXETINE | 6 | Active GEE/MSM prediction | Target 50 mg | 50 | 13.965 | 1.584 | 10.859 | 17.070 | 79 | 79 |
| Ordinal IPTW | 29060_003 | 29060_003_PAROXETINE | 7 | Observed placebo | Observed placebo | NA | 13.486 | 1.137 | 11.257 | 15.715 | 82 | 82 |
| Ordinal IPTW | 29060_003 | 29060_003_PAROXETINE | 7 | Active GEE/MSM prediction | Target 10 mg | 10 | 17.683 | 7.819 | 2.357 | 33.008 | 74 | 74 |
| Ordinal IPTW | 29060_003 | 29060_003_PAROXETINE | 7 | Active GEE/MSM prediction | Target 20 mg | 20 | 17.522 | 1.430 | 14.719 | 20.324 | 74 | 74 |
| Ordinal IPTW | 29060_003 | 29060_003_PAROXETINE | 7 | Active GEE/MSM prediction | Target 30 mg | 30 | 13.853 | 1.167 | 11.566 | 16.141 | 74 | 74 |
| Ordinal IPTW | 29060_003 | 29060_003_PAROXETINE | 7 | Active GEE/MSM prediction | Target 40 mg | 40 | 15.719 | 1.457 | 12.864 | 18.574 | 74 | 74 |
| Ordinal IPTW | 29060_003 | 29060_003_PAROXETINE | 7 | Active GEE/MSM prediction | Target 50 mg | 50 | 13.850 | 1.950 | 10.028 | 17.671 | 74 | 74 |
| Ordinal IPTW | 29060_003 | 29060_003_PAROXETINE | 8 | Observed placebo | Observed placebo | NA | 12.693 | NA | NA | NA | 1 | 1 |
| Ordinal IPTW | 29060_003 | 29060_003_PAROXETINE | 8 | Active GEE/MSM prediction | Target 10 mg | 10 | 21.045 | 10.849 | -0.218 | 42.309 | 6 | 6 |
| Ordinal IPTW | 29060_003 | 29060_003_PAROXETINE | 8 | Active GEE/MSM prediction | Target 20 mg | 20 | 19.633 | 2.067 | 15.582 | 23.683 | 6 | 6 |
| Ordinal IPTW | 29060_003 | 29060_003_PAROXETINE | 8 | Active GEE/MSM prediction | Target 30 mg | 30 | 14.760 | 1.630 | 11.565 | 17.956 | 6 | 6 |
| Ordinal IPTW | 29060_003 | 29060_003_PAROXETINE | 8 | Active GEE/MSM prediction | Target 40 mg | 40 | 16.833 | 2.139 | 12.641 | 21.025 | 6 | 6 |
| Ordinal IPTW | 29060_003 | 29060_003_PAROXETINE | 8 | Active GEE/MSM prediction | Target 50 mg | 50 | 13.671 | 2.603 | 8.569 | 18.772 | 6 | 6 |

<small><em>Active GEE/MSM predicted means and observed placebo means by
week</em></small>

## References

1.  J M Robins, M A Hernán, B Brumback. Marginal structural models and
    causal inference in epidemiology. Epidemiology. 2000
    Sep;11(5):550-60. doi: 10.1097/00001648-200009000-00011.
2.  Ilya Lipkovich, David H Adams, Craig Mallinckrodt, Doug Faries,
    David Baron, John P Houston. Evaluating dose response from flexible
    dose clinical trials. BMC Psychiatry. 2008 Jan 7;8:3. doi:
    10.1186/1471-244X-8-3.

## Appendix

### Predicted probabilities under the ordinal dose model

In the primary dose-weight model, dose was modeled as an ordered
categorical variable with levels $20 < 30 < 40 < 50$. For each
cumulative dose threshold $c \in \{20,30,40\}$, the ordinal logistic
model was written as

``` math
\log
\left[
\frac{
\Pr(D_{it} \le c \mid H_{it})
}{
1 - \Pr(D_{it} \le c \mid H_{it})
}
\right]
=
\alpha_{0c} + X_{it}\alpha.
```

Let

``` math
\eta_{it,c} = \alpha_{0c} + X_{it}\alpha.
```

Then

``` math
\log
\left[
\frac{
\Pr(D_{it} \le c \mid H_{it})
}{
1 - \Pr(D_{it} \le c \mid H_{it})
}
\right]
=
\eta_{it,c}.
```

Exponentiating both sides gives

``` math
\frac{
\Pr(D_{it} \le c \mid H_{it})
}{
1 - \Pr(D_{it} \le c \mid H_{it})
}
=
\exp(\eta_{it,c}).
```

Let

``` math
F_c =
\Pr(D_{it} \le c \mid H_{it}).
```

Then

``` math
\frac{F_c}{1-F_c}
=
\exp(\eta_{it,c}).
```

Solving for $F_c$ gives

``` math
F_c
=
\frac{\exp(\eta_{it,c})}{1+\exp(\eta_{it,c})}
=
\frac{1}{1+\exp(-\eta_{it,c})}.
```

Therefore, the ordinal model gives the cumulative predicted
probabilities

``` math
F_{20}
=
\Pr(D_{it} \le 20 \mid H_{it}),
```

``` math
F_{30}
=
\Pr(D_{it} \le 30 \mid H_{it}),
```

and

``` math
F_{40}
=
\Pr(D_{it} \le 40 \mid H_{it}).
```

The predicted probabilities for the individual dose categories are then
obtained as

``` math
\Pr(D_{it}=20 \mid H_{it})
=
F_{20},
```

``` math
\Pr(D_{it}=30 \mid H_{it})
=
F_{30} - F_{20},
```

``` math
\Pr(D_{it}=40 \mid H_{it})
=
F_{40} - F_{30},
```

and

``` math
\Pr(D_{it}=50 \mid H_{it})
=
1 - F_{40}.
```

### Alternative multinomial dose model for stabilized dose weights

As a sensitivity analysis to the ordinal dose model, we considered a
multinomial logistic regression for the dose-weight models. Unlike the
ordinal model, which treats the dose categories as ordered and assumes
proportional odds, the multinomial model treats dose as a nominal
categorical variable. This allows the associations between patient
history and dose assignment to differ across dose levels.

Using 20 mg as the reference category, the denominator multinomial dose
model was specified for each non-reference dose $d \in \{30,40,50\}$ as

``` math
\begin{aligned}
\log
\left[
\frac{
\Pr(D_{it}=d \mid H_{it})
}{
\Pr(D_{it}=20 \mid H_{it})
}
\right]
&=
\alpha_{0d}
+ \alpha_{1d} f(\mathrm{week}_{it})
+ \alpha_{2d} \Delta Y_{it}
+ \alpha_{3d} S_{it}
+ \alpha_{4d} d_{it-1} \\
&\quad
+ \alpha_{5d} Y_{i0}.
\end{aligned}
```

The fitted multinomial model provides predicted probabilities for all
dose categories. For each patient-visit row, the **denominator
probability** was the fitted probability corresponding to the dose
actually received:

``` math
\hat{p}_{it}^{D,\mathrm{mult}}
=
\widehat{\Pr}
\left(
D_{it}=d_{it}
\mid
H_{it}
\right).
```

The **numerator dose model** had the same multinomial structure but
excluded the time-varying confounders:

``` math
\begin{aligned}
\log
\left[
\frac{
\Pr(D_{it}=d \mid H_{it}^{*})
}{
\Pr(D_{it}=20 \mid H_{it}^{*})
}
\right]
&=
\beta_{0d}
+ \beta_{1d} f(\mathrm{week}_{it})
+ \beta_{2d} d_{it-1}
+ \beta_{3d} Y_{i0}.
\end{aligned}
```

The fitted probability corresponding to the observed dose $d_{it}$ was
extracted as

``` math
\hat{q}_{it}^{D,\mathrm{mult}}
=
\widehat{\Pr}
\left(
D_{it}=d_{it}
\mid
H_{it}^{*}
\right).
```

The visit-specific and cumulative stabilized dose weights were then
constructed analogously to the primary ordinal-model analysis:

``` math
\mathrm{SW}_{it}^{D,\mathrm{mult}}
=
\frac{
\hat{q}_{it}^{D,\mathrm{mult}}
}{
\hat{p}_{it}^{D,\mathrm{mult}}
}.
```

``` math
\mathrm{cSW}_{it}^{D,\mathrm{mult}}
=
\prod_{s=1}^{t}
\mathrm{SW}_{is}^{D,\mathrm{mult}}.
```

The multinomial model was used to assess the robustness of the results
to relaxing the proportional-odds assumption imposed by the ordinal dose
model. Whereas the ordinal model assumes that covariates shift the
probability toward higher or lower dose categories similarly across dose
thresholds, the multinomial model allows separate covariate effects for
each non-reference dose level relative to 20 mg.

#### Predicted probabilities under the multinomial dose model

For a multinomial logit model with $K$ possible dose categories, one
category is chosen as the reference. Let category $K$ denote the
reference category. For each non-reference category $k = 1,\ldots,K-1$,
the model can be written as

``` math
\log
\left[
\frac{
\Pr(Y_i = k \mid X_i)
}{
\Pr(Y_i = K \mid X_i)
}
\right]
=
\beta_k^{\top}X_i.
```

Exponentiating both sides gives

``` math
\Pr(Y_i = k \mid X_i)
=
\Pr(Y_i = K \mid X_i)
\exp(\beta_k^{\top}X_i),
\qquad
k = 1,\ldots,K-1.
```

Because the category probabilities must sum to one,

``` math
\begin{aligned}
\Pr(Y_i = K \mid X_i)
&=
1 - \sum_{j=1}^{K-1} \Pr(Y_i = j \mid X_i) \\
&=
1 - \sum_{j=1}^{K-1}
\Pr(Y_i = K \mid X_i)\exp(\beta_j^{\top}X_i).
\end{aligned}
```

Therefore,

``` math
\begin{aligned}
\Pr(Y_i = K \mid X_i)
+
\Pr(Y_i = K \mid X_i)
\sum_{j=1}^{K-1}\exp(\beta_j^{\top}X_i)
&=
1, \\
\Pr(Y_i = K \mid X_i)
\left[
1 + \sum_{j=1}^{K-1}\exp(\beta_j^{\top}X_i)
\right]
&=
1.
\end{aligned}
```

Thus, for the reference category,

``` math
\Pr(Y_i = K \mid X_i)
=
\frac{
1
}{
1 + \sum_{j=1}^{K-1}\exp(\beta_j^{\top}X_i)
}.
```

For each non-reference category $k = 1,\ldots,K-1$, we obtain

``` math
\Pr(Y_i = k \mid X_i)
=
\frac{
\exp(\beta_k^{\top}X_i)
}{
1 + \sum_{j=1}^{K-1}\exp(\beta_j^{\top}X_i)
},
\qquad
1 \le k < K.
```

Thus, the multinomial model estimates log-odds relative to the reference
dose category. It also relies on the independence of irrelevant
alternatives assumption, meaning that the odds comparing two dose
categories are assumed not to depend on the presence or characteristics
of the other dose categories.

#### Reasults using multiomial model for the IPTW step

| Dose_level_vs_reference | Term | Estimate | SE | z_value | p_value |
|:---|:---|---:|---:|---:|---:|
| 20 | (Intercept) | 3.579 | 2.678 | 1.336 | 0.181 |
| 20 | rms::rcs(visit, visit_df)visit | -0.070 | 0.436 | -0.162 | 0.872 |
| 20 | rms::rcs(visit, visit_df)visit’ | -0.065 | 0.422 | -0.155 | 0.877 |
| 20 | delta_outcome | -0.003 | 0.053 | -0.057 | 0.955 |
| 20 | side.effects | -0.298 | 0.190 | -1.567 | 0.117 |
| 20 | dose_lag1_f20 | 15.727 | 8.908 | 1.765 | 0.077 |
| 20 | dose_lag1_f30 | 0.002 | 2.700 | 0.001 | 0.999 |
| 20 | dose_lag1_f40 | 169.165 | 1.374 | 123.084 | 0.000 |
| 20 | dose_lag1_f50 | 262.122 | 1.362 | 192.523 | 0.000 |
| 20 | outcome_0 | -0.050 | 0.074 | -0.673 | 0.501 |
| 20 | delta_outcome:dose_lag1_f20 | -0.033 | 0.066 | -0.503 | 0.615 |
| 20 | delta_outcome:dose_lag1_f30 | 0.017 | 0.076 | 0.225 | 0.822 |
| 20 | delta_outcome:dose_lag1_f40 | 0.146 | 0.134 | 1.089 | 0.276 |
| 20 | delta_outcome:dose_lag1_f50 | 14.596 | 0.143 | 102.019 | 0.000 |
| 20 | side.effects:dose_lag1_f20 | -1.283 | 0.915 | -1.402 | 0.161 |
| 20 | side.effects:dose_lag1_f30 | 0.376 | 0.333 | 1.130 | 0.259 |
| 20 | side.effects:dose_lag1_f40 | -16.714 | 0.202 | -82.801 | 0.000 |
| 20 | side.effects:dose_lag1_f50 | -31.934 | 0.200 | -159.715 | 0.000 |
| 30 | (Intercept) | 5.361 | 3.045 | 1.761 | 0.078 |
| 30 | rms::rcs(visit, visit_df)visit | -0.116 | 0.450 | -0.257 | 0.797 |
| 30 | rms::rcs(visit, visit_df)visit’ | -0.041 | 0.439 | -0.094 | 0.925 |
| 30 | delta_outcome | -0.021 | 0.073 | -0.290 | 0.772 |
| 30 | side.effects | -0.679 | 0.298 | -2.280 | 0.023 |
| 30 | dose_lag1_f20 | 15.008 | 9.004 | 1.667 | 0.096 |
| 30 | dose_lag1_f30 | 1.799 | 2.938 | 0.613 | 0.540 |
| 30 | dose_lag1_f40 | 168.442 | 1.145 | 147.123 | 0.000 |
| 30 | dose_lag1_f50 | 261.660 | 1.132 | 231.077 | 0.000 |
| 30 | outcome_0 | -0.065 | 0.078 | -0.835 | 0.404 |
| 30 | delta_outcome:dose_lag1_f20 | -0.028 | 0.085 | -0.328 | 0.743 |
| 30 | delta_outcome:dose_lag1_f30 | 0.043 | 0.091 | 0.473 | 0.636 |
| 30 | delta_outcome:dose_lag1_f40 | 0.167 | 0.144 | 1.164 | 0.245 |
| 30 | delta_outcome:dose_lag1_f50 | 14.635 | 0.142 | 103.350 | 0.000 |
| 30 | side.effects:dose_lag1_f20 | -1.140 | 0.946 | -1.205 | 0.228 |
| 30 | side.effects:dose_lag1_f30 | 0.325 | 0.403 | 0.807 | 0.419 |
| 30 | side.effects:dose_lag1_f40 | -16.425 | 0.249 | -65.835 | 0.000 |
| 30 | side.effects:dose_lag1_f50 | -31.833 | 0.210 | -151.386 | 0.000 |
| 40 | (Intercept) | 5.930 | 3.353 | 1.769 | 0.077 |
| 40 | rms::rcs(visit, visit_df)visit | -0.027 | 0.463 | -0.059 | 0.953 |
| 40 | rms::rcs(visit, visit_df)visit’ | -0.248 | 0.453 | -0.548 | 0.584 |
| 40 | delta_outcome | -0.038 | 0.096 | -0.400 | 0.689 |
| 40 | side.effects | -1.098 | 0.508 | -2.163 | 0.031 |
| 40 | dose_lag1_f20 | 14.251 | 9.091 | 1.568 | 0.117 |
| 40 | dose_lag1_f30 | 0.984 | 3.211 | 0.307 | 0.759 |
| 40 | dose_lag1_f40 | 169.982 | 1.169 | 145.416 | 0.000 |
| 40 | dose_lag1_f50 | 261.223 | 1.192 | 219.231 | 0.000 |
| 40 | outcome_0 | -0.068 | 0.080 | -0.848 | 0.396 |
| 40 | delta_outcome:dose_lag1_f20 | 0.053 | 0.107 | 0.495 | 0.621 |
| 40 | delta_outcome:dose_lag1_f30 | 0.088 | 0.113 | 0.783 | 0.433 |
| 40 | delta_outcome:dose_lag1_f40 | 0.169 | 0.157 | 1.076 | 0.282 |
| 40 | delta_outcome:dose_lag1_f50 | 14.708 | 0.144 | 102.360 | 0.000 |
| 40 | side.effects:dose_lag1_f20 | -0.916 | 1.033 | -0.887 | 0.375 |
| 40 | side.effects:dose_lag1_f30 | 0.540 | 0.581 | 0.929 | 0.353 |
| 40 | side.effects:dose_lag1_f40 | -16.309 | 0.446 | -36.606 | 0.000 |
| 40 | side.effects:dose_lag1_f50 | -31.560 | 0.349 | -90.370 | 0.000 |
| 50 | (Intercept) | 6.446 | 3.617 | 1.782 | 0.075 |
| 50 | rms::rcs(visit, visit_df)visit | -0.202 | 0.470 | -0.430 | 0.667 |
| 50 | rms::rcs(visit, visit_df)visit’ | 0.049 | 0.461 | 0.105 | 0.916 |
| 50 | delta_outcome | -0.276 | 0.535 | -0.516 | 0.606 |
| 50 | side.effects | -42.757 | 0.222 | -192.312 | 0.000 |
| 50 | dose_lag1_f20 | 14.183 | 9.183 | 1.544 | 0.122 |
| 50 | dose_lag1_f30 | 0.831 | 3.469 | 0.239 | 0.811 |
| 50 | dose_lag1_f40 | 168.259 | 1.498 | 112.332 | 0.000 |
| 50 | dose_lag1_f50 | 261.979 | 1.469 | 178.382 | 0.000 |
| 50 | outcome_0 | -0.028 | 0.082 | -0.347 | 0.729 |
| 50 | delta_outcome:dose_lag1_f20 | 0.246 | 0.537 | 0.458 | 0.647 |
| 50 | delta_outcome:dose_lag1_f30 | 0.297 | 0.538 | 0.552 | 0.581 |
| 50 | delta_outcome:dose_lag1_f40 | 0.399 | 0.549 | 0.727 | 0.467 |
| 50 | delta_outcome:dose_lag1_f50 | 14.849 | 0.397 | 37.377 | 0.000 |
| 50 | side.effects:dose_lag1_f20 | 40.493 | 0.721 | 56.145 | 0.000 |
| 50 | side.effects:dose_lag1_f30 | 41.783 | 0.324 | 128.899 | 0.000 |
| 50 | side.effects:dose_lag1_f40 | 24.765 | 0.308 | 80.398 | 0.000 |
| 50 | side.effects:dose_lag1_f50 | 9.929 | 0.219 | 45.337 | 0.000 |

<small><em>Multinomial dose-weight denominator model:
29060_002_PAROXETINE</em></small>

| Dose_level_vs_reference | Term | Estimate | SE | z_value | p_value |
|:---|:---|---:|---:|---:|---:|
| 65 | (Intercept) | -4.624 | 3.093 | -1.495 | 0.135 |
| 65 | rms::rcs(visit, visit_df)visit | -0.458 | 0.452 | -1.013 | 0.311 |
| 65 | rms::rcs(visit, visit_df)visit’ | 0.867 | 0.612 | 1.416 | 0.157 |
| 65 | delta_outcome | -0.048 | 0.082 | -0.588 | 0.557 |
| 65 | side.effects | -0.205 | 0.239 | -0.856 | 0.392 |
| 65 | dose_lag1_f65 | 11.604 | 5.122 | 2.265 | 0.023 |
| 65 | dose_lag1_f80 | 0.052 | 0.275 | 0.188 | 0.851 |
| 65 | dose_lag1_f145 | 1.534 | 2.545 | 0.603 | 0.547 |
| 65 | dose_lag1_f210 | 5.050 | 3.069 | 1.646 | 0.100 |
| 65 | dose_lag1_f275 | -2.036 | 0.067 | -30.213 | 0.000 |
| 65 | outcome_0 | 0.113 | 0.084 | 1.341 | 0.180 |
| 65 | delta_outcome:dose_lag1_f65 | -0.027 | 0.102 | -0.266 | 0.790 |
| 65 | delta_outcome:dose_lag1_f80 | 0.392 | 5.535 | 0.071 | 0.944 |
| 65 | delta_outcome:dose_lag1_f145 | -0.120 | 0.105 | -1.142 | 0.253 |
| 65 | delta_outcome:dose_lag1_f210 | 0.153 | 0.134 | 1.141 | 0.254 |
| 65 | delta_outcome:dose_lag1_f275 | -0.141 | 5.646 | -0.025 | 0.980 |
| 65 | side.effects:dose_lag1_f65 | -0.633 | 0.530 | -1.193 | 0.233 |
| 65 | side.effects:dose_lag1_f80 | -1.389 | 2.622 | -0.530 | 0.596 |
| 65 | side.effects:dose_lag1_f145 | 0.113 | 0.320 | 0.354 | 0.723 |
| 65 | side.effects:dose_lag1_f210 | -0.712 | 0.451 | -1.579 | 0.114 |
| 65 | side.effects:dose_lag1_f275 | -1.361 | 0.073 | -18.634 | 0.000 |
| 80 | (Intercept) | -1.812 | 2.285 | -0.793 | 0.428 |
| 80 | rms::rcs(visit, visit_df)visit | -0.751 | 0.420 | -1.788 | 0.074 |
| 80 | rms::rcs(visit, visit_df)visit’ | 1.151 | 0.560 | 2.055 | 0.040 |
| 80 | delta_outcome | 0.000 | 0.051 | 0.006 | 0.995 |
| 80 | side.effects | -0.267 | 0.138 | -1.938 | 0.053 |
| 80 | dose_lag1_f65 | 9.177 | 5.666 | 1.619 | 0.105 |
| 80 | dose_lag1_f80 | 23.053 | 1.276 | 18.070 | 0.000 |
| 80 | dose_lag1_f145 | 1.097 | 1.702 | 0.645 | 0.519 |
| 80 | dose_lag1_f210 | 3.216 | 2.682 | 1.199 | 0.230 |
| 80 | dose_lag1_f275 | 2.746 | 2.475 | 1.110 | 0.267 |
| 80 | outcome_0 | 0.076 | 0.069 | 1.098 | 0.272 |
| 80 | delta_outcome:dose_lag1_f65 | -0.102 | 0.150 | -0.680 | 0.497 |
| 80 | delta_outcome:dose_lag1_f80 | 0.106 | 0.113 | 0.945 | 0.345 |
| 80 | delta_outcome:dose_lag1_f145 | -0.012 | 0.068 | -0.170 | 0.865 |
| 80 | delta_outcome:dose_lag1_f210 | 0.038 | 0.116 | 0.331 | 0.741 |
| 80 | delta_outcome:dose_lag1_f275 | -0.029 | 0.160 | -0.180 | 0.857 |
| 80 | side.effects:dose_lag1_f65 | -0.808 | 0.615 | -1.314 | 0.189 |
| 80 | side.effects:dose_lag1_f80 | -2.221 | 0.217 | -10.237 | 0.000 |
| 80 | side.effects:dose_lag1_f145 | -0.082 | 0.230 | -0.357 | 0.721 |
| 80 | side.effects:dose_lag1_f210 | -0.411 | 0.404 | -1.018 | 0.309 |
| 80 | side.effects:dose_lag1_f275 | -0.560 | 0.499 | -1.123 | 0.262 |
| 145 | (Intercept) | 1.761 | 1.094 | 1.610 | 0.107 |
| 145 | rms::rcs(visit, visit_df)visit | -0.324 | 0.207 | -1.566 | 0.117 |
| 145 | rms::rcs(visit, visit_df)visit’ | 0.604 | 0.271 | 2.228 | 0.026 |
| 145 | delta_outcome | 0.065 | 0.025 | 2.645 | 0.008 |
| 145 | side.effects | -0.384 | 0.070 | -5.474 | 0.000 |
| 145 | dose_lag1_f65 | 8.227 | 4.980 | 1.652 | 0.099 |
| 145 | dose_lag1_f80 | 22.170 | 0.939 | 23.615 | 0.000 |
| 145 | dose_lag1_f145 | 3.095 | 0.852 | 3.634 | 0.000 |
| 145 | dose_lag1_f210 | 4.568 | 1.779 | 2.568 | 0.010 |
| 145 | dose_lag1_f275 | 1.865 | 1.385 | 1.346 | 0.178 |
| 145 | outcome_0 | -0.009 | 0.033 | -0.281 | 0.779 |
| 145 | delta_outcome:dose_lag1_f65 | -0.150 | 0.087 | -1.719 | 0.086 |
| 145 | delta_outcome:dose_lag1_f80 | 0.013 | 0.119 | 0.110 | 0.912 |
| 145 | delta_outcome:dose_lag1_f145 | -0.079 | 0.029 | -2.679 | 0.007 |
| 145 | delta_outcome:dose_lag1_f210 | 0.086 | 0.065 | 1.307 | 0.191 |
| 145 | delta_outcome:dose_lag1_f275 | -0.073 | 0.062 | -1.178 | 0.239 |
| 145 | side.effects:dose_lag1_f65 | -0.626 | 0.508 | -1.231 | 0.219 |
| 145 | side.effects:dose_lag1_f80 | -2.485 | 0.223 | -11.141 | 0.000 |
| 145 | side.effects:dose_lag1_f145 | -0.094 | 0.108 | -0.875 | 0.382 |
| 145 | side.effects:dose_lag1_f210 | -0.730 | 0.283 | -2.582 | 0.010 |
| 145 | side.effects:dose_lag1_f275 | -0.110 | 0.180 | -0.609 | 0.542 |
| 210 | (Intercept) | 2.741 | 1.303 | 2.102 | 0.036 |
| 210 | rms::rcs(visit, visit_df)visit | -0.537 | 0.252 | -2.137 | 0.033 |
| 210 | rms::rcs(visit, visit_df)visit’ | 0.861 | 0.331 | 2.599 | 0.009 |
| 210 | delta_outcome | 0.075 | 0.027 | 2.744 | 0.006 |
| 210 | side.effects | -0.579 | 0.090 | -6.446 | 0.000 |
| 210 | dose_lag1_f65 | 9.866 | 5.232 | 1.886 | 0.059 |
| 210 | dose_lag1_f80 | 22.608 | 0.805 | 28.082 | 0.000 |
| 210 | dose_lag1_f145 | 1.916 | 0.920 | 2.083 | 0.037 |
| 210 | dose_lag1_f210 | 6.761 | 1.766 | 3.828 | 0.000 |
| 210 | dose_lag1_f275 | 0.272 | 1.582 | 0.172 | 0.864 |
| 210 | outcome_0 | -0.023 | 0.040 | -0.572 | 0.568 |
| 210 | delta_outcome:dose_lag1_f65 | -0.241 | 0.138 | -1.748 | 0.080 |
| 210 | delta_outcome:dose_lag1_f80 | -0.017 | 0.113 | -0.152 | 0.879 |
| 210 | delta_outcome:dose_lag1_f145 | -0.048 | 0.034 | -1.405 | 0.160 |
| 210 | delta_outcome:dose_lag1_f210 | 0.036 | 0.065 | 0.559 | 0.576 |
| 210 | delta_outcome:dose_lag1_f275 | -0.100 | 0.071 | -1.411 | 0.158 |
| 210 | side.effects:dose_lag1_f65 | -0.780 | 0.560 | -1.393 | 0.164 |
| 210 | side.effects:dose_lag1_f80 | -2.241 | 0.201 | -11.164 | 0.000 |
| 210 | side.effects:dose_lag1_f145 | -0.226 | 0.153 | -1.474 | 0.141 |
| 210 | side.effects:dose_lag1_f210 | -0.773 | 0.290 | -2.664 | 0.008 |
| 210 | side.effects:dose_lag1_f275 | 0.209 | 0.212 | 0.983 | 0.325 |
| 275 | (Intercept) | 2.623 | 1.702 | 1.541 | 0.123 |
| 275 | rms::rcs(visit, visit_df)visit | -0.032 | 0.341 | -0.094 | 0.925 |
| 275 | rms::rcs(visit, visit_df)visit’ | 0.282 | 0.443 | 0.637 | 0.524 |
| 275 | delta_outcome | 0.010 | 0.043 | 0.239 | 0.811 |
| 275 | side.effects | -0.851 | 0.157 | -5.431 | 0.000 |
| 275 | dose_lag1_f65 | 7.962 | 5.131 | 1.552 | 0.121 |
| 275 | dose_lag1_f80 | 23.293 | 1.047 | 22.244 | 0.000 |
| 275 | dose_lag1_f145 | 1.019 | 1.074 | 0.949 | 0.342 |
| 275 | dose_lag1_f210 | 4.017 | 1.943 | 2.068 | 0.039 |
| 275 | dose_lag1_f275 | 2.916 | 1.463 | 1.993 | 0.046 |
| 275 | outcome_0 | -0.036 | 0.054 | -0.662 | 0.508 |
| 275 | delta_outcome:dose_lag1_f65 | -0.080 | 0.109 | -0.734 | 0.463 |
| 275 | delta_outcome:dose_lag1_f80 | -0.023 | 0.142 | -0.159 | 0.874 |
| 275 | delta_outcome:dose_lag1_f145 | 0.030 | 0.057 | 0.520 | 0.603 |
| 275 | delta_outcome:dose_lag1_f210 | 0.086 | 0.088 | 0.986 | 0.324 |
| 275 | delta_outcome:dose_lag1_f275 | -0.041 | 0.077 | -0.538 | 0.591 |
| 275 | side.effects:dose_lag1_f65 | -0.488 | 0.569 | -0.857 | 0.392 |
| 275 | side.effects:dose_lag1_f80 | -3.176 | 0.886 | -3.584 | 0.000 |
| 275 | side.effects:dose_lag1_f145 | -0.568 | 0.357 | -1.592 | 0.111 |
| 275 | side.effects:dose_lag1_f210 | -0.684 | 0.393 | -1.743 | 0.081 |
| 275 | side.effects:dose_lag1_f275 | -0.015 | 0.257 | -0.058 | 0.954 |

<small><em>Multinomial dose-weight denominator model:
29060_003_IMIPRAMINE</em></small>

| Dose_level_vs_reference | Term | Estimate | SE | z_value | p_value |
|:---|:---|---:|---:|---:|---:|
| 20 | (Intercept) | -2.546 | 2.751 | -0.925 | 0.355 |
| 20 | rms::rcs(visit, visit_df)visit | 0.561 | 0.371 | 1.513 | 0.130 |
| 20 | rms::rcs(visit, visit_df)visit’ | -1.200 | 0.599 | -2.002 | 0.045 |
| 20 | delta_outcome | -0.007 | 0.070 | -0.106 | 0.916 |
| 20 | side.effects | -0.069 | 0.254 | -0.271 | 0.787 |
| 20 | dose_lag1_f20 | 3.729 | 2.513 | 1.484 | 0.138 |
| 20 | dose_lag1_f30 | 123.291 | 1.563 | 78.884 | 0.000 |
| 20 | dose_lag1_f40 | 1.951 | 2.703 | 0.722 | 0.470 |
| 20 | dose_lag1_f50 | 3.323 | 3.052 | 1.089 | 0.276 |
| 20 | outcome_0 | 0.051 | 0.056 | 0.915 | 0.360 |
| 20 | delta_outcome:dose_lag1_f20 | 0.030 | 0.079 | 0.381 | 0.703 |
| 20 | delta_outcome:dose_lag1_f30 | 0.004 | 0.097 | 0.045 | 0.964 |
| 20 | delta_outcome:dose_lag1_f40 | 0.045 | 0.096 | 0.466 | 0.641 |
| 20 | delta_outcome:dose_lag1_f50 | 0.056 | 0.097 | 0.573 | 0.567 |
| 20 | side.effects:dose_lag1_f20 | 0.019 | 0.293 | 0.063 | 0.949 |
| 20 | side.effects:dose_lag1_f30 | -12.063 | 0.211 | -57.289 | 0.000 |
| 20 | side.effects:dose_lag1_f40 | 0.171 | 0.336 | 0.509 | 0.611 |
| 20 | side.effects:dose_lag1_f50 | -0.172 | 0.353 | -0.487 | 0.626 |
| 30 | (Intercept) | 0.207 | 2.242 | 0.092 | 0.927 |
| 30 | rms::rcs(visit, visit_df)visit | 0.655 | 0.380 | 1.723 | 0.085 |
| 30 | rms::rcs(visit, visit_df)visit’ | -1.485 | 0.615 | -2.415 | 0.016 |
| 30 | delta_outcome | 0.130 | 0.077 | 1.680 | 0.093 |
| 30 | side.effects | -0.503 | 0.214 | -2.347 | 0.019 |
| 30 | dose_lag1_f20 | 1.733 | 1.982 | 0.874 | 0.382 |
| 30 | dose_lag1_f30 | 123.798 | 0.931 | 132.986 | 0.000 |
| 30 | dose_lag1_f40 | 0.734 | 2.202 | 0.333 | 0.739 |
| 30 | dose_lag1_f50 | 1.706 | 2.627 | 0.649 | 0.516 |
| 30 | outcome_0 | 0.017 | 0.058 | 0.294 | 0.769 |
| 30 | delta_outcome:dose_lag1_f20 | -0.086 | 0.087 | -0.992 | 0.321 |
| 30 | delta_outcome:dose_lag1_f30 | -0.105 | 0.102 | -1.029 | 0.303 |
| 30 | delta_outcome:dose_lag1_f40 | -0.137 | 0.102 | -1.344 | 0.179 |
| 30 | delta_outcome:dose_lag1_f50 | -0.109 | 0.103 | -1.053 | 0.292 |
| 30 | side.effects:dose_lag1_f20 | 0.260 | 0.261 | 0.994 | 0.320 |
| 30 | side.effects:dose_lag1_f30 | -11.919 | 0.193 | -61.742 | 0.000 |
| 30 | side.effects:dose_lag1_f40 | 0.469 | 0.310 | 1.514 | 0.130 |
| 30 | side.effects:dose_lag1_f50 | 0.231 | 0.328 | 0.703 | 0.482 |
| 40 | (Intercept) | 0.865 | 2.427 | 0.356 | 0.721 |
| 40 | rms::rcs(visit, visit_df)visit | 0.329 | 0.392 | 0.839 | 0.402 |
| 40 | rms::rcs(visit, visit_df)visit’ | -0.798 | 0.635 | -1.258 | 0.208 |
| 40 | delta_outcome | -0.090 | 0.111 | -0.810 | 0.418 |
| 40 | side.effects | -0.547 | 0.241 | -2.264 | 0.024 |
| 40 | dose_lag1_f20 | 0.618 | 2.113 | 0.292 | 0.770 |
| 40 | dose_lag1_f30 | 121.074 | 1.056 | 114.706 | 0.000 |
| 40 | dose_lag1_f40 | 1.656 | 2.276 | 0.728 | 0.467 |
| 40 | dose_lag1_f50 | 1.282 | 2.746 | 0.467 | 0.641 |
| 40 | outcome_0 | 0.070 | 0.061 | 1.152 | 0.249 |
| 40 | delta_outcome:dose_lag1_f20 | 0.111 | 0.120 | 0.930 | 0.352 |
| 40 | delta_outcome:dose_lag1_f30 | 0.118 | 0.133 | 0.885 | 0.376 |
| 40 | delta_outcome:dose_lag1_f40 | 0.127 | 0.129 | 0.984 | 0.325 |
| 40 | delta_outcome:dose_lag1_f50 | 0.112 | 0.135 | 0.828 | 0.408 |
| 40 | side.effects:dose_lag1_f20 | 0.109 | 0.289 | 0.376 | 0.707 |
| 40 | side.effects:dose_lag1_f30 | -11.985 | 0.211 | -56.801 | 0.000 |
| 40 | side.effects:dose_lag1_f40 | 0.157 | 0.326 | 0.481 | 0.631 |
| 40 | side.effects:dose_lag1_f50 | -0.118 | 0.366 | -0.323 | 0.747 |
| 50 | (Intercept) | -0.123 | 2.483 | -0.049 | 0.961 |
| 50 | rms::rcs(visit, visit_df)visit | 0.638 | 0.421 | 1.515 | 0.130 |
| 50 | rms::rcs(visit, visit_df)visit’ | -1.353 | 0.676 | -2.002 | 0.045 |
| 50 | delta_outcome | 0.060 | 0.118 | 0.510 | 0.610 |
| 50 | side.effects | -0.860 | 0.294 | -2.923 | 0.003 |
| 50 | dose_lag1_f20 | 0.814 | 2.092 | 0.389 | 0.697 |
| 50 | dose_lag1_f30 | 121.851 | 1.011 | 120.486 | 0.000 |
| 50 | dose_lag1_f40 | -0.881 | 2.360 | -0.373 | 0.709 |
| 50 | dose_lag1_f50 | 3.037 | 2.686 | 1.131 | 0.258 |
| 50 | outcome_0 | 0.091 | 0.066 | 1.376 | 0.169 |
| 50 | delta_outcome:dose_lag1_f20 | -0.030 | 0.126 | -0.237 | 0.813 |
| 50 | delta_outcome:dose_lag1_f30 | -0.052 | 0.140 | -0.373 | 0.709 |
| 50 | delta_outcome:dose_lag1_f40 | 0.049 | 0.142 | 0.348 | 0.728 |
| 50 | delta_outcome:dose_lag1_f50 | -0.058 | 0.138 | -0.421 | 0.674 |
| 50 | side.effects:dose_lag1_f20 | 0.077 | 0.352 | 0.219 | 0.826 |
| 50 | side.effects:dose_lag1_f30 | -12.123 | 0.304 | -39.855 | 0.000 |
| 50 | side.effects:dose_lag1_f40 | 0.264 | 0.394 | 0.669 | 0.503 |
| 50 | side.effects:dose_lag1_f50 | -0.049 | 0.399 | -0.123 | 0.902 |

<small><em>Multinomial dose-weight denominator model:
29060_003_PAROXETINE</em></small>

| Dose_level_vs_reference | Term | Estimate | SE | z_value | p_value |
|:---|:---|---:|---:|---:|---:|
| 20 | (Intercept) | -0.243 | 1.670 | -0.145 | 0.884 |
| 20 | rms::rcs(visit, visit_df)visit | 0.239 | 0.318 | 0.751 | 0.452 |
| 20 | rms::rcs(visit, visit_df)visit’ | -0.292 | 0.327 | -0.894 | 0.371 |
| 20 | dose_lag1_f20 | 3.073 | 0.615 | 4.997 | 0.000 |
| 20 | dose_lag1_f30 | 2.936 | 0.758 | 3.875 | 0.000 |
| 20 | dose_lag1_f40 | 3.567 | 1.113 | 3.204 | 0.001 |
| 20 | dose_lag1_f50 | 3.047 | 0.861 | 3.539 | 0.000 |
| 20 | outcome_0 | -0.022 | 0.050 | -0.441 | 0.659 |
| 30 | (Intercept) | -0.630 | 1.755 | -0.359 | 0.720 |
| 30 | rms::rcs(visit, visit_df)visit | 0.207 | 0.327 | 0.635 | 0.526 |
| 30 | rms::rcs(visit, visit_df)visit’ | -0.285 | 0.336 | -0.847 | 0.397 |
| 30 | dose_lag1_f20 | 2.855 | 0.733 | 3.895 | 0.000 |
| 30 | dose_lag1_f30 | 4.041 | 0.841 | 4.805 | 0.000 |
| 30 | dose_lag1_f40 | 4.364 | 1.175 | 3.716 | 0.000 |
| 30 | dose_lag1_f50 | 3.345 | 0.951 | 3.518 | 0.000 |
| 30 | outcome_0 | -0.033 | 0.052 | -0.635 | 0.525 |
| 40 | (Intercept) | -2.804 | 1.890 | -1.484 | 0.138 |
| 40 | rms::rcs(visit, visit_df)visit | 0.390 | 0.333 | 1.169 | 0.242 |
| 40 | rms::rcs(visit, visit_df)visit’ | -0.548 | 0.346 | -1.584 | 0.113 |
| 40 | dose_lag1_f20 | 3.291 | 0.899 | 3.660 | 0.000 |
| 40 | dose_lag1_f30 | 3.954 | 0.999 | 3.958 | 0.000 |
| 40 | dose_lag1_f40 | 5.858 | 1.272 | 4.607 | 0.000 |
| 40 | dose_lag1_f50 | 4.318 | 1.076 | 4.015 | 0.000 |
| 40 | outcome_0 | 0.007 | 0.053 | 0.135 | 0.892 |
| 50 | (Intercept) | -2.946 | 2.009 | -1.467 | 0.142 |
| 50 | rms::rcs(visit, visit_df)visit | 0.222 | 0.332 | 0.670 | 0.503 |
| 50 | rms::rcs(visit, visit_df)visit’ | -0.358 | 0.344 | -1.041 | 0.298 |
| 50 | dose_lag1_f20 | 4.111 | 1.139 | 3.611 | 0.000 |
| 50 | dose_lag1_f30 | 4.588 | 1.221 | 3.757 | 0.000 |
| 50 | dose_lag1_f40 | 5.248 | 1.470 | 3.569 | 0.000 |
| 50 | dose_lag1_f50 | 5.902 | 1.273 | 4.637 | 0.000 |
| 50 | outcome_0 | 0.002 | 0.053 | 0.038 | 0.970 |

<small><em>Multinomial dose-weight numerator model:
29060_002_PAROXETINE</em></small>

| Dose_level_vs_reference | Term | Estimate | SE | z_value | p_value |
|:---|:---|---:|---:|---:|---:|
| 65 | (Intercept) | -5.008 | 2.061 | -2.429 | 0.015 |
| 65 | rms::rcs(visit, visit_df)visit | -0.476 | 0.417 | -1.142 | 0.254 |
| 65 | rms::rcs(visit, visit_df)visit’ | 0.697 | 0.566 | 1.232 | 0.218 |
| 65 | dose_lag1_f65 | 5.666 | 0.904 | 6.270 | 0.000 |
| 65 | dose_lag1_f80 | -9.696 | 0.000 | -1032841.491 | 0.000 |
| 65 | dose_lag1_f145 | 1.679 | 0.887 | 1.892 | 0.058 |
| 65 | dose_lag1_f210 | 1.493 | 1.024 | 1.457 | 0.145 |
| 65 | dose_lag1_f275 | -12.714 | 0.000 | -7651531.275 | 0.000 |
| 65 | outcome_0 | 0.064 | 0.064 | 1.011 | 0.312 |
| 80 | (Intercept) | -3.979 | 1.836 | -2.167 | 0.030 |
| 80 | rms::rcs(visit, visit_df)visit | -0.724 | 0.390 | -1.854 | 0.064 |
| 80 | rms::rcs(visit, visit_df)visit’ | 1.091 | 0.509 | 2.144 | 0.032 |
| 80 | dose_lag1_f65 | 1.593 | 1.200 | 1.327 | 0.184 |
| 80 | dose_lag1_f80 | 3.884 | 0.699 | 5.553 | 0.000 |
| 80 | dose_lag1_f145 | 0.598 | 0.675 | 0.886 | 0.376 |
| 80 | dose_lag1_f210 | 0.434 | 0.848 | 0.512 | 0.609 |
| 80 | dose_lag1_f275 | 0.287 | 1.122 | 0.256 | 0.798 |
| 80 | outcome_0 | 0.080 | 0.057 | 1.396 | 0.163 |
| 145 | (Intercept) | -0.965 | 0.757 | -1.275 | 0.202 |
| 145 | rms::rcs(visit, visit_df)visit | -0.159 | 0.169 | -0.942 | 0.346 |
| 145 | rms::rcs(visit, visit_df)visit’ | 0.349 | 0.224 | 1.558 | 0.119 |
| 145 | dose_lag1_f65 | 1.163 | 0.694 | 1.674 | 0.094 |
| 145 | dose_lag1_f80 | 0.728 | 0.695 | 1.048 | 0.295 |
| 145 | dose_lag1_f145 | 1.900 | 0.245 | 7.762 | 0.000 |
| 145 | dose_lag1_f210 | 0.634 | 0.336 | 1.885 | 0.059 |
| 145 | dose_lag1_f275 | 0.605 | 0.423 | 1.431 | 0.152 |
| 145 | outcome_0 | 0.008 | 0.024 | 0.319 | 0.750 |
| 210 | (Intercept) | -1.452 | 0.886 | -1.639 | 0.101 |
| 210 | rms::rcs(visit, visit_df)visit | -0.274 | 0.192 | -1.424 | 0.155 |
| 210 | rms::rcs(visit, visit_df)visit’ | 0.499 | 0.255 | 1.960 | 0.050 |
| 210 | dose_lag1_f65 | 1.130 | 0.794 | 1.424 | 0.154 |
| 210 | dose_lag1_f80 | 1.785 | 0.624 | 2.859 | 0.004 |
| 210 | dose_lag1_f145 | 0.736 | 0.328 | 2.245 | 0.025 |
| 210 | dose_lag1_f210 | 2.278 | 0.308 | 7.396 | 0.000 |
| 210 | dose_lag1_f275 | 0.224 | 0.557 | 0.402 | 0.688 |
| 210 | outcome_0 | 0.016 | 0.028 | 0.569 | 0.569 |
| 275 | (Intercept) | -1.924 | 1.154 | -1.667 | 0.095 |
| 275 | rms::rcs(visit, visit_df)visit | 0.226 | 0.264 | 0.856 | 0.392 |
| 275 | rms::rcs(visit, visit_df)visit’ | -0.244 | 0.350 | -0.698 | 0.485 |
| 275 | dose_lag1_f65 | 1.804 | 0.813 | 2.218 | 0.027 |
| 275 | dose_lag1_f80 | 1.812 | 0.726 | 2.495 | 0.013 |
| 275 | dose_lag1_f145 | 0.560 | 0.434 | 1.291 | 0.197 |
| 275 | dose_lag1_f210 | 0.270 | 0.558 | 0.485 | 0.628 |
| 275 | dose_lag1_f275 | 2.417 | 0.429 | 5.640 | 0.000 |
| 275 | outcome_0 | -0.024 | 0.036 | -0.672 | 0.502 |

<small><em>Multinomial dose-weight numerator model:
29060_003_IMIPRAMINE</em></small>

| Dose_level_vs_reference | Term | Estimate | SE | z_value | p_value |
|:---|:---|---:|---:|---:|---:|
| 20 | (Intercept) | -3.267 | 1.511 | -2.162 | 0.031 |
| 20 | rms::rcs(visit, visit_df)visit | 0.400 | 0.338 | 1.185 | 0.236 |
| 20 | rms::rcs(visit, visit_df)visit’ | -0.833 | 0.550 | -1.513 | 0.130 |
| 20 | dose_lag1_f20 | 4.034 | 0.718 | 5.622 | 0.000 |
| 20 | dose_lag1_f30 | 3.789 | 0.849 | 4.462 | 0.000 |
| 20 | dose_lag1_f40 | 3.430 | 0.852 | 4.025 | 0.000 |
| 20 | dose_lag1_f50 | 2.646 | 0.816 | 3.243 | 0.001 |
| 20 | outcome_0 | 0.067 | 0.047 | 1.408 | 0.159 |
| 30 | (Intercept) | -2.149 | 1.483 | -1.449 | 0.147 |
| 30 | rms::rcs(visit, visit_df)visit | 0.552 | 0.342 | 1.612 | 0.107 |
| 30 | rms::rcs(visit, visit_df)visit’ | -1.173 | 0.560 | -2.095 | 0.036 |
| 30 | dose_lag1_f20 | 2.452 | 0.650 | 3.773 | 0.000 |
| 30 | dose_lag1_f30 | 4.151 | 0.778 | 5.335 | 0.000 |
| 30 | dose_lag1_f40 | 2.400 | 0.799 | 3.004 | 0.003 |
| 30 | dose_lag1_f50 | 1.892 | 0.759 | 2.494 | 0.013 |
| 30 | outcome_0 | 0.036 | 0.048 | 0.749 | 0.454 |
| 40 | (Intercept) | -3.093 | 1.563 | -1.979 | 0.048 |
| 40 | rms::rcs(visit, visit_df)visit | 0.197 | 0.350 | 0.563 | 0.573 |
| 40 | rms::rcs(visit, visit_df)visit’ | -0.491 | 0.572 | -0.858 | 0.391 |
| 40 | dose_lag1_f20 | 2.388 | 0.734 | 3.252 | 0.001 |
| 40 | dose_lag1_f30 | 3.046 | 0.861 | 3.537 | 0.000 |
| 40 | dose_lag1_f40 | 4.106 | 0.843 | 4.872 | 0.000 |
| 40 | dose_lag1_f50 | 1.810 | 0.849 | 2.132 | 0.033 |
| 40 | outcome_0 | 0.075 | 0.050 | 1.519 | 0.129 |
| 50 | (Intercept) | -3.854 | 1.617 | -2.383 | 0.017 |
| 50 | rms::rcs(visit, visit_df)visit | 0.582 | 0.365 | 1.595 | 0.111 |
| 50 | rms::rcs(visit, visit_df)visit’ | -1.157 | 0.595 | -1.943 | 0.052 |
| 50 | dose_lag1_f20 | 1.948 | 0.711 | 2.738 | 0.006 |
| 50 | dose_lag1_f30 | 2.476 | 0.847 | 2.924 | 0.003 |
| 50 | dose_lag1_f40 | 1.898 | 0.870 | 2.182 | 0.029 |
| 50 | dose_lag1_f50 | 2.771 | 0.781 | 3.547 | 0.000 |
| 50 | outcome_0 | 0.083 | 0.051 | 1.628 | 0.104 |

<small><em>Multinomial dose-weight numerator model:
29060_003_PAROXETINE</em></small>

| arm_name | n | n_patients | mean_SW_treatment_multinom | sd_SW_treatment_multinom | min_SW_treatment_multinom | p1_SW_treatment_multinom | p50_SW_treatment_multinom | p99_SW_treatment_multinom | max_SW_treatment_multinom | mean_cSW_treatment_multinom | sd_cSW_treatment_multinom | min_cSW_treatment_multinom | p1_cSW_treatment_multinom | p50_cSW_treatment_multinom | p99_cSW_treatment_multinom | max_cSW_treatment_multinom | ESS_cSW_treatment_multinom |
|:---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 29060_002_PAROXETINE | 564 | 165 | 0.939 | 0.871 | 0.015 | 0.134 | 0.742 | 4.840 | 11.128 | 0.812 | 0.888 | 0.002 | 0.018 | 0.586 | 4.540 | 8.115 | 256.910 |
| 29060_003_IMIPRAMINE | 684 | 229 | 0.983 | 1.706 | 0.073 | 0.178 | 0.731 | 7.444 | 29.497 | 1.055 | 2.901 | 0.011 | 0.025 | 0.585 | 15.287 | 34.215 | 80.052 |
| 29060_003_PAROXETINE | 758 | 237 | 1.100 | 3.176 | 0.081 | 0.213 | 0.808 | 5.688 | 83.302 | 1.370 | 6.762 | 0.020 | 0.045 | 0.703 | 11.451 | 133.670 | 29.935 |

<small><em>Summary of multinomial stabilized treatment
weights</em></small>

| arm_name | n | n_patients | mean_SW_total_multinom_trunc | sd_SW_total_multinom_trunc | min_SW_total_multinom_trunc | p1_SW_total_multinom_trunc | p50_SW_total_multinom_trunc | p99_SW_total_multinom_trunc | max_SW_total_multinom_trunc | ESS_SW_total_multinom_trunc |
|:---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 29060_002_PAROXETINE | 525 | 163 | 0.802 | 0.814 | 0.019 | 0.020 | 0.591 | 4.806 | 4.828 | 258.888 |
| 29060_003_IMIPRAMINE | 682 | 229 | 0.949 | 1.920 | 0.024 | 0.024 | 0.599 | 14.464 | 15.067 | 134.129 |
| 29060_003_PAROXETINE | 756 | 237 | 1.004 | 1.512 | 0.050 | 0.051 | 0.713 | 10.997 | 11.738 | 231.635 |

<small><em>Summary of truncated total weights using multinomial
IPTW</em></small>

| Term                                   | Estimate | Std.err |   Wald | Pr(\>\|W\|) |
|:---------------------------------------|---------:|--------:|-------:|------------:|
| (Intercept)                            |  -28.421 |   8.664 | 10.760 |       0.001 |
| rms::rcs(visit, 3)visit                |    5.622 |   4.728 |  1.414 |       0.234 |
| rms::rcs(visit, 3)visit’               |   -3.455 |   8.180 |  0.178 |       0.673 |
| outcome_0                              |    0.776 |   0.252 |  9.467 |       0.002 |
| dose_lag1_f20                          |    1.226 |   4.966 |  0.061 |       0.805 |
| dose_lag1_f30                          |   11.914 |   5.055 |  5.554 |       0.018 |
| dose_lag1_f40                          |    7.678 |   4.872 |  2.483 |       0.115 |
| dose_lag1_f50                          |    3.597 |   4.952 |  0.528 |       0.468 |
| dose_lag2_f10                          |   -3.537 |   2.816 |  1.577 |       0.209 |
| dose_lag2_f20                          |   -2.366 |   1.639 |  2.085 |       0.149 |
| dose_lag2_f30                          |   -1.081 |   1.825 |  0.351 |       0.554 |
| dose_lag2_f40                          |   -0.610 |   2.238 |  0.074 |       0.785 |
| dose_lag2_f50                          |   -0.014 |   1.884 |  0.000 |       0.994 |
| dose_lag3_f10                          |   -4.300 |   2.202 |  3.812 |       0.051 |
| dose_lag3_f20                          |   -0.087 |   1.728 |  0.003 |       0.960 |
| dose_lag3_f30                          |    1.551 |   2.073 |  0.560 |       0.454 |
| dose_lag3_f40                          |   -1.481 |   2.165 |  0.468 |       0.494 |
| dose_lag3_f50                          |    0.757 |   2.210 |  0.117 |       0.732 |
| avg_dose_before_lag3                   |   -0.031 |   0.041 |  0.580 |       0.446 |
| rms::rcs(visit, 3)visit:outcome_0      |    0.116 |   0.134 |  0.744 |       0.388 |
| rms::rcs(visit, 3)visit’:outcome_0     |   -0.244 |   0.229 |  1.138 |       0.286 |
| rms::rcs(visit, 3)visit:dose_lag1_f20  |   -1.740 |   2.675 |  0.423 |       0.515 |
| rms::rcs(visit, 3)visit’:dose_lag1_f20 |    1.903 |   4.629 |  0.169 |       0.681 |
| rms::rcs(visit, 3)visit:dose_lag1_f30  |   -6.741 |   2.734 |  6.080 |       0.014 |
| rms::rcs(visit, 3)visit’:dose_lag1_f30 |    9.117 |   4.925 |  3.426 |       0.064 |
| rms::rcs(visit, 3)visit:dose_lag1_f40  |   -5.116 |   2.849 |  3.226 |       0.072 |
| rms::rcs(visit, 3)visit’:dose_lag1_f40 |    7.138 |   4.948 |  2.081 |       0.149 |
| rms::rcs(visit, 3)visit:dose_lag1_f50  |   -3.635 |   2.718 |  1.788 |       0.181 |
| rms::rcs(visit, 3)visit’:dose_lag1_f50 |    3.900 |   4.747 |  0.675 |       0.411 |

<small><em>Weighted MSM using multinomial IPTW:
29060_002_PAROXETINE</em></small>

| Term                                    | Estimate | Std.err |   Wald | Pr(\>\|W\|) |
|:----------------------------------------|---------:|--------:|-------:|------------:|
| (Intercept)                             |  -25.017 |   7.816 | 10.245 |       0.001 |
| rms::rcs(visit, 3)visit                 |    3.590 |   3.982 |  0.813 |       0.367 |
| rms::rcs(visit, 3)visit’                |   -0.853 |   5.291 |  0.026 |       0.872 |
| outcome_0                               |    1.079 |   0.277 | 15.204 |       0.000 |
| dose_lag1_f65                           |   -0.010 |   4.894 |  0.000 |       0.998 |
| dose_lag1_f80                           |    5.138 |   7.348 |  0.489 |       0.484 |
| dose_lag1_f145                          |   -4.779 |   2.844 |  2.823 |       0.093 |
| dose_lag1_f210                          |   -2.595 |   2.982 |  0.757 |       0.384 |
| dose_lag1_f275                          |    0.658 |   4.102 |  0.026 |       0.873 |
| dose_lag2_f20                           |    2.084 |   1.115 |  3.490 |       0.062 |
| dose_lag2_f65                           |   -2.218 |   2.667 |  0.692 |       0.405 |
| dose_lag2_f80                           |    2.002 |   1.813 |  1.219 |       0.270 |
| dose_lag2_f145                          |   -0.449 |   1.098 |  0.167 |       0.682 |
| dose_lag2_f210                          |   -2.203 |   1.389 |  2.516 |       0.113 |
| dose_lag2_f275                          |   -2.327 |   1.329 |  3.068 |       0.080 |
| dose_lag3_f20                           |    0.055 |   2.122 |  0.001 |       0.979 |
| dose_lag3_f65                           |   -0.547 |   2.449 |  0.050 |       0.823 |
| dose_lag3_f80                           |    1.282 |   2.424 |  0.280 |       0.597 |
| dose_lag3_f145                          |    1.056 |   1.384 |  0.582 |       0.445 |
| dose_lag3_f210                          |    1.957 |   1.337 |  2.141 |       0.143 |
| dose_lag3_f275                          |    2.418 |   1.843 |  1.721 |       0.190 |
| avg_dose_before_lag3                    |    0.000 |   0.007 |  0.001 |       0.975 |
| rms::rcs(visit, 3)visit:outcome_0       |   -0.086 |   0.138 |  0.388 |       0.534 |
| rms::rcs(visit, 3)visit’:outcome_0      |    0.072 |   0.179 |  0.160 |       0.689 |
| rms::rcs(visit, 3)visit:dose_lag1_f65   |    1.830 |   2.898 |  0.399 |       0.528 |
| rms::rcs(visit, 3)visit’:dose_lag1_f65  |   -2.978 |   4.068 |  0.536 |       0.464 |
| rms::rcs(visit, 3)visit:dose_lag1_f80   |   -0.648 |   3.212 |  0.041 |       0.840 |
| rms::rcs(visit, 3)visit’:dose_lag1_f80  |   -0.934 |   4.072 |  0.053 |       0.819 |
| rms::rcs(visit, 3)visit:dose_lag1_f145  |    2.994 |   1.398 |  4.586 |       0.032 |
| rms::rcs(visit, 3)visit’:dose_lag1_f145 |   -3.658 |   1.978 |  3.422 |       0.064 |
| rms::rcs(visit, 3)visit:dose_lag1_f210  |    1.884 |   1.465 |  1.654 |       0.198 |
| rms::rcs(visit, 3)visit’:dose_lag1_f210 |   -3.188 |   1.807 |  3.112 |       0.078 |
| rms::rcs(visit, 3)visit:dose_lag1_f275  |    0.218 |   1.886 |  0.013 |       0.908 |
| rms::rcs(visit, 3)visit’:dose_lag1_f275 |   -0.986 |   2.724 |  0.131 |       0.717 |

<small><em>Weighted MSM using multinomial IPTW:
29060_003_IMIPRAMINE</em></small>

| Term                                   | Estimate | Std.err |   Wald | Pr(\>\|W\|) |
|:---------------------------------------|---------:|--------:|-------:|------------:|
| (Intercept)                            |  -27.439 |   9.806 |  7.830 |       0.005 |
| rms::rcs(visit, 3)visit                |    4.838 |   4.257 |  1.292 |       0.256 |
| rms::rcs(visit, 3)visit’               |   -6.541 |   7.181 |  0.830 |       0.362 |
| outcome_0                              |    1.074 |   0.288 | 13.866 |       0.000 |
| dose_lag1_f20                          |    0.222 |   6.643 |  0.001 |       0.973 |
| dose_lag1_f30                          |    0.196 |   6.702 |  0.001 |       0.977 |
| dose_lag1_f40                          |    0.103 |   7.529 |  0.000 |       0.989 |
| dose_lag1_f50                          |   -4.256 |   7.301 |  0.340 |       0.560 |
| dose_lag2_f10                          |   -0.835 |   2.126 |  0.154 |       0.694 |
| dose_lag2_f20                          |    0.580 |   1.449 |  0.160 |       0.689 |
| dose_lag2_f30                          |    0.328 |   1.305 |  0.063 |       0.802 |
| dose_lag2_f40                          |    1.551 |   1.470 |  1.113 |       0.291 |
| dose_lag2_f50                          |    1.634 |   1.412 |  1.340 |       0.247 |
| dose_lag3_f10                          |   -0.120 |   1.948 |  0.004 |       0.951 |
| dose_lag3_f20                          |   -0.691 |   1.185 |  0.340 |       0.560 |
| dose_lag3_f30                          |   -0.826 |   1.072 |  0.594 |       0.441 |
| dose_lag3_f40                          |   -0.015 |   1.542 |  0.000 |       0.992 |
| dose_lag3_f50                          |   -2.474 |   1.673 |  2.186 |       0.139 |
| avg_dose_before_lag3                   |    0.041 |   0.036 |  1.316 |       0.251 |
| rms::rcs(visit, 3)visit:outcome_0      |   -0.073 |   0.116 |  0.392 |       0.531 |
| rms::rcs(visit, 3)visit’:outcome_0     |    0.109 |   0.173 |  0.397 |       0.529 |
| rms::rcs(visit, 3)visit:dose_lag1_f20  |   -0.428 |   3.117 |  0.019 |       0.891 |
| rms::rcs(visit, 3)visit’:dose_lag1_f20 |    3.589 |   5.918 |  0.368 |       0.544 |
| rms::rcs(visit, 3)visit:dose_lag1_f30  |   -0.472 |   3.130 |  0.023 |       0.880 |
| rms::rcs(visit, 3)visit’:dose_lag1_f30 |    2.034 |   5.966 |  0.116 |       0.733 |
| rms::rcs(visit, 3)visit:dose_lag1_f40  |   -0.182 |   3.442 |  0.003 |       0.958 |
| rms::rcs(visit, 3)visit’:dose_lag1_f40 |    1.286 |   6.354 |  0.041 |       0.840 |
| rms::rcs(visit, 3)visit:dose_lag1_f50  |    1.889 |   3.332 |  0.321 |       0.571 |
| rms::rcs(visit, 3)visit’:dose_lag1_f50 |   -1.802 |   6.269 |  0.083 |       0.774 |

<small><em>Weighted MSM using multinomial IPTW:
29060_003_PAROXETINE</em></small>

### Predicted trajectories using multinomial IPTW

![](README_files/figure-gfm/appendix-multinomial-dose-model-profile-grid-1.png)<!-- -->

![](README_files/figure-gfm/appendix-multinomial-prediction-grid-1.png)<!-- -->

##### Comparison of ordinal and multinomial IPTW results

| weight_model | arm_name | trial_name | week | dose_strategy | strategy_dose | active_n_patients | active_n_visit_rows | placebo_n_patients | placebo_n_visit_rows | predicted_active_improvement | SE_active_prediction | observed_placebo_improvement | SE_observed_placebo | dose_vs_placebo_difference | SE | lower_95 | upper_95 |
|:---|:---|:---|---:|:---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| Ordinal IPTW | 29060_002_PAROXETINE | 29060_002 | 1 | Target 10 mg | 10 | 86 | 86 | 73 | 73 | 2.091 | 2.791 | 3.753 | 0.839 | -1.662 | 2.915 | -7.375 | 4.051 |
| Ordinal IPTW | 29060_002_PAROXETINE | 29060_002 | 2 | Target 10 mg | 10 | 91 | 91 | 99 | 99 | 5.231 | 2.126 | 4.322 | 0.827 | 0.909 | 2.281 | -3.561 | 5.380 |
| Ordinal IPTW | 29060_002_PAROXETINE | 29060_002 | 3 | Target 10 mg | 10 | 89 | 89 | 89 | 89 | 6.881 | 3.192 | 7.746 | 0.933 | -0.865 | 3.326 | -7.384 | 5.653 |
| Ordinal IPTW | 29060_002_PAROXETINE | 29060_002 | 4 | Target 10 mg | 10 | 81 | 81 | 86 | 86 | 10.936 | 3.068 | 9.505 | 1.155 | 1.431 | 3.278 | -4.993 | 7.856 |
| Ordinal IPTW | 29060_002_PAROXETINE | 29060_002 | 5 | Target 10 mg | 10 | 50 | 50 | 37 | 37 | 12.562 | 2.622 | 9.757 | 1.674 | 2.805 | 3.111 | -3.292 | 8.902 |
| Ordinal IPTW | 29060_002_PAROXETINE | 29060_002 | 6 | Target 10 mg | 10 | 55 | 55 | 49 | 49 | 12.569 | 2.248 | 10.814 | 1.163 | 1.754 | 2.531 | -3.206 | 6.715 |
| Ordinal IPTW | 29060_002_PAROXETINE | 29060_002 | 7 | Target 10 mg | 10 | 62 | 62 | 45 | 45 | 11.765 | 2.638 | 12.941 | 1.436 | -1.176 | 3.003 | -7.063 | 4.710 |
| Ordinal IPTW | 29060_002_PAROXETINE | 29060_002 | 8 | Target 10 mg | 10 | 11 | 11 | 6 | 6 | 10.827 | 3.697 | 19.861 | 4.748 | -9.035 | 6.017 | -20.828 | 2.758 |
| Ordinal IPTW | 29060_002_PAROXETINE | 29060_002 | 1 | Target 20 mg | 20 | 86 | 86 | 73 | 73 | 1.365 | 1.561 | 3.753 | 0.839 | -2.388 | 1.772 | -5.862 | 1.085 |
| Ordinal IPTW | 29060_002_PAROXETINE | 29060_002 | 2 | Target 20 mg | 20 | 91 | 91 | 99 | 99 | 5.151 | 1.409 | 4.322 | 0.827 | 0.830 | 1.633 | -2.372 | 4.031 |
| Ordinal IPTW | 29060_002_PAROXETINE | 29060_002 | 3 | Target 20 mg | 20 | 89 | 89 | 89 | 89 | 9.388 | 1.649 | 7.746 | 0.933 | 1.642 | 1.895 | -2.072 | 5.356 |
| Ordinal IPTW | 29060_002_PAROXETINE | 29060_002 | 4 | Target 20 mg | 20 | 81 | 81 | 86 | 86 | 12.415 | 1.506 | 9.505 | 1.155 | 2.910 | 1.898 | -0.810 | 6.630 |
| Ordinal IPTW | 29060_002_PAROXETINE | 29060_002 | 5 | Target 20 mg | 20 | 50 | 50 | 37 | 37 | 13.775 | 1.318 | 9.757 | 1.674 | 4.018 | 2.130 | -0.158 | 8.193 |
| Ordinal IPTW | 29060_002_PAROXETINE | 29060_002 | 6 | Target 20 mg | 20 | 55 | 55 | 49 | 49 | 14.025 | 1.264 | 10.814 | 1.163 | 3.211 | 1.718 | -0.156 | 6.577 |
| Ordinal IPTW | 29060_002_PAROXETINE | 29060_002 | 7 | Target 20 mg | 20 | 62 | 62 | 45 | 45 | 13.719 | 1.604 | 12.941 | 1.436 | 0.778 | 2.153 | -3.441 | 4.997 |
| Ordinal IPTW | 29060_002_PAROXETINE | 29060_002 | 8 | Target 20 mg | 20 | 11 | 11 | 6 | 6 | 13.321 | 2.222 | 19.861 | 4.748 | -6.540 | 5.242 | -16.814 | 3.733 |
| Ordinal IPTW | 29060_002_PAROXETINE | 29060_002 | 1 | Target 30 mg | 30 | 86 | 86 | 73 | 73 | 5.586 | 1.403 | 3.753 | 0.839 | 1.833 | 1.635 | -1.371 | 5.037 |
| Ordinal IPTW | 29060_002_PAROXETINE | 29060_002 | 2 | Target 30 mg | 30 | 91 | 91 | 99 | 99 | 7.872 | 1.428 | 4.322 | 0.827 | 3.551 | 1.650 | 0.317 | 6.785 |
| Ordinal IPTW | 29060_002_PAROXETINE | 29060_002 | 3 | Target 30 mg | 30 | 89 | 89 | 89 | 89 | 9.481 | 1.790 | 7.746 | 0.933 | 1.735 | 2.019 | -2.222 | 5.692 |
| Ordinal IPTW | 29060_002_PAROXETINE | 29060_002 | 4 | Target 30 mg | 30 | 81 | 81 | 86 | 86 | 11.466 | 1.726 | 9.505 | 1.155 | 1.962 | 2.076 | -2.108 | 6.031 |
| Ordinal IPTW | 29060_002_PAROXETINE | 29060_002 | 5 | Target 30 mg | 30 | 50 | 50 | 37 | 37 | 13.120 | 1.839 | 9.757 | 1.674 | 3.363 | 2.487 | -1.511 | 8.236 |
| Ordinal IPTW | 29060_002_PAROXETINE | 29060_002 | 6 | Target 30 mg | 30 | 55 | 55 | 49 | 49 | 14.552 | 2.433 | 10.814 | 1.163 | 3.738 | 2.696 | -1.547 | 9.023 |
| Ordinal IPTW | 29060_002_PAROXETINE | 29060_002 | 7 | Target 30 mg | 30 | 62 | 62 | 45 | 45 | 15.874 | 3.456 | 12.941 | 1.436 | 2.933 | 3.742 | -4.401 | 10.267 |
| Ordinal IPTW | 29060_002_PAROXETINE | 29060_002 | 8 | Target 30 mg | 30 | 11 | 11 | 6 | 6 | 17.178 | 4.657 | 19.861 | 4.748 | -2.684 | 6.650 | -15.718 | 10.351 |
| Ordinal IPTW | 29060_002_PAROXETINE | 29060_002 | 1 | Target 40 mg | 40 | 86 | 86 | 73 | 73 | 3.868 | 1.189 | 3.753 | 0.839 | 0.115 | 1.455 | -2.738 | 2.967 |
| Ordinal IPTW | 29060_002_PAROXETINE | 29060_002 | 2 | Target 40 mg | 40 | 91 | 91 | 99 | 99 | 5.308 | 1.720 | 4.322 | 0.827 | 0.986 | 1.908 | -2.753 | 4.726 |
| Ordinal IPTW | 29060_002_PAROXETINE | 29060_002 | 3 | Target 40 mg | 40 | 89 | 89 | 89 | 89 | 4.746 | 1.938 | 7.746 | 0.933 | -3.000 | 2.151 | -7.215 | 1.216 |
| Ordinal IPTW | 29060_002_PAROXETINE | 29060_002 | 4 | Target 40 mg | 40 | 81 | 81 | 86 | 86 | 6.813 | 1.798 | 9.505 | 1.155 | -2.692 | 2.137 | -6.880 | 1.496 |
| Ordinal IPTW | 29060_002_PAROXETINE | 29060_002 | 5 | Target 40 mg | 40 | 50 | 50 | 37 | 37 | 8.736 | 1.465 | 9.757 | 1.674 | -1.021 | 2.224 | -5.380 | 3.338 |
| Ordinal IPTW | 29060_002_PAROXETINE | 29060_002 | 6 | Target 40 mg | 40 | 55 | 55 | 49 | 49 | 10.563 | 1.200 | 10.814 | 1.163 | -0.251 | 1.671 | -3.526 | 3.025 |
| Ordinal IPTW | 29060_002_PAROXETINE | 29060_002 | 7 | Target 40 mg | 40 | 62 | 62 | 45 | 45 | 12.343 | 1.499 | 12.941 | 1.436 | -0.599 | 2.075 | -4.666 | 3.469 |
| Ordinal IPTW | 29060_002_PAROXETINE | 29060_002 | 8 | Target 40 mg | 40 | 11 | 11 | 6 | 6 | 14.114 | 2.233 | 19.861 | 4.748 | -5.747 | 5.246 | -16.030 | 4.535 |
| Ordinal IPTW | 29060_002_PAROXETINE | 29060_002 | 1 | Target 50 mg | 50 | 86 | 86 | 73 | 73 | 0.879 | 1.427 | 3.753 | 0.839 | -2.874 | 1.655 | -6.118 | 0.370 |
| Ordinal IPTW | 29060_002_PAROXETINE | 29060_002 | 2 | Target 50 mg | 50 | 91 | 91 | 99 | 99 | 6.137 | 1.991 | 4.322 | 0.827 | 1.816 | 2.156 | -2.410 | 6.042 |
| Ordinal IPTW | 29060_002_PAROXETINE | 29060_002 | 3 | Target 50 mg | 50 | 89 | 89 | 89 | 89 | 7.472 | 2.005 | 7.746 | 0.933 | -0.274 | 2.211 | -4.608 | 4.061 |
| Ordinal IPTW | 29060_002_PAROXETINE | 29060_002 | 4 | Target 50 mg | 50 | 81 | 81 | 86 | 86 | 9.828 | 1.933 | 9.505 | 1.155 | 0.323 | 2.252 | -4.090 | 4.737 |
| Ordinal IPTW | 29060_002_PAROXETINE | 29060_002 | 5 | Target 50 mg | 50 | 50 | 50 | 37 | 37 | 10.586 | 1.778 | 9.757 | 1.674 | 0.829 | 2.442 | -3.957 | 5.615 |
| Ordinal IPTW | 29060_002_PAROXETINE | 29060_002 | 6 | Target 50 mg | 50 | 55 | 55 | 49 | 49 | 10.279 | 1.670 | 10.814 | 1.163 | -0.535 | 2.035 | -4.524 | 3.454 |
| Ordinal IPTW | 29060_002_PAROXETINE | 29060_002 | 7 | Target 50 mg | 50 | 62 | 62 | 45 | 45 | 9.440 | 1.845 | 12.941 | 1.436 | -3.501 | 2.338 | -8.084 | 1.081 |
| Ordinal IPTW | 29060_002_PAROXETINE | 29060_002 | 8 | Target 50 mg | 50 | 11 | 11 | 6 | 6 | 8.512 | 2.306 | 19.861 | 4.748 | -11.349 | 5.278 | -21.694 | -1.005 |
| Ordinal IPTW | 29060_003_IMIPRAMINE | 29060_003 | 1 | Target 20 mg | 20 | 120 | 120 | 125 | 125 | 4.922 | 1.969 | 3.113 | 0.645 | 1.808 | 2.072 | -2.252 | 5.869 |
| Ordinal IPTW | 29060_003_IMIPRAMINE | 29060_003 | 2 | Target 20 mg | 20 | 144 | 144 | 144 | 144 | 7.686 | 2.107 | 4.042 | 0.681 | 3.644 | 2.214 | -0.696 | 7.984 |
| Ordinal IPTW | 29060_003_IMIPRAMINE | 29060_003 | 3 | Target 20 mg | 20 | 132 | 132 | 144 | 144 | 10.118 | 1.506 | 5.359 | 0.671 | 4.759 | 1.648 | 1.528 | 7.990 |
| Ordinal IPTW | 29060_003_IMIPRAMINE | 29060_003 | 4 | Target 20 mg | 20 | 108 | 108 | 134 | 134 | 11.980 | 1.443 | 7.574 | 0.824 | 4.405 | 1.662 | 1.149 | 7.662 |
| Ordinal IPTW | 29060_003_IMIPRAMINE | 29060_003 | 5 | Target 20 mg | 20 | 53 | 53 | 58 | 58 | 14.393 | 1.379 | 8.501 | 1.149 | 5.892 | 1.795 | 2.374 | 9.410 |
| Ordinal IPTW | 29060_003_IMIPRAMINE | 29060_003 | 6 | Target 20 mg | 20 | 73 | 73 | 68 | 68 | 17.083 | 1.855 | 10.208 | 1.163 | 6.875 | 2.190 | 2.583 | 11.167 |
| Ordinal IPTW | 29060_003_IMIPRAMINE | 29060_003 | 7 | Target 20 mg | 20 | 50 | 50 | 82 | 82 | 19.819 | 2.730 | 13.486 | 1.137 | 6.333 | 2.958 | 0.536 | 12.130 |
| Ordinal IPTW | 29060_003_IMIPRAMINE | 29060_003 | 8 | Target 20 mg | 20 | 2 | 2 | 1 | 1 | 22.555 | 3.733 | 12.693 | NA | 9.862 | NA | NA | NA |
| Ordinal IPTW | 29060_003_IMIPRAMINE | 29060_003 | 1 | Target 65 mg | 65 | 120 | 120 | 125 | 125 | 6.121 | 2.102 | 3.113 | 0.645 | 3.008 | 2.199 | -1.302 | 7.317 |
| Ordinal IPTW | 29060_003_IMIPRAMINE | 29060_003 | 2 | Target 65 mg | 65 | 144 | 144 | 144 | 144 | 8.566 | 1.803 | 4.042 | 0.681 | 4.524 | 1.927 | 0.747 | 8.301 |
| Ordinal IPTW | 29060_003_IMIPRAMINE | 29060_003 | 3 | Target 65 mg | 65 | 132 | 132 | 144 | 144 | 9.660 | 2.635 | 5.359 | 0.671 | 4.301 | 2.719 | -1.029 | 9.630 |
| Ordinal IPTW | 29060_003_IMIPRAMINE | 29060_003 | 4 | Target 65 mg | 65 | 108 | 108 | 134 | 134 | 11.154 | 2.676 | 7.574 | 0.824 | 3.580 | 2.800 | -1.908 | 9.068 |
| Ordinal IPTW | 29060_003_IMIPRAMINE | 29060_003 | 5 | Target 65 mg | 65 | 53 | 53 | 58 | 58 | 11.242 | 2.451 | 8.501 | 1.149 | 2.741 | 2.707 | -2.565 | 8.046 |
| Ordinal IPTW | 29060_003_IMIPRAMINE | 29060_003 | 6 | Target 65 mg | 65 | 73 | 73 | 68 | 68 | 10.628 | 3.111 | 10.208 | 1.163 | 0.420 | 3.321 | -6.090 | 6.929 |
| Ordinal IPTW | 29060_003_IMIPRAMINE | 29060_003 | 7 | Target 65 mg | 65 | 50 | 50 | 82 | 82 | 9.896 | 4.574 | 13.486 | 1.137 | -3.590 | 4.713 | -12.828 | 5.648 |
| Ordinal IPTW | 29060_003_IMIPRAMINE | 29060_003 | 8 | Target 65 mg | 65 | 2 | 2 | 1 | 1 | 9.165 | 6.314 | 12.693 | NA | -3.528 | NA | NA | NA |
| Ordinal IPTW | 29060_003_IMIPRAMINE | 29060_003 | 1 | Target 80 mg | 80 | 120 | 120 | 125 | 125 | 7.145 | 5.084 | 3.113 | 0.645 | 4.032 | 5.125 | -6.013 | 14.077 |
| Ordinal IPTW | 29060_003_IMIPRAMINE | 29060_003 | 2 | Target 80 mg | 80 | 144 | 144 | 144 | 144 | 9.185 | 2.445 | 4.042 | 0.681 | 5.143 | 2.538 | 0.169 | 10.118 |
| Ordinal IPTW | 29060_003_IMIPRAMINE | 29060_003 | 3 | Target 80 mg | 80 | 132 | 132 | 144 | 144 | 12.242 | 3.107 | 5.359 | 0.671 | 6.883 | 3.179 | 0.653 | 13.114 |
| Ordinal IPTW | 29060_003_IMIPRAMINE | 29060_003 | 4 | Target 80 mg | 80 | 108 | 108 | 134 | 134 | 13.753 | 3.461 | 7.574 | 0.824 | 6.179 | 3.558 | -0.794 | 13.152 |
| Ordinal IPTW | 29060_003_IMIPRAMINE | 29060_003 | 5 | Target 80 mg | 80 | 53 | 53 | 58 | 58 | 15.025 | 3.287 | 8.501 | 1.149 | 6.524 | 3.482 | -0.301 | 13.349 |
| Ordinal IPTW | 29060_003_IMIPRAMINE | 29060_003 | 6 | Target 80 mg | 80 | 73 | 73 | 68 | 68 | 16.178 | 3.753 | 10.208 | 1.163 | 5.970 | 3.929 | -1.732 | 13.671 |
| Ordinal IPTW | 29060_003_IMIPRAMINE | 29060_003 | 7 | Target 80 mg | 80 | 50 | 50 | 82 | 82 | 17.310 | 5.156 | 13.486 | 1.137 | 3.824 | 5.280 | -6.524 | 14.173 |
| Ordinal IPTW | 29060_003_IMIPRAMINE | 29060_003 | 8 | Target 80 mg | 80 | 2 | 2 | 1 | 1 | 18.443 | 6.975 | 12.693 | NA | 5.750 | NA | NA | NA |
| Ordinal IPTW | 29060_003_IMIPRAMINE | 29060_003 | 1 | Target 145 mg | 145 | 120 | 120 | 125 | 125 | 2.548 | 0.943 | 3.113 | 0.645 | -0.565 | 1.142 | -2.804 | 1.674 |
| Ordinal IPTW | 29060_003_IMIPRAMINE | 29060_003 | 2 | Target 145 mg | 145 | 144 | 144 | 144 | 144 | 6.369 | 0.858 | 4.042 | 0.681 | 2.327 | 1.095 | 0.181 | 4.473 |
| Ordinal IPTW | 29060_003_IMIPRAMINE | 29060_003 | 3 | Target 145 mg | 145 | 132 | 132 | 144 | 144 | 12.020 | 1.536 | 5.359 | 0.671 | 6.661 | 1.677 | 3.375 | 9.947 |
| Ordinal IPTW | 29060_003_IMIPRAMINE | 29060_003 | 4 | Target 145 mg | 145 | 108 | 108 | 134 | 134 | 14.301 | 1.387 | 7.574 | 0.824 | 6.727 | 1.613 | 3.566 | 9.889 |
| Ordinal IPTW | 29060_003_IMIPRAMINE | 29060_003 | 5 | Target 145 mg | 145 | 53 | 53 | 58 | 58 | 15.673 | 1.365 | 8.501 | 1.149 | 7.171 | 1.784 | 3.674 | 10.669 |
| Ordinal IPTW | 29060_003_IMIPRAMINE | 29060_003 | 6 | Target 145 mg | 145 | 73 | 73 | 68 | 68 | 16.589 | 1.767 | 10.208 | 1.163 | 6.381 | 2.116 | 2.234 | 10.527 |
| Ordinal IPTW | 29060_003_IMIPRAMINE | 29060_003 | 7 | Target 145 mg | 145 | 50 | 50 | 82 | 82 | 17.429 | 2.465 | 13.486 | 1.137 | 3.943 | 2.715 | -1.378 | 9.263 |
| Ordinal IPTW | 29060_003_IMIPRAMINE | 29060_003 | 8 | Target 145 mg | 145 | 2 | 2 | 1 | 1 | 18.269 | 3.274 | 12.693 | NA | 5.576 | NA | NA | NA |
| Ordinal IPTW | 29060_003_IMIPRAMINE | 29060_003 | 1 | Target 210 mg | 210 | 120 | 120 | 125 | 125 | 4.307 | 0.974 | 3.113 | 0.645 | 1.193 | 1.168 | -1.096 | 3.483 |
| Ordinal IPTW | 29060_003_IMIPRAMINE | 29060_003 | 2 | Target 210 mg | 210 | 144 | 144 | 144 | 144 | 5.879 | 1.287 | 4.042 | 0.681 | 1.837 | 1.456 | -1.016 | 4.690 |
| Ordinal IPTW | 29060_003_IMIPRAMINE | 29060_003 | 3 | Target 210 mg | 210 | 132 | 132 | 144 | 144 | 9.800 | 1.668 | 5.359 | 0.671 | 4.441 | 1.798 | 0.918 | 7.964 |
| Ordinal IPTW | 29060_003_IMIPRAMINE | 29060_003 | 4 | Target 210 mg | 210 | 108 | 108 | 134 | 134 | 11.430 | 1.503 | 7.574 | 0.824 | 3.856 | 1.714 | 0.496 | 7.216 |
| Ordinal IPTW | 29060_003_IMIPRAMINE | 29060_003 | 5 | Target 210 mg | 210 | 53 | 53 | 58 | 58 | 12.509 | 1.328 | 8.501 | 1.149 | 4.008 | 1.756 | 0.565 | 7.450 |
| Ordinal IPTW | 29060_003_IMIPRAMINE | 29060_003 | 6 | Target 210 mg | 210 | 73 | 73 | 68 | 68 | 13.313 | 1.243 | 10.208 | 1.163 | 3.105 | 1.702 | -0.232 | 6.441 |
| Ordinal IPTW | 29060_003_IMIPRAMINE | 29060_003 | 7 | Target 210 mg | 210 | 50 | 50 | 82 | 82 | 14.070 | 1.394 | 13.486 | 1.137 | 0.584 | 1.799 | -2.941 | 4.110 |
| Ordinal IPTW | 29060_003_IMIPRAMINE | 29060_003 | 8 | Target 210 mg | 210 | 2 | 2 | 1 | 1 | 14.828 | 1.733 | 12.693 | NA | 2.135 | NA | NA | NA |
| Ordinal IPTW | 29060_003_IMIPRAMINE | 29060_003 | 1 | Target 275 mg | 275 | 120 | 120 | 125 | 125 | 6.633 | 2.165 | 3.113 | 0.645 | 3.520 | 2.259 | -0.907 | 7.946 |
| Ordinal IPTW | 29060_003_IMIPRAMINE | 29060_003 | 2 | Target 275 mg | 275 | 144 | 144 | 144 | 144 | 7.471 | 2.628 | 4.042 | 0.681 | 3.429 | 2.715 | -1.892 | 8.750 |
| Ordinal IPTW | 29060_003_IMIPRAMINE | 29060_003 | 3 | Target 275 mg | 275 | 132 | 132 | 144 | 144 | 12.307 | 3.726 | 5.359 | 0.671 | 6.948 | 3.786 | -0.472 | 14.368 |
| Ordinal IPTW | 29060_003_IMIPRAMINE | 29060_003 | 4 | Target 275 mg | 275 | 108 | 108 | 134 | 134 | 13.723 | 3.603 | 7.574 | 0.824 | 6.149 | 3.696 | -1.094 | 13.392 |
| Ordinal IPTW | 29060_003_IMIPRAMINE | 29060_003 | 5 | Target 275 mg | 275 | 53 | 53 | 58 | 58 | 14.990 | 3.231 | 8.501 | 1.149 | 6.488 | 3.429 | -0.232 | 13.209 |
| Ordinal IPTW | 29060_003_IMIPRAMINE | 29060_003 | 6 | Target 275 mg | 275 | 73 | 73 | 68 | 68 | 16.181 | 3.200 | 10.208 | 1.163 | 5.973 | 3.404 | -0.700 | 12.646 |
| Ordinal IPTW | 29060_003_IMIPRAMINE | 29060_003 | 7 | Target 275 mg | 275 | 50 | 50 | 82 | 82 | 17.360 | 3.793 | 13.486 | 1.137 | 3.874 | 3.960 | -3.887 | 11.635 |
| Ordinal IPTW | 29060_003_IMIPRAMINE | 29060_003 | 8 | Target 275 mg | 275 | 2 | 2 | 1 | 1 | 18.539 | 4.800 | 12.693 | NA | 5.846 | NA | NA | NA |
| Ordinal IPTW | 29060_003_PAROXETINE | 29060_003 | 1 | Target 10 mg | 10 | 126 | 126 | 125 | 125 | 9.634 | 4.035 | 3.113 | 0.645 | 6.521 | 4.086 | -1.488 | 14.529 |
| Ordinal IPTW | 29060_003_PAROXETINE | 29060_003 | 2 | Target 10 mg | 10 | 153 | 153 | 144 | 144 | 7.487 | 2.897 | 4.042 | 0.681 | 3.445 | 2.976 | -2.387 | 9.277 |
| Ordinal IPTW | 29060_003_PAROXETINE | 29060_003 | 3 | Target 10 mg | 10 | 145 | 145 | 144 | 144 | 9.286 | 1.549 | 5.359 | 0.671 | 3.927 | 1.689 | 0.617 | 7.236 |
| Ordinal IPTW | 29060_003_PAROXETINE | 29060_003 | 4 | Target 10 mg | 10 | 115 | 115 | 134 | 134 | 9.727 | 1.221 | 7.574 | 0.824 | 2.153 | 1.473 | -0.734 | 5.040 |
| Ordinal IPTW | 29060_003_PAROXETINE | 29060_003 | 5 | Target 10 mg | 10 | 58 | 58 | 58 | 58 | 11.589 | 2.488 | 8.501 | 1.149 | 3.088 | 2.740 | -2.283 | 8.459 |
| Ordinal IPTW | 29060_003_PAROXETINE | 29060_003 | 6 | Target 10 mg | 10 | 79 | 79 | 68 | 68 | 14.399 | 4.902 | 10.208 | 1.163 | 4.191 | 5.038 | -5.683 | 14.066 |
| Ordinal IPTW | 29060_003_PAROXETINE | 29060_003 | 7 | Target 10 mg | 10 | 74 | 74 | 82 | 82 | 17.683 | 7.819 | 13.486 | 1.137 | 4.197 | 7.902 | -11.290 | 19.683 |
| Ordinal IPTW | 29060_003_PAROXETINE | 29060_003 | 8 | Target 10 mg | 10 | 6 | 6 | 1 | 1 | 21.045 | 10.849 | 12.693 | NA | 8.353 | NA | NA | NA |
| Ordinal IPTW | 29060_003_PAROXETINE | 29060_003 | 1 | Target 20 mg | 20 | 126 | 126 | 125 | 125 | 3.708 | 0.969 | 3.113 | 0.645 | 0.595 | 1.164 | -1.686 | 2.876 |
| Ordinal IPTW | 29060_003_PAROXETINE | 29060_003 | 2 | Target 20 mg | 20 | 153 | 153 | 144 | 144 | 6.245 | 1.052 | 4.042 | 0.681 | 2.203 | 1.253 | -0.252 | 4.658 |
| Ordinal IPTW | 29060_003_PAROXETINE | 29060_003 | 3 | Target 20 mg | 20 | 145 | 145 | 144 | 144 | 8.550 | 1.218 | 5.359 | 0.671 | 3.191 | 1.391 | 0.465 | 5.917 |
| Ordinal IPTW | 29060_003_PAROXETINE | 29060_003 | 4 | Target 20 mg | 20 | 115 | 115 | 134 | 134 | 10.966 | 1.151 | 7.574 | 0.824 | 3.392 | 1.416 | 0.617 | 6.167 |
| Ordinal IPTW | 29060_003_PAROXETINE | 29060_003 | 5 | Target 20 mg | 20 | 58 | 58 | 58 | 58 | 13.234 | 0.985 | 8.501 | 1.149 | 4.732 | 1.513 | 1.766 | 7.699 |
| Ordinal IPTW | 29060_003_PAROXETINE | 29060_003 | 6 | Target 20 mg | 20 | 79 | 79 | 68 | 68 | 15.402 | 1.002 | 10.208 | 1.163 | 5.194 | 1.535 | 2.185 | 8.204 |
| Ordinal IPTW | 29060_003_PAROXETINE | 29060_003 | 7 | Target 20 mg | 20 | 74 | 74 | 82 | 82 | 17.522 | 1.430 | 13.486 | 1.137 | 4.035 | 1.827 | 0.455 | 7.616 |
| Ordinal IPTW | 29060_003_PAROXETINE | 29060_003 | 8 | Target 20 mg | 20 | 6 | 6 | 1 | 1 | 19.633 | 2.067 | 12.693 | NA | 6.940 | NA | NA | NA |
| Ordinal IPTW | 29060_003_PAROXETINE | 29060_003 | 1 | Target 30 mg | 30 | 126 | 126 | 125 | 125 | 3.853 | 1.006 | 3.113 | 0.645 | 0.740 | 1.195 | -1.602 | 3.082 |
| Ordinal IPTW | 29060_003_PAROXETINE | 29060_003 | 2 | Target 30 mg | 30 | 153 | 153 | 144 | 144 | 6.552 | 0.833 | 4.042 | 0.681 | 2.510 | 1.076 | 0.402 | 4.618 |
| Ordinal IPTW | 29060_003_PAROXETINE | 29060_003 | 3 | Target 30 mg | 30 | 145 | 145 | 144 | 144 | 8.888 | 1.100 | 5.359 | 0.671 | 3.529 | 1.289 | 1.004 | 6.055 |
| Ordinal IPTW | 29060_003_PAROXETINE | 29060_003 | 4 | Target 30 mg | 30 | 115 | 115 | 134 | 134 | 10.568 | 1.063 | 7.574 | 0.824 | 2.994 | 1.345 | 0.359 | 5.630 |
| Ordinal IPTW | 29060_003_PAROXETINE | 29060_003 | 5 | Target 30 mg | 30 | 58 | 58 | 58 | 58 | 11.872 | 0.951 | 8.501 | 1.149 | 3.371 | 1.491 | 0.448 | 6.294 |
| Ordinal IPTW | 29060_003_PAROXETINE | 29060_003 | 6 | Target 30 mg | 30 | 79 | 79 | 68 | 68 | 12.926 | 0.911 | 10.208 | 1.163 | 2.717 | 1.478 | -0.179 | 5.614 |
| Ordinal IPTW | 29060_003_PAROXETINE | 29060_003 | 7 | Target 30 mg | 30 | 74 | 74 | 82 | 82 | 13.853 | 1.167 | 13.486 | 1.137 | 0.367 | 1.630 | -2.827 | 3.561 |
| Ordinal IPTW | 29060_003_PAROXETINE | 29060_003 | 8 | Target 30 mg | 30 | 6 | 6 | 1 | 1 | 14.760 | 1.630 | 12.693 | NA | 2.068 | NA | NA | NA |
| Ordinal IPTW | 29060_003_PAROXETINE | 29060_003 | 1 | Target 40 mg | 40 | 126 | 126 | 125 | 125 | 5.088 | 1.805 | 3.113 | 0.645 | 1.975 | 1.916 | -1.781 | 5.731 |
| Ordinal IPTW | 29060_003_PAROXETINE | 29060_003 | 2 | Target 40 mg | 40 | 153 | 153 | 144 | 144 | 8.559 | 1.311 | 4.042 | 0.681 | 4.517 | 1.477 | 1.622 | 7.411 |
| Ordinal IPTW | 29060_003_PAROXETINE | 29060_003 | 3 | Target 40 mg | 40 | 145 | 145 | 144 | 144 | 10.862 | 1.643 | 5.359 | 0.671 | 5.503 | 1.775 | 2.024 | 8.982 |
| Ordinal IPTW | 29060_003_PAROXETINE | 29060_003 | 4 | Target 40 mg | 40 | 115 | 115 | 134 | 134 | 12.207 | 1.564 | 7.574 | 0.824 | 4.633 | 1.768 | 1.168 | 8.099 |
| Ordinal IPTW | 29060_003_PAROXETINE | 29060_003 | 5 | Target 40 mg | 40 | 58 | 58 | 58 | 58 | 13.440 | 1.317 | 8.501 | 1.149 | 4.939 | 1.748 | 1.514 | 8.364 |
| Ordinal IPTW | 29060_003_PAROXETINE | 29060_003 | 6 | Target 40 mg | 40 | 79 | 79 | 68 | 68 | 14.598 | 1.137 | 10.208 | 1.163 | 4.390 | 1.627 | 1.202 | 7.579 |
| Ordinal IPTW | 29060_003_PAROXETINE | 29060_003 | 7 | Target 40 mg | 40 | 74 | 74 | 82 | 82 | 15.719 | 1.457 | 13.486 | 1.137 | 2.233 | 1.848 | -1.389 | 5.855 |
| Ordinal IPTW | 29060_003_PAROXETINE | 29060_003 | 8 | Target 40 mg | 40 | 6 | 6 | 1 | 1 | 16.833 | 2.139 | 12.693 | NA | 4.140 | NA | NA | NA |
| Ordinal IPTW | 29060_003_PAROXETINE | 29060_003 | 1 | Target 50 mg | 50 | 126 | 126 | 125 | 125 | 1.417 | 2.253 | 3.113 | 0.645 | -1.696 | 2.344 | -6.290 | 2.898 |
| Ordinal IPTW | 29060_003_PAROXETINE | 29060_003 | 2 | Target 50 mg | 50 | 153 | 153 | 144 | 144 | 7.583 | 1.249 | 4.042 | 0.681 | 3.541 | 1.422 | 0.753 | 6.329 |
| Ordinal IPTW | 29060_003_PAROXETINE | 29060_003 | 3 | Target 50 mg | 50 | 145 | 145 | 144 | 144 | 10.464 | 1.979 | 5.359 | 0.671 | 5.105 | 2.090 | 1.009 | 9.201 |
| Ordinal IPTW | 29060_003_PAROXETINE | 29060_003 | 4 | Target 50 mg | 50 | 115 | 115 | 134 | 134 | 12.656 | 1.806 | 7.574 | 0.824 | 5.082 | 1.985 | 1.191 | 8.973 |
| Ordinal IPTW | 29060_003_PAROXETINE | 29060_003 | 5 | Target 50 mg | 50 | 58 | 58 | 58 | 58 | 13.695 | 1.620 | 8.501 | 1.149 | 5.194 | 1.986 | 1.302 | 9.086 |
| Ordinal IPTW | 29060_003_PAROXETINE | 29060_003 | 6 | Target 50 mg | 50 | 79 | 79 | 68 | 68 | 13.965 | 1.584 | 10.208 | 1.163 | 3.757 | 1.965 | -0.096 | 7.609 |
| Ordinal IPTW | 29060_003_PAROXETINE | 29060_003 | 7 | Target 50 mg | 50 | 74 | 74 | 82 | 82 | 13.850 | 1.950 | 13.486 | 1.137 | 0.363 | 2.257 | -4.061 | 4.788 |
| Ordinal IPTW | 29060_003_PAROXETINE | 29060_003 | 8 | Target 50 mg | 50 | 6 | 6 | 1 | 1 | 13.671 | 2.603 | 12.693 | NA | 0.978 | NA | NA | NA |
| Multinomial IPTW | 29060_002_PAROXETINE | 29060_002 | 1 | Target 10 mg | 10 | 86 | 86 | 73 | 73 | 1.214 | 2.379 | 3.753 | 0.839 | -2.539 | 2.523 | -7.483 | 2.406 |
| Multinomial IPTW | 29060_002_PAROXETINE | 29060_002 | 2 | Target 10 mg | 10 | 91 | 91 | 99 | 99 | 6.141 | 2.255 | 4.322 | 0.827 | 1.819 | 2.401 | -2.888 | 6.526 |
| Multinomial IPTW | 29060_002_PAROXETINE | 29060_002 | 3 | Target 10 mg | 10 | 89 | 89 | 89 | 89 | 8.320 | 3.378 | 7.746 | 0.933 | 0.574 | 3.505 | -6.295 | 7.443 |
| Multinomial IPTW | 29060_002_PAROXETINE | 29060_002 | 4 | Target 10 mg | 10 | 81 | 81 | 86 | 86 | 12.189 | 3.271 | 9.505 | 1.155 | 2.684 | 3.469 | -4.115 | 9.484 |
| Multinomial IPTW | 29060_002_PAROXETINE | 29060_002 | 5 | Target 10 mg | 10 | 50 | 50 | 37 | 37 | 13.552 | 2.638 | 9.757 | 1.674 | 3.795 | 3.124 | -2.329 | 9.918 |
| Multinomial IPTW | 29060_002_PAROXETINE | 29060_002 | 6 | Target 10 mg | 10 | 55 | 55 | 49 | 49 | 13.244 | 2.291 | 10.814 | 1.163 | 2.430 | 2.569 | -2.605 | 7.465 |
| Multinomial IPTW | 29060_002_PAROXETINE | 29060_002 | 7 | Target 10 mg | 10 | 62 | 62 | 45 | 45 | 12.101 | 3.276 | 12.941 | 1.436 | -0.841 | 3.577 | -7.851 | 6.170 |
| Multinomial IPTW | 29060_002_PAROXETINE | 29060_002 | 8 | Target 10 mg | 10 | 11 | 11 | 6 | 6 | 10.818 | 5.007 | 19.861 | 4.748 | -9.043 | 6.900 | -22.567 | 4.481 |
| Multinomial IPTW | 29060_002_PAROXETINE | 29060_002 | 1 | Target 20 mg | 20 | 86 | 86 | 73 | 73 | 0.700 | 1.586 | 3.753 | 0.839 | -3.053 | 1.794 | -6.570 | 0.463 |
| Multinomial IPTW | 29060_002_PAROXETINE | 29060_002 | 2 | Target 20 mg | 20 | 91 | 91 | 99 | 99 | 5.111 | 1.547 | 4.322 | 0.827 | 0.789 | 1.754 | -2.648 | 4.226 |
| Multinomial IPTW | 29060_002_PAROXETINE | 29060_002 | 3 | Target 20 mg | 20 | 89 | 89 | 89 | 89 | 9.819 | 1.798 | 7.746 | 0.933 | 2.073 | 2.026 | -1.897 | 6.043 |
| Multinomial IPTW | 29060_002_PAROXETINE | 29060_002 | 4 | Target 20 mg | 20 | 81 | 81 | 86 | 86 | 12.873 | 1.686 | 9.505 | 1.155 | 3.368 | 2.043 | -0.637 | 7.374 |
| Multinomial IPTW | 29060_002_PAROXETINE | 29060_002 | 5 | Target 20 mg | 20 | 50 | 50 | 37 | 37 | 13.897 | 1.468 | 9.757 | 1.674 | 4.140 | 2.226 | -0.223 | 8.503 |
| Multinomial IPTW | 29060_002_PAROXETINE | 29060_002 | 6 | Target 20 mg | 20 | 55 | 55 | 49 | 49 | 13.567 | 1.423 | 10.814 | 1.163 | 2.753 | 1.838 | -0.849 | 6.356 |
| Multinomial IPTW | 29060_002_PAROXETINE | 29060_002 | 7 | Target 20 mg | 20 | 62 | 62 | 45 | 45 | 12.561 | 1.867 | 12.941 | 1.436 | -0.381 | 2.355 | -4.997 | 4.236 |
| Multinomial IPTW | 29060_002_PAROXETINE | 29060_002 | 8 | Target 20 mg | 20 | 11 | 11 | 6 | 6 | 11.441 | 2.626 | 19.861 | 4.748 | -8.420 | 5.425 | -19.054 | 2.214 |
| Multinomial IPTW | 29060_002_PAROXETINE | 29060_002 | 1 | Target 30 mg | 30 | 86 | 86 | 73 | 73 | 6.387 | 1.561 | 3.753 | 0.839 | 2.634 | 1.772 | -0.839 | 6.107 |
| Multinomial IPTW | 29060_002_PAROXETINE | 29060_002 | 2 | Target 30 mg | 30 | 91 | 91 | 99 | 99 | 7.282 | 1.458 | 4.322 | 0.827 | 2.960 | 1.676 | -0.325 | 6.245 |
| Multinomial IPTW | 29060_002_PAROXETINE | 29060_002 | 3 | Target 30 mg | 30 | 89 | 89 | 89 | 89 | 9.716 | 1.788 | 7.746 | 0.933 | 1.970 | 2.017 | -1.984 | 5.923 |
| Multinomial IPTW | 29060_002_PAROXETINE | 29060_002 | 4 | Target 30 mg | 30 | 81 | 81 | 86 | 86 | 11.275 | 1.686 | 9.505 | 1.155 | 1.771 | 2.044 | -2.235 | 5.777 |
| Multinomial IPTW | 29060_002_PAROXETINE | 29060_002 | 5 | Target 30 mg | 30 | 50 | 50 | 37 | 37 | 12.608 | 1.833 | 9.757 | 1.674 | 2.851 | 2.482 | -2.015 | 7.716 |
| Multinomial IPTW | 29060_002_PAROXETINE | 29060_002 | 6 | Target 30 mg | 30 | 55 | 55 | 49 | 49 | 13.789 | 2.569 | 10.814 | 1.163 | 2.975 | 2.820 | -2.553 | 8.503 |
| Multinomial IPTW | 29060_002_PAROXETINE | 29060_002 | 7 | Target 30 mg | 30 | 62 | 62 | 45 | 45 | 14.895 | 3.762 | 12.941 | 1.436 | 1.954 | 4.027 | -5.939 | 9.846 |
| Multinomial IPTW | 29060_002_PAROXETINE | 29060_002 | 8 | Target 30 mg | 30 | 11 | 11 | 6 | 6 | 15.988 | 5.123 | 19.861 | 4.748 | -3.873 | 6.985 | -17.563 | 9.816 |
| Multinomial IPTW | 29060_002_PAROXETINE | 29060_002 | 1 | Target 40 mg | 40 | 86 | 86 | 73 | 73 | 3.776 | 1.277 | 3.753 | 0.839 | 0.023 | 1.528 | -2.971 | 3.017 |
| Multinomial IPTW | 29060_002_PAROXETINE | 29060_002 | 2 | Target 40 mg | 40 | 91 | 91 | 99 | 99 | 6.711 | 1.739 | 4.322 | 0.827 | 2.389 | 1.926 | -1.385 | 6.163 |
| Multinomial IPTW | 29060_002_PAROXETINE | 29060_002 | 3 | Target 40 mg | 40 | 89 | 89 | 89 | 89 | 7.040 | 2.625 | 7.746 | 0.933 | -0.707 | 2.786 | -6.166 | 4.753 |
| Multinomial IPTW | 29060_002_PAROXETINE | 29060_002 | 4 | Target 40 mg | 40 | 81 | 81 | 86 | 86 | 9.263 | 2.571 | 9.505 | 1.155 | -0.242 | 2.819 | -5.767 | 5.282 |
| Multinomial IPTW | 29060_002_PAROXETINE | 29060_002 | 5 | Target 40 mg | 40 | 50 | 50 | 37 | 37 | 10.764 | 1.980 | 9.757 | 1.674 | 1.007 | 2.593 | -4.075 | 6.088 |
| Multinomial IPTW | 29060_002_PAROXETINE | 29060_002 | 6 | Target 40 mg | 40 | 55 | 55 | 49 | 49 | 11.785 | 1.333 | 10.814 | 1.163 | 0.970 | 1.769 | -2.496 | 4.437 |
| Multinomial IPTW | 29060_002_PAROXETINE | 29060_002 | 7 | Target 40 mg | 40 | 62 | 62 | 45 | 45 | 12.565 | 1.780 | 12.941 | 1.436 | -0.377 | 2.287 | -4.859 | 4.105 |
| Multinomial IPTW | 29060_002_PAROXETINE | 29060_002 | 8 | Target 40 mg | 40 | 11 | 11 | 6 | 6 | 13.305 | 3.028 | 19.861 | 4.748 | -6.557 | 5.631 | -17.594 | 4.480 |
| Multinomial IPTW | 29060_002_PAROXETINE | 29060_002 | 1 | Target 50 mg | 50 | 86 | 86 | 73 | 73 | 1.175 | 1.489 | 3.753 | 0.839 | -2.577 | 1.709 | -5.927 | 0.772 |
| Multinomial IPTW | 29060_002_PAROXETINE | 29060_002 | 2 | Target 50 mg | 50 | 91 | 91 | 99 | 99 | 6.098 | 1.843 | 4.322 | 0.827 | 1.777 | 2.020 | -2.183 | 5.736 |
| Multinomial IPTW | 29060_002_PAROXETINE | 29060_002 | 3 | Target 50 mg | 50 | 89 | 89 | 89 | 89 | 9.202 | 3.193 | 7.746 | 0.933 | 1.456 | 3.326 | -5.064 | 7.975 |
| Multinomial IPTW | 29060_002_PAROXETINE | 29060_002 | 4 | Target 50 mg | 50 | 81 | 81 | 86 | 86 | 11.332 | 2.956 | 9.505 | 1.155 | 1.827 | 3.173 | -4.393 | 8.047 |
| Multinomial IPTW | 29060_002_PAROXETINE | 29060_002 | 5 | Target 50 mg | 50 | 50 | 50 | 37 | 37 | 11.931 | 2.713 | 9.757 | 1.674 | 2.173 | 3.188 | -4.075 | 8.422 |
| Multinomial IPTW | 29060_002_PAROXETINE | 29060_002 | 6 | Target 50 mg | 50 | 55 | 55 | 49 | 49 | 11.509 | 2.511 | 10.814 | 1.163 | 0.695 | 2.768 | -4.729 | 6.119 |
| Multinomial IPTW | 29060_002_PAROXETINE | 29060_002 | 7 | Target 50 mg | 50 | 62 | 62 | 45 | 45 | 10.577 | 2.485 | 12.941 | 1.436 | -2.364 | 2.870 | -7.989 | 3.260 |
| Multinomial IPTW | 29060_002_PAROXETINE | 29060_002 | 8 | Target 50 mg | 50 | 11 | 11 | 6 | 6 | 9.560 | 2.683 | 19.861 | 4.748 | -10.302 | 5.453 | -20.989 | 0.386 |
| Multinomial IPTW | 29060_003_IMIPRAMINE | 29060_003 | 1 | Target 20 mg | 20 | 120 | 120 | 125 | 125 | 4.805 | 1.627 | 3.113 | 0.645 | 1.692 | 1.750 | -1.739 | 5.122 |
| Multinomial IPTW | 29060_003_IMIPRAMINE | 29060_003 | 2 | Target 20 mg | 20 | 144 | 144 | 144 | 144 | 8.250 | 1.786 | 4.042 | 0.681 | 4.208 | 1.911 | 0.462 | 7.955 |
| Multinomial IPTW | 29060_003_IMIPRAMINE | 29060_003 | 3 | Target 20 mg | 20 | 132 | 132 | 144 | 144 | 9.920 | 1.152 | 5.359 | 0.671 | 4.561 | 1.334 | 1.947 | 7.175 |
| Multinomial IPTW | 29060_003_IMIPRAMINE | 29060_003 | 4 | Target 20 mg | 20 | 108 | 108 | 134 | 134 | 11.958 | 1.100 | 7.574 | 0.824 | 4.384 | 1.374 | 1.692 | 7.077 |
| Multinomial IPTW | 29060_003_IMIPRAMINE | 29060_003 | 5 | Target 20 mg | 20 | 53 | 53 | 58 | 58 | 14.328 | 1.086 | 8.501 | 1.149 | 5.826 | 1.581 | 2.728 | 8.925 |
| Multinomial IPTW | 29060_003_IMIPRAMINE | 29060_003 | 6 | Target 20 mg | 20 | 73 | 73 | 68 | 68 | 16.863 | 1.583 | 10.208 | 1.163 | 6.655 | 1.964 | 2.805 | 10.505 |
| Multinomial IPTW | 29060_003_IMIPRAMINE | 29060_003 | 7 | Target 20 mg | 20 | 50 | 50 | 82 | 82 | 19.425 | 2.388 | 13.486 | 1.137 | 5.939 | 2.645 | 0.756 | 11.123 |
| Multinomial IPTW | 29060_003_IMIPRAMINE | 29060_003 | 8 | Target 20 mg | 20 | 2 | 2 | 1 | 1 | 21.988 | 3.279 | 12.693 | NA | 9.295 | NA | NA | NA |
| Multinomial IPTW | 29060_003_IMIPRAMINE | 29060_003 | 1 | Target 65 mg | 65 | 120 | 120 | 125 | 125 | 6.624 | 1.984 | 3.113 | 0.645 | 3.511 | 2.086 | -0.577 | 7.599 |
| Multinomial IPTW | 29060_003_IMIPRAMINE | 29060_003 | 2 | Target 65 mg | 65 | 144 | 144 | 144 | 144 | 7.478 | 2.502 | 4.042 | 0.681 | 3.436 | 2.593 | -1.647 | 8.518 |
| Multinomial IPTW | 29060_003_IMIPRAMINE | 29060_003 | 3 | Target 65 mg | 65 | 132 | 132 | 144 | 144 | 9.552 | 2.930 | 5.359 | 0.671 | 4.193 | 3.006 | -1.698 | 10.084 |
| Multinomial IPTW | 29060_003_IMIPRAMINE | 29060_003 | 4 | Target 65 mg | 65 | 108 | 108 | 134 | 134 | 11.355 | 3.020 | 7.574 | 0.824 | 3.781 | 3.131 | -2.355 | 9.917 |
| Multinomial IPTW | 29060_003_IMIPRAMINE | 29060_003 | 5 | Target 65 mg | 65 | 53 | 53 | 58 | 58 | 12.536 | 2.339 | 8.501 | 1.149 | 4.035 | 2.606 | -1.073 | 9.143 |
| Multinomial IPTW | 29060_003_IMIPRAMINE | 29060_003 | 6 | Target 65 mg | 65 | 73 | 73 | 68 | 68 | 13.407 | 2.313 | 10.208 | 1.163 | 3.199 | 2.589 | -1.876 | 8.273 |
| Multinomial IPTW | 29060_003_IMIPRAMINE | 29060_003 | 7 | Target 65 mg | 65 | 50 | 50 | 82 | 82 | 14.226 | 3.564 | 13.486 | 1.137 | 0.740 | 3.741 | -6.592 | 8.072 |
| Multinomial IPTW | 29060_003_IMIPRAMINE | 29060_003 | 8 | Target 65 mg | 65 | 2 | 2 | 1 | 1 | 15.045 | 5.274 | 12.693 | NA | 2.352 | NA | NA | NA |
| Multinomial IPTW | 29060_003_IMIPRAMINE | 29060_003 | 1 | Target 80 mg | 80 | 120 | 120 | 125 | 125 | 9.295 | 4.052 | 3.113 | 0.645 | 6.181 | 4.103 | -1.861 | 14.223 |
| Multinomial IPTW | 29060_003_IMIPRAMINE | 29060_003 | 2 | Target 80 mg | 80 | 144 | 144 | 144 | 144 | 11.972 | 1.767 | 4.042 | 0.681 | 7.930 | 1.894 | 4.218 | 11.642 |
| Multinomial IPTW | 29060_003_IMIPRAMINE | 29060_003 | 3 | Target 80 mg | 80 | 132 | 132 | 144 | 144 | 13.973 | 3.194 | 5.359 | 0.671 | 8.614 | 3.264 | 2.217 | 15.011 |
| Multinomial IPTW | 29060_003_IMIPRAMINE | 29060_003 | 4 | Target 80 mg | 80 | 108 | 108 | 134 | 134 | 14.715 | 3.457 | 7.574 | 0.824 | 7.141 | 3.554 | 0.175 | 14.107 |
| Multinomial IPTW | 29060_003_IMIPRAMINE | 29060_003 | 5 | Target 80 mg | 80 | 53 | 53 | 58 | 58 | 15.490 | 3.279 | 8.501 | 1.149 | 6.988 | 3.474 | 0.179 | 13.798 |
| Multinomial IPTW | 29060_003_IMIPRAMINE | 29060_003 | 6 | Target 80 mg | 80 | 73 | 73 | 68 | 68 | 16.281 | 3.540 | 10.208 | 1.163 | 6.073 | 3.726 | -1.231 | 13.376 |
| Multinomial IPTW | 29060_003_IMIPRAMINE | 29060_003 | 7 | Target 80 mg | 80 | 50 | 50 | 82 | 82 | 17.075 | 4.555 | 13.486 | 1.137 | 3.589 | 4.695 | -5.614 | 12.791 |
| Multinomial IPTW | 29060_003_IMIPRAMINE | 29060_003 | 8 | Target 80 mg | 80 | 2 | 2 | 1 | 1 | 17.869 | 5.974 | 12.693 | NA | 5.176 | NA | NA | NA |
| Multinomial IPTW | 29060_003_IMIPRAMINE | 29060_003 | 1 | Target 145 mg | 145 | 120 | 120 | 125 | 125 | 3.021 | 0.850 | 3.113 | 0.645 | -0.093 | 1.067 | -2.183 | 1.998 |
| Multinomial IPTW | 29060_003_IMIPRAMINE | 29060_003 | 2 | Target 145 mg | 145 | 144 | 144 | 144 | 144 | 6.781 | 0.876 | 4.042 | 0.681 | 2.739 | 1.109 | 0.565 | 4.913 |
| Multinomial IPTW | 29060_003_IMIPRAMINE | 29060_003 | 3 | Target 145 mg | 145 | 132 | 132 | 144 | 144 | 11.450 | 1.218 | 5.359 | 0.671 | 6.091 | 1.391 | 3.365 | 8.818 |
| Multinomial IPTW | 29060_003_IMIPRAMINE | 29060_003 | 4 | Target 145 mg | 145 | 108 | 108 | 134 | 134 | 13.946 | 1.201 | 7.574 | 0.824 | 6.372 | 1.457 | 3.517 | 9.227 |
| Multinomial IPTW | 29060_003_IMIPRAMINE | 29060_003 | 5 | Target 145 mg | 145 | 53 | 53 | 58 | 58 | 15.603 | 1.296 | 8.501 | 1.149 | 7.101 | 1.732 | 3.707 | 10.496 |
| Multinomial IPTW | 29060_003_IMIPRAMINE | 29060_003 | 6 | Target 145 mg | 145 | 73 | 73 | 68 | 68 | 16.840 | 1.605 | 10.208 | 1.163 | 6.632 | 1.982 | 2.748 | 10.516 |
| Multinomial IPTW | 29060_003_IMIPRAMINE | 29060_003 | 7 | Target 145 mg | 145 | 50 | 50 | 82 | 82 | 18.007 | 2.071 | 13.486 | 1.137 | 4.521 | 2.363 | -0.111 | 9.153 |
| Multinomial IPTW | 29060_003_IMIPRAMINE | 29060_003 | 8 | Target 145 mg | 145 | 2 | 2 | 1 | 1 | 19.174 | 2.612 | 12.693 | NA | 6.481 | NA | NA | NA |
| Multinomial IPTW | 29060_003_IMIPRAMINE | 29060_003 | 1 | Target 210 mg | 210 | 120 | 120 | 125 | 125 | 4.094 | 1.002 | 3.113 | 0.645 | 0.981 | 1.191 | -1.354 | 3.316 |
| Multinomial IPTW | 29060_003_IMIPRAMINE | 29060_003 | 2 | Target 210 mg | 210 | 144 | 144 | 144 | 144 | 5.009 | 1.079 | 4.042 | 0.681 | 0.967 | 1.276 | -1.533 | 3.467 |
| Multinomial IPTW | 29060_003_IMIPRAMINE | 29060_003 | 3 | Target 210 mg | 210 | 132 | 132 | 144 | 144 | 9.615 | 1.676 | 5.359 | 0.671 | 4.256 | 1.806 | 0.717 | 7.795 |
| Multinomial IPTW | 29060_003_IMIPRAMINE | 29060_003 | 4 | Target 210 mg | 210 | 108 | 108 | 134 | 134 | 11.326 | 1.660 | 7.574 | 0.824 | 3.752 | 1.853 | 0.121 | 7.384 |
| Multinomial IPTW | 29060_003_IMIPRAMINE | 29060_003 | 5 | Target 210 mg | 210 | 53 | 53 | 58 | 58 | 12.349 | 1.737 | 8.501 | 1.149 | 3.848 | 2.083 | -0.235 | 7.930 |
| Multinomial IPTW | 29060_003_IMIPRAMINE | 29060_003 | 6 | Target 210 mg | 210 | 73 | 73 | 68 | 68 | 13.027 | 2.037 | 10.208 | 1.163 | 2.819 | 2.346 | -1.778 | 7.416 |
| Multinomial IPTW | 29060_003_IMIPRAMINE | 29060_003 | 7 | Target 210 mg | 210 | 50 | 50 | 82 | 82 | 13.648 | 2.547 | 13.486 | 1.137 | 0.162 | 2.789 | -5.304 | 5.628 |
| Multinomial IPTW | 29060_003_IMIPRAMINE | 29060_003 | 8 | Target 210 mg | 210 | 2 | 2 | 1 | 1 | 14.269 | 3.171 | 12.693 | NA | 1.576 | NA | NA | NA |
| Multinomial IPTW | 29060_003_IMIPRAMINE | 29060_003 | 1 | Target 275 mg | 275 | 120 | 120 | 125 | 125 | 5.681 | 1.971 | 3.113 | 0.645 | 2.567 | 2.074 | -1.497 | 6.631 |
| Multinomial IPTW | 29060_003_IMIPRAMINE | 29060_003 | 2 | Target 275 mg | 275 | 144 | 144 | 144 | 144 | 4.894 | 1.774 | 4.042 | 0.681 | 0.851 | 1.900 | -2.872 | 4.575 |
| Multinomial IPTW | 29060_003_IMIPRAMINE | 29060_003 | 3 | Target 275 mg | 275 | 132 | 132 | 144 | 144 | 8.925 | 3.175 | 5.359 | 0.671 | 3.566 | 3.246 | -2.795 | 9.927 |
| Multinomial IPTW | 29060_003_IMIPRAMINE | 29060_003 | 4 | Target 275 mg | 275 | 108 | 108 | 134 | 134 | 10.498 | 3.083 | 7.574 | 0.824 | 2.923 | 3.191 | -3.331 | 9.177 |
| Multinomial IPTW | 29060_003_IMIPRAMINE | 29060_003 | 5 | Target 275 mg | 275 | 53 | 53 | 58 | 58 | 12.086 | 2.980 | 8.501 | 1.149 | 3.584 | 3.194 | -2.676 | 9.845 |
| Multinomial IPTW | 29060_003_IMIPRAMINE | 29060_003 | 6 | Target 275 mg | 275 | 73 | 73 | 68 | 68 | 13.681 | 3.137 | 10.208 | 1.163 | 3.473 | 3.346 | -3.084 | 10.030 |
| Multinomial IPTW | 29060_003_IMIPRAMINE | 29060_003 | 7 | Target 275 mg | 275 | 50 | 50 | 82 | 82 | 15.279 | 3.625 | 13.486 | 1.137 | 1.793 | 3.800 | -5.654 | 9.240 |
| Multinomial IPTW | 29060_003_IMIPRAMINE | 29060_003 | 8 | Target 275 mg | 275 | 2 | 2 | 1 | 1 | 16.876 | 4.338 | 12.693 | NA | 4.183 | NA | NA | NA |
| Multinomial IPTW | 29060_003_PAROXETINE | 29060_003 | 1 | Target 10 mg | 10 | 126 | 126 | 125 | 125 | 3.947 | 3.745 | 3.113 | 0.645 | 0.834 | 3.800 | -6.614 | 8.282 |
| Multinomial IPTW | 29060_003_PAROXETINE | 29060_003 | 2 | Target 10 mg | 10 | 153 | 153 | 144 | 144 | 5.925 | 2.765 | 4.042 | 0.681 | 1.883 | 2.848 | -3.699 | 7.464 |
| Multinomial IPTW | 29060_003_PAROXETINE | 29060_003 | 3 | Target 10 mg | 10 | 145 | 145 | 144 | 144 | 8.423 | 1.149 | 5.359 | 0.671 | 3.064 | 1.331 | 0.455 | 5.672 |
| Multinomial IPTW | 29060_003_PAROXETINE | 29060_003 | 4 | Target 10 mg | 10 | 115 | 115 | 134 | 134 | 9.562 | 0.927 | 7.574 | 0.824 | 1.988 | 1.240 | -0.443 | 4.419 |
| Multinomial IPTW | 29060_003_PAROXETINE | 29060_003 | 5 | Target 10 mg | 10 | 58 | 58 | 58 | 58 | 9.788 | 1.906 | 8.501 | 1.149 | 1.286 | 2.225 | -3.075 | 5.648 |
| Multinomial IPTW | 29060_003_PAROXETINE | 29060_003 | 6 | Target 10 mg | 10 | 79 | 79 | 68 | 68 | 9.405 | 4.274 | 10.208 | 1.163 | -0.803 | 4.430 | -9.485 | 7.879 |
| Multinomial IPTW | 29060_003_PAROXETINE | 29060_003 | 7 | Target 10 mg | 10 | 74 | 74 | 82 | 82 | 8.718 | 7.229 | 13.486 | 1.137 | -4.768 | 7.318 | -19.112 | 9.576 |
| Multinomial IPTW | 29060_003_PAROXETINE | 29060_003 | 8 | Target 10 mg | 10 | 6 | 6 | 1 | 1 | 7.980 | 10.312 | 12.693 | NA | -4.713 | NA | NA | NA |
| Multinomial IPTW | 29060_003_PAROXETINE | 29060_003 | 1 | Target 20 mg | 20 | 126 | 126 | 125 | 125 | 3.742 | 0.950 | 3.113 | 0.645 | 0.628 | 1.148 | -1.622 | 2.878 |
| Multinomial IPTW | 29060_003_PAROXETINE | 29060_003 | 2 | Target 20 mg | 20 | 153 | 153 | 144 | 144 | 6.806 | 1.109 | 4.042 | 0.681 | 2.764 | 1.301 | 0.214 | 5.314 |
| Multinomial IPTW | 29060_003_PAROXETINE | 29060_003 | 3 | Target 20 mg | 20 | 145 | 145 | 144 | 144 | 9.417 | 1.343 | 5.359 | 0.671 | 4.058 | 1.501 | 1.115 | 7.000 |
| Multinomial IPTW | 29060_003_PAROXETINE | 29060_003 | 4 | Target 20 mg | 20 | 115 | 115 | 134 | 134 | 11.873 | 1.255 | 7.574 | 0.824 | 4.299 | 1.501 | 1.356 | 7.241 |
| Multinomial IPTW | 29060_003_PAROXETINE | 29060_003 | 5 | Target 20 mg | 20 | 58 | 58 | 58 | 58 | 14.313 | 1.075 | 8.501 | 1.149 | 5.811 | 1.574 | 2.727 | 8.896 |
| Multinomial IPTW | 29060_003_PAROXETINE | 29060_003 | 6 | Target 20 mg | 20 | 79 | 79 | 68 | 68 | 16.742 | 1.021 | 10.208 | 1.163 | 6.534 | 1.547 | 3.501 | 9.567 |
| Multinomial IPTW | 29060_003_PAROXETINE | 29060_003 | 7 | Target 20 mg | 20 | 74 | 74 | 82 | 82 | 19.166 | 1.341 | 13.486 | 1.137 | 5.680 | 1.759 | 2.233 | 9.127 |
| Multinomial IPTW | 29060_003_PAROXETINE | 29060_003 | 8 | Target 20 mg | 20 | 6 | 6 | 1 | 1 | 21.590 | 1.901 | 12.693 | NA | 8.897 | NA | NA | NA |
| Multinomial IPTW | 29060_003_PAROXETINE | 29060_003 | 1 | Target 30 mg | 30 | 126 | 126 | 125 | 125 | 3.671 | 1.062 | 3.113 | 0.645 | 0.558 | 1.242 | -1.877 | 2.993 |
| Multinomial IPTW | 29060_003_PAROXETINE | 29060_003 | 2 | Target 30 mg | 30 | 153 | 153 | 144 | 144 | 6.396 | 0.829 | 4.042 | 0.681 | 2.354 | 1.072 | 0.252 | 4.456 |
| Multinomial IPTW | 29060_003_PAROXETINE | 29060_003 | 3 | Target 30 mg | 30 | 145 | 145 | 144 | 144 | 8.939 | 1.169 | 5.359 | 0.671 | 3.580 | 1.348 | 0.939 | 6.221 |
| Multinomial IPTW | 29060_003_PAROXETINE | 29060_003 | 4 | Target 30 mg | 30 | 115 | 115 | 134 | 134 | 10.594 | 1.121 | 7.574 | 0.824 | 3.020 | 1.391 | 0.293 | 5.747 |
| Multinomial IPTW | 29060_003_PAROXETINE | 29060_003 | 5 | Target 30 mg | 30 | 58 | 58 | 58 | 58 | 11.845 | 1.006 | 8.501 | 1.149 | 3.344 | 1.527 | 0.351 | 6.337 |
| Multinomial IPTW | 29060_003_PAROXETINE | 29060_003 | 6 | Target 30 mg | 30 | 79 | 79 | 68 | 68 | 12.827 | 0.977 | 10.208 | 1.163 | 2.618 | 1.519 | -0.358 | 5.595 |
| Multinomial IPTW | 29060_003_PAROXETINE | 29060_003 | 7 | Target 30 mg | 30 | 74 | 74 | 82 | 82 | 13.673 | 1.246 | 13.486 | 1.137 | 0.187 | 1.687 | -3.120 | 3.494 |
| Multinomial IPTW | 29060_003_PAROXETINE | 29060_003 | 8 | Target 30 mg | 30 | 6 | 6 | 1 | 1 | 14.497 | 1.724 | 12.693 | NA | 1.804 | NA | NA | NA |
| Multinomial IPTW | 29060_003_PAROXETINE | 29060_003 | 1 | Target 40 mg | 40 | 126 | 126 | 125 | 125 | 3.867 | 2.270 | 3.113 | 0.645 | 0.754 | 2.360 | -3.870 | 5.379 |
| Multinomial IPTW | 29060_003_PAROXETINE | 29060_003 | 2 | Target 40 mg | 40 | 153 | 153 | 144 | 144 | 8.085 | 1.459 | 4.042 | 0.681 | 4.042 | 1.610 | 0.887 | 7.198 |
| Multinomial IPTW | 29060_003_PAROXETINE | 29060_003 | 3 | Target 40 mg | 40 | 145 | 145 | 144 | 144 | 11.996 | 1.887 | 5.359 | 0.671 | 6.637 | 2.003 | 2.711 | 10.563 |
| Multinomial IPTW | 29060_003_PAROXETINE | 29060_003 | 4 | Target 40 mg | 40 | 115 | 115 | 134 | 134 | 13.578 | 1.910 | 7.574 | 0.824 | 6.004 | 2.080 | 1.926 | 10.082 |
| Multinomial IPTW | 29060_003_PAROXETINE | 29060_003 | 5 | Target 40 mg | 40 | 58 | 58 | 58 | 58 | 14.568 | 1.647 | 8.501 | 1.149 | 6.067 | 2.008 | 2.131 | 10.003 |
| Multinomial IPTW | 29060_003_PAROXETINE | 29060_003 | 6 | Target 40 mg | 40 | 79 | 79 | 68 | 68 | 15.165 | 1.358 | 10.208 | 1.163 | 4.957 | 1.788 | 1.452 | 8.461 |
| Multinomial IPTW | 29060_003_PAROXETINE | 29060_003 | 7 | Target 40 mg | 40 | 74 | 74 | 82 | 82 | 15.564 | 1.578 | 13.486 | 1.137 | 2.078 | 1.945 | -1.734 | 5.889 |
| Multinomial IPTW | 29060_003_PAROXETINE | 29060_003 | 8 | Target 40 mg | 40 | 6 | 6 | 1 | 1 | 15.930 | 2.284 | 12.693 | NA | 3.237 | NA | NA | NA |
| Multinomial IPTW | 29060_003_PAROXETINE | 29060_003 | 1 | Target 50 mg | 50 | 126 | 126 | 125 | 125 | 1.580 | 2.219 | 3.113 | 0.645 | -1.533 | 2.311 | -6.062 | 2.997 |
| Multinomial IPTW | 29060_003_PAROXETINE | 29060_003 | 2 | Target 50 mg | 50 | 153 | 153 | 144 | 144 | 7.867 | 1.096 | 4.042 | 0.681 | 3.825 | 1.290 | 1.296 | 6.353 |
| Multinomial IPTW | 29060_003_PAROXETINE | 29060_003 | 3 | Target 50 mg | 50 | 145 | 145 | 144 | 144 | 11.204 | 2.117 | 5.359 | 0.671 | 5.845 | 2.221 | 1.492 | 10.199 |
| Multinomial IPTW | 29060_003_PAROXETINE | 29060_003 | 4 | Target 50 mg | 50 | 115 | 115 | 134 | 134 | 13.357 | 1.956 | 7.574 | 0.824 | 5.783 | 2.122 | 1.623 | 9.942 |
| Multinomial IPTW | 29060_003_PAROXETINE | 29060_003 | 5 | Target 50 mg | 50 | 58 | 58 | 58 | 58 | 14.146 | 1.620 | 8.501 | 1.149 | 5.644 | 1.986 | 1.752 | 9.537 |
| Multinomial IPTW | 29060_003_PAROXETINE | 29060_003 | 6 | Target 50 mg | 50 | 79 | 79 | 68 | 68 | 14.025 | 1.498 | 10.208 | 1.163 | 3.817 | 1.897 | 0.100 | 7.534 |
| Multinomial IPTW | 29060_003_PAROXETINE | 29060_003 | 7 | Target 50 mg | 50 | 74 | 74 | 82 | 82 | 13.451 | 2.062 | 13.486 | 1.137 | -0.036 | 2.355 | -4.652 | 4.580 |
| Multinomial IPTW | 29060_003_PAROXETINE | 29060_003 | 8 | Target 50 mg | 50 | 6 | 6 | 1 | 1 | 12.800 | 3.023 | 12.693 | NA | 0.107 | NA | NA | NA |

<small><em>Comparison of dose-versus-placebo predictions using ordinal
versus multinomial treatment weights</em></small>
