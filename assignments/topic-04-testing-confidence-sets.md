# ORF 524 Practice Module 4: Testing and Confidence Sets

**Status:** Ready for Fall 2026  
**Last updated:** August 30, 2026  
**Primary chapter:** [Week 4: Hypothesis testing and confidence sets](../lectures/week-04.md)

## Purpose

This ungraded module develops tests as calibrated decision rules and confidence sets as inverted test families. Problems 1--4 form the core bank. Problems 5--7 are optional extensions on combining p values, prediction, and exact two sample inference.

<a id="problem-map"></a>

## Problem map

| Problem | Bank | Main task |
|---|---|---|
| [1. Exact Normal inference with unknown variance](#problem-1) | Core | Derive Student tests, power, local power, and confidence sets by inversion. |
| [2. Neyman--Pearson, discreteness, and MLR](#problem-2) | Core | Use likelihood ratios, randomization, and monotone likelihood ratio to study optimal tests. |
| [3. P-values as random variables](#problem-3) | Core | Establish exact and valid p value behavior under simple, composite, and discrete nulls. |
| [4. Test inversion with parameter dependent support](#problem-4) | Core | Invert exact tests in a Uniform model with moving support. |
| [5. Fisher's combination under independence](#problem-5) | Further | Derive Fisher's combination rule and audit its independence assumption. |
| [6. A prediction interval is not a confidence interval](#problem-6) | Further | Construct a prediction interval and distinguish its target from a confidence interval. |
| [7. Two sample Student inference](#problem-7) | Further | Derive the pooled Student statistic, likelihood ratio test, power, and confidence interval. |

## How to work with this module

For each test, identify the null and alternative parameter sets, test statistic, rejection region, and power function before making an optimality claim. Separate calibration from comparison: first verify level over the entire null, then ask whether the test is most powerful, UMP, or UMP only within a restricted class. For confidence sets, write the test family being inverted and verify coverage as an event identity.

AI may be used on these ungraded problems under the syllabus policy and the staged protocol in [`AGENTS.md`](AGENTS.md). Begin with a genuine attempt, request the least revealing useful hint, and verify every generated probability statement under the correct sampling law. The closed book midterms remain fully unaided.

For released exam practice, use the [question level first half guide](../exams/README.md#question-level-guide-for-the-first-half). Confirm that the listed material has been covered, attempt the selected question before opening its solution, and treat historical exam instructions as superseded by the current syllabus.

## Core practice

<a id="problem-1"></a>

### Problem 1. Exact Normal inference with unknown variance

**Chapter connections:** [Sections 2.3--2.4: Two sided and Student inference](../lectures/week-04.md#23-the-other-one-sided-direction-and-the-two-sided-test), [Sections 4.4--4.5: Two sided restricted optimality](../lectures/week-04.md#44-why-two-sided-ump-tests-usually-fail), and [Sections 5.2--5.3: Test inversion and the Normal interval](../lectures/week-04.md#52-inverting-tests)

Let $X_1,\ldots,X_n$ be i.i.d. $\mathcal N(\mu,\sigma^2)$, where $n\geq2$ and both parameters are unknown, and let $S_n^2=(n-1)^{-1}\sum_i(X_i-\overline X_n)^2$.

1. Show that

$$
T_n(\mu)=\frac{\sqrt n(\overline X_n-\mu)}{S_n}
\sim t_{n-1}.
$$

2. Construct a level $\alpha$ test of

$$
H_0:\mu\leq\mu_0
\quad\text{versus}\quad
H_1:\mu>\mu_0.
$$

   Derive its power in terms of a noncentral $t$ c.d.f. and justify the supremum over the composite null. Show that power tends to one for every fixed $\mu>\mu_0$. Then, for fixed $\sigma$ and $h>0$, find the limiting power under the local alternatives $\mu_n=\mu_0+h/\sqrt n$. Treat this last calculation as a bridge to the Week 5 limit theorems.

3. Construct the equal tail level $\alpha$ test of $H_0:\mu=\mu_0$ against $H_1:\mu\neq\mu_0$. Derive its power in terms of a noncentral $t$ c.d.f. and describe its shape as a function of $\mu$ for fixed $\sigma$. Invert the family of tests to obtain a confidence interval.

4. Is the two sided test unrestricted UMP? State the relevant restricted optimality conclusion and explain the distinction.

<!-- Source lineage: Fall 2024 pset4, “Normal testing, unknown variance”; shortened to one one-sided test, one two-sided test, power, and inversion. -->

[Back to the problem map](#problem-map)

<a id="problem-2"></a>

### Problem 2. Neyman--Pearson, discreteness, and MLR

**Chapter connections:** [Section 1: Tests as decision rules](../lectures/week-04.md#1-tests-as-decision-rules), [Section 3.2: Neyman--Pearson lemma](../lectures/week-04.md#32-neyman--pearson-lemma), and [Section 4.2: Monotone likelihood ratio](../lectures/week-04.md#42-monotone-likelihood-ratio)

Let $X_1,\ldots,X_n$ be i.i.d. $\mathsf{Bernoulli}(p)$ and $T=\sum_iX_i$.

1. For $0<p_0<p_1<1$, show that the likelihood ratio in favor of $p_1$ is strictly increasing in $T$.

2. Choose an integer $c$ satisfying

$$
\mathbb{P}_{p_0}(T>c)\leq\alpha\leq\mathbb{P}_{p_0}(T\geq c).
$$

   Find the randomization probability $\gamma$ at $T=c$ that makes the test have exact size $\alpha$.

3. Use the Neyman--Pearson lemma to show that this test is most powerful for $p=p_0$ versus $p=p_1$. Then use monotone likelihood ratio to state the UMP conclusion for

$$
H_0:p\leq p_0
\quad\text{versus}\quad
H_1:p>p_0.
$$

4. Explain why “most powerful against $p_1$” and “UMP against every $p>p_0$” are distinct claims, even though the same rejection rule works here.

<!-- Source lineage: Fall 2024 pset4, “MLR property” and “UMP inference”; corrected to distinguish simple-alternative Neyman--Pearson optimality from composite-alternative UMP optimality and to include randomization. -->

[Back to the problem map](#problem-map)

<a id="problem-3"></a>

### Problem 3. P-values as random variables

**Chapter connection:** [Section 2.5: P-values](../lectures/week-04.md#25-p-values)

1. Let $T$ have continuous c.d.f. $F_0$ under a simple null and suppose large $T$ is evidence against the null. Show that $p(T)=1-F_0(T)$ is $\mathsf{Uniform}(0,1)$ under the null.

2. In the known variance Normal model, test $H_0:\mu\leq\mu_0$ against $H_1:\mu>\mu_0$ with

$$
p(X)=1-\Phi\mkern-3mu\left(\frac{\sqrt n(\overline X_n-\mu_0)}{\sigma}\right).
$$

   Show that this p-value is valid over the composite null and exactly Uniform at $\mu=\mu_0$.

3. Let $T\sim\mathsf{Binomial}(n,p_0)$ under a simple null and define

$$
p(T)=\mathbb{P}_{p_0}(T'\geq T),
$$

   where $T'$ is an independent null copy. Show that $p(T)$ is valid but generally not exactly Uniform.

4. For each setting, explain why the p-value is not the posterior probability that the null is true.

<!-- Source lineage: Fall 2024 pset4, “p-value”; extended from a simple continuous null to composite and discrete validity. -->

[Back to the problem map](#problem-map)

<a id="problem-4"></a>

### Problem 4. Test inversion with parameter dependent support

**Chapter connections:** [Section 5.1: Coverage](../lectures/week-04.md#51-coverage) and [Section 5.2: Inverting tests](../lectures/week-04.md#52-inverting-tests)

Let $X_1,\ldots,X_n$ be i.i.d. $\mathsf{Uniform}(0,\theta)$, where $\theta>0$, and let $M_n=\max_iX_i$. For $0<\alpha<1$, define

$$
q_L=(\alpha/2)^{1/n},
\qquad
q_U=(1-\alpha/2)^{1/n}.
$$

1. Show that $M_n/\theta$ has c.d.f. $u^n$ on $[0,1]$.

2. For every candidate $\theta_0$, construct an exact level $\alpha$ test of $H_0:\theta=\theta_0$ that rejects when $M_n/\theta_0<q_L$ or $M_n/\theta_0>q_U$.

3. Invert these tests and show that the resulting confidence set is

$$
C(X)
=\left[\frac{M_n}{q_U},\frac{M_n}{q_L}\right]\cap(0,\infty).
$$

   Verify its coverage directly.

4. On the event $M_n>0$, the interval does not contain the likelihood maximizer $M_n$. Explain why this is not a contradiction. What is the inverted set on the probability zero sample $M_n=0$, and why does that edge case not affect coverage?

<!-- Source lineage: Uniform likelihood exercises and legacy confidence-set problems; new inversion problem emphasizing moving support and exact coverage. -->

[Back to the problem map](#problem-map)

## Further practice

<a id="problem-5"></a>

### Problem 5. Fisher's combination under independence

**Chapter connection:** [Section 2.5: P-values](../lectures/week-04.md#25-p-values)

Suppose $p_1,\ldots,p_K$ are independent p-values that are exactly $\mathsf{Uniform}(0,1)$ under a common simple null. Show that

$$
-2\sum_{k=1}^K\log p_k\sim\chi^2_{2K}.
$$

Identify which conclusion can fail if the p-values are dependent, discrete, or merely valid rather than exactly Uniform.

<!-- Source lineage: Fall 2024 pset4, “p-value,” Fisher-combination extension. -->

[Back to the problem map](#problem-map)

<a id="problem-6"></a>

### Problem 6. A prediction interval is not a confidence interval

**Chapter connection:** [Section 5.1: Coverage](../lectures/week-04.md#51-coverage)

Let $X_1,\ldots,X_n,X_{n+1}$ be i.i.d. Exponential with mean $\theta>0$.

1. Show that $2n\overline X_n/\theta\sim\chi^2_{2n}$ and

$$
\frac{X_{n+1}}{\overline X_n}\sim F_{2,2n}.
$$

2. Construct an equal tail prediction interval with coverage $1-\alpha$ for the future observation $X_{n+1}$.

3. Explain precisely how its coverage statement differs from a confidence interval for $\theta$.

<!-- Source lineage: Fall 2024 pset4, “Confidence interval for prediction.” -->

[Back to the problem map](#problem-map)

<a id="problem-7"></a>

### Problem 7. Two sample Student inference

**Chapter connections:** [Section 2.4: Exact Student inference with unknown variance](../lectures/week-04.md#24-exact-student-inference-with-unknown-variance), [Section 4.5: Unbiased tests and two sided optimality](../lectures/week-04.md#45-unbiased-tests-and-two-sided-optimality), and [Section 5.2: Inverting tests](../lectures/week-04.md#52-inverting-tests)

Let $X_1,\ldots,X_m$ and $Y_1,\ldots,Y_n$ be independent samples from $\mathcal N(\mu_X,\sigma^2)$ and $\mathcal N(\mu_Y,\sigma^2)$, respectively, where $m,n\geq2$ and the common variance is unknown. Write $N=m+n$, $\Delta=\mu_X-\mu_Y$, $D=\overline X_m-\overline Y_n$, and

$$
S_p^2
=\frac{(m-1)S_X^2+(n-1)S_Y^2}{N-2}.
$$

1. Derive the exact distribution of $D$.

2. Prove that

$$
\frac{(N-2)S_p^2}{\sigma^2}\sim\chi^2_{N-2}
$$

   and that this quantity is independent of $D$.

3. For a candidate difference $\Delta_0$, show that

$$
T(\Delta_0)
=\frac{D-\Delta_0}{S_p\sqrt{1/m+1/n}}
$$

   has a central $t_{N-2}$ distribution under $H_0:\Delta=\Delta_0$ and a noncentral $t_{N-2}$ distribution away from the null. Identify its noncentrality parameter.

4. For $H_0:\Delta=0$, derive the restricted and unrestricted maximum likelihood estimators and show that the likelihood ratio is

$$
\Lambda(X,Y)
=\left\lbrace1+\frac{T(0)^2}{N-2}\right\rbrace^{-N/2}.
$$

5. Construct the equal tail level $\alpha$ likelihood ratio test. Write its exact power in terms of a noncentral $t$ c.d.f. and prove consistency against every fixed $\Delta\neq0$ along sequences for which $m,n\to\infty$ and $m/N\to\rho\in(0,1)$. State the standard restricted optimality conclusion and explain why it is not an unrestricted UMP conclusion.

6. Invert the test family for $H_0:\Delta=\Delta_0$ to construct an exact confidence interval with coverage $1-\alpha$ for $\Delta$.

<!-- Source lineage: Fall 2019 through Fall 2023 assignment banks, “Two sided Normal testing, differences in means”; restored with the exact sampling argument, likelihood ratio calculation, power, consistency, restricted optimality, and inversion. -->

[Back to the problem map](#problem-map)

## Completion check

After the core, you should be able to explain:

1. how a test statistic, critical value, rejection region, and power function describe a testing procedure;

2. the difference among level, size, power, and the two error probabilities;

3. when randomization is needed to attain exact size;

4. why fixed and local alternative power describe different operating characteristics;

5. why a valid p-value need not be exactly Uniform;

6. how Neyman--Pearson optimality differs from a UMP conclusion and why MLR can bridge them;

7. why two sided UMP tests generally fail and what UMPU means; and

8. why test inversion automatically produces a coverage statement.
