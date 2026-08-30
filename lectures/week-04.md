# ORF 524: Statistical Theory and Methods

## Week 4: Hypothesis testing and confidence sets

**Status:** Ready for Fall 2026  
**Last updated:** August 30, 2026  
**Meetings:** Monday, September 21, and Wednesday, September 23, 2026

<!-- Legacy sources: 11 ORF 524/lectures/ORF524-Fall2025-Lecture-Week-04.tex and
ORF524-Fall2025-Lecture-Week-04_Lecture.pdf; 09 ORF 524/lectures/
ORF524-Fall2024-Lecture-Week-04_Notes.pdf. -->

## Central question

How can a decision rule control false rejections, retain power against relevant alternatives, and produce confidence sets with calibrated repeated sampling coverage?

## Learning goals

By the end of the week, you should be able to:

1. formulate simple and composite hypotheses and distinguish level, size, power, and both error probabilities for possibly randomized tests;
2. construct exact Normal and Student tests from pivotal quantities and interpret valid p-values;
3. derive the Neyman--Pearson most powerful test for simple hypotheses;
4. use monotone likelihood ratio for one sided UMP tests and unbiasedness for two sided UMPU tests;
5. compare power under fixed and local alternatives; and
6. construct confidence sets by inverting tests and interpret their coverage correctly.

## In-class route

This is one continuous route through the week, not a division into class meetings. Proceed in order, stop wherever the discussion ends, and resume at the next unfinished stop. The linked sections contain the full student facing exposition; the route identifies the concepts and activities to foreground live.

| Stop | Live focus | Mode |
|---|---|---|
| **W4.1** | [Tests as decision rules](#w4-stop-1) | Discuss + Checkpoint 1 |
| **W4.2** | [Size, power, and exact Normal calibration](#w4-stop-2) | Discuss + board work |
| **W4.3** | [Student inference and p-values](#w4-stop-3) | Discuss + AI interpretation audit + Checkpoint 2 |
| **W4.4** | [Neyman--Pearson and most powerful tests](#w4-stop-4) | Discuss proof + board example + AI proof audit |
| **W4.5** | [MLR, UMP tests, and two sided restricted optimality](#w4-stop-5) | Discuss + Checkpoint 3 |
| **W4.6** | [Confidence sets by test inversion](#w4-stop-6) | Board work + AI interpretation audit + Checkpoint 4 |

Normal sampling distribution and power curve diagrams belong in Stops W4.2--W4.3. The geometric likelihood ratio argument is the live core of Stop W4.4; extended two sided calculations can be developed selectively within Stop W4.5.

## How to use this chapter

**Prepare:** Read Sections 1 and 2 and review Normal and Student $t$ sampling distributions. Attempt Checkpoint 1 before class.

**In class:** Follow the continuous route above through calibration, power, p-values, Neyman--Pearson, monotone likelihood ratio, two sided restricted optimality, and confidence sets by test inversion. Prepared responses supply an alternative for every live prompt if AI is unavailable.

**Review:** Re derive the Normal examples without looking at the formulas. Use the AI dialogues to audit interpretations and proof steps; the formal definitions in the chapter remain authoritative.

**Practice:** The matching assignment complements this chapter with ungraded problems for consolidating the material, extending selected ideas, and preparing for cumulative review.

**Prerequisites:** Weeks 1--3; likelihood ratios; decision rules and risk; sufficient statistics; Normal, chi squared, and Student $t$ distributions.

## Full chapter map

1. [Tests as decision rules](#1-tests-as-decision-rules)
2. [Calibration, power, and p-values](#2-calibration-power-and-p-values)
3. [Most powerful tests](#3-most-powerful-tests)
4. [Monotone likelihood ratio and restricted optimality](#4-monotone-likelihood-ratio-and-restricted-optimality)
5. [Confidence sets and test inversion](#5-confidence-sets-and-test-inversion)
6. [Recap](#6-recap)

<a id="w4-stop-1"></a>

## 1. Tests as decision rules

Let $X\sim \mathbb{P}_\theta$ with $\theta\in\Theta$, and partition the parameter space into disjoint sets $\Theta_0$ and $\Theta_1$. We consider

$$
H_0:\theta\in\Theta_0
\qquad\text{versus}\qquad
H_1:\theta\in\Theta_1.
$$

Equivalently, testing chooses between two subsets of the maintained statistical model $\mathcal{P}=\lbrace\mathbb{P}_\theta:\theta\in\Theta\rbrace$. The parametric notation is convenient, but the same decision structure applies when the hypotheses are sets of probability laws in a nonparametric model.

A **nonrandomized test** is a statistic $\phi(X)\in\lbrace0,1\rbrace$, where $\phi=1$ means reject $H_0$. Equivalently, there is a rejection region $C$ such that

$$
\phi(X)=\mathbf 1\lbrace X\in C\rbrace.
$$

The set $C$ is the **rejection region** or **critical region**, and its complement is the **nonrejection region**. Many tests are represented by a real valued **test statistic** $T(X)$ and a **critical value** $c$, for example $C=\lbrace x:T(x)>c\rbrace$. Choosing the statistic determines how the sample is ordered from less to more incompatible with the null; choosing the critical value calibrates how much null rejection is allowed.

A **randomized test** takes values in $[0,1]$; conditional on $X=x$, the null is rejected with probability $\phi(x)$. Randomization is sometimes necessary to attain an exact level in a discrete model.

The **power function** is

$$
\beta_\phi(\theta)=\mathbb E_\theta[\phi(X)].
$$

Thus:

- for $\theta\in\Theta_0$, $\beta_\phi(\theta)$ is the probability of type I error;
- for $\theta\in\Theta_1$, $1-\beta_\phi(\theta)$ is the probability of type II error.

The two hypotheses are **simple** if each contains one distribution and **composite** otherwise.

### 1.1 Running Normal example: direction determines the rejection region

Let $X_1,\ldots,X_n$ be i.i.d. $\mathcal N(\mu,\sigma^2)$ with $\sigma^2$ known. Before calibrating a critical value or proving optimality, the direction of the alternative suggests three natural rules based on $\overline X_n$:

| Hypotheses | Evidence against $H_0$ | Candidate rejection region |
|---|---|---|
| $H_0:\mu\leq\mu_0$ versus $H_1:\mu>\mu_0$ | $\overline X_n$ is large | $\lbrace\overline X_n>c\rbrace$ |
| $H_0:\mu\geq\mu_0$ versus $H_1:\mu<\mu_0$ | $\overline X_n$ is small | $\lbrace\overline X_n<c\rbrace$ |
| $H_0:\mu=\mu_0$ versus $H_1:\mu\neq\mu_0$ | $\overline X_n$ is far from $\mu_0$ | $\lbrace\lvert\overline X_n-\mu_0\rvert>c\rbrace$ |

This first step is only structural. It does not yet explain why $\overline X_n$ should be used rather than another statistic, how $c$ should be chosen, or whether the resulting rule is optimal. Calibration answers the second question; Neyman--Pearson and monotone likelihood ratio later answer the first and third in this model. [Week 9](week-09.md#6-wald-score-and-likelihood-ratio-inference) shows how likelihood ratio, score, and Wald constructions generate test statistics beyond settings with an exact finite sample solution.

The two possible decisions and the two possible states produce the following error table. “Do not reject” is deliberately different from “accept”: a test may fail to find sufficient evidence against $H_0$ without establishing that $H_0$ is true.

| True state | Do not reject $H_0$ | Reject $H_0$ |
|---|---:|---:|
| $H_0$ is true | Correct decision | Type I error |
| $H_1$ is true | Type II error | Correct decision |

Reducing the rejection region tends to reduce type I error but increase type II error; enlarging it tends to do the reverse. The Neyman approach handles this tradeoff by first controlling type I error and then seeking high power against alternatives.

### Checkpoint 1

1. Is a test a statistic, a decision rule, or both?
2. Why is the type I error probability a function of $\theta$ under a composite null?
3. What role can randomization play in a discrete model?
4. Can both type I and type II error probabilities be made uniformly smaller without changing the experiment?

<a id="w4-stop-2"></a>

## 2. Calibration, power, and p-values

### 2.1 Level and size

The **size** of a test is

$$
\sup_{\theta\in\Theta_0}\beta_\phi(\theta).
$$

The test has **level $\alpha$** if

$$
\sup_{\theta\in\Theta_0}\mathbb E_\theta[\phi(X)]\leq\alpha.
$$

A test can have level $\alpha$ while having size strictly below $\alpha$. Calibrating rejection probability at one convenient null value does not control a composite null unless that value is least favorable or another argument controls the supremum.

Calibration requires the null distribution of the chosen statistic. A **nuisance parameter** is an unknown component of the model that is not the target of the hypothesis but may affect that distribution. In favorable cases the exact finite sample null law is known or can be made nuisance free by a pivot. Later in the course, asymptotic approximations provide calibration when no convenient exact law is available.

### 2.2 Exact Normal test with known variance

A **pivotal quantity** is a function of the data and parameter whose sampling distribution does not depend on the unknown parameter. A pivot can calibrate a test by comparing its null value with fixed quantiles.

A pivot need not be a statistic because it may contain the unknown parameter, as $\sqrt n(\overline X_n-\mu)/\sigma$ does. An ancillary statistic, by contrast, is a function of the data alone whose distribution is parameter free. A data only pivot is ancillary, but the terms are not interchangeable in general.

Let $X_1,\ldots,X_n$ be i.i.d. $\mathcal N(\mu,\sigma^2)$ with $\sigma^2$ known. For

$$
H_0:\mu\leq\mu_0
\qquad\text{versus}\qquad
H_1:\mu>\mu_0,
$$

define

$$
Z=\frac{\sqrt n(\overline X_n-\mu_0)}{\sigma}.
$$

The level $\alpha$ test rejects when

$$
Z>z_{1-\alpha},
$$

where $z_q$ is the $q$ quantile of $\mathcal N(0,1)$. At the boundary $\mu=\mu_0$, rejection probability is exactly $\alpha$; for $\mu<\mu_0$, it is smaller. The power function is

$$
\beta(\mu)
=1-\Phi\mkern-3mu\left(
z_{1-\alpha}-\frac{\sqrt n(\mu-\mu_0)}{\sigma}
\right).
$$

This formula displays the roles of effect size, noise, sample size, and level.

For every fixed $\mu>\mu_0$, power tends to one as $n\to\infty$. A more demanding comparison uses a **local alternative** approaching the null boundary at the standard error rate,

$$
\mu_n=\mu_0+\frac{h}{\sqrt n}.
$$

Under this sequence, $Z\sim\mathcal N(h/\sigma,1)$ and

$$
\beta_n(\mu_n)
=1-\Phi\mkern-3mu\left(z_{1-\alpha}-\frac h\sigma\right).
$$

The nondegenerate limit records which $n^{-1/2}$ departures the experiment can reliably detect; consistency against fixed alternatives alone does not reveal this resolution.

### 2.3 The other one sided direction and the two sided test

The same pivot produces the other two tests from the running example. Writing

$$
Z=\frac{\sqrt n(\overline X_n-\mu_0)}{\sigma},
$$

we obtain:

| Hypotheses | Level $\alpha$ rejection rule |
|---|---|
| $H_0:\mu\leq\mu_0$ versus $H_1:\mu>\mu_0$ | $Z>z_{1-\alpha}$ |
| $H_0:\mu\geq\mu_0$ versus $H_1:\mu<\mu_0$ | $Z<z_{\alpha}$ |
| $H_0:\mu=\mu_0$ versus $H_1:\mu\neq\mu_0$ | $\lvert Z\rvert>z_{1-\alpha/2}$ |

For the two sided rule, put $\delta=\sqrt n(\mu-\mu_0)/\sigma$ and $z=z_{1-\alpha/2}$. Since $Z\sim\mathcal N(\delta,1)$ under $\mu$, its power is

$$
\begin{aligned}
\beta(\mu)
&=\mathbb{P}_\mu(Z<-z)+\mathbb{P}_\mu(Z>z)\\
&=\Phi(-z-\delta)+1-\Phi(z-\delta).
\end{aligned}
$$

Thus $\beta(\mu_0)=\alpha$, while $\beta(\mu)\to1$ as $\mu\to-\infty$ or $\mu\to+\infty$. The one sided power curve is monotone in the direction of the alternative; the two sided curve has its minimum at the null and rises in both directions. These shapes are operating characteristics of the procedures, not pictures of posterior uncertainty about $\mu$.

> [!IMPORTANT]
> **Board work 1 — Normal power and sample size**
>
> Derive the power function and verify that its supremum over $\mu\leq\mu_0$ occurs at $\mu=\mu_0$. Then determine the sample size needed to attain power $1-\gamma$ at a specified alternative $\mu_1>\mu_0$.
>
> Draw two linked pictures while doing the calculation: first, the sampling distributions of $\overline X_n$ under the boundary null and several alternatives with the common rejection cutoff; second, the resulting power curve as a function of $\mu$. Use the pictures to explain why the null boundary determines size and why power increases under the one sided alternative. Then sketch the two sided power curve, mark $\beta(\mu_0)=\alpha$, and explain why it rises in both directions.

<a id="w4-stop-3"></a>

### 2.4 Exact Student inference with unknown variance

Now suppose both $\mu$ and $\sigma^2$ are unknown and define

$$
S_n^2=\frac1{n-1}\sum_{i=1}^n(X_i-\overline X_n)^2,
\qquad
T_n(\mu_0)=\frac{\sqrt n(\overline X_n-\mu_0)}{S_n}.
$$

Under $\mu=\mu_0$, the Normal sampling theorem gives

$$
T_n(\mu_0)\sim t_{n-1}.
$$

Consequently, rejecting when $T_n(\mu_0)>t_{n-1,1-\alpha}$ gives an exact level $\alpha$ test of $H_0:\mu\leq\mu_0$ against $H_1:\mu>\mu_0$. The boundary rejection probability is exactly $\alpha$ for every $\sigma>0$, and it decreases when $\mu<\mu_0$. For the point null $H_0:\mu=\mu_0$, the equal tail test rejects when

$$
|T_n(\mu_0)|>t_{n-1,1-\alpha/2}.
$$

Under the Normal family theorem this two sided test is UMPU, not unrestricted UMP. Inverting it gives the exact Student interval

$$
\left[
\overline X_n-t_{n-1,1-\alpha/2}\frac{S_n}{\sqrt n},
\overline X_n+t_{n-1,1-\alpha/2}\frac{S_n}{\sqrt n}
\right].
$$

### 2.5 P-values

A **valid p-value** is a statistic $p(X)\in[0,1]$ satisfying

$$
\sup_{\theta\in\Theta_0}
\mathbb{P}_\theta\lbrace p(X)\leq u\rbrace
\leq u
\qquad\text{for every }u\in[0,1].
$$

Rejecting whenever $p(X)\leq\alpha$ therefore defines a level $\alpha$ test.

A p-value can also be constructed from a nested family of tests. If $\phi_a$ is the level $a$ test and the rejection regions expand as $a$ increases, then

$$
p(x)=\inf\lbrace a\in[0,1]:\phi_a(x)=1\rbrace,
$$

with the convention that the infimum is one if the set is empty, is the smallest nominal level at which the observed sample would lead to rejection. In the usual upper tail construction based on a statistic $T$, this becomes

$$
p(x)=\sup_{\theta\in\Theta_0}\mathbb{P}_\theta\lbrace T(X)\geq T(x)\rbrace,
$$

with the tail convention adjusted to the test and to possible randomization. Under a continuous boundary null this p-value is often exactly $\mathsf{Uniform}(0,1)$; with a discrete statistic, the attainable tail probabilities are typically conservative rather than exactly Uniform.

For the one sided Normal test above,

$$
p(x)=1-\Phi\mkern-3mu\left(
\frac{\sqrt n(\overline x_n-\mu_0)}{\sigma}
\right).
$$

The p-value is not $\mathbb{P}(H_0\mid X)$, the probability that the null is true, or the probability that the result occurred by chance. It measures how small a level is needed for this test to reject at the observed data.

The nested test picture makes this operational: draw the null distribution of $T$, mark the observed value $T(x)$, and move the upper tail critical value until it reaches that observation. The corresponding tail area is the p-value. More extreme values move farther into the tail and therefore enter the rejection region at smaller levels.

> [!TIP]
> **AI interaction 1 — Interpret this p-value**
>
> Copy the prompt below into an AI interface and audit the response.

```text
An iid Normal sample with known variance produces a one sided p-value of 0.03
for testing H0: mu <= 0 against H1: mu > 0.

Give three mathematically correct interpretations and three common incorrect
interpretations. For every probability statement, identify what is treated as
random and under which distribution. Do not interpret the p-value as a
posterior probability or as the probability that the result occurred by chance.
```

**Audit question:** Does the response condition on the observed data while simultaneously treating the same data as random?

### Checkpoint 2

1. Does a smaller p-value imply a larger effect size?
2. If $p=0.03$, which tests in the nested family $\lbrace\phi_\alpha\rbrace$ reject?
3. Why can a discrete p-value be conservative rather than exactly $\mathsf{Uniform}(0,1)$ under the null?
4. What feature of a local alternative power calculation is invisible under a fixed alternative?

<a id="w4-stop-4"></a>

## 3. Most powerful tests

### 3.1 The optimization problem

For simple hypotheses

$$
H_0:\theta=\theta_0
\qquad\text{versus}\qquad
H_1:\theta=\theta_1,
$$

a **most powerful level $\alpha$ test** maximizes

$$
\mathbb E_{\theta_1}[\phi(X)]
$$

subject to

$$
\mathbb E_{\theta_0}[\phi(X)]\leq\alpha.
$$

This is constrained optimization: spend a limited type I error budget on sample points that are most diagnostic of the alternative.

### 3.2 Neyman--Pearson lemma

Suppose the laws under the two simple hypotheses have densities $f_0$ and $f_1$ with respect to a common measure. A test of the form

$$
\phi^\star(x)=
\begin{cases}
1, & f_1(x)>k f_0(x),\\
\gamma, & f_1(x)=k f_0(x),\\
0, & f_1(x)<k f_0(x),
\end{cases}
$$

with $k\geq0$ and $\gamma\in[0,1]$ chosen to attain level $\alpha$, is most powerful against $\theta_1$. Conversely, under the standard conditions, a most powerful test has this form up to null sets.

The statistic $f_1(X)/f_0(X)$ is the likelihood ratio in favor of the alternative. Randomization on the boundary may be needed to spend the type I error budget exactly.

### Proof

Let $\nu$ be the common dominating measure. For any competing level $\alpha$ test $\phi$, compare it pointwise with $\phi^\star$.

- On $\lbrace f_1>kf_0\rbrace$, we have $\phi^\star=1\geq\phi$, and both factors below are nonnegative.
- On $\lbrace f_1<kf_0\rbrace$, we have $\phi^\star=0\leq\phi$, and both factors below are nonpositive.
- On $\lbrace f_1=kf_0\rbrace$, the second factor is zero, regardless of the boundary randomization.

Therefore

$$
(\phi^\star-\phi)(f_1-kf_0)\geq0
$$

pointwise almost everywhere with respect to $\nu$. Integrating with respect to $\nu$ gives

$$
\begin{aligned}
0
&\leq \int (\phi^\star-\phi)(f_1-kf_0)\mkern3mu d\nu\\
&=\mathbb E_{\theta_1}[\phi^\star(X)-\phi(X)]
-k\mathbb E_{\theta_0}[\phi^\star(X)-\phi(X)].
\end{aligned}
$$

If $k>0$, exact calibration of $\phi^\star$ and the level constraint on $\phi$ imply

$$
\mathbb E_{\theta_0}[\phi^\star(X)-\phi(X)]
=\alpha-\mathbb E_{\theta_0}[\phi(X)]
\geq0.
$$

Consequently,

$$
\mathbb E_{\theta_1}[\phi^\star(X)-\phi(X)]
\geq
k\mathbb E_{\theta_0}[\phi^\star(X)-\phi(X)]
\geq0,
$$

which proves that $\phi^\star$ is most powerful. If $k=0$, the integrated pointwise inequality directly gives $\mathbb E_{\theta_1}[\phi^\star-\phi]\geq0$, so the conclusion still follows. Randomization on $\lbrace f_1=kf_0\rbrace$ causes no difficulty because the pointwise product is zero there.

The same calculation gives the usual converse and uniqueness statement. Suppose $k>0$ and another level $\alpha$ test $\widetilde\phi$ is also most powerful. Equality of alternative power forces equality throughout the display above: $\widetilde\phi$ must also use the full null budget, and the nonnegative pointwise product must vanish almost everywhere with respect to $\nu$. Hence $\widetilde\phi=\phi^\star$ away from the boundary set $\lbrace f_1=kf_0\rbrace$; on that set, different randomizations are possible provided they preserve size.

> [!IMPORTANT]
> **Board work 2 — Neyman--Pearson in the Normal model**
>
> Derive the most powerful test for
>
> $$
> \mathcal N(\mu_0,\sigma^2)
> \quad\text{versus}\quad
> \mathcal N(\mu_1,\sigma^2),
> \qquad \mu_1>\mu_0,
> $$
>
> from an i.i.d. sample. Expand the joint log likelihood ratio, reduce it to an increasing function of $\overline X_n$, and calibrate the cutoff under $\mu_0$. Explain geometrically why observations with larger $\overline X_n$ receive the type I error budget when $\mu_1>\mu_0$.

> [!TIP]
> **AI interaction 2 — Audit a Neyman--Pearson proof**
>
> Copy the prompt below into an AI interface and audit the response.

```text
Independently prove the Neyman--Pearson lemma for two simple hypotheses,
including boundary randomization and the case k = 0. Then compare your proof
with the proof in the course chapter. Identify any missing assumption or
incorrect inequality, and explain why maximizing the likelihood under each
hypothesis separately is not the same optimization problem.
```

**Audit question:** Did the response use the level constraint in the correct direction?

### 3.3 Sufficiency and the likelihood ratio bridge

Suppose $T(X)$ is sufficient and the factorization theorem gives

$$
f_\theta(x)=g_\theta(T(x))h(x).
$$

Then the Neyman--Pearson likelihood ratio satisfies

$$
\frac{f_{\theta_1}(x)}{f_{\theta_0}(x)}
=\frac{g_{\theta_1}(T(x))}{g_{\theta_0}(T(x))},
$$

so the most powerful rule can be expressed through the sufficient statistic. This explains why sufficient reductions recur in optimal testing, but sufficiency alone does not determine which values of $T$ should lead to rejection; the ordering comes from the likelihood ratio.

There is also a direct link to the generalized likelihood ratio statistic. If the maintained parameter space consists only of $\lbrace\theta_0,\theta_1\rbrace$, then

$$
\Lambda(x)
=\frac{f_{\theta_0}(x)}{\max\lbrace f_{\theta_0}(x),f_{\theta_1}(x)\rbrace}.
$$

Where $f_{\theta_1}(x)>f_{\theta_0}(x)$, the statistic satisfies $\Lambda(x)=f_{\theta_0}(x)/f_{\theta_1}(x)$, so smaller values correspond exactly to larger Neyman--Pearson ratios. In particular, a small $\Lambda$ rejection region whose cutoff is below one is an upper likelihood ratio region. The Neyman--Pearson lemma therefore supplies an exact optimality foundation for likelihood ratio ordering in the simple versus simple case. It does not by itself establish finite sample optimality for arbitrary composite likelihood ratio tests.

<a id="w4-stop-5"></a>

## 4. Monotone likelihood ratio and restricted optimality

### 4.1 Uniformly most powerful tests

A level $\alpha$ test $\phi^\star$ is **uniformly most powerful** (UMP) against a composite alternative $\Theta_1$ if

$$
\beta_{\phi^\star}(\theta)
\geq
\beta_\phi(\theta)
$$

for every $\theta\in\Theta_1$ and every other level $\alpha$ test $\phi$.

UMP tests are exceptional: the most powerful rejection region for one alternative may not be most powerful for another.

### 4.2 Monotone likelihood ratio

A one parameter family has **monotone likelihood ratio** in a statistic $T(X)$ if, whenever $\theta_1>\theta_0$,

$$
\frac{f_X(x;\theta_1)}{f_X(x;\theta_0)}
$$

is a nondecreasing function of $T(x)$.

The **Karlin--Rubin theorem** implies, under its conditions, that for

$$
H_0:\theta\leq\theta_0
\qquad\text{versus}\qquad
H_1:\theta>\theta_0,
$$

a level $\alpha$ test that rejects for large $T(X)$ is UMP. The theorem converts a collection of Neyman--Pearson problems into a common rejection region.

Regular one parameter exponential families often have monotone likelihood ratio in their natural statistic.

### 4.3 Normal model: from Neyman--Pearson to UMP

Return to i.i.d. $\mathcal N(\mu,\sigma^2)$ observations with known $\sigma^2$. For $\mu_1>\mu_0$,

$$
\log\frac{f_{\mu_1}(x)}{f_{\mu_0}(x)}
=\frac{n(\mu_1-\mu_0)}{\sigma^2}
\left(\overline x_n-\frac{\mu_1+\mu_0}{2}\right).
$$

The likelihood ratio is therefore strictly increasing in $\overline X_n$. Neyman--Pearson says that the most powerful level $\alpha$ test of $\mu=\mu_0$ against any fixed $\mu_1>\mu_0$ rejects for large $\overline X_n$. Once the cutoff is calibrated under $\mu_0$, it does not depend on the chosen $\mu_1$:

$$
\overline X_n>\mu_0+\frac{\sigma}{\sqrt n}z_{1-\alpha}.
$$

The same rule is most powerful for every $\mu_1>\mu_0$, and its rejection probability is increasing in $\mu$. Hence its size over $\mu\leq\mu_0$ is attained at the boundary and it is UMP for

$$
H_0:\mu\leq\mu_0
\qquad\text{versus}\qquad
H_1:\mu>\mu_0.
$$

This calculation displays the full bridge: sufficiency reduces the data to $\overline X_n$, monotone likelihood ratio orders its values, boundary calibration controls the composite null, and the common Neyman--Pearson region supplies uniform optimality.

### 4.4 Why two sided UMP tests usually fail

For

$$
H_0:\theta=\theta_0
\qquad\text{versus}\qquad
H_1:\theta\neq\theta_0,
$$

alternatives on opposite sides favor opposite tails. A single level $\alpha$ rejection region generally cannot be most powerful in both directions. Optimality may be recovered only after restricting the class—for example, to unbiased, invariant, or equal tail tests.

### 4.5 Unbiased tests and two sided optimality

A level $\alpha$ test with power function $\beta_\phi$ is **unbiased** if

$$
\beta_\phi(\theta)\geq\alpha
\qquad\text{for every }\theta\in\Theta_1.
$$

Thus the test rejects at least as often under every alternative as the nominal rejection probability; because its null rejection probabilities are at most $\alpha$, it also satisfies

$$
\inf_{\theta\in\Theta_1}\beta_\phi(\theta)
\geq\alpha
\geq\sup_{\theta\in\Theta_0}\beta_\phi(\theta).
$$

A **uniformly most powerful unbiased** (UMPU) level $\alpha$ test has at least as much power at every alternative as every other unbiased level $\alpha$ test. This is optimality within a restricted class, not unrestricted UMP optimality.

For i.i.d. $\mathcal N(\mu,\sigma^2)$ observations with known $\sigma$, the test of

$$
H_0:\mu=\mu_0
\qquad\text{versus}\qquad
H_1:\mu\neq\mu_0
$$

that rejects when

$$
\left|
\frac{\sqrt n(\overline X_n-\mu_0)}{\sigma}
\right|>z_{1-\alpha/2}
$$

is UMPU. The equal tail construction is not an arbitrary convention here: within the unbiased class, it balances the competing one sided alternatives. More general UMPU conclusions require their own exponential family and regularity conditions.

### Checkpoint 3

1. What does monotonicity of the likelihood ratio contribute beyond sufficiency of $T$?
2. Why does a two sided alternative create competing power objectives?
3. Does failure of a UMP test mean that no reasonable test exists?
4. How does unbiasedness restrict the class of competing tests?
5. Why does UMPU optimality not imply unrestricted UMP optimality?

<a id="w4-stop-6"></a>

## 5. Confidence sets and test inversion

Point estimation reports one data dependent value, and hypothesis testing evaluates a researcher specified division of the parameter space. Confidence set estimation instead reports a data dependent subset of parameter values. Test inversion will show that the last two procedures are dual descriptions rather than unrelated forms of inference.

### 5.1 Coverage

A confidence set is a random set $C(X)\subseteq\Theta$.

For a scalar parameter, an **interval estimator** is a random interval $I(X)=[L(X),U(X)]$, where both endpoints are statistics. The word “random” refers to the procedure before the sample is observed; the realized endpoints are fixed numbers after observing $X=x$.

The coverage function of a confidence set is

$$
c(\theta)=\mathbb{P}_\theta\lbrace\theta\in C(X)\rbrace.
$$

It has confidence coefficient at least $1-\alpha$ if

$$
\inf_{\theta\in\Theta}c(\theta)\geq1-\alpha.
$$

This is a repeated sampling statement about the random procedure $C(X)$. After observing $X=x$, the set $C(x)$ is fixed. In the frequentist model, the parameter is not assigned a posterior probability of belonging to that realized set.

### 5.2 Inverting tests

For every candidate value $\theta_0$, suppose $\phi_{\theta_0}(X)$ is a nonrandomized level $\alpha$ test of

$$
H_0:\theta=\theta_0.
$$

Define

$$
C(X)=\lbrace\theta_0:\phi_{\theta_0}(X)=0\rbrace.
$$

Then

$$
\mathbb{P}_\theta\lbrace\theta\in C(X)\rbrace
=\mathbb{P}_\theta\lbrace\phi_\theta(X)=0\rbrace
\geq1-\alpha.
$$

Thus a family of tests generates a confidence set by collecting the parameter values not rejected. Conversely, a confidence set generates tests by rejecting values outside the set.

Explicitly, if $C(X)$ has coverage at least $1-\alpha$, then for every $\theta_0$ the set

$$
A(\theta_0)=\lbrace x:\theta_0\in C(x)\rbrace
$$

is the nonrejection region of a level $\alpha$ test of $H_0:\theta=\theta_0$. The two constructions are inverse descriptions of the same event:

$$
\lbrace\theta_0\in C(X)\rbrace=\lbrace X\in A(\theta_0)\rbrace.
$$

### 5.3 Normal confidence interval

For i.i.d. $\mathcal N(\mu,\sigma^2)$ observations with known $\sigma$, invert the two sided tests

$$
\left|
\frac{\sqrt n(\overline X_n-\mu_0)}{\sigma}
\right|>z_{1-\alpha/2}.
$$

The resulting exact interval is

$$
C(X)=
\left[
\overline X_n-z_{1-\alpha/2}\frac{\sigma}{\sqrt n},
\overline X_n+z_{1-\alpha/2}\frac{\sigma}{\sqrt n}
\right].
$$

Its coverage can be checked directly:

$$
\begin{aligned}
\mathbb{P}_\mu\lbrace\mu\in C(X)\rbrace
&=\mathbb{P}_\mu\mkern-3mu\Bigl\lbrace
-z_{1-\alpha/2}
\leq \frac{\sqrt n(\overline X_n-\mu)}{\sigma}
\leq z_{1-\alpha/2}
\Bigr\rbrace\\
&=1-\alpha.
\end{aligned}
$$

With unknown variance under Gaussian sampling, replacing $\sigma$ by $S_n$ and using the exact $t_{n-1}$ distribution yields the usual Student interval.

> [!IMPORTANT]
> **Board work 3 — Confidence sets by test inversion**
>
> Invert the one sided Normal test from Section 2 to obtain a one sided confidence bound. Then show directly that the coverage event and the nonrejection event are identical.

> [!TIP]
> **AI interaction 3 — What does 95% confidence mean?**
>
> Copy the prompt below into an AI interface and audit the response.

```text
A student reports: "Given the data, there is a 95% probability that mu lies in
this confidence interval."

Evaluate the statement in the frequentist Normal model. Write the exact
coverage probability before observing the data, identify what becomes fixed
after observation, and contrast coverage with a Bayesian credible probability
without claiming that either concept is universally superior.
```

### Checkpoint 4

1. Can two confidence procedures have the same confidence coefficient but different expected length?
2. What happens to coverage when a family of pointwise tests fails to control a composite null?
3. Does a value inside a 95% confidence set necessarily have p-value $0.95$?
4. Why can inversion of asymmetric tests produce asymmetric confidence sets?

## 6. Recap

- A test is a decision rule; a test statistic orders the evidence, a critical value determines the rejection region, and the power function records rejection probability across the parameter space.
- Level controls the supremum of type I error probability over the null; size is the achieved supremum.
- A valid p-value indexes a nested family of level controlled tests; it is not a posterior null probability.
- Normal pivots give exact $z$ tests when variance is known and exact Student tests when it is unknown; the two sided Student test is UMPU rather than unrestricted UMP.
- Fixed alternative consistency and local alternative power answer different questions about a test's sensitivity.
- Neyman--Pearson allocates a type I error budget to observations with the largest likelihood ratio.
- Sufficiency can reduce the Neyman--Pearson ratio to a function of $T$, while monotone likelihood ratio determines which values of $T$ are most favorable to larger parameter values.
- Monotone likelihood ratio can produce one sided UMP tests, but two sided UMP tests often do not exist without additional restrictions.
- Unbiasedness supplies one such restriction; UMPU optimality is optimality only within the resulting class.
- Confidence sets are random procedures with repeated sampling coverage and can be constructed by inverting tests.

## Notation introduced this week

- $\Theta_0$ and $\Theta_1$: null and alternative parameter sets.
- $\phi(X)$: possibly randomized test.
- $T(X)$ and $c$: test statistic and critical value.
- $\beta_\phi(\theta)$: power function.
- $\alpha$: nominal test level.
- $p(X)$: valid p-value.
- $\Lambda(X)$: generalized likelihood ratio statistic in the simple versus simple bridge.
- UMPU: uniformly most powerful among unbiased level controlled tests.
- $z_q$: $q$ quantile of the standard Normal distribution.
- $t_{\nu,q}$: $q$ quantile of the Student $t$ distribution with $\nu$ degrees of freedom.
- $C(X)$: confidence set.
- $L(X)$ and $U(X)$: endpoints of an interval estimator.
- $c(\theta)$: coverage function.

## References

- Jun Shao, *Mathematical Statistics*, Springer, 2003.
- Robert W. Keener, *Theoretical Statistics: Topics for a Core Course*, Springer, 2010.
- Erich Lehmann and Joseph Romano, *Testing Statistical Hypotheses*, 3rd ed., Springer, 2005.
- Aad van der Vaart, *Asymptotic Statistics*, Cambridge University Press, 1998.
