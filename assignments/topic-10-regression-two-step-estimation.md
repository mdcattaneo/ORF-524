# ORF 524 Practice Module 10: Regression and Two Step Estimation

**Status:** Ready for Fall 2026  
**Last updated:** August 30, 2026  
**Primary chapter:** [Week 10: Regression applications and two step estimation](../lectures/week-10.md)

## Purpose

This ungraded module applies M- and Z-estimation to regression and then makes the effect of an estimated first step explicit. Problems 1--4 form the core bank. Problems 5--6 provide exact Gaussian regression practice and a control variate interpretation of nuisance estimation.

<a id="problem-map"></a>

## Problem map

| Problem | Bank | Main task |
|---|---|---|
| [1. OLS without assuming a linear conditional mean](#problem-1) | Core | Separate the linear projection target from the conditional mean and derive robust OLS inference. |
| [2. Binary response: likelihood, NLS, and weighting](#problem-2) | Core | Compare likelihood, nonlinear least squares, and fixed or feasible weighting. |
| [3. The first step adjustment in a scalar two step estimator](#problem-3) | Core | Derive the first step adjustment using sequential and stacked equations. |
| [4. Why feasible weighting can be first order free](#problem-4) | Core | Prove when weight estimation has no first order effect. |
| [5. Exact Normal regression geometry](#problem-5) | Further | Derive exact Gaussian regression results from projection geometry. |
| [6. A mean zero moment as a control variate](#problem-6) | Further | Optimize an influence function adjustment using a mean zero moment. |

## How to work with this module

For each procedure, state the population target, conditional restrictions, sample criterion or estimating equation, derivative matrix, moment variance, and influence function. In the weighting problems, keep fixed infeasible weights, estimated and frozen feasible weights, continuously updated equations, and parameter weighted criteria distinct. In the two step problems, compare the known nuisance calculation with the adjustment induced by estimating the nuisance.

AI may be used on these ungraded problems under the syllabus policy and the staged protocol in [the assignment AI guide](AGENTS.md). Begin with a genuine attempt, request the least revealing useful hint, and verify every generated conditional expectation argument, derivative, matrix sign, and efficiency claim. The closed book midterms remain fully unaided.

For released exam practice, use the [question level Week 10 guide](../exams/README.md#question-level-guide-for-week-10). Attempt only the listed parts after confirming their prerequisites, and translate legacy weighting and expansion language into the distinctions and coordinatewise mean value argument used in the current chapter.

## Core practice

<a id="problem-1"></a>

### Problem 1. OLS without assuming a linear conditional mean

**Chapter connections:** [Section 1.1: The population target](../lectures/week-10.md#11-the-population-target) and [Section 1.2: Sample OLS and its first order distribution](../lectures/week-10.md#12-sample-ols-and-its-first-order-distribution)

Let $(Y_i,X_i)$ be i.i.d., with $Y_i\in\mathbb R$, $X_i\in\mathbb R^d$, $\mathbb E[Y_i^2+\Vert X_i\Vert^2]<\infty$, and $H=\mathbb E[X_iX_i']$ positive definite. Define

$$
\beta_0=\arg\min_\beta\mathbb E[(Y_i-X_i'\beta)^2],
\qquad
\varepsilon_i=Y_i-X_i'\beta_0.
$$

1. Derive $\beta_0=H^{-1}\mathbb E[X_iY_i]$ and prove the population risk decomposition

$$
\mathbb E[(Y_i-X_i'\beta)^2]
=\mathbb E[\varepsilon_i^2]
+(\beta-\beta_0)'H(\beta-\beta_0).
$$

2. Let $\widehat\beta_n$ be OLS and $\widehat H_n=\mathbb{E}_n[XX']$. Show that $\widehat H_n$ is nonsingular with probability approaching one and, on that event, derive the exact identity

$$
\widehat\beta_n-\beta_0
=\widehat H_n^{-1}\mathbb{E}_n[X\varepsilon],
$$

   For the distributional conclusions, additionally assume $\mathbb E[\Vert X_i\varepsilon_i\Vert^2]<\infty$. Derive its asymptotic linear representation and limiting distribution. Check the sign.

3. Construct a heteroskedasticity robust estimator of the asymptotic variance. Show how the formula simplifies under $\mathbb E[\varepsilon_i^2\mid X_i]=\sigma^2$.

4. Let $Z\sim\mathcal N(0,1)$, $X=(1,Z)'$, and $Y=Z^2+U$ with $\mathbb E[U\mid Z]=0$. Find $\beta_0$ and explain why OLS is correctly targeting the best linear predictor although $\mathbb E[Y\mid Z]$ is not linear in $Z$.

<!-- Source lineage: legacy linear-regression and multiple-linear-prediction problems; rewritten to separate the projection target from conditional-mean specification. -->

[Back to the problem map](#problem-map)

<a id="problem-2"></a>

### Problem 2. Binary response: likelihood, NLS, and weighting

**Chapter connections:** [Section 2.2: Binary response: likelihood versus squared loss](../lectures/week-10.md#22-binary-response-likelihood-versus-squared-loss) and [Section 2.3: Infeasible, feasible, and continuously updated reweighting](../lectures/week-10.md#23-infeasible-feasible-and-continuously-updated-reweighting)

Let $Y\in\lbrace0,1\rbrace$ and suppose

$$
\mathbb{P}(Y=1\mid X)=p(X'\beta_0),
$$

where $p$ is known, smooth, strictly between zero and one, and $\beta_0$ is identified. Write $p_\beta=p(X'\beta)$ and $\dot p_\beta=\dot p(X'\beta)$.

1. Derive the conditional likelihood score for $\beta$.

2. Derive the nonlinear least squares estimating equation. Under correct specification, show that its derivative and moment variance matrices are

$$
\Psi_{\beta,\mathrm{NLS}}(\beta_0)
=-\mathbb E[XX'\dot p_{\beta_0}^2],
$$

   and

$$
B_{\mathrm{NLS}}
=\mathbb E[XX'\dot p_{\beta_0}^2p_{\beta_0}(1-p_{\beta_0})].
$$

3. Let $v_0(X)=p_{\beta_0}(1-p_{\beta_0})$. Form the infeasible WNLS criterion using the fixed weight $1/v_0(X)$, derive its estimating equation, and show that its sandwich variance collapses to the inverse Bernoulli information matrix.

4. Starting from a consistent preliminary NLS estimator $\widetilde\beta_n$, construct feasible WNLS by estimating $v_0(X)$ and freezing the resulting weights in a second step criterion. Explain why this is a two step estimator and state conditions under which it is first order equivalent to infeasible WNLS.

5. Form the continuously updated estimating equation using $1/\lbrace p_\beta(1-p_\beta)\rbrace$ and show that it is exactly the likelihood score equation. Differentiate the naive criterion $\mathbb{E}_n[(Y-p_\beta)^2/\lbrace p_\beta(1-p_\beta)\rbrace]$ and explain why minimizing it is a different procedure.

6. Explain why MLE and NLS estimate the same parameter under correct specification but can estimate different pseudo true parameters under misspecification.

<!-- Source lineage: Fall 2024 pset7, “Estimation, probit,” and the legacy binary-choice slides; generalized to an arbitrary smooth link and focused on target and efficiency. -->

[Back to the problem map](#problem-map)

<a id="problem-3"></a>

### Problem 3. The first step adjustment in a scalar two step estimator

**Chapter connections:** [Section 3.2: The first step contribution](../lectures/week-10.md#32-the-first-step-contribution) and [Section 3.3: Stacking gives the same answer](../lectures/week-10.md#33-stacking-gives-the-same-answer)

Let $\beta_0$ and $\gamma_0$ be scalar and satisfy

$$
\mathbb E[m(W,\beta_0,\gamma_0)]=0,
\qquad
\mathbb E[g(W,\gamma_0)]=0.
$$

Let

$$
a=\mathbb E[\partial_\beta m(W,\beta_0,\gamma_0)],
\quad
b=\mathbb E[\partial_\gamma m(W,\beta_0,\gamma_0)],
\quad
c=\mathbb E[\partial_\gamma g(W,\gamma_0)],
$$

with $a\neq0$ and $c\neq0$. Assume consistency, finite second moments, and the local regularity conditions needed below.

1. If $\gamma_0$ were known, derive the influence function of the root $\widetilde\beta_n$ of $\mathbb{E}_n m(W,\beta,\gamma_0)=0$.

2. Let $\widehat\gamma_n$ solve $\mathbb{E}_n g(W,\gamma)=0$, and let $\widehat\beta_n$ solve $\mathbb{E}_n m(W,\beta,\widehat\gamma_n)=0$. Derive the influence function of $\widehat\beta_n$ sequentially and by stacking the equations.

3. Show that its asymptotic variance is

$$
\frac1{a^2}
\left\lbrace
\mathbb V[m]
+\frac{b^2}{c^2}\mathbb V[g]
-2\frac bc\mathrm{Cov}(m,g)
\right\rbrace,
$$

   with all moments evaluated at the truth. Construct a plug in variance estimator.

4. What happens when $b=0$? Interpret the condition as orthogonality rather than independence.

<!-- Source lineage: Fall 2024 pset7 and Fall 2025 Problem Set 4, “Two-step procedures and efficiency”; streamlined to expose the first step influence adjustment. -->

[Back to the problem map](#problem-map)

<a id="problem-4"></a>

### Problem 4. Why feasible weighting can be first order free

**Chapter connections:** [Section 4.1: Finite dimensional Neyman orthogonality](../lectures/week-10.md#41-finite-dimensional-neyman-orthogonality) and [Section 4.2: Why feasible weighting can be first order equivalent](../lectures/week-10.md#42-why-feasible-weighting-can-be-first-order-equivalent)

Suppose

$$
\mathbb E[Y-G(X'\beta_0)\mid X]=0,
\qquad
\mathbb E[\lbrace Y-G(X'\beta_0)\rbrace^2\mid X]=v(X,\gamma_0)>0.
$$

Let $q(X,\beta)=\nabla_\beta G(X'\beta)$. A preliminary estimator $\widehat\gamma_n$ is consistent at the $\sqrt n$ rate, and $\widehat\beta_n$ solves

$$
0
=\mathbb{E}_n\left[
\frac{q(X,\beta)}{v(X,\widehat\gamma_n)}
\lbrace Y-G(X'\beta)\rbrace
\right].
$$

1. Verify that the population equation is zero at $\beta_0$ for every admissible weighting function measurable with respect to $X$.

2. Differentiate the population moment with respect to $\gamma$ and prove that the derivative is zero at $(\beta_0,\gamma_0)$.

3. Derive the influence function and show that estimating $\gamma_0$ has no first order effect. Identify the asymptotic variance when the conditional variance model is correct.

4. Specialize to $G(u)=p(u)$ and Bernoulli $Y$. Explain how a preliminary NLS estimate produces the frozen feasible weights in Section 2.3 and why the resulting second step estimator is first order equivalent to infeasible WNLS and the correctly specified MLE.

<!-- Source lineage: the feasible weighted-NLS binary-choice example in the legacy Week 9 slides; rewritten as an explicit orthogonality calculation. -->

[Back to the problem map](#problem-map)

## Further practice

<a id="problem-5"></a>

### Problem 5. Exact Normal regression geometry

**Chapter connections:** [Week 10, Section 1: Best linear prediction and OLS](../lectures/week-10.md#1-best-linear-prediction-and-ols) and [Week 6, Section 2.3: Studentization](../lectures/week-06.md#23-studentization)

Condition on a nonrandom full column rank $n\times d$ matrix $X$, where $n>d$, and suppose

$$
Y=X\beta+\varepsilon,
\qquad
\varepsilon\sim\mathcal N(0,\sigma^2I_n).
$$

Let $P=X(X'X)^{-1}X'$ and $M=I_n-P$.

1. Prove that $P$ and $M$ are symmetric idempotent projections, with ranks $d$ and $n-d$.

2. Derive the conditional distributions of $PY$ and $MY$ and prove their independence.

3. Show that $\widehat\beta\mid X\sim\mathcal N(\beta,\sigma^2(X'X)^{-1})$ and $Y'MY/\sigma^2\sim\chi^2_{n-d}$, independently.

4. Construct the exact Student statistic for a scalar contrast $r'\beta$.

<!-- Source lineage: Fall 2024 pset7 and Fall 2025 Problem Set 4, “Classical Normal linear regression.” -->

[Back to the problem map](#problem-map)

<a id="problem-6"></a>

### Problem 6. A mean zero moment as a control variate

**Chapter connections:** [Section 3.2: The first step contribution](../lectures/week-10.md#32-the-first-step-contribution) and [Section 4.1: Finite dimensional Neyman orthogonality](../lectures/week-10.md#41-finite-dimensional-neyman-orthogonality)

Suppose $\beta_0$ is scalar and satisfies $\mathbb E[m(W,\beta_0)]=0$, with $a=\mathbb E[\partial_\beta m(W,\beta_0)]\neq0$. Also suppose a scalar auxiliary function $g(W)$ has known mean zero. For a constant $\lambda$, define $\widehat\beta_\lambda$ by

$$
\mathbb{E}_n\lbrace m(W,\widehat\beta_\lambda)-\lambda g(W)\rbrace=0.
$$

1. Derive the influence function and asymptotic variance of $\widehat\beta_\lambda$.

2. Find the variance minimizing $\lambda^\star$ when $\mathbb V[g]>0$.

3. Show that the optimized adjustment cannot increase asymptotic variance and identify the equality case.

4. Explain how this calculation clarifies the apparently paradoxical possibility that estimating or using an auxiliary nuisance equation can improve precision even when the nuisance value is known.

<!-- Source lineage: the efficiency comparison in the legacy two step estimation problem; reframed as a transparent control variate calculation. -->

[Back to the problem map](#problem-map)

## Completion check

After the core, you should be able to explain:

1. what OLS targets without a linear conditional mean;

2. how NLS, infeasible WNLS, feasible WNLS, continuously updated reweighting, and Bernoulli MLE are related;

3. where the first step term enters a two step influence function;

4. why consistency alone does not remove that term; and

5. how orthogonality makes feasible weighting first order equivalent to infeasible weighting.
