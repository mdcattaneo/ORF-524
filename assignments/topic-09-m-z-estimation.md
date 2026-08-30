# ORF 524 Practice Module 9: M- and Z-Estimation

**Status:** Ready for Fall 2026  
**Last updated:** August 30, 2026  
**Primary chapter:** [Week 9: M- and Z-estimation](../lectures/week-09.md)

## Purpose

This ungraded module carries the first half foundations into the second half. Problems 1--4 form the core bank: representation, consistency, asymptotic linearity, sandwich variance, and inference. Problems 5--7 extend the framework to one step estimation, misspecified likelihood, and a cumulative curved Gaussian calculation.

<a id="problem-map"></a>

## Problem map

| Problem | Bank | Main task |
|---|---|---|
| [1. The representation determines the theorem](#problem-1) | Core | Choose M and Z representations and identify their population targets. |
| [2. Uniform convergence and argmin consistency](#problem-2) | Core | Prove argmin consistency and show why pointwise convergence is insufficient. |
| [3. Nonlinear least squares and sandwich variance](#problem-3) | Core | Derive identification, asymptotic linearity, and feasible sandwich inference. |
| [4. Wald, score, and likelihood ratio in a Poisson model](#problem-4) | Core | Derive and compare three large sample likelihood tests. |
| [5. One Newton step can recover first order efficiency](#problem-5) | Further | Show how one Newton step upgrades a preliminary estimator. |
| [6. Misspecified likelihood needs the sandwich](#problem-6) | Further | Derive robust variance when a likelihood model is misspecified. |
| [7. A curved Gaussian likelihood laboratory](#problem-7) | Further | Work through moments, sufficiency, maximum likelihood, asymptotic variance, and feasible inference. |

## How to work with this module

For each problem, separate the estimator's representation, population target, identification condition, global consistency argument, and local distributional argument. Before applying an M- or Z-estimator theorem, state the precise uniform LLN, derivative or Hessian condition, and CLT being used.

AI may be used on these ungraded problems under the syllabus policy and the staged protocol in [the assignment AI guide](AGENTS.md). Begin with a genuine attempt, request the least revealing useful hint, and verify every generated regularity condition and matrix calculation. The closed book midterms remain fully unaided.

For released exam practice, use the [question level Week 9 guide](../exams/README.md#question-level-guide-for-week-9). Attempt only the listed parts after confirming their prerequisites, and translate legacy notation and proof language into the current chapter's notation and mean value/ULLN argument.

## Core practice

<a id="problem-1"></a>

### Problem 1. The representation determines the theorem

**Chapter connections:** [Section 1: Two general estimator representations](../lectures/week-09.md#1-two-general-estimator-representations) and [Section 2: Population identification](../lectures/week-09.md#2-population-identification)

For each procedure below, give an M-estimator representation, a Z-estimator representation when one is valid, and the corresponding population target. State any conditions needed to connect the two representations.

1. The sample mean.

2. The empirical $\tau$ quantile, using the check loss from Module 1.

3. Ordinary least squares for observations $W=(Y,X)$.

4. The Poisson MLE for an i.i.d. sample with mean parameter $\theta\in[0,\infty)$, using the degenerate distribution at zero as the boundary model.

Then explain why the $\mathsf{Uniform}(0,\theta)$ MLE is naturally an M-estimator but is not characterized by an ordinary interior score equation. Which representation would you choose for its consistency proof?

<!-- New bridge problem built from the Week 1 loss examples, Week 2 likelihood examples, and legacy likelihood/regression exercises. -->

[Back to the problem map](#problem-map)

<a id="problem-2"></a>

### Problem 2. Uniform convergence and argmin consistency

**Chapter connections:** [Section 2: Population identification](../lectures/week-09.md#2-population-identification) and [Section 3: Consistency](../lectures/week-09.md#3-consistency)

Let $\Theta$ be a parameter space, let $Q$ have a well separated minimizer $\theta_0$, and suppose

$$
\sup_{\theta\in\Theta}|Q_n(\theta)-Q(\theta)|\to_{\mathbb{P}}0.
$$

Let $\widehat\theta_n$ satisfy, for $r_n\geq0$,

$$
Q_n(\widehat\theta_n)
\leq\inf_{\theta\in\Theta}Q_n(\theta)+r_n,
\qquad r_n=o_{\mathbb{P}}(1).
$$

1. Prove that $\widehat\theta_n\to_{\mathbb{P}}\theta_0$ using the separation gap

$$
\eta_\varepsilon
=\inf_{\Vert\theta-\theta_0\Vert\geq\varepsilon}
\lbrace Q(\theta)-Q(\theta_0)\rbrace.
$$

2. On $\Theta=[0,1]$, let $Q(\theta)=\theta^2$ and

$$
Q_n(\theta)=\theta^2
-2\mathbf 1\mkern-3mu\left\lbrace
\left|\theta-\left(1-\frac1n\right)\right|\leq\frac1{n^2}
\right\rbrace.
$$

   Show that $Q_n(\theta)\to Q(\theta)$ for every fixed $\theta$, but minimizers of $Q_n$ do not converge to the unique minimizer of $Q$.

3. State a usable i.i.d. uniform LLN based on compactness, continuity in the parameter, and an integrable envelope. Verify its conditions for

$$
\rho((Y,X),\beta)=(Y-X'\beta)^2
$$

   on a compact set $\mathcal B\subset\lbrace\beta:\Vert\beta\Vert\leq K\rbrace$ when $\mathbb E[Y^2+\Vert X\Vert^2]<\infty$.

4. Identify the exact quantifier error in replacing uniform convergence by pointwise convergence in part 1. Why does allowing $r_n=o_{\mathbb{P}}(1)$ not create the same problem?

5. State the parallel consistency theorem for a Z-estimator based on $\sup_\theta\Vert\Psi_n(\theta)-\Psi(\theta)\Vert$ and a well separated population zero.

<!-- New bridge problem based on the Week 9 argmin theorem and its moving-spike counterexample; fills a gap in the legacy application-focused sets. -->

[Back to the problem map](#problem-map)

<a id="problem-3"></a>

### Problem 3. Nonlinear least squares and sandwich variance

**Chapter connections:** [Section 4.1: Smooth Z-estimator mean value argument](../lectures/week-09.md#41-smooth-z-estimator-mean-value-argument), [Section 5.1: Estimating the variance](../lectures/week-09.md#51-estimating-the-variance), and [Section 5.3: Best linear predictor and least squares](../lectures/week-09.md#53-best-linear-predictor-and-least-squares)

Let $(Y_i,X_i)$ be i.i.d., with $X_i\in\mathbb R^d$, and suppose

$$
\mathbb E[Y_i\mid X_i]=\mu(X_i'\beta_0),
\qquad
\mathbb V(Y_i\mid X_i)=\sigma^2(X_i),
$$

where the link $\mu$ is known and twice continuously differentiable. Define

$$
\widehat\beta_n
\in\arg\min_{\beta\in\Theta}
\frac1n\sum_{i=1}^n\lbrace Y_i-\mu(X_i'\beta)\rbrace^2.
$$

1. Decompose the population criterion into irreducible conditional variance plus approximation error. Give a condition that identifies $\beta_0$ as its unique minimizer.

2. Derive the first order estimating equation and express $\widehat\beta_n$ as a Z-estimator with

$$
\psi((Y,X),\beta)
=X\dot\mu(X'\beta)\lbrace Y-\mu(X'\beta)\rbrace.
$$

3. Assuming consistency and the local regularity conditions of Week 9, derive

$$
\sqrt n(\widehat\beta_n-\beta_0)
\rightsquigarrow\mathcal N(0,H^{-1}BH^{-1}),
$$

   identifying $H$ and $B$ explicitly under heteroskedasticity.

4. Construct a consistent plug in sandwich estimator. Show how the asymptotic variance simplifies under $\sigma^2(X)=\sigma^2$.

<!-- Source lineage: Fall 2024 pset6, “Non-linear least squares”; shortened to identification, Z-representation, asymptotic normality, and robust variance. -->

[Back to the problem map](#problem-map)

<a id="problem-4"></a>

### Problem 4. Wald, score, and likelihood ratio in a Poisson model

**Chapter connections:** [Section 4: Asymptotic linearity and normality](../lectures/week-09.md#4-asymptotic-linearity-and-normality), [Section 5.1: Estimating the variance](../lectures/week-09.md#51-estimating-the-variance), and [Section 6: Wald, score, and likelihood ratio inference](../lectures/week-09.md#6-wald-score-and-likelihood-ratio-inference)

Let $X_1,\ldots,X_n$ be i.i.d. $\mathsf{Poisson}(\theta)$ with the extended parameter space $\theta\in[0,\infty)$, and consider $H_0:\theta=\theta_0$ for a fixed $\theta_0>0$.

1. Derive the MLE, its exact asymptotic linear representation, its asymptotic distribution, and a consistent estimator of its asymptotic variance.

2. On the event $\overline X_n>0$, derive the scalar statistics

$$
W_n=\frac{n(\overline X_n-\theta_0)^2}{\overline X_n},
\qquad
LM_n=\frac{n(\overline X_n-\theta_0)^2}{\theta_0},
$$

   and

$$
LR_n=2n\left\lbrace
\overline X_n\log\frac{\overline X_n}{\theta_0}
-(\overline X_n-\theta_0)
\right\rbrace,
$$

   with $0\log0=0$. Define $W_n$ arbitrarily when $\overline X_n=0$ and explain why this event is asymptotically negligible under $H_0$.

3. Under $H_0$, use a local expansion to show that all three equal

$$
\frac{n(\overline X_n-\theta_0)^2}{\theta_0}+o_{\mathbb{P}}(1)
$$

   and converge to $\chi_1^2$.

4. Explain why first order equivalence does not imply numerical equality in finite samples. Identify where the unrestricted estimate, restricted estimate, and criterion values enter the three tests.

<!-- Source lineage: Fall 2024 pset6, “MLE inference”; one regular model selected to expose the common quadratic approximation cleanly. -->

[Back to the problem map](#problem-map)

## Further practice

<a id="problem-5"></a>

### Problem 5. One Newton step can recover first order efficiency

**Chapter connection:** [Section 4.1: Smooth Z-estimator mean value argument](../lectures/week-09.md#41-smooth-z-estimator-mean-value-argument)

Let $\Psi_n(\theta)=\mathbb{E}_n\psi(W,\theta)$ be scalar, let $\theta_0$ be its identified zero, and suppose a preliminary estimator satisfies $\sqrt n(\overline\theta_n-\theta_0)=O_{\mathbb{P}}(1)$. Define

$$
\widetilde\theta_n
=\overline\theta_n
-\Psi_{\theta,n}(\overline\theta_n)^{-1}\Psi_n(\overline\theta_n),
$$

where $\Psi_{\theta,n}(\theta)=\partial\Psi_n(\theta)/\partial\theta$. Define $\Psi_\theta(\theta)=\partial\Psi(\theta)/\partial\theta$. Under uniform local derivative convergence to a nonzero $\Psi_\theta(\theta_0)$ and a CLT for $\sqrt n\Psi_n(\theta_0)$, derive

$$
\sqrt n(\widetilde\theta_n-\theta_0)
=-[\Psi_\theta(\theta_0)]^{-1}\sqrt n\Psi_n(\theta_0)+o_{\mathbb{P}}(1).
$$

Explain why $\sqrt n$ consistency of the preliminary estimator matters and why pointwise derivative convergence is insufficient at the random intermediate value.

<!-- Source lineage: Fall 2024 pset6, “MLE updating”; generalized from likelihood to a scalar Z-estimator. -->

[Back to the problem map](#problem-map)

<a id="problem-6"></a>

### Problem 6. Misspecified likelihood needs the sandwich

**Chapter connections:** [Section 5.1: Estimating the variance](../lectures/week-09.md#51-estimating-the-variance) and [Section 5.4: Regular maximum likelihood](../lectures/week-09.md#54-regular-maximum-likelihood)

Let $W_1,\ldots,W_n$ be i.i.d. nonnegative integer valued variables with mean $\mu_0>0$ and variance $\tau^2<\infty$. A researcher fits a Poisson likelihood even though the true distribution need not be Poisson.

1. Show that the pseudo true population maximizer is $\mu_0$ and the Poisson quasi-MLE is $\overline W_n$.

2. Using the Poisson score $\psi(W,\theta)=W/\theta-1$, compute $\Psi_\theta(\mu_0)$ and $B$ and recover asymptotic variance $\tau^2$.

3. Show that the information equality holds exactly when $\tau^2=\mu_0$. How does imposing Poisson information distort the standard error under overdispersion?

<!-- Source lineage: legacy likelihood-inference problems and Week 9's misspecification discussion; new sandwich audit. -->

[Back to the problem map](#problem-map)

<a id="problem-7"></a>

### Problem 7. A curved Gaussian likelihood laboratory

**Chapter connections:** [Section 3: Consistency](../lectures/week-09.md#3-consistency), [Section 4: Asymptotic linearity and normality](../lectures/week-09.md#4-asymptotic-linearity-and-normality), [Section 5.4: Regular maximum likelihood](../lectures/week-09.md#54-regular-maximum-likelihood), and [Section 6: Wald, score, and likelihood ratio inference](../lectures/week-09.md#6-wald-score-and-likelihood-ratio-inference)

Let $X_1,\ldots,X_n$ be i.i.d. $\mathcal N(\theta,\theta^2)$ with $\theta>0$, and define

$$
S_n=\mathbb E_n[X],
\qquad
T_n=\mathbb E_n[X^2].
$$

1. Write $X=\theta(1+Z)$ for $Z\sim\mathcal N(0,1)$ and derive $\mathbb E[X^r]$ for $r=1,2,3,4$.

2. Write the likelihood and prove that $(S_n,T_n)$ is sufficient. Use the likelihood ratio criterion to determine whether it is minimal sufficient.

3. On the probability one event $T_n>0$, derive the maximum likelihood estimator

$$
\widehat\theta_n
=\frac{\sqrt{S_n^2+4T_n}-S_n}{2}.
$$

   Verify that this stationary point is the unique global maximum over $\theta>0$. What happens on the null sample $T_n=0$?

4. Prove consistency using a law of large numbers and continuous mapping.

5. Derive the joint central limit theorem for $(S_n,T_n)$, including every covariance calculation, and use the delta method to prove

$$
\sqrt n(\widehat\theta_n-\theta)
\rightsquigarrow
\mathcal N\mkern-3mu\left(0,\frac{\theta^2}{3}\right).
$$

6. Derive the score and Fisher information and verify that the asymptotic variance in part 5 attains the information bound.

7. Construct a feasible standard error, a direct Wald interval, and a positivity preserving interval obtained on the $\log\theta$ scale. State the pointwise asymptotic coverage conclusion for each interval.

<!-- Source lineage: Fall 2019 through Fall 2022 assignment banks, “Curved Gaussian distribution”; restored as a cumulative likelihood laboratory and corrected to use E[X^3]=4 theta^3. -->

[Back to the problem map](#problem-map)

## Completion check

After the core, you should be able to separate:

1. representation from identification;

2. a usable uniform LLN from an abstract uniform convergence assumption;

3. global consistency conditions from local distributional conditions;

4. the derivative matrix $\Psi_\theta(\theta_0)$ from the moment variance $B$; and

5. first order equivalence of tests from finite sample equality.
