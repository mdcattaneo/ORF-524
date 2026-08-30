# ORF 524 Practice Module 2: Models, Likelihood, and Sufficiency

**Status:** Ready for Fall 2026  
**Last updated:** August 30, 2026  
**Primary chapter:** [Week 2: Models, identification, likelihood, and sufficiency](../lectures/week-02.md)

## Purpose

This ungraded module asks what can be learned before choosing an estimator. Problems 1--4 form the core bank. Problems 5--7 are optional extensions on mixed distributions, nonunique maximum likelihood, and exact Gaussian sampling calculations.

<a id="problem-map"></a>

## Problem map

| Problem | Bank | Main task |
|---|---|---|
| [1. Identification in a threshold model](#problem-1) | Core | Separate observable identification from latent parameters and population likelihood. |
| [2. Normalizations in fixed effects models](#problem-2) | Core | Identify fixed effects using normalizations and contrasts. |
| [3. Uniform likelihood with moving support](#problem-3) | Core | Handle moving support in likelihood, maximum likelihood, and minimal sufficiency. |
| [4. How much does sufficiency reduce?](#problem-4) | Core | Compare the information reduction achieved by sufficient statistics across models. |
| [5. A censored observation needs a mixed measure](#problem-5) | Further | Construct a likelihood for an observation with discrete and continuous components. |
| [6. Must an MLE be a function of a sufficient statistic?](#problem-6) | Further | Qualify the relationship between maximum likelihood estimators and sufficient statistics. |
| [7. Exact Gaussian sampling calculations](#problem-7) | Further | Derive the Normal, chi squared, independence, and Student laws from Gaussian geometry. |

## How to work with this module

For each problem, identify the observable model, parameter or target, support, and dominating measure when relevant. Separate identification from estimation, keep parameter dependent support visible, and state the criterion used to establish sufficiency or minimal sufficiency.

AI may be used on these ungraded problems under the syllabus policy and the staged protocol in [`AGENTS.md`](AGENTS.md). Begin with a genuine attempt, request the least revealing useful hint, and verify every generated claim line by line. The closed book midterms remain fully unaided.

For released exam practice, use the [question level first half guide](../exams/README.md#question-level-guide-for-the-first-half). Confirm that the listed material has been covered, attempt the selected question before opening its solution, and treat historical exam instructions as superseded by the current syllabus.

## Core practice

<a id="problem-1"></a>

### Problem 1. Identification in a threshold model

**Chapter connections:** [Section 2: Identification before estimation](../lectures/week-02.md#2-identification-before-estimation) and [Section 3.2: Why population likelihood favors the true law](../lectures/week-02.md#32-why-population-likelihood-favors-the-true-law)

Let latent variables $Z_i$ be i.i.d. $\mathcal N(\mu,\sigma^2)$, with $\sigma>0$, and suppose only

$$
Y_i=\mathbf 1\lbrace Z_i\leq\kappa\rbrace
$$

is observed. The parameter is $\theta=(\kappa,\mu,\sigma)$.

1. Derive the observable distribution of $Y_i$ and show that it depends on $\theta$ only through

$$
\delta(\theta)=\frac{\kappa-\mu}{\sigma}.
$$

2. Characterize the observational equivalence class of a given $\theta$. Is the full parameter identified? Is $\delta(\theta)$ identified?

3. Determine what is identified in each case: (i) $(\mu,\sigma)$ known; (ii) only $\kappa$ known; (iii) only $\sigma$ known; and (iv) the normalization $\mu=0$, $\sigma=1$ imposed.

4. Let $p_0$ be the true success probability. Write the population Bernoulli log likelihood as a function of a candidate probability $p\in(0,1)$ and show that it is uniquely maximized at $p=p_0$. Explain why this identifies $\delta$ but not necessarily $\theta$.

<!-- Source lineage: Fall 2024 pset2, “Identification, probit” and “Probit model”; corrected and expanded to separate the identified index from the redundant parameterization and population likelihood. -->

[Back to the problem map](#problem-map)

<a id="problem-2"></a>

### Problem 2. Normalizations in fixed effects models

**Chapter connection:** [Section 2: Identification before estimation](../lectures/week-02.md#2-identification-before-estimation)

Suppose independent observations satisfy

$$
X_{it}\sim\mathcal N(\mu+\alpha_i+\beta_t,\sigma^2),
\qquad i=1,\ldots,N,\quad t=1,\ldots,T,
$$

with all components of $(\mu,\alpha_1,\ldots,\alpha_N,\beta_1,\ldots,\beta_T,\sigma^2)$ unknown.

1. Exhibit a two parameter family of transformations of $(\mu,\alpha,\beta)$ that leaves every observable distribution unchanged. Conclude that the unrestricted parameterization is not identified.

2. Show that the restrictions

$$
\sum_{i=1}^N\alpha_i=0,
\qquad
\sum_{t=1}^T\beta_t=0
$$

   identify all mean parameters. Give formulas for them in terms of the population cell means $m_{it}=\mathbb E[X_{it}]$.

3. Explain why $\sigma^2$ is identified without those location normalizations.

4. Give two nontrivial targets involving the fixed effects that are identified even in the unrestricted parameterization, and give one target that is not.

<!-- Source lineage: Fall 2024 pset2, “Identification, fixed effects models”; the two-way invariance and identified contrasts are made explicit. -->

[Back to the problem map](#problem-map)

<a id="problem-3"></a>

### Problem 3. Uniform likelihood with moving support

**Chapter connections:** [Section 3.3: Uniform model](../lectures/week-02.md#33-uniform-model-the-support-matters), [Section 5.3: Factorization theorem](../lectures/week-02.md#53-factorization-theorem), and [Section 5.4: Minimal sufficiency](../lectures/week-02.md#54-minimal-sufficiency)

Let $X_1,\ldots,X_n$ be i.i.d. $\mathsf{Uniform}(0,\theta)$, with $\theta>0$, and let $M_n=X_{(n)}=\max_iX_i$.

1. Write the complete likelihood, including its support indicator. On the event $M_n>0$, derive the MLE of $\theta$ and explain why a score equation does not locate it. What happens on the null sample $X_1=\cdots=X_n=0$ when the parameter space is $(0,\infty)$?

2. Derive the c.d.f. and density of $M_n$.

3. Use factorization to prove that $M_n$ is sufficient. Use the likelihood ratio criterion to prove that it is minimal sufficient.

4. Let $\theta_0$ be the true value. Evaluate

$$
Q(\theta)=\mathbb E_{\theta_0}[\log f(X_1;\theta)]
$$

   for $\theta<\theta_0$ and $\theta\geq\theta_0$, using the extended value $-\infty$ where appropriate. Show that $Q$ is uniquely maximized at $\theta_0$. What roles are played by support and identification?

<!-- Source lineage: Fall 2024 pset2, “Uniform distribution”; rewritten to integrate likelihood, parameter-dependent support, population likelihood, and minimal sufficiency. -->

[Back to the problem map](#problem-map)

<a id="problem-4"></a>

### Problem 4. How much does sufficiency reduce?

**Chapter connections:** [Section 4.3: Exponential families](../lectures/week-02.md#43-exponential-families) and [Section 5.4: Minimal sufficiency](../lectures/week-02.md#54-minimal-sufficiency)

For each model below, let $X_1,\ldots,X_n$ be i.i.d. Identify the family, determine whether it is a regular one parameter exponential family in the stated parameter, and find a minimal sufficient statistic.

1. $f(x;\theta)=\theta^{-1}e^{-x/\theta}\mathbf 1\lbrace x>0\rbrace$, $\theta>0$.

2. $f(x;\theta)=e^{-(x-\theta)}\mathbf 1\lbrace x>\theta\rbrace$, $\theta\in\mathbb R$.

3. $f(x;\theta)=\tfrac12e^{-|x-\theta|}$, $\theta\in\mathbb R$.

For part 3, show that the vector of order statistics is minimal sufficient. In particular, explain why there is generally no fixed dimensional reduction as $n$ grows.

<!-- Source lineage: Fall 2024 pset2, “Sufficiency”; reduced to three contrasting cases and strengthened with likelihood-ratio proofs of minimality. -->

[Back to the problem map](#problem-map)

## Further practice

<a id="problem-5"></a>

### Problem 5. A censored observation needs a mixed measure

**Chapter connection:** This extension develops the common dominating measure requirement used in [Section 3.1: Sample likelihood](../lectures/week-02.md#31-sample-likelihood).

Let $X\sim\mathcal N(\mu,\sigma^2)$, with $\sigma>0$, and observe $Y=\max\lbrace X,c\rbrace$ for known $c$.

1. Find the probability mass at $c$ and the density on $(c,\infty)$.

2. Let $\nu=\delta_c+\lambda|_{(c,\infty)}$, where $\delta_c$ is point mass and $\lambda$ is Lebesgue measure. Write one density of $Y$ with respect to $\nu$.

3. Write the likelihood for an i.i.d. sample $Y_1,\ldots,Y_n$. Explain what information would be lost by pretending that the distribution had only a Lebesgue density.

<!-- Source lineage: Fall 2024 pset2, “Left censoring”; refocused on the dominating measure and likelihood. -->

[Back to the problem map](#problem-map)

<a id="problem-6"></a>

### Problem 6. Must an MLE be a function of a sufficient statistic?

**Chapter connections:** [Section 3.1: Sample likelihood](../lectures/week-02.md#31-sample-likelihood) and [Section 5.3: Factorization theorem](../lectures/week-02.md#53-factorization-theorem)

Suppose a dominated model admits the factorization

$$
f(x;\theta)=g(T(x),\theta)h(x),
$$

with $h(x)>0$ on the common sample space.

1. Show that the set $\arg\max_\theta f(x;\theta)$ depends on $x$ only through $T(x)$.

2. If the maximizer is unique, show that the MLE is a function of $T$.

3. If maximizers are nonunique, explain why a measurable MLE can be chosen as a function of $T$ when an appropriate measurable selection exists, but why an arbitrary selection based on the full sample need not be a function of $T$.

<!-- Source lineage: legacy “MLE and sufficiency” exercises; qualifications added for nonexistence, ties, and measurable selection. -->

[Back to the problem map](#problem-map)

<a id="problem-7"></a>

### Problem 7. Exact Gaussian sampling calculations

**Chapter connections:** [Section 4.1: Finite sample statistics](../lectures/week-02.md#41-finite-sample-statistics) and [Week 4, Section 2.4: Exact Student inference with unknown variance](../lectures/week-04.md#24-exact-student-inference-with-unknown-variance)

Let $X_1,\ldots,X_n$ be i.i.d. $\mathcal N(\mu,\sigma^2)$, where $n\geq2$, and define

$$
\overline X_n=\frac1n\sum_{i=1}^nX_i,
\qquad
S_n^2=\frac1{n-1}\sum_{i=1}^n(X_i-\overline X_n)^2.
$$

Do not cite the standard Gaussian sampling results without proof in this problem.

1. Derive $\mathbb E[\overline X_n]$, $\mathbb V(\overline X_n)$, and the exact distribution of $\overline X_n$.

2. Prove the identity

$$
\sum_{i=1}^n(X_i-\overline X_n)^2
=\sum_{i=1}^n(X_i-\mu)^2-n(\overline X_n-\mu)^2,
$$

   and use it to show that $S_n^2$ is unbiased for $\sigma^2$.

3. Let $Z=(X-\mu\mathbf 1_n)/\sigma$ and choose an orthogonal matrix $Q$ whose first row is $\mathbf 1_n'/\sqrt n$. For $U=QZ$, identify $U_1$ and $\sum_{j=2}^nU_j^2$ in terms of $\overline X_n$ and $S_n^2$. Use the joint distribution of $U$ to prove

$$
\frac{(n-1)S_n^2}{\sigma^2}\sim\chi^2_{n-1}
$$

   and the independence of $\overline X_n$ and $S_n^2$.

4. Use the chi squared moments to derive $\mathbb E[S_n^2]$ again and compute $\mathbb V(S_n^2)$.

5. Derive the exact Student law

$$
\frac{\sqrt n(\overline X_n-\mu)}{S_n}\sim t_{n-1},
$$

   and use it to construct an exact equal tail confidence interval with coverage $1-\alpha$ for $\mu$.

<!-- Source lineage: Fall 2019 through Fall 2025 assignment banks, “Gaussian model”; restored as a complete finite sample calculation and corrected to prove independence through an orthogonal Gaussian decomposition. -->

[Back to the problem map](#problem-map)

## Completion check

After the core, you should be able to explain:

1. why consistency cannot repair a failure of identification;

2. why a normalization can identify parameters without changing the observable model;

3. why moving support affects likelihood calculus but not factorization; and

4. why low dimension and usefulness for estimation do not imply sufficiency.
