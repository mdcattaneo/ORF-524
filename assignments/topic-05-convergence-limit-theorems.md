# ORF 524 Practice Module 5: Convergence and Limit Theorems

**Status:** Ready for Fall 2026  
**Last updated:** August 30, 2026  
**Primary chapter:** [Week 5: Convergence and limit theorems](../lectures/week-05.md)

## Purpose

This ungraded module develops the probability tools behind large sample statistics. Problems 1--4 form the core bank. Problems 5--6 provide additional practice with mapping theorems and with failures of standard LLNs and CLTs.

<a id="problem-map"></a>

## Problem map

| Problem | Bank | Main task |
|---|---|---|
| [1. Convergence modes and stochastic orders](#problem-1) | Core | Compare convergence modes, prove implications, and manipulate stochastic orders. |
| [2. Laws of large numbers and Studentization](#problem-2) | Core | Establish laws of large numbers, sample variance consistency, and Studentization. |
| [3. Triangular array and multivariate CLTs](#problem-3) | Core | Verify central limit theorems beyond the scalar identically distributed setting. |
| [4. A non Gaussian boundary limit](#problem-4) | Core | Derive a non Gaussian limit for an extreme order statistic. |
| [5. Continuous mapping and Slutsky under stress](#problem-5) | Further | Test mapping and Slutsky arguments under discontinuity and missing joint convergence. |
| [6. Diagnose failed limit theorems](#problem-6) | Further | Diagnose failures caused by dependence, heavy tails, and dominant summands. |

## How to work with this module

For every proposed limit, name the convergence mode, centering, normalization, limiting object, and theorem before calculating. Verify the theorem's dependence and moment assumptions rather than citing “the LLN” or “the CLT” generically, and distinguish a stochastic order statement from a distributional approximation.

AI may be used on these ungraded problems under the syllabus policy and the staged protocol in [`AGENTS.md`](AGENTS.md). Begin with a genuine attempt, ask for a theorem level hint before requesting a proof, and check every generated counterexample and limit directly from its definition. The closed book midterms remain fully unaided.

For released exam practice, use the [question level first half guide](../exams/README.md#question-level-guide-for-the-first-half). Confirm that the listed material has been covered, attempt the selected question before opening its solution, and treat historical exam instructions as superseded by the current syllabus.

## Core practice

<a id="problem-1"></a>

### Problem 1. Convergence modes and stochastic orders

**Chapter connections:** [Section 1: Modes of convergence](../lectures/week-05.md#1-modes-of-convergence) and [Section 3: Stochastic orders](../lectures/week-05.md#3-stochastic-orders)

1. Prove that $X_n\overset{L^p}{\longrightarrow}X$ implies $X_n\to_{\mathbb{P}}X$ for every $p>0$.

2. Prove that $X_n\rightsquigarrow X$ implies $X_n=O_{\mathbb{P}}(1)$.

3. Give fully specified counterexamples showing that:

   - convergence in distribution need not imply convergence in probability; and
   - convergence in probability need not imply convergence of expectations.

4. For deterministic positive sequences $a_n,b_n$, prove

$$
A_n=O_{\mathbb{P}}(a_n),\quad B_n=o_{\mathbb{P}}(b_n)
\quad\Longrightarrow\quad
A_nB_n=o_{\mathbb{P}}(a_nb_n),
$$

   and

$$
A_n=O_{\mathbb{P}}(a_n),\quad C_n=O_{\mathbb{P}}(b_n)
\quad\Longrightarrow\quad
A_n+C_n=O_{\mathbb{P}}(a_n+b_n).
$$

<!-- Source lineage: Fall 2024 pset5, “Modes of convergence”; extended with explicit stochastic-order algebra and two false converses. -->

[Back to the problem map](#problem-map)

<a id="problem-2"></a>

### Problem 2. Laws of large numbers and Studentization

**Chapter connections:** [Sections 4.1--4.2: Classical LLN, CLT, and Studentization](../lectures/week-05.md#41-laws-of-large-numbers) and [Section 4.4: Beyond i.i.d. laws of large numbers](../lectures/week-05.md#44-beyond-iid-laws-of-large-numbers)

1. Let $X_1,X_2,\ldots$ be independent with means $\mu_i$ and finite variances. Suppose

$$
\frac1n\sum_{i=1}^n\mu_i\to\mu,
\qquad
\frac1{n^2}\sum_{i=1}^n\mathbb V[X_i]\to0.
$$

   Prove that $n^{-1}\sum_{i=1}^nX_i\to_{\mathbb{P}}\mu$.

2. Suppose instead that the variables have common mean $\mu$, common finite variance, and

$$
\mathrm{Cov}(X_i,X_j)=\rho(|i-j|),
\qquad |\rho(h)|\to0.
$$

   Prove the same conclusion without assuming independence.

3. For each $n$, let $X_{n1},\ldots,X_{nk_n}$ be independent with finite variances. Let $v_n>0$. Show that

$$
\frac1{v_n^2}\sum_{i=1}^{k_n}\mathbb V[X_{ni}]\to0
$$

   implies

$$
\frac1{v_n}\sum_{i=1}^{k_n}
\lbrace X_{ni}-\mathbb E[X_{ni}]\rbrace
\to_{\mathbb{P}}0.
$$

4. Let $X_1,X_2,\ldots$ be i.i.d. with mean $\mu$ and variance $0<\sigma^2<\infty$, and define

$$
S_n^2=\frac1{n-1}\sum_{i=1}^n(X_i-\overline X_n)^2.
$$

   Derive the sample variance decomposition using explicit sums, prove $S_n^2\to_{\mathbb{P}}\sigma^2$, and combine this result with the classical CLT to prove

$$
\frac{\sqrt n(\overline X_n-\mu)}{S_n}
\rightsquigarrow\mathcal N(0,1).
$$

5. For each part, identify exactly where independence, covariance decay, finite second moments, or the chosen normalization enters the proof.

<!-- Source lineage: Fall 2024 pset5, “Non-iid LLNs”; reorganized to separate LLN conditions from CLT conditions and augmented with the sample-variance/Studentization bridge. -->

[Back to the problem map](#problem-map)

<a id="problem-3"></a>

### Problem 3. Triangular array and multivariate CLTs

**Chapter connections:** [Section 4.3: Cramer--Wold and the multivariate CLT](../lectures/week-05.md#43-cramer--wold-and-the-multivariate-clt) and [Section 4.5: Triangular array CLT](../lectures/week-05.md#45-reference-triangular-array-central-limit-theorem)

1. For each $n$, let $X_{n1},\ldots,X_{nn}$ be independent centered Bernoulli variables

$$
X_{ni}=B_{ni}-p_{ni},
\qquad B_{ni}\sim\mathsf{Bernoulli}(p_{ni}),
$$

   with $p_{ni}\in[\varepsilon,1-\varepsilon]$ for a fixed $\varepsilon\in(0,1/2)$. Let $s_n^2=\sum_i\mathbb V[X_{ni}]$. Verify the Lindeberg condition and conclude that

$$
\frac{\sum_iX_{ni}}{s_n}\rightsquigarrow\mathcal N(0,1).
$$

2. Explain why bounded summands alone would not verify Lindeberg if $s_n$ failed to diverge.

3. Let $Y_1,Y_2,\ldots$ be i.i.d. vectors in $\mathbb R^d$ with finite second moments, mean $\mu$, and covariance matrix $\Sigma$. Use the scalar CLT and Cramer--Wold to prove

$$
\sqrt n(\overline Y_n-\mu)\rightsquigarrow\mathcal N_d(0,\Sigma),
$$

   allowing $\Sigma$ to be singular.

4. Explain how the LLN and CLT conclusions differ in centering, normalization, and information supplied about the approximation error.

<!-- Source lineage: Fall 2024 pset5 and the legacy multivariate-limit material; retains the targeted Lindeberg verification. -->

[Back to the problem map](#problem-map)

<a id="problem-4"></a>

### Problem 4. A non Gaussian boundary limit

**Chapter connection:** [Section 4.6: A non Gaussian order statistic limit](../lectures/week-05.md#46-a-non-gaussian-order-statistic-limit)

Let $X_1,\ldots,X_n$ be i.i.d. $\mathsf{Uniform}(0,\theta)$, where $\theta>0$, and let $M_n=\max_iX_i$.

1. Derive the distribution function of $M_n$ and prove $M_n\to_{\mathbb{P}}\theta$.

2. Show that

$$
\frac{n(\theta-M_n)}{\theta}
\rightsquigarrow \mathsf{Exponential}(1).
$$

3. Deduce that $\theta-M_n=O_{\mathbb{P}}(n^{-1})$. Explain why a Normal approximation scaled by $\sqrt n$ is inappropriate.

4. For $\alpha\in(0,1)$, show that

$$
C_n(X)=\left[M_n,\frac{M_n}{\alpha^{1/n}}\right]
$$

   has exact coverage $1-\alpha$ for $\theta$. Contrast this exact pivot calculation with the limiting distribution calculation in part 2.

<!-- Source lineage: legacy Uniform-model likelihood and extreme-value exercises; moved here to contrast Gaussian and non-Gaussian limits. -->

[Back to the problem map](#problem-map)

## Further practice

<a id="problem-5"></a>

### Problem 5. Continuous mapping and Slutsky under stress

**Chapter connection:** [Section 2: Continuous mapping and Slutsky](../lectures/week-05.md#2-continuous-mapping-and-slutsky)

1. Suppose $X_n\rightsquigarrow X$ and $Y_n\to_{\mathbb{P}}c$. Derive the limits of $X_n+Y_n$, $X_nY_n$, and, when $c\neq0$, $X_n/Y_n$.

2. Let $Z\sim\mathcal N(0,1)$ and set $X_n=Z/n$. Give a function $g$ discontinuous at zero for which $g(X_n)$ does not converge in probability to $g(0)$.

3. Construct a sequence of pairs $(U_n,V_n)$ such that both marginal distributions are $\mathcal N(0,1)$ for every $n$, but $U_n+V_n$ has no limiting distribution. Explain why marginal convergence is insufficient.

<!-- Source lineage: Week 5 mapping-theorem checkpoints; expanded into explicit counterexamples. -->

[Back to the problem map](#problem-map)

<a id="problem-6"></a>

### Problem 6. Diagnose failed limit theorems

**Chapter connections:** [Section 4.4: Beyond i.i.d. laws of large numbers](../lectures/week-05.md#44-beyond-iid-laws-of-large-numbers) and [Section 4.5: Triangular array CLT](../lectures/week-05.md#45-reference-triangular-array-central-limit-theorem)

For each construction, determine whether the stated LLN or CLT conclusion holds and identify the failed assumption.

1. Let $X_i=Z+\varepsilon_i$, where $Z$ is nondegenerate with mean zero and the $\varepsilon_i$ are i.i.d. mean zero variables independent of $Z$. Does $\overline X_n\to_{\mathbb{P}}0$?

2. Let the $X_i$ be i.i.d. standard Cauchy. Does $\overline X_n\to_{\mathbb{P}}0$?

3. In row $n$, let $X_{n1}$ be Rademacher and $X_{ni}=0$ for $i>1$. Does the standardized row sum satisfy a Normal CLT?

<!-- Source lineage: legacy probability-limit counterexamples; selected to expose dependence, tail, and no-dominant-term requirements. -->

[Back to the problem map](#problem-map)

## Completion check

After the core, you should be able to:

1. name the convergence mode and normalization in any proposed limit;

2. distinguish stochastic magnitude from limiting shape;

3. verify the variance calculation behind an LLN and the consistency calculation behind Studentization;

4. verify a Lindeberg condition rather than merely naming it;

5. use Cramer--Wold for multivariate limits; and

6. recognize a valid approximation with a rate other than $\sqrt n$ and a non-Gaussian limit.
