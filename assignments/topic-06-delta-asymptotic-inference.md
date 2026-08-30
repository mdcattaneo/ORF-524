# ORF 524 Practice Module 6: Delta Method and Asymptotic Inference

**Status:** Ready for Fall 2026  
**Last updated:** August 30, 2026  
**Primary chapter:** [Week 6: Delta method and asymptotic inference](../lectures/week-06.md)

## Purpose

This ungraded module turns the limit theorems from Week 5 into complete statistical approximations. Problems 1--4 form the core bank. Problems 5--7 provide variance stabilization, second order limits, and higher order bias and MSE calculations. The chapter's uniform LLN toolkit is applied directly in [Practice Module 9](topic-09-m-z-estimation.md), where uniform convergence supports estimator consistency.

<a id="problem-map"></a>

## Problem map

| Problem | Bank | Main task |
|---|---|---|
| [1. First and second order delta methods](#problem-1) | Core | Apply delta methods at regular points and where the first derivative vanishes. |
| [2. A complete asymptotic chain for sample variance](#problem-2) | Core | Derive consistency, asymptotic linearity, a central limit theorem, and Studentization. |
| [3. Pointwise level without uniform level](#problem-3) | Core | Distinguish pointwise asymptotic validity from uniform validity. |
| [4. A transformed, studentized confidence interval](#problem-4) | Core | Build a transformed feasible interval and map it back to the original parameter. |
| [5. Variance stabilization for a Bernoulli proportion](#problem-5) | Further | Derive and assess a variance stabilizing transformation. |
| [6. A general second order delta method](#problem-6) | Further | Prove the multivariate second order delta method. |
| [7. Higher order bias and MSE calculations](#problem-7) | Further | Derive a full transformed mean expansion through order $n^{-2}$. |

## How to work with this module

For each problem, identify the centering, rate, limiting distribution, and assumptions before applying a transformation. Construct feasible inference as a full chain: establish the limit, identify its variance, prove variance consistency, Studentize, and classify the resulting validity claim as exact, pointwise asymptotic, or uniform asymptotic.

AI may be used on these ungraded problems under the syllabus policy and the staged protocol in [the assignment AI guide](AGENTS.md). Begin with a genuine attempt, request the least revealing useful hint, and verify every generated expansion, variance calculation, and coverage statement. The closed book midterms remain fully unaided.

For released exam practice, use the [question level first half guide](../exams/README.md#question-level-guide-for-the-first-half). Confirm that the listed material has been covered, attempt the selected question before opening its solution, and treat historical exam instructions as superseded by the current syllabus.

## Core practice

<a id="problem-1"></a>

### Problem 1. First and second order delta methods

**Chapter connections:** [Section 1: Delta method](../lectures/week-06.md#1-delta-method) and [Second order delta method](../lectures/week-06.md#second-order-delta-method)

Let $X_1,\ldots,X_n$ be i.i.d. $\mathsf{Bernoulli}(p)$, $p\in(0,1)$, let $\widehat p_n=\overline X_n$, and define $g(p)=p(1-p)$.

1. If $p\neq1/2$, derive the limiting distribution of

$$
\sqrt n\lbrace g(\widehat p_n)-g(p)\rbrace.
$$

2. If $p=1/2$, show that the limit under $\sqrt n$ scaling is degenerate and derive the nondegenerate limit of

$$
n\lbrace g(\widehat p_n)-g(1/2)\rbrace.
$$

3. For fixed $p\neq1/2$, construct a consistent plug in estimator of the first order asymptotic variance and an asymptotic confidence interval for $g(p)$.

4. Explain why the first order approximation in part 3 is pointwise and deteriorates as $p$ approaches $1/2$.

<!-- Source lineage: Fall 2024 pset5, “Delta method”; retains the zero-derivative boundary and adds feasible inference plus a uniformity audit. -->

[Back to the problem map](#problem-map)

<a id="problem-2"></a>

### Problem 2. A complete asymptotic chain for sample variance

**Chapter connections:** [Section 2.1: Consistency and rate](../lectures/week-06.md#21-consistency-and-rate), [Bias that disappears at first order](../lectures/week-06.md#example-bias-that-disappears-at-first-order), and [Section 2.3: Studentization](../lectures/week-06.md#23-studentization)

Let $X_1,\ldots,X_n$ be i.i.d. with $\mathbb E[|X_1|^4]<\infty$, mean $\mu$, and variance $\sigma^2>0$. Define

$$
\widehat\sigma_n^2
=\frac1n\sum_{i=1}^n(X_i-\overline X_n)^2,
\qquad
S_n^2=\frac1{n-1}\sum_{i=1}^n(X_i-\overline X_n)^2.
$$

1. Prove that both estimators are consistent for $\sigma^2$.

2. Establish the asymptotic linear representation

$$
\sqrt n(\widehat\sigma_n^2-\sigma^2)
=\frac1{\sqrt n}\sum_{i=1}^n
\lbrace(X_i-\mu)^2-\sigma^2\rbrace+o_{\mathbb{P}}(1).
$$

3. Derive the limiting distribution and identify its variance

$$
V=\mathbb E[(X_1-\mu)^4]-\sigma^4.
$$

4. Construct a consistent estimator $\widehat V_n$ of $V$ from centered sample moments.

5. Show that $S_n^2$ and $\widehat\sigma_n^2$ are first order asymptotically equivalent. If $V>0$, construct a studentized confidence interval for $\sigma^2$ and classify its coverage as exact, pointwise asymptotic, or uniform asymptotic.

<!-- Source lineage: Fall 2024 pset5, “Sample variance”; preserves consistency, asymptotic linearity, standard errors, and inference. -->

[Back to the problem map](#problem-map)

<a id="problem-3"></a>

### Problem 3. Pointwise level without uniform level

**Chapter connections:** [Section 2.4: Asymptotic size and power](../lectures/week-06.md#24-asymptotic-size-and-power) and [Section 2.5: Pointwise versus uniform coverage](../lectures/week-06.md#25-pointwise-versus-uniform-coverage)

For each $n$, consider the statistical experiment

$$
X_n\sim\mathsf{Bernoulli}(e^{-n\theta}),
\qquad \theta\in(0,1],
$$

and the test $\phi_n=X_n$ of the null consisting of the entire parameter space.

1. Show that for every fixed $\theta>0$, $\mathbb E_\theta[\phi_n]\to0$.

2. Show that $\sup_{\theta\in(0,1]}\mathbb E_\theta[\phi_n]=1$ for every $n$.

3. Identify a sequence $\theta_n$ along which rejection probability does not approach zero.

4. State precisely the difference between pointwise asymptotic level and uniform asymptotic level. Which one does this test possess at level zero?

<!-- New bridge problem motivated by Week 6's pointwise-versus-uniform validity discussion. -->

[Back to the problem map](#problem-map)

<a id="problem-4"></a>

### Problem 4. A transformed, studentized confidence interval

**Chapter connections:** [Section 1: Delta method](../lectures/week-06.md#1-delta-method) and [Section 2.3: Studentization](../lectures/week-06.md#23-studentization)

Let $X_1,\ldots,X_n$ be i.i.d. with mean $\mu>0$, variance $0<\sigma^2<\infty$, and consistent sample variance $S_n^2$. Define $\widehat\eta_n=\log\overline X_n$ on $\lbrace\overline X_n>0\rbrace$ and extend it arbitrarily on the asymptotically negligible complement.

1. Derive the limiting distribution of $\sqrt n(\widehat\eta_n-\log\mu)$.

2. Construct a consistent estimator of its asymptotic variance.

3. Form an asymptotic confidence interval for $\log\mu$ and transform it into an interval for $\mu$.

4. State the coverage conclusion carefully. Why need the interval not be symmetric on the original $\mu$ scale?

<!-- Source lineage: legacy variance-transformation and studentization exercises; rewritten as a complete approximation chain. -->

[Back to the problem map](#problem-map)

## Further practice

<a id="problem-5"></a>

### Problem 5. Variance stabilization for a Bernoulli proportion

**Chapter connections:** [Section 1: Delta method](../lectures/week-06.md#1-delta-method) and [Section 2.5: Pointwise versus uniform coverage](../lectures/week-06.md#25-pointwise-versus-uniform-coverage)

Let $g(p)=2\arcsin\sqrt p$ for $p\in(0,1)$.

1. Use the delta method to show that

$$
\sqrt n\lbrace g(\widehat p_n)-g(p)\rbrace\rightsquigarrow\mathcal N(0,1).
$$

2. Construct a transformed asymptotic interval for $p$ by first forming an interval on the $g$ scale, intersecting it with $[0,\pi]$, and then applying $g^{-1}$.

3. Explain why the argument is pointwise for interior $p$ and what fails near $0$ or $1$.

<!-- Source lineage: legacy variance-stabilizing transformation exercises. -->

[Back to the problem map](#problem-map)

<a id="problem-6"></a>

### Problem 6. A general second order delta method

**Chapter connection:** [Second order delta method](../lectures/week-06.md#second-order-delta-method)

Suppose

$$
\sqrt n(T_n-\theta)\rightsquigarrow\mathcal N(0,V),
\qquad V>0,
$$

and $g$ is twice continuously differentiable near $\theta$, with $g'(\theta)=0$ and $g''(\theta)\neq0$.

1. Show that $\sqrt n\lbrace g(T_n)-g(\theta)\rbrace\to_{\mathbb{P}}0$.

2. Under a second order Taylor remainder condition, derive the limit of $n\lbrace g(T_n)-g(\theta)\rbrace$.

3. Explain why the resulting limit is generally neither centered nor Normal.

4. Recover Problem 1, part 2, as a special case.

<!-- Source lineage: the second-order delta-method extension in the redesigned Week 6 chapter. -->

[Back to the problem map](#problem-map)

<a id="problem-7"></a>

### Problem 7. Higher order bias and MSE calculations

**Chapter connections:** [Section 1: Delta method](../lectures/week-06.md#1-delta-method), [Second order delta method](../lectures/week-06.md#second-order-delta-method), and [Section 2.2: Leading asymptotic bias, variance, and efficiency](../lectures/week-06.md#22-leading-asymptotic-bias-variance-and-efficiency)

Let $X_1,\ldots,X_n$ be i.i.d. with mean $\mu$, central moments

$$
m_r=\mathbb E[(X_1-\mu)^r],
$$

and $m_2=\sigma^2>0$. Assume $\mathbb E[|X_1-\mu|^8]<\infty$. Let $g:\mathbb R\to\mathbb R$ be five times continuously differentiable with bounded fourth and fifth derivatives, and let $\delta_n=\overline X_n-\mu$.

1. Derive the exact formulas

$$
\mathbb E[\delta_n^2]=\frac{m_2}{n},
\qquad
\mathbb E[\delta_n^3]=\frac{m_3}{n^2},
$$

and

$$
\mathbb E[\delta_n^4]
=\frac{3m_2^2}{n^2}+\frac{m_4-3m_2^2}{n^3}.
$$

2. Use a fourth order Taylor expansion and justify its expected remainder to prove

$$
\begin{aligned}
\mathbb E[g(\overline X_n)]
={}&g(\mu)+\frac{g^{(2)}(\mu)m_2}{2n}\\
&+\frac1{n^2}\left\lbrace
\frac{g^{(3)}(\mu)m_3}{6}
+\frac{g^{(4)}(\mu)m_2^2}{8}
\right\rbrace
+o(n^{-2}).
\end{aligned}
$$

3. Derive the squared bias through order $n^{-2}$.

4. Expand the squared estimation error directly and prove

$$
\begin{aligned}
\mathbb E[\lbrace g(\overline X_n)-g(\mu)\rbrace^2]
={}&\frac{g^{(1)}(\mu)^2m_2}{n}\\
&+\frac1{n^2}\left\lbrace
g^{(1)}(\mu)g^{(2)}(\mu)m_3
+g^{(1)}(\mu)g^{(3)}(\mu)m_2^2
+\frac34g^{(2)}(\mu)^2m_2^2
\right\rbrace
+o(n^{-2}).
\end{aligned}
$$

5. Compare the regular case $g^{(1)}(\mu)\neq0$ with the case $g^{(1)}(\mu)=0$ and $g^{(2)}(\mu)\neq0$. Reconcile the leading MSE in the second case with the second order delta method limit.

6. Specialize every expansion to $g(x)=x^2$. Give the exact MSE of $\overline X_n^2$ as an estimator of $\mu^2$, and interpret the special case $\mu=0$.

<!-- Source lineage: Fall 2019 through Fall 2023 assignment banks, “Asymptotic MSE”; restored with sufficient remainder conditions and corrected to retain every term of order n^{-2}. -->

[Back to the problem map](#problem-map)

## Completion check

After the core, you should be able to reconstruct this dependency chain:

1. consistency locates the probability limit;

2. a rate and expansion identify the leading random term;

3. a CLT or other limit theorem gives its limiting distribution;

4. the delta method handles smooth transformations;

5. variance consistency makes studentization feasible; and

6. a separate uniform argument is needed for worst case validity; and

7. discretization, finite dimensional convergence, and oscillation control are the three ingredients in the chapter's uniform LLN proof.
