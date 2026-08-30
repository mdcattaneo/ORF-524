# ORF 524: Statistical Theory and Methods

## Week 5: Convergence and limit theorems

**Status:** Ready for Fall 2026  
**Last updated:** August 30, 2026  
**Meetings:** Monday, September 28, and Wednesday, September 30, 2026

<!-- Legacy sources: 11 ORF 524/lectures/ORF524-Fall2025-Lecture-Week-05.tex and
ORF524-Fall2025-Lecture-Week-05_Lecture.pdf; 09 ORF 524/lectures/
ORF524-Fall2024-Lecture-Week-05_Notes.pdf. The original asymptotics chapter is now taught across
Weeks 5 and 6. -->

## Central question

Which convergence concepts and limit theorems justify replacing random sample quantities by stable population approximations?

## Learning goals

By the end of the week, you should be able to:

1. distinguish pointwise, almost sure, in probability, in distribution, and $L^p$ convergence;
2. apply continuous mapping and Slutsky's theorem with their actual assumptions;
3. manipulate $O_{\mathbb{P}}$ and $o_{\mathbb{P}}$ without confusing magnitude with distribution;
4. state and apply i.i.d., independent but not identically distributed, and triangular array laws of large numbers and central limit theorems;
5. use Cramer--Wold to derive the multivariate i.i.d. CLT; and
6. recognize when a rate differs from $\sqrt n$ or a limiting distribution is non-Gaussian.

## In-class route

This is one continuous route through the week, not a division into class meetings. Proceed in order, stop wherever the discussion ends, and resume at the next unfinished stop. The linked sections contain the full student facing exposition; the route identifies the concepts and activities to foreground live.

| Stop | Live focus | Mode |
|---|---|---|
| **W5.1** | [Modes of convergence and false converses](#w5-stop-1) | Discuss + counterexample + AI implication audit + Checkpoint 1 |
| **W5.2** | [Continuous mapping and Slutsky](#w5-stop-2) | Board work + Checkpoint 2 |
| **W5.3** | [Stochastic orders](#w5-stop-3) | Discuss + algebra audit + Checkpoint 3 |
| **W5.4** | [LLNs and the classical CLT](#w5-stop-4) | Board work + assumption comparison |
| **W5.5** | [Multivariate, triangular array, and non-Gaussian limits](#w5-stop-5) | Board work + boundary comparison + Checkpoint 4 |

Within Stop W5.4, the Chebyshev LLN, sample variance consistency, and exact versus asymptotic Studentization are the main live calculations. Berry--Esseen and the Lyapunov to Lindeberg bridge may be used selectively before Stop W5.5's multivariate and boundary examples.

## How to use this chapter

**Prepare:** Read Sections 1 and 2, classify every convergence statement by mode and limit, and attempt Checkpoint 1.

**In class:** Follow the continuous route above through convergence modes, mapping results, stochastic orders, LLNs, CLTs, multivariate limits, and a non-Gaussian boundary example. Prepared responses supply an alternative for every live prompt if AI is unavailable.

**Review:** For each limit theorem, identify the centering, normalization, dependence and moment conditions, and exact conclusion. Week 6 will turn these probability limits into statistical inference.

**Practice:** The matching assignment complements this chapter with ungraded problems for consolidating the material, extending selected ideas, and preparing for cumulative review.

**Prerequisites:** Weeks 1--4; probability inequalities; expectations, variance, and random vectors.

## Full chapter map

1. [Modes of convergence](#1-modes-of-convergence)
2. [Continuous mapping and Slutsky](#2-continuous-mapping-and-slutsky)
3. [Stochastic orders](#3-stochastic-orders)
4. [Laws of large numbers and central limit theorems](#4-laws-of-large-numbers-and-central-limit-theorems)
5. [Recap](#5-recap)

<a id="w5-stop-1"></a>

## 1. Modes of convergence

Let $X_n$ and $X$ be random elements defined on a common probability space when the definition requires it.

### 1.1 Convergence in probability

We say $X_n$ **converges in probability** to $X$, written $X_n\to_{\mathbb{P}}X$, if for every $\varepsilon>0$,

$$
\mathbb{P}(\Vert X_n-X\Vert>\varepsilon)\to0.
$$

For a constant limit $c$, this says the random variable concentrates in every fixed neighborhood of $c$.

### 1.2 Convergence in distribution

For real valued random variables, $X_n$ **converges in distribution** to $X$, written $X_n\rightsquigarrow X$, if

$$
F_{X_n}(x)\to F_X(x)
$$

at every continuity point of $F_X$. Equivalently, in a metric space,

$$
\mathbb E[f(X_n)]\to\mathbb E[f(X)]
$$

for every bounded continuous $f$.

On a separable metric space, it is enough to verify this expectation criterion for every bounded Lipschitz $f$. For random vectors in $\mathbb R^d$, Lévy's continuity theorem provides another route: pointwise convergence of characteristic functions to a characteristic function continuous at the origin implies convergence in distribution.

Convergence in distribution concerns laws and does not require $X_n$ and $X$ to be defined on the same probability space.

### 1.3 Almost sure and $L^p$ convergence

We have **almost sure convergence**, $X_n\overset{a.s.}{\to}X$, if

$$
\mathbb{P}\left(\lim_{n\to\infty}X_n=X\right)=1.
$$

For $p>0$, we have **$L^p$ convergence**, $X_n\overset{L^p}{\to}X$, if

$$
\mathbb E[\Vert X_n-X\Vert^p]\to0.
$$

If all variables are defined as functions on a common sample space, **pointwise convergence** means

$$
X_n(\omega)\to X(\omega)
\qquad\text{for every }\omega\in\Omega.
$$

Pointwise convergence is stronger than almost sure convergence, which permits failure on a null set. Because random variables are identified up to almost sure equality in most statistical arguments, almost sure convergence is usually the more natural notion.

The main implications are

$$
X_n(\omega)\to X(\omega)\ \text{for every }\omega
\quad\Longrightarrow\quad
X_n\overset{a.s.}{\to}X
\quad\Longrightarrow\quad
X_n\to_{\mathbb{P}}X
\quad\Longrightarrow\quad
X_n\rightsquigarrow X,
$$

and

$$
X_n\overset{L^p}{\to}X
\quad\Longrightarrow\quad
X_n\to_{\mathbb{P}}X.
$$

The reverse implications generally fail. An important exception is

$$
X_n\rightsquigarrow c
\quad\Longrightarrow\quad
X_n\to_{\mathbb{P}}c
$$

when the limit is constant.

### Counterexample: identical laws do not imply convergence in probability

Let $X$ be nondegenerate and symmetric about zero, and define $X_n=(-1)^nX$. Every $X_n$ has the same distribution as $X$, so $X_n\rightsquigarrow X$ in the distributional sense. But $X_n-X$ does not converge to zero in probability. Convergence in distribution does not track the joint coupling between $X_n$ and $X$.

### Checkpoint 1

For each claim, decide whether it is true and supply a theorem or counterexample.

1. Convergence in distribution implies convergence in probability.
2. Convergence in probability implies convergence in distribution.
3. Convergence in probability implies convergence of expectations.
4. $L^2$ convergence implies $L^1$ convergence on a probability space.
5. Almost sure convergence implies convergence in probability.
6. Pointwise convergence and almost sure convergence are equivalent for random variables regarded only up to almost sure equality.

> [!TIP]
> **AI interaction 1 — Build the implication diagram**
>
> Copy the prompt below into an AI interface and audit the response.

```text
Construct an implication diagram for pointwise, almost sure, in probability,
in distribution, and Lp convergence. For every implication you assert, state
its conditions. For every tempting reverse implication that fails, give a
fully specified counterexample and verify the failure. Treat a constant weak
limit as a special case.
```

**Audit question:** Does the response mistake equality of marginal distributions for convergence in probability on a common probability space?

<a id="w5-stop-2"></a>

## 2. Continuous mapping and Slutsky

### 2.1 Continuous mapping theorem

If $X_n\to_{\mathbb{P}}X$ and $g$ is continuous at $X$ almost surely, then

$$
g(X_n)\to_{\mathbb{P}}g(X).
$$

If $X_n\rightsquigarrow X$ and $g$ is continuous at $X$ almost surely, then

$$
g(X_n)\rightsquigarrow g(X).
$$

Continuity is needed only on a set containing the limit with probability one. A discontinuity that the limiting variable hits with positive probability can invalidate the conclusion.

### 2.2 Slutsky's theorem

If

$$
X_n\rightsquigarrow X,
\qquad
Y_n\to_{\mathbb{P}}c,
$$

where $c$ is constant, then jointly

$$
(X_n,Y_n)\rightsquigarrow(X,c).
$$

Consequently,

$$
X_n+Y_n\rightsquigarrow X+c,
\qquad
X_nY_n\rightsquigarrow cX,
$$

and, when $c\neq0$,

$$
\frac{X_n}{Y_n}\rightsquigarrow\frac{X}{c}.
$$

Slutsky permits a random quantity converging to a constant to be replaced asymptotically. It does not permit replacing a random quantity by a nondegenerate random limit without establishing joint convergence.

> [!IMPORTANT]
> **Board work 1 — Slutsky from joint convergence**
>
> Starting only from a joint version of continuous mapping, derive the three Slutsky conclusions. Then construct an example showing why marginal convergence of $X_n$ and $Y_n$ is insufficient to determine the distribution of $X_n+Y_n$ when both limits are nondegenerate.

### Checkpoint 2

Suppose $\sqrt n(T_n-\theta)\rightsquigarrow\mathcal N(0,V)$ and $\widehat V_n\to_{\mathbb{P}}V>0$. Derive the limit of

$$
\frac{\sqrt n(T_n-\theta)}{\sqrt{\widehat V_n}}.
$$

Name every theorem used.

<a id="w5-stop-3"></a>

## 3. Stochastic orders

A sequence $X_n$ is **bounded in probability**, written $X_n=O_{\mathbb{P}}(1)$, if for every $\varepsilon>0$ there are $M<\infty$ and $N$ such that

$$
\mathbb{P}(\Vert X_n\Vert>M)<\varepsilon
\quad\text{for all }n\geq N.
$$

For positive constants $a_n$,

$$
X_n=O_{\mathbb{P}}(a_n)
\quad\Longleftrightarrow\quad
X_n/a_n=O_{\mathbb{P}}(1),
$$

and

$$
X_n=o_{\mathbb{P}}(a_n)
\quad\Longleftrightarrow\quad
X_n/a_n\to_{\mathbb{P}}0.
$$

More generally, if $Y_n>0$ with probability approaching one, then

$$
X_n=O_{\mathbb{P}}(Y_n)
\quad\Longleftrightarrow\quad
X_n/Y_n=O_{\mathbb{P}}(1),
$$

and

$$
X_n=o_{\mathbb{P}}(Y_n)
\quad\Longleftrightarrow\quad
X_n/Y_n\to_{\mathbb{P}}0.
$$

Random normalizations require care when $Y_n$ can vanish or change sign. Deterministic positive sequences are the default unless stated otherwise.

Useful algebra includes

$$
O_{\mathbb{P}}(a_n)O_{\mathbb{P}}(b_n)=O_{\mathbb{P}}(a_nb_n),
$$

$$
o_{\mathbb{P}}(a_n)O_{\mathbb{P}}(b_n)=o_{\mathbb{P}}(a_nb_n),
$$

$$
O_{\mathbb{P}}(a_n)+O_{\mathbb{P}}(b_n)=O_{\mathbb{P}}(a_n+b_n).
$$

If $X_n\rightsquigarrow X$, then $X_n=O_{\mathbb{P}}(1)$. If $\sqrt n(T_n-\theta)=O_{\mathbb{P}}(1)$, then

$$
T_n-\theta=O_{\mathbb{P}}(n^{-1/2}).
$$

Stochastic order describes magnitude, not a limiting distribution.

### Checkpoint 3

Let $X_n=O_{\mathbb{P}}(1)$, $Y_n=o_{\mathbb{P}}(1)$, and $a_n\to0$.

1. Classify $X_nY_n$.
2. Classify $a_nX_n$.
3. Does $X_n=O_{\mathbb{P}}(1)$ imply $\mathbb E[X_n^2]=O(1)$?
4. Does $T_n-\theta=O_{\mathbb{P}}(n^{-1/2})$ imply asymptotic normality?

<a id="w5-stop-4"></a>

## 4. Laws of large numbers and central limit theorems

### 4.1 Laws of large numbers

If $X_1,X_2,\ldots$ are i.i.d. and $\mathbb E[|X_1|]<\infty$, the weak law of large numbers gives

$$
\overline X_n\to_{\mathbb{P}}\mathbb E[X_1].
$$

With finite variance, Chebyshev provides the short proof

$$
\mathbb{P}(|\overline X_n-\mu|>\varepsilon)
\leq
\frac{\sigma^2}{n\varepsilon^2}
\to0.
$$

The finite first moment result is stronger than this proof and requires a different argument. A law of large numbers supplies concentration around a population quantity; it does not supply the shape of the approximation error.

### 4.2 Classical central limit theorem

If $X_1,X_2,\ldots$ are i.i.d. with mean $\mu$ and variance $0<\sigma^2<\infty$, then

$$
\frac{\sqrt n(\overline X_n-\mu)}{\sigma}
\rightsquigarrow
\mathcal N(0,1).
$$

The LLN and CLT answer different questions:

- the LLN gives $\overline X_n-\mu=o_{\mathbb{P}}(1)$;
- the CLT gives the scale $n^{-1/2}$ and a limiting distribution.

If $S_n^2\to_{\mathbb{P}}\sigma^2$, Slutsky gives the self normalized result

$$
\frac{\sqrt n(\overline X_n-\mu)}{S_n}
\rightsquigarrow
\mathcal N(0,1).
$$

This conclusion is asymptotic outside the Gaussian model; it should not be confused with the exact Student $t$ distribution under Gaussian sampling.

For the usual sample variance, consistency follows from an explicit algebraic reduction:

$$
S_n^2
=\frac{1}{n-1}\sum_{i=1}^n(X_i-\overline X_n)^2
=\frac{n}{n-1}
\left\lbrace
\frac1n\sum_{i=1}^n(X_i-\mu)^2-(\overline X_n-\mu)^2
\right\rbrace
\to_{\mathbb{P}}\sigma^2.
$$

The LLN handles the explicit sample average $n^{-1}\sum_{i=1}^n(X_i-\mu)^2$, while the LLN and continuous mapping make the squared sample mean error vanish. This calculation is the missing bridge between the numerator CLT and Studentization.

> [!IMPORTANT]
> **Board work 2 — Studentization dependency map**
>
> Write the self normalized CLT as an explicit dependency graph:
>
> 1. the CLT for the numerator;
> 2. consistency of $S_n^2$;
> 3. continuous mapping for $S_n$ and $1/S_n$; and
> 4. Slutsky for the ratio.
>
> Identify which assumptions are used at each step.
>
> Then derive the displayed identity for $S_n^2$ and state exactly where a finite second moment is used.

<a id="w5-stop-5"></a>

### 4.3 Cramer--Wold and the multivariate CLT

The **Cramer--Wold device** says that random vectors $Y_n,Y\in\mathbb R^d$ satisfy $Y_n\rightsquigarrow Y$ if and only if

$$
a'Y_n\rightsquigarrow a'Y
\qquad\text{for every }a\in\mathbb R^d.
$$

Let $X_1,X_2,\ldots$ now be i.i.d. random vectors with $\mathbb E[\Vert X_1\Vert^2]<\infty$, mean $\mu$, and covariance matrix $\Sigma$. Then

$$
\sqrt n(\overline X_n-\mu)
\rightsquigarrow
\mathcal N_d(0,\Sigma).
$$

**Proof map.** For every fixed $a$, the scalar variables $a'X_i$ are i.i.d. with mean $a'\mu$ and variance $a'\Sigma a$. The scalar CLT gives

$$
a'\sqrt n(\overline X_n-\mu)
\rightsquigarrow
\mathcal N(0,a'\Sigma a).
$$

These are exactly the projections of $\mathcal N_d(0,\Sigma)$, so Cramer--Wold completes the argument. The covariance matrix may be singular; directions with $a'\Sigma a=0$ simply have a degenerate limiting projection.

### 4.4 Beyond i.i.d. laws of large numbers

Independence and identical distributions are convenient sufficient conditions, not the definition of an LLN. If $X_1,X_2,\ldots$ are independent with means $\mu_i$ and finite variances, and

$$
\frac1n\sum_{i=1}^n\mu_i\to\mu,
\qquad
\frac1{n^2}\sum_{i=1}^n\mathbb V[X_i]\to0,
$$

then Chebyshev's inequality gives

$$
\frac1n\sum_{i=1}^nX_i\to_{\mathbb{P}}\mu.
$$

Independence can also be weakened. For example, if the variables have common mean $\mu$, common finite variance, and

$$
\mathrm{Cov}(X_i,X_j)=\rho(|i-j|),
\qquad \rho(h)\to0,
$$

then $\mathbb V(\overline X_n)\to0$ and hence $\overline X_n\to_{\mathbb{P}}\mu$. The proof is a variance calculation followed by Chebyshev.

As a simpler special case, pairwise uncorrelated variables with common mean $\mu$ and common finite variance $\sigma^2$ satisfy

$$
\mathbb V(\overline X_n)=\frac{\sigma^2}{n},
$$

so the same Chebyshev proof works without mutual independence.

A useful row wise version makes the normalization explicit. Let $X_{n1},\ldots,X_{nk_n}$ be independent within each row, with finite variances. If $v_n>0$ and

$$
\frac1{v_n^2}\sum_{i=1}^{k_n}\mathbb V[X_{ni}]\to0,
$$

then the triangular array LLN gives

$$
\frac1{v_n}\sum_{i=1}^{k_n}
\lbrace X_{ni}-\mathbb E[X_{ni}]\rbrace
\to_{\mathbb{P}}0.
$$

Truncation arguments yield versions that do not assume finite variances. The conditions of the LLN must match the dependence, tail behavior, and normalization in the application.

### 4.5 Reference: triangular array central limit theorem

For independent mean zero variables $X_{n1},\ldots,X_{nk_n}$, let

$$
s_n^2=\sum_{i=1}^{k_n}\mathbb V[X_{ni}].
$$

If $s_n^2>0$ and, for every $\varepsilon>0$,

$$
\frac1{s_n^2}
\sum_{i=1}^{k_n}
\mathbb E\left[
X_{ni}^2\mathbf 1\lbrace|X_{ni}|>\varepsilon s_n\rbrace
\right]
\to0,
$$

then the Lindeberg--Feller theorem gives

$$
\frac{\sum_iX_{ni}}{s_n}
\rightsquigarrow\mathcal N(0,1).
$$

The Lindeberg condition prevents a single large term from dominating the normalized sum. This result is reference material for later arguments involving observations that are not identically distributed and for local asymptotic arguments.

A convenient sufficient condition is **Lyapunov's condition**: for some $\delta>0$,

$$
\frac1{s_n^{2+\delta}}
\sum_{i=1}^{k_n}\mathbb E[|X_{ni}|^{2+\delta}]
\to0.
$$

Indeed, on $|X_{ni}|>\varepsilon s_n$,

$$
X_{ni}^2
\leq
\frac{|X_{ni}|^{2+\delta}}{(\varepsilon s_n)^\delta},
$$

so Lyapunov implies Lindeberg after summing and dividing by $s_n^2$.

For i.i.d. summands, the **Berry--Esseen theorem** makes the classical CLT quantitative. If $\rho=\mathbb E[|X_1-\mu|^3]<\infty$, then for a universal constant $C$,

$$
\sup_{t\in\mathbb R}
\left|
\mathbb{P}\mkern-3mu\left(
\frac{\sqrt n(\overline X_n-\mu)}{\sigma}\leq t
\right)-\Phi(t)
\right|
\leq
\frac{C\rho}{\sigma^3\sqrt n}.
$$

This is a uniform CDF error bound, not a relative tail approximation; it also requires a third absolute moment, which the ordinary finite variance CLT does not.

### 4.6 A non-Gaussian order statistic limit

Asymptotic distributions need not be Normal or use the $\sqrt n$ rate. If $X_1,\ldots,X_n$ are i.i.d. $\mathsf{Uniform}(0,\theta)$ and $M_n=\max_iX_i$, then for every fixed $t\geq0$ and all $n>t$,

$$
\mathbb{P}_\theta\mkern-3mu\left\lbrace
\frac{n(\theta-M_n)}{\theta}>t
\right\rbrace
=\mathbb{P}_\theta\mkern-3mu\left\lbrace M_n<\theta\left(1-\frac tn\right)\right\rbrace
=\left(1-\frac tn\right)^n
\longrightarrow e^{-t}.
$$

Therefore

$$
\frac{n(\theta-M_n)}{\theta}
\rightsquigarrow\mathsf{Exponential}(1).
$$

The maximum converges at rate $n$, its error is one sided, and its limit is non-Gaussian. This is the asymptotic signature of the parameter dependent boundary already seen in Weeks 2--4.

> [!IMPORTANT]
> **Board work 3 — Multivariate Gaussian and boundary limits**
>
> 1. Derive the multivariate i.i.d. CLT by applying the scalar CLT to every fixed projection and then invoking Cramer--Wold.
> 2. Explain how the argument accommodates a singular covariance matrix.
> 3. Contrast the centering, rate, support, and limiting law with the i.i.d. Uniform maximum.

### Checkpoint 4

1. Why does a finite mean suffice for the i.i.d. weak LLN but not for the classical finite variance CLT stated here?
2. Give a setting in which an exact Normal result is available and one in which Normality is only an approximation.
3. What feature of a triangular array is absent from a single i.i.d. sequence?
4. Which variance condition proves the LLN for independent but not identically distributed observations above?
5. Why does an LLN for a triangular array not by itself imply a triangular array CLT?
6. How does Cramer--Wold reduce a multivariate CLT to scalar CLTs?
7. Which features of the Uniform maximum show that a Normal approximation scaled by $\sqrt n$ is not universal?
8. How does Lyapunov's condition imply Lindeberg's condition, and which one is closer to the no dominant term idea?
9. What does Berry--Esseen add to the CLT, and what extra moment does it require?

## 5. Recap

- Pointwise, almost sure, in probability, in distribution, and $L^p$ convergence answer different questions and are not interchangeable.
- Continuous mapping transforms convergent sequences; Slutsky permits consistent random quantities to replace constants.
- $O_{\mathbb{P}}$ and $o_{\mathbb{P}}$ encode stochastic magnitude, not limiting shape.
- LLNs justify concentration; CLTs justify normalized distributional approximations.
- Theorems for observations that are not identically distributed and for triangular arrays require conditions matched to dependence, tails, and normalization.
- Cramer--Wold converts scalar projection limits into a multivariate limit.
- Boundary estimators can have faster rates and non-Gaussian limits, as the Uniform maximum shows.
- A theorem's assumptions are part of the conclusion; “by the LLN” or “by the CLT” is not a complete argument.

## Notation introduced this week

- $\to_{\mathbb{P}}$, $\rightsquigarrow$, $\overset{a.s.}{\to}$, and $\overset{L^p}{\to}$: convergence in probability, distribution, almost surely, and in $L^p$.
- $O_{\mathbb{P}}(a_n)$ and $o_{\mathbb{P}}(a_n)$: stochastic order and negligible stochastic order.
- $\Sigma$: covariance matrix in a multivariate CLT.
- $s_n$: row specific standard deviation in a triangular array CLT.

## References

- Jun Shao, *Mathematical Statistics*, Springer, 2003.
- Robert W. Keener, *Theoretical Statistics: Topics for a Core Course*, Springer, 2010.
- Aad van der Vaart, *Asymptotic Statistics*, Cambridge University Press, 1998.
- Martin J. Wainwright, *High-Dimensional Statistics: A Non-Asymptotic Viewpoint*, Cambridge University Press, 2019.
