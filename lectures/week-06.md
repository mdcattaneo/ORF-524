# ORF 524: Statistical Theory and Methods

## Week 6: Delta method and asymptotic inference

**Status:** Ready for Fall 2026  
**Last updated:** August 30, 2026  
**Meetings:** Monday, October 5, and Wednesday, October 7, 2026

<!-- Legacy sources: 11 ORF 524/lectures/ORF524-Fall2025-Lecture-Week-05.tex and
ORF524-Fall2025-Lecture-Week-05_Lecture.pdf; 09 ORF 524/lectures/
ORF524-Fall2024-Lecture-Week-05_Notes.pdf. This chapter is the inference half of the original
asymptotic-approximations week. -->

## Central question

How do a limiting distribution, a smooth transformation, and a consistent standard error combine to justify large sample statistical inference?

## Learning goals

By the end of the week, you should be able to:

1. derive scalar and multivariate transformed limits with the delta method;
2. diagnose a vanishing first derivative and obtain an appropriate second order limit;
3. distinguish consistency, convergence rate, asymptotic distribution, and moment convergence;
4. organize and execute a complete approximation chain to construct studentized asymptotic tests and confidence intervals;
5. distinguish exact, pointwise asymptotic, and uniform asymptotic validity; and
6. prove a usable uniform law of large numbers and explain how it supports consistency arguments for M- and Z-estimators.

## In-class route

This is one continuous route through the week, not a division into class meetings. Proceed in order, stop wherever the discussion ends, and resume at the next unfinished stop. The linked sections contain the full student facing exposition; the route identifies the concepts and activities to foreground live.

| Stop | Live focus | Mode |
|---|---|---|
| **W6.1** | [First order delta method at a general rate](#w6-stop-1) | Board work |
| **W6.2** | [Vanishing derivative and second order limits](#w6-stop-2) | AI diagnosis + Checkpoint 1 |
| **W6.3** | [Consistency, rates, bias, variance, and efficiency](#w6-stop-3) | Discuss + claim classification |
| **W6.4** | [Studentization, tests, and confidence intervals](#w6-stop-4) | Board work |
| **W6.5** | [Size, power, and pointwise versus uniform validity](#w6-stop-5) | AI interval audit + Checkpoint 2 |
| **W6.6** | [Bridge to general estimation](#w6-stop-6) | Discuss + transfer Checkpoint 3 |
| **W6.7** | [A uniform law of large numbers](#w6-stop-7) | Discuss proof architecture |

The general rate Taylor expansion and stochastic order bookkeeping are the board spine of W6.1. Stops W6.4--W6.5 should preserve the entire dependency chain from a limiting distribution through a consistent standard error to the precisely qualified inferential claim. W6.7 proves a uniform LLN from compactness, continuity, and an integrable envelope; this theorem becomes a working tool in Weeks 9 and 10.

## How to use this chapter

**Prepare:** Review Week 5's CLTs, Slutsky theorem, and stochastic orders. Read Section 1 and attempt Checkpoint 1.

**In class:** Follow the continuous route above through first and second order delta methods, rates, studentization, asymptotic testing and confidence intervals, pointwise versus uniform validity, and the uniform LLN. Prepared responses supply an alternative for every live prompt if AI is unavailable.

**Review:** Reconstruct each inference result as a dependency chain. “By asymptotic normality” is incomplete unless the centering, scaling, variance, standard error consistency, and validity claim have all been identified.

**Practice:** The matching assignment complements this chapter with ungraded problems for consolidating the material, extending selected ideas, and preparing for cumulative review.

**Prerequisites:** Weeks 1--5; Taylor expansion; continuous mapping; Slutsky; stochastic orders; scalar and multivariate CLTs; exact finite sample inference in a regular model.

## Full chapter map

1. [Delta method](#1-delta-method)
2. [Consistency, rates, and asymptotic inference](#2-consistency-rates-and-asymptotic-inference)
3. [Bridge to general estimation](#3-bridge-to-general-estimation)
4. [A uniform law of large numbers](#4-a-uniform-law-of-large-numbers)
5. [Recap](#5-recap)

<a id="w6-stop-1"></a>

## 1. Delta method

Suppose

$$
r_n(T_n-\theta)\rightsquigarrow Z,
\qquad r_n\to\infty,
$$

and $g$ is differentiable at $\theta$. Then

$$
r_n\lbrace g(T_n)-g(\theta)\rbrace
\rightsquigarrow
g'(\theta)Z.
$$

In the multivariate case, if

$$
\sqrt n(T_n-\theta)\rightsquigarrow\mathcal N(0,V)
$$

and $g:\mathbb R^k\to\mathbb R^m$ is differentiable with Jacobian $G(\theta)$, then

$$
\sqrt n\lbrace g(T_n)-g(\theta)\rbrace
\rightsquigarrow
\mathcal N\mkern-3mu\left(0,G(\theta)VG(\theta)'\right).
$$

### Proof

Differentiability gives

$$
g(T_n)-g(\theta)
=G(\theta)(T_n-\theta)+r(T_n-\theta),
$$

where

$$
\frac{\Vert r(h)\Vert}{\Vert h\Vert}\to0
\quad\text{as }h\to0.
$$

The assumed limiting distribution implies $T_n\to_{\mathbb{P}}\theta$ and $r_n(T_n-\theta)=O_{\mathbb{P}}(1)$. Therefore

$$
\Vert r_n r(T_n-\theta)\Vert
=\frac{\Vert r(T_n-\theta)\Vert}{\Vert T_n-\theta\Vert}
\Vert r_n(T_n-\theta)\Vert
=o_{\mathbb{P}}(1).
$$

The ratio is defined as zero on the event $T_n=\theta$, consistent with $r(0)=0$.

Slutsky then yields the result.

> [!IMPORTANT]
> **Board work 1 — The delta method at a general rate**
>
> Start from $r_n(T_n-\theta)\rightsquigarrow Z$ rather than assuming $r_n=\sqrt n$. Write the Taylor expansion with an explicit remainder, show separately that $T_n\to_{\mathbb{P}}\theta$ and $r_n(T_n-\theta)=O_{\mathbb{P}}(1)$, and then prove that the scaled remainder is $o_{\mathbb{P}}(1)$. Mark the precise step at which differentiability, rather than mere continuity, is needed.

<a id="w6-stop-2"></a>

### Failure mode: zero derivative

If $g'(\theta)=0$, the first order delta method gives a degenerate limit. It does not follow that the transformed estimator has no sampling variation. A second order expansion may yield a different rate and a non-Gaussian limit.

### Second order delta method

Suppose $r_n(T_n-\theta)\rightsquigarrow Z$, where $r_n\to\infty$, and $g$ is twice differentiable at $\theta$ with $g'(\theta)=0$. Then

$$
r_n^2\lbrace g(T_n)-g(\theta)\rbrace
\rightsquigarrow
\frac12 g''(\theta)Z^2.
$$

Indeed, the second order expansion gives

$$
g(T_n)-g(\theta)
=\frac12g''(\theta)(T_n-\theta)^2
+o_{\mathbb{P}}\mkern-3mu\left((T_n-\theta)^2\right).
$$

Multiplication by $r_n^2$ yields

$$
r_n^2\lbrace g(T_n)-g(\theta)\rbrace
=\frac12g''(\theta)\lbrace r_n(T_n-\theta)\rbrace^2+o_{\mathbb{P}}(1),
$$

because $r_n(T_n-\theta)=O_{\mathbb{P}}(1)$. The continuous mapping theorem and Slutsky give the result. If $Z\sim\mathcal N(0,V)$, then $Z^2\sim V\chi_1^2$, so the limit is generally neither centered nor Normal.

For example, if $\sqrt n T_n\rightsquigarrow Z$ and $g(t)=t^2$ at $\theta=0$, then

$$
nT_n^2\rightsquigarrow Z^2,
$$

not a nondegenerate Normal limit scaled by $\sqrt n$.

> [!TIP]
> **AI interaction 1 — The derivative vanished**
>
> Copy the prompt below into an AI interface and audit the response.

```text
Suppose sqrt(n)(T_n-theta) converges in distribution to N(0,V), and g is
twice continuously differentiable near theta with g'(theta)=0 and
g''(theta) not equal to zero.

Audit the claim that sqrt(n)(g(T_n)-g(theta)) is asymptotically Normal with
variance zero. Explain what that statement does and does not say, derive the
appropriate second order scaling, and identify the resulting limit. State any
additional remainder conditions you use.
```

### Checkpoint 1

Let $\sqrt n(\widehat p_n-p)\rightsquigarrow\mathcal N(0,p(1-p))$, with the logit defined on $\lbrace\widehat p_n\in(0,1)\rbrace$ and extended arbitrarily on its asymptotically negligible complement.

1. Derive the limit for $\log\lbrace\widehat p_n/(1-\widehat p_n)\rbrace$ when $p\in(0,1)$.
2. What fails as $p$ approaches zero or one?
3. Is the resulting approximation pointwise or uniform over $p\in(0,1)$?

<a id="w6-stop-3"></a>

## 2. Consistency, rates, and asymptotic inference

### 2.1 Consistency and rate

An estimator $\widehat\theta_n$ is **consistent** for $\theta_0$ if

$$
\widehat\theta_n\to_{\mathbb{P}}\theta_0.
$$

It is consistent at rate $r_n$ if

$$
r_n(\widehat\theta_n-\theta_0)=O_{\mathbb{P}}(1).
$$

Two estimators $T_n$ and $S_n$ are **first order asymptotically equivalent at rate $r_n$** if

$$
r_n(T_n-S_n)\to_{\mathbb{P}}0.
$$

If $r_n(T_n-\theta_0)\rightsquigarrow Z$, then Slutsky implies $r_n(S_n-\theta_0)\rightsquigarrow Z$. Thus an $o_{\mathbb{P}}(r_n^{-1})$ difference does not alter the first order limiting distribution.

Consistency identifies the probability limit. A rate describes the magnitude of the error. An asymptotic distribution describes the limiting shape after centering and scaling. None of these statements alone supplies the others.

Convergence in distribution does not generally imply convergence of moments. Consequently, asymptotic bias or variance cannot be read from a limiting distribution without conditions such as uniform integrability.

### 2.2 Leading asymptotic bias, variance, and efficiency

An integrable estimator is **asymptotically unbiased** for $\theta_0$ if

$$
\mathbb E[\widehat\theta_n]-\theta_0\to0.
$$

This is a statement about moments and is neither implied by consistency nor sufficient for it without additional conditions.

Suppose

$$
r_n(\widehat\theta_n-\theta_0)\rightsquigarrow Z
$$

and the relevant moments converge as well—for example, under suitable uniform integrability conditions. For a scalar estimator, the first order asymptotic bias, variance, and MSE on the $r_n$ scale are

$$
b=\mathbb E[Z],
\qquad
V=\mathbb V[Z],
\qquad
M=\mathbb E[Z^2]=V+b^2.
$$

Equivalently, when moment convergence holds,

$$
r_n\lbrace\mathbb E[\widehat\theta_n]-\theta_0\rbrace\to b,
$$

and

$$
r_n^2\mathrm{MSE}(\widehat\theta_n)\to M.
$$

For two estimators of the same scalar target with the same rate and centered limiting distributions, we define the asymptotic relative efficiency of estimator 1 to estimator 2 here as

$$
\mathrm{ARE}(1,2)=\frac{V_2}{V_1}.
$$

Under this convention, values above one favor estimator 1. The convention must be stated, because some sources reverse the ratio.

If the limiting distributions are not centered, the corresponding asymptotic MSE comparison replaces $V_j$ by $M_j=\mathbb E[Z_j^2]$. Estimators must be compared at the same rate and for the same target.

Weak convergence alone does not justify any of these moment conclusions.

### Example: bias that disappears at first order

Let

$$
\widehat\sigma_n^2=\frac1n\sum_{i=1}^n(X_i-\overline X_n)^2,
\qquad
S_n^2=\frac1{n-1}\sum_{i=1}^n(X_i-\overline X_n)^2.
$$

When $X_1,\ldots,X_n$ are i.i.d. with variance $\sigma^2<\infty$,

$$
\mathbb E[\widehat\sigma_n^2]=\frac{n-1}{n}\sigma^2,
\qquad
\mathbb E[S_n^2]=\sigma^2.
$$

Thus $\widehat\sigma_n^2$ has finite sample bias $-\sigma^2/n$, whereas $S_n^2$ is unbiased. Nevertheless,

$$
S_n^2-\widehat\sigma_n^2
=\frac{\widehat\sigma_n^2}{n-1}
=O_{\mathbb{P}}(n^{-1}),
$$

provided $\widehat\sigma_n^2=O_{\mathbb{P}}(1)$. Consequently, whenever either estimator has a limiting distribution under $\sqrt n$ scaling, the other has the same limit. In particular, if $\mathbb E[|X_1|^4]<\infty$ and the limiting variance is positive, they have the same first order asymptotic variance and ARE equal to one. First order asymptotic comparisons can therefore miss meaningful finite sample differences.

<a id="w6-stop-4"></a>

### 2.3 Studentization

Suppose

$$
\sqrt n(\widehat\theta_n-\theta_0)
\rightsquigarrow
\mathcal N(0,V),
$$

and $\widehat V_n\to_{\mathbb{P}}V>0$. Then

$$
\frac{\sqrt n(\widehat\theta_n-\theta_0)}{\sqrt{\widehat V_n}}
\rightsquigarrow
\mathcal N(0,1).
$$

An asymptotic two sided level $\alpha$ test of $H_0:\theta=\theta_0$ rejects when

$$
\left|
\frac{\sqrt n(\widehat\theta_n-\theta_0)}{\sqrt{\widehat V_n}}
\right|>z_{1-\alpha/2}.
$$

The corresponding confidence interval is

$$
\left[
\widehat\theta_n-z_{1-\alpha/2}\sqrt{\frac{\widehat V_n}{n}},
\widehat\theta_n+z_{1-\alpha/2}\sqrt{\frac{\widehat V_n}{n}}
\right].
$$

These are pointwise asymptotic statements unless a uniform theorem has been established.

### Worked example: a ratio of means

Let $(X_i,Y_i)'$, $i=1,\ldots,n$, be i.i.d. with mean $\mu=(\mu_X,\mu_Y)'$, covariance matrix $\Sigma$, finite second moments, and $\mu_Y\neq0$. The multivariate CLT gives

$$
\sqrt n
\left\lbrace
\begin{pmatrix}\overline X_n\\ \overline Y_n\end{pmatrix} -
\begin{pmatrix}\mu_X\\ \mu_Y\end{pmatrix}
\right\rbrace
\rightsquigarrow
\mathcal N(0,\Sigma).
$$

For $g(a,b)=a/b$, define $\theta=\mu_X/\mu_Y$ and

$$
Dg(\mu) =
\begin{pmatrix}
1/\mu_Y & -\mu_X/\mu_Y^2
\end{pmatrix}.
$$

The multivariate delta method yields

$$
\sqrt n
\left(
\frac{\overline X_n}{\overline Y_n}-\theta
\right)
\rightsquigarrow
\mathcal N(0,V),
\qquad
V=Dg(\mu)\Sigma Dg(\mu)'.
$$

Let $\widehat\Sigma_n$ be the sample covariance matrix and let

$$
\widehat D_n =
\begin{pmatrix}
1/\overline Y_n & -\overline X_n/\overline Y_n^2
\end{pmatrix},
\qquad
\widehat V_n=\widehat D_n\widehat\Sigma_n\widehat D_n'.
$$

On the event $\overline Y_n\neq0$, which has probability approaching one, the law of large numbers and continuous mapping give $\widehat V_n\to_{\mathbb{P}}V$. Define the estimators arbitrarily on the asymptotically negligible complementary event. If $V>0$, Studentization therefore gives

$$
\frac{\sqrt n\left(\overline X_n/\overline Y_n-\theta\right)}{\sqrt{\widehat V_n}}
\rightsquigarrow
\mathcal N(0,1),
$$

and a pointwise asymptotic $(1-\alpha)$ confidence interval is

$$
\left[
\frac{\overline X_n}{\overline Y_n}
-z_{1-\alpha/2}\sqrt{\frac{\widehat V_n}{n}},
\frac{\overline X_n}{\overline Y_n}
+z_{1-\alpha/2}\sqrt{\frac{\widehat V_n}{n}}
\right].
$$

The qualification $\mu_Y\neq0$ is substantive. This argument is not uniform over parameter values whose denominator can approach zero.

> [!IMPORTANT]
> **Board work 2 — From a limit theorem to an interval**
>
> Use the ratio of means example to build the inference chain without skipping a link:
>
> 1. begin with a linear or asymptotically Normal representation;
> 2. identify the asymptotic variance $V$;
> 3. prove $\widehat V_n\to_{\mathbb{P}}V$;
> 4. apply continuous mapping to $\sqrt{\widehat V_n}$ and its reciprocal;
> 5. apply Slutsky to obtain the pivot; and
> 6. invert the pivot to obtain the confidence interval.
>
> At each step, state whether the conclusion is exact, pointwise asymptotic, or uniform asymptotic.

<a id="w6-stop-5"></a>

### 2.4 Asymptotic size and power

For a sequence of tests $\phi_n$ with power functions

$$
\beta_n(\theta)=\mathbb E_\theta[\phi_n(X)],
$$

**pointwise asymptotic level $\alpha$** means

$$
\limsup_{n\to\infty}\beta_n(\theta)\leq\alpha
\qquad\text{for every fixed }\theta\in\Theta_0.
$$

The stronger **uniform asymptotic level $\alpha$** condition is

$$
\limsup_{n\to\infty}
\sup_{\theta\in\Theta_0}\beta_n(\theta)
\leq\alpha.
$$

When the limit exists, $\lim_n\sup_{\theta\in\Theta_0}\beta_n(\theta)$ is the limiting size. The test sequence is **pointwise consistent** against $\Theta_1$ if $\beta_n(\theta)\to1$ for every fixed $\theta\in\Theta_1$. Pointwise size control and pointwise consistency do not imply uniform performance over parameter sequences approaching a boundary.

### 2.5 Pointwise versus uniform coverage

Pointwise asymptotic exactness means that for each fixed $\theta$,

$$
\mathbb{P}_\theta\lbrace\theta\in C_n(X)\rbrace\to1-\alpha.
$$

Uniform asymptotic validity requires the stronger worst case statement

$$
\liminf_{n\to\infty}
\inf_{\theta\in\Theta}
\mathbb{P}_\theta\lbrace\theta\in C_n(X)\rbrace
\geq1-\alpha.
$$

The still stronger condition

$$
\sup_{\theta\in\Theta}
\left|
\mathbb{P}_\theta\lbrace\theta\in C_n(X)\rbrace-(1-\alpha)
\right|
\to0
$$

gives uniform asymptotic exactness rather than merely preventing undercoverage.

A pointwise theorem allows the approximation to deteriorate along sequences $\theta_n$ that change with $n$. This distinction will matter for weak identification, boundaries, nonparametric smoothing, and semiparametric procedures.

> [!TIP]
> **AI interaction 2 — Audit an asymptotic interval**
>
> Copy the prompt below into an AI interface and audit the response.

```text
A student writes: "Because sqrt(n)(theta_hat-theta) is asymptotically Normal,
theta_hat plus or minus 1.96 times its standard error is an exact 95 percent
confidence interval for every n and uniformly over theta."

Identify every unjustified step. State the additional requirements for
studentization, distinguish exact from asymptotic coverage, and distinguish
pointwise from uniform validity. Rewrite the conclusion as a correct theorem.
```

### Checkpoint 2

1. What must be proved before replacing $V$ by $\widehat V_n$?
2. Does $\widehat\theta_n-\theta_0=O_{\mathbb{P}}(n^{-1/2})$ imply consistency?
3. Does consistency imply a $\sqrt n$ rate?
4. Does pointwise asymptotic level control the limiting supremum of rejection probability over a composite null?
5. Does a pointwise asymptotic confidence interval control worst case coverage over $\Theta$?
6. What additional condition is needed before reading asymptotic bias or variance from a weak limit?

<a id="w6-stop-6"></a>

## 3. Bridge to general estimation

A reusable asymptotic inference argument has five logically separate steps:

1. show that the estimator is consistent for an identified population target;
2. obtain a first order expansion at an explicit rate;
3. apply an appropriate limit theorem to the leading random term;
4. estimate the asymptotic variance consistently; and
5. state whether the resulting test or confidence set is exact, pointwise asymptotic, or uniformly valid.

Week 9 develops general M- and Z-estimation theorems that automate parts of this chain. The global argument establishing consistency and the local argument establishing asymptotic normality remain distinct.

### Checkpoint 3

For an estimator you encountered in Weeks 2--4, identify which of the five steps above are exact, which require a limit theorem, and which would need additional work before claiming uniform validity.

<a id="w6-stop-7"></a>

## 4. A uniform law of large numbers

The ordinary law of large numbers gives convergence at each fixed parameter value. Estimation by optimization evaluates a sample criterion at a random, data dependent parameter, so pointwise convergence is generally insufficient. The following theorem gives one practical route to uniform convergence.

Let $W_1,\ldots,W_n$ be i.i.d., let $(\Theta,d)$ be a compact metric space, and let $q(w,\theta)$ be measurable in $w$ for each $\theta$. Define

$$
Q_n(\theta)=\frac1n\sum_{i=1}^n q(W_i,\theta),
\qquad
Q(\theta)=\mathbb E[q(W,\theta)].
$$

Suppose:

1. for almost every $w$, the map $\theta\mapsto q(w,\theta)$ is continuous on $\Theta$;
2. there is a measurable envelope $F$ such that

$$
\sup_{\theta\in\Theta}|q(W,\theta)|\leq F(W),
\qquad
\mathbb E[F(W)]<\infty;
$$

3. the suprema used below are measurable.

Then

$$
\sup_{\theta\in\Theta}|Q_n(\theta)-Q(\theta)|
\to_{\mathbb{P}}0.
$$

### Proof: three ingredients

The argument upgrades pointwise laws of large numbers through discretization, finite dimensional convergence, and oscillation control.

#### Step 1: Discretization

Fix $\delta>0$. Compactness of $\Theta$ gives a finite net of radius $\delta$, with points $\theta_1,\ldots,\theta_J$: for every $\theta\in\Theta$, some $\theta_j$ satisfies $d(\theta,\theta_j)<\delta$. Define the global oscillation modulus

$$
\omega_\delta(w) =
\sup_{\substack{\theta,\vartheta\in\Theta\\ d(\theta,\vartheta)<\delta}}
|q(w,\theta)-q(w,\vartheta)|.
$$

For a center $\theta_j$ within $\delta$ of $\theta$,

$$
|Q_n(\theta)-Q(\theta)|
\leq
|Q_n(\theta_j)-Q(\theta_j)|
+\frac1n\sum_{i=1}^n\omega_\delta(W_i)
+\mathbb E[\omega_\delta(W)].
$$

Consequently,

$$
\sup_{\theta\in\Theta}|Q_n(\theta)-Q(\theta)|
\leq
\max_{1\leq j\leq J}|Q_n(\theta_j)-Q(\theta_j)|
+\frac1n\sum_{i=1}^n\omega_\delta(W_i)
+\mathbb E[\omega_\delta(W)].
$$

This inequality separates the problem into behavior on a finite grid and behavior between grid points.

#### Step 2: Finite dimensional convergence

For fixed $\delta$, the ordinary law of large numbers applies at each of the finitely many centers:

$$
Q_n(\theta_j)\to_{\mathbb{P}}Q(\theta_j),
\qquad j=1,\ldots,J.
$$

Because $J$ is finite,

$$
\max_{1\leq j\leq J}|Q_n(\theta_j)-Q(\theta_j)|
\to_{\mathbb{P}}0.
$$

This is the only step that uses pointwise convergence directly.

#### Step 3: Oscillation control

For almost every $w$, continuity of $q(w,\cdot)$ on compact $\Theta$ implies uniform continuity. Hence $\omega_\delta(W)\to0$ almost surely as $\delta\downarrow0$. Also,

$$
0\leq\omega_\delta(W)\leq2F(W),
$$

so dominated convergence gives $\mathbb E[\omega_\delta(W)]\to0$. For each fixed $\delta$, the ordinary law of large numbers also gives

$$
\frac1n\sum_{i=1}^n\omega_\delta(W_i)
\to_{\mathbb{P}}
\mathbb E[\omega_\delta(W)].
$$

Choose $\delta$ so that the expected oscillation is arbitrarily small, and then let $n\to\infty$. The discretization bound, finite dimensional convergence, and oscillation control together yield

$$
\sup_{\theta\in\Theta}|Q_n(\theta)-Q(\theta)|
\to_{\mathbb{P}}0.
$$

This theorem is intentionally limited to a compact, continuously indexed class with an integrable envelope. Week 9 applies it to criterion functions and then combines uniform convergence with identification and separation to prove consistency of M- and Z-estimators. Week 10 reuses the same global logic for regression and two step procedures. More general empirical process machinery is not required here.

## 5. Recap

- The delta method is a Taylor expansion combined with stochastic order control.
- A zero first derivative can change both the rate and the limiting distribution.
- Consistency, rate, asymptotic distribution, and convergence of moments are distinct claims.
- Asymptotic bias, variance, MSE, and relative efficiency require moment control in addition to weak convergence.
- Studentization requires both a limiting distribution and a consistent variance estimator.
- Pointwise asymptotic size or coverage is weaker than uniform validity.
- Exact finite sample inference and asymptotic inference must be labeled separately.
- These tools provide the probability foundation for the M- and Z-estimation framework in Week 9.
- Compactness, continuity, and an integrable envelope convert finitely many ordinary LLNs and local oscillation bounds into a uniform LLN.

## Notation introduced this week

- $r_n$: convergence rate sequence.
- $G(\theta)$: Jacobian of a transformation.
- $b$, $V$, and $M$: first order asymptotic bias, variance, and MSE on the chosen scale.
- $\mathrm{ARE}(1,2)$: asymptotic relative efficiency under the stated variance ratio convention.
- $\widehat V_n$: consistent estimator of the asymptotic variance $V$.
- $\beta_n(\theta)$: power function of test number $n$ in an asymptotic sequence.
- $Q_n(\theta)$ and $Q(\theta)$: sample and population criterion functions in the uniform LLN toolkit.

## References

- Jun Shao, *Mathematical Statistics*, Springer, 2003.
- Robert W. Keener, *Theoretical Statistics: Topics for a Core Course*, Springer, 2010.
- Aad van der Vaart, *Asymptotic Statistics*, Cambridge University Press, 1998.
