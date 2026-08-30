# ORF 524: Statistical Theory and Methods

## Week 9: M- and Z-estimation

**Status:** Ready for Fall 2026  
**Last updated:** August 30, 2026  
**Meetings:** Monday, October 26, and Wednesday, October 28, 2026

<!-- Legacy sources: 11 ORF 524/lectures/ORF524-Fall2025-Lecture-Week-06.tex and
ORF524-Fall2025-Lecture-Week-06_Lecture.pdf; 09 ORF 524/lectures/
ORF524-Fall2024-Lecture-Week-06_Notes.pdf. -->

## Central question

How can one argument establish consistency, asymptotic normality, standard errors, and inference for many estimators defined by optimization or estimating equations?

## Learning goals

By the end of the week, you should be able to:

1. represent familiar procedures as M-estimators, Z-estimators, or both;
2. distinguish population identification from sample optimization or root finding;
3. verify a usable uniform LLN and prove consistency using uniform convergence and separation;
4. derive asymptotic linear representations for smooth M- and Z-estimators;
5. construct and interpret the sandwich variance estimator; and
6. compare Wald, score, and likelihood ratio tests under regularity and identify failure modes.

## In-class route

This is one continuous route through the week, not a division into class meetings. Proceed in order, stop wherever the discussion ends, and resume at the next unfinished stop. The linked sections contain the full student facing exposition; the route identifies the concepts and activities to foreground live.

| Stop | Live focus | Mode |
|---|---|---|
| **W9.1** | [M- and Z-estimator representations](#w9-stop-1) | Discuss + AI representation audit + Checkpoint 1 |
| **W9.2** | [Population identification](#w9-stop-2) | Discuss + Checkpoint 2 |
| **W9.3** | [Uniform convergence and consistency](#w9-stop-3) | Board work + AI quantifier audit + Checkpoint 3 |
| **W9.4** | [Local linearization and influence functions](#w9-stop-4) | Board work + AI mean value audit + Checkpoint 4 |
| **W9.5** | [Sandwich variance and canonical examples](#w9-stop-5) | Discuss + board work + Checkpoint 5 |
| **W9.6** | [Wald, score, and likelihood ratio inference](#w9-stop-6) | Geometric comparison + AI dialogue + Checkpoint 6 |
| **W9.7** | [Synthesis and bridge to applications](#w9-stop-7) | Transfer + Checkpoint 7 |

Stop W9.3 should recall the Week 6 uniform LLN without reproving it, keep the wandering minimizer counterexample visible, and then use uniform convergence in the argmin theorem. Stop W9.4 is the principal handwritten board derivation: the random intermediate derivative, asymptotic linear representation, and influence function form should be developed as one chain.

## How to use this chapter

**Prepare:** Review stochastic orders, the uniform LLN and pointwise versus uniform distinction from Week 6, the multivariate CLT, and the mean value theorem. Read Sections 1 and 2 and attempt Checkpoint 1.

**In class:** Follow the continuous route above through estimator representations, identification, consistency, local linearization, sandwich variance, and the three classical large sample tests. Prepared responses supply an alternative for every live prompt if AI is unavailable.

**Review:** For every theorem, separate the global conditions that locate the estimator from the local conditions that determine its distribution. Re derive the sample mean and regression examples from the general formulas.

**Practice:** The matching assignment complements this chapter with ungraded problems for consolidating the material, extending selected ideas, and preparing for cumulative review.

**Prerequisites:** Weeks 1--6; identification; maximum likelihood and method of moments; LLN, CLT, Slutsky, stochastic orders, matrix differentiation, and quadratic forms.

## Full chapter map

1. [Two general estimator representations](#1-two-general-estimator-representations)
2. [Population identification](#2-population-identification)
3. [Consistency](#3-consistency)
4. [Asymptotic linearity and normality](#4-asymptotic-linearity-and-normality)
5. [Sandwich variance and examples](#5-sandwich-variance-and-examples)
6. [Wald, score, and likelihood ratio inference](#6-wald-score-and-likelihood-ratio-inference)
7. [Synthesis and bridge to applications](#7-synthesis-and-bridge-to-applications)
8. [Recap](#8-recap)

<a id="w9-stop-1"></a>

## 1. Two general estimator representations

Let $W_1,\ldots,W_n$ be i.i.d. observations from $\mathbb{P}_0$, and write

$$
\mathbb{E}_n f=\frac1n\sum_{i=1}^n f(W_i).
$$

The parameter $\theta_0$ belongs to $\Theta\subseteq\mathbb R^k$.

### 1.1 M-estimators

Given a criterion $\rho(w,\theta)$, define

$$
Q_n(\theta)=\mathbb{E}_n\rho(W,\theta).
$$

An **M-estimator** approximately minimizes the sample criterion: for some $r_n\geq0$ with $r_n=o_{\mathbb{P}}(1)$,

$$
Q_n(\widehat\theta_n)
\leq
\inf_{\theta\in\Theta}Q_n(\theta)+r_n.
$$

Exact minimization is a special case. Maximization problems are included by changing the sign of the criterion.

Examples include:

- maximum likelihood, with $\rho(w,\theta)=-\log f(w;\theta)$;
- least squares, with $\rho((y,x),\beta)=(y-x'\beta)^2$; and
- quantile estimation, with the check loss $\rho_\tau(u)=u\lbrace\tau-\mathbf 1(u<0)\rbrace$.

### 1.2 Z-estimators

Given an estimating function $\psi(w,\theta)\in\mathbb R^k$, define

$$
\Psi_n(\theta)=\mathbb{E}_n\psi(W,\theta).
$$

A **Z-estimator** approximately solves

$$
\Psi_n(\widehat\theta_n)=o_{\mathbb{P}}(n^{-1/2})
$$

for asymptotic linearity, or at least $o_{\mathbb{P}}(1)$ for a consistency argument.

Examples include:

- method of moments;
- the sample mean, with $\psi(w,\theta)=w-\theta$;
- least squares normal equations, with $\psi((y,x),\beta)=x(y-x'\beta)$; and
- the likelihood score in a regular interior problem.

### 1.3 Relationship between the representations

If $Q_n$ is differentiable and $\widehat\theta_n$ is an interior minimizer, then it may also be a Z-estimator with

$$
\psi(w,\theta)=\nabla_\theta\rho(w,\theta).
$$

The representations are not interchangeable without conditions. Boundary optima, nonunique solutions, and nonsmooth criteria may not satisfy an ordinary gradient equation. Conversely, an estimating equation need not be the derivative of any scalar criterion.

### 1.4 Method of moments as Z-estimation and extremum estimation

Suppose the parameter is identified by

$$
\mathbb E[g(W)]=m(\theta_0).
$$

The method of moments estimator is a Z-estimator with

$$
\psi(W,\theta)=g(W)-m(\theta).
$$

In the broader extremum estimation sense, it can also minimize

$$
Q_n(\theta)
=\lbrace\mathbb{E}_n g-m(\theta)\rbrace'
\Omega_n
\lbrace\mathbb{E}_n g-m(\theta)\rbrace,
$$

where $\Omega_n$ is positive semidefinite. In a just identified problem, a zero of the sample moment attains the minimum zero. In an overidentified problem, an exact zero may not exist and the weighting matrix affects the estimator. Notice that this quadratic criterion need not have the separable form $\mathbb{E}_n\rho(W,\theta)$ used in Section 1.1.

Moment and criterion methods can require less distributional structure than a full parametric likelihood. When the likelihood is correctly specified, its additional restrictions can improve efficiency; under misspecification, they can instead change the target. Robustness and precision are therefore properties to establish in a particular model, not automatic labels attached to an estimator's name.

### Checkpoint 1

Classify each procedure as M, Z, both under conditions, or neither from the information given.

1. sample mean;
2. sample median;
3. maximum likelihood estimator in $\mathsf{Uniform}(0,\theta)$;
4. ordinary least squares;
5. a root of an arbitrary vector estimating equation; and
6. a just identified method of moments estimator.

> [!TIP]
> **AI interaction 1 — Same estimator, different theorem**
>
> Copy the prompt below into an AI interface and audit the response.

```text
Give two representations of the sample mean: one as an M-estimator and one as
a Z-estimator. Then compare the assumptions a generic argmin consistency
theorem and a generic estimating equation theorem would require. Explain why
two algebraically equivalent representations of the same estimator need not
satisfy the assumptions of the same theorem.
```

**Audit question:** Does the response infer theorem assumptions from the estimator's name rather than from its chosen representation?

<a id="w9-stop-2"></a>

## 2. Population identification

### 2.1 M-identification

The population criterion is

$$
Q(\theta)=\mathbb{P}_0\rho(W,\theta)=\mathbb E_{\mathbb{P}_0}[\rho(W,\theta)].
$$

For M-estimation, the target must be identified as a population minimizer. A useful condition is that $\theta_0$ is the unique, well separated minimizer:

$$
Q(\theta_0)<Q(\theta)
\quad\text{for }\theta\neq\theta_0,
$$

and, for every $\varepsilon>0$,

$$
\inf_{\theta:\Vert\theta-\theta_0\Vert\geq\varepsilon}
\lbrace Q(\theta)-Q(\theta_0)\rbrace>0.
$$

Uniqueness alone need not provide separation on a noncompact or irregular parameter space.

### 2.2 Z-identification

The population estimating equation is

$$
\Psi(\theta)=\mathbb{P}_0\psi(W,\theta).
$$

For Z-estimation, require

$$
\Psi(\theta_0)=0
$$

and separation of the root:

$$
\inf_{\theta:\Vert\theta-\theta_0\Vert\geq\varepsilon}
\Vert\Psi(\theta)\Vert>0
$$

for every $\varepsilon>0$.

These are population identification conditions. A sample minimizer or sample root does not establish that the corresponding population target is identified.

### Checkpoint 2

1. Why is a unique minimizer not automatically well separated?
2. Can $Q$ identify $\theta_0$ even when the statistical model's full parameterization is redundant?
3. Can $\Psi(\theta)=0$ have multiple roots even when each component moment is correctly specified?

<a id="w9-stop-3"></a>

## 3. Consistency

### 3.1 Uniform convergence is the bridge

Pointwise convergence

$$
Q_n(\theta)\to_{\mathbb{P}}Q(\theta)
\quad\text{for each fixed }\theta
$$

does not control the random point $\widehat\theta_n$ selected by optimization. The standard bridge is a uniform law of large numbers:

$$
\sup_{\theta\in\Theta}|Q_n(\theta)-Q(\theta)|\to_{\mathbb{P}}0.
$$

Uniform convergence prevents the sample criterion from developing a deep, narrow, moving minimum that is invisible at every fixed $\theta$.

### 3.2 Reusing the Week 6 uniform LLN

[Week 6 proved a usable uniform LLN](week-06.md#4-a-uniform-law-of-large-numbers) for i.i.d. observations on a compact parameter space under almost sure continuity, an integrable envelope, and measurability. Its three proof ingredients were discretization, finite dimensional convergence, and oscillation control. Here we treat that theorem as an available tool: verify its conditions for the chosen criterion, obtain the uniform convergence displayed above, and then combine it with population separation.

### 3.3 Counterexample: pointwise convergence with wandering minimizers

On $\Theta=[0,1]$, let $Q(\theta)=\theta^2$ and

$$
Q_n(\theta)
=\theta^2
-2\mathbf 1\left\lbrace
\left|\theta-\left(1-\frac1n\right)\right|\leq\frac1{n^2}
\right\rbrace.
$$

For every fixed $\theta$, eventually the moving spike misses $\theta$, so $Q_n(\theta)\to Q(\theta)$. Yet a minimizer of $Q_n$ lies near $1$, not near the unique population minimizer $0$.

### 3.4 Argmin consistency theorem

Suppose:

1. $\theta_0$ is a well separated minimizer of $Q$;
2. $\sup_{\theta\in\Theta}|Q_n(\theta)-Q(\theta)|\to_{\mathbb{P}}0$; and
3. $Q_n(\widehat\theta_n)\leq\inf_{\theta\in\Theta}Q_n(\theta)+r_n$ for some $r_n\geq0$ with $r_n=o_{\mathbb{P}}(1)$.

Then

$$
\widehat\theta_n\to_{\mathbb{P}}\theta_0.
$$

### Proof

Fix $\varepsilon>0$ and define the separation gap

$$
\eta_\varepsilon
=\inf_{\Vert\theta-\theta_0\Vert\geq\varepsilon}
\lbrace Q(\theta)-Q(\theta_0)\rbrace>0.
$$

On the event

$$
\sup_{\theta\in\Theta}|Q_n(\theta)-Q(\theta)|<\eta_\varepsilon/3
$$

and when the optimization error is below $\eta_\varepsilon/3$, every $\theta$ outside the $\varepsilon$ ball satisfies

$$
Q_n(\theta)
>Q(\theta)-\eta_\varepsilon/3
\geq Q(\theta_0)+2\eta_\varepsilon/3
>Q_n(\theta_0)+\eta_\varepsilon/3.
$$

An approximate minimizer cannot lie there. The probability of the uniform convergence event tends to one, and $r_n=o_{\mathbb{P}}(1)$, so

$$
\mathbb{P}(\Vert\widehat\theta_n-\theta_0\Vert\geq\varepsilon)\to0.
$$

### 3.5 Z-estimator consistency

The parallel result replaces criteria by estimating equations. If

$$
\sup_{\theta\in\Theta}
\Vert\Psi_n(\theta)-\Psi(\theta)\Vert
\to_{\mathbb{P}}0,
$$

$\theta_0$ is a well separated zero of $\Psi$, and $\Vert\Psi_n(\widehat\theta_n)\Vert=o_{\mathbb{P}}(1)$, then

$$
\widehat\theta_n\to_{\mathbb{P}}\theta_0.
$$

> [!IMPORTANT]
> **Board work 1 — Z-estimator consistency**
>
> Adapt the argmin proof to the Z-estimator theorem. Identify the event on which every point outside an $\varepsilon$ ball has sample moment norm bounded away from zero.

> [!TIP]
> **AI interaction 2 — Find the quantifier error**
>
> Copy the prompt below into an AI interface and audit the response.

```text
A proof of M-estimator consistency says: "For each theta, Q_n(theta)
converges to Q(theta). Therefore Q_n(theta_hat) converges to Q(theta_hat), so
the minimizer theta_hat converges to the unique minimizer of Q."

Identify the first invalid step. Explain the mismatch of quantifiers, use a
moving spike counterexample, and state a corrected argmin consistency theorem
with separation, uniform convergence, and approximate optimization.
```

### Checkpoint 3

1. Which condition is global: separation, uniform convergence, or nonsingularity of a derivative?
2. Does consistency require an exact sample minimizer?
3. Can computational error be asymptotically negligible?
4. Why does compactness often appear in sufficient conditions for uniform convergence and separation, even though it is not itself the conclusion?

<a id="w9-stop-4"></a>

## 4. Asymptotic linearity and normality

Consistency localizes the estimator. A mean value argument then determines its first order distribution.

### 4.1 Smooth Z-estimator mean value argument

Suppose $\widehat\theta_n\to_{\mathbb{P}}\theta_0$ and $\Psi_n(\widehat\theta_n)=o_{\mathbb{P}}(n^{-1/2})$. Define the derivative of the population estimating equation by

$$
\Psi_\theta(\theta)
:=
\frac{\partial}{\partial\theta'}\Psi(\theta) =
\frac{\partial}{\partial\theta'}
\mathbb E_{\mathbb{P}_0}[\psi(W,\theta)].
$$

When differentiation may be moved through the expectation,

$$
\Psi_\theta(\theta) =
\mathbb E_{\mathbb{P}_0}
\left[
\frac{\partial}{\partial\theta'}\psi(W,\theta)
\right].
$$

Define the corresponding sample derivative by

$$
\Psi_{\theta,n}(\theta)
:=
\frac{\partial}{\partial\theta'}\Psi_n(\theta) =
\mathbb{E}_n
\left[
\frac{\partial}{\partial\theta'}\psi(W,\theta)
\right].
$$

Assume $\theta_0$ is interior and the estimating equation is continuously differentiable on a convex neighborhood of $\theta_0$. Apply the scalar mean value theorem to each coordinate $\Psi_{n,j}$. For every $j=1,\ldots,k$, there is an intermediate point $\widetilde\theta_{n,j}$ on the segment from $\theta_0$ to $\widehat\theta_n$ such that

$$
\Psi_{n,j}(\widehat\theta_n)-\Psi_{n,j}(\theta_0) =
\Psi_{\theta,n,j}(\widetilde\theta_{n,j})
(\widehat\theta_n-\theta_0),
$$

where $\Psi_{\theta,n,j}$ denotes the $j$-th row of the derivative matrix. The intermediate point may differ across coordinates. Stacking these $k$ equations gives

$$
\Psi_n(\widehat\theta_n)-\Psi_n(\theta_0) =
\widetilde\Psi_{\theta,n}
(\widehat\theta_n-\theta_0),
$$

where row $j$ of $\widetilde\Psi_{\theta,n}$ is $\Psi_{\theta,n,j}(\widetilde\theta_{n,j})$. There need not be one common intermediate point, but the rowwise identity is all the argument requires.

The remaining issue is to show that the derivative evaluated along this random segment is consistent. Suppose that for some neighborhood $N_\delta=\lbrace\theta:\Vert\theta-\theta_0\Vert\leq\delta\rbrace$,

$$
\sup_{\theta\in N_\delta}
\Vert\Psi_{\theta,n}(\theta)-\Psi_\theta(\theta)\Vert
\to_{\mathbb{P}}0,
$$

and that $\Psi_\theta(\theta)$ is continuous at $\theta_0$. The displayed convergence is a local uniform LLN. A convenient sufficient condition is that $\partial\psi(W,\theta)/\partial\theta'$ be almost surely continuous in $\theta$ on $N_\delta$ and dominated there by an integrable envelope. Consistency places every intermediate point inside $N_\delta$ with probability approaching one. The ULLN controls all the rows at once, and continuity then gives

$$
\widetilde\Psi_{\theta,n}
\to_{\mathbb{P}}
\Psi_\theta(\theta_0).
$$

Assume $\Psi_\theta(\theta_0)$ is nonsingular. Then $\widetilde\Psi_{\theta,n}$ is invertible with probability approaching one, and the mean value identity and approximate root condition yield

$$
\sqrt n(\widehat\theta_n-\theta_0)
=-\widetilde\Psi_{\theta,n}^{-1}\sqrt n\Psi_n(\theta_0)
+o_{\mathbb{P}}(1).
$$

Slutsky's theorem now gives the **asymptotic linear representation**

$$
\sqrt n(\widehat\theta_n-\theta_0)
=-[\Psi_\theta(\theta_0)]^{-1}
\sqrt n\Psi_n(\theta_0)
+o_{\mathbb{P}}(1).
$$

If

$$
\sqrt n\Psi_n(\theta_0)
=\frac1{\sqrt n}\sum_{i=1}^n\psi(W_i,\theta_0)
\rightsquigarrow
\mathcal N(0,B),
$$

where

$$
B=\mathbb V_{\mathbb{P}_0}[\psi(W,\theta_0)],
$$

then

$$
\sqrt n(\widehat\theta_n-\theta_0)
\rightsquigarrow
\mathcal N(0,V),
$$

with

$$
V
=[\Psi_\theta(\theta_0)]^{-1}
B
[\Psi_\theta(\theta_0)]^{-\prime}.
$$

### Proof map

The argument has four logically distinct inputs:

1. consistency places the random segment from $\theta_0$ to $\widehat\theta_n$ in a fixed local neighborhood;
2. a local uniform LLN makes the sample derivative consistent uniformly along that segment;
3. nonsingularity of $\Psi_\theta(\theta_0)$ permits solving for estimation error; and
4. a CLT determines the distribution of the empirical moment at $\theta_0$.

Removing any one input changes or destroys the conclusion.

### 4.2 Smooth M-estimator theorem

The corresponding result for a smooth M-estimator follows by treating the gradient of the criterion as an estimating equation. Suppose $\theta_0$ is an interior population minimizer with $\mathbb E[\nabla_\theta\rho(W,\theta_0)]=0$, and $\widehat\theta_n\to_{\mathbb{P}}\theta_0$ is an interior approximate sample minimizer satisfying

$$
\mathbb{E}_n\nabla_\theta\rho(W,\widehat\theta_n)
=o_{\mathbb{P}}(n^{-1/2}).
$$

Let

$$
H=\mathbb E[\nabla_{\theta\theta'}^2\rho(W,\theta_0)],
\qquad
\Sigma=\mathbb V[\nabla_\theta\rho(W,\theta_0)].
$$

Write

$$
H_n(\theta)=\mathbb{E}_n\nabla_{\theta\theta'}^2\rho(W,\theta),
\qquad
H(\theta)=\mathbb E[\nabla_{\theta\theta'}^2\rho(W,\theta)].
$$

The required Hessian consistency is obtained from the local uniform LLN

$$
\sup_{\theta\in N_\delta}
\Vert H_n(\theta)-H(\theta)\Vert
\to_{\mathbb{P}}0.
$$

As above, almost sure continuity of the Hessian in $\theta$ and an integrable local envelope are convenient sufficient conditions. Applying the mean value theorem coordinate by coordinate to the criterion gradient evaluates different rows of the sample Hessian at possibly different intermediate points. Consistency puts all those points inside $N_\delta$, and the uniform LLN plus continuity makes the resulting rowwise matrix converge to $H=H(\theta_0)$. If $H$ is nonsingular and a CLT holds for the criterion gradient, then

$$
\sqrt n(\widehat\theta_n-\theta_0)
=-H^{-1}\frac1{\sqrt n}
\sum_{i=1}^n\nabla_\theta\rho(W_i,\theta_0)
+o_{\mathbb{P}}(1)
$$

and therefore

$$
\sqrt n(\widehat\theta_n-\theta_0)
\rightsquigarrow
\mathcal N(0,H^{-1}\Sigma H^{-\prime}).
$$

Twice continuous differentiability with a dominated local Hessian supplies convenient sufficient conditions, but it is stronger than necessary. The sample median, for example, is an asymptotically normal M-estimator at the $\sqrt n$ rate even though the absolute loss criterion is nonsmooth. For maximum likelihood, $\rho(W,\theta)=-\log f(W;\theta)$, so this same argument applies to the score and its derivative, the negative sample Hessian.

### 4.3 Influence function form

Define

$$
\varphi(W)
=-[\Psi_\theta(\theta_0)]^{-1}\psi(W,\theta_0).
$$

Then

$$
\widehat\theta_n-\theta_0
=\frac1n\sum_{i=1}^n\varphi(W_i)+o_{\mathbb{P}}(n^{-1/2}).
$$

The function $\varphi$ describes the first order contribution of one observation. This representation will become central in the second half of the course.

### Failure modes

- If $\Psi_\theta(\theta_0)$ is singular or nearly singular, the equation does not determine a stable local approximation at the $\sqrt n$ rate.
- If the estimating equation is nonsmooth, the mean value argument may not apply.
- If the estimator lies on a boundary, the limit may be truncated or otherwise non-Gaussian.
- If $\psi$ has infinite variance, the classical CLT and sandwich formula may fail.
- If identification weakens with $n$, pointwise regular theory may be misleading.

> [!IMPORTANT]
> **Board work 2 — From consistency to asymptotic linearity**
>
> Apply the mean value theorem coordinate by coordinate to a smooth Z-estimating equation between $\theta_0$ and $\widehat\theta_n$. Allow the intermediate point to differ across rows, explain how consistency and a local uniform LLN control all the resulting derivatives, and then recover the asymptotic linear and sandwich variance formulas.

> [!TIP]
> **AI interaction 3 — Audit the linearization**
>
> Copy the prompt below into an AI interface and audit the response.

```text
Derive the asymptotic linear representation of a smooth Z-estimator using the
mean value theorem, not a second order Taylor remainder. For a vector parameter,
apply the scalar theorem coordinate by coordinate and allow a different
intermediate point for each row of the derivative matrix. Explain how
consistency and a local uniform LLN control all those intermediate derivatives.
Then label the steps using nonsingularity and the CLT. Do not hide the argument
under "standard regularity."
```

### Checkpoint 4

1. Why is a pointwise LLN for the derivative or Hessian not enough along the random segment between $\widehat\theta_n$ and $\theta_0$?
2. What does the minus sign in $-[\Psi_\theta(\theta_0)]^{-1}$ do to the asymptotic variance?
3. Does asymptotic normality imply that the estimator is unbiased in finite samples?
4. How do consistency, a local uniform LLN, and continuity together make the rowwise intermediate derivative matrix converge to $\Psi_\theta(\theta_0)$?

<a id="w9-stop-5"></a>

## 5. Sandwich variance and examples

### 5.1 Estimating the variance

Recall the sample derivative

$$
\Psi_{\theta,n}(\theta)
:=
\frac{\partial}{\partial\theta'}\Psi_n(\theta) =
\mathbb{E}_n
\left[
\frac{\partial}{\partial\theta'}\psi(W,\theta)
\right],
$$

and let

$$
\widehat B_n
=\mathbb{E}_n
\left[
\psi(W,\widehat\theta_n)
\psi(W,\widehat\theta_n)'
\right].
$$

Under consistency and the required laws of large numbers,

$$
\Psi_{\theta,n}(\widehat\theta_n)
\to_{\mathbb{P}}
\Psi_\theta(\theta_0),
\qquad
\widehat B_n\to_{\mathbb{P}}B.
$$

Therefore

$$
\widehat V_n
=[\Psi_{\theta,n}(\widehat\theta_n)]^{-1}
\widehat B_n
[\Psi_{\theta,n}(\widehat\theta_n)]^{-\prime}
\to_{\mathbb{P}}V.
$$

This is the **sandwich variance**. The estimated covariance matrix of $\widehat\theta_n$ itself is $\widehat V_n/n$.

### 5.2 Sample mean

For $\psi(W,\theta)=W-\theta$,

$$
\Psi_\theta(\mu)=-1,
\qquad
B=\mathbb V[W]=\sigma^2,
$$

so

$$
V=\sigma^2.
$$

Thus

$$
\sqrt n(\overline W_n-\mu)
\rightsquigarrow\mathcal N(0,\sigma^2).
$$

### 5.3 Best linear predictor and least squares

Let $W=(Y,X)$ and define $\beta_0$ by

$$
\mathbb E[X(Y-X'\beta_0)]=0.
$$

The consistency argument now makes the Section 3 toolkit concrete. Let

$$
Q_n(\beta)=\mathbb{E}_n(Y-X'\beta)^2,
\qquad
Q(\beta)=\mathbb E[(Y-X'\beta)^2],
$$

and suppose optimization is over a compact set $\mathcal B\subset\lbrace\beta:\Vert\beta\Vert\leq K\rbrace$ containing $\beta_0$. The criterion is continuous in $\beta$, and

$$
(Y-X'\beta)^2
\leq2Y^2+2K^2\Vert X\Vert^2.
$$

Thus $\mathbb E[Y^2+\Vert X\Vert^2]<\infty$ supplies an integrable envelope, and the Week 6 uniform LLN gives $\sup_{\beta\in\mathcal B}|Q_n(\beta)-Q(\beta)|\to_{\mathbb{P}}0$. If $\mathbb E[XX']$ is positive definite, $\beta_0$ is the unique population minimizer; continuity and compactness give separation. The argmin theorem therefore yields consistency of the least squares estimator. If $\beta_0$ is in the interior of $\mathcal B$, the estimator is also interior with probability approaching one, and its first order condition is then valid.

The sample estimator solves

$$
\mathbb{E}_n[X(Y-X'\widehat\beta_n)]=0.
$$

Writing $\varepsilon=Y-X'\beta_0$,

$$
\Psi_\beta(\beta_0)=-\mathbb E[XX'],
\qquad
B=\mathbb E[XX'\varepsilon^2].
$$

If $\mathbb E[XX']$ is nonsingular,

$$
V
=\mathbb E[XX']^{-1}
\mathbb E[XX'\varepsilon^2]
\mathbb E[XX']^{-1}.
$$

This formula does not require homoskedasticity. Under $\mathbb E[\varepsilon^2\mid X]=\sigma^2$, it simplifies to

$$
V=\sigma^2\mathbb E[XX']^{-1}.
$$

The target is the best linear predictor coefficient under squared loss. It need not describe the full conditional expectation unless that function is linear.

### 5.4 Regular maximum likelihood

Applying the argmin theorem in Section 3 to negative log likelihood gives consistency when the population likelihood uniquely identifies $\theta_0$, the sample log likelihood converges uniformly, and the optimizer is asymptotically exact. Pointwise likelihood convergence alone is insufficient.

For a correctly specified regular likelihood, take $\psi$ to be the score. The information equality gives

$$
B=\mathcal I_1(\theta_0),
\qquad
\Psi_\theta(\theta_0)=-\mathcal I_1(\theta_0),
$$

and therefore

$$
V=\mathcal I_1(\theta_0)^{-1}.
$$

Writing $s(W,\theta)=\nabla_\theta\log f(W;\theta)$, two familiar plug in estimators of Fisher information are

$$
\widehat{\mathcal I}_{\mathrm{OPG}}
=\mathbb{E}_n[s(W,\widehat\theta_n)s(W,\widehat\theta_n)'],
$$

and

$$
\widehat{\mathcal I}_{\mathrm{H}}
=-\mathbb{E}_n[\nabla_{\theta\theta'}^2\log f(W;\widehat\theta_n)].
$$

These are the outer product of gradients and negative Hessian estimators. Under correct specification and regularity they converge to the same information matrix, although they need not be numerically equal in finite samples.

Under misspecification, the expected Hessian and score variance need not agree; the sandwich form is then essential for the pseudo true population minimizer, subject to the corresponding conditions. Replacing both components by a single information estimator would discard that distinction.

> [!IMPORTANT]
> **Board work 3 — Least squares sandwich variance**
>
> Derive the least squares sandwich variance. Identify exactly where conditional homoskedasticity would simplify the formula and why it is not needed for consistency of the best linear predictor target.

### Checkpoint 5

1. Is $\widehat V_n$ the variance of $\widehat\theta_n$ or of $\sqrt n(\widehat\theta_n-\theta_0)$?
2. Why do both inverse derivative factors appear when $\Psi_\theta(\theta_0)$ is not symmetric?
3. What equality makes the regular MLE sandwich collapse to inverse information?
4. What target does least squares estimate when the conditional mean is nonlinear?
5. What envelope verifies the uniform LLN for squared loss on a compact coefficient set, and which separate condition identifies the population minimizer?

<a id="w9-stop-6"></a>

## 6. Wald, score, and likelihood ratio inference

### 6.1 Smooth restrictions and Wald inference

Consider

$$
H_0:r(\theta_0)=0,
$$

where $r:\mathbb R^k\to\mathbb R^q$ is differentiable and its Jacobian $R(\theta_0)$ has full row rank. By the delta method,

$$
\sqrt n\mkern3mu r(\widehat\theta_n)
\rightsquigarrow
\mathcal N\mkern-3mu\left(0,R(\theta_0)VR(\theta_0)'
\right)
$$

under the null. The Wald statistic is

$$
W_n
=n\mkern3mu r(\widehat\theta_n)'
\left[
R(\widehat\theta_n)\widehat V_nR(\widehat\theta_n)'
\right]^{-1}
r(\widehat\theta_n),
$$

and

$$
W_n\rightsquigarrow\chi_q^2.
$$

Wald inference measures how far the unrestricted estimator is from satisfying the restriction.

### 6.2 Score type moment tests and likelihood ratio tests

For a point null $H_0:\theta=\theta_0$, suppose a $k$ dimensional estimating equation satisfies

$$
\sqrt n\mkern3mu\Psi_n(\theta_0)
\rightsquigarrow
\mathcal N(0,B),
$$

with $B$ positive definite and $\widehat B_n\to_{\mathbb{P}}B$. Then the score type moment statistic

$$
LM_n^Z
=n\mkern3mu\Psi_n(\theta_0)'
\widehat B_n^{-1}
\Psi_n(\theta_0)
\rightsquigarrow\chi_k^2.
$$

When $\Psi_n$ is the gradient of an M-criterion, this is the corresponding gradient or Lagrange multiplier test. Nuisance parameters and lower dimensional restrictions require evaluating the restricted problem and projecting onto the $q$ tested directions; the resulting degrees of freedom are $q$, not automatically the full parameter dimension.

In a regular likelihood model, let $\widehat\theta_n$ maximize the unrestricted likelihood and $\widetilde\theta_n$ maximize it subject to $r(\theta)=0$.

The **likelihood ratio statistic** is

$$
LR_n
=2\lbrace\ell_n(\widehat\theta_n)-\ell_n(\widetilde\theta_n)\rbrace.
$$

It measures the loss in maximized fit caused by imposing the null.

The **score**, or **Lagrange multiplier**, test evaluates the likelihood slope at the restricted estimator and forms the appropriate information standardized quadratic form. For a point null $H_0:\theta=\theta_0$, it is

$$
LM_n
=s_n(\theta_0)'
\mathcal I_n(\theta_0)^{-1}
s_n(\theta_0).
$$

Under regularity and a $q$ dimensional restriction,

$$
W_n\rightsquigarrow\chi_q^2,
\qquad
LM_n\rightsquigarrow\chi_q^2,
\qquad
LR_n\rightsquigarrow\chi_q^2.
$$

The three tests are first order asymptotically equivalent under the null, but they are not numerically equal in finite samples and may behave differently under reparameterization or weak regularity.

### 6.3 Geometric comparison

- **Wald:** start at the unrestricted estimate and measure distance to the null.
- **Score/LM:** stay on the null and measure the local slope toward the unrestricted model.
- **Likelihood ratio:** compare the optimized criterion values under the unrestricted and restricted models.

Their equivalence comes from the same local quadratic approximation to the log likelihood. Boundary parameters, unidentified nuisance parameters, singular information, or nonsmooth likelihoods can change the limiting distribution.

### 6.4 Normal mean: exact coincidence

Let $X_1,\ldots,X_n$ be i.i.d. $\mathcal N(\mu,\sigma^2)$ with known $\sigma$, and test $H_0:\mu=\mu_0$. The unrestricted MLE is $\widehat\mu=\overline X_n$, while the restricted MLE is $\widetilde\mu=\mu_0$. Direct calculation gives

$$
W_n
=LM_n
=LR_n
=\frac{n(\overline X_n-\mu_0)^2}{\sigma^2}.
$$

Under the null this statistic has an exact $\chi_1^2$ distribution. The exact equality of the three statistics comes from the quadratic Gaussian log likelihood; in a general regular model, only their first order asymptotic equivalence is guaranteed.

A raw difference between unrestricted and restricted values of an arbitrary M-estimation criterion is not automatically chi squared. Its scaling depends on the relationship between criterion curvature and gradient variance. The likelihood ratio statistic obtains its standard scaling from the likelihood structure and information equality. General criterion difference tests require their own normalization and regularity theorem.

> [!TIP]
> **AI interaction 4 — Three tests, one quadratic approximation**
>
> Copy the prompt below into an AI interface and audit the response.

```text
Compare Wald, score/LM, and likelihood ratio tests for a smooth q dimensional
restriction in a regular likelihood model. State where each statistic is
evaluated, explain their common chi squared limit using a local quadratic
approximation, and give two settings in which that equivalence can fail.
Do not claim that the three statistics are equal in finite samples.
```

### Checkpoint 6

1. Which tests require the unrestricted estimator, the restricted estimator, or both?
2. Why does the Wald test depend visibly on the parameterization in finite samples?
3. What determines the chi squared degrees of freedom?
4. Why might a boundary null produce a nonstandard limiting distribution?
5. Why do Wald, score, and LR coincide exactly in the known variance Normal mean example?
6. Why is a generic difference in M-estimation criterion values not automatically chi squared?

<a id="w9-stop-7"></a>

## 7. Synthesis and bridge to applications

M- and Z-estimation assemble the modeling, decision, and asymptotic tools from Weeks 1--6 into a general framework. A complete application should identify the statistical model and target, state the relevant assumptions, and distinguish an exact conclusion from an asymptotic one.

Use these questions to connect the chapters:

1. **Bernoulli model:** Trace the sample mean from sufficient and complete statistic, through UMVU and Cramer--Rao optimality, to exact and asymptotic testing.
2. **Uniform model:** Explain how parameter dependent support affects the MLE, sufficiency, regular score identities, and estimator comparison.
3. **Normal mean:** Compare exact finite sample inference with CLT based inference and identify when Student $t$ is exact.
4. **Best linear prediction:** Begin with conditional expectation under squared loss, restrict the predictor class to linear functions, and derive the Z-estimator and its sandwich variance.
5. **Theorem audit:** For a proposed estimator, separate identification, consistency, rate, asymptotic distribution, standard error consistency, and confidence set validity.

### Checkpoint 7

Choose one of the four models above and draw a dependency graph in which every arrow is labeled by a theorem and its assumptions. Identify one arrow that fails after removing a substantive assumption.

## 8. Recap

- M-estimators optimize sample criteria; Z-estimators solve sample equations.
- Representation determines which theorem assumptions must be checked.
- Method of moments is naturally a Z-procedure and can also be written as a weighted extremum problem; overidentification makes the weight consequential.
- Population identification and sample optimization are distinct.
- Uniform convergence plus separation converts population identification into consistency.
- Compactness, continuity, and an integrable envelope provide one practical route to the required uniform LLN.
- Consistency is global; smooth M- and Z-estimator asymptotic normality follows from a local mean value argument, a uniform LLN for the derivative or Hessian, and a CLT.
- The asymptotic linear representation exposes the first order contribution of each observation.
- Sandwich variance remains valid beyond information equality special cases when its assumptions hold.
- Wald, score, and likelihood ratio tests use different sample calculations but share a regular local quadratic limit.
- A generic criterion difference is not automatically a likelihood ratio statistic and need not have a chi squared limit without problem specific normalization.
- Boundaries, nonsmoothness, weak identification, and singular derivatives are substantive failure modes rather than technical footnotes.

## Notation introduced this week

- $\mathbb{E}_n$: empirical average.
- $Q_n(\theta)$ and $Q(\theta)$: sample and population M-criteria.
- $\Psi_n(\theta)$ and $\Psi(\theta)$: sample and population estimating equations.
- $\psi(W,\theta)$: estimating function.
- $\Psi_\theta(\theta)=\partial\Psi(\theta)/\partial\theta'$: derivative matrix of the population estimating equation.
- $\Psi_{\theta,n}(\theta)=\partial\Psi_n(\theta)/\partial\theta'$: sample derivative matrix.
- $B$: variance of the estimating function.
- $V=[\Psi_\theta(\theta_0)]^{-1}B[\Psi_\theta(\theta_0)]^{-\prime}$: sandwich asymptotic variance.
- $H$ and $\Sigma$: M-criterion Hessian and gradient variance.
- $\Omega_n$: weighting matrix for a quadratic moment criterion.
- $\varphi(W)$: influence function in the asymptotic linear representation.
- $\widehat{\mathcal I}_{\mathrm{OPG}}$ and $\widehat{\mathcal I}_{\mathrm{H}}$: outer product and negative Hessian information estimators.
- $r(\theta)$ and $R(\theta)$: restriction and its Jacobian.
- $W_n$, $LM_n^Z$, $LM_n$, and $LR_n$: Wald, moment score, likelihood score, and likelihood ratio statistics.

## References

- Whitney Newey and Daniel McFadden, “Large Sample Estimation and Hypothesis Testing,” in *Handbook of Econometrics*, vol. 4, 1994.
- Aad van der Vaart, *Asymptotic Statistics*, Cambridge University Press, 1998.
- Jun Shao, *Mathematical Statistics*, Springer, 2003.
- Robert W. Keener, *Theoretical Statistics: Topics for a Core Course*, Springer, 2010.
