# ORF 524: Statistical Theory and Methods

## Week 2: Models, identification, likelihood, and sufficiency

**Status:** Ready for Fall 2026  
**Last updated:** August 30, 2026  
**Meeting:** Wednesday, September 9, 2026

<!-- Legacy sources: 11 ORF 524/lectures/ORF524-Fall2025-Lecture-Week-02.tex and
ORF524-Fall2025-Lecture-Week-02_Lecture.pdf; 09 ORF 524/lectures/
ORF524-Fall2024-Lecture-Week-02_Notes.pdf. -->

## Central question

Before trying to estimate a parameter, how do we determine whether the data can distinguish it and which features of the sample contain information about it?

## Learning goals

By the end of the week, you should be able to:

1. formulate a statistical experiment and distinguish a model, a parameterization, and a target;
2. diagnose whether a parameter or target is identified;
3. construct and interpret a likelihood, including when the support depends on the parameter;
4. explain how population likelihood, Kullback--Leibler divergence, and identification are related;
5. establish sufficiency using conditional distributions or factorization; and
6. distinguish sufficiency from minimal sufficiency, low dimensionality, and usefulness for estimation.

## In-class route

This is one continuous route through the week, not a division into class meetings. Proceed in order, stop wherever the discussion ends, and resume at the next unfinished stop. The linked sections contain the full student facing exposition; the route identifies the concepts and activities to foreground live.

| Stop | Live focus | Mode |
|---|---|---|
| **W2.1** | [Experiments, model size, and random sampling](#w2-stop-1) | Discuss + Checkpoint 1 |
| **W2.2** | [Identification, observational equivalence, and normalization](#w2-stop-2) | Discuss + AI dialogue + Checkpoint 2 |
| **W2.3** | [Sample likelihood, population likelihood, and KL divergence](#w2-stop-3) | Discuss + board work |
| **W2.4** | [The Uniform boundary MLE](#w2-stop-4) | AI audit + Checkpoint 3 |
| **W2.5** | [Statistics and distribution family toolkit](#w2-stop-5) | Brief orientation + selected examples + Checkpoint 4 |
| **W2.6** | [Sufficiency from the conditional distribution definition](#w2-stop-6) | Discuss + Bernoulli example |
| **W2.7** | [Factorization in Bernoulli, Uniform, and exponential families](#w2-stop-7) | Board work |
| **W2.8** | [Minimal sufficiency and false shortcuts](#w2-stop-8) | Discuss + Checkpoint 5; AI dialogue if time |

Section 4 is a reference toolkit rather than a separate extended lecture. Stop W2.5 orients students to the identities and family forms that later arguments will call on; select examples according to need instead of reading every formula. Score, Fisher information, Cramér--Rao bounds, ancillarity, and completeness are intentionally absent from this route and enter in Week 3.

## How to use this chapter

**Prepare:** Read Section 1 and the reference toolkit in Section 4. Attempt Checkpoint 1 before class. You should be comfortable multiplying densities under independence and manipulating sample means and order statistics.

**In class:** Follow the continuous route above. We will focus on three arcs: identification, likelihood with parameter dependent support, and sufficiency. AI dialogues 2 and 3 are the primary live interactions; dialogue 4 may be used if time permits. Prepared responses supply an alternative for every live prompt if AI is unavailable.

**Review:** Work through the Bernoulli sufficiency calculation and the family recognition exercises. Score, Fisher information, Cramer--Rao bounds, ancillarity, and completeness are deliberately deferred to Week 3, where they will be tied to estimator comparison and optimality.

**Practice:** The matching assignment complements this chapter with ungraded problems for consolidating the material, extending selected ideas, and preparing for cumulative review.

**Prerequisites:** Week 1; joint, marginal, and conditional distributions; densities with respect to a dominating measure; independence; and elementary optimization.

## Full chapter map

1. [Statistical experiments and models](#1-statistical-experiments-and-models)
2. [Identification before estimation](#2-identification-before-estimation)
3. [Likelihood and its boundary cases](#3-likelihood-and-its-boundary-cases)
4. [Reference toolkit: statistics and distribution families](#4-reference-toolkit-statistics-and-distribution-families)
5. [Sufficiency and information reduction](#5-sufficiency-and-information-reduction)
6. [Recap](#6-recap)

<a id="w2-stop-1"></a>

## 1. Statistical experiments and models

### 1.1 The experiment

A statistical experiment consists of an observation space $(\mathcal X,\mathcal A)$ and a collection of probability laws

$$
\mathcal P=\lbrace\mathbb{P}_\theta:\theta\in\Theta\rbrace.
$$

Before observing the data, $X$ is a random element satisfying $X\sim\mathbb{P}_{\theta_0}$. After observation, $x=X(\omega)$ is the realized data. The model $\mathcal P$ is an assumption about possible data generating laws, not something observed directly.

A **parameterization** is the map $\theta\mapsto\mathbb{P}_\theta$. A target may use the full index $\theta$ or only a feature

$$
\beta(\theta)=\psi(\mathbb{P}_\theta).
$$

These objects must remain distinct:

- $X$: the random data;
- $x$: the realized data;
- $\mathcal P$: the collection of possible laws;
- $\theta_0$: an index for the law generating the data; and
- $\beta(\theta_0)$: the feature of that law we want to learn.

### 1.2 Model size

A **parametric model** has a finite dimensional parameter. A **nonparametric model** has an infinite dimensional parameter, such as an unrestricted distribution or regression function. A **semiparametric model** combines a finite dimensional target with an infinite dimensional nuisance component.

For a scalar observation, one may consider nested classes such as

$$
\mathcal P_{\mathrm{Gaussian}}
\subset
\mathcal P_{\mathrm{symmetric}}
\subset
\mathcal P_{\mathrm{all}}.
$$

A smaller model supplies more structure and may permit stronger conclusions. If the restriction is wrong, those conclusions need not describe the data generating law.

### Checkpoint 1

For each setting, identify the model, target, and nuisance component. Then state whether the target is obviously identified or requires further analysis.

1. $X_i\sim\mathcal N(\mu,1)$, with $\mu$ unknown.
2. $X_i\sim\mathcal N(\mu,\sigma^2)$, with interest only in $\mu$.
3. $\mathbb E[Y\mid X=x]=m(x)$, with $m$ unrestricted.
4. $\mathbb E[Y\mid X=x]=\beta x+g(x)$, with $\beta$ the target and $g$ unrestricted.

> [!TIP]
> **Prepare — AI interaction 1 — More assumptions, more information?**
>
> Copy the prompt below into an AI interface and audit the response.

```text
Compare three models for a scalar observation:
(A) all probability laws with finite variance,
(B) all laws symmetric about an unknown center, and
(C) all Gaussian laws.

For each model, describe the parameter, what is restricted, and what remains
unknown. Explain one statistical benefit and one risk of moving from A to C.
Do not equate a model assumption with an empirical fact, and distinguish a
mean, median, and symmetry point.
```

**Audit question:** Does the response say precisely which conclusions require finite moments or a unique center?

### 1.3 Random sampling

The default sampling assumption in the first half of the course is

$$
X_i\overset{\mathrm{iid}}{\sim}\mathbb{P}_{\theta_0},
\qquad
(X_1,\ldots,X_n)\sim\mathbb{P}_{\theta_0}^{\otimes n}.
$$

If densities or probability mass functions exist, independence gives

$$
f_X(x;\theta)=\prod_{i=1}^n f(x_i;\theta).
$$

Without independence, the joint density instead obeys the chain rule

$$
f_X(x;\theta)
=f_{X_1}(x_1;\theta)
\prod_{i=2}^n
f_{X_i\mid X_1,\ldots,X_{i-1}}
(x_i\mid x_1,\ldots,x_{i-1};\theta).
$$

Identical marginal distributions alone do not justify the product likelihood.

<a id="w2-stop-2"></a>

## 2. Identification before estimation

### 2.1 Definition

A parameterization is **identified** if it is injective:

$$
\mathbb{P}_{\theta_1}=\mathbb{P}_{\theta_2}
\quad\Longrightarrow\quad
\theta_1=\theta_2.
$$

The full parameter may fail to be identified while a target remains identified. The target $\beta(\theta)$ is identified if

$$
\mathbb{P}_{\theta_1}=\mathbb{P}_{\theta_2}
\quad\Longrightarrow\quad
\beta(\theta_1)=\beta(\theta_2).
$$

Identification is a population property. Increasing $n$ cannot distinguish two parameter values that generate exactly the same observable law.

### 2.2 An unrestricted location decomposition

Suppose

$$
X=\mu+\varepsilon,
$$

where both $\mu\in\mathbb R$ and the distribution $\mathbb{P}_\varepsilon$ are unknown. For any constant $c$,

$$
(\mu,\mathbb{P}_\varepsilon)
\quad\text{and}\quad
(\mu+c,\mathbb{P}_{\varepsilon-c})
$$

generate the same distribution of $X$. Thus the decomposition is not identified.

If the model additionally requires $\mathbb E[\varepsilon]=0$ and the mean exists, then

$$
\mu=\mathbb E[X],
$$

so the normalization identifies $\mu$.

### 2.3 A hidden identification failure

Consider

$$
\mathbb E[Y\mid X=x]=\beta x+g(x),
$$

with $g$ unrestricted. For any $c\in\mathbb R$,

$$
\beta x+g(x)
=(\beta+c)x+\lbrace g(x)-cx\rbrace.
$$

The observable regression function is unchanged, so $\beta$ is not identified. Calling $\beta$ the target and $g$ a nuisance component does not solve the problem; the model needs a restriction that separates their roles.

> [!TIP]
> **AI interaction 2 — Identification is not consistency**
>
> Copy the prompt below into an AI interface and audit the response.

First decide individually whether the claim in the prompt is correct.

```text
In the model X = mu + epsilon, both mu and the distribution of epsilon are
unknown. A student claims that mu is identified because, with enough data, the
sample mean converges to E[X].

Evaluate the claim. Give two observationally equivalent parameter values,
state a normalization that identifies mu, and distinguish identification from
consistency. Do not assume E[epsilon]=0 unless you label it as an added model
restriction.
```

**Takeaway:** Consistency describes an estimator under an identified target and suitable sampling conditions. It cannot repair observational equivalence in the model.

### Checkpoint 2

1. Can two different identified parameterizations describe the same collection of probability laws?
2. Give an example of a full parameter that is not identified but a target of that parameter that is.
3. What kind of restriction could identify $\beta$ in the decomposition $m(x)=\beta x+g(x)$?

<a id="w2-stop-3"></a>

## 3. Likelihood and its boundary cases

### 3.1 Sample likelihood

Suppose the model is dominated by a common measure and has joint density $f_X(x;\theta)$. For observed data $x$, define

$$
L_n(\theta;x)=f_X(x;\theta),
\qquad
\ell_n(\theta;x)=\log L_n(\theta;x).
$$

The likelihood is a function of $\theta$ with $x$ held fixed. It is not a probability distribution for $\theta$ and, for continuous data, it is not the probability of observing exactly $x$.

A **maximum likelihood estimator** is a measurable selection

$$
\widehat\theta_{\mathrm{ML}}
\in\underset{\theta\in\Theta}{\arg\max}\mkern5mu\ell_n(\theta;X).
$$

The set valued notation allows for nonunique maximizers.

### 3.2 Why population likelihood favors the true law

Let $f_0=f(\cdot;\theta_0)$. Under the conditions needed for the expectations below,

$$
\mathbb E_{\theta_0}[\log f(X;\theta)]
\leq
\mathbb E_{\theta_0}[\log f(X;\theta_0)].
$$

The gap is the Kullback--Leibler divergence

$$
D_{\mathrm{KL}}(\mathbb{P}_{\theta_0}\Vert\mathbb{P}_\theta)
=\mathbb E_{\theta_0}
\left[
\log\frac{f(X;\theta_0)}{f(X;\theta)}
\right]
\geq0.
$$

If the true law is not absolutely continuous with respect to the candidate law, the divergence is $+\infty$ under the extended definition. Equality of population objectives identifies the true distribution; identification of the parameterization is what turns equality of distributions into equality of their parameter values.

**Proof map.** On the set where the density ratio is defined, apply Jensen's inequality to $f(X;\theta)/f(X;\theta_0)$. Its expectation under $\mathbb{P}_{\theta_0}$ is at most one, with a strict loss when the candidate assigns mass differently. Track zero density sets separately rather than silently dividing by zero.

> [!IMPORTANT]
> **Board work 1 — Population likelihood and identification**
>
> Complete the population likelihood proof and identify the role of:
>
> 1. the dominating measure;
> 2. integrability of the log ratio;
> 3. absolute continuity; and
> 4. identification.

<a id="w2-stop-4"></a>

### 3.3 Uniform model: the support matters

Let $X_1,\ldots,X_n$ be i.i.d. $\mathsf{Uniform}(0,\theta)$, with $\theta>0$. The likelihood is

$$
L_n(\theta;x)
=\theta^{-n}
\mathbf 1\lbrace0\leq x_{(1)},\ x_{(n)}\leq\theta\rbrace.
$$

On the event $x_{(n)}>0$, feasible values satisfy $\theta\geq x_{(n)}$, and the factor $\theta^{-n}$ decreases with $\theta$. Therefore

$$
\widehat\theta_{\mathrm{ML}}=X_{(n)}.
$$

The support depends on $\theta$, and the maximum occurs at a boundary. Differentiating the smooth part and setting it equal to zero would miss the estimator entirely. Under every $\theta>0$, the event $X_{(n)}>0$ has probability one. On the null sample $x_1=\cdots=x_n=0$, the likelihood has a supremum as $\theta\downarrow0$ but no maximizer in $(0,\infty)$; an estimator may be defined arbitrarily there without changing its sampling properties.

> [!TIP]
> **AI interaction 3 — A boundary MLE**
>
> Copy the prompt below into an AI interface and audit the response.

```text
Derive the MLE of theta from an iid Uniform(0,theta) sample. Write the complete
likelihood, including its indicator. Explain why the sample maximum appears and
why a score equation does not characterize the maximizer. Then identify the
step that would fail in an attempted proof that every score has mean zero.
Do not discard the parameter dependent support.
```

**Audit question:** Does the response differentiate an expression whose support changes with $\theta$ as if its indicator were constant?

### Checkpoint 3

Classify each claim and repair it if necessary.

1. A likelihood is a probability distribution for the parameter.
2. Every MLE is an interior solution of a first order condition.
3. Maximizing the average log likelihood and the log likelihood gives the same argmax.
4. Population likelihood uniquely identifies $\theta_0$ without an identification assumption.

<a id="w2-stop-5"></a>

## 4. Reference toolkit: statistics and distribution families

This section collects material needed later in the course. It is preparation and reference material, not a separate live class arc.

### 4.1 Finite sample statistics

A **statistic** is a measurable map $T:\mathcal X^n\to\mathbb R^k$ whose definition does not involve an unknown parameter. Under i.i.d. sampling with mean $\mu$ and variance $\sigma^2<\infty$,

$$
\overline X_n=\frac1n\sum_{i=1}^nX_i,
\qquad
S_n^2=\frac1{n-1}\sum_{i=1}^n(X_i-\overline X_n)^2
$$

satisfy

$$
\mathbb E[\overline X_n]=\mu,
\qquad
\mathbb V[\overline X_n]=\frac{\sigma^2}{n},
\qquad
\mathbb E[S_n^2]=\sigma^2.
$$

If the sample is Gaussian, then exactly

$$
\overline X_n\sim\mathcal N\mkern-3mu\left(\mu,\frac{\sigma^2}{n}\right),
\qquad
\frac{(n-1)S_n^2}{\sigma^2}\sim\chi^2_{n-1},
$$

and $\overline X_n$ is independent of $S_n^2$.

For order statistics $X_{(1)}\leq\cdots\leq X_{(n)}$ from an i.i.d. sample with c.d.f. $F$,

$$
\mathbb{P}(X_{(n)}\leq x)=F(x)^n.
$$

### 4.2 Location scale families

Given a standard density $f_0$, a location scale family has density

$$
f(x;\mu,\sigma)
=\frac1\sigma f_0\mkern-3mu\left(\frac{x-\mu}{\sigma}\right),
\qquad \mu\in\mathbb R,\quad \sigma>0.
$$

The location and scale parameters equal the mean and standard deviation only when the standardized variable has mean zero and variance one.

### 4.3 Exponential families

A $k$ parameter exponential family can be written

$$
f(x;\eta)=h(x)\exp\lbrace\eta't(x)-A(\eta)\rbrace.
$$

In a regular exponential family, the support does not depend on $\eta$ and the natural parameter space is open. Differentiating the normalization identity under suitable regularity conditions gives

$$
\mathbb E_\eta[t(X)]=\nabla_\eta A(\eta),
\qquad
\mathbb V_\eta[t(X)]=\nabla_\eta^2A(\eta).
$$

For example, if $X\sim\mathsf{Poisson}(\theta)$, then

$$
\eta=\log\theta,
\qquad
t(x)=x,
\qquad
A(\eta)=e^\eta.
$$

The Uniform model is not a regular one parameter exponential family in $\theta$ because its support depends on the parameter.

### Checkpoint 4

1. Which identities for $\overline X_n$ and $S_n^2$ require independence, identical distributions, or Gaussianity?
2. When do $\mu$ and $\sigma$ in a location scale family equal the mean and standard deviation?
3. Put Bernoulli, Poisson, and $\mathcal N(\mu,1)$ into exponential family form.
4. Explain precisely why $\mathsf{Uniform}(0,\theta)$ is not regular in this sense.

<a id="w2-stop-6"></a>

## 5. Sufficiency and information reduction

### 5.1 Definition

A statistic $T(X)$ is **sufficient** for $\theta$ if a version of the conditional distribution of the full sample $X$ given $T(X)$ does not depend on $\theta$.

Sufficiency is relative to a model. It says that, after $T(X)$ is observed, the remaining variation in the sample contains no additional information about $\theta$ within that model.

### 5.2 Bernoulli example from the definition

Let $X_1,\ldots,X_n$ be i.i.d. $\mathsf{Bernoulli}(\theta)$ and define

$$
T(X)=\sum_{i=1}^nX_i.
$$

For any binary vector $x$ with $\sum_i x_i=t$,

$$
\mathbb{P}_\theta(X=x\mid T=t)
=\frac{\theta^t(1-\theta)^{n-t}}
{\binom nt\theta^t(1-\theta)^{n-t}}
=\binom nt^{-1}.
$$

The conditional probability is independent of $\theta$. Once the number of successes is known, their ordering contains no further information about the Bernoulli success probability.

<a id="w2-stop-7"></a>

### 5.3 Factorization theorem

For a dominated model, the **Neyman--Fisher factorization theorem** states that $T(X)$ is sufficient if and only if there exist nonnegative functions $g$ and $h$ such that

$$
f_X(x;\theta)=g(T(x),\theta)h(x)
$$

for all relevant $x$ and $\theta$. All dependence on $\theta$ must enter through $T(x)$.

For an i.i.d. exponential family,

$$
\prod_{i=1}^n f(X_i;\eta)
=\left\lbrace\prod_{i=1}^n h(X_i)\right\rbrace
\exp\left\lbrace
\eta'\sum_{i=1}^n t(X_i)-nA(\eta)
\right\rbrace,
$$

so $\sum_i t(X_i)$ is sufficient for $\eta$.

This conclusion establishes sufficiency, not completeness. [Week 3](week-03.md#33-bernoulli-and-exponential-family-completeness) adds the distinct full rank exponential family completeness theorem: under standard common support and integrability conditions, the same natural statistic is complete when the image of the natural parameter map contains a nonempty open subset of $\mathbb R^k$. The Bernoulli model illustrates how to verify that condition. Exponential family form by itself does not justify a completeness claim.

> [!IMPORTANT]
> **Board work 2 — Sufficiency in Bernoulli and Uniform models**
>
> Use factorization to establish sufficiency in:
>
> 1. the Bernoulli model, where $T=\sum_iX_i$; and
> 2. the Uniform model, where $T=X_{(n)}$.
>
> Explain why the parameter dependent support is an obstacle for regular likelihood calculus but not for the factorization theorem.

<a id="w2-stop-8"></a>

### 5.4 Minimal sufficiency

A sufficient statistic $T$ is **minimal sufficient** if, for every sufficient statistic $S$, there is a measurable function $r$ such that $T=r(S)$ almost surely. Thus every sufficient statistic retains enough information to reconstruct $T$.

Under standard dominated model conditions, the likelihood ratio criterion says that

$$
T(x)=T(y)
\quad\Longleftrightarrow\quad
\frac{f_X(x;\theta)}{f_X(y;\theta)}
\text{ is independent of }\theta,
$$

for sample points where the ratio is meaningful.

> [!TIP]
> **Optional — AI interaction 4 — Audit a sufficiency claim**
>
> Copy the prompt below into an AI interface and audit the response.

```text
A student argues: "The sample mean is sufficient for every iid model because
it summarizes all observations in one number and is a natural estimator of the
population mean."

Find every error in the argument. Explain why dimensionality reduction and
usefulness for estimation do not imply sufficiency. Give one iid model where
the sample mean is sufficient and one where it is not, justifying each
conclusion by factorization or by the conditional distribution definition.
```

### Checkpoint 5

Classify each statement as true or false and justify the answer.

1. Every measurable one-to-one transformation of a sufficient statistic, with measurable inverse on its range, is sufficient.
2. Every low dimensional statistic is sufficient.
3. Every maximum likelihood estimator is sufficient.
4. Every sufficient statistic is minimal sufficient.
5. Sufficiency may change when the statistical model changes.

## 6. Recap

- A statistical model is a collection of possible observable laws; a parameterization labels them.
- Identification is logically prior to estimation and cannot be repaired by a larger sample.
- Likelihood compares parameter values with the observed data held fixed.
- Population likelihood identifies a distribution through Kullback--Leibler divergence; parameter identification is a separate step.
- Parameter dependent support can create boundary MLEs and invalidate routine differentiation.
- Sufficiency is model relative, lossless data reduction; it is not implied by low dimensionality or usefulness for estimation.
- Exponential family factorization produces a natural sufficient statistic; completeness requires an additional full rank open set condition.
- Completeness and information bounds enter in Week 3, when we ask whether an estimator is optimal.

## Notation introduced this week

- $(\mathcal X,\mathcal A,\mathcal P)$: observation space and statistical model.
- $\theta\mapsto\mathbb{P}_\theta$: parameterization.
- $\beta(\theta)$: target parameter.
- $L_n(\theta;x)$ and $\ell_n(\theta;x)$: likelihood and log likelihood.
- $D_{\mathrm{KL}}(\mathbb{P}\Vert\mathbb Q)$: Kullback--Leibler divergence.
- $X_{(j)}$: $j$-th order statistic.
- $T(X)$: statistic, and in Section 5 a sufficient statistic.
- $\eta$, $t(x)$, and $A(\eta)$: natural parameter, statistic, and log partition function of an exponential family.

## References

- Jun Shao, *Mathematical Statistics*, Springer, 2003.
- Robert W. Keener, *Theoretical Statistics: Topics for a Core Course*, Springer, 2010.
- Erich Lehmann and George Casella, *Theory of Point Estimation*, 2nd ed., Springer, 1998.
- Thomas A. Severini, *Elements of Distribution Theory*, Cambridge University Press, 2005.
- Aad van der Vaart, *Asymptotic Statistics*, Cambridge University Press, 1998.
