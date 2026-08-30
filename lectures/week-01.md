# ORF 524: Statistical Theory and Methods

## Week 1: Statistical questions and expectations

**Status:** Ready for Fall 2026  
**Last updated:** August 30, 2026  
**Meeting:** Wednesday, September 2, 2026

<!-- Legacy sources: 11 ORF 524/lectures/ORF524-Fall2025-Lecture-Week-01.tex and
ORF524-Fall2025-Lecture-Week-01_Lecture.pdf; 09 ORF 524/lectures/
ORF524-Fall2024-Lecture-Week-01_Notes.pdf. -->

## Central question

How do the statistical model, the loss function, and the allowed class of procedures determine what we should estimate or predict?

## Learning goals

By the end of the week, you should be able to:

1. distinguish a statistical model, a parameter, an estimator, and a statistic;
2. interpret expectation and conditional expectation as solutions to prediction problems;
3. select and apply the main inequalities used throughout the course; and
4. explain how the choice of loss function determines the optimal predictor.

## In-class route

This is one continuous route through the week, not a division into class meetings. Proceed in order, stop wherever the discussion ends, and resume at the next unfinished stop. The linked sections contain the full student facing exposition; the route identifies the concepts and activities to foreground live.

| Stop           | Live focus                                                      | Mode                                    |
| -------------- | --------------------------------------------------------------- | --------------------------------------- |
| **W1.1** | [Model, target, estimator, and statistical question](#w1-stop-1) | Discuss + Checkpoint 1                  |
| **W1.2** | [Expectation, moments, and existence conditions](#w1-stop-2)     | Discuss + Checkpoint 2                  |
| **W1.3** | [Loss determines the optimal constant predictor](#w1-stop-3)     | Board work + AI assisted discussion     |
| **W1.4** | [Moment inequalities and their assumptions](#w1-stop-4)          | Board work + Checkpoint 3 + assumption audit |
| **W1.5** | [Conditional expectation as projection](#w1-stop-5)              | Board work + AI dialogue + Checkpoint 4 |

Sections or dialogues labeled **Prepare** or **Review** remain part of the chapter but are not planned live stops. Within Stop W1.4, Markov and Chebyshev are the main derivation; Jensen, Hölder, Cauchy--Schwarz, and Minkowski form a compact toolkit to discuss selectively rather than prove in full.

## How to use this chapter

**Prepare:** Read Sections 1 and 2, attempt Checkpoints 1 and 2, and review the statements of the inequalities in Section 4. Record every moment or measurability condition you see.

**In class:** Follow the continuous route above. We will focus on loss dependent targets, prove Markov and derive Chebyshev, and develop the orthogonal projection argument for conditional expectation. AI interactions 2 and 4 are the primary live interactions. A prepared response supplies an alternative if a live AI system is unavailable.

**Review:** Revisit Checkpoints 3 and 4. AI dialogues 1 and 3 provide additional terminology and assumption audit practice.

**Practice:** The matching assignment complements this chapter with ungraded problems for consolidating the material, extending selected ideas, and preparing for cumulative review.

**Prerequisites:** Basic probability spaces, random variables and distributions, integration and expectation, conditional probability, and elementary convexity. Measure theoretic conditional expectation is useful but not required at the start of the week.

## Full chapter map

1. [From data to statistical questions](#1-from-data-to-statistical-questions)
2. [Expectation and moments](#2-expectation-and-moments)
3. [Expectation as a best constant predictor](#3-expectation-as-a-best-constant-predictor)
4. [Probability and moment inequalities](#4-probability-and-moment-inequalities)
5. [Conditional expectation as a best predictor](#5-conditional-expectation-as-a-best-predictor)
6. [Recap](#6-recap)

<a id="w1-stop-1"></a>

## 1. From data to statistical questions

Let $X$ denote an observable random variable with distribution

$$
X\sim \mathbb{P}_{\theta_0},
\qquad
\mathbb{P}_{\theta_0}\in
\mathcal P=\lbrace\mathbb{P}_\theta:\theta\in\Theta\rbrace.
$$

The set $\mathcal P$ is a **statistical model** and $\theta_0$ indexes the distribution that generated the data. Often we care not about all of $\theta_0$, but about a feature

$$
\beta_0=\beta(\theta_0)=\psi(\mathbb{P}_{\theta_0}).
$$

Given observations $X_1,\ldots,X_n$, an **estimator** of $\beta_0$ is a measurable function of the sample,

$$
\widehat\beta_n=\widehat\beta_n(X_1,\ldots,X_n).
$$

A **statistic** is any measurable function of the sample whose definition does not involve an unknown parameter. An estimator is therefore a statistic used for point estimation. A statistic may instead be used for testing, confidence sets, diagnostics, or data reduction.

The course will repeatedly return to three questions:

- **Point estimation:** Which estimator should we use, and according to what criterion?
- **Hypothesis testing:** How should evidence against a claim be measured and calibrated?
- **Uncertainty quantification:** How should we construct and interpret confidence sets?

The same observation can be embedded in different models, and the model determines which questions are meaningful. For example,

$$
\mathcal P_1=\lbrace\mathcal N(\theta,1):\theta\in\mathbb R\rbrace
\quad\text{and}\quad
\mathcal P_2=\lbrace\mathcal N(\mu,\sigma^2):\mu\in\mathbb R,\ \sigma^2>0\rbrace
$$

both describe a real valued observation, but only the second treats the variance as unknown. Even within one model, different functionals define different targets. Two examples are

$$
\mu(\mathbb{P})=\mathbb E_{\mathbb{P}}[X]
\qquad\text{and}\qquad
F_{\mathbb{P}}(t)=\mathbb E_{\mathbb{P}}[\mathbf 1\lbrace X\leq t\rbrace],\quad t\in\mathbb R.
$$

The first target is a number; the second is an entire function. Neither the model nor the data alone selects which target matters.

Possible properties of an estimator include finite sample unbiasedness,

$$
\mathbb E_\theta[\widehat\beta_n]=\beta(\theta),
$$

and consistency,

$$
\widehat\beta_n\to_{\mathbb{P}}\beta(\theta_0).
$$

Neither property, by itself, determines a unique estimator. Statistical procedures must be compared using an explicit objective or loss.

### Checkpoint 1

Suppose $X_1,\ldots,X_n$ are i.i.d. with unknown mean $\mu$ and variance $\sigma^2$.

1. Give two statistics computed from the sample.
2. Which one would naturally estimate $\mu$?
3. Is an estimator still a random variable after the data have been observed?

> [!TIP]
> **Review — AI interaction 1 — What exactly is unknown?**
>
> Copy the prompt below into an AI interface and audit the response.

First answer the following prompt individually. We will then ask an AI system and evaluate its answer.

```text
In the model X_1,...,X_n iid ~ N(mu,sigma^2), distinguish carefully among:
(i) the model, (ii) the parameter, (iii) the true parameter value,
(iv) an estimator, (v) an estimate, and (vi) a statistic.
Give one example of each and identify any terms that overlap.
Do not discuss asymptotics.
```

**Question for the class:** Which distinctions in the answer are mathematical, and which are merely terminological?

<a id="w1-stop-2"></a>

## 2. Expectation and moments

Let $X$ be a random variable on $(\Omega,\mathcal F,\mathbb{P})$. If $X$ is integrable, its **expected value** is

$$
\mathbb E[X]
=\int_\Omega X(\omega)\mkern3mu d\mathbb{P}(\omega)
=\int_{\mathbb R}x\mkern3mu dF_X(x).
$$

All properties of the integral yield corresponding properties of expectation. In particular, for integrable $X$ and $Y$ and constants $a,b$,

$$
\mathbb E[aX+bY]=a\mathbb E[X]+b\mathbb E[Y].
$$

If $X\leq Y$ almost surely, then $\mathbb E[X]\leq\mathbb E[Y]$ whenever the expectations are well defined.

If $X\in L^2$, its **variance** and **standard deviation** are

$$
\mathbb V[X]
=\mathbb E\mkern-3mu\left[(X-\mathbb E[X])^2\right]
=\mathbb E[X^2]-\mathbb E[X]^2,
\qquad
\mathrm{sd}(X)=\sqrt{\mathbb V[X]}.
$$

For constants $a,b$,

$$
\mathbb V[a+bX]=b^2\mathbb V[X].
$$

More generally, when the expectations exist:

- the raw moment of order $k$ is $\mu_k'=\mathbb E[X^k]$;
- the central moment of order $k$ is $\mu_k=\mathbb E[(X-\mathbb E[X])^k]$.

Moments describe some, but not necessarily all, aspects of a distribution. Existence must never be assumed silently.

### Checkpoint 2

Give an example of a distribution for which $\mathbb E[X]$ does not exist. Which expressions in this section then become undefined?

<a id="w1-stop-3"></a>

## 3. Expectation as a best constant predictor

Let $Y$ be square integrable. If we must predict $Y$ using a constant $a$ and assess error with squared loss, then

$$
\mathbb E[Y]
=\underset{a\in\mathbb R}{\arg\min}\mkern5mu
\mathbb E[(Y-a)^2].
$$

The identity behind this result is

$$
\mathbb E[(Y-a)^2]
=\mathbb V[Y]+(\mathbb E[Y]-a)^2.
$$

Thus $a=\mathbb E[Y]$ is the unique minimizer whenever $Y\in L^2$.

> [!IMPORTANT]
> **Board work 1 — Derive the squared loss decomposition**
>
> Develop the proof live on the iPad or board. Starting from $m=\mathbb E[Y]$, derive the displayed decomposition, identify why the cross term vanishes, state where square integrability is used, and explain why the minimizer is unique. The written derivation is intentionally reserved for class.

The answer changes when the loss changes. Under absolute loss, the set of minimizers is the set of medians of $Y$:

$$
\underset{a\in\mathbb R}{\arg\min}\mkern5mu
\mathbb E[|Y-a|]
=\lbrace\text{medians of }Y\rbrace.
$$

The minimizer need not be unique. This simple observation anticipates a major theme of the course: an “optimal” statistical procedure is meaningful only after the criterion has been specified.

> [!TIP]
> **AI interaction 2 — Why absolute loss selects a median**
>
> Copy the prompt below into an AI interface and audit the response.

Use the following prompt to compare the squared and absolute loss targets without developing a second full proof on the board.

```text
Let Y be an integrable real random variable. Compare the optimal constant predictor
of Y under squared error loss and absolute error loss.

1. State the relevant integrability conditions.
2. Give a proof map—not a complete proof—for why the absolute loss minimizers satisfy
   F_Y(a-) <= 1/2 <= F_Y(a). Explicitly account for atoms in the distribution.
3. Give a discrete distribution for which the absolute loss minimizer is not unique.
4. Explain why changing the loss changes the statistical target.
```

**Audit questions:** Did the answer distinguish a median from a unique median? Did it impose a second moment for the squared loss problem? Did it incorrectly assume differentiability when $Y$ has an atom?

<a id="w1-stop-4"></a>

## 4. Probability and moment inequalities

These inequalities translate information about moments into bounds on probabilities and other moments. They will be used repeatedly in finite sample and asymptotic arguments.

**Prepare:** For each inequality, identify its input assumptions and the form of its output. In class we will prove Markov and derive Chebyshev; the remaining inequalities form a reference toolkit for later weeks.

### 4.1 Markov and Chebyshev

**Markov's inequality.** If $g(X)\geq0$ almost surely and $\mathbb E[g(X)]<\infty$, then for every $\varepsilon>0$,

$$
\mathbb{P}\lbrace g(X)\geq\varepsilon\rbrace
\leq \frac{\mathbb E[g(X)]}{\varepsilon}.
$$

**Chebyshev's inequality.** If $X$ has mean $\mu$ and finite variance $\sigma^2$, then for every $\varepsilon>0$,

$$
\mathbb{P}\lbrace|X-\mu|\geq\varepsilon\rbrace
\leq\frac{\sigma^2}{\varepsilon^2}.
$$

> [!IMPORTANT]
> **Board work 2 — Markov and Chebyshev inequalities**
>
> 1. Prove Markov's inequality by comparing $g(X)$ with $\varepsilon\mathbf 1\lbrace g(X)\geq\varepsilon\rbrace$.
> 2. Derive Chebyshev's inequality as one application of Markov's inequality.
> 3. Determine when Markov's bound can hold with equality.

### 4.2 Jensen

If $g$ is convex and the relevant expectations exist, then

$$
g(\mathbb E[X])\leq\mathbb E[g(X)].
$$

For concave $g$, the inequality is reversed:

$$
\mathbb E[g(X)]\leq g(\mathbb E[X]).
$$

Strict convexity yields a strict inequality unless $X$ is constant almost surely, subject to the usual domain conditions.

### 4.3 $L^p$ norms and related inequalities

For $p\geq1$, define

$$
\Vert X\Vert_p=\left(\mathbb E[|X|^p]\right)^{1/p}.
$$

**Hölder's inequality.** If $p,q>1$, $p^{-1}+q^{-1}=1$, $X\in L^p$, and $Y\in L^q$, then

$$
\mathbb E[|XY|]\leq\Vert X\Vert_p\Vert Y\Vert_q.
$$

**Cauchy–Schwarz inequality.** If $X,Y\in L^2$, taking $p=q=2$ gives

$$
|\mathbb E[XY]|
\leq \mathbb E[|XY|]
\leq \sqrt{\mathbb E[X^2]}\sqrt{\mathbb E[Y^2]}.
$$

Applied to centered variables, it gives the covariance inequality

$$
|\mathrm{Cov}(X,Y)|^2
\leq \mathbb V[X]\mathbb V[Y].
$$

**Minkowski's inequality.** If $p\geq1$ and $X,Y\in L^p$, then

$$
\Vert X+Y\Vert_p\leq\Vert X\Vert_p+\Vert Y\Vert_p.
$$

### Checkpoint 3

For each goal below, select an inequality and state the conditions needed to apply it.

1. Bound $\mathbb{P}(|X-\mathbb E[X]|\geq t)$ using a variance.
2. Bound $|\mathbb E[XY]|$ using second moments.
3. Compare $g(\mathbb E[X])$ and $\mathbb E[g(X)]$.
4. Bound the $L^p$ norm of a sum.

> [!TIP]
> **Review — AI interaction 3 — Find the missing assumption**
>
> Copy the prompt below into an AI interface and audit the response.

```text
For the theorem selected by the instructor from Markov, Jensen,
Cauchy--Schwarz, and Minkowski, construct a short plausible looking proof that
omits exactly one essential assumption without announcing which assumption was
omitted. Keep the conclusion mathematically standard. We will identify and
repair the gap.
```

The AI output is an object to audit, not a source of the theorem. A prepared flawed proof will be available if the generated example does not isolate a single clear error.

<a id="w1-stop-5"></a>

## 5. Conditional expectation as a best predictor

Let $Y\in L^2$ and let $X$ be another random variable. Among all square integrable measurable functions of $X$, the conditional expectation minimizes mean squared error:

$$
\mathbb E[Y\mid X]
=\underset{g:\mkern3mu g(X)\in L^2}{\arg\min}\mkern5mu
\mathbb E[(Y-g(X))^2].
$$

The minimizer is unique up to almost sure equality. Its minimum risk is

$$
\mathbb E\mkern-3mu\left[(Y-\mathbb E[Y\mid X])^2\right]
=\mathbb E[\mathbb V(Y\mid X)].
$$

The key orthogonality property is

$$
\mathbb E\mkern-3mu\left[(Y-\mathbb E[Y\mid X])h(X)\right]=0
$$

for every square integrable $h(X)$. Consequently,

$$
\mathbb E[(Y-g(X))^2]
=\mathbb E[(Y-\mathbb E[Y\mid X])^2]
+\mathbb E[(\mathbb E[Y\mid X]-g(X))^2].
$$

This is a Pythagorean identity in $L^2$. Conditional expectation is the orthogonal projection of $Y$ onto the space of square integrable functions of $X$.

**Proof map.** Set $m(X)=\mathbb E[Y\mid X]$ and write

$$
Y-g(X)=\lbrace Y-m(X)\rbrace+\lbrace m(X)-g(X)\rbrace.
$$

After expanding the square, conditional expectation makes the cross term zero. The remaining excess risk is $\mathbb E[(m(X)-g(X))^2]$, which is nonnegative and vanishes only up to almost sure equality.

The unrestricted square integrable class and the affine class can be written as

$$
\mathcal G_2=\lbrace g:\mathbb E[g(X)^2]<\infty\rbrace,
\qquad
\mathcal G_L=\lbrace g\in\mathcal G_2:g(x)=a+bx,\ a,b\in\mathbb R\rbrace.
$$

Over $\mathcal G_2$, the optimizer is $\mathbb E[Y\mid X]$. Over the smaller class $\mathcal G_L$, the optimizer is the best linear predictor and need not equal the conditional mean. Statistical and machine learning methods often differ precisely in the function class over which they optimize.

> [!IMPORTANT]
> **Board work 3 — Conditional expectation Pythagorean identity**
>
> Derive the Pythagorean identity. Then show that the risk of the optimal predictor is $\mathbb E[\mathbb V(Y\mid X)]$.

### Checkpoint 4

Let $Y=X^2+\varepsilon$, where $\mathbb E[\varepsilon\mid X]=0$.

1. What is the best unrestricted predictor of $Y$ given $X$ under squared loss?
2. Would the best linear predictor generally be the same?
3. What additional condition would make $X^2$ an affine function of $X$ almost surely?

> [!TIP]
> **AI interaction 4 — Projection, not magic**
>
> Copy the prompt below into an AI interface and audit the response.

```text
Explain why E[Y|X] is the best mean squared error predictor of Y using two
different arguments:

(i) an orthogonal projection argument in L2, and
(ii) when a regular conditional distribution is available, a pointwise
argument that conditions on X=x.

State all moment and measurability requirements, and do not treat E[Y|X=x] as
canonically defined at every x. Then explain what changes if the predictor is
restricted to be affine in X. Keep the three optimization problems clearly
separated.
```

**Audit question:** Does the response confuse the conditional expectation with the best linear predictor?

## 6. Recap

- A statistical model lists the probability laws under consideration; a parameter indexes or extracts features of those laws.
- An estimator is a statistic used to learn an unknown target from data.
- Expectation and conditional expectation solve squared loss prediction problems.
- Optimality depends on the loss and on the allowed class of predictors.
- Markov, Chebyshev, Jensen, Hölder, Cauchy–Schwarz, and Minkowski convert assumptions into useful bounds.
- Every theorem comes with conditions. Identifying those conditions is part of using the theorem.

## Notation introduced this week

- $\mathcal P=\lbrace\mathbb{P}_\theta:\theta\in\Theta\rbrace$: statistical model.
- $\beta(\theta)=\psi(\mathbb{P}_\theta)$: target parameter or functional.
- $\widehat\beta_n$: estimator based on $X_1,\ldots,X_n$.
- $\mathbb E$, $\mathbb V$, and $\mathrm{Cov}$: expectation, variance, and covariance.
- $\Vert X\Vert_p$: $L^p$ norm.
- $\mathbb E[Y\mid X]$: conditional expectation, defined up to almost sure equality.

## References

- Jun Shao, *Mathematical Statistics*, Springer, 2003.
- Robert W. Keener, *Theoretical Statistics: Topics for a Core Course*, Springer, 2010.
- Aad van der Vaart, *Asymptotic Statistics*, Cambridge University Press, 1998.
