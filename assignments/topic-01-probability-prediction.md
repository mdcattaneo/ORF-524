# ORF 524 Practice Module 1: Probability, Loss, and Prediction

**Status:** Ready for Fall 2026  
**Last updated:** August 30, 2026  
**Primary chapter:** [Week 1: Statistical questions and expectations](../lectures/week-01.md)

## Purpose

This ungraded module develops the main idea of Week 1: a statistical target is determined jointly by the probability model, the loss function, and the class of procedures under consideration. Problems 1--4 form the core bank. Problems 5--6 are an optional bank, not additional required work; they review foundations and expose an important failure caused by missing integrability.

<a id="problem-map"></a>

## Problem map

| Problem | Bank | Main task |
|---|---|---|
| [1. Loss determines the target](#problem-1) | Core | Compare mean, median, and quantile targets under different loss functions. |
| [2. Inequalities and their assumptions](#problem-2) | Core | Prove moment inequalities and audit their assumptions and equality cases. |
| [3. Conditional expectation as projection](#problem-3) | Core | Use conditional expectation as an optimal predictor and derive variance decompositions. |
| [4. The predictor class matters](#problem-4) | Core | Compare unrestricted and affine prediction under squared loss. |
| [5. An integrability boundary for iterated expectations](#problem-5) | Further | Diagnose why iterated expectation can fail without integrability. |
| [6. Continuity of probability](#problem-6) | Further | Prove continuity of probability and countable subadditivity. |

## How to work with this module

For each problem, first write down the assumptions and the object being optimized or bounded. Give a mathematical argument, not only a named theorem or a numerical answer. When a minimizer is not unique, report the entire set of minimizers.

AI may be used on these ungraded problems under the syllabus policy and the staged protocol in [`AGENTS.md`](AGENTS.md). Begin with a genuine attempt, request the least revealing useful hint, and verify every generated claim line by line. The closed book midterms remain fully unaided.

For released exam practice, use the [question level first half guide](../exams/README.md#question-level-guide-for-the-first-half). Confirm that the listed material has been covered, attempt the selected question before opening its solution, and treat historical exam instructions as superseded by the current syllabus.

## Core practice

<a id="problem-1"></a>

### Problem 1. Loss determines the target

**Chapter connection:** [Section 3: Expectation as a best constant predictor](../lectures/week-01.md#3-expectation-as-a-best-constant-predictor)

Let $Y$ be a real valued random variable. For $\tau\in(0,1)$, define the check loss

$$
\rho_\tau(u)=u\bigl(\tau-\mathbf 1\lbrace u<0\rbrace\bigr) =
\begin{cases}
\tau u, & u\geq 0,\\
(1-\tau)(-u), & u<0.
\end{cases}
$$

1. Suppose $Y\in L^2$. Prove that

$$
\mathbb E[(Y-a)^2]
=\mathbb V[Y]+(\mathbb E[Y]-a)^2,
$$

   and identify the unique optimal constant predictor under squared loss.

2. Suppose $Y\in L^1$. Show that the set of optimal constant predictors under absolute loss is

$$
\left\lbrace a\in\mathbb R:
F_Y(a^-)\leq\frac12\leq F_Y(a)\right\rbrace.
$$

3. Still assuming $Y\in L^1$, show more generally that

$$
\underset{a\in\mathbb R}{\arg\min}\mkern5mu
\mathbb E[\rho_\tau(Y-a)]
=\lbrace a\in\mathbb R:F_Y(a^-)\leq\tau\leq F_Y(a)\rbrace.
$$

   You may establish and use, for $b>a$,

$$
\mathbb E[\rho_\tau(Y-b)]-
\mathbb E[\rho_\tau(Y-a)]
=\int_a^b\lbrace F_Y(t)-\tau\rbrace\mkern3mu dt.
$$

4. Let $\mathbb{P}(Y=0)=\mathbb{P}(Y=2)=1/2$. Find the optimal constant predictor under squared loss, absolute loss, $\rho_{1/4}$ loss, and $\rho_{3/4}$ loss. Explain in one sentence what this example says about the phrase “the optimal predictor.”

<!-- Source lineage: Fall 2024 pset2, “Median,” and Fall 2025 Problem Set Weeks 1--2, “Median”; substantially expanded to unify squared, absolute, and check loss. -->

[Back to the problem map](#problem-map)

<a id="problem-2"></a>

### Problem 2. Inequalities and their assumptions

**Chapter connection:** [Section 4: Probability and moment inequalities](../lectures/week-01.md#4-probability-and-moment-inequalities)

Let $X$ be a real valued random variable.

1. If $1\leq p<q$ and $\mathbb E[|X|^q]<\infty$, prove that

$$
\Vert X\Vert_p\leq\Vert X\Vert_q.
$$

   State precisely where you use the fact that $\mathbb{P}$ is a probability measure.

2. Suppose $X\in L^2$, let $\mu=\mathbb E[X]$ and $\sigma^2=\mathbb V[X]$, and let $m$ be any median of $X$. Prove the mean--median inequality

$$
|\mu-m|\leq\sigma.
$$

   Identify each inequality or optimality property used in your argument.

3. Let $W\geq0$ almost surely with $\mathbb E[W]<\infty$. Prove Markov's inequality

$$
\mathbb{P}(W\geq t)\leq\frac{\mathbb E[W]}{t},
\qquad t>0,
$$

   and characterize all distributions of $W$ for which equality holds at a fixed $t$.

4. Let $X_1,\ldots,X_n$ be independent with common mean $\mu$ and common finite variance $\sigma^2$. Use Markov's inequality to derive

$$
\mathbb{P}\mkern-3mu\left(\left|\frac1n\sum_{i=1}^nX_i-\mu\right|\geq\varepsilon\right)
\leq\frac{\sigma^2}{n\varepsilon^2},
\qquad \varepsilon>0.
$$

   Which stated assumption can be weakened if the variables are pairwise uncorrelated rather than independent?

<!-- Source lineage: Fall 2024 pset1 and Fall 2025 Problem Set Weeks 1--2, “Inequalities”; curated and extended with an equality case and an assumption audit. -->

[Back to the problem map](#problem-map)

<a id="problem-3"></a>

### Problem 3. Conditional expectation as projection

**Chapter connection:** [Section 5: Conditional expectation as a best predictor](../lectures/week-01.md#5-conditional-expectation-as-a-best-predictor)

Let $X$ and $Y$ be random variables with $Y\in L^2$, and write

$$
m(X)=\mathbb E[Y\mid X].
$$

1. Prove that

$$
\mathbb V(Y\mid X)
=\mathbb E[Y^2\mid X]-m(X)^2
$$

   almost surely, and then prove the law of total variance

$$
\mathbb V[Y]
=\mathbb E[\mathbb V(Y\mid X)]+\mathbb V[m(X)].
$$

2. For every measurable $g$ such that $g(X)\in L^2$, prove the Pythagorean identity

$$
\mathbb E[(Y-g(X))^2]
=\mathbb E[(Y-m(X))^2]
+\mathbb E[(m(X)-g(X))^2].
$$

   Conclude that $m(X)$ is the unique mean squared error minimizer up to almost sure equality.

3. Suppose now that

$$
\mathbb E[Y\mid X]=X,
\qquad
\mathbb E[Y^2]=\mathbb E[X^2]<\infty.
$$

   Prove that $X=Y$ almost surely.

4. Show that neither condition in part 3 can simply be dropped. Give one square integrable example satisfying $\mathbb E[Y\mid X]=X$ but not the equal second moment condition, and another satisfying the equal second moment condition but not $\mathbb E[Y\mid X]=X$. In each case, ensure that $\mathbb{P}(X\neq Y)>0$.

<!-- Source lineage: Fall 2024 pset2 and Fall 2025 Problem Set Weeks 1--2, “Conditional Variance, Properties” and “Perfect Prediction”; reorganized around the projection identity. -->

[Back to the problem map](#problem-map)

<a id="problem-4"></a>

### Problem 4. The predictor class matters

**Chapter connection:** [Section 5: Conditional expectation as a best predictor](../lectures/week-01.md#5-conditional-expectation-as-a-best-predictor)

Let $X\sim\mathsf{Uniform}(-1,1)$ and

$$
Y=X^2+\varepsilon,
$$

where $\varepsilon\in L^2$ satisfies

$$
\mathbb E[\varepsilon\mid X]=0,
\qquad
\mathbb E[\varepsilon^2\mid X]=\sigma^2
$$

almost surely for some $\sigma^2\geq0$.

1. Find the best mean squared error predictor of $Y$ among all square integrable measurable functions of $X$. Compute its mean squared error.

2. Find the best mean squared error predictor of $Y$ in the affine class

$$
\mathcal G_{\mathrm{aff}}=\lbrace a+bX:a,b\in\mathbb R\rbrace.
$$

   Derive the normal equations rather than quoting a regression formula.

3. Compute the minimum mean squared error over $\mathcal G_{\mathrm{aff}}$ and the excess risk caused by the affine restriction.

4. Both parts use squared loss. Explain why their answers differ, and identify exactly which feature of the optimization problem changed.

<!-- Source lineage: Fall 2024 pset2 and Fall 2025 Problem Set Weeks 1--2, “Prediction I” and “Prediction II”; new example designed to make unrestricted and affine prediction differ. -->

[Back to the problem map](#problem-map)

## Further practice

<a id="problem-5"></a>

### Problem 5. An integrability boundary for iterated expectations

**Chapter connections:** [Section 2: Expectation and moments](../lectures/week-01.md#2-expectation-and-moments) and [Section 5: Conditional expectation as a best predictor](../lectures/week-01.md#5-conditional-expectation-as-a-best-predictor)

Let $U,V\stackrel{\mathrm{iid}}{\sim}\mathcal N(0,1)$ and define $Z=U/V$.

1. By transforming $(U,V)$ to $(Z,V)$, show that $Z$ has density

$$
f_Z(z)=\frac{1}{\pi(1+z^2)},
\qquad z\in\mathbb R.
$$

2. Show that $\mathbb E[|Z|]=\infty$ and explain why $\mathbb E[Z]$ is not zero but undefined.

3. For every $v\neq0$, the conditional distribution of $Z$ given $V=v$ has mean zero. Explain why this pointwise calculation does not justify

$$
\mathbb E[Z]=\mathbb E[\mathbb E[Z\mid V]].
$$

   State the missing assumption.

<!-- Source lineage: Fall 2024 pset2 and Fall 2025 Problem Set Weeks 1--2, “LIE”; corrected to distinguish a pointwise conditional-distribution calculation from an L1 conditional expectation. -->

[Back to the problem map](#problem-map)

<a id="problem-6"></a>

### Problem 6. Continuity of probability

**Chapter connection:** This foundational extension supports the probability toolkit used in [Section 4](../lectures/week-01.md#4-probability-and-moment-inequalities), although continuity of probability is developed here rather than in the chapter.

Let $(A_n)_{n\geq1}$ be events in $(\Omega,\mathcal F,\mathbb{P})$.

1. If $A_n\uparrow A$, prove from countable additivity that

$$
\mathbb{P}(A_n)\longrightarrow\mathbb{P}(A).
$$

2. If $A_n\downarrow A$, prove the same conclusion. Identify where finiteness of $\mathbb{P}(\Omega)$ enters the argument.

3. Use continuity from below to prove that

$$
\mathbb{P}\mkern-3mu\left(\bigcup_{n=1}^{\infty}A_n\right)
\leq\sum_{n=1}^{\infty}\mathbb{P}(A_n).
$$

<!-- Source lineage: Fall 2024 pset1, “Probability continuity properties”; shortened to the foundational results used repeatedly later in the course. -->

[Back to the problem map](#problem-map)

## Completion check

After finishing the core, you should be able to answer all four questions without calculation:

1. Why is “best predictor” incomplete unless a loss and predictor class are specified?

2. What moment condition is needed for each inequality you used?

3. Why does the cross term vanish in the conditional expectation projection identity?

4. What can fail when an expectation is manipulated without first checking integrability?
