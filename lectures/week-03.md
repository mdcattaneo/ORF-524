# ORF 524: Statistical Theory and Methods

## Week 3: Decision theory and point estimation

**Status:** Ready for Fall 2026  
**Last updated:** August 30, 2026  
**Meetings:** Monday, September 14, and Wednesday, September 16, 2026

<!-- Legacy sources: 11 ORF 524/lectures/ORF524-Fall2025-Lecture-Week-02.tex and
ORF524-Fall2025-Lecture-Week-03.tex, with the corresponding _Lecture.pdf board notes;
09 ORF 524/lectures/ORF524-Fall2024-Lecture-Week-02_Notes.pdf and Week-03_Notes.pdf. -->

## Central question

What does it mean for one estimator to be better than another, and when can sufficiency or information bounds identify an optimal estimator?

## Learning goals

By the end of the week, you should be able to:

1. formulate a statistical decision problem and compare rules using dominance, admissibility, Bayes risk, and minimax risk;
2. decompose mean squared error and distinguish unrestricted from restricted optimality;
3. use Rao--Blackwellization, completeness, and Lehmann--Scheffe to construct UMVU estimators, and state the structural roles of Bahadur and Basu;
4. derive score and Fisher information identities under explicit regularity conditions;
5. prove and apply scalar and matrix Cramer--Rao lower bounds, including nuisance parameters, and use bound attainment as a second route to a UMVU conclusion; and
6. distinguish estimator construction by moments or likelihood from an optimality argument.

## In-class route

This is one continuous route through the week, not a division into class meetings. Proceed in order, stop wherever the discussion ends, and resume at the next unfinished stop. The linked sections contain the full student facing exposition; the route identifies the concepts and activities to foreground live.

| Stop | Live focus | Mode |
|---|---|---|
| **W3.1** | [Decision rules, risk, admissibility, Bayes, and minimax](#w3-stop-1) | Discuss + Checkpoint 1 + AI dialogue |
| **W3.2** | [Bias, variance, MSE, and restricted optimality](#w3-stop-2) | Discuss + Checkpoint 2 |
| **W3.3** | [Rao--Blackwell, completeness, UMVU, and ancillarity](#w3-stop-3) | Discuss + Rao--Blackwell proof + complete sufficiency construction + AI implication audit + Checkpoint 3 |
| **W3.4** | [Score and Fisher information identities](#w3-stop-4) | Board work + regularity audit |
| **W3.5** | [Cramer--Rao bounds and the second UMVU route](#w3-stop-5) | Board work + Bernoulli dual certification + AI differentiation audit + Checkpoint 4 |
| **W3.6** | [Moment identification and construction by moments and likelihood](#w3-stop-6) | Discuss + moment identification examples + Uniform board work + AI construction audit + Checkpoint 5 |

The score identity, information equality, and Cramer--Rao argument are separate route stops so their regularity conditions and moving support failure mode remain visible. The Bernoulli sample mean then reunites Stops W3.3 and W3.5 by certifying the same UMVU estimator through complete sufficiency and through attainment of the Cramer--Rao bound. Bahadur and Basu provide the structural close to Stop W3.3; their longer consequences may be abbreviated if needed.

## How to use this chapter

**Prepare:** Review sufficiency and exponential families from Week 2. Read Sections 1 and 2 and attempt Checkpoints 1 and 2.

**In class:** Follow the continuous route above through decision theoretic comparisons, Rao--Blackwell--Lehmann--Scheffe optimality, score and information identities, Cramer--Rao bounds, and estimator construction. Prepared responses supply an alternative for every live prompt if AI is unavailable.

**Review:** Complete the Bernoulli and Uniform calculations and audit every use of differentiation under the integral. The AI prompts are designed to test implications among concepts, not to supply the definitions.

**Practice:** The matching assignment complements this chapter with ungraded problems for consolidating the material, extending selected ideas, and preparing for cumulative review.

**Prerequisites:** Weeks 1--2; conditional expectation; variance decomposition; likelihood; sufficiency; elementary differentiation and Cauchy--Schwarz.

## Full chapter map

1. [Decision rules, loss, and risk](#1-decision-rules-loss-and-risk)
2. [Bias, variance, and restricted optimality](#2-bias-variance-and-restricted-optimality)
3. [Rao--Blackwell, completeness, and UMVU estimation](#3-rao--blackwell-completeness-and-umvu-estimation)
4. [Score, information, Cramer--Rao, and two UMVU routes](#4-score-information-and-cramer--rao)
5. [Constructing estimators](#5-constructing-estimators)
6. [Recap](#6-recap)

<a id="w3-stop-1"></a>

## 1. Decision rules, loss, and risk

### 1.1 Statistical decision problem

Let $X\sim\mathbb{P}_\theta$, where $\theta\in\Theta$. A statistical decision problem specifies:

- an **action space** $\mathcal A$;
- a **decision rule** $\delta:\mathcal X\to\mathcal A$; and
- a **loss function** $L:\Theta\times\mathcal A\to[0,\infty)$.

The **risk function** of $\delta$ is

$$
R(\theta,\delta)
=\mathbb E_\theta[L(\theta,\delta(X))].
$$

Risk is not a single number until a value of $\theta$ or another method of aggregating over $\theta$ has been specified. Two decision rules may therefore have crossing risk functions.

Point estimation takes $\mathcal A=\Theta$ or a space containing the target. Under squared loss,

$$
L(\theta,a)=(a-\theta)^2.
$$

Testing takes $\mathcal A=\lbrace0,1\rbrace$, with the action interpreted as retaining or rejecting a null hypothesis. Week 4 develops that decision problem in detail.

### 1.2 Dominance and admissibility

A rule $\delta_1$ **dominates** $\delta_2$ if

$$
R(\theta,\delta_1)\leq R(\theta,\delta_2)
\quad\text{for every }\theta,
$$

with strict inequality for at least one $\theta$. A rule is **admissible** within a class $\mathcal D$ if no other rule in $\mathcal D$ dominates it.

Admissibility is a minimal requirement: it rules out uniform improvement. It does not say that the risk is small, that the rule is unique, or that it minimizes risk at every parameter value.

### Why unrestricted uniform optimality usually fails

For the Normal location model $X\sim\mathcal N(\theta,1)$, consider every constant rule $\delta_a(X)=a$. Its risk under squared loss is

$$
R(\theta,\delta_a)=(a-\theta)^2,
$$

which is zero at $\theta=a$. A rule that dominated every constant rule would need zero risk at every $\theta$. Because the Normal location laws are mutually absolutely continuous, no statistic can equal every possible $\theta$ almost surely. Pointwise minimization by the unknown action $a=\theta$ therefore does not produce an implementable decision rule.

### 1.3 Bayes and minimax criteria

A prior distribution $\Pi$ on $\Theta$ turns the risk function into the **Bayes risk**

$$
r(\Pi,\delta)=\int_\Theta R(\theta,\delta)\mkern3mu d\Pi(\theta).
$$

A **Bayes rule** minimizes this integrated risk. A **minimax rule** instead minimizes the worst case risk,

$$
\inf_\delta\sup_{\theta\in\Theta}R(\theta,\delta).
$$

These answer different questions: Bayes risk averages with respect to a specified distribution on the parameter, whereas minimax risk protects against its least favorable value.

For $X\sim\mathcal N(\theta,1)$ under squared loss and the prior $\theta\sim\mathcal N(0,v)$, the posterior mean and Bayes rule are

$$
\delta_v(X)=\frac{v}{1+v}X,
$$

and the Bayes risk is the posterior variance $v/(1+v)$. The identity rule $\delta_{\mathrm{id}}(X)=X$ has constant risk one. For every rule $\delta$ and every $v>0$,

$$
\sup_\theta R(\theta,\delta)
\geq r(\Pi_v,\delta)
\geq r(\Pi_v,\delta_v)
=\frac{v}{1+v}.
$$

Letting $v\to\infty$ shows that every rule has worst case risk at least one, so $X$ is minimax. The increasingly diffuse priors provide lower bounds; no improper prior is needed for the argument.

### Checkpoint 1

1. Give an example in which two estimators have crossing squared error risk functions.
2. Can an admissible rule perform badly for some parameter values?
3. Why is the action $a=\theta$ not an estimator?
4. How do Bayes and minimax criteria turn a risk function into a scalar objective?

> [!TIP]
> **AI interaction 1 — Does admissible mean good?**
>
> Copy the prompt below into an AI interface and audit the response.

```text
In the model X ~ N(theta,1), compare the constant estimator delta_a(X)=a
with the estimator delta(X)=X under squared error loss. A student says that
any admissible estimator must have lower risk than every inadmissible estimator
at every theta.

Evaluate the claim using risk functions. Distinguish admissibility, dominance,
pointwise risk minimization, and uniform optimality. Do not use "admissible" as
a synonym for "best."
```

**Audit question:** Does the response compare whole risk functions or only risk at one convenient parameter value?

<a id="w3-stop-2"></a>

## 2. Bias, variance, and restricted optimality

Let $T(X)$ estimate a scalar target $\tau(\theta)$. Its bias and mean squared error are

$$
\mathrm{Bias}_\theta(T)
=\mathbb E_\theta[T]-\tau(\theta),
$$

$$
\mathrm{MSE}_\theta(T)
=\mathbb E_\theta[(T-\tau(\theta))^2]
=\mathbb V_\theta[T]
+\mathrm{Bias}_\theta(T)^2.
$$

The identity follows by adding and subtracting $\mathbb E_\theta[T]$. An estimator is **unbiased** for $\tau(\theta)$ if

$$
\mathbb E_\theta[T]=\tau(\theta)
\quad\text{for every }\theta\in\Theta.
$$

Within the class of unbiased estimators with finite variance, minimizing MSE is equivalent to minimizing variance.

An estimator $T^\star$ is **uniformly minimum variance unbiased** (UMVU) for $\tau(\theta)$ if it is unbiased and

$$
\mathbb V_\theta[T^\star]
\leq
\mathbb V_\theta[U]
$$

for every $\theta$ and every other unbiased estimator $U$ with finite variance.

Unbiasedness is a restriction, not an automatic virtue. A biased estimator can have smaller MSE, and an unbiased estimator need not exist for every target.

### Checkpoint 2

Suppose $T$ has bias $b(\theta)$ and variance $v(\theta)$.

1. When can a biased estimator have lower MSE than an unbiased estimator?
2. Does minimum variance among unbiased estimators imply minimum MSE among all estimators?
3. Which part of the definition of UMVU is uniform in $\theta$?

<a id="w3-stop-3"></a>

## 3. Rao--Blackwell, completeness, and UMVU estimation

### 3.1 Rao--Blackwellization

Let $U(X)$ be unbiased for $\tau(\theta)$ with finite variance, and let $T(X)$ be sufficient for $\theta$. Define

$$
U^\star(T)=\mathbb E_\theta[U(X)\mid T(X)].
$$

Sufficiency ensures that a version of this conditional expectation can be chosen as a statistic that does not depend on the unknown parameter.

**Rao--Blackwell theorem.** The statistic $U^\star(T)$ is unbiased for $\tau(\theta)$ and, for every $\theta$,

$$
\mathbb V_\theta[U^\star(T)]\leq\mathbb V_\theta[U].
$$

Equality holds at a given $\theta$ if and only if $U=U^\star(T)$ almost surely under $\mathbb{P}_\theta$.

**Proof.** The tower property gives

$$
\mathbb E_\theta[U^\star(T)]
=\mathbb E_\theta\mkern-3mu\left[\mathbb E_\theta(U\mid T)\right]
=\mathbb E_\theta[U]
=\tau(\theta),
$$

so $U^\star(T)$ is unbiased. Let $R=U-U^\star(T)$. Because $\mathbb E_\theta[R\mid T]=0$, for every square integrable function $h(T)$,

$$
\mathbb E_\theta[R h(T)]
=\mathbb E_\theta\mkern-3mu\left[
\mathbb E_\theta(R\mid T)h(T)
\right]
=0.
$$

Taking $h(T)=U^\star(T)-\tau(\theta)$ and expanding $U-\tau(\theta)=R+\lbrace U^\star(T)-\tau(\theta)\rbrace$ yields the orthogonal decomposition

$$
\mathbb E_\theta[(U-\tau(\theta))^2]
=\mathbb E_\theta[(U-U^\star(T))^2]
+\mathbb E_\theta[(U^\star(T)-\tau(\theta))^2].
$$

Both estimators are unbiased, so this is

$$
\mathbb V_\theta[U]
=\mathbb E_\theta[(U-U^\star(T))^2]
+\mathbb V_\theta[U^\star(T)].
$$

The middle term is nonnegative, proving the variance inequality. It equals zero exactly when $U=U^\star(T)$ almost surely under $\mathbb{P}_\theta$. This is the same projection identity used for conditional expectation in Week 1; sufficiency is what makes the projected rule a parameter free statistic.

### 3.2 Completeness

A statistic $T$ is **complete** for $\lbrace\mathbb{P}_\theta:\theta\in\Theta\rbrace$ if, for every measurable $g$ for which the expectations exist,

$$
\mathbb E_\theta[g(T)]=0
\quad\text{for every }\theta
$$

implies

$$
g(T)=0
\quad \mathbb{P}_\theta\text{-almost surely for every }\theta.
$$

Completeness says that a function of $T$ cannot be unbiased for zero throughout the model unless it is the zero function. It is a property of the family of distributions of $T$, not of a realized value.

### 3.3 Bernoulli and exponential family completeness

If $T=\sum_iX_i\sim\mathsf{Binomial}(n,\theta)$, $\theta\in(0,1)$, and $\mathbb E_\theta[g(T)]=0$ for every $\theta$, then

$$
\sum_{t=0}^n g(t)\binom nt\theta^t(1-\theta)^{n-t}=0.
$$

Divide by $(1-\theta)^n$ and set $z=\theta/(1-\theta)>0$:

$$
\sum_{t=0}^n g(t)\binom nt z^t=0
\quad\text{for every }z>0.
$$

A polynomial that vanishes on an interval has every coefficient equal to zero. Hence $g(t)=0$ for $t=0,\ldots,n$, proving completeness.

#### Full rank exponential family theorem

Suppose $X_1,\ldots,X_n$ are i.i.d. from a common support exponential family

$$
f(x;\eta)=h(x)\exp\lbrace\eta't(x)-A(\eta)\rbrace,
\qquad
\eta\in\mathcal H\subseteq\mathbb R^k,
$$

under the usual integrability conditions. If $\mathcal H$ contains a nonempty open subset of $\mathbb R^k$, then

$$
T_n=\sum_{i=1}^n t(X_i)
$$

is complete and sufficient for $\eta$. Sufficiency follows from factorization. For completeness, the condition $\mathbb E_\eta[g(T_n)]=0$ says that the Laplace transform of a signed measure induced by $g(T_n)$ vanishes on an open set; uniqueness of the Laplace transform forces that signed measure, and hence $g(T_n)$, to vanish. The Bernoulli polynomial proof above is the finite support version of this argument.

To apply the theorem:

1. write the model in natural exponential family form with parameter independent support;
2. identify the natural statistic $t(X)$ and its sample sum $T_n$;
3. determine the image $\mathcal H$ of the natural parameter map; and
4. verify that $\mathcal H$ contains a nonempty open subset of the full ambient space $\mathbb R^k$.

If the natural parameter image lies only on a lower dimensional curve or boundary, the theorem is silent; failure of its open set condition is not itself a proof of incompleteness.

For the Bernoulli model, with $x\in\lbrace0,1\rbrace$,

$$
p^x(1-p)^{1-x}
=\exp\lbrace x\eta-A(\eta)\rbrace,
\qquad
\eta=\log\frac{p}{1-p},
\qquad
A(\eta)=\log(1+e^\eta).
$$

As $p$ ranges over $(0,1)$, $\eta$ ranges over all of $\mathbb R$. Thus the open set condition holds and $\sum_iX_i$ is complete and sufficient. This theorem gives a quick verification; the preceding polynomial proof shows why the conclusion is true in this model.

### 3.4 Lehmann--Scheffe theorem

Suppose $T$ is complete and sufficient for $\theta$. If $h(T)$ is unbiased for $\tau(\theta)$ and has finite variance, then $h(T)$ is the unique unbiased function of $T$ and is UMVU for $\tau(\theta)$.

**Proof.** Rao--Blackwellizing any unbiased estimator $U$ produces an unbiased function $U^\star(T)$. Since both $U^\star(T)$ and $h(T)$ are unbiased,

$$
\mathbb E_\theta[U^\star(T)-h(T)]=0
$$

for every $\theta$. Completeness gives $U^\star(T)=h(T)$ almost surely. Rao--Blackwell then gives

$$
\mathbb V_\theta[h(T)]
=\mathbb V_\theta[U^\star(T)]
\leq\mathbb V_\theta[U].
$$

### 3.5 Ancillarity as a contrast

A statistic $A(X)$ is **ancillary** if its sampling distribution does not depend on $\theta$. A sufficient statistic concentrates the model's information about $\theta$; an ancillary statistic has a parameter free distribution. Neither property by itself implies the other.

For $X_i\sim\mathcal N(\mu,1)$, the sample variance $S_n^2$ is ancillary for $\mu$ because $(n-1)S_n^2\sim\chi^2_{n-1}$.

### 3.6 Structural consequences: Bahadur and Basu

Two theorems explain why complete sufficient statistics have an especially strong structural role.

- **Bahadur's theorem:** in the usual dominated setting, a complete sufficient statistic is minimal sufficient, with statistics identified up to common null sets.
- **Basu's theorem:** a complete sufficient statistic is independent of every ancillary statistic.

**Proof map for Basu.** For an ancillary statistic $A$ and a measurable set $B$, sufficiency makes

$$
g_B(T)=\mathbb{P}_\theta(A\in B\mid T)
$$

a function that can be chosen independently of $\theta$. Ancillarity makes $\mathbb E_\theta[g_B(T)]=\mathbb{P}_\theta(A\in B)$ constant in $\theta$. Completeness therefore forces $g_B(T)$ to equal that constant almost surely, which is precisely independence.

In the model $X_i\sim\mathcal N(\mu,1)$, $\overline X_n$ is complete and sufficient for $\mu$, whereas $S_n^2$ is ancillary. Basu's theorem therefore gives $\overline X_n\perp S_n^2$. This fixed variance example isolates the theorem's logic. The analogous independence result when both $\mu$ and $\sigma^2$ are unknown is established by the Normal sampling theorem, not by treating $S_n^2$ as ancillary for the two parameter model.

> [!TIP]
> **AI interaction 2 — Repair the implication diagram**
>
> Copy the prompt below into an AI interface and audit the response.

```text
Compare sufficient, minimal sufficient, complete, ancillary, unbiased, and
UMVU. For each concept, give its definition and one valid implication. Then
evaluate these claims:

1. completeness implies sufficiency;
2. sufficiency implies completeness;
3. a complete sufficient statistic is minimal sufficient;
4. an unbiased function of a complete sufficient statistic is UMVU;
5. a complete sufficient statistic is independent of every ancillary statistic.

Use the iid Bernoulli model where possible. State any theorem conditions.
```

### Checkpoint 3

1. Why is sufficiency needed in the Rao--Blackwell theorem?
2. Why does completeness provide uniqueness rather than existence?
3. What does Bahadur's theorem add to completeness plus sufficiency?
4. Why does Basu's theorem need both completeness and sufficiency?
5. Can a statistic be both ancillary and nonconstant?

<a id="w3-stop-4"></a>

## 4. Score, information, and Cramer--Rao

### 4.1 Regular likelihood identities

For a scalar parameter, define the score

$$
s_\theta(X)=\frac{\partial}{\partial\theta}\log f_X(X;\theta).
$$

The standard identities require conditions such as:

- a parameter independent support;
- differentiability in $\theta$;
- permission to interchange differentiation and integration; and
- finite score moments.

Under these conditions,

$$
\mathbb E_\theta[s_\theta(X)]=0.
$$

Indeed,

$$
\mathbb E_\theta[s_\theta(X)]
=\int \frac{\partial_\theta f_X(x;\theta)}{f_X(x;\theta)}
f_X(x;\theta)\mkern3mu dx
=\frac{\partial}{\partial\theta}\int f_X(x;\theta)\mkern3mu dx
=0.
$$

The **Fisher information** is

$$
\mathcal I_n(\theta)
=\mathbb E_\theta[s_\theta(X)^2].
$$

Under the additional conditions needed to differentiate the score,

$$
\mathcal I_n(\theta)
=-\mathbb E_\theta
\left[
\frac{\partial^2}{\partial\theta^2}\log f_X(X;\theta)
\right].
$$

To see the second identity rather than memorize it, differentiate the mean zero score equation:

$$
0
=\frac{\partial}{\partial\theta}
\int s_\theta(x)f_X(x;\theta)\mkern3mu dx
=\mathbb E_\theta[\partial_\theta s_\theta(X)]
+\mathbb E_\theta[s_\theta(X)^2].
$$

The product rule uses $\partial_\theta f_X=s_\theta f_X$. Thus expected negative curvature equals score variance only when both differentiations under the integral are justified. Parameter dependent support is a standard way for the calculation to fail.

For an i.i.d. sample, information adds:

$$
\mathcal I_n(\theta)=n\mathcal I_1(\theta).
$$

> [!IMPORTANT]
> **Board work 1 — Score and information identities**
>
> Starting from $\int f_X(x;\theta)\mkern3mu dx=1$, derive the mean zero score identity and then the information identity. At each differentiation, state the interchange and support conditions being used. Finally, attempt the score calculation for an i.i.d. $\mathsf{Uniform}(0,\theta)$ sample and identify why the regular derivation fails.

<a id="w3-stop-5"></a>

### 4.2 Cramer--Rao lower bound

Let $U(X)$ have finite variance, and suppose its expectation is differentiable with the derivative passing under the integral. Then

$$
\frac{d}{d\theta}\mathbb E_\theta[U]
=\mathbb E_\theta[U s_\theta(X)]
=\mathrm{Cov}_\theta(U,s_\theta),
$$

where the last equality uses the mean zero score. Cauchy--Schwarz gives

$$
\left(
\frac{d}{d\theta}\mathbb E_\theta[U]
\right)^2
\leq
\mathbb V_\theta[U]\mkern3mu\mathcal I_n(\theta).
$$

Therefore

$$
\mathbb V_\theta[U]
\geq
\frac{
\left(\dfrac{d}{d\theta}\mathbb E_\theta[U]\right)^2
}{\mathcal I_n(\theta)}.
$$

If $U$ is unbiased for a differentiable target $\tau(\theta)$, then

$$
\mathbb V_\theta[U]
\geq
\frac{\tau'(\theta)^2}{\mathcal I_n(\theta)}.
$$

Equality in the bound requires $U-\mathbb E_\theta[U]$ to be proportional to the score almost surely. The bound may be unattainable, and it is not valid when its regularity conditions fail.

### 4.3 Matrix information and nuisance parameters

For $\theta\in\mathbb R^p$, define the score vector and information matrix by

$$
s_\theta(X)=\nabla_\theta\log f_X(X;\theta),
\qquad
\mathcal I_n(\theta)=\mathbb E_\theta[s_\theta(X)s_\theta(X)'].
$$

Under the multivariate versions of the regularity conditions above, let $U\in\mathbb R^m$ be unbiased for $\tau(\theta)\in\mathbb R^m$, and let $D\tau(\theta)$ be its $m\times p$ Jacobian. If $\mathcal I_n(\theta)$ is nonsingular, then

$$
\mathrm{Cov}_\theta(U)
\succeq
D\tau(\theta)\mathcal I_n(\theta)^{-1}D\tau(\theta)',
$$

where $A\succeq B$ means that $A-B$ is positive semidefinite. For a scalar target this becomes

$$
\mathbb V_\theta[U]
\geq
\nabla\tau(\theta)'\mathcal I_n(\theta)^{-1}\nabla\tau(\theta).
$$

This form automatically accounts for nuisance parameters through the inverse of the full information matrix. Using only one diagonal entry of $\mathcal I_n(\theta)$ is generally invalid unless additional orthogonality makes that reduction legitimate.

### 4.4 Constructing a UMVU estimator: two routes

The UMVU definition compares one candidate with the entire class of unbiased estimators, so a direct variance comparison is rarely practical. The preceding results supply two general certification routes.

| Route | Construction and verification | Conclusion | Main limitation |
|---|---|---|---|
| Complete sufficiency | Find a complete sufficient statistic $T$, then construct an unbiased function $h(T)$ of the target | Lehmann--Scheffe makes $h(T)$ the unique UMVU estimator | A usable complete sufficient statistic may not exist or may be hard to identify |
| Cramer--Rao bound | Verify regularity, derive the variance bound for the target, then find an unbiased estimator attaining it for every parameter value | The attaining estimator is UMVU within the class covered by the bound | The bound may be invalid, unattainable, or not sharp |

The Bernoulli model demonstrates both routes with the same estimator. Let $X_1,\ldots,X_n$ be i.i.d. $\mathsf{Bernoulli}(\theta)$, $\theta\in(0,1)$, and put $T=\sum_iX_i$.

**Route 1: complete sufficiency and unbiasedness.** Factorization makes $T$ sufficient, and the polynomial argument in Section 3.3 makes it complete. Because

$$
\mathbb E_\theta\left[\frac Tn\right]=\theta,
$$

Lehmann--Scheffe implies that $T/n=\overline X_n$ is the unique UMVU estimator of $\theta$.

**Route 2: attainment of the Cramer--Rao bound.** The sample score and information are

$$
s_\theta(X)
=\frac{T-n\theta}{\theta(1-\theta)},
\qquad
\mathcal I_n(\theta)=\frac{n}{\theta(1-\theta)}.
$$

Every unbiased estimator of $\theta$ in this finite support model satisfies the regular differentiation step, so its variance is at least

$$
\frac1{\mathcal I_n(\theta)}
=\frac{\theta(1-\theta)}n.
$$

But

$$
\mathbb V_\theta(\overline X_n)
=\frac{\theta(1-\theta)}n,
$$

so $\overline X_n$ attains the bound for every $\theta$ and is therefore UMVU. The completeness route establishes uniqueness through model structure; the bound route establishes optimality through a universal variance benchmark. Neither route should be invoked unless its own assumptions have been verified.

> [!IMPORTANT]
> **Board work 2 — Two UMVU routes in the Bernoulli model**
>
> For $X_i\overset{\mathrm{iid}}{\sim}\mathsf{Bernoulli}(\theta)$:
>
> 1. reconstruct the complete sufficiency and unbiasedness argument for $\overline X_n$;
> 2. derive the sample score and show $\mathcal I_n(\theta)=n/[\theta(1-\theta)]$;
> 3. compute the Cramer--Rao bound for unbiased estimators of $\theta$;
> 4. verify that $\overline X_n$ attains it; and
> 5. state what each route proves and the assumptions that the other route does not use.
>
> End with a compact decision rule: try complete sufficiency when model structure supplies it; try the Cramer--Rao route when regular information calculations produce an attainable bound.

> [!TIP]
> **AI interaction 3 — Find the illegal differentiation**
>
> Copy the prompt below into an AI interface and audit the response.

```text
A proposed proof says that the score has mean zero in every differentiable
parametric model because the derivative of the integral of a density is zero.

Audit the proof. State conditions that justify each equality, then apply the
argument to an iid Uniform(0,theta) sample and identify the first invalid step.
Explain why differentiability of the smooth part of the density is not enough.
```

### Checkpoint 4

Classify each statement as always true, true under additional conditions, or false.

1. Every score has expectation zero.
2. Fisher information is the expected negative Hessian.
3. Every unbiased estimator attains the Cramer--Rao bound.
4. Attaining the Cramer--Rao bound proves UMVU within the bound's class.
5. Failure of the regular Cramer--Rao calculation proves that no useful lower bound exists.
6. With nuisance parameters, the reciprocal of the target parameter's diagonal information entry is always the correct variance bound.
7. In the Bernoulli example, what distinct fact is supplied by complete sufficiency and by Cramer--Rao bound attainment?

<a id="w3-stop-6"></a>

## 5. Constructing estimators

Optimality criteria compare estimators after they have been proposed. Two broad construction principles are method of moments and maximum likelihood.

### 5.1 Method of moments

Sample averages provide the empirical analogues of population moments. For any chosen function $g$, replace

$$
\mathbb E_\theta[g(X)]
\qquad\text{by}\qquad
\frac1n\sum_{i=1}^n g(X_i).
$$

Method of moments chooses functions $g$ whose population expectations identify the parameter and then equates those expectations to their sample averages.

Suppose a $k$ dimensional parameter satisfies population moment equations

$$
\mathbb E_\theta[g(X)]=m(\theta),
$$

where $m$ is one-to-one on the parameter space. The method of moments estimator replaces population moments by sample moments and solves

$$
\frac1n\sum_{i=1}^n g(X_i)=m(\widehat\theta_{\mathrm{MM}}).
$$

The method requires the moments to exist, the sample equation to have an admissible solution, and the population equation to identify the target.

For example, if $X\sim\mathcal N(\mu,\sigma^2)$, the choices $g_1(x)=x$ and $g_2(x)=x^2$ give

$$
\mathbb E[X]=\mu,
\qquad
\mathbb E[X^2]=\mu^2+\sigma^2.
$$

Solving the sample moment equations yields

$$
\widehat\mu=\overline X_n,
\qquad
\widehat\sigma^2
=\frac1n\sum_{i=1}^nX_i^2-\overline X_n^2
=\frac1n\sum_{i=1}^n(X_i-\overline X_n)^2.
$$

This is the sample variance with denominator $n$; the usual version with denominator $n-1$ is a finite sample unbiasedness correction, not the direct method of moments estimator.

#### Identification by the selected moments

Identification of the full probability model does not guarantee that a chosen moment condition identifies the parameter. This is the [Week 2 identification question](week-02.md#2-identification-before-estimation) applied to the reduced map $\theta\mapsto\mathbb E_\theta[g(X)]$. Consider

$$
X\sim\mathcal N(\theta,\theta^2),
\qquad
\theta\in\mathbb R\setminus\lbrace0\rbrace.
$$

The model is identified because $\mathbb E_\theta[X]=\theta$. If one uses only $g(x)=x^2$, however,

$$
\mathbb E_\theta[X^2]=2\theta^2,
$$

so the population moment identifies $|\theta|$ but not its sign. The sample equation has the two solutions

$$
\widehat\theta
=\pm\sqrt{\frac1{2n}\sum_{i=1}^nX_i^2}.
$$

Using $g(x)=x$ instead gives $\mathbb E_\theta[X]=\theta$ and the uniquely identifying equation $\widehat\theta=\overline X_n$. The relevant requirement is injectivity of the selected moment map, not merely identification of the underlying parametric family.

#### Different moments, different estimators

Even when several moments identify the same parameter, they need not produce the same estimator. If $X\sim\mathsf{Uniform}(0,\theta)$ with $\theta>0$, then

$$
\mathbb E_\theta[X]=\frac\theta2,
\qquad
\mathbb E_\theta[X^2]=\frac{\theta^2}{3}.
$$

The first and second raw moments therefore give, respectively,

$$
\widehat\theta_1=2\overline X_n,
\qquad
\widehat\theta_2=\sqrt{\frac3n\sum_{i=1}^nX_i^2}.
$$

Both population moment maps identify $\theta$ on $(0,\infty)$, but the estimators differ in finite samples and their bias, variance, and robustness require separate comparison. More generally, the empirical moment can fall outside the range of the population moment map, leaving no admissible sample solution, and moments must exist before the method is defined. Even knowledge of every moment need not identify a distribution without additional determinacy conditions.

#### Moment restrictions beyond parametric models

The moment principle does not require a completely specified parametric distribution. If

$$
X=\mu+\varepsilon,
\qquad
\mathbb E[\varepsilon]=0,
$$

while the remainder of the error law is unrestricted, the moment restriction identifies $\mu=\mathbb E[X]$ and yields $\widehat\mu=\overline X_n$. Likewise, for an otherwise unrestricted distribution with c.d.f. $F$,

$$
F(x)=\mathbb E[\mathbf 1\lbrace X\leq x\rbrace]
$$

has the empirical analogue

$$
\widehat F_n(x)=\frac1n\sum_{i=1}^n\mathbf 1\lbrace X_i\leq x\rbrace.
$$

These examples connect method of moments to semiparametric and nonparametric plug in estimation: the model restriction must identify the target, but it need not identify a finite dimensional distribution family.

### 5.2 Maximum likelihood

The maximum likelihood estimator satisfies

$$
\widehat\theta_{\mathrm{ML}}
\in\arg\max_{\theta\in\Theta}\ell_n(\theta).
$$

An interior differentiable maximizer satisfies a score equation, but boundary and nonsmooth maxima need not. For a transformed parameter $\tau=r(\theta)$, define the induced profile likelihood

$$
L^\star(\tau)=\sup_{\theta:\mkern3mu r(\theta)=\tau}L(\theta).
$$

The **invariance property of maximum likelihood** gives

$$
\widehat{r(\theta)}_{\mathrm{ML}}=r(\widehat\theta_{\mathrm{ML}}).
$$

For a one-to-one $r$, this is ordinary reparameterization. For a many to one map, the profile definition and possible nonuniqueness of maximizers must be kept explicit.

Neither construction method, by itself, proves unbiasedness, minimum risk, consistency, or asymptotic normality. Those are separate claims requiring separate arguments.

> [!IMPORTANT]
> **Board work 3 — Three estimators in the Uniform model**
>
> For $X_i\overset{\mathrm{iid}}{\sim}\mathsf{Uniform}(0,\theta)$, compare
>
> $$
> 2\overline X_n,
> \qquad
> X_{(n)},
> \qquad
> \frac{n+1}{n}X_{(n)}.
> $$
>
> Identify which is method of moments, which is maximum likelihood on the probability one event $X_{(n)}>0$, and which are unbiased. Then compute their variances or MSEs. Explain why the regular Cramer--Rao calculation from Section 4 is not the right benchmark for this model.

> [!TIP]
> **AI interaction 4 — Construction is not optimality**
>
> Copy the prompt below into an AI interface and audit the response.

```text
Compare method of moments and maximum likelihood as ways to construct an
estimator. For each method, list what follows from the definition and what
requires a separate theorem. For method of moments, distinguish identification
of the full model from injectivity of the selected moment map, using the model
Normal(theta,theta^2) and the choices g(x)=x and g(x)=x^2. Then use the
Uniform(0,theta) model to show why the
claims "the MLE solves the score equation," "the MLE is unbiased," and "the
MLE attains the Cramer--Rao bound" are not automatic.
```

### Checkpoint 5

1. Can a model be identified while a selected moment condition fails to identify its parameter?
2. Can a method of moments equation have multiple solutions or no admissible solution?
3. Can an MLE fail to exist even when the likelihood is bounded above?
4. Does invariance of the MLE imply unbiasedness after transformation?
5. Which Week 5 results will be needed to establish consistency of these estimators?

## 6. Recap

- Statistical procedures can be compared only after specifying actions, loss, and a parameter indexed risk function.
- Admissibility means absence of a uniform improvement; it does not mean uniformly smallest risk.
- Bayes risk averages risk under a prior; minimax risk controls its worst case value.
- Under squared loss, MSE equals variance plus squared bias.
- Rao--Blackwellization uses sufficiency to improve an unbiased estimator.
- Completeness supplies uniqueness; together with sufficiency it yields the Lehmann--Scheffe theorem.
- Bahadur links complete sufficiency to minimal sufficiency; Basu links it to ancillaries through independence.
- Score and information identities are regular model results with substantive assumptions.
- Scalar and matrix Cramer--Rao bounds follow from score identities and covariance inequalities; their regularity conditions and treatment of nuisance parameters matter.
- A UMVU estimator can be certified by constructing an unbiased function of a complete sufficient statistic or by showing that an unbiased estimator attains a valid Cramer--Rao bound for every parameter value.
- A method of moments construction requires the selected population moment map to identify the target; different identifying moments can produce different estimators.
- Method of moments and maximum likelihood construct candidates. Their statistical properties must be proved separately.

## Notation introduced this week

- $\mathcal A$: action space.
- $\delta$: decision rule.
- $L(\theta,a)$ and $R(\theta,\delta)$: loss and risk.
- $r(\Pi,\delta)$: Bayes risk under prior $\Pi$.
- $\tau(\theta)$: scalar target.
- $\mathrm{Bias}_\theta$: bias.
- $\mathrm{MSE}_\theta$: mean squared error.
- $s_\theta(X)$: score.
- $\mathcal I_n(\theta)$: scalar Fisher information or the Fisher information matrix in the sample.
- $\widehat\theta_{\mathrm{MM}}$ and $\widehat\theta_{\mathrm{ML}}$: method of moments and maximum likelihood estimators.

## References

- Jun Shao, *Mathematical Statistics*, Springer, 2003.
- Robert W. Keener, *Theoretical Statistics: Topics for a Core Course*, Springer, 2010.
- Peter J. Bickel and Kjell A. Doksum, *Mathematical Statistics: Basic Ideas and Selected Topics*, vol. I, 2nd ed., CRC Press, 2015.
- Erich Lehmann and George Casella, *Theory of Point Estimation*, 2nd ed., Springer, 1998.
