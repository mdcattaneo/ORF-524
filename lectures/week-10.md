# ORF 524: Statistical Theory and Methods

## Week 10: Regression applications and two step estimation

**Status:** Ready for Fall 2026  
**Last updated:** August 30, 2026  
**Meetings:** Monday, November 2, and Wednesday, November 4, 2026

<!-- Legacy sources: 11 ORF 524/lectures/ORF524-Fall2025-Lecture-Week-09.tex and
ORF524-Fall2025-Lecture-Week-09_Lecture.pdf; 09 ORF 524/lectures/
ORF524-Fall2024-Lecture-Week-09_Notes.pdf; 09 ORF 524/assignments/pset7.tex; and
11 ORF 524 Problem Set 4. -->

## Central question

How does the general M- and Z-estimation theory from Week 9 translate into regression procedures, and how does estimating a nuisance parameter in a first step change the second step distribution?

## Learning goals

By the end of the week, you should be able to:

1. distinguish the best linear predictor target from a correctly specified linear conditional mean;
2. derive the OLS asymptotic linear representation and heteroskedastic sandwich variance;
3. analyze nonlinear least squares and binary response likelihood as M- and Z-procedures;
4. compare unweighted, infeasibly weighted, feasibly weighted, and continuously updated estimating equations;
5. derive the influence function of a finite dimensional two step Z-estimator; and
6. recognize orthogonality as the condition that removes the first step effect at first order.

## In-class route

This is one continuous route through the week, not a division into class meetings. Proceed in order, stop wherever the discussion ends, and resume at the next unfinished stop. The linked sections contain the full student facing exposition; the route identifies the concepts and activities to foreground live.

| Stop | Live focus | Mode |
|---|---|---|
| **W10.1** | [Best linear prediction, OLS, and robust variance](#w10-stop-1) | Board work + AI target audit + Checkpoint 1 |
| **W10.2** | [Nonlinear regression and binary response](#w10-stop-2) | Board work + efficiency comparison + two step motivation + Checkpoint 2 |
| **W10.3** | [Sequential two step influence adjustment](#w10-stop-3) | Board work + first step adjustment |
| **W10.4** | [Stacked two step representation](#w10-stop-4) | Board work + AI nuisance audit + Checkpoint 3 |
| **W10.5** | [Orthogonality and feasible weighting](#w10-stop-5) | Discuss + transfer + Checkpoint 4 |

Stop W10.2 should end with the construction of feasible WNLS from a preliminary NLS estimator, making the need for a general two step theory concrete. Stops W10.3--W10.4 then present the same first step adjustment sequentially and through a stacked system; their sign conventions and derivative blocks should be compared explicitly. Stop W10.5 is the finite dimensional bridge to the orthogonal scores and cross fitting of Week 14.

## How to use this chapter

**Prepare:** Review the M/Z consistency and linearization theorems in Week 9. Re derive the least squares sandwich variance without assuming a linear conditional mean.

**In class:** Follow the continuous route above through linear and nonlinear regression, binary response, the two step influence adjustment, stacked estimation, feasible weighting, and orthogonality. Prepared responses supply an alternative for every live prompt if AI is unavailable.

**Review:** For every estimator, write down the population target, sample equation, derivative matrix, moment variance, and influence function. Do not describe an estimator as efficient until the model and competitor class have been specified.

**Practice:** The matching assignment complements this chapter with ungraded problems for consolidating the material, extending selected ideas, and preparing for cumulative review.

**Prerequisites:** Weeks 1--6 and Week 9; conditional expectation; best linear prediction; M- and Z-estimation; multivariate CLT; delta method; matrix differentiation; sandwich variance.

## Full chapter map

1. [Best linear prediction and OLS](#1-best-linear-prediction-and-ols)
2. [Nonlinear regression and binary response](#2-nonlinear-regression-and-binary-response)
3. [Finite dimensional two step Z-estimation](#3-finite-dimensional-two-step-z-estimation)
4. [Orthogonality and feasible weighting](#4-orthogonality-and-feasible-weighting)
5. [Recap](#5-recap)

<a id="w10-stop-1"></a>

## 1. Best linear prediction and OLS

### 1.1 The population target

Let $W=(Y,X)$ be i.i.d. across observations, where $Y\in\mathbb R$ and $X\in\mathbb R^d$ includes an intercept when one is desired. Suppose $\mathbb E[Y^2+\Vert X\Vert^2]<\infty$ and

$$
H=\mathbb E[XX']
$$

is positive definite. The **best linear predictor coefficient** is

$$
\beta_0
=\arg\min_{\beta\in\mathbb R^d}
\mathbb E[(Y-X'\beta)^2]
=H^{-1}\mathbb E[XY].
$$

Writing $\varepsilon=Y-X'\beta_0$, the population normal equation is

$$
\mathbb E[X\varepsilon]=0.
$$

This orthogonality identifies the best linear approximation even when

$$
\mathbb E[Y\mid X]\neq X'\beta_0.
$$

The stronger condition $\mathbb E[\varepsilon\mid X]=0$ says that the conditional mean is correctly specified as linear. It is not required to define or consistently estimate the best linear predictor target.

### 1.2 Sample OLS and its first order distribution

The OLS estimator is both an M-estimator and a Z-estimator:

$$
\widehat\beta_n
\in\arg\min_\beta\mathbb{E}_n(Y-X'\beta)^2,
\qquad
0=\mathbb{E}_n[X(Y-X'\widehat\beta_n)].
$$

Let

$$
\widehat H_n=\mathbb{E}_n[XX'].
$$

Whenever $\widehat H_n$ is nonsingular,

$$
\widehat\beta_n=\widehat H_n^{-1}\mathbb{E}_n[XY].
$$

The derivative of the estimating equation is $\Psi_\beta(\beta_0)=-H$. Therefore the Week 9 Z-estimator theorem gives the correctly signed representation

$$
\sqrt n(\widehat\beta_n-\beta_0)
=H^{-1}\frac1{\sqrt n}\sum_{i=1}^nX_i\varepsilon_i
+o_{\mathbb{P}}(1).
$$

If $\mathbb E[\Vert X\varepsilon\Vert^2]<\infty$, then

$$
\sqrt n(\widehat\beta_n-\beta_0)
\rightsquigarrow
\mathcal N(0,V),
$$

where

$$
V=H^{-1}\Sigma H^{-1},
\qquad
\Sigma=\mathbb E[XX'\varepsilon^2].
$$

With residuals $\widehat\varepsilon_i=Y_i-X_i'\widehat\beta_n$, the heteroskedasticity robust estimator is

$$
\widehat V_n
=\widehat H_n^{-1}\widehat\Sigma_n\widehat H_n^{-1},
\qquad
\widehat\Sigma_n
=\mathbb{E}_n[XX'\widehat\varepsilon^2].
$$

Fixed dimensional leverage corrections such as HC1--HC3 alter the finite sample formula but not its first order target under the usual conditions.

If the stronger restriction

$$
\mathbb E[\varepsilon^2\mid X]=\sigma^2
$$

holds, then $\Sigma=\sigma^2H$ and $V=\sigma^2H^{-1}$. This is a simplification under an additional model restriction, not a general ordering between two variance matrices.

> [!IMPORTANT]
> **Board work 1 — OLS without a correct linear conditional mean**
>
> Derive the best linear predictor normal equation, the OLS influence function, and the heteroskedastic sandwich variance without assuming $\mathbb E[Y\mid X]=X'\beta_0$. Then identify the two additional conditional restrictions that respectively justify a linear conditional mean interpretation and the homoskedastic variance simplification.

### Checkpoint 1

1. What parameter does OLS estimate when the conditional mean is nonlinear?
2. Which moment condition defines the best linear predictor?
3. Why is the influence function $H^{-1}X\varepsilon$ rather than $-H^{-1}X\varepsilon$?
4. Which assumption simplifies the sandwich variance to $\sigma^2H^{-1}$?

> [!TIP]
> **AI interaction 1 — Audit the regression target**
>
> Copy the prompt below into an AI interface and audit the response.

```text
An analyst fits OLS and reports that the regression consistently estimates
E[Y|X=x] because the OLS residual is orthogonal to X.

Audit the claim. Distinguish the best linear predictor target from a correctly
specified linear conditional mean. Derive the population normal equations and
give a nonlinear conditional mean example for which OLS still has a well defined
population target. State the assumptions needed for the robust sandwich limit.
```

**Audit question:** Does the response replace unconditional orthogonality $\mathbb E[X\varepsilon]=0$ by conditional mean independence without justification?

<a id="w10-stop-2"></a>

## 2. Nonlinear regression and binary response

### 2.1 Nonlinear least squares

Let $G:\mathbb R\to\mathbb R$ be known and smooth. The nonlinear least squares target is

$$
\beta_0
\in\arg\min_{\beta\in\Theta}
\mathbb E[\lbrace Y-G(X'\beta)\rbrace^2].
$$

Define

$$
q(X,\beta)
=\nabla_\beta G(X'\beta)
=\dot G(X'\beta)X.
$$

An interior sample minimizer solves

$$
0=\mathbb{E}_n\psi(W,\widehat\beta_n),
\qquad
\psi(W,\beta)=q(X,\beta)\lbrace Y-G(X'\beta)\rbrace.
$$

If the conditional mean is correctly specified,

$$
\mathbb E[Y\mid X]=G(X'\beta_0),
$$

then, writing $\varepsilon=Y-G(X'\beta_0)$,

$$
\Psi_\beta(\beta_0)
=\mathbb E[\nabla_{\beta'}\psi(W,\beta_0)]
=-\mathbb E[q(X,\beta_0)q(X,\beta_0)'],
$$

because the term multiplying $\varepsilon$ has conditional mean zero. Also,

$$
B
=\mathbb E[q(X,\beta_0)q(X,\beta_0)'\varepsilon^2].
$$

The sandwich variance is $[\Psi_\beta(\beta_0)]^{-1}B[\Psi_\beta(\beta_0)]^{-\prime}$. Under misspecification, nonlinear least squares still targets a best approximation when the population minimizer is identified, but the derivative $\Psi_\beta(\beta_0)$ also contains the expected residual times the second derivative of $G$. The correct specification simplification cannot then be used automatically.

### 2.2 Binary response: likelihood versus squared loss

Suppose $Y\in\lbrace0,1\rbrace$ and

$$
\mathbb{P}(Y=1\mid X)=p(X'\beta_0),
$$

where $p$ is a known smooth link such as the logistic or Normal c.d.f. Put $p_\beta(X)=p(X'\beta)$ and $\dot p_\beta(X)=\dot p(X'\beta)$. The conditional likelihood score is

$$
s(W,\beta)
=\frac{X\dot p_\beta(X)}{p_\beta(X)\lbrace1-p_\beta(X)\rbrace}
\lbrace Y-p_\beta(X)\rbrace.
$$

Nonlinear least squares instead uses

$$
\psi_{\mathrm{NLS}}(W,\beta)
=X\dot p_\beta(X)\lbrace Y-p_\beta(X)\rbrace.
$$

Under correct specification and identification, both equations target $\beta_0$. They generally have different asymptotic variances. Because

$$
\mathbb V(Y\mid X)=p_{\beta_0}(X)\lbrace1-p_{\beta_0}(X)\rbrace,
$$

multiplying the NLS estimating equation at each candidate $\beta$ by the inverse candidate variance $\lbrace p_\beta(X)[1-p_\beta(X)]\rbrace^{-1}$ produces the likelihood score. Thus the Bernoulli MLE is the optimally weighted member of this conditional moment family under the correct model. Probabilities bounded away from zero and one, together with appropriate moments, are a convenient sufficient condition for the regularity needed in this comparison.

Under misspecification, the likelihood and squared loss criteria can define different pseudo true parameters. “MLE versus NLS” is then not merely a variance comparison unless the targets have first been aligned.

### 2.3 Infeasible, feasible, and continuously updated reweighting

The preceding comparison assumes knowledge of the variance that determines the efficient weight. To expose the resulting feasibility problem, write

$$
v_0(X)
=p_{\beta_0}(X)\lbrace1-p_{\beta_0}(X)\rbrace,
\qquad
q_\beta(X)=X\dot p_\beta(X).
$$

#### Infeasible weighted nonlinear least squares

If $v_0(X)$ were known, the nonlinear analogue of GLS would be the infeasible weighted nonlinear least squares estimator (IWNLS), also called infeasible nonlinear GLS:

$$
\widehat\beta_{\mathrm{IWNLS}}
\in
\arg\min_\beta
\mathbb{E}_n
\left[
\frac{\lbrace Y-p_\beta(X)\rbrace^2}{v_0(X)}
\right].
$$

Because the weight is fixed with respect to the candidate $\beta$, an interior minimizer solves

$$
0 =
\mathbb{E}_n
\left[
\frac{q_\beta(X)}{v_0(X)}
\lbrace Y-p_\beta(X)\rbrace
\right].
$$

At $\beta_0$, this estimating function equals the Bernoulli score. If

$$
\mathcal I(\beta_0) =
\mathbb E
\left[
\frac{q_{\beta_0}(X)q_{\beta_0}(X)'}{v_0(X)}
\right],
$$

then its derivative matrix is $-\mathcal I(\beta_0)$ and its variance is $\mathcal I(\beta_0)$. Its asymptotic variance is therefore $\mathcal I(\beta_0)^{-1}$, the same efficiency bound attained by the correctly specified Bernoulli MLE. The estimator is infeasible precisely because $v_0(X)$ depends on the unknown $\beta_0$.

#### Feasible weighted nonlinear least squares

Start with a consistent preliminary estimator $\widetilde\beta_n$, such as NLS, and form the plug in conditional variance

$$
\widetilde v_n(X) =
p_{\widetilde\beta_n}(X)
\lbrace1-p_{\widetilde\beta_n}(X)\rbrace.
$$

The feasible weighted nonlinear least squares estimator (FWNLS), or nonlinear FGLS, freezes these estimated weights in the second step:

$$
\widehat\beta_{\mathrm{FWNLS}}
\in
\arg\min_\beta
\mathbb{E}_n
\left[
\frac{\lbrace Y-p_\beta(X)\rbrace^2}{\widetilde v_n(X)}
\right].
$$

Its second step equation is

$$
0 =
\mathbb{E}_n
\left[
\frac{q_\beta(X)}{\widetilde v_n(X)}
\lbrace Y-p_\beta(X)\rbrace
\right].
$$

This is a genuine two step estimator: the first step estimates the variance weights, and the second step estimates the regression parameter using those weights. Under correct specification and the regularity conditions developed in Sections 3 and 4, estimating the weights has no first order effect, so FWNLS is first order equivalent to infeasible WNLS and hence to the Bernoulli MLE.

#### Continuously updated reweighting

Instead of stopping after one feasible reweighting step, one can iterate. Given $\widehat\beta^{(t)}$, form

$$
v^{(t)}(X) =
p_{\widehat\beta^{(t)}}(X)
\lbrace1-p_{\widehat\beta^{(t)}}(X)\rbrace,
$$

and compute

$$
\widehat\beta^{(t+1)}
\in
\arg\min_\beta
\mathbb{E}_n
\left[
\frac{\lbrace Y-p_\beta(X)\rbrace^2}{v^{(t)}(X)}
\right],
$$

treating $v^{(t)}$ as fixed during that update. At a fixed point, the parameter and its implied variance weights agree, and the estimating equation is

$$
0 =
\mathbb{E}_n
\left[
\frac{q_\beta(X)}{p_\beta(X)\lbrace1-p_\beta(X)\rbrace}
\lbrace Y-p_\beta(X)\rbrace
\right].
$$

This continuously updated estimating equation is exactly the Bernoulli likelihood score equation. An interior Bernoulli MLE therefore solves it; when the score root uniquely determines the likelihood maximizer, the continuously updated estimator is the MLE. The iteration estimates the regression parameter and updates its implied conditional variance weights jointly.

There is an important distinction between continuously updating an estimating equation and minimizing the naive parameter weighted residual criterion

$$
\mathbb{E}_n
\left[
\frac{\lbrace Y-p_\beta(X)\rbrace^2}
{p_\beta(X)\lbrace1-p_\beta(X)\rbrace}
\right].
$$

Differentiating this criterion also differentiates the weight and produces an additional term. Its first order condition is therefore not the likelihood score equation and need not target $\beta_0$. The continuously updated procedure above is defined by the weighted moment equation—or, equivalently here, by maximizing the Bernoulli likelihood.

> [!IMPORTANT]
> **Board work 2 — From efficiency to two step estimation**
>
> For logistic $p(u)=e^u/(1+e^u)$:
>
> 1. simplify the likelihood score;
> 2. derive the NLS estimating equation;
> 3. derive infeasible WNLS and show that it has inverse information asymptotic variance;
> 4. construct feasible WNLS from a preliminary NLS estimator and identify its two steps;
> 5. show that continuously updating the variance weight in the estimating equation reproduces the likelihood score; and
> 6. explain why differentiating a parameter weighted least squares objective is a different procedure.

### Checkpoint 2

1. Which term in $\Psi_\beta(\beta_0)$ disappears under a correctly specified nonlinear conditional mean?
2. Why can two consistent estimating equations for the same target have different variances?
3. Why is variance weighted NLS infeasible even though the Bernoulli conditional variance has a known functional form?
4. Which two estimators constitute feasible WNLS?
5. Why does the continuously updated estimating equation equal the MLE score, while the naive parameter weighted least squares objective does not?
6. Under misspecification, must MLE and NLS target the same parameter?

<a id="w10-stop-3"></a>

## 3. Finite dimensional two step Z-estimation

### 3.1 Setup and identification

Let $\beta_0\in\mathbb R^p$ be the parameter of interest and $\gamma_0\in\mathbb R^q$ a nuisance parameter. Suppose they are identified by

$$
\mathbb E[m(W,\beta_0,\gamma_0)]=0,
\qquad
\mathbb E[g(W,\gamma_0)]=0,
$$

where $m\in\mathbb R^p$ and $g\in\mathbb R^q$. A two step estimator first solves

$$
\mathbb{E}_n g(W,\widehat\gamma_n)=o_{\mathbb{P}}(n^{-1/2})
$$

and then solves

$$
\mathbb{E}_n m(W,\widehat\beta_n,\widehat\gamma_n)
=o_{\mathbb{P}}(n^{-1/2}).
$$

Consistency still requires the global identification and uniform convergence arguments from Week 9. The calculations below begin after consistency has localized both estimators.

### 3.2 The first step contribution

Define the population derivative matrices

$$
M_\beta
=\mathbb E[\nabla_{\beta'}m(W,\beta_0,\gamma_0)],
\qquad
M_\gamma
=\mathbb E[\nabla_{\gamma'}m(W,\beta_0,\gamma_0)],
$$

and

$$
G_\gamma
=\mathbb E[\nabla_{\gamma'}g(W,\gamma_0)].
$$

Suppose $M_\beta$ and $G_\gamma$ are nonsingular and the local mean value and uniform derivative conditions from Week 9 hold. The first step has representation

$$
\sqrt n(\widehat\gamma_n-\gamma_0)
=-G_\gamma^{-1}
\frac1{\sqrt n}\sum_{i=1}^n g(W_i,\gamma_0)
+o_{\mathbb{P}}(1).
$$

Applying the coordinatewise mean value argument to the second step equation gives

$$
0
=\frac1{\sqrt n}\sum_{i=1}^n m(W_i,\beta_0,\gamma_0)
+M_\beta\sqrt n(\widehat\beta_n-\beta_0)
+M_\gamma\sqrt n(\widehat\gamma_n-\gamma_0)
+o_{\mathbb{P}}(1).
$$

Substitution yields

$$
\sqrt n(\widehat\beta_n-\beta_0)
=\frac1{\sqrt n}\sum_{i=1}^n\varphi_\beta(W_i)
+o_{\mathbb{P}}(1),
$$

where

$$
\varphi_\beta(W)
=-M_\beta^{-1}
\left[
m(W,\beta_0,\gamma_0)
-M_\gamma G_\gamma^{-1}g(W,\gamma_0)
\right].
$$

The second term inside the brackets is the first step adjustment. Omitting it treats $\widehat\gamma_n$ as known and generally gives the wrong variance.

If $\mathbb E[\Vert\varphi_\beta(W)\Vert^2]<\infty$, then

$$
\sqrt n(\widehat\beta_n-\beta_0)
\rightsquigarrow
\mathcal N(0,V_\beta),
\qquad
V_\beta=\mathbb E[\varphi_\beta(W)\varphi_\beta(W)'].
$$

A plug in estimator replaces the parameters and derivative matrices in $\varphi_\beta$ and takes the empirical second moment of the estimated influence functions.

> [!IMPORTANT]
> **Board work 3 — Sequential two step adjustment**
>
> Apply the coordinatewise mean value theorem to the first and second step equations, substitute the first step influence representation into the second equation, and verify every matrix sign and dimension. Identify the nuisance sensitivity matrix whose vanishing removes the first step contribution.

<a id="w10-stop-4"></a>

### 3.3 Stacking gives the same answer

Define

$$
h(W,\beta,\gamma) =
\begin{pmatrix}
m(W,\beta,\gamma)\\
g(W,\gamma)
\end{pmatrix}.
$$

The Jacobian is block triangular:

$$
\mathbb E[\nabla_{(\beta',\gamma')}h] =
\begin{pmatrix}
M_\beta&M_\gamma\\
0&G_\gamma
\end{pmatrix}.
$$

Applying the ordinary multivariate Z-estimator theorem to this stacked system and taking the first block produces the same $\varphi_\beta$. The sequential derivation is useful because it exposes exactly where nuisance estimation enters.

> [!IMPORTANT]
> **Board work 4 — Stacked two step influence function**
>
> Invert the block triangular Jacobian and recover the two step influence function. Then identify the variance formula when $m$ and $g$ are correlated.

> [!TIP]
> **AI interaction 2 — The nuisance was estimated, not known**
>
> Copy the prompt below into an AI interface and audit the response.

```text
Derive the influence function of a two step Z-estimator. Use one estimating
equation for beta that depends on gamma and a second equation for gamma.
Keep all matrix dimensions explicit and verify the signs by comparing the
sequential derivation with inversion of the stacked block triangular Jacobian.

Then audit the claim that replacing gamma by a consistent estimator never
changes first order variance. State the condition under which that claim is true.
```

**Audit question:** Does the response discard $M_\gamma\sqrt n(\widehat\gamma_n-\gamma_0)$ merely because $\widehat\gamma_n\to_{\mathbb{P}}\gamma_0$?

### Checkpoint 3

1. Why is consistency of $\widehat\gamma_n$ alone insufficient for ignoring the first step?
2. Which matrix measures sensitivity of the second step moment to the nuisance parameter?
3. How does correlation between $m$ and $g$ enter the variance?
4. When would using an independent first step sample change the covariance calculation?

<a id="w10-stop-5"></a>

## 4. Orthogonality and feasible weighting

### 4.1 Finite dimensional Neyman orthogonality

The second step moment is **Neyman orthogonal** to $\gamma$ at the truth if

$$
M_\gamma
=\left.
\nabla_{\gamma'}\mathbb E[m(W,\beta_0,\gamma)]
\right|_{\gamma=\gamma_0}
=0.
$$

Then

$$
\varphi_\beta(W)=-M_\beta^{-1}m(W,\beta_0,\gamma_0),
$$

so estimating a nuisance parameter at the $\sqrt n$ rate has no first order effect on the second step estimator. Orthogonality is local insensitivity of a population moment, not statistical independence and not a claim that the nuisance is unimportant in finite samples.

For slower nonparametric first steps, $M_\gamma=0$ removes the linear nuisance error, but a second order remainder must still be negligible. A typical requirement has the form

$$
\sqrt n\mkern3mu\Vert\widehat\gamma_n-\gamma_0\Vert^2=o_{\mathbb{P}}(1),
$$

or, with several nuisances, a product rate condition. Week 14 combines this idea with cross fitting.

### 4.2 Why feasible weighting can be first order equivalent

The binary response FWNLS estimator in Section 2.3 is a special case of the following argument. Consider a conditional mean model

$$
\mathbb E[Y-G(X'\beta_0)\mid X]=0
$$

and a weighted estimating equation

$$
m(W,\beta,\gamma)
=w(X,\gamma)q(X,\beta)\lbrace Y-G(X'\beta)\rbrace.
$$

If $w$ is smooth and the model is correctly specified, then

$$
\mathbb E[\nabla_{\gamma'}m(W,\beta_0,\gamma_0)]
=\mathbb E\mkern-3mu\left[
\nabla_{\gamma'}\lbrace w(X,\gamma_0)q(X,\beta_0)\rbrace
\mathbb E[\varepsilon\mid X]
\right]
=0.
$$

Thus estimating the weights from a regular preliminary model has no first order effect. This is why feasible weighted nonlinear least squares can share the first order distribution of its infeasible counterpart. The conclusion depends on the conditional moment restriction; it is not automatic for every plug in weight.

### Checkpoint 4

1. Does orthogonality mean that $\widehat\beta_n$ and $\widehat\gamma_n$ are independent?
2. Which derivative must vanish?
3. Why does conditional mean zero make estimated weighting orthogonal in the example?
4. What additional remainder control is needed when the nuisance converges slower than root $n$?

## 5. Recap

- OLS always needs a population target; the best linear predictor does not require a linear conditional mean.
- The OLS influence function is $H^{-1}X\varepsilon$, and its general variance has sandwich form.
- Nonlinear least squares is a generic M/Z procedure; correct specification simplifies its derivative matrix.
- In Bernoulli regression, infeasible WNLS uses the true conditional variance and attains the correctly specified MLE's first order efficiency.
- Feasible WNLS estimates and freezes those variance weights in a preliminary step, providing the motivating example for two step estimation.
- Continuously updating the estimating equation weight produces the Bernoulli score, but differentiating a parameter weighted residual criterion produces an additional term.
- A two step estimator inherits a first step influence adjustment through $M_\gamma G_\gamma^{-1}g$.
- Stacking and sequential mean value arguments give the same two step result.
- Neyman orthogonality sets $M_\gamma=0$ and removes the linear first step contribution.
- Feasible weighting is a concrete orthogonality example and a bridge to nonparametric nuisance estimation in Week 14.

## Notation introduced this week

- $H=\mathbb E[XX']$: population second moment matrix for OLS.
- $\varepsilon=Y-X'\beta_0$ or $Y-G(X'\beta_0)$: population residual.
- $q(X,\beta)=\nabla_\beta G(X'\beta)$: nonlinear regression gradient.
- $v_0(X)=p_{\beta_0}(X)\lbrace1-p_{\beta_0}(X)\rbrace$: Bernoulli conditional variance at the truth.
- $\widehat\beta_{\mathrm{IWNLS}}$ and $\widehat\beta_{\mathrm{FWNLS}}$: infeasible and feasible weighted nonlinear least squares estimators.
- $M_\beta$ and $M_\gamma$: derivatives of the second step population moment.
- $g(W,\gamma)$ and $G_\gamma$: first step moment and its population derivative.
- $\varphi_\beta$: influence function of the second step estimator.
- Neyman orthogonality: $M_\gamma=0$ at the true parameters.

## References

- Whitney Newey and Daniel McFadden, “Large Sample Estimation and Hypothesis Testing,” in *Handbook of Econometrics*, vol. 4, 1994.
- Aad van der Vaart, *Asymptotic Statistics*, Cambridge University Press, 1998.
- Jun Shao, *Mathematical Statistics*, Springer, 2003.
- Robert W. Keener, *Theoretical Statistics: Topics for a Core Course*, Springer, 2010.
