# ORF 524: Statistical Theory and Methods

## Week 14: Semiparametric methods

**Status:** Ready for Fall 2026  
**Last updated:** August 30, 2026  
**Meetings:** Monday, November 30, and Wednesday, December 2, 2026

<!-- Legacy sources: 11 ORF 524/lectures/ORF524-Fall2025-Lecture-Week-12.tex and
ORF524-Fall2025-Lecture-Week-12_Lecture.pdf; ORF524-Fall2025-Lecture-Week-13.tex;
09 ORF 524/lectures/ORF524-Fall2024-Lecture-Week-12_Notes.pdf and Week-13_Notes.pdf; and
09 ORF 524/assignments/pset8.tex and pset9.tex. -->

## Central question

How do plug in estimators in canonical semiparametric models reveal the general two step structure, and when does estimating an infinite dimensional first step still permit $\sqrt n$ rate inference for a finite dimensional target?

## Learning goals

By the end of the week, you should be able to:

1. identify the finite dimensional target, infinite dimensional nuisance, plug in estimator, and second step in three canonical semiparametric models;
2. derive the residualized identification formula and plug in estimator for the partially linear model;
3. derive direct and integration by parts representations for average derivatives;
4. work from the integrated squared density plug in estimator through its V-statistic, U-statistic, Hoeffding decomposition, bias, variance, and asymptotic linear representation;
5. formulate a general two step semiparametric estimator through a moment condition;
6. state the Fréchet linearization, stochastic equicontinuity, and first step representation conditions that produce the adjustment term $\alpha_0$; and
7. map the canonical models into the general theory and distinguish a zero first step adjustment from a nonzero one.

## In-class route

This is one continuous route through the week, not a division into class meetings. Proceed in order, stop wherever the discussion ends, and resume at the next unfinished stop. The linked sections contain the full student facing exposition; the route identifies the concepts and activities to foreground live.

| Stop | Live focus | Mode |
|---|---|---|
| **W14.1** | [Plug in motivation and canonical semiparametric models](#w14-stop-1) | Discuss + Checkpoint 1 |
| **W14.2** | [Partially linear regression as a two step estimator](#w14-stop-2) | Board work 1 |
| **W14.3** | [Average derivatives by plug in and integration by parts](#w14-stop-3) | Board work 2 + Checkpoint 2 |
| **W14.4** | [Integrated squared density fully worked](#w14-stop-4) | Board work 3 + Board work 4 + AI rate audit + Checkpoint 3 |
| **W14.5** | [General two step semiparametric theory](#w14-stop-5) | Board work 5 + Checkpoint 4 |
| **W14.6** | [The examples inside the general theory](#w14-stop-6) | Synthesis + AI proof map audit + Optional cross fitting + Checkpoint 5 |

Stops W14.1--W14.4 form the model driven arc: plug in estimation first, then three canonical examples, with integrated squared density developed completely. Stops W14.5--W14.6 abstract the common two step structure only after the examples have supplied the objects and proof problems that the general theory must organize.

## How to use this chapter

**Prepare:** Review two step estimation from Week 10, kernel density and derivative estimation from Week 11, and series regression from Week 12. For each of the three models in Section 1, try to identify what is estimated first and what finite dimensional quantity is computed second.

**In class:** Follow the continuous route from concrete plug in estimators to the general theory. The partially linear and average derivative examples establish the recurring two step construction. Integrated squared density is the fully worked bridge: its Hoeffding decomposition makes the first order adjustment and higher order remainder visible before the abstract notation is introduced.

**Review:** For every estimator, write down the target, nuisance function, first step estimate, second step formula, leading linear term, deterministic bias, and higher order remainder. Then identify which of the three general theory conditions controls each part.

**Practice:** The matching assignment complements this chapter with ungraded problems for consolidating the live material and extending it to influence function calculations, missing outcomes, orthogonality, and cross fitting.

**Prerequisites:** Weeks 1--6 and 9--12; conditional expectation; CLT and Slutsky; delta method; M- and Z-estimation; OLS residualization; kernel and series estimators; stochastic orders.

## Full chapter map

1. [Plug in estimation and canonical models](#1-plug-in-estimation-and-canonical-models)
2. [The partially linear model](#2-the-partially-linear-model)
3. [Average derivatives](#3-average-derivatives)
4. [Integrated squared density](#4-integrated-squared-density)
5. [General two step semiparametric theory](#5-general-two-step-semiparametric-theory)
6. [The examples inside the general theory](#6-the-examples-inside-the-general-theory)
7. [Recap](#7-recap)

<a id="w14-stop-1"></a>

## 1. Plug in estimation and canonical models

### 1.1 Functionals of nonparametric objects

A semiparametric problem combines a finite dimensional target with an infinite dimensional nuisance function. In the simplest plug in notation,

$$
\theta_0=T(\gamma_0),
\qquad
\widehat\theta=T(\widehat\gamma).
$$

The first step estimates the function $\gamma_0$. The second step evaluates the target map $T$ at that estimate. A more general estimator chooses $\widehat\theta$ so that

$$
\mathbb{E}_n[m(W,\widehat\theta,\widehat\gamma)]=0,
$$

or approximately minimizes the norm of this sample moment. This is a two step estimator even when the final formula is presented as one compact expression.

The central difficulty is not merely that $\widehat\gamma$ converges more slowly than a parametric estimator. The relevant question is how its error enters the low dimensional target. A slowly converging first step may have a negligible first order effect, while a seemingly accurate first step may still create nonnegligible bias or dependence.

### 1.2 Three canonical models

| Model | Finite dimensional target | Infinite dimensional nuisance | Plug in and second step |
|---|---|---|---|
| Partially linear regression | Regression coefficient $\beta_0$ | Conditional means of $X$ and $Y$ given $Z$ | Estimate both conditional means, form residuals, then regress residualized $Y$ on residualized $X$ |
| Average derivative | $\theta_w=\mathbb E[w(X)\nabla\mu_0(X)]$ | Regression function $\mu_0$ and, after integration by parts, density $f_0$ | Estimate a regression derivative directly or estimate a density derivative and average |
| Integrated squared density | $\theta_0=\int f_0(x)^2\mkern3mu dx$ | Density $f_0$ | Estimate $f_0$, average fitted density values or integrate the squared fit, then analyze the resulting pairwise statistic |

These examples are not illustrations added after the theory. They generate the theory. Partially linear regression shows residualization and a first step whose local effect can vanish. Average derivatives show that changing the representation of the target can change which nuisance derivative must be estimated. Integrated squared density shows explicitly how a plug in estimator produces a first order projection and a higher order pairwise remainder.

### 1.3 The target can be regular when the nuisance is not

The nuisance functions in the table are estimated nonparametrically and generally do not converge at the $\sqrt n$ rate. Nevertheless, the finite dimensional target may admit

$$
\sqrt n(\widehat\theta-\theta_0) =
\frac1{\sqrt n}\sum_{i=1}^n\varphi_0(W_i)
+o_{\mathbb{P}}(1).
$$

The examples will reveal where $\varphi_0$ comes from. Only after that calculation will Section 5 state the abstract conditions that produce it.

### Checkpoint 1

1. In each model, what is estimated first and what is computed second?
2. Why is the convergence rate of the whole nuisance function not by itself the convergence rate of the target?
3. Which example naturally produces a pairwise statistic?
4. Which example uses a change of representation through integration by parts?

<a id="w14-stop-2"></a>

## 2. The partially linear model

### 2.1 Model and residualized identification

Let $W=(Y,X,Z)$, where $X\in\mathbb R^p$, and suppose

$$
Y=X'\beta_0+g_0(Z)+\varepsilon,
\qquad
\mathbb E[\varepsilon\mid X,Z]=0.
$$

Define the conditional means

$$
\mu_{X,0}(Z)=\mathbb E[X\mid Z],
\qquad
\mu_{Y,0}(Z)=\mathbb E[Y\mid Z],
$$

and the residuals

$$
\widetilde X=X-\mu_{X,0}(Z),
\qquad
\widetilde Y=Y-\mu_{Y,0}(Z).
$$

Taking conditional expectations given $Z$ gives

$$
\mu_{Y,0}(Z)=\mu_{X,0}(Z)'\beta_0+g_0(Z).
$$

Subtracting this equation from the model yields

$$
\widetilde Y=\widetilde X'\beta_0+\varepsilon.
$$

If

$$
Q_0=\mathbb E[\widetilde X\widetilde X']
$$

is nonsingular, then

$$
\beta_0 =
Q_0^{-1}\mathbb E[\widetilde X\widetilde Y] =
Q_0^{-1}\mathbb E[\widetilde XY].
$$

The second equality uses $\mathbb E[\widetilde X\mu_{Y,0}(Z)]=0$. This is the population residualization formula developed in the handwritten notes.

### 2.2 The plug in estimator

Estimate the two conditional means and form

$$
\widehat{\widetilde X}_i=X_i-\widehat\mu_X(Z_i),
\qquad
\widehat{\widetilde Y}_i=Y_i-\widehat\mu_Y(Z_i).
$$

The plug in estimator is

$$
\widehat\beta =
\left\lbrace
\mathbb{E}_n[\widehat{\widetilde X}\widehat{\widetilde X}']
\right\rbrace^{-1}
\mathbb{E}_n[\widehat{\widetilde X}\widehat{\widetilde Y}].
$$

The first step estimates $\mu_{X,0}$ and $\mu_{Y,0}$. The second step is an ordinary least squares regression using the two residuals. The handwritten construction records two natural implementations.

| First step | Conditional mean fit | Second step |
|---|---|---|
| Local polynomial | Fit each component of $X$ and the outcome $Y$ locally on $Z$ | Regress the local residual for $Y$ on the local residuals for $X$ |
| Series | Regress each component of $X$ and $Y$ on a common basis $p_k(Z)$ | Regress the two series residuals |

For a common series basis matrix $P$, let $M_P=I-P(P'P)^{-1}P'$. Then

$$
\widehat\beta=(X'M_PX)^{-1}X'M_PY.
$$

This is the Frisch--Waugh--Lovell coefficient from a joint regression of $Y$ on $X$ and the series terms. When $k$ grows, the formula remains exact but fixed dimension OLS theory does not by itself justify inference.

### 2.3 The moment and its leading term

Let $\gamma=(\mu_X,\mu_Y)$ and define

$$
m(W,\beta,\gamma) =
\lbrace X-\mu_X(Z)\rbrace
\left[
Y-\mu_Y(Z)
-\lbrace X-\mu_X(Z)\rbrace'\beta
\right].
$$

At $(\beta_0,\gamma_0)$,

$$
m(W,\beta_0,\gamma_0)=\widetilde X\varepsilon.
$$

Under conditions that make the first step remainder negligible,

$$
\sqrt n(\widehat\beta-\beta_0) =
Q_0^{-1}\frac1{\sqrt n}\sum_{i=1}^n
\widetilde X_i\varepsilon_i
+o_{\mathbb{P}}(1).
$$

This conclusion is not obtained by pretending that the fitted conditional means were known. Section 6 will use the general theory to explain why the linear first step adjustment is zero for this residualized moment.

> [!IMPORTANT]
> **Board work 1 — Partially linear residualization**
>
> 1. Start from $Y=X'\beta_0+g_0(Z)+\varepsilon$ and derive $\mu_{Y,0}(Z)=\mu_{X,0}(Z)'\beta_0+g_0(Z)$.
> 2. Subtract conditional means and derive the residualized regression.
> 3. Prove both population formulas for $\beta_0$.
> 4. Replace the two conditional means by local polynomial or series fits and write the plug in estimator.
> 5. For a common series basis, derive the Frisch--Waugh--Lovell representation and identify why growing dimension inference needs more than the fixed dimension OLS theorem.

<a id="w14-stop-3"></a>

## 3. Average derivatives

### 3.1 Direct plug in estimation

Let $X\in\mathbb R^d$ have density $f_0$, and let

$$
\mu_0(x)=\mathbb E[Y\mid X=x].
$$

For a scalar weight $w$, consider the vector target

$$
\theta_w=\mathbb E[w(X)\nabla\mu_0(X)].
$$

A direct plug in estimator is

$$
\widehat\theta_w^{\mathrm{direct}} =
\mathbb{E}_n[w(X)\nabla\widehat\mu(X)].
$$

A local polynomial estimator supplies derivative coefficients at each evaluation point. A differentiable series fit $\widehat\mu_k(x)=p_k(x)'\widehat\beta_k$ supplies

$$
\nabla\widehat\mu_k(x) =
\nabla p_k(x)'\widehat\beta_k.
$$

This is a genuine two step estimator: first estimate the regression function and its derivative, then average that derivative with the chosen weight.

### 3.2 Density weighting and integration by parts

The handwritten notes specialize to $w=f_0$. For coordinate $j$,

$$
\theta_{f,j} =
\int f_0(x)^2\partial_j\mu_0(x)\mkern3mu dx.
$$

If the boundary term $\mu_0(x)f_0(x)^2$ vanishes on the relevant faces or tails and the derivatives are integrable, then

$$
\theta_{f,j} =
\left[\mu_0(x)f_0(x)^2\right]_{\partial_j\mathcal X} -
\int\mu_0(x)\partial_j\lbrace f_0(x)^2\rbrace\mkern3mu dx =
-2\int\mu_0(x)f_0(x)\partial_j f_0(x)\mkern3mu dx.
$$

Since $\mu_0(X)=\mathbb E[Y\mid X]$,

$$
\theta_f=-2\mathbb E[Y\nabla f_0(X)].
$$

The representation has changed the first step. Instead of estimating $\nabla\mu_0$, one can estimate $\nabla f_0$. The boundary condition is part of the identification argument and cannot be suppressed without comment.

### 3.3 Density derivative plug in estimator

For a differentiable symmetric kernel $K$, use the leave one out KDE

$$
\widehat f_{-i,h}(x)=\frac1{(n-1)h^d}\sum_{j\ne i}K\mkern-3mu\left(\frac{X_j-x}{h}\right).
$$

Differentiating with respect to $x$ and using symmetry gives

$$
\nabla\widehat f_{-i,h}(X_i) =
-\frac1{(n-1)h^{d+1}}
\sum_{j\ne i}
\nabla K\mkern-3mu\left(\frac{X_j-X_i}{h}\right) =
\frac1{(n-1)h^{d+1}}
\sum_{j\ne i}
\nabla K\mkern-3mu\left(\frac{X_i-X_j}{h}\right).
$$

The density weighted estimator is

$$
\widehat\theta_{f,h} =
-\frac2n\sum_{i=1}^n
Y_i\nabla\widehat f_{-i,h}(X_i).
$$

The direct regression derivative estimator and this density derivative estimator target the same parameter under the integration by parts conditions, but their first steps, smoothing biases, and stochastic remainders differ.

> [!IMPORTANT]
> **Board work 2 — Average derivative representations**
>
> 1. Write $\theta_w$ as an ordinary integral.
> 2. Specialize to $w=f_0$ and display the boundary term.
> 3. Perform integration by parts coordinate by coordinate and derive $\theta_f=-2\mathbb E[Y\nabla f_0(X)]$.
> 4. Write the direct regression derivative plug in estimator and the density derivative plug in estimator.
> 5. Identify the first step and second step in each representation.

### Checkpoint 2

1. What boundary condition makes the density weighted identity valid?
2. Why does the factor two appear?
3. Which nuisance derivative is estimated in the direct representation and which is estimated after integration by parts?
4. Why can two estimators of the same target require different remainder conditions?

<a id="w14-stop-4"></a>

## 4. Integrated squared density

### 4.1 Target and two natural plug in estimators

Let $X_1,\ldots,X_n$ be i.i.d. with Lebesgue density $f_0$ on $\mathbb R^d$. The target is

$$
\theta_0 =
\mathbb E[f_0(X)] =
\int f_0(x)^2\mkern3mu dx.
$$

With a symmetric kernel $K$, define

$$
\widehat f_h(x) =
\frac1{nh^d}\sum_{j=1}^n
K\mkern-3mu\left(\frac{x-X_j}{h}\right).
$$

The first plug in estimator averages fitted density values:

$$
\widehat\theta_{\mathrm{PI}} =
\mathbb{E}_n[\widehat f_h(X)] =
\frac1{n^2h^d}
\sum_{i=1}^n\sum_{j=1}^n
K\mkern-3mu\left(\frac{X_i-X_j}{h}\right).
$$

The second integrates the squared density estimate:

$$
\widehat\theta_{\mathrm{ID}} =
\int\widehat f_h(x)^2\mkern3mu dx =
\frac1{n^2h^d}
\sum_{i=1}^n\sum_{j=1}^n
L\mkern-3mu\left(\frac{X_i-X_j}{h}\right),
$$

where

$$
L(v)=\int K(u)K(u+v)\mkern3mu du.
$$

Thus both natural plug in estimators are V-statistics. Their pairwise kernels differ: the first uses $K$, while the second uses the convolution kernel $L$. The remaining derivation applies to either one, so let $G$ denote the relevant symmetric kernel, with $\int G=1$, and define

$$
G_h(u)=h^{-d}G(u/h).
$$

### 4.2 From the V-statistic to the U-statistic

Write

$$
V_{n,h} =
\frac1{n^2}\sum_{i=1}^n\sum_{j=1}^n
G_h(X_i-X_j).
$$

Separate diagonal and off diagonal pairs. Since $G_h(0)=h^{-d}G(0)$,

$$
V_{n,h} =
\frac{G(0)}{nh^d}
+
\frac1{n^2}\sum_{i\ne j}G_h(X_i-X_j).
$$

Define the order two U-statistic

$$
U_{n,h} =
\binom{n}{2}^{-1}
\sum_{1\leq i<j\leq n}
G_h(X_i-X_j).
$$

Symmetry gives the exact relation

$$
V_{n,h} =
\frac{n-1}{n}U_{n,h}
+
\frac{G(0)}{nh^d}.
$$

Removing the diagonal is also a leave one out correction:

$$
U_{n,h} =
\frac1n\sum_{i=1}^n
\widehat f_{-i,h}^{G}(X_i),
\qquad
\widehat f_{-i,h}^{G}(X_i) =
\frac1{n-1}\sum_{j\ne i}G_h(X_i-X_j).
$$

This identity is the algebraic bridge from the plug in estimator to the U-statistic analyzed below.

### 4.3 Hoeffding decomposition

Let

$$
\theta_h =
\mathbb E[G_h(X_1-X_2)] =
\int f_0(x)(G_h*f_0)(x)\mkern3mu dx,
$$

and define the conditional projection

$$
g_h(x) =
\mathbb E[G_h(x-X_2)] =
(G_h*f_0)(x).
$$

The linear projection kernel is

$$
\ell_h(x)=2\lbrace g_h(x)-\theta_h\rbrace.
$$

Define the degenerate kernel

$$
q_h(x_1,x_2) =
G_h(x_1-x_2)
-g_h(x_1)
-g_h(x_2)
+\theta_h.
$$

Then

$$
U_{n,h}-\theta_h =
\frac1n\sum_{i=1}^n\ell_h(X_i)
+
Q_{n,h},
$$

where

$$
Q_{n,h} =
\binom{n}{2}^{-1}
\sum_{1\leq i<j\leq n}
q_h(X_i,X_j).
$$

The kernel is degenerate:

$$
\mathbb E[q_h(X_1,X_2)\mid X_1]=0,
\qquad
\mathbb E[q_h(X_1,X_2)\mid X_2]=0.
$$

Consequently, the linear projection and degenerate remainder are uncorrelated. This is the exact decomposition drawn in the handwritten notes.

> [!IMPORTANT]
> **Board work 3 — Plug in, diagonal removal, and Hoeffding projection**
>
> 1. Derive both plug in V-statistics and identify why their kernels differ.
> 2. Split the generic V-statistic into diagonal and off diagonal terms.
> 3. Derive the exact relation between $V_{n,h}$ and $U_{n,h}$.
> 4. Show that $U_{n,h}$ is the average leave one out fitted density.
> 5. Add and subtract the two conditional projections to derive the Hoeffding decomposition.
> 6. Verify degeneracy of $q_h$ and orthogonality of the linear and degenerate terms.

### 4.4 Bias, variance, and asymptotic linearity

Suppose $G$ has order $p$ and $f_0$ is sufficiently smooth and integrable. Then

$$
\theta_h-\theta_0=O(h^p).
$$

This is deterministic smoothing bias, not a stochastic order. If $g_h\to f_0$ in $L^2(\mathbb{P})$, then

$$
\ell_h(X)
\to
2\lbrace f_0(X)-\theta_0\rbrace
$$

in $L^2(\mathbb{P})$. Thus the linear part has variance

$$
\mathbb V\mkern-3mu\left(
\frac1n\sum_{i=1}^n\ell_h(X_i)
\right) =
\frac1n\mathbb V[\ell_h(X)]
\to
\frac4n\mathbb V[f_0(X)].
$$

For the degenerate part,

$$
\mathbb V[Q_{n,h}] =
\frac{2}{n(n-1)}
\mathbb E[q_h(X_1,X_2)^2] =
O\mkern-3mu\left(\frac1{n^2h^d}\right).
$$

The three pieces therefore give the representative mean squared error expansion

$$
\mathbb E[(U_{n,h}-\theta_0)^2]
\lesssim
h^{2p}
+
\frac1n
+
\frac1{n^2h^d}.
$$

On the $\sqrt n$ scale,

$$
\sqrt n(U_{n,h}-\theta_0) =
\frac1{\sqrt n}\sum_{i=1}^n
2\lbrace f_0(X_i)-\theta_0\rbrace
+
o_{\mathbb{P}}(1)
$$

provided, among other regularity conditions,

$$
h\to0,
\qquad
nh^d\to\infty,
\qquad
\sqrt n h^p\to0.
$$

The first condition makes the projection approach its limit, the second removes the degenerate remainder, and the third removes smoothing bias. If $h=n^{-a}$, the last two conditions require

$$
\frac1{2p}<a<\frac1d,
$$

which is feasible when $p>d/2$. These are transparent sufficient conditions for this estimator, not a characterization of the weakest smoothness under which the functional is root $n$ estimable.

The influence function and asymptotic variance are

$$
\varphi_0(X)=2\lbrace f_0(X)-\theta_0\rbrace,
\qquad
V_0=4\mathbb V[f_0(X)].
$$

A feasible estimate uses

$$
\widehat\varphi_i =
2\lbrace
\widehat f_{-i,h}^{G}(X_i)-U_{n,h}
\rbrace,
\qquad
\widehat V=\mathbb{E}_n[\widehat\varphi\widehat\varphi'].
$$

For a scalar target, the standard error is

$$
\mathrm{se}(U_{n,h}) =
\sqrt{\widehat V/n}.
$$

### 4.5 Why the diagonal matters

For the original V-statistic, diagonal negligibility at the $\sqrt n$ scale requires

$$
\sqrt n h^d\to\infty.
$$

With $h=n^{-a}$, this adds $a<1/(2d)$. Compatibility with $\sqrt n h^p\to0$ then requires $p>d$. The diagonal free U-statistic needs only $p>d/2$ under this representative argument. The exact diagonal term therefore changes the bandwidth region rather than merely changing notation.

> [!IMPORTANT]
> **Board work 4 — Rates and the asymptotic linear representation**
>
> 1. Derive the smoothing bias order for $\theta_h-\theta_0$.
> 2. Derive the variance orders of the linear and degenerate terms separately.
> 3. Assemble the mean squared error expansion.
> 4. Put every term on the $\sqrt n$ scale and derive the bandwidth interval.
> 5. Compare the U-statistic and V-statistic bandwidth requirements.
> 6. Construct the projected influence values, variance estimate, and standard error.

> [!TIP]
> **AI interaction 1 — Audit an integrated squared density proof**
>
> Copy the prompt below into an AI interface and audit the response.

```text
An analyst estimates integrated squared density by averaging a KDE at the same
observations used to construct it. The proof calls the smoothing bias
O_{\mathbb{P}}(h^p), ignores diagonal terms, applies an ordinary iid CLT directly
to the pairwise average, and reports an influence function variance as the
standard error.

Repair the proof from the beginning. Derive the V statistic, remove or account
for its diagonal, give the Hoeffding decomposition of the resulting U statistic,
and state separate conditions for smoothing bias, the linear projection, and
the degenerate remainder. End with the correct standard error.
```

**Audit question:** Does the response recover all three terms in the handwritten asymptotic linear representation rather than merely quote a final limit?

### Checkpoint 3

1. Why are both natural plug in estimators V-statistics?
2. What exact term disappears when the diagonal is removed?
3. Which term in the Hoeffding decomposition becomes the influence function?
4. Why are $nh^d\to\infty$ and $\sqrt n h^p\to0$ logically different conditions?
5. Why can the U-statistic permit a wider bandwidth region than the V-statistic?

<a id="w14-stop-5"></a>

## 5. General two step semiparametric theory

### 5.1 Setup

Let $W_1,\ldots,W_n$ be i.i.d. observations. Let $\beta_0\in\mathcal B\subset\mathbb R^q$ be the finite dimensional target and $\gamma_0\in\Gamma$ an infinite dimensional nuisance. Write

$$
\beta_0
\in
\arg\min_{\beta\in\mathcal B}
\Vert
\mathbb E[m(W,\beta,\gamma_0)]
\Vert,
$$

and define the two step estimator

$$
\widehat\beta
\in
\arg\min_{\beta\in\mathcal B}
\Vert
\mathbb{E}_n[m(W,\beta,\widehat\gamma)]
\Vert.
$$

The derivation below specializes to the exactly identified case: $m$ has the same dimension as $\beta$, the population moment has a unique zero at $\beta_0$, and the sample moment is solved up to a negligible error. The derivative in the finite dimensional parameter is

$$
H_0 =
\left.
\frac{\partial}{\partial\beta'}
\mathbb E[m(W,\beta,\gamma_0)]
\right|_{\beta=\beta_0},
$$

which is assumed nonsingular.

### 5.2 First linearize in the finite dimensional parameter

A Taylor expansion in $\beta$ gives the proof map

$$
\sqrt n(\widehat\beta-\beta_0) =
-\widehat H^{-1}
\frac1{\sqrt n}\sum_{i=1}^n
m(W_i,\beta_0,\widehat\gamma)
+
o_{\mathbb{P}}(1),
$$

with $\widehat H\to_{\mathbb{P}}H_0$. The remaining problem is entirely in

$$
\frac1{\sqrt n}\sum_{i=1}^n
m(W_i,\beta_0,\widehat\gamma).
$$

Add and subtract $m(W_i,\beta_0,\gamma_0)$. The first step contribution must now be linearized in a function space.

### 5.3 Three key conditions for the first step

Let $D(W,\beta_0,\gamma_0)[v]$ denote the derivative in nuisance direction $v$.

**Condition 1: Fréchet linearization.** For $\gamma$ near $\gamma_0$,

$$
\Vert
m(W,\beta_0,\gamma) -
m(W,\beta_0,\gamma_0) -
D(W,\beta_0,\gamma_0)[\gamma-\gamma_0]
\Vert
\leq
b(W)\Vert\gamma-\gamma_0\Vert_\Gamma^2,
$$

where $\mathbb E[b(W)]<\infty$, together with

$$
\sqrt n
\Vert\widehat\gamma-\gamma_0\Vert_\Gamma^2 =
o_{\mathbb{P}}(1).
$$

This condition removes the nonlinear first step remainder after multiplication by $\sqrt n$.

**Condition 2: Stochastic equicontinuity.** The centered empirical linear term satisfies

$$
\sqrt n\left\lbrace
\mathbb{E}_n
D(W,\beta_0,\gamma_0)
[\widehat\gamma-\gamma_0] -
\mathbb E
D(W,\beta_0,\gamma_0)
[\widehat\gamma-\gamma_0]
\right\rbrace =
o_{\mathbb{P}}(1).
$$

This condition controls evaluation of a random fitted nuisance on the same empirical process.

**Condition 3: First step linear representation.** There is a mean zero function $\alpha_0(W)$ such that

$$
\sqrt n\mkern3mu
\mathbb E
D(W,\beta_0,\gamma_0)
[\widehat\gamma-\gamma_0] =
\frac1{\sqrt n}\sum_{i=1}^n\alpha_0(W_i)
+
o_{\mathbb{P}}(1).
$$

This is the adjustment generated by estimating the nuisance. It may be nonzero, as in integrated squared density or inverse weighting, or zero because the population moment is locally insensitive to the nuisance, as in the residualized partially linear score.

### 5.4 General asymptotic linear representation

Combining the three conditions gives

$$
\sqrt n(\widehat\beta-\beta_0) =
-\frac1{\sqrt n}\sum_{i=1}^n
H_0^{-1}
\lbrace
m(W_i,\beta_0,\gamma_0)
+
\alpha_0(W_i)
\rbrace
+
o_{\mathbb{P}}(1).
$$

Thus the influence function is

$$
\varphi_0(W) =
-H_0^{-1}
\lbrace
m(W,\beta_0,\gamma_0)
+
\alpha_0(W)
\rbrace.
$$

If $\mathbb E[\varphi_0\varphi_0']$ is finite and nonsingular in the relevant directions, a multivariate CLT supplies the limiting distribution. A consistent estimate of that covariance matrix must account for both the original moment and the first step adjustment.

> [!IMPORTANT]
> **Board work 5 — General two step linearization**
>
> 1. Linearize the sample moment in $\beta$ and isolate $m(W_i,\beta_0,\widehat\gamma)$.
> 2. Add and subtract $m(W_i,\beta_0,\gamma_0)$.
> 3. Apply the Fréchet derivative and bound the quadratic remainder.
> 4. Split the derivative term into its centered empirical part and population mean.
> 5. Apply stochastic equicontinuity and the first step linear representation.
> 6. Derive the adjustment $\alpha_0$ and the final influence function.

### Checkpoint 4

1. Which expansion concerns $\beta$ and which concerns $\gamma$?
2. What does the quadratic rate condition remove?
3. What dependence problem is handled by stochastic equicontinuity?
4. Where does $\alpha_0$ enter, and when can it be zero?
5. Why would treating $\widehat\gamma$ as known generally omit part of the variance?

<a id="w14-stop-6"></a>

## 6. The examples inside the general theory

### 6.1 Three examples in common notation

| Example | Moment or functional | First step effect |
|---|---|---|
| Integrated squared density | $m(X,\theta,f)=f(X)-\theta$ | Nonzero; estimating $f$ contributes a second copy of $f_0(X)-\theta_0$ |
| Partially linear regression | Residualized moment using $\gamma=(\mu_X,\mu_Y)$ | Zero at first order for the symmetric residualized moment |
| Average derivatives | $m(W,\theta,\mu)=w(X)\nabla\mu(X)-\theta$, or the density derivative representation when $w=f_0$ | Generally nonzero in either raw plug in representation, but the adjustment depends on whether the first step estimates $\nabla\mu_0$ or $\nabla f_0$ |

For integrated squared density,

$$
m(X,\theta_0,f_0)=f_0(X)-\theta_0.
$$

The Hoeffding calculation showed that the first step contributes

$$
\alpha_0(X)=f_0(X)-\theta_0.
$$

Since $H_0=-1$, the general influence function becomes

$$
\varphi_0(X) =
2\lbrace f_0(X)-\theta_0\rbrace,
$$

exactly as obtained from the first projection of the U-statistic.

For the partially linear moment, perturb $\mu_{X,0}$ by $a(Z)$ and $\mu_{Y,0}$ by $b(Z)$. Its directional derivative at the truth is

$$
-a(Z)\varepsilon
+
\widetilde X
\lbrace
a(Z)'\beta_0-b(Z)
\rbrace.
$$

Its expectation is zero because $\mathbb E[\varepsilon\mid Z]=0$ and $\mathbb E[\widetilde X\mid Z]=0$. Thus $\alpha_0=0$ for this moment: the first step still needs to be accurate enough for the nonlinear remainder to vanish, but it does not contribute a separate linear adjustment.

For average derivatives, the direct moment $m(W,\theta,\mu)=w(X)\nabla\mu(X)-\theta$ changes linearly when $\mu$ is perturbed, so its first step adjustment is generally nonzero. After integration by parts with $w=f_0$, the density derivative moment $m(W,\theta,f)=-2Y\nabla f(X)-\theta$ also has a generally nonzero first step adjustment. The two adjustments need not have the same form because the nuisance and its estimated derivative have changed. Neither raw plug in moment is Neyman orthogonal without an additional construction.

As an optional fourth comparison, let $\gamma_0(x)$ be a positive conditional mean used as an inverse weight and consider

$$
m(W,\beta,\gamma) =
\frac{a(W,\beta)}{\gamma(X)}.
$$

The reciprocal expansion gives

$$
D(W,\beta_0,\gamma_0)[v] =
-\frac{a(W,\beta_0)}{\gamma_0(X)^2}v(X).
$$

This derivative is generally nonzero, so the asymptotic linear representation of $\widehat\gamma$ generates a nonzero $\alpha_0$. This is the inverse weighting example developed in the second handwritten record; it does not require adding missing outcome methods to the live chapter.

### 6.2 Reading the integrated squared density proof through the three conditions

| Concrete term | General theory role |
|---|---|
| $\theta_h-\theta_0$ | Deterministic approximation or smoothing bias |
| $n^{-1}\sum_i\ell_h(X_i)$ | Original moment plus first step linear adjustment |
| $Q_{n,h}$ | Higher order remainder generated by using the same estimated density in the second step |
| $nh^d\to\infty$ | Condition making the higher order pairwise remainder negligible |
| $\sqrt n h^p\to0$ | Condition making smoothing bias negligible |
| $2\lbrace f_0(X)-\theta_0\rbrace$ | Final influence function $-H_0^{-1}(m+\alpha_0)$ |

The abstract notation is therefore a compression of the worked example, not a replacement for it.

> [!TIP]
> **AI interaction 2 — Audit a general two step proof**
>
> Copy the prompt below into an AI interface and audit the response.

```text
A proof Taylor expands a two step estimator only in the finite dimensional
parameter beta. It then replaces the fitted nuisance gamma_hat by gamma_0
because gamma_hat is consistent and applies a central limit theorem to the
original moment.

Identify the missing steps. State a Frechet linearization in the nuisance,
separate the centered empirical derivative term from its population mean,
explain the roles of stochastic equicontinuity and the first step adjustment
alpha, and show how the final influence function combines the original moment
with alpha. Use integrated squared density and partially linear regression to
illustrate nonzero and zero adjustments.
```

**Audit question:** Does the response explain why consistency alone neither removes the first step adjustment nor controls the higher order remainder?

### 6.3 Optional modern extension: orthogonality and cross fitting

The partially linear moment is now often called orthogonal because its population derivative in every nuisance direction is zero. Orthogonality sets $\alpha_0=0$, but the quadratic remainder still needs a product or squared rate condition.

Cross fitting trains nuisance estimates outside the observations used to evaluate the moment. It can simplify control of the centered empirical derivative term and reduce overfitting dependence. It does not create orthogonality, guarantee nuisance convergence rates, remove deterministic bias, or repair identification. These ideas extend the slide theory but are not part of the main live derivation.

### Checkpoint 5

1. Why does integrated squared density have a nonzero $\alpha_0$?
2. Why is $\alpha_0=0$ for the residualized partially linear moment?
3. Why is the first step adjustment generally nonzero for both raw average derivative representations even though the adjustments differ?
4. What feature of inverse weighting makes its first step derivative nonzero?
5. Which pieces of the Hoeffding decomposition correspond to the general theory remainder terms?
6. What can cross fitting simplify, and what does it leave unchanged?

## 7. Recap

- Semiparametric estimation begins with concrete plug in and two step constructions.
- Partially linear regression estimates two conditional means before a residual on residual regression.
- Average derivatives can be estimated directly through a regression derivative or indirectly through a density derivative after integration by parts.
- Integrated squared density turns a density plug in estimator into a V-statistic, a diagonal free U-statistic, and a Hoeffding decomposition.
- The integrated squared density influence function is the limit of the linear Hoeffding projection.
- Smoothing bias, the linear projection, and the degenerate remainder require separate conditions.
- The general theory first linearizes in the finite dimensional parameter and then in the infinite dimensional nuisance.
- Fréchet linearization, stochastic equicontinuity, and a first step linear representation produce the adjustment term $\alpha_0$.
- The first step adjustment is generally nonzero for integrated squared density, the raw average derivative plug in representations, and inverse weighting, but zero for the residualized partially linear moment.
- Orthogonality and cross fitting are useful extensions, not substitutes for the general two step proof.

## Notation introduced this week

- $\beta_0$ or $\theta_0$: finite dimensional target.
- $\gamma_0$: infinite dimensional nuisance.
- $\mu_{X,0}$ and $\mu_{Y,0}$: conditional means in the partially linear model.
- $\widetilde X$ and $\widetilde Y$: residualized variables.
- $Q_0=\mathbb E[\widetilde X\widetilde X']$: residual Gram matrix in the partially linear model.
- $\theta_w$: weighted average derivative.
- $V_{n,h}$ and $U_{n,h}$: V-statistic and U-statistic versions of an integrated squared density plug in estimator.
- $\ell_h$: linear Hoeffding projection; $q_h$ and $Q_{n,h}$: degenerate kernel and remainder.
- $D[\mkern3mu\cdot\mkern3mu]$: Fréchet derivative in a nuisance direction.
- $H_0$: derivative of the population moment in the finite dimensional parameter.
- $\alpha_0$: first step linear adjustment.
- $\varphi_0=-H_0^{-1}(m+\alpha_0)$: influence function of the final estimator.

## References

- Peter J. Bickel and Kjell A. Doksum, *Mathematical Statistics: Basic Ideas and Selected Topics*, vol. II, CRC Press, 2016.
- Whitney Newey and Daniel McFadden, “Large Sample Estimation and Hypothesis Testing,” in *Handbook of Econometrics*, vol. 4, 1994.
- Aad van der Vaart, *Asymptotic Statistics*, Cambridge University Press, 1998.
