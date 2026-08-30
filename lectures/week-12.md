# ORF 524: Statistical Theory and Methods

## Week 12: Series based nonparametric methods

**Status:** Ready for Fall 2026  
**Last updated:** August 30, 2026  
**Meetings:** Monday, November 16, and Wednesday, November 18, 2026

<!-- Legacy sources: 11 ORF 524/lectures/ORF524-Fall2025-Lecture-Week-11.tex and
ORF524-Fall2025-Lecture-Week-11_Lecture.pdf; 09 ORF 524/lectures/
ORF524-Fall2024-Lecture-Week-11_Notes.pdf; 09 ORF 524/assignments/pset8.tex and pset9.tex; and
earlier series/partitioning exercises. -->

## Central question

How does a growing collection of basis functions approximate an unknown regression function, and how do least squares logic, approximation bias, and out of sample risk govern inference and dimension choice?

## Learning goals

By the end of the week, you should be able to:

1. describe series estimation through nested approximation spaces and concrete polynomial, spline, trigonometric, and partition bases;
2. define the population series projection and its least squares estimator;
3. decompose error into approximation error and coefficient estimation error, both algebraically and geometrically;
4. derive integrated and usable supremum norm rate bounds;
5. derive the fixed $k$ OLS limit, construct feasible pointwise sandwich standard errors, and separate the CLT from approximation bias;
6. compare series undersmoothing with RBC based on a richer auxiliary approximation;
7. distinguish training, validation, test, and leave one out prediction error; and
8. derive the linear smoother shortcut for leave one out cross validation.

## In-class route

This is one continuous route through the week, not a division into class meetings. Proceed in order, stop wherever the discussion ends, and resume at the next unfinished stop. The linked sections contain the full student facing exposition; the route identifies the concepts and activities to foreground live.

| Stop | Live focus | Mode |
|---|---|---|
| **W12.1** | [Series projection and least squares](#w12-stop-1) | Discuss + basis map + Board work 1 + Checkpoint 1 |
| **W12.2** | [Approximation and estimation error](#w12-stop-2) | Board work 2 |
| **W12.3** | [Rates, dimension choice, and the curse of dimensionality](#w12-stop-3) | Discuss + Checkpoint 2 |
| **W12.4** | [Pointwise inference and approximation bias](#w12-stop-4) | Discuss + CLT and bias + feasible OLS inference + RBC + AI inference audit |
| **W12.5** | [Prediction risk and leave one out cross validation](#w12-stop-5) | Board work 3 + AI objective audit |
| **W12.6** | [Kernel and series methods compared](#w12-stop-6) | Synthesis + Checkpoint 3 |

Stops W12.1--W12.3 should repeatedly distinguish the finite $k$ projection from the underlying regression function. Stop W12.4 begins with fixed $k$ OLS inference, then separates the CLT from approximation bias before developing feasible and RBC inference. Stop W12.5 treats cross validation as a tuning method for out of sample risk, not as a guarantee of valid inference or structural model selection.

## How to use this chapter

**Prepare:** Review Week 10 OLS and Week 11 bias and variance tradeoffs. Write indicator bin, polynomial, and spline bases as vectors and attempt Checkpoint 1.

**In class:** Follow the continuous route above through nested approximation spaces, projection, estimation, rates, pointwise OLS inference, approximation bias, RBC, basis choice, prediction risk, cross validation, the leave one out formula, and comparison with kernel methods. Prepared responses supply an alternative for every live prompt if AI is unavailable.

**Review:** For every series claim, state whether the target is the finite $k$ projection $\mu_k$ or the underlying function $\mu$, and identify the basis normalization, approximation norm, growth of $k$, inferential bias condition, variance estimator, and whether bias correction changes the studentization.

**Practice:** The matching assignment complements this chapter with ungraded problems for consolidating the material, extending selected ideas, and preparing for cumulative review.

**Prerequisites:** Weeks 1--6 and 9--11; linear projection; OLS and sandwich variance; stochastic orders; pointwise asymptotic inference; nonparametric bias and variance tradeoffs.

## Full chapter map

1. [Series projection and least squares](#1-series-projection-and-least-squares)
2. [Approximation and estimation error](#2-approximation-and-estimation-error)
3. [Rates, dimension choice, and the curse of dimensionality](#3-rates-dimension-choice-and-the-curse-of-dimensionality)
4. [Pointwise inference](#4-pointwise-inference)
5. [Cross validation and out of sample risk](#5-cross-validation-and-out-of-sample-risk)
6. [Kernel and series methods compared](#6-kernel-and-series-methods-compared)
7. [Recap](#7-recap)

<a id="w12-stop-1"></a>

## 1. Series projection and least squares

### 1.1 Approximation spaces and the regression target

Let $(Y_i,X_i)$ be i.i.d. with

$$
Y_i=\mu(X_i)+\varepsilon_i,
\qquad
\mu(x)=\mathbb E[Y_i\mid X_i=x],
\qquad
\mathbb E[\varepsilon_i\mid X_i]=0.
$$

The conditional mean has a population least squares characterization over all square integrable functions:

$$
\mu
\in
\arg\min_{g\in L^2(\mathbb{P}_X)}
\mathbb E[\lbrace Y-g(X)\rbrace^2].
$$

Series methods replace this unrestricted problem with a sequence of finite dimensional approximation spaces

$$
\mathcal G_1\subset\mathcal G_2\subset\cdots,
\qquad
\mathcal G_k=\lbrace x\mapsto p_k(x)'\beta:\beta\in\mathbb R^k\rbrace.
$$

The dimension $k$ is a smoothing parameter. A small space is stable but may be too rigid; a large space can approximate more functions but requires estimating more coefficients and may produce high leverage.

Let

$$
p_k(x)=(p_{1k}(x),\ldots,p_{kk}(x))'
$$

be a $k$ dimensional basis. The basis is a coordinate system for $\mathcal G_k$; changing coordinates without changing the span does not change the population projection or the fitted values.

### 1.2 Concrete bases

| Basis | Representative functions | Geometry and use |
|---|---|---|
| Algebraic polynomial | $1,x,x^2,\ldots$ | Global support; simple, but raw powers can be poorly conditioned. |
| Orthogonal polynomial | Legendre or another basis orthogonal under a reference measure | Same polynomial span with better numerical and theoretical normalization. |
| Spline | Piecewise polynomials joined at knots | Local support and flexible smoothness control through degree and knots. |
| Trigonometric | $1,\cos(2\pi x),\sin(2\pi x),\ldots$ | Global support and natural periodic approximation. |
| Haar or partition | Indicators of cells, or multiscale differences of cell indicators | Local support; a fixed partition regression estimator is series least squares. |

For example, let $I_{1J},\ldots,I_{JJ}$ partition the support of a scalar $X$ and take $p_J(x)=(1\lbrace x\in I_{1J}\rbrace,\ldots,1\lbrace x\in I_{JJ}\rbrace)'$. The fitted value within each occupied cell is the sample mean of $Y_i$ in that cell. Thus a familiar partition estimator is exactly least squares on indicator regressors.

The handwritten construction also suggests an important comparison. A fixed polynomial, spline, or partition space is selected before seeing the outcomes. A regression tree instead chooses a partition from the data. Both regularize by restricting the function class, but data adaptive selection creates an additional stochastic problem and is outside this week's core theory.

### 1.3 Population projection

The **population series projection** is the best approximation to $\mu$ in $\mathcal G_k$:

$$
\beta_k
\in
\arg\min_{\beta\in\mathbb R^k}
\mathbb E[\lbrace \mu(X)-p_k(X)'\beta\rbrace^2].
$$

Write

$$
Q_k=\mathbb E[p_k(X)p_k(X)'].
$$

If $Q_k$ is nonsingular, then

$$
\beta_k
=Q_k^{-1}\mathbb E[p_k(X)\mu(X)]
=Q_k^{-1}\mathbb E[p_k(X)Y].
$$

Define

$$
\mu_k(x)=p_k(x)'\beta_k,
\qquad
r_k(x)=\mu(x)-\mu_k(x).
$$

The projection normal equation is

$$
\mathbb E[p_k(X)r_k(X)]=0.
$$

This is an unconditional projection statement. It does not say that $\mathbb E[r_k(X)\mid X]=0$, since $r_k(X)$ is already a deterministic function of $X$.

Projection orthogonality gives a Pythagorean identity. For every $g_k\in\mathcal G_k$,

$$
\Vert \mu-g_k\Vert_{L^2(\mathbb{P}_X)}^2 =
\Vert \mu-\mu_k\Vert_{L^2(\mathbb{P}_X)}^2
+
\Vert \mu_k-g_k\Vert_{L^2(\mathbb{P}_X)}^2.
$$

Consequently, enlarging a nested space cannot increase population approximation error, although it can increase sample estimation error.

> [!IMPORTANT]
> **Board work 1 — Population series projection**
>
> 1. Begin with the unrestricted least squares characterization of $\mu$ and restrict the minimization to $\mathcal G_k$.
> 2. Derive $\mathbb E[p_k(X)r_k(X)]=0$ and the Pythagorean identity.
> 3. Use cell indicators to show that partition regression returns within cell sample means.
> 4. Locate polynomial, spline, Haar, and fixed partition estimators in the nested space picture, and explain what changes when a partition is selected from the data.

### 1.4 Sample least squares

The series least squares estimator is

$$
\widehat\beta_k
\in
\arg\min_{\beta\in\mathbb R^k}
\mathbb{E}_n\lbrace Y-p_k(X)'\beta\rbrace^2,
$$

with

$$
\widehat Q_k=\mathbb{E}_n[p_k(X)p_k(X)'].
$$

When $\widehat Q_k$ is nonsingular,

$$
\widehat\beta_k
=\widehat Q_k^{-1}\mathbb{E}_n[p_k(X)Y],
\qquad
\widehat\mu_k(x)=p_k(x)'\widehat\beta_k.
$$

The estimator is ordinary least squares with a dimension $k=k_n$ that may grow with $n$. Treating it as fixed dimension regression can be a useful heuristic, but growing dimension probability bounds must justify the approximation.

### Checkpoint 1

1. Why are both $\mu_k$ and $\mu$ relevant targets?
2. Which equation defines the population projection?
3. Why does a fixed $k$ generally leave nonvanishing approximation bias?
4. How is a partition estimator a series estimator, and why is a data selected partition a harder object?

<a id="w12-stop-2"></a>

## 2. Approximation and estimation error

### 2.1 Exact decomposition

Let

$$
u_{ik}=Y_i-p_k(X_i)'\beta_k
=\varepsilon_i+r_k(X_i).
$$

The sample normal equation gives

$$
\widehat\beta_k-\beta_k =
\widehat Q_k^{-1}\mathbb{E}_n[p_k(X)u_k].
$$

Therefore

$$
\widehat\mu_k(x)-\mu(x) =
p_k(x)'(\widehat\beta_k-\beta_k)-r_k(x).
$$

This exact decomposition separates:

- **estimation error**, caused by fitting $k$ coefficients from a finite sample; and
- **approximation error**, caused by representing $\mu$ in a $k$ dimensional space.

Increasing $k$ reduces approximation error but increases estimation variability and can destabilize the empirical Gram matrix.

Because $\widehat\mu_k-\mu_k\in\mathcal G_k$ and $r_k$ is orthogonal to $\mathcal G_k$, the function decomposition also gives the exact population norm identity

$$
\Vert\widehat\mu_k-\mu\Vert_{L^2(\mathbb{P}_X)}^2 =
\Vert\widehat\mu_k-\mu_k\Vert_{L^2(\mathbb{P}_X)}^2
+
\Vert r_k\Vert_{L^2(\mathbb{P}_X)}^2.
$$

This identity is exact even though $\widehat\mu_k$ is random: conditional on the fitted coefficients, $\widehat\mu_k-\mu_k$ remains an element of $\mathcal G_k$. It is the geometric source of the estimation plus approximation decomposition used in rate calculations.

### 2.2 Basis normalization

Series bases can be reparameterized without changing their span. It is therefore useful to normalize the basis so that eigenvalues of

$$
Q_k=\mathbb E[p_k(X)p_k(X)']
$$

remain bounded above and away from zero. The sample analogue must also be well behaved. One representative requirement is

$$
\Vert\widehat Q_k-Q_k\Vert_{\mathrm{op}}=o_{\mathbb{P}}(1).
$$

This is a growing dimension uniform law, not an automatic consequence of a fixed dimension LLN.

Put

$$
\zeta_k=\sup_{x\in\mathcal X}\Vert p_k(x)\Vert.
$$

The size of $\zeta_k$ depends on the basis and enters pointwise and supremum norm control.

> [!IMPORTANT]
> **Board work 2 — Series estimation error decomposition**
>
> Starting from the sample normal equation:
>
> 1. derive the exact coefficient and function error decompositions;
> 2. prove the population norm identity and identify where projection orthogonality is used;
> 3. explain why a sample Gram matrix result must be uniform in the growing dimension; and
> 4. compare indicator bin and globally supported polynomial bases.

<a id="w12-stop-3"></a>

## 3. Rates, dimension choice, and the curse of dimensionality

### 3.1 Integrated rate and its source

Let

$$
c_k=\Vert r_k\Vert_{L^2(\mathbb{P}_X)}.
$$

Under bounded Gram matrix eigenvalues, suitable moments, and a growth condition on $k$,

$$
\Vert\widehat\beta_k-\beta_k\Vert =
O_{\mathbb{P}}\mkern-3mu\left(\sqrt{\frac{k}{n}}\right).
$$

The dimension factor is not merely a slogan. If $Q_k$ is normalized and $\mathbb E[u_k^2\Vert p_k(X)\Vert^2]\lesssim k$, then

$$
\mathbb E\left[
\left\Vert
\mathbb{E}_n[p_k(X)u_k]
\right\Vert^2
\right] =
\frac1n
\mathbb E[u_k^2\Vert p_k(X)\Vert^2]
\lesssim
\frac{k}{n}.
$$

Combining this score bound with stability of $\widehat Q_k^{-1}$ produces the coefficient rate. This calculation also shows why controlling the Gram matrix and the score are separate tasks.

Consequently,

$$
\Vert\widehat\mu_k-\mu\Vert_{L^2(\mathbb{P}_X)} =
O_{\mathbb{P}}\mkern-3mu\left(
\sqrt{\frac{k}{n}}+c_k
\right),
$$

and the squared integrated error has order

$$
O_{\mathbb{P}}\mkern-3mu\left(
\frac{k}{n}+c_k^2
\right).
$$

Suppose $c_k\lesssim k^{-\alpha}$ for some approximation exponent $\alpha>0$. Balancing terms gives

$$
k_{\mathrm{MSE}}\asymp n^{1/(2\alpha+1)}
$$

and squared integrated error order

$$
n^{-2\alpha/(2\alpha+1)}.
$$

The exponent $\alpha$ encodes both smoothness and how well the chosen basis approximates the function class.

For a function with smoothness $s$ in $d$ dimensions, a common approximation relation is $c_k\lesssim k^{-s/d}$. Then

$$
k_{\mathrm{MSE}}\asymp n^{d/(2s+d)},
\qquad
\Vert\widehat\mu_k-\mu\Vert_{L^2(\mathbb{P}_X)}^2 =
O_{\mathbb{P}}\mkern-3mu\left(n^{-2s/(2s+d)}\right).
$$

The denominator $2s+d$ makes the dimensional cost explicit: holding smoothness fixed, approximation becomes harder as $d$ grows.

### 3.2 A usable supremum norm bound

Let

$$
c_{\infty,k} =
\sup_{x\in\mathcal X}|r_k(x)|.
$$

The deterministic inequality

$$
\sup_{x\in\mathcal X}
|p_k(x)'(\widehat\beta_k-\beta_k)|
\leq
\zeta_k\Vert\widehat\beta_k-\beta_k\Vert
$$

gives the representative bound

$$
\sup_{x\in\mathcal X}
|\widehat\mu_k(x)-\mu(x)| =
O_{\mathbb{P}}\mkern-3mu\left(
\zeta_k\sqrt{\frac{k}{n}}+c_{\infty,k}
\right).
$$

This bound is deliberately transparent and may be conservative. Sharper basis specific results may replace $\zeta_k\sqrt{k/n}$ by a smaller stochastic term. The essential correction is that a supremum norm error rate is on the square root scale; it is not $\zeta_k\lbrace k/n+k^{-2\alpha}\rbrace$ without a square root or further justification.

### 3.3 Dimension and tensor products

For a $d$ dimensional regressor, tensor product bases can require many terms to achieve a given approximation resolution. If smoothness is held fixed, the effective approximation exponent often deteriorates with $d$. This is the series version of the curse of dimensionality.

Partitioning illustrates the tradeoff. With $J$ approximately equal probability bins in one dimension and a piecewise constant regression estimate, a typical integrated MSE has order

$$
\frac{J}{n}+\frac1{J^2}
$$

for a once differentiable regression function. Balancing gives

$$
J\asymp n^{1/3},
\qquad
\mathrm{IMSE}\asymp n^{-2/3}.
$$

Partitioning here is a fixed basis series method. Adaptive decision tree algorithms are not part of this week's core. The comparison is nevertheless useful: fixed cells make the approximation space explicit, whereas an outcome selected tree makes the space random and adds selection error to the approximation and coefficient estimation tradeoff.

### Checkpoint 2

1. Which term grows with $k$ and which falls with $k$?
2. Why must the stochastic term in a norm rate be on a square root scale?
3. What does $\zeta_k$ measure?
4. Why does a tensor product basis reveal the curse of dimensionality?

<a id="w12-stop-4"></a>

## 4. Pointwise inference

### 4.1 Fixed $k$ OLS inference for the projection

Recall that

$$
u_k=Y-p_k(X)'\beta_k
=\varepsilon+r_k(X),
\qquad
\mathbb E[p_k(X)u_k]=0.
$$

The residual $u_k$ need not have conditional mean zero because it contains approximation error, but its population OLS normal equation is exactly zero. For fixed $k$, the sample normal equation gives

$$
\sqrt n(\widehat\beta_k-\beta_k) =
\widehat Q_k^{-1}
\frac1{\sqrt n}
\sum_{i=1}^n
p_k(X_i)u_{ik}.
$$

If $Q_k$ is nonsingular, the score has finite second moments, and the usual fixed dimension CLT conditions hold, then

$$
\sqrt n(\widehat\beta_k-\beta_k)
\rightsquigarrow
\mathcal N\lbrace
0,Q_k^{-1}\Omega_kQ_k^{-1}
\rbrace,
$$

where

$$
\Omega_k =
\mathbb E[p_k(X)p_k(X)'u_k^2].
$$

For one fixed evaluation point $x$, define

$$
V_k(x) =
p_k(x)'Q_k^{-1}\Omega_kQ_k^{-1}p_k(x).
$$

Provided $V_k(x)>0$,

$$
\frac{\widehat\mu_k(x)-\mu_k(x)}
{\sqrt{V_k(x)/n}}
\rightsquigarrow
\mathcal N(0,1).
$$

This is ordinary heteroskedasticity robust OLS inference for the population series projection. With fixed $k$, it is an inference statement about $\mu_k(x)$, not automatically about $\mu(x)$.

### 4.2 CLT plus approximation bias

The exact function error decomposition gives

$$
\frac{\widehat\mu_k(x)-\mu(x)}
{\sqrt{V_k(x)/n}} =
\underbrace{
\frac{p_k(x)'(\widehat\beta_k-\beta_k)}
{\sqrt{V_k(x)/n}}
}_{\text{OLS random component}} -
\underbrace{
\frac{r_k(x)}{\sqrt{V_k(x)/n}}
}_{\text{approximation bias in standard error units}}.
$$

This is the series version of the KDE and local polynomial CLT plus bias decomposition. If

$$
-\frac{r_k(x)}{\sqrt{V_k(x)/n}}
\to
\delta(x),
$$

then the statistic centered at $\mu(x)$ converges to $\mathcal N\lbrace\delta(x),1\rbrace$. Conventional inference requires

$$
\frac{|r_k(x)|}{\sqrt{V_k(x)/n}}
\to0.
$$

For fixed $k$, $r_k(x)$ generally remains fixed while the standard error shrinks, so the OLS interval is normally an interval for $\mu_k(x)$. When $k\to\infty$, choosing a larger space can reduce $|r_k(x)|$, but the variance and leverage also rise. Series undersmoothing means choosing $k$ large enough that approximation bias is negligible relative to the standard error. A dimension selected for prediction or integrated MSE need not satisfy this pointwise condition.

### 4.3 Feasible sandwich inference

Let

$$
\widehat u_{ik}=Y_i-p_k(X_i)'\widehat\beta_k,
\qquad
\widehat\Omega_k =
\mathbb{E}_n[p_k(X)p_k(X)'\widehat u_k^2].
$$

The OLS sandwich estimator is

$$
\widehat V_k(x) =
p_k(x)'\widehat Q_k^{-1}
\widehat\Omega_k
\widehat Q_k^{-1}p_k(x),
$$

and the feasible standard error is

$$
\widehat s_k(x)=
\sqrt{\widehat V_k(x)/n}.
$$

For fixed $k$, this is the familiar heteroskedasticity robust OLS standard error. For growing $k$, consistency requires additional Gram matrix, moment, and maximum leverage conditions. Residual leverage adjustments may improve finite sample performance, but they do not remove approximation bias.

The conventional interval

$$
\left[
\widehat\mu_k(x)
\pm
z_{1-\alpha/2}\widehat s_k(x)
\right]
$$

is therefore feasible for $\mu_k(x)$ under fixed $k$ OLS conditions. It is feasible for $\mu(x)$ only after the growing dimension CLT and the approximation bias condition have both been established.

| Inferential object | Centering target | Main extra issue |
|---|---|---|
| Fixed $k$ projection | $\mu_k(x)$ | Ordinary finite dimension OLS conditions and sandwich consistency can suffice. |
| Growing $k$ projection | $\mu_k(x)$ | Gram matrix stability, leverage control, residual variance consistency, and a triangular array CLT are required. |
| Underlying function | $\mu(x)$ | All growing $k$ conditions plus approximation bias negligible relative to the standard error. |

### 4.4 RBC inference in practice

Robust bias correction (RBC) is available when the chosen basis admits a leading approximation bias expansion that can be estimated by a richer auxiliary approximation. Suppose

$$
\mu_k(x)-\mu(x)=B_k(x)+R_k(x),
$$

where $B_k(x)$ is the leading bias and $R_k(x)$ is a smaller remainder. Use a higher order or richer series fit to construct

$$
\widehat B_k(x)=
\sum_{i=1}^n b_{i,k}(x)Y_i.
$$

The robust bias corrected center is

$$
\widehat\mu_{\mathrm{bc},k}(x) =
\widehat\mu_k(x)-\widehat B_k(x).
$$

Since series least squares is linear in $Y$, let $P$ be the $n\times k$ matrix with rows $p_k(X_i)'$ and write

$$
\widehat\mu_k(x)=
\sum_{i=1}^n a_{i,k}(x)Y_i,
\qquad
a_{i,k}(x)=
p_k(x)'(P'P)^{-1}p_k(X_i).
$$

Then

$$
\widehat\mu_{\mathrm{bc},k}(x)=
\sum_{i=1}^n
\lbrace a_{i,k}(x)-b_{i,k}(x)\rbrace Y_i.
$$

Let $\widehat u_{i,+}$ be residuals from the richer auxiliary fit. The RBC sandwich variance is

$$
\widehat s_{\mathrm{rbc},k}^2(x) =
\sum_{i=1}^n
\lbrace a_{i,k}(x)-b_{i,k}(x)\rbrace^2
\widehat u_{i,+}^2.
$$

This combined weight formula automatically includes the variance of the estimated bias and its covariance with the original series estimate. Using $\widehat s_k(x)$ after subtracting $\widehat B_k(x)$ is generally invalid.

For partitioning or piecewise polynomial series, a concrete implementation fits one additional polynomial degree on the same partition, or a related richer partition space, and uses that auxiliary fit to estimate the leading Taylor approximation term. Other bases require their own verified approximation expansion. A generic increase in $k$ is not automatically an RBC construction.

If the remaining bias $R_k(x)$ is negligible relative to $\widehat s_{\mathrm{rbc},k}(x)$, the main and auxiliary Gram matrices are stable, no combined weight dominates, and the RBC variance estimator is consistent, then

$$
\frac{\widehat\mu_{\mathrm{bc},k}(x)-\mu(x)}
{\widehat s_{\mathrm{rbc},k}(x)}
\rightsquigarrow
\mathcal N(0,1).
$$

The practical interval is

$$
\left[
\widehat\mu_{\mathrm{bc},k}(x)
\pm
z_{1-\alpha/2}\widehat s_{\mathrm{rbc},k}(x)
\right].
$$

RBC can permit an estimation oriented series dimension while correcting its leading approximation bias, but only when that leading bias is estimable and the remaining bias is small. It is not a generic consequence of fitting more basis terms.

> [!TIP]
> **AI interaction 1 — Fixed regression dressed as nonparametrics**
>
> Copy the prompt below into an AI interface and audit the response.

```text
A proposed proof treats k as fixed, applies the ordinary OLS CLT, and reports
a confidence interval for the unknown regression function mu(x). It then
subtracts a bias estimate from a richer series fit but keeps the original
OLS standard error.

Repair the argument. First derive fixed k sandwich inference for the projection
mu_k(x). Then decompose inference for mu(x) into the OLS random component and
approximation bias. State the growing dimension conditions that must supplement
the fixed k proof. Finally, describe an RBC center based on a richer auxiliary
series fit and derive its variance from the combined linear weights.
```

**Audit question:** Does the response distinguish projection inference from function inference and studentize the bias corrected center with its own variance?

<a id="w12-stop-5"></a>

## 5. Cross validation and out of sample risk

### 5.1 Prediction risk

Let $\mathcal D_n$ denote the training sample, and let $(Y^\circ,X^\circ)$ be an independent test observation. The conditional prediction risk of an estimator $\widehat\mu_c$ with tuning parameter $c$ is

$$
R(c\mid\mathcal D_n) =
\mathbb E[
\lbrace Y^\circ-\widehat\mu_c(X^\circ)\rbrace^2
\mid\mathcal D_n].
$$

Because $Y^\circ=\mu(X^\circ)+\varepsilon^\circ$,

$$
R(c\mid\mathcal D_n) =
\mathbb E[(\varepsilon^\circ)^2]
+
\mathbb E[
\lbrace\widehat\mu_c(X^\circ)-\mu(X^\circ)\rbrace^2
\mid\mathcal D_n].
$$

The irreducible noise does not depend on $c$. Minimizing prediction risk therefore targets integrated squared estimation error.

Training error reuses the same outcomes to fit and evaluate the model and is generally optimistic. A validation sample estimates out of sample risk for tuning. A test sample should be reserved for final evaluation rather than repeatedly reused for tuning.

### 5.2 Leave one out cross validation

For each $i$, fit the estimator without observation $i$ and call it $\widehat\mu_{-i,c}$. Leave one out cross validation uses

$$
\mathrm{CV}(c) =
\frac1n\sum_{i=1}^n
\lbrace Y_i-\widehat\mu_{-i,c}(X_i)\rbrace^2.
$$

The held out residual compares $Y_i$ with the leave one out prediction. It is not the difference between a leave one out fit and a full sample fit.

For series least squares at dimension $k$, let $P_k$ be the $n\times k$ design matrix and

$$
H_k=P_k(P_k'P_k)^{-1}P_k'
$$

the hat matrix. Put $\widehat Y=H_kY$ and let $H_{k,ii}$ denote leverage. When $H_{k,ii}\neq1$,

$$
Y_i-\widehat\mu_{-i,k}(X_i) =
\frac{Y_i-\widehat Y_i}{1-H_{k,ii}}.
$$

Therefore

$$
\mathrm{CV}(k) =
\frac1n\sum_{i=1}^n
\left(
\frac{Y_i-\widehat Y_i}{1-H_{k,ii}}
\right)^2.
$$

To see the algebra, delete row $p_k(X_i)'$ from the normal equations and apply the Sherman--Morrison inverse identity to $P_k'P_k-p_k(X_i)p_k(X_i)'$. Multiplying the resulting coefficient update by $p_k(X_i)'$ yields

$$
\widehat\mu_{-i,k}(X_i) =
\widehat Y_i-
\frac{H_{k,ii}}{1-H_{k,ii}}
(Y_i-\widehat Y_i),
$$

which rearranges to the displayed residual formula. The leverage correction is large when an observation is difficult to predict without itself. The identity avoids refitting the model $n$ times.

### 5.3 Regression and density objectives

The validation objective must match the statistical target. Regression uses held out squared prediction error. Density estimation cannot use an unobserved response, so integrated squared error is estimated differently. Since

$$
\int\lbrace\widehat f_h(x)-f(x)\rbrace^2\mkern3mu dx =
\int\widehat f_h(x)^2\mkern3mu dx
-2\mathbb E[\widehat f_h(X)]
+\int f(x)^2\mkern3mu dx,
$$

the last term does not depend on $h$. A leave one out estimate of the relevant part is

$$
\mathrm{CV}_{\mathrm{density}}(h) =
\int\widehat f_h(x)^2\mkern3mu dx
-\frac2n\sum_{i=1}^n\widehat f_{-i,h}(X_i).
$$

| Problem | Held out contribution | Risk targeted |
|---|---|---|
| Regression | $\lbrace Y_i-\widehat\mu_{-i,c}(X_i)\rbrace^2$ | Prediction risk, equivalently integrated regression error up to irreducible noise |
| Density | $\int\widehat f_h^2-(2/n)\sum_i\widehat f_{-i,h}(X_i)$ | Integrated squared density error up to $\int f^2$ |

The common principle is out of sample evaluation; the formula changes with the loss and target.

### 5.4 What cross validation does not guarantee

Cross validation targets the risk represented by its held out loss and candidate set. It does not automatically:

- select a dimension appropriate for pointwise inference;
- prove a minimax theorem;
- protect against distribution shift;
- make a misspecified loss scientifically relevant; or
- turn the selected estimator into an independent test of performance.

> [!IMPORTANT]
> **Board work 3 — Leave one out cross validation**
>
> Derive the leave one out residual identity from the Sherman--Morrison formula. Then compare:
>
> 1. training residual;
> 2. leave one out residual;
> 3. validation error; and
> 4. test error.

> [!TIP]
> **AI interaction 2 — Audit a cross validation objective**
>
> Copy the prompt below into an AI interface and audit the response.

```text
An analyst defines CV(c) as the average squared difference between the
leave one out fitted value and the full sample fitted value.

Explain why this is not regression prediction error. Write the correct
leave one out criterion, derive its linear smoother shortcut, and state what
risk it estimates and what it does not establish about inference.
```

**Audit question:** Does the response retain the held out outcome $Y_i$ in the residual?

<a id="w12-stop-6"></a>

## 6. Kernel and series methods compared

| Feature | Kernel/local polynomial | Series least squares |
|---|---|---|
| Approximation | Local around each evaluation point | Global in a chosen basis |
| Main tuning parameter | Bandwidth $h$ | Dimension $k$ |
| More flexible fit | Smaller $h$ | Larger $k$ |
| Variance pressure | Effective local sample size | Number and leverage of regressors |
| Boundary behavior | Local polynomial reproduction matters | Basis and knot/partition construction matter |
| Computation | Refit locally across $x$ | One global least squares fit |
| Conventional pointwise inference | Choose smaller $h$ so smoothing bias is negligible | Choose larger $k$ so approximation bias is negligible |
| RBC | Use a higher degree local fit to estimate leading smoothing bias | Use a richer auxiliary series fit to estimate leading approximation bias |
| CV use | Select bandwidth | Select dimension |

Neither method dominates universally. Smoothness, dimension, evaluation loss, boundary behavior, computational constraints, and inferential target determine the useful choice.

Week 14 uses kernel or series fits as possible first step nuisance estimators. The key question there will not be whether the entire function is estimated at a $\sqrt n$ rate, but whether its error enters the final low dimensional estimator at first order.

### Checkpoint 3

1. Why does smaller $h$ correspond to greater kernel flexibility while larger $k$ corresponds to greater series flexibility?
2. What is the analogue of bandwidth bias in a series estimator?
3. Why can a prediction optimal tuning parameter fail to be inference optimal?
4. Which features of an estimator matter when it is used only as a nuisance first step?
5. Why does fixed $k$ OLS inference normally target $\mu_k(x)$ rather than $\mu(x)$?
6. Which weights determine the standard error after a series bias correction?

## 7. Recap

- Series regression is least squares with a basis dimension that may grow.
- The population projection $\mu_k$ and underlying function $\mu$ are distinct targets.
- Error decomposes into coefficient estimation error and approximation error.
- Integrated squared error typically balances $k/n$ against squared approximation error.
- Supremum norm control requires basis growth factors and a correctly scaled stochastic term.
- Fixed $k$ pointwise inference is ordinary heteroskedasticity robust OLS inference for the projection $\mu_k(x)$.
- Pointwise inference for $\mu(x)$ combines the OLS CLT with approximation bias measured in standard error units.
- Feasible series inference uses the OLS sandwich; growing dimension validity also requires Gram matrix, leverage, moment, and residual variance controls.
- Series RBC uses a verified leading approximation expansion, a richer auxiliary fit, and the variance of the combined main and correction weights.
- Training, validation, test, and leave one out errors serve different roles.
- Linear smoothers admit a fast leave one out residual formula.
- Cross validation targets out of sample loss, not automatic inferential validity.
- Series and kernel methods express the same regularization principle through different tuning parameters.

## Notation introduced this week

- $p_k(x)$: $k$ dimensional series basis.
- $Q_k=\mathbb E[p_k(X)p_k(X)']$: population Gram matrix.
- $\mu_k=p_k'\beta_k$: population series projection.
- $r_k=\mu-\mu_k$: approximation error.
- $\zeta_k=\sup_x\Vert p_k(x)\Vert$: basis size factor.
- $V_k(x)$ and $\widehat s_k(x)$: pointwise projection variance constant and feasible OLS standard error.
- $B_k(x)$ and $\widehat B_k(x)$: leading series approximation bias and its auxiliary estimate.
- $a_{i,k}(x)-b_{i,k}(x)$: combined linear weight for the series RBC center.
- $H_k$: series least squares hat matrix.
- $H_{k,ii}$: leverage.
- $\mathrm{CV}(k)$: leave one out cross validation criterion.

## References

- Paul P. B. Eggermont and Vincent N. LaRiccia, *Maximum Penalized Likelihood Estimation, Volume II: Regression*, Springer, 2009.
- László Györfi et al., *A Distribution-Free Theory of Nonparametric Regression*, Springer, 2002.
- Martin J. Wainwright, *High-Dimensional Statistics: A Non-Asymptotic Viewpoint*, Cambridge University Press, 2019.
