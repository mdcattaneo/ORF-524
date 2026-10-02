# ORF 524 Practice Module 10: Solutions

**Status:** Released solutions  
**Last updated:** October 2, 2026  
**Practice module:** [Module 10](../topic-10-regression-two-step-estimation.md)  
**Primary chapter:** [Week 10](../../lectures/week-10.md)

Attempt each problem before consulting its solution. After studying a solution, close it and reconstruct the argument, including the assumptions used at each step.

<a id="solution-map"></a>

## Solution map

| Problem | Solution |
|---|---|
| [1. OLS without assuming a linear conditional mean](../topic-10-regression-two-step-estimation.md#problem-1) | [Solution](#problem-1) |
| [2. Binary response: likelihood, NLS, and weighting](../topic-10-regression-two-step-estimation.md#problem-2) | [Solution](#problem-2) |
| [3. The first step adjustment in a scalar two step estimator](../topic-10-regression-two-step-estimation.md#problem-3) | [Solution](#problem-3) |
| [4. Why feasible weighting can be first order free](../topic-10-regression-two-step-estimation.md#problem-4) | [Solution](#problem-4) |
| [5. Exact Normal regression geometry](../topic-10-regression-two-step-estimation.md#problem-5) | [Solution](#problem-5) |
| [6. A mean zero moment as a control variate](../topic-10-regression-two-step-estimation.md#problem-6) | [Solution](#problem-6) |

<a id="problem-1"></a>

## Problem 1. OLS without assuming a linear conditional mean

[Problem statement](../topic-10-regression-two-step-estimation.md#problem-1)

Differentiating the population criterion gives

$$
\mathbb E[X(Y-X'\beta_0)]=0,
$$

so $H\beta_0=\mathbb E[XY]$ and $\beta_0=H^{-1}\mathbb E[XY]$. For $d=\beta-\beta_0$,

$$
Y-X'\beta=\varepsilon-X'd.
$$

Consequently,

$$
\mathbb E[(Y-X'\beta)^2]
=\mathbb E[\varepsilon^2]
-2d'\mathbb E[X\varepsilon]
+d'Hd
=\mathbb E[\varepsilon^2]+d'Hd.
$$

Because $\widehat H_n\to_{\mathbb{P}}H$ and $H$ is positive definite, the smallest eigenvalue of $\widehat H_n$ is positive with probability approaching one. On that event, the sample normal equation gives

$$
0=\mathbb{E}_n[X(Y-X'\widehat\beta_n)]
=\mathbb{E}_n[X\varepsilon]
-\widehat H_n(\widehat\beta_n-\beta_0),
$$

and hence

$$
\widehat\beta_n-\beta_0
=\widehat H_n^{-1}\mathbb{E}_n[X\varepsilon].
$$

The additional condition $\mathbb E[\Vert X\varepsilon\Vert^2]<\infty$ makes the covariance matrix below finite and supplies the ordinary multivariate CLT (together with i.i.d. sampling). Thus,

$$
\sqrt n(\widehat\beta_n-\beta_0)
=H^{-1}\frac1{\sqrt n}\sum_{i=1}^nX_i\varepsilon_i
+o_{\mathbb{P}}(1)
\rightsquigarrow
\mathcal N(0,H^{-1}\Sigma H^{-1}),
$$

where $\Sigma=\mathbb E[XX'\varepsilon^2]$. The sign is positive because the derivative of $X(Y-X'\beta)$ is $-XX'$ and the generic representation has the additional leading minus sign.

With $\widehat\varepsilon_i=Y_i-X_i'\widehat\beta_n$, define

$$
\widehat\Sigma_n=\mathbb{E}_n[XX'\widehat\varepsilon^2],
\qquad
\widehat V_n=\widehat H_n^{-1}\widehat\Sigma_n\widehat H_n^{-1}.
$$

Under the corresponding plug in LLN conditions, $\widehat V_n\to_{\mathbb{P}}H^{-1}\Sigma H^{-1}$. Under $\mathbb E[\varepsilon^2\mid X]=\sigma^2$, $\Sigma=\sigma^2H$ and the asymptotic variance becomes $\sigma^2H^{-1}$.

In the example, $H=I_2$ because $\mathbb E[Z]=0$ and $\mathbb E[Z^2]=1$. Also,

$$
\mathbb E[XY] =
\begin{pmatrix}
\mathbb E[Z^2+U]\\
\mathbb E[Z(Z^2+U)]
\end{pmatrix} =
\begin{pmatrix}
1\\
0
\end{pmatrix},
$$

using $\mathbb E[U\mid Z]=0$ and $\mathbb E[Z^3]=0$. Thus $\beta_0=(1,0)'$, while $\mathbb E[Y\mid Z]=Z^2$. The linear projection is well defined and correctly estimated even though it does not equal the nonlinear conditional mean.

[Back to the solution map](#solution-map)

<a id="problem-2"></a>

## Problem 2. Binary response: likelihood, NLS, and weighting

[Problem statement](../topic-10-regression-two-step-estimation.md#problem-2)

The conditional log likelihood contribution is

$$
\ell(W,\beta)
=Y\log p_\beta+(1-Y)\log(1-p_\beta).
$$

Differentiation gives

$$
s(W,\beta)
=\frac{X\dot p_\beta}{p_\beta(1-p_\beta)}(Y-p_\beta).
$$

The squared loss criterion has gradient equal to minus twice

$$
\psi_{\mathrm{NLS}}(W,\beta)
=X\dot p_\beta(Y-p_\beta).
$$

Its derivative is

$$
\nabla_{\beta'}\psi_{\mathrm{NLS}}
=XX'\ddot p_\beta(Y-p_\beta)-XX'\dot p_\beta^2.
$$

At the correctly specified truth, the conditional expectation of the first term is zero. Therefore

$$
\Psi_{\beta,\mathrm{NLS}}(\beta_0)
=-\mathbb E[XX'\dot p_{\beta_0}^2].
$$

Because $\mathbb V(Y\mid X)=p_{\beta_0}(1-p_{\beta_0})$,

$$
B_{\mathrm{NLS}}
=\mathbb E[XX'\dot p_{\beta_0}^2p_{\beta_0}(1-p_{\beta_0})].
$$

Put $v_0(X)=p_{\beta_0}(1-p_{\beta_0})$. The infeasible WNLS criterion uses $v_0(X)^{-1}$ as a fixed weight, so its estimating function is

$$
\psi_{\mathrm{IWNLS}}(W,\beta)
=\frac{X\dot p_\beta}{v_0(X)}(Y-p_\beta).
$$

At $\beta_0$, its derivative and variance matrices are respectively $-\mathcal I(\beta_0)$ and $\mathcal I(\beta_0)$, where

$$
\mathcal I(\beta_0)
=\mathbb E\mkern-3mu\left[
\frac{XX'\dot p_{\beta_0}^2}{p_{\beta_0}(1-p_{\beta_0})}
\right].
$$

Thus infeasible WNLS has asymptotic variance $\mathcal I(\beta_0)^{-1}$, the same as the correctly specified MLE.

Given a consistent preliminary NLS estimator $\widetilde\beta_n$, define $\widetilde v_n(X)=p_{\widetilde\beta_n}(1-p_{\widetilde\beta_n})$ and minimize the squared residual criterion with $1/\widetilde v_n(X)$ held fixed. This produces the feasible second step equation

$$
0
=\mathbb{E}_n
\left[
\frac{X\dot p_\beta}{\widetilde v_n(X)}
(Y-p_\beta)
\right].
$$

The first step estimates the weights and the second estimates $\beta_0$. Under correct specification, smoothness, weights bounded away from zero, the required plug in LLNs, and the rate conditions used in the two step theorem, the population moment is orthogonal to the weight parameter. Feasible and infeasible WNLS therefore have the same first order representation.

Continuously updating can be implemented by repeatedly estimating $v^{(t)}(X)=p_{\widehat\beta^{(t)}}(1-p_{\widehat\beta^{(t)}})$ and minimizing weighted squared residuals with that iteration's weights held fixed. At a fixed point, the estimating equation weight gives

$$
\frac{X\dot p_\beta}{p_\beta(1-p_\beta)}(Y-p_\beta),
$$

which is exactly the likelihood score. By contrast, if $v_\beta=p_\beta(1-p_\beta)$, then

$$
\nabla_\beta
\mathbb{E}_n
\left[
\frac{(Y-p_\beta)^2}{v_\beta}
\right] =
\mathbb{E}_n
\left[
-\frac{2X\dot p_\beta(Y-p_\beta)}{v_\beta}
-\frac{(Y-p_\beta)^2\nabla_\beta v_\beta}{v_\beta^2}
\right].
$$

The second term comes from differentiating the weight, so this naive parameter weighted criterion does not generate the likelihood score and need not identify $\beta_0$.

Indeed, at $\beta_0$ the conditional expectation of the displayed gradient is

$$
-\frac{\nabla_\beta v_{\beta_0}(X)}{v_0(X)},
$$

which is not generally zero. Under correct specification, MLE, NLS, infeasible WNLS, feasible WNLS, and the continuously updated score equation identify $\beta_0$; the naive parameter weighted residual criterion is excluded from this list. Under misspecification, likelihood minimizes expected log loss while NLS minimizes expected squared loss. Their first order conditions and pseudo true minimizers generally differ, so a direct variance ranking would compare estimators of different targets.

[Back to the solution map](#solution-map)

<a id="problem-3"></a>

## Problem 3. The first step adjustment in a scalar two step estimator

[Problem statement](../topic-10-regression-two-step-estimation.md#problem-3)

If $\gamma_0$ is known, the coordinatewise mean value argument for the sample moment gives

$$
\sqrt n(\widetilde\beta_n-\beta_0)
=-\frac1a\frac1{\sqrt n}\sum_{i=1}^n m_i+o_{\mathbb{P}}(1),
$$

where $m_i=m(W_i,\beta_0,\gamma_0)$. Thus the influence function is $-m_i/a$.

For the feasible estimator, the first step satisfies

$$
\sqrt n(\widehat\gamma_n-\gamma_0)
=-\frac1c\frac1{\sqrt n}\sum_{i=1}^n g_i+o_{\mathbb{P}}(1).
$$

The corresponding second step mean value representation is

$$
0
=\frac1{\sqrt n}\sum_i m_i
+a\sqrt n(\widehat\beta_n-\beta_0)
+b\sqrt n(\widehat\gamma_n-\gamma_0)
+o_{\mathbb{P}}(1).
$$

Substitution gives

$$
\sqrt n(\widehat\beta_n-\beta_0)
=-\frac1a\frac1{\sqrt n}\sum_{i=1}^n
\left(m_i-\frac bcg_i\right)
+o_{\mathbb{P}}(1).
$$

The stacked Jacobian and its inverse are

$$
J=
\begin{pmatrix}
a&b\\
0&c
\end{pmatrix},
\qquad
J^{-1}=
\begin{pmatrix}
a^{-1}&-a^{-1}bc^{-1}\\
0&c^{-1}
\end{pmatrix}.
$$

The first coordinate of $-J^{-1}(m_i,g_i)'$ is the same influence function.

Its variance is

$$
V_\beta
=\frac1{a^2}
\mathbb V\mkern-3mu\left(m-\frac bcg\right)
=\frac1{a^2}
\left\lbrace
\mathbb V[m]
+\frac{b^2}{c^2}\mathbb V[g]
-2\frac bc\mathrm{Cov}(m,g)
\right\rbrace.
$$

Let $\widehat a,\widehat b,\widehat c$ be empirical plug in derivative estimates and put

$$
\widehat\varphi_i
=-\widehat a^{-1}
\left\lbrace
m(W_i,\widehat\beta_n,\widehat\gamma_n)
-\widehat b\widehat c^{-1}g(W_i,\widehat\gamma_n)
\right\rbrace.
$$

Then $\widehat V_\beta=\mathbb E_n[\widehat\varphi^2]$ is consistent under the corresponding plug in LLNs, and the standard error of $\widehat\beta_n$ is $\sqrt{\widehat V_\beta/n}$.

If $b=0$, the feasible and known nuisance influence functions agree. This is local insensitivity of the population moment to $\gamma$, not independence between the estimators or between $m$ and $g$.

[Back to the solution map](#solution-map)

<a id="problem-4"></a>

## Problem 4. Why feasible weighting can be first order free

[Problem statement](../topic-10-regression-two-step-estimation.md#problem-4)

Write $\varepsilon=Y-G(X'\beta_0)$. For any integrable vector $r(X)$ measurable with respect to $X$,

$$
\mathbb E[r(X)\varepsilon]
=\mathbb E[r(X)\mathbb E(\varepsilon\mid X)]
=0.
$$

This proves the population moment condition for every admissible weight.

For

$$
m(W,\beta,\gamma)
=\frac{q(X,\beta)}{v(X,\gamma)}\lbrace Y-G(X'\beta)\rbrace,
$$

differentiation with respect to $\gamma$ at the truth gives

$$
\mathbb E[\nabla_{\gamma'}m(W,\beta_0,\gamma_0)]
=\mathbb E\mkern-3mu\left[
-\frac{q(X,\beta_0)\nabla_{\gamma'}v(X,\gamma_0)}{v(X,\gamma_0)^2}
\mathbb E(\varepsilon\mid X)
\right]
=0.
$$

The moment is therefore orthogonal to the weight parameter. Put

$$
J=\mathbb E\mkern-3mu\left[
\frac{q(X,\beta_0)q(X,\beta_0)'}{v(X,\gamma_0)}
\right].
$$

Correct conditional mean and variance specification give $\Psi_\beta(\beta_0)=-J$ and

$$
B
=\mathbb E\mkern-3mu\left[
\frac{q(X,\beta_0)q(X,\beta_0)'\varepsilon^2}{v(X,\gamma_0)^2}
\right]
=J.
$$

Thus

$$
\varphi_\beta(W)
=J^{-1}\frac{q(X,\beta_0)}{v(X,\gamma_0)}\varepsilon,
\qquad
V_\beta=J^{-1}.
$$

There is no first step influence term because the nuisance derivative is zero.

For Bernoulli regression, take $G(u)=p(u)$, $q(X,\beta)=X\dot p(X'\beta)$, and $v(X,\gamma)=p(X'\gamma)\lbrace1-p(X'\gamma)\rbrace$. A preliminary consistent NLS estimator supplies the frozen feasible weights from Section 2.3. At $\beta=\gamma=\beta_0$, the weighted moment is the likelihood score. Because the weighting nuisance is orthogonal, using the preliminary estimate has no first order cost under the stated regularity conditions, so feasible WNLS, infeasible WNLS, and the correctly specified MLE share the same first order distribution.

[Back to the solution map](#solution-map)

<a id="problem-5"></a>

## Problem 5. Exact Normal regression geometry

[Problem statement](../topic-10-regression-two-step-estimation.md#problem-5)

Direct calculation gives

$$
P'=P,
\qquad
P^2=P,
\qquad
M'=M,
\qquad
M^2=M,
\qquad
PM=0.
$$

The range of $P$ is the $d$ dimensional column space of $X$, so $\mathrm{rank}(P)=d$ and $\mathrm{rank}(M)=n-d$.

Conditional on $X$,

$$
PY\sim\mathcal N(X\beta,\sigma^2P),
\qquad
MY\sim\mathcal N(0,\sigma^2M).
$$

Their conditional covariance is $\sigma^2PM=0$. They are jointly Gaussian, so they are independent. Also,

$$
\widehat\beta
=(X'X)^{-1}X'Y
\sim
\mathcal N\lbrace\beta,\sigma^2(X'X)^{-1}\rbrace.
$$

Because $n>d$, the residual degrees of freedom are positive. Cochran's theorem, or diagonalization of the projection $M$ with rank $n-d$, gives

$$
\frac{Y'MY}{\sigma^2}
=\frac{\varepsilon'M\varepsilon}{\sigma^2}
\sim\chi^2_{n-d},
$$

independently of $\widehat\beta$. With $S^2=Y'MY/(n-d)$, for any nonzero $r\in\mathbb R^d$,

$$
\frac{r'\widehat\beta-r'\beta}
{S\sqrt{r'(X'X)^{-1}r}}
\sim t_{n-d}.
$$

[Back to the solution map](#solution-map)

<a id="problem-6"></a>

## Problem 6. A mean zero moment as a control variate

[Problem statement](../topic-10-regression-two-step-estimation.md#problem-6)

The adjusted sample equation has derivative limit $a$, so

$$
\sqrt n(\widehat\beta_\lambda-\beta_0)
=-\frac1a\frac1{\sqrt n}\sum_{i=1}^n
\lbrace m(W_i,\beta_0)-\lambda g(W_i)\rbrace
+o_{\mathbb{P}}(1).
$$

Its asymptotic variance is

$$
V(\lambda)
=\frac1{a^2}
\left\lbrace
\mathbb V[m]
+\lambda^2\mathbb V[g]
-2\lambda\mathrm{Cov}(m,g)
\right\rbrace.
$$

When $\mathbb V[g]>0$, differentiation gives

$$
\lambda^\star
=\frac{\mathrm{Cov}(m,g)}{\mathbb V[g]}.
$$

The minimized variance is

$$
V(\lambda^\star)
=\frac1{a^2}
\left\lbrace
\mathbb V[m]
-\frac{\mathrm{Cov}(m,g)^2}{\mathbb V[g]}
\right\rbrace
\leq
\frac{\mathbb V[m]}{a^2},
$$

where nonnegativity follows from Cauchy--Schwarz. Equality holds exactly when $\mathrm{Cov}(m,g)=0$. The auxiliary moment acts as a control variate. Consequently, a procedure that estimates or uses an auxiliary nuisance equation can outperform a procedure that knows the nuisance value but discards the auxiliary mean zero information.

[Back to the solution map](#solution-map)
