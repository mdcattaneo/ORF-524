# ORF 524 Practice Module 14: Semiparametric Methods

**Status:** Ready for Fall 2026  
**Last updated:** August 30, 2026  
**Primary chapter:** [Week 14: Semiparametric methods](../lectures/week-14.md)

## Purpose

This ungraded module complements the restored model first route. Problems 2, 3, and 5 deepen the integrated squared density, partially linear, and average derivative calculations developed live. Problems 1, 4, and 6 deliberately extend the chapter to pathwise influence function calculations, missing outcomes, and cross fitting; these extensions may go beyond the material developed in class. In this module, the Core label identifies serious recommended practice and does not imply that every part was developed in the live route.

<a id="problem-map"></a>

## Problem map

| Problem | Bank | Main task |
|---|---|---|
| [1. Influence functions and transformations](#problem-1) | Core | Derive influence functions for means, variances, and smooth transformations. |
| [2. Integrated squared density and the Hoeffding decomposition](#problem-2) | Core | Work from plug in estimators to U-statistics, MSE expansions, and feasible inference. |
| [3. Orthogonal estimation in the partially linear model](#problem-3) | Core | Derive the orthogonal score, remainder conditions, and cross fitted estimator. |
| [4. Missing outcomes: IPW, first step adjustment, and augmentation](#problem-4) | Core | Compare inverse probability weighting with known, estimated, and augmented propensity scores. |
| [5. Weighted average derivatives](#problem-5) | Further | Derive equivalent representations and a U-statistic estimator for weighted average derivatives. |
| [6. Cross fitting is not a rate assumption](#problem-6) | Further | Separate the role of cross fitting from nuisance convergence rate requirements. |

## Core practice

<a id="problem-1"></a>

### Problem 1. Influence functions and transformations

**Chapter connections:** [Section 5.1: Setup](../lectures/week-14.md#51-setup) and [Section 5.4: General asymptotic linear representation](../lectures/week-14.md#54-general-asymptotic-linear-representation)

Let $W_1,\ldots,W_n$ be i.i.d. from $\mathbb{P}$.

1. For $\mu(\mathbb{P})=\mathbb E_{\mathbb{P}}[W]$, use a smooth submodel with score $s$ to show that

$$
\varphi_\mu(W)=W-\mu(\mathbb{P})
$$

   represents the pathwise derivative.

2. For $\sigma^2(\mathbb{P})=\mathbb E_{\mathbb{P}}[(W-\mu(\mathbb{P}))^2]$, derive an influence function without treating the mean as known.

3. If $\tau(\mathbb{P})=a\lbrace\mu(\mathbb{P})\rbrace$ for a continuously differentiable function $a$, derive the influence function of $\tau$.

4. Explain how the answers reproduce the delta method variance formulas from Week 6.

5. State the moment conditions needed for $\sqrt n$ rate Gaussian inference in each case.

<!-- New bridge from Week 6 to the influence function language of legacy Weeks 12--13. -->

[Back to the problem map](#problem-map)

<a id="problem-2"></a>

### Problem 2. Integrated squared density and the Hoeffding decomposition

**Chapter connections:** [Section 4: Integrated squared density](../lectures/week-14.md#4-integrated-squared-density) and [Section 6.2: Reading the integrated squared density proof through the three conditions](../lectures/week-14.md#62-reading-the-integrated-squared-density-proof-through-the-three-conditions)

Let $X_1,\ldots,X_n$ be i.i.d. with density $f$ on $\mathbb R^d$ and target $\theta=\int f^2$. Let $K$ be a bounded symmetric kernel and define

$$
K_h(u)=h^{-d}K(u/h),
\qquad
\widehat f_h(x)=\frac1n\sum_{i=1}^nK_h(x-X_i).
$$

Let

$$
L(t)=\int_{\mathbb R^d}K(u)K(u+t)du,
\qquad
L_h(t)=h^{-d}L(t/h),
$$

and assume that $L$ is a kernel of order $p$, meaning that its relevant mixed moments below total order $p$ vanish.

1. Starting from the plug in estimator, prove the exact V-statistic identity

$$
\widehat\theta_h^V
=\int_{\mathbb R^d}\widehat f_h(x)^2dx
=\frac1{n^2}\sum_{i=1}^n\sum_{j=1}^nL_h(X_i-X_j).
$$

2. Delete the diagonal and define

$$
\widehat\theta_h =
\binom{n}{2}^{-1}\sum_{i<j}L_h(X_i-X_j).
$$

Show that

$$
\theta_h=\mathbb E[\widehat\theta_h]
=\int f(x)(L_h*f)(x)dx.
$$

3. Under appropriate smoothness and integrability conditions, prove $\theta_h-\theta=O(h^p)$.

4. Define $g_h(x)=(L_h*f)(x)$, construct

$$
q_h(x_1,x_2) =
L_h(x_1-x_2)-g_h(x_1)-g_h(x_2)+\theta_h,
$$

verify that $q_h$ is degenerate, and derive the Hoeffding decomposition

$$
\widehat\theta_h-\theta_h
=\frac2n\sum_i\lbrace g_h(X_i)-\theta_h\rbrace+R_{n,h}.
$$

5. Prove, under suitable boundedness and smoothness conditions,

$$
\mathbb E\left[
\left\lbrace
2\lbrace g_h(X)-\theta_h\rbrace
-2\lbrace f(X)-\theta\rbrace
\right\rbrace^2
\right]
=O(h^{2p})
$$

and $\mathbb E[q_h(X_1,X_2)^2]=O(h^{-d})$. Use degeneracy to deduce $\mathbb E[R_{n,h}^2]=O(n^{-2}h^{-d})$.

6. Derive the representative MSE order

$$
h^{2p}+\frac1n+\frac1{n^2h^d}.
$$

State bandwidth conditions that make both the degenerate remainder and smoothing bias negligible at the $\sqrt n$ scale, and state when they are compatible for $h=n^{-a}$.

7. Derive the limiting influence function and asymptotic variance. Construct a leave one out influence function variance estimator and give the standard error on the correct scale.

8. Return to the plug in V-statistic and prove

$$
\widehat\theta_h^V =
\frac{n-1}{n}\widehat\theta_h
+
\frac{L(0)}{nh^d}.
$$

Why does the V-statistic require a stronger bandwidth condition, and what stronger smoothness order restriction makes its bias and diagonal conditions compatible?

<!-- Source lineage: legacy Fall 2025 Weeks 12--13 and Fall 2024 pset9, “Integrated squared density estimation”; restores the plug in derivation, first projection approximation, degenerate variance, MSE expansion, bandwidth audit, and feasible inference while correcting the diagonal and standard error conditions. -->

[Back to the problem map](#problem-map)

<a id="problem-3"></a>

### Problem 3. Orthogonal estimation in the partially linear model

**Chapter connections:** [Section 2: The partially linear model](../lectures/week-14.md#2-the-partially-linear-model), [Section 6.1: Three examples in common notation](../lectures/week-14.md#61-three-examples-in-common-notation), and [Section 6.3: Optional modern extension](../lectures/week-14.md#63-optional-modern-extension-orthogonality-and-cross-fitting)

Suppose

$$
Y=D\theta_0+g_0(X)+\varepsilon,
\qquad
\mathbb E[\varepsilon\mid D,X]=0,
$$

and define

$$
m_0(X)=\mathbb E[D\mid X],
\qquad
\ell_0(X)=\mathbb E[Y\mid X].
$$

1. Put $V=D-m_0(X)$. Derive both $\theta_0=\mathbb E[VY]/\mathbb E[VD]$ and $\theta_0=\mathbb E[V\lbrace Y-\ell_0(X)\rbrace]/\mathbb E[V^2]$, and prove that the denominators agree.

2. For

$$
\psi(W;\theta,m,\ell)
=\lbrace D-m(X)\rbrace[Y-\ell(X)-\theta\lbrace D-m(X)\rbrace],
$$

   verify the population moment condition at the truth.

3. Compute the directional derivatives of the population moment with respect to $m$ and $\ell$ and show that both vanish at the truth.

4. Derive the influence function and its heteroskedastic asymptotic variance.

5. Let $\delta m=\widehat m-m_0$ and $\delta\ell=\widehat\ell-\ell_0$. Derive the exact population nuisance remainder at $\theta_0$ and give sufficient $L^2$ rate conditions for it to be $o_{\mathbb{P}}(n^{-1/2})$.

6. Describe how local polynomial and series first steps implement residualization. For a common series matrix $P$, derive $\widehat\theta=D'M_PY/(D'M_PD)$. Then write a $K$ fold cross fitting algorithm, the residual on residual estimator, and an influence function standard error.

7. Explain separately the roles of orthogonality, cross fitting, nuisance rates, and the residual variation condition.

<!-- Source lineage: legacy Fall 2025 Week 12 partially linear model; augmented with modern orthogonal score and cross fitting analysis. -->

[Back to the problem map](#problem-map)

<a id="problem-4"></a>

### Problem 4. Missing outcomes: IPW, first step adjustment, and augmentation

**Chapter connections:** This extension applies [Section 5.3: Three key conditions for the first step](../lectures/week-14.md#53-three-key-conditions-for-the-first-step), [Section 5.4: General asymptotic linear representation](../lectures/week-14.md#54-general-asymptotic-linear-representation), and [Section 6.3: Optional modern extension](../lectures/week-14.md#63-optional-modern-extension-orthogonality-and-cross-fitting) to a model not developed in the live route.

**Extension beyond the live route:** The chapter does not develop missing outcome methods. Parts 1--5 restore the detailed inverse probability weighting exercise from the earlier assignments. Parts 6--8 then connect that calculation to the chapter's [optional orthogonality and cross fitting analysis](../lectures/week-14.md#63-optional-modern-extension-orthogonality-and-cross-fitting). Use the staged AI protocol in [`AGENTS.md`](AGENTS.md) for targeted hints if needed, and check every generated step from the stated missingness, moment, and positivity assumptions.

Let $Y^\star$ be an outcome, $T\in\lbrace0,1\rbrace$ indicate that it is observed, let $G=g(Y^\star)$ for a fixed measurable function $g$, and let the observed data be $W=(X,T,TG)$. The target is $\theta_0=\mathbb E[G]$.

First assume $G$ and $T$ are independent and put $p_0=\mathbb{P}(T=1)>0$.

1. Show that $\theta_0=\mathbb E[TG]/p_0$. Define

$$
\widehat p=\mathbb{E}_n[T],
\qquad
\widetilde\theta=\frac{\mathbb{E}_n[TG]}{\widehat p}.
$$

Show that $\widetilde\theta$ is the observed outcome sample average and give sufficient conditions for consistency.

2. Derive the asymptotic linear representation

$$
\sqrt n(\widetilde\theta-\theta_0)
=\frac1{\sqrt n}\sum_{i=1}^n
\frac{T_i}{p_0}(G_i-\theta_0)+o_{\mathbb{P}}(1).
$$

Identify its asymptotic variance, construct a consistent variance estimator, and give a feasible confidence interval.

3. Compare $\widetilde\theta$ with the oracle estimator $\mathbb{E}_n[TG]/p_0$ that treats $p_0$ as known. Derive both influence functions and explain why estimating $p_0$ changes the first order variance.

Now assume missing at random and positivity:

$$
\mathbb E[G\mid T,X]=m_0(X),
\qquad
p_0(X)=\mathbb{P}(T=1\mid X)\geq c>0.
$$

4. Show that

$$
\theta_0=\mathbb E\mkern-3mu\left[\frac{TG}{p_0(X)}\right].
$$

For

$$
\widehat\theta_{\mathrm{IPW}}
=\mathbb{E}_n\mkern-3mu\left[\frac{TG}{\widehat p(X)}\right],
$$

prove consistency when $\Vert\widehat p-p_0\Vert_\infty=o_{\mathbb{P}}(1)$, along with suitable positivity and moment conditions.

5. Suppose $p_0(x)=p(x;\gamma_0)$, where $p(x;\gamma)$ is known and smooth, and

$$
\sqrt n(\widehat\gamma-\gamma_0)
=\frac1{\sqrt n}\sum_{i=1}^n\varphi_\gamma(W_i)+o_{\mathbb{P}}(1).
$$

Derive the first step adjusted influence function of $\widehat\theta_{\mathrm{IPW}}$, its asymptotic variance, a consistent plug in variance estimator, and a feasible confidence interval.

6. Show that the regression and inverse probability identities

$$
\theta_0=\mathbb E[m_0(X)]
=\mathbb E\mkern-3mu\left[\frac{TG}{p_0(X)}\right]
$$

identify the same target. Consider the augmented moment

$$
\psi(W;\theta,m,p)
=m(X)-\theta+\frac{T}{p(X)}\lbrace G-m(X)\rbrace.
$$

Derive its expectation at $\theta_0$ for generic $m$ and $p$.

7. Prove that the augmented moment is correct if either $m=m_0$ or $p=p_0$. Explain why this is a population bias statement rather than finite sample immunity to poor estimation. Verify Neyman orthogonality at $(m_0,p_0)$ and identify the product remainder.

8. Give the augmented influence function, a cross fitted estimator, and its standard error. State product rate and moment conditions sufficient for asymptotic linearity. Explain separately what fails without positivity for identification, variance, and numerical stability.

<!-- Source lineage: Fall 2024 pset8 and Fall 2025 Problem Set 4, “Missing data, IPW”; restores the full constant probability, conditional probability, parametric first step, variance, and feasible inference progression, then retains the augmented orthogonal extension. -->

[Back to the problem map](#problem-map)

## Further practice

<a id="problem-5"></a>

### Problem 5. Weighted average derivatives

**Chapter connection:** [Section 3: Average derivatives](../lectures/week-14.md#3-average-derivatives)

Let $X$ be scalar with density $f$, let $\mu(x)=\mathbb E[Y\mid X=x]$, and, for a differentiable weight $w$, define

$$
\theta_w=\mathbb E[w(X)\mu'(X)].
$$

1. Write the parameter as an ordinary integral.

2. State boundary and smoothness conditions and derive the general integration by parts identity

$$
\theta_w
=-\int \mu(x)\lbrace w'(x)f(x)+w(x)f'(x)\rbrace\mkern3mu dx.
$$

3. Now take $w=f$. Show that the identity simplifies to

$$
\theta_f=-2\mathbb E[Yf'(X)].
$$

4. Display the boundary term before setting it to zero.

5. Construct a leave one out kernel estimator of the density weighted right hand side. If $K$ is symmetric, show that it equals the U-statistic with pairwise kernel

$$
u_h(W_i,W_j) =
-\frac1{h^2}(Y_i-Y_j)
K'\mkern-3mu\left(\frac{X_i-X_j}{h}\right).
$$

6. Explain how its $\sqrt n$ rate proof would separate smoothing bias, first projection, and degenerate remainder.

7. What changes for a vector regressor and a directional derivative?

<!-- Source lineage: legacy Fall 2025 Week 12 average derivative example. -->

[Back to the problem map](#problem-map)

<a id="problem-6"></a>

### Problem 6. Cross fitting is not a rate assumption

**Chapter connections:** [Section 5.3: Three key conditions for the first step](../lectures/week-14.md#53-three-key-conditions-for-the-first-step) and [Section 6.3: Optional modern extension](../lectures/week-14.md#63-optional-modern-extension-orthogonality-and-cross-fitting)

Suppose a cross fitted orthogonal estimator has remainder bounded by

$$
O_{\mathbb{P}}(r_{m,n}r_{\ell,n}+r_{m,n}^2),
$$

where $r_{m,n}=\Vert\widehat m-m_0\Vert_{L^2}$ and $r_{\ell,n}=\Vert\widehat\ell-\ell_0\Vert_{L^2}$.

1. If $r_{m,n}=n^{-a}$ and $r_{\ell,n}=n^{-b}$, characterize the pairs $(a,b)$ that make the displayed remainder $o(n^{-1/2})$.

2. Does $a=b=1/4$ suffice without a little $o$ improvement? Explain.

3. Give an asymmetric pair of rates that works and one that fails.

4. Explain why cross fitting by itself implies none of these inequalities.

5. Give two additional conditions, unrelated to nuisance rates, that are needed for valid partially linear inference.

<!-- New diagnostic problem addressing a common misinterpretation of cross fitting. -->

[Back to the problem map](#problem-map)

## Completion check

After the core, you should be able to:

1. derive simple influence functions by pathwise differentiation;

2. decompose a kernel U-statistic and state separate bias and remainder conditions;

3. prove orthogonality of the partially linear score;

4. recognize and bound a nuisance product remainder;

5. implement the logic of cross fitting without overstating what it guarantees; and

6. derive both the first step adjustment for IPW and the product remainder for augmented IPW.
