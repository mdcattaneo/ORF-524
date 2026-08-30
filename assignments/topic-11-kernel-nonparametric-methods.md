# ORF 524 Practice Module 11: Kernel Based Nonparametric Methods

**Status:** Ready for Fall 2026  
**Last updated:** August 30, 2026  
**Primary chapter:** [Week 11: Kernel based nonparametric methods](../lectures/week-11.md)

## Purpose

This ungraded module develops local smoothing from exact expectation and variance calculations through bandwidth choice and inference. Problems 1--4 form the core bank. Problems 5--6 examine positivity and the quantifiers behind uniform inference.

<a id="problem-map"></a>

## Problem map

| Problem | Bank | Main task |
|---|---|---|
| [1. A fourth order kernel density estimator](#problem-1) | Core | Compute kernel moments, bias, variance, MSE, and optimal bandwidth rates. |
| [2. Bandwidth choice for pointwise inference](#problem-2) | Core | Compare centered limits, undersmoothing, shifted limits, and robust bias correction. |
| [3. Local polynomial regression and boundary adaptation](#problem-3) | Core | Derive local polynomial weights, reproduction, boundary behavior, bias, and variance. |
| [4. Density and local polynomial rates under dimension and derivatives](#problem-4) | Core | Compare how dimension, derivatives, and smoothness determine rates. |
| [5. Higher order kernels and positivity](#problem-5) | Further | Analyze signed higher order kernels and a positive part modification. |
| [6. Three different meanings of uniform](#problem-6) | Further | Distinguish pointwise, simultaneous, and uniform validity through their quantifiers. |

## Core practice

<a id="problem-1"></a>

### Problem 1. A fourth order kernel density estimator

**Chapter connections:** [Section 1.2: The exact expectation is a convolution](../lectures/week-11.md#12-the-exact-expectation-is-a-convolution), [Section 2.1: Pointwise variance](../lectures/week-11.md#21-pointwise-variance), and [Section 2.2: MSE and the smoothing tradeoff](../lectures/week-11.md#22-mse-and-the-smoothing-tradeoff)

Let $X_1,\ldots,X_n$ be i.i.d. with density $f$, and define

$$
K(u)=\frac{45}{32}\left(1-\frac{7u^2}{3}\right)(1-u^2)
\mathbf 1\lbrace|u|\leq1\rbrace.
$$

Let

$$
\widehat f_h(x)=\frac1{nh}\sum_{i=1}^nK\mkern-3mu\left(\frac{X_i-x}{h}\right).
$$

Assume $f$ has four bounded continuous derivatives in a neighborhood of an interior point $x$.

1. Verify that $\int K=1$, $\int u^jK(u)\mkern3mu du=0$ for $j=1,2,3$, and $\mu_4(K)=-1/21$.

2. Show that $R(K)=\int K(u)^2du=5/4$.

3. Derive the leading bias and variance at $x$.

4. Show that the pointwise MSE has order $h^8+(nh)^{-1}$, and obtain the MSE optimal bandwidth and MSE orders.

5. Compare those orders with a second order kernel. State the extra smoothness needed to benefit from the fourth order kernel.

<!-- Source lineage: Fall 2024 pset8, “Kernel density estimation with a fourth order kernel”; corrected to retain the negative optimal-MSE exponent. -->

[Back to the problem map](#problem-map)

<a id="problem-2"></a>

### Problem 2. Bandwidth choice for pointwise inference

**Chapter connections:** [Section 3.3: Three bandwidth regimes](../lectures/week-11.md#33-three-bandwidth-regimes), [Section 3.4: Feasible pointwise inference](../lectures/week-11.md#34-feasible-pointwise-inference), and [Section 3.6: Robust bias correction with an MSE optimal bandwidth](../lectures/week-11.md#36-robust-bias-correction-with-an-mse-optimal-bandwidth)

Suppose $f(x)>0$ at a fixed interior point, the order $p$ bias expansion holds, and

$$
\mathbb V\lbrace\widehat f_h(x)\rbrace
=\frac{R(K)f(x)}{nh}\lbrace1+o(1)\rbrace.
$$

Write $W_{i,h}(x)=h^{-1}K((X_i-x)/h)$. Assume the boundedness, moment, and bandwidth conditions needed for the calculations below.

1. With $s_h^2(x)=\mathbb V\lbrace\widehat f_h(x)\rbrace$ and $b_h(x)=\mathbb E[\widehat f_h(x)]-f(x)$, add and subtract $\mathbb E[\widehat f_h(x)]$ to decompose $(\widehat f_h(x)-f(x))/s_h(x)$ into a random component and a signal to noise ratio. Explain why consistency does not make the second component negligible automatically.

2. Express the random component as $\sum_{i=1}^n\xi_{i,n}(x)$. Verify that the summands have mean zero and total variance one. If $\mathbb E|W_{i,h}(x)-\mathbb EW_{i,h}(x)|^{2+\eta}=O(h^{-1-\eta})$ for some $\eta>0$, verify Lyapunov's condition when $nh\to\infty$ and obtain the centered CLT.

3. Construct a ratio consistent estimator of $s_h^2(x)$ from the empirical variance of $W_{i,h}(x)$, keeping track of the variance of the mean factor. State a sufficient fourth moment condition for ratio consistency.

4. Let $h=n^{-a}$. Find conditions on $a$ for $h\to0$, $nh\to\infty$, and $\sqrt{nh}h^p\to0$.

5. Show that the MSE optimal exponent fails the last condition. More precisely, if $h\sim c n^{-1/(2p+1)}$ and $B_f(x)\neq0$, derive the shifted limiting distribution of the statistic centered at $f(x)$.

6. Construct an undersmoothed pointwise interval. If the signal to noise ratio instead converges to a finite $\delta$, derive the limiting coverage of the same interval.

7. For a second order KDE, use a twice differentiable pilot kernel $L$ and define $\widehat f_b^{(2)}(x)=(nb^3)^{-1}\sum_iL^{(2)}((X_i-x)/b)$. Form the bias corrected estimator $\widehat f_h(x)-h^2\mu_2(K)\widehat f_b^{(2)}(x)/2$, express it as a sample average, and construct its robust standard error from the variance of the resulting summand. Explain why this can support inference with $h\asymp n^{-1/5}$ under additional smoothness and remainder conditions, and why retaining the original KDE standard error is invalid.

<!-- Source lineage: legacy Week 10 KDE inference slides; repairs the missing averaging factor and separates undersmoothing from robust bias correction. -->

[Back to the problem map](#problem-map)

<a id="problem-3"></a>

### Problem 3. Local polynomial regression and boundary adaptation

**Chapter connections:** [Section 4.1: Regression target and estimator](../lectures/week-11.md#41-regression-target-and-estimator), [Section 4.3: Polynomial reproduction and boundary behavior](../lectures/week-11.md#43-polynomial-reproduction-and-boundary-behavior), and [Section 4.6: Feasible inference from the OLS sandwich](../lectures/week-11.md#46-feasible-inference-from-the-ols-sandwich)

Let $Y_i=\mu(X_i)+\varepsilon_i$ with $\mathbb E[\varepsilon_i\mid X_i]=0$. For fixed $x$, define the local polynomial estimator of degree $q$ as in Week 11.

1. Derive its weighted least squares solution and express $\widehat\mu(x)$ as $\sum_i\ell_i(x)Y_i$.

2. Prove that if the local design matrix is nonsingular, local polynomial regression of degree $q$ reproduces every polynomial of degree at most $q$.

3. Deduce for local linear regression that

$$
\sum_i\ell_i(x)=1,
\qquad
\sum_i\ell_i(x)(X_i-x)=0.
$$

4. Consider a design supported on $[0,\infty)$ and the boundary point $x=0$. If $\mu(t)=a+bt$, show at the population smoothing level that a local constant estimator generally has bias of order $h$, while a local linear estimator reproduces $\mu(0)=a$.

5. Explain why negative local linear weights are compatible with boundary bias reduction.

6. Let $\mathbf X=(X_1,\ldots,X_n)$ and $\sigma^2(t)=\mathbb V(Y_i\mid X_i=t)$. Derive the exact finite sample conditional bias, variance, and MSE of $\widehat\mu(x)$. Then state the leading one dimensional interior local linear bias and variance under the regularity conditions in Week 11 and obtain the MSE optimal bandwidth order.

<!-- Source lineage: Fall 2024 pset8, “Linear smoothers, cross-validation”; focuses the core on representation and boundary behavior. -->

[Back to the problem map](#problem-map)

<a id="problem-4"></a>

### Problem 4. Density and local polynomial rates under dimension and derivatives

**Chapter connections:** [Section 5.1: A common pointwise bias and variance template](../lectures/week-11.md#51-a-common-pointwise-bias--variance-template) and [Section 5.2: Dimension and derivative order](../lectures/week-11.md#52-dimension-and-derivative-order)

Suppose a KDE for a density derivative of total order $s$ has bias order $h^p$, while a local polynomial estimator of degree $q$ for a regression derivative of total order $\nu$ has bias order $h^{r_{q,\nu}}$.

1. Using variance order $(nh^{d+2s})^{-1}$, derive the KDE's pointwise MSE, MSE optimal bandwidth, and optimal MSE orders.

2. Using conditional variance order $(nh^{d+2\nu})^{-1}$, derive the corresponding orders for the local polynomial estimator.

3. Specialize to level estimation with $s=\nu=0$ and $p=r_{q,\nu}=2$. Compare the KDE and local polynomial bandwidth and MSE rates in one dimension, and explain why equal rates do not imply equal constants.

4. Compare the effect of adding one covariate with the effect of increasing the derivative order by one. Explain the smoothness and moment cancellation considerations that determine $p$ and $r_{q,\nu}$.

5. Using Section 2.4 of Week 11, explain why these calculations provide estimator specific upper bound rate heuristics rather than a complete minimax theorem. State the additional lower bound claim that minimaxity requires, and explain why a pointwise calculation does not establish the $L_\infty$ minimax rate.

<!-- Source lineage: legacy multivariate and derivative-estimation slides; corrects the derivative-bias exponent and aligns density and regression smoothing. -->

[Back to the problem map](#problem-map)

## Further practice

<a id="problem-5"></a>

### Problem 5. Higher order kernels and positivity

**Chapter connections:** [Section 1.2: The exact expectation is a convolution](../lectures/week-11.md#12-the-exact-expectation-is-a-convolution) and [Section 2.2: MSE and the smoothing tradeoff](../lectures/week-11.md#22-mse-and-the-smoothing-tradeoff)

For the fourth order kernel in Problem 1:

1. Find its minimum on $[-1,1]$ and show that it takes negative values.

2. Define $\widetilde f_h(x)=\max\lbrace\widehat f_h(x),0\rbrace$. Prove pointwise that

$$
|\widetilde f_h(x)-f(x)|
\leq
|\widehat f_h(x)-f(x)|.
$$

3. Does $\widetilde f_h$ necessarily integrate to one? Explain why renormalizing it is a separate modification whose risk must be analyzed.

<!-- Source lineage: Fall 2024 pset8 positivity extension. -->

[Back to the problem map](#problem-map)

<a id="problem-6"></a>

### Problem 6. Three different meanings of uniform

**Chapter connections:** [Section 2.3: Pointwise, integrated, and supremum norm rates](../lectures/week-11.md#23-pointwise-integrated-and-supremum-norm-rates), [Section 2.4: What minimaxity means](../lectures/week-11.md#24-what-minimaxity-means), and [Section 3.7: Pointwise intervals are not simultaneous bands](../lectures/week-11.md#37-pointwise-intervals-are-not-simultaneous-bands)

For each statement below, write the correct order of quantifiers and explain whether it follows from a pointwise CLT:

1. coverage at one fixed $x$ for one fixed $f$;

2. simultaneous coverage for every $x\in\mathcal X$ for one fixed $f$;

3. pointwise coverage at $x$, uniformly over $f\in\mathcal F$; and

4. simultaneous coverage uniformly over both $x\in\mathcal X$ and $f\in\mathcal F$.

Explain why replacing the worst case infimum of coverage probabilities by a supremum makes a uniform validity claim nearly vacuous.

<!-- Source lineage: legacy uniform-inference slides; rewritten to repair the coverage quantifier. -->

[Back to the problem map](#problem-map)

## Completion check

After the core, you should be able to:

1. derive KDE bias and variance rather than quote their orders;

2. distinguish MSE and root MSE rates;

3. explain what a minimax rate asserts under pointwise, squared $L_2$, or $L_\infty$ loss and why an estimator specific upper bound is only half of the claim;

4. use the random plus signal to noise decomposition to explain why estimation optimal smoothing need not be inference valid without undersmoothing or a correctly studentized robust bias correction;

5. derive local polynomial linear smoother weights and their exact finite sample conditional bias and variance;

6. explain local linear boundary adaptation and the leading asymptotic bias--variance orders; and

7. compare how dimension, derivative order, and uniformity affect density and regression smoothing.
