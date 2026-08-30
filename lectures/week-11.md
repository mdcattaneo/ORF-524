# ORF 524: Statistical Theory and Methods

## Week 11: Kernel based nonparametric methods

**Status:** Ready for Fall 2026  
**Last updated:** August 30, 2026  
**Meetings:** Monday, November 9, and Wednesday, November 11, 2026

<!-- Legacy sources: 11 ORF 524/lectures/ORF524-Fall2025-Lecture-Week-10.tex and
ORF524-Fall2025-Lecture-Week-10_Lecture.pdf; 09 ORF 524/lectures/
ORF524-Fall2024-Lecture-Week-10_Notes.pdf and Week-11_Notes.pdf; 09 ORF 524/assignments/pset8.tex;
and earlier kernel and local polynomial problem sets. -->

## Central question

How can local averaging estimate an unknown density or regression function, and how do bandwidth, smoothness, dimension, and derivative order determine bias, variance, and valid inference?

## Learning goals

By the end of the week, you should be able to:

1. define kernel density and local polynomial regression estimators;
2. derive leading interior bias and variance terms for kernel density estimation;
3. balance smoothing bias and sampling variance to obtain rates and bandwidth orders, and interpret minimax statements under pointwise, $L_2$, and $L_\infty$ loss;
4. decompose a standardized estimation error into a random component and a bias to standard error signal to noise ratio, and compare undersmoothing with robust bias correction based on an MSE optimal bandwidth;
5. derive the exact conditional bias and variance of a local polynomial estimator, connect them to their leading asymptotic orders, and explain boundary adaptation;
6. conduct feasible local polynomial inference by combining a conditional CLT with OLS sandwich variance estimation, and explain how RBC changes both the center and the standard error; and
7. compare how dimension and derivative order affect KDE and local polynomial rates, while distinguishing pointwise inference, simultaneous bands, and uniform validity over a function class.

## In-class route

This is one continuous route through the week, not a division into class meetings. Proceed in order, stop wherever the discussion ends, and resume at the next unfinished stop. The linked sections contain the full student facing exposition; the route identifies the concepts and activities to foreground live.

| Stop | Live focus | Mode |
|---|---|---|
| **W11.1** | [From the empirical CDF to kernel density estimation](#w11-stop-1) | Discuss + exact expectation + Checkpoint 1 |
| **W11.2** | [Bias, variance, rates, bandwidth, and minimaxity](#w11-stop-2) | Board work + minimax interpretation + AI bandwidth audit |
| **W11.3** | [Pointwise inference, smoothing bias, and robust bias correction](#w11-stop-3) | Board work + inference audit + robust bias correction + Checkpoint 2 |
| **W11.4** | [Local polynomial regression, boundary behavior, and inference](#w11-stop-4) | Board work 3 + Board work 4 + AI residual audit + feasible inference + RBC |
| **W11.5** | [KDE and local polynomial rates across dimension and derivatives](#w11-stop-5) | Discuss + comparison table + Checkpoint 3 |

The exact convolution, leading bias and variance, and bandwidth balance form one continuous board arc across Stops W11.1--W11.2. Begin Stop W11.3 with the generic random plus signal to noise decomposition before specializing it to the KDE and its bandwidth regimes. In Stop W11.4, keep the weighted least squares form, polynomial reproduction, boundary comparison, conditional bias and variance, OLS inference, and RBC together as one regression smoothing arc.

## How to use this chapter

**Prepare:** Review Week 6 studentization and Week 10 weighted least squares. Before class, derive the expectation of a kernel density estimator by a change of variables and attempt Checkpoint 1.

**In class:** Follow the continuous route above through density estimation, bandwidth choice, rates, pointwise inference, local polynomial regression, boundary behavior, feasible regression inference, RBC, multivariate and derivative extensions, and forms of uniform inference. Prepared responses supply an alternative for every live prompt if AI is unavailable.

**Review:** For every smoothing claim, identify the evaluation point, smoothness class, kernel moments, bandwidth sequence, norm, and whether the claim is an expansion, a stochastic rate, or a limit theorem.

**Practice:** The matching assignment complements this chapter with ungraded problems for consolidating the material, extending selected ideas, and preparing for cumulative review.

**Prerequisites:** Weeks 1--6 and 9--10; Taylor expansion; conditional expectation; LLN, triangular array CLT, stochastic orders, studentization, and weighted least squares.

## Full chapter map

1. [Local smoothing and kernel density estimation](#1-local-smoothing-and-kernel-density-estimation)
2. [Bias, variance, rates, and bandwidth](#2-bias-variance-rates-and-bandwidth)
3. [Pointwise inference and its bias problem](#3-pointwise-inference-and-its-bias-problem)
4. [Local polynomial regression](#4-local-polynomial-regression)
5. [Comparing KDE and local polynomial smoothing](#5-comparing-kde-and-local-polynomial-smoothing)
6. [Recap](#6-recap)

<a id="w11-stop-1"></a>

## 1. Local smoothing and kernel density estimation

### 1.1 The nonparametric target

Let $X_1,\ldots,X_n$ be i.i.d. real valued observations with Lebesgue density $f$. Unlike a finite dimensional parametric model, the target is the function $x\mapsto f(x)$. Estimating the value $f(x)$ at one fixed point, estimating the entire curve in an integrated norm, and constructing a simultaneous band are different statistical problems.

A bandwidth $h=h_n>0$ defines a local neighborhood around the evaluation point. The kernel density estimator is

$$
\widehat f_h(x)
=\frac1{nh}\sum_{i=1}^n
K\mkern-3mu\left(\frac{X_i-x}{h}\right).
$$

One route to this estimator starts from the empirical CDF

$$
\widehat F_n(x)=\frac1n\sum_{i=1}^n\mathbf 1\lbrace X_i\leq x\rbrace.
$$

It is a step function, so differentiating it directly does not produce an ordinary density. A symmetric difference quotient instead gives

$$
\frac{\widehat F_n(x+h)-\widehat F_n(x-h)}{2h}
=\frac1{2nh}\sum_{i=1}^n\mathbf 1\lbrace x-h<X_i\leq x+h\rbrace,
$$

which is the KDE with the uniform kernel $K(u)=\tfrac12\mathbf 1\lbrace|u|\leq1\rbrace$ up to an immaterial endpoint convention. General kernels replace the hard neighborhood boundary by a weighted local average. This derivation makes the bandwidth's role concrete: it determines both the neighborhood and the effective local sample size.

The kernel $K:\mathbb R\to\mathbb R$ is integrable and satisfies $\int K(u)\mkern3mu du=1$. A kernel is of order $p$ if

$$
\int u^jK(u)\mkern3mu du=0,\qquad j=1,\ldots,p-1,
$$

and

$$
\mu_p(K)=\int u^pK(u)\mkern3mu du\neq0.
$$

Write

$$
R(K)=\int K(u)^2\mkern3mu du.
$$

A nonnegative symmetric kernel is ordinarily of order two. Higher order kernels reduce leading bias for sufficiently smooth densities but generally take negative values, so the resulting estimate need not itself be nonnegative.

### 1.2 The exact expectation is a convolution

At a point where the relevant integrals exist,

$$
\mathbb E[\widehat f_h(x)]
=\int K(u)f(x+hu)\mkern3mu du.
$$

This is an exact identity. A Taylor expansion converts it into a leading approximation only after smoothness and an interior point condition have been imposed.

If $f$ has $p$ continuous derivatives near an interior point $x$, the Taylor remainder is controlled, and $K$ has order $p$, then

$$
\mathbb E[\widehat f_h(x)]-f(x)
=h^p\frac{\mu_p(K)}{p!}f^{(p)}(x)+o(h^p).
$$

Define

$$
B_f(x)=\frac{\mu_p(K)}{p!}f^{(p)}(x).
$$

Then the leading bias is $h^pB_f(x)$. This expansion is generally false at a support boundary for an uncorrected symmetric KDE because the local kernel mass is truncated.

### Checkpoint 1

1. Which step is exact: the convolution identity or the Taylor bias formula?
2. Which kernel moments remove the lower order Taylor terms?
3. Why can a fourth order kernel be negative?
4. What goes wrong if $x$ is at the boundary of a compact support?

<a id="w11-stop-2"></a>

## 2. Bias, variance, rates, and bandwidth

### 2.1 Pointwise variance

Put

$$
W_{i,h}(x)=\frac1hK\mkern-3mu\left(\frac{X_i-x}{h}\right).
$$

Because $\widehat f_h(x)=\mathbb{E}_nW_{i,h}(x)$,

$$
\mathbb V[\widehat f_h(x)]
=\frac1n\mathbb V[W_{i,h}(x)].
$$

If $h\to0$, $f$ is continuous at $x$, and $R(K)<\infty$, then

$$
\mathbb V[\widehat f_h(x)]
=\frac{R(K)f(x)}{nh}
+o\mkern-3mu\left(\frac1{nh}\right).
$$

The term $1/(nh)$ is the effective local sample size cost: only observations within a neighborhood of width proportional to $h$ contribute materially.

A direct sample variance estimator is

$$
\widehat{\mathbb V}\lbrace\widehat f_h(x)\rbrace =
\frac1{n(n-1)}
\sum_{i=1}^n
\lbrace W_{i,h}(x)-\overline W_h(x)\rbrace^2,
$$

where $\overline W_h(x)=\widehat f_h(x)$. The factor $1/n$ outside the sample variance is essential because the target is the variance of a sample average.

### 2.2 MSE and the smoothing tradeoff

Combining the leading variance and squared bias gives

$$
\mathrm{MSE}\lbrace\widehat f_h(x)\rbrace =
\frac{R(K)f(x)}{nh}
+h^{2p}B_f(x)^2
+o\mkern-3mu\left(\frac1{nh}+h^{2p}\right).
$$

Ignoring constants, balance

$$
\frac1{nh}\asymp h^{2p}.
$$

This gives

$$
h_{\mathrm{MSE}}\asymp n^{-1/(2p+1)}
$$

and

$$
\mathrm{MSE}\lbrace\widehat f_{h_{\mathrm{MSE}}}(x)\rbrace
\asymp n^{-2p/(2p+1)}.
$$

The root MSE rate is $n^{-p/(2p+1)}$, slower than $n^{-1/2}$ for every fixed finite $p$. Increasing smoothness permits faster nonparametric rates, but never changes the fact that the target is an unknown function rather than a fixed dimensional parameter.

For a nonnegative weight $w$ such that the displayed integrals and remainders are finite, integrated MSE has the expansion

$$
\int \mathbb E[(\widehat f_h(x)-f(x))^2]w(x)\mkern3mu dx =
\frac{R(K)}{nh}\int f(x)w(x)\mkern3mu dx
+h^{2p}\int B_f(x)^2w(x)\mkern3mu dx
+o\mkern-3mu\left(\frac1{nh}+h^{2p}\right)
$$

under conditions that justify integrating the pointwise expansions. It therefore has the same leading rate as pointwise MSE, but different constants and a different loss. With $w\equiv1$, the first leading integral equals one for a density.

### 2.3 Pointwise, integrated, and supremum norm rates

Under appropriate smoothness, tail, and empirical process conditions,

$$
|\widehat f_h(x)-f(x)| =
O_{\mathbb{P}}\mkern-3mu\left(
h^p+\frac1{\sqrt{nh}}
\right)
$$

at a fixed interior point, while a typical compact set supremum norm bound is

$$
\sup_{x\in\mathcal X}
|\widehat f_h(x)-f(x)| =
O_{\mathbb{P}}\mkern-3mu\left(
h^p+\sqrt{\frac{\log n}{nh}}
\right).
$$

The logarithm is the price of controlling many evaluation points simultaneously. These displays are rates, not automatically asymptotic distributions or confidence band theorems.

The three losses should not be conflated:

- pointwise loss asks about one fixed $x$;
- integrated loss averages squared error across $x$; and
- supremum norm loss is governed by the worst evaluation point and typically pays a logarithmic factor.

### 2.4 What minimaxity means

A convergence rate for one estimator is not yet a statement of optimality. Minimax analysis first fixes a function class and a loss, then asks for the best worst case risk among all estimators. For example, let $\mathcal H(\beta,L)$ be a regular class of densities with Hölder smoothness $\beta>0$ and radius $L$, and let $x_0$ be a fixed interior point. The pointwise minimax squared risk is

$$
R^\star_{n,x_0} =
\inf_{\widetilde f_n}
\sup_{f\in\mathcal H(\beta,L)}
\mathbb E_f
\left[
\lbrace\widetilde f_n(x_0)-f(x_0)\rbrace^2
\right],
$$

where the infimum ranges over every estimator based on the sample, not only KDEs. Analogously, define

$$
R^\star_{n,2} =
\inf_{\widetilde f_n}
\sup_{f\in\mathcal H(\beta,L)}
\mathbb E_f
\Vert\widetilde f_n-f\Vert_2^2
$$

and

$$
R^\star_{n,\infty} =
\inf_{\widetilde f_n}
\sup_{f\in\mathcal H(\beta,L)}
\mathbb E_f
\Vert\widetilde f_n-f\Vert_\infty.
$$

Under standard regularity and nondegeneracy conditions in one dimension, the benchmark rates are

| Loss convention | Minimax risk order | Corresponding estimation error scale |
|---|---:|---:|
| Pointwise squared error at fixed $x_0$ | $n^{-2\beta/(2\beta+1)}$ | $n^{-\beta/(2\beta+1)}$ |
| Squared $L_2$ loss | $n^{-2\beta/(2\beta+1)}$ | $n^{-\beta/(2\beta+1)}$ |
| Expected $L_\infty$ loss | $(\log n/n)^{\beta/(2\beta+1)}$ | same as the risk convention shown |

If squared $L_\infty$ loss is used, square the final rate. Likewise, if “$L_2$ rate” refers to $\mathbb E_f\Vert\widetilde f_n-f\Vert_2$ rather than expected squared $L_2$ loss, take the square root of the displayed squared risk rate. Stating the loss convention prevents an apparent factor of two disagreement in exponents.

The notation $R_n^\star\asymp r_n$ contains two logically distinct claims:

1. **Upper bound:** there exists an estimator whose worst case risk over $\mathcal H(\beta,L)$ is at most a constant times $r_n$.
2. **Lower bound:** every estimator has worst case risk at least a constant times $r_n$.

The bias--variance calculations above can establish the upper bound order for a suitably chosen KDE, with kernel order and boundary treatment matched to the smoothness class. Balancing $h^{2\beta}$ with $(nh)^{-1}$ gives the pointwise and squared $L_2$ rate. Balancing $h^\beta$ with $\sqrt{\log n/(nh)}$ gives the $L_\infty$ rate and explains its logarithmic penalty. A separate information theoretic lower bound argument is needed before calling either rate minimax.

Minimaxity is a worst case estimation statement, not a confidence coverage statement. It also assumes that the smoothness class is specified; selecting a rate without knowing $\beta$ is the distinct problem of adaptation.

> [!IMPORTANT]
> **Board work 1 — KDE bias, variance, and optimal bandwidth**
>
> Starting from the exact convolution and variance identities:
>
> 1. derive the leading bias and variance;
> 2. balance them to obtain the MSE optimal bandwidth order;
> 3. calculate the resulting MSE and root MSE rates; and
> 4. identify every assumption that would fail at a boundary point; and
> 5. distinguish the KDE upper bound calculation from the lower bound required for a minimax claim under pointwise, squared $L_2$, and $L_\infty$ loss.

> [!TIP]
> **AI interaction 1 — Audit a bandwidth calculation**
>
> Copy the prompt below into an AI interface and audit the response.

```text
A proposed solution claims that an order p KDE has bias h^(p-1),
variance 1/(n h^2), optimal bandwidth n^(-1/(2p)), and root n rate.

Find the first invalid step, derive the correct pointwise orders, and state
the interior smoothness and kernel moment assumptions. Then explain why a
higher order kernel may cease to be nonnegative.
```

**Audit question:** Does the response distinguish MSE from root MSE and retain the factor $h$ created by the change of variables?

<a id="w11-stop-3"></a>

## 3. Pointwise inference and its bias problem

### 3.1 Random component plus signal to noise ratio

Start with any scalar estimator $\widehat\theta_n$ of a target $\theta$. Define its exact standard deviation and bias by

$$
s_n^2=\mathbb V(\widehat\theta_n)>0,
\qquad
b_n=\mathbb E[\widehat\theta_n]-\theta.
$$

The denominator below is the standard deviation $s_n=\sqrt{\mathbb V(\widehat\theta_n)}$, not the variance. Adding and subtracting $\mathbb E[\widehat\theta_n]$ gives the exact decomposition

$$
\frac{\widehat\theta_n-\theta}{s_n} =
\underbrace{
\frac{\widehat\theta_n-\mathbb E[\widehat\theta_n]}{s_n}
}_{\text{random component}}
+
\underbrace{
\frac{\mathbb E[\widehat\theta_n]-\theta}{s_n}
}_{\text{signal to noise ratio}}
=Z_n+\delta_n.
$$

This identity is algebra, not an asymptotic approximation. The random component $Z_n$ has mean zero and variance one, although it need not be Gaussian in finite samples. The deterministic term

$$
\delta_n=\frac{b_n}{s_n}
$$

measures the estimator's centering error in standard error units. In this inference calculation, the “signal” is smoothing bias and the “noise” is sampling standard deviation.

Suppose a central limit theorem gives $Z_n\rightsquigarrow\mathcal N(0,1)$. If $\delta_n\to\delta\in\mathbb R$, then Slutsky's theorem gives

$$
\frac{\widehat\theta_n-\theta}{s_n}
\rightsquigarrow
\mathcal N(\delta,1).
$$

Thus a standard Normal approximation centered at the target requires the stronger condition $\delta_n\to0$. Consistency alone does not ensure this: consistency can hold when both $b_n\to0$ and $s_n\to0$ even though $b_n/s_n$ converges to a nonzero constant or diverges.

For pointwise KDE inference, take

$$
\widehat\theta_n=\widehat f_h(x),
\qquad
\theta=f(x).
$$

At a fixed interior point with $f(x)>0$, the expansions from Sections 1--2 give

$$
b_h(x)
=\mathbb E[\widehat f_h(x)]-f(x)
=h^pB_f(x)+o(h^p)
$$

and

$$
s_h^2(x)
=\mathbb V\lbrace\widehat f_h(x)\rbrace
=\frac{R(K)f(x)}{nh}\lbrace1+o(1)\rbrace.
$$

Consequently,

$$
\delta_n(x)
=\frac{b_h(x)}{s_h(x)}
=\frac{B_f(x)}{\sqrt{R(K)f(x)}}
\sqrt{n h^{2p+1}}
+o\mkern-3mu\left(\sqrt{n h^{2p+1}}\right).
$$

The bandwidth therefore controls not only estimation error but also the centering of the inferential statistic.

### 3.2 Why the random component is asymptotically Gaussian

Let

$$
W_{i,h}(x)
=\frac1hK\mkern-3mu\left(\frac{X_i-x}{h}\right),
\qquad
\widehat f_h(x)=\frac1n\sum_{i=1}^nW_{i,h}(x).
$$

The random component can be written as a triangular array sum

$$
Z_n(x)
=\frac{\widehat f_h(x)-\mathbb E[\widehat f_h(x)]}
{\sqrt{\mathbb V\lbrace\widehat f_h(x)\rbrace}}
=\sum_{i=1}^n\xi_{i,n}(x),
$$

where

$$
\xi_{i,n}(x) =
\frac{W_{i,h}(x)-\mathbb E[W_{i,h}(x)]}
{\sqrt{n\mathbb V\lbrace W_{i,h}(x)\rbrace}}.
$$

For every $n$, the summands are independent, have mean zero, and their variances sum to one. A convenient sufficient route to the Lindeberg--Feller CLT is Lyapunov's condition. If, for some $\eta>0$, the kernel and density conditions imply

$$
\mathbb E\mkern-3mu\left[
|W_{i,h}(x)-\mathbb E W_{i,h}(x)|^{2+\eta}
\right]
=O(h^{-1-\eta})
$$

and $\mathbb V\lbrace W_{i,h}(x)\rbrace\asymp h^{-1}$, then

$$
\sum_{i=1}^n
\mathbb E[|\xi_{i,n}(x)|^{2+\eta}]
=O\mkern-3mu\left((nh)^{-\eta/2}\right)
\to0
$$

whenever $nh\to\infty$. Bounded compactly supported $K$, local boundedness and continuity of $f$, $f(x)>0$, and the required kernel moments provide one simple sufficient set of conditions. Hence

$$
Z_n(x)\rightsquigarrow\mathcal N(0,1).
$$

Using the leading variance constant, the same result can be written as

$$
\frac{\sqrt{nh}\lbrace
\widehat f_h(x)-\mathbb E[\widehat f_h(x)]
\rbrace}
{\sqrt{R(K)f(x)}}
\rightsquigarrow
\mathcal N(0,1).
$$

Writing $r_h(x)=b_h(x)-h^pB_f(x)=o(h^p)$, if the bias remainder also satisfies $\sqrt{nh}\mkern3mu r_h(x)\to0$, adding and subtracting the leading bias gives the familiar bias centered form

$$
\frac{\sqrt{nh}\lbrace
\widehat f_h(x)-f(x)-h^pB_f(x)
\rbrace}
{\sqrt{R(K)f(x)}}
\rightsquigarrow
\mathcal N(0,1).
$$

This proof reveals why $nh$ is the effective sample size. If $nh$ does not diverge, the local average need not accumulate enough small contributions for this Gaussian approximation.

### 3.3 Three bandwidth regimes

Combining the centered CLT with the signal to noise calculation produces three regimes when $B_f(x)\neq0$:

1. **Undersmoothing:** If $nh\to\infty$ and $nh^{2p+1}\to0$, then $\delta_n(x)\to0$ and

$$
\frac{\widehat f_h(x)-f(x)}{s_h(x)}
\rightsquigarrow\mathcal N(0,1).
$$

2. **MSE optimal smoothing:** If $h\sim c n^{-1/(2p+1)}$ for $c>0$, then

$$
\delta_n(x)
\to
\delta(x)
=\frac{c^{p+1/2}B_f(x)}{\sqrt{R(K)f(x)}},
$$

   so the statistic centered at $f(x)$ converges to $\mathcal N\lbrace\delta(x),1\rbrace$ rather than $\mathcal N(0,1)$.

3. **Oversmoothing for conventional inference:** If $nh^{2p+1}\to\infty$, then the smoothing bias dominates the standard error in magnitude and conventional bias ignoring inference fails.

If $B_f(x)=0$, the displayed leading bias vanishes and the next nonzero bias term determines the regimes. One must not impose an undersmoothing condition based on a leading term that is identically zero.

> [!IMPORTANT]
> **Board work 2 — Random fluctuation and smoothing signal**
>
> 1. add and subtract $\mathbb E[\widehat\theta_n]$ to obtain the exact random plus signal to noise decomposition;
> 2. specialize to $\widehat\theta_n=\widehat f_h(x)$ and derive the order $\delta_n(x)\asymp\sqrt{nh}\mkern3mu h^p$;
> 3. express the centered random component as a triangular array sum and verify the Lyapunov order $O\lbrace(nh)^{-\eta/2}\rbrace$; and
> 4. compare undersmoothing, MSE optimal smoothing, and more smoothing, explaining why consistency does not by itself validate inference.

### 3.4 Feasible pointwise inference

Let

$$
\widehat s_h^2(x)
=\frac1{n(n-1)}
\sum_{i=1}^n
\lbrace W_{i,h}(x)-\overline W_h(x)\rbrace^2,
\qquad
\overline W_h(x)=\widehat f_h(x).
$$

Under conditions ensuring relative variance consistency,

$$
\frac{\widehat s_h(x)}{s_h(x)}\to_{\mathbb{P}}1.
$$

The feasible statistic retains the same decomposition:

$$
\frac{\widehat f_h(x)-f(x)}{\widehat s_h(x)} =
\frac{s_h(x)}{\widehat s_h(x)}Z_n(x)
+
\frac{b_h(x)}{\widehat s_h(x)}.
$$

If $\delta_n(x)\to\delta(x)$, then this statistic converges to $\mathcal N\lbrace\delta(x),1\rbrace$. Therefore the conventional interval

$$
\mathcal C_n(x) =
\left[
\widehat f_h(x)-z_{1-\alpha/2}\widehat s_h(x),
\widehat f_h(x)+z_{1-\alpha/2}\widehat s_h(x)
\right]
$$

has limiting coverage

$$
\Phi\lbrace z_{1-\alpha/2}-\delta(x)\rbrace -
\Phi\lbrace-z_{1-\alpha/2}-\delta(x)\rbrace.
$$

This equals the nominal level $1-\alpha$ when $\delta(x)=0$. A nonzero limiting signal to noise ratio shifts the distribution and produces asymptotic coverage distortion.

### 3.5 Bias correction as a change of center

A different strategy estimates the leading bias and subtracts it. If $\widehat B_f(x)$ estimates $B_f(x)$, define

$$
\widehat f_{\mathrm{bc}}(x) =
\widehat f_h(x)-h^p\widehat B_f(x).
$$

If the bias were known exactly, this change of center would remove the leading signal to noise term. An estimated bias is random, however, so the distribution of $\widehat f_{\mathrm{bc}}(x)$ includes the estimation noise in $\widehat B_f(x)$ and its covariance with $\widehat f_h(x)$. Simply subtracting an estimated bias while retaining the original standard error $\widehat s_h(x)$ is generally invalid.

### 3.6 Robust bias correction with an MSE optimal bandwidth

Consider a second order KDE, for which

$$
B_f(x)=\frac{\mu_2(K)}2f''(x).
$$

Using a twice differentiable pilot kernel $L$ and pilot bandwidth $b$, estimate

$$
f''(x)
\quad\text{by}\quad
\widehat f_b^{(2)}(x) =
\frac1{nb^3}
\sum_{i=1}^n
L^{(2)}\mkern-3mu\left(\frac{X_i-x}{b}\right).
$$

The bias corrected estimator is

$$
\widehat f_{\mathrm{bc}}(x) =
\widehat f_h(x)
-h^2\frac{\mu_2(K)}2
\widehat f_b^{(2)}(x).
$$

It can be written as a single sample average

$$
\widehat f_{\mathrm{bc}}(x)
=\frac1n\sum_{i=1}^nU_{i,h,b}(x),
$$

where

$$
U_{i,h,b}(x) =
\frac1hK\mkern-3mu\left(\frac{X_i-x}{h}\right)
-h^2\frac{\mu_2(K)}{2b^3}
L^{(2)}\mkern-3mu\left(\frac{X_i-x}{b}\right).
$$

This representation displays the correct variance:

$$
\mathbb V\lbrace\widehat f_{\mathrm{bc}}(x)\rbrace =
\mathbb V\lbrace\widehat f_h(x)\rbrace
+h^4\frac{\mu_2(K)^2}{4}
\mathbb V\lbrace\widehat f_b^{(2)}(x)\rbrace
-h^2\mu_2(K)
\mathrm{Cov}\lbrace\widehat f_h(x),\widehat f_b^{(2)}(x)\rbrace.
$$

Equivalently, estimate this variance by taking the empirical variance of $U_{i,h,b}(x)$ and dividing by $n$. The resulting standard error $\widehat s_{\mathrm{rbc}}(x)$ accounts for both the randomness of the estimated bias and its covariance with the original KDE.

Suppose $h\asymp n^{-1/5}$, the pointwise MSE optimal order for a second order KDE, and $b$ is of the same order as $h$. With enough smoothness, the leading $h^2$ bias is removed and the remaining bias is of smaller order than the corrected standard error. Then

$$
\frac{\widehat f_{\mathrm{bc}}(x)-f(x)}
{\widehat s_{\mathrm{rbc}}(x)}
\rightsquigarrow
\mathcal N(0,1),
$$

which motivates the robust bias corrected interval

$$
\left[
\widehat f_{\mathrm{bc}}(x)
\pm
z_{1-\alpha/2}\widehat s_{\mathrm{rbc}}(x)
\right].
$$

The word **robust** refers to studentizing the bias corrected estimator with its own full variability, not to robustness against outliers or arbitrary misspecification. Robust bias correction can support inference with an MSE optimal bandwidth, but the conclusion still requires adequate smoothness, appropriate relative behavior of $h$ and $b$, consistent variance estimation, and negligible remaining bias. It does not make every data selected MSE bandwidth automatically inference valid.

### 3.7 Pointwise intervals are not simultaneous bands

For one fixed $x$, an interval $\mathcal C_n(x)$ may satisfy

$$
\mathbb{P}_f\lbrace f(x)\in\mathcal C_n(x)\rbrace\to1-\alpha.
$$

A simultaneous band over $\mathcal X$ instead requires

$$
\mathbb{P}_f\lbrace f(x)\in\mathcal C_n(x)\ \text{for every }x\in\mathcal X\rbrace
\to1-\alpha.
$$

Uniform validity over a function class $\mathcal F$ is a second, distinct uniformity:

$$
\liminf_{n\to\infty}
\inf_{f\in\mathcal F}
\mathbb{P}_f\lbrace f(x)\in\mathcal C_n(x)\rbrace
\geq1-\alpha.
$$

Both forms may be combined, but neither follows from a pointwise CLT. Worst case coverage is governed by an infimum, not a supremum.

### Checkpoint 2

1. Which part of the random plus signal to noise decomposition is exact, and which parts require asymptotic approximations?
2. Why do $b_n\to0$ and $s_n\to0$ not imply $b_n/s_n\to0$?
3. Which condition produces the Gaussian random component, and which additional condition removes the KDE's smoothing signal?
4. Why does an MSE optimal bandwidth generally produce a shifted limiting distribution for a statistic centered at $f(x)$?
5. In the second order KDE example, what are the corrected center and corrected standard error, and why is the original KDE standard error insufficient?
6. What is the difference among pointwise coverage, simultaneous coverage over $x$, and uniform validity over $f\in\mathcal F$?

<a id="w11-stop-4"></a>

## 4. Local polynomial regression

### 4.1 Regression target and estimator

Let $(Y_i,X_i)$ be i.i.d. and write

$$
Y_i=\mu(X_i)+\varepsilon_i,
\qquad
\mu(x)=\mathbb E[Y_i\mid X_i=x],
\qquad
\mathbb E[\varepsilon_i\mid X_i]=0.
$$

For an integer $q\geq0$, let $K$ be a nonnegative local polynomial weighting kernel and define

$$
r_q(u)=(1,u,\ldots,u^q)'.
$$

At an evaluation point $x$, the local polynomial coefficient solves

$$
\widehat b(x)
\in
\arg\min_{b\in\mathbb R^{q+1}}
\sum_{i=1}^n
\left[
Y_i-r_q\mkern-3mu\left(\frac{X_i-x}{h}\right)'b
\right]^2
K\mkern-3mu\left(\frac{X_i-x}{h}\right).
$$

The regression function estimator is

$$
\widehat\mu(x)=e_0'\widehat b(x),
$$

and, for $\nu\leq q$, the derivative estimator is

$$
\widehat\mu^{(\nu)}(x)
=\frac{\nu!}{h^\nu}e_\nu'\widehat b(x).
$$

Here $e_\nu$ is the coordinate vector selecting the coefficient on $u^\nu$.

The case $q=0$ is local constant or Nadaraya--Watson regression. The case $q=1$ is local linear regression.

### 4.2 Weighted least squares and linear smoother form

Let $R_x$ have rows $r_q((X_i-x)/h)'$, let $W_x$ be diagonal with entries $K((X_i-x)/h)$, and let $Y=(Y_1,\ldots,Y_n)'$. Whenever the local design matrix is nonsingular,

$$
\widehat b(x)
=(R_x'W_xR_x)^{-1}R_x'W_xY.
$$

Thus

$$
\widehat\mu(x)
=\sum_{i=1}^n\ell_i(x)Y_i,
$$

where

$$
\ell_i(x) =
e_0'(R_x'W_xR_x)^{-1}
r_q\mkern-3mu\left(\frac{X_i-x}{h}\right)
K\mkern-3mu\left(\frac{X_i-x}{h}\right).
$$

Conditional on the design points, the estimator is linear in the outcomes. Its weights generally need not be nonnegative.

### 4.3 Polynomial reproduction and boundary behavior

If $\mu$ is a polynomial of degree at most $q$, the noiseless local polynomial fit reproduces $\mu(x)$ exactly whenever the local design matrix is nonsingular. This reproduction property implies moment conditions for the weights. In particular, local linear regression satisfies

$$
\sum_i\ell_i(x)=1,
\qquad
\sum_i\ell_i(x)(X_i-x)=0.
$$

The second condition cancels the leading local slope contribution even when the neighborhood is asymmetric at a boundary. This is why local linear regression has substantially better boundary behavior than local constant regression.

For local linear estimation of a twice differentiable regression function, a representative interior or boundary expansion has squared bias of order $h^4$ and variance of order $1/(nh)$ under regularity conditions. Constants and parity effects depend on the polynomial degree, derivative, kernel, evaluation point, and design density.

> [!IMPORTANT]
> **Board work 3 — Local polynomial weights and boundary adaptation**
>
> 1. derive the linear smoother weights from weighted least squares;
> 2. prove the two local linear reproduction identities;
> 3. explain geometrically why the identities remain useful at a support boundary; and
> 4. compare local constant and local linear bias without assuming symmetric design points.

> [!TIP]
> **AI interaction 2 — The wrong residual restriction**
>
> Copy the prompt below into an AI interface and audit the response.

```text
A regression slide writes Y=mu(X)+epsilon and then assumes E[Y|X]=0.
Audit and repair the model. Derive local polynomial weighted least squares,
show that the fitted value is a linear smoother, and explain which weight
moment makes local linear regression boundary adaptive.
```

**Audit question:** Does the response impose the zero conditional mean on $\varepsilon$, rather than on $Y$?

### 4.4 Conditional bias and variance: exact formulas and asymptotic orders

Let $\mathbf X=(X_1,\ldots,X_n)$ and define $\sigma^2(t)=\mathbb V(Y_i\mid X_i=t)$. Because $\widehat\mu(x)=\sum_i\ell_i(x)Y_i$ and the weights are fixed conditional on $\mathbf X$, the finite sample conditional bias is exactly

$$
\mathrm{Bias}\lbrace\widehat\mu(x)\mid\mathbf X\rbrace =
\sum_{i=1}^n\ell_i(x)\mu(X_i)-\mu(x),
$$

and conditional independence gives the exact conditional variance

$$
\mathbb V\lbrace\widehat\mu(x)\mid\mathbf X\rbrace =
\sum_{i=1}^n\ell_i(x)^2\sigma^2(X_i).
$$

Consequently, the finite sample conditional MSE is the squared conditional bias plus the conditional variance. These identities require no Taylor expansion and make the two sources of error transparent: how well the local polynomial approximates $\mu$ and how much conditional noise the smoother weights transmit.

For a one dimensional local linear estimator at a fixed interior point, suppose that $\mu$ is twice continuously differentiable, the design density $f_X(x)$ is positive and continuous, $\sigma^2(x)$ is continuous, the kernel is symmetric with $\int K=1$, and $h\to0$ with $nh\to\infty$. Then the corresponding conditional expansions are

$$
\mathrm{Bias}\lbrace\widehat\mu(x)\mid\mathbf X\rbrace =
\frac{h^2\mu_2(K)}2\mu''(x)
+o_{\mathbb{P}}(h^2)
$$

and

$$
\mathbb V\lbrace\widehat\mu(x)\mid\mathbf X\rbrace =
\frac{\sigma^2(x)R(K)}{nhf_X(x)}
\lbrace1+o_{\mathbb{P}}(1)\rbrace.
$$

Thus the pointwise conditional MSE has leading order $h^4+(nh)^{-1}$, giving $h_{\mathrm{MSE}}\asymp n^{-1/5}$ and optimal MSE order $n^{-4/5}$. At a boundary, local linear regression retains these orders under standard regularity conditions, although the constants change. The exact formulas above remain valid at every fixed $n$; the asymptotic formulas summarize their leading behavior as the local design fills in.

### 4.5 Pointwise inference: conditional CLT plus smoothing bias

Write

$$
\mu_{\mathbf X,h}(x)=
\mathbb E[\widehat\mu_h(x)\mid\mathbf X],
\qquad
s_{\mathbf X,h}^2(x)=
\mathbb V[\widehat\mu_h(x)\mid\mathbf X]
=\sum_{i=1}^n\ell_i(x)^2\sigma^2(X_i).
$$

Adding and subtracting the conditional mean gives the exact local polynomial analogue of the KDE decomposition:

$$
\frac{\widehat\mu_h(x)-\mu(x)}{s_{\mathbf X,h}(x)} =
\underbrace{
\frac{\sum_{i=1}^n\ell_i(x)\varepsilon_i}{s_{\mathbf X,h}(x)}
}_{\text{random component}}
+
\underbrace{
\frac{\mu_{\mathbf X,h}(x)-\mu(x)}{s_{\mathbf X,h}(x)}
}_{\text{smoothing bias in standard error units}}.
$$

Conditional on $\mathbf X$, the random component is a weighted sum of independent mean zero errors. A conditional Lindeberg condition, for example

$$
\sum_{i=1}^n
\mathbb E\left[
\left.
\left|
\frac{\ell_i(x)\varepsilon_i}{s_{\mathbf X,h}(x)}
\right|^{2+\eta}
\right|\mathbf X
\right]
\to_{\mathbb{P}}0
$$

for some $\eta>0$, gives a conditional standard Normal limit. This condition is plausible when the effective local sample size diverges, no local observation has dominant leverage, and conditional error moments are controlled.

The least squares interpretation is clearest if $h$ is held fixed. Define the population weighted projection coefficient

$$
b_{q,h}(x)
\in
\arg\min_{b\in\mathbb R^{q+1}}
\mathbb E\left[
\left\lbrace
Y-r_q\mkern-3mu\left(\frac{X-x}{h}\right)'b
\right\rbrace^2
K\mkern-3mu\left(\frac{X-x}{h}\right)
\right]
$$

and put $\mu_{q,h}(x)=e_0'b_{q,h}(x)$. Ordinary weighted least squares asymptotics apply to $\widehat b(x)$ around $b_{q,h}(x)$ under fixed dimension regularity conditions. The resulting interval targets the finite bandwidth projection $\mu_{q,h}(x)$, not automatically the regression function $\mu(x)$. Letting $h\to0$ turns the problem into a triangular array and makes the difference $\mu_{q,h}(x)-\mu(x)$ the smoothing bias that must be controlled.

For one dimensional local linear level estimation at an interior point,

$$
\mu_{\mathbf X,h}(x)-\mu(x)=O_{\mathbb{P}}(h^2),
\qquad
s_{\mathbf X,h}(x)\asymp_{\mathbb{P}}(nh)^{-1/2}.
$$

The bias to standard error ratio therefore has order $\sqrt{nh}\mkern3mu h^2$. Conventional inference centered at $\mu(x)$ requires $nh^5\to0$. The MSE optimal order $h\asymp n^{-1/5}$ generally leaves a nonzero first order bias and hence does not justify the conventional interval.

### 4.6 Feasible inference from the OLS sandwich

The handwritten OLS calculation applies directly. Let $\widehat\varepsilon_i(x)=Y_i-r_q((X_i-x)/h)'\widehat b(x)$ and define

$$
\widehat\Sigma_x=
\mathrm{diag}\lbrace
\widehat\varepsilon_1(x)^2,\ldots,
\widehat\varepsilon_n(x)^2
\rbrace.
$$

The same heteroskedasticity robust weighted least squares sandwich applies to the local fit. Define

$$
\widehat V_b(x) =
(R_x'W_xR_x)^{-1}
R_x'W_x\widehat\Sigma_xW_xR_x
(R_x'W_xR_x)^{-1}.
$$

Thus

$$
\widehat s_{\mathrm{lp}}^2(x) =
e_0'
\widehat V_b(x)
e_0 =
\sum_{i=1}^n\ell_i(x)^2\widehat\varepsilon_i(x)^2.
$$

For a derivative, premultiply and postmultiply by $(\nu!/h^\nu)e_\nu$ instead of $e_0$. With fixed $h$, this is the Eicker--Huber--White sandwich for the population weighted projection, including local approximation error in its regression residual. Along a sequence $h\to0$, the approximation error vanishes locally and the residuals estimate the conditional noise. Under vanishing maximum leverage and the same local moment and design conditions used for the CLT, the feasible standard error is ratio consistent. HC leverage adjustments can be used when the effective local sample is modest, but they do not remove smoothing bias.

For local linear level estimation, if $nh\to\infty$, $nh^5\to0$, and the CLT and variance consistency conditions hold, the feasible undersmoothed interval is

$$
\left[
\widehat\mu_h(x)
\pm
z_{1-\alpha/2}\widehat s_{\mathrm{lp}}(x)
\right].
$$

### 4.7 RBC inference for local polynomial regression

RBC estimates the leading smoothing bias, subtracts it, and studentizes the corrected estimator with its own variance. For a simple interior local linear illustration, fit a local quadratic with bandwidth $b$ and write its second derivative estimator as

$$
\widehat\mu_{2,b}^{(2)}(x)=
\sum_{i=1}^n\ell_{i,2,b}^{(2)}(x)Y_i.
$$

The leading bias corrected center is

$$
\widehat\mu_{\mathrm{bc}}(x) =
\widehat\mu_{1,h}(x) -
h^2\frac{\mu_2(K)}2
\widehat\mu_{2,b}^{(2)}(x) =
\sum_{i=1}^n\ell_i^{\mathrm{bc}}(x)Y_i,
$$

where

$$
\ell_i^{\mathrm{bc}}(x) =
\ell_{i,1,h}(x) -
h^2\frac{\mu_2(K)}2
\ell_{i,2,b}^{(2)}(x).
$$

The interior constant $\mu_2(K)/2$ is replaced by the appropriate local polynomial moment matrix expression at a boundary or for other degrees and derivatives. Let $\widehat\varepsilon_{i,+}(x)$ be residuals from the auxiliary higher degree fit. The robust bias corrected variance estimator is

$$
\widehat s_{\mathrm{rbc}}^2(x) =
\sum_{i=1}^n
\ell_i^{\mathrm{bc}}(x)^2
\widehat\varepsilon_{i,+}(x)^2.
$$

Because it uses the combined weights, this variance includes the variability of the estimated bias and its covariance with the original local linear estimator. Subtracting the estimated bias while retaining $\widehat s_{\mathrm{lp}}(x)$ would omit both terms and is not RBC inference.

If $h$ has its MSE optimal order, $h/b\to\rho\in(0,\infty)$, the regression function is sufficiently smooth, the remaining bias is negligible relative to $\widehat s_{\mathrm{rbc}}(x)$, and the conditional CLT and variance consistency conditions hold, then

$$
\frac{\widehat\mu_{\mathrm{bc}}(x)-\mu(x)}
{\widehat s_{\mathrm{rbc}}(x)}
\rightsquigarrow
\mathcal N(0,1).
$$

The corresponding pointwise interval is

$$
\left[
\widehat\mu_{\mathrm{bc}}(x)
\pm
z_{1-\alpha/2}\widehat s_{\mathrm{rbc}}(x)
\right].
$$

> [!IMPORTANT]
> **Board work 4 — Local polynomial inference from OLS to RBC**
>
> 1. Decompose the local polynomial statistic into its conditional random component and smoothing bias in standard error units.
> 2. Explain what fixed $h$ OLS inference targets and why $h\to0$ introduces a triangular array and a centering problem.
> 3. Derive the conditional sandwich variance and its residual based feasible estimator.
> 4. For local linear regression, derive the order $\sqrt{nh}\mkern3mu h^2$ and compare undersmoothing with MSE optimal smoothing.
> 5. Express the RBC estimator through combined local linear and local quadratic weights and derive its robust standard error.
> 6. Identify the variance terms omitted if the original local polynomial standard error is retained after bias correction.

<a id="w11-stop-5"></a>

## 5. Comparing KDE and local polynomial smoothing

### 5.1 A common pointwise bias--variance template

For $X\in\mathbb R^d$ and an isotropic bandwidth $h$, both density estimation and regression smoothing use about $nh^d$ observations in a local neighborhood. Let $p$ be the leading KDE bias exponent and let $r_\mu$ be the first nonzero local polynomial bias exponent after polynomial reproduction and any kernel symmetry are taken into account. For level estimation, the main pointwise orders are:

| Quantity | KDE for $f(x)$ | Local polynomial for $\mu(x)$ |
|---|---:|---:|
| Leading bias | $O(h^p)$ | $O_{\mathbb{P}}(h^{r_\mu})$ |
| Leading variance | $O((nh^d)^{-1})$ | $O_{\mathbb{P}}((nh^d)^{-1})$ conditionally on the design |
| Pointwise MSE | $O(h^{2p}+(nh^d)^{-1})$ | $O_{\mathbb{P}}(h^{2r_\mu}+(nh^d)^{-1})$ |
| MSE optimal bandwidth | $n^{-1/(2p+d)}$ | $n^{-1/(2r_\mu+d)}$ |
| Optimal MSE | $n^{-2p/(2p+d)}$ | $n^{-2r_\mu/(2r_\mu+d)}$ |

For a second order KDE and local linear level estimation, $p=r_\mu=2$. Both therefore have bandwidth order $n^{-1/(4+d)}$ and MSE order $n^{-4/(4+d)}$; in one dimension these become $n^{-1/5}$ and $n^{-4/5}$. The rates coincide, but the constants do not: KDE variance depends on $f(x)$ and the kernel, whereas regression variance also depends on the conditional variance $\sigma^2(x)$, the design density $f_X(x)$, and the local polynomial equivalent kernel. For local polynomials, $r_\mu$ depends on the polynomial degree, derivative, evaluation point, and parity; local linear level estimation has $r_\mu=2$ both in the interior and at a regular boundary.

### 5.2 Dimension and derivative order

For a KDE estimating a density derivative of total order $s$, sufficient smoothness and a kernel of order $p$ give

$$
\mathrm{Bias}=O(h^p),
\qquad
\mathbb V=O\mkern-3mu\left(\frac1{nh^{d+2s}}\right),
$$

and hence

$$
\mathrm{MSE}
=O\mkern-3mu\left(h^{2p}+\frac1{nh^{d+2s}}\right),
\qquad
h_{\mathrm{MSE},s}\asymp n^{-1/(2p+d+2s)}.
$$

For a local polynomial of degree $q$ estimating a regression derivative of total order $\nu$, let $r_{q,\nu}$ denote its leading bias exponent. Then

$$
\mathrm{Bias}=O_{\mathbb{P}}(h^{r_{q,\nu}}),
\qquad
\mathbb V(\widehat\mu^{(\nu)}\mid\mathbf X)
=O_{\mathbb{P}}\mkern-3mu\left(\frac1{nh^{d+2\nu}}\right),
$$

so

$$
\mathrm{MSE}
=O_{\mathbb{P}}\mkern-3mu\left(h^{2r_{q,\nu}}+\frac1{nh^{d+2\nu}}\right),
\qquad
h_{\mathrm{MSE},\nu}
\asymp n^{-1/(2r_{q,\nu}+d+2\nu)}.
$$

Generically $r_{q,\nu}=q+1-\nu$; at an interior point, kernel symmetry and parity may eliminate that term and increase the exponent by one. Density derivatives require the corresponding additional derivatives of $f$, while regression derivatives require additional smoothness of $\mu$. In both columns, adding dimensions shrinks the effective local sample, and estimating one more derivative adds two powers of $h$ to the variance denominator. Local polynomials estimate regression derivatives directly through fitted coefficients rather than by differentiating a noisy fitted curve after the fact.

The table and formulas are pointwise. Supremum norm rates generally acquire a logarithmic cost, and simultaneous confidence bands require a uniform stochastic approximation rather than merely a pointwise CLT.

### Checkpoint 3

1. What is the common effective local sample size for KDE and local polynomial level estimation in $d$ dimensions?
2. For second order KDE and local linear regression, what are the MSE optimal bandwidth and MSE orders when $d=1$?
3. Why does estimating an $s$-th density derivative or a $\nu$-th regression derivative add twice the derivative order to the variance exponent?
4. Which additional theorem is needed to turn a supremum norm rate into a confidence band?
5. With fixed $h$, what population object does ordinary weighted least squares inference target?
6. Why must an RBC local polynomial interval use the combined main and bias correction weights in its standard error?

## 6. Recap

- A KDE is a local average whose expectation is exactly a convolution.
- Interior Taylor expansion gives bias order $h^p$ for a kernel of order $p$.
- Pointwise variance is of order $1/(nh)$, not $1/(nh^2)$.
- Estimation optimal bandwidth balances squared bias and variance.
- Minimaxity compares the best worst case risk over a specified function class and requires matching upper and lower bounds for a stated loss.
- A standardized estimation error equals a centered random component plus bias measured in standard error units; consistency does not make this signal to noise ratio negligible.
- Conventional pointwise inference requires bias control: undersmoothing makes bias negligible, while robust bias correction changes both the center and the standard error and can accommodate an MSE optimal bandwidth under additional conditions.
- A pointwise interval, a simultaneous band over $x$, and uniform validity over a function class are different claims.
- Local polynomial regression is weighted least squares and a linear smoother; conditional on the design, its finite sample bias and variance follow exactly from its weights.
- Polynomial reproduction explains the boundary advantage of local linear regression, while Taylor expansion and local design regularity turn the exact conditional formulas into leading asymptotic orders.
- Local polynomial inference repeats the KDE logic: a conditional CLT controls the random component, while smoothing bias must be negligible or explicitly corrected.
- Fixed bandwidth OLS inference targets a local weighted projection; feasible inference uses the weighted least squares sandwich, and RBC uses the variance of the combined main and bias correction weights.
- KDE and local polynomial regression share the same bias--variance template, but their constants and targets differ.
- Dimension and derivative order slow both density and regression smoothing.

## Notation introduced this week

- $K$: kernel function; $\mu_p(K)$: $p$-th kernel moment.
- $R(K)=\int K^2$: kernel roughness constant.
- $h$: bandwidth.
- $\widehat f_h(x)$: kernel density estimator.
- $\mathcal H(\beta,L)$ and $R_n^\star$: a Hölder type smoothness class and its minimax risk under a specified loss.
- $b_h(x)$, $s_h(x)$, and $\delta_n(x)=b_h(x)/s_h(x)$: pointwise bias, standard deviation, and bias to standard error ratio.
- $L$ and $b$: pilot kernel and pilot bandwidth used to estimate the leading bias.
- $\widehat f_{\mathrm{bc}}(x)$ and $\widehat s_{\mathrm{rbc}}(x)$: bias corrected center and robust standard error including bias estimation variability.
- $\mu(x)=\mathbb E[Y\mid X=x]$, $\sigma^2(x)=\mathbb V(Y\mid X=x)$, and $f_X(x)$: conditional mean, conditional variance, and design density.
- $q$: local polynomial degree; $\nu$: regression derivative order.
- $r_\mu$ and $r_{q,\nu}$: leading local polynomial bias exponents for level and derivative estimation.
- $\ell_i(x)$: local polynomial linear smoother weight.
- $\widehat s_{\mathrm{lp}}(x)$: feasible local polynomial OLS standard error.
- $\ell_i^{\mathrm{bc}}(x)$ and $\widehat s_{\mathrm{rbc}}(x)$: combined RBC weight and its robust standard error.
- $\mathcal X$: evaluation set; $\mathcal F$: function class.

## References

- Paul P. B. Eggermont and Vincent N. LaRiccia, *Maximum Penalized Likelihood Estimation, Volume I: Density Estimation*, Springer, 2001.
- Paul P. B. Eggermont and Vincent N. LaRiccia, *Maximum Penalized Likelihood Estimation, Volume II: Regression*, Springer, 2009.
- László Györfi et al., *A Distribution-Free Theory of Nonparametric Regression*, Springer, 2002.
- Martin J. Wainwright, *High-Dimensional Statistics: A Non-Asymptotic Viewpoint*, Cambridge University Press, 2019.
