# ORF 524 Practice Module 11: Solutions

**Status:** Released solutions  
**Last updated:** October 2, 2026  
**Practice module:** [Module 11](../topic-11-kernel-nonparametric-methods.md)  
**Primary chapter:** [Week 11](../../lectures/week-11.md)

Attempt each problem before consulting its solution. After studying a solution, close it and reconstruct the argument, including the assumptions used at each step.

<a id="solution-map"></a>

## Solution map

| Problem | Solution |
|---|---|
| [1. A fourth order kernel density estimator](../topic-11-kernel-nonparametric-methods.md#problem-1) | [Solution](#problem-1) |
| [2. Bandwidth choice for pointwise inference](../topic-11-kernel-nonparametric-methods.md#problem-2) | [Solution](#problem-2) |
| [3. Local polynomial regression and boundary adaptation](../topic-11-kernel-nonparametric-methods.md#problem-3) | [Solution](#problem-3) |
| [4. Density and local polynomial rates under dimension and derivatives](../topic-11-kernel-nonparametric-methods.md#problem-4) | [Solution](#problem-4) |
| [5. Higher order kernels and positivity](../topic-11-kernel-nonparametric-methods.md#problem-5) | [Solution](#problem-5) |
| [6. Three different meanings of uniform](../topic-11-kernel-nonparametric-methods.md#problem-6) | [Solution](#problem-6) |

<a id="problem-1"></a>

## Problem 1. A fourth order kernel density estimator

[Problem statement](../topic-11-kernel-nonparametric-methods.md#problem-1)

The kernel is even, so its first and third moments vanish. Expanding on $[-1,1]$ gives

$$
K(u)=\frac1{32}(45-150u^2+105u^4).
$$

Direct integration yields

$$
\int_{-1}^1K(u)\mkern3mu du=1,
\qquad
\int_{-1}^1u^2K(u)\mkern3mu du=0,
$$

and

$$
\mu_4(K)
=\int_{-1}^1u^4K(u)\mkern3mu du
=-\frac1{21}.
$$

Squaring the polynomial and integrating gives $R(K)=5/4$.

The convolution identity and a fourth order Taylor expansion give

$$
\mathbb E[\widehat f_h(x)]-f(x) =
h^4\frac{\mu_4(K)}{4!}f^{(4)}(x)+o(h^4) =
-\frac{h^4}{504}f^{(4)}(x)+o(h^4).
$$

Also,

$$
\mathbb V[\widehat f_h(x)] =
\frac{R(K)f(x)}{nh}+o\lbrace(nh)^{-1}\rbrace =
\frac{5f(x)}{4nh}+o\lbrace(nh)^{-1}\rbrace.
$$

Thus the MSE order is $h^8+(nh)^{-1}$. Balancing terms gives

$$
h_{\mathrm{MSE}}\asymp n^{-1/9},
\qquad
\mathrm{MSE}_{\mathrm{opt}}\asymp n^{-8/9}.
$$

More precisely, if $f^{(4)}(x)\neq0$, write $B_f(x)=\mu_4(K)f^{(4)}(x)/4!=-f^{(4)}(x)/504$. Minimizing the leading approximation $B_f(x)^2h^8+R(K)f(x)/(nh)$ gives

$$
h_{\mathrm{MSE}}
=\left\lbrace\frac{R(K)f(x)}{8B_f(x)^2n}\right\rbrace^{1/9}
=\left\lbrace\frac{5\mkern3mu504^2 f(x)}{32 f^{(4)}(x)^2n}\right\rbrace^{1/9}.
$$

If $f^{(4)}(x)=0$, this leading term optimizer is not informative and a higher order expansion is needed.

For a second order kernel, the corresponding orders are $h_{\mathrm{MSE}}\asymp n^{-1/5}$ and $\mathrm{MSE}_{\mathrm{opt}}\asymp n^{-4/5}$. The fourth order improvement requires four derivatives and control of the Taylor remainder; without that smoothness the formal cancellation does not deliver fourth order bias.

[Back to the solution map](#solution-map)

<a id="problem-2"></a>

## Problem 2. Bandwidth choice for pointwise inference

[Problem statement](../topic-11-kernel-nonparametric-methods.md#problem-2)

Let

$$
s_h^2(x)=\mathbb V\lbrace\widehat f_h(x)\rbrace,
\qquad
b_h(x)=\mathbb E[\widehat f_h(x)]-f(x).
$$

Adding and subtracting the exact expectation gives

$$
\frac{\widehat f_h(x)-f(x)}{s_h(x)} =
\frac{\widehat f_h(x)-\mathbb E[\widehat f_h(x)]}{s_h(x)}
+
\frac{b_h(x)}{s_h(x)}.
$$

The first term is the random component and the second is the bias to standard error signal to noise ratio. Consistency requires both the unscaled random error and bias to vanish, but it does not require their ratio to vanish.

Let

$$
\overline W_h=\frac1n\sum_iW_{i,h}=\widehat f_h(x).
$$

Because $s_h^2(x)=\mathbb V(W_{i,h})/n$, the random component is

$$
\sum_{i=1}^n\xi_{i,n}(x),
\qquad
\xi_{i,n}(x) =
\frac{W_{i,h}(x)-\mathbb E W_{i,h}(x)}
{\sqrt{n\mathbb V\lbrace W_{i,h}(x)\rbrace}}.
$$

The summands have mean zero and their variances sum to one. Moreover,

$$
\sum_{i=1}^n\mathbb E|\xi_{i,n}(x)|^{2+\eta} =
\frac{n\mkern3mu O(h^{-1-\eta})}
{\lbrace n\mathbb V(W_{i,h})\rbrace^{1+\eta/2}}
=O\lbrace(nh)^{-\eta/2}\rbrace
\to0,
$$

where $\mathbb V(W_{i,h})\asymp h^{-1}$. Lyapunov's theorem therefore gives

$$
\frac{\widehat f_h(x)-\mathbb E[\widehat f_h(x)]}{s_h(x)}
\rightsquigarrow
\mathcal N(0,1).
$$

A consistent estimator of the variance of the sample mean is

$$
\widehat{\mathbb V}\lbrace\widehat f_h(x)\rbrace =
\frac1{n(n-1)}
\sum_{i=1}^n(W_{i,h}-\overline W_h)^2.
$$

The sample variance estimates $\mathbb V(W_{i,h})$; division by $n$ converts it to the variance of $\mathbb E_nW_{i,h}$.

For example, suppose $\mathbb E|W_{i,h}-\mathbb EW_{i,h}|^4=O(h^{-3})$. Since $\mathbb V(W_{i,h})\asymp h^{-1}$,

$$
\frac{
\mathbb V\mkern-3mu\left(
\mathbb{E}_n[(W_{i,h}-\mathbb EW_{i,h})^2]
\right)
}{\mathbb V(W_{i,h})^2}
=O\mkern-3mu\left(\frac1{nh}\right)
\to0.
$$

Chebyshev's inequality and the negligible sample centering correction then give

$$
\frac{\widehat s_h^2(x)}{s_h^2(x)}
\to_{\mathbb{P}}1.
$$

For $h=n^{-a}$, $h\to0$ requires $a>0$, and $nh\to\infty$ requires $a<1$. Moreover,

$$
\sqrt{nh}h^p =
n^{1/2-a(p+1/2)}
\to0
$$

if and only if

$$
a>\frac1{2p+1}.
$$

Thus a simple admissible range is $1/(2p+1)<a<1$. The MSE optimal exponent is exactly $1/(2p+1)$, for which the standardized bias has constant order rather than vanishing.

Indeed, the bias and variance expansions imply

$$
\frac{b_h(x)}{s_h(x)} =
\frac{B_f(x)}{\sqrt{R(K)f(x)}}
\sqrt{nh^{2p+1}}
+o\mkern-3mu\left(\sqrt{nh^{2p+1}}\right).
$$

If $h\sim c n^{-1/(2p+1)}$ and $B_f(x)\neq0$, then

$$
\frac{b_h(x)}{s_h(x)}
\to
\delta(x) =
\frac{c^{p+1/2}B_f(x)}{\sqrt{R(K)f(x)}}.
$$

Consequently,

$$
\frac{\widehat f_h(x)-f(x)}{s_h(x)}
\rightsquigarrow
\mathcal N\lbrace\delta(x),1\rbrace,
$$

not $\mathcal N(0,1)$. If $B_f(x)=0$, the next nonzero bias term must be used instead.

With $\widehat s_h^2$ denoting the variance estimator above, an undersmoothed interval is

$$
\left[
\widehat f_h(x)-z_{1-\alpha/2}\widehat s_h,\mkern3mu
\widehat f_h(x)+z_{1-\alpha/2}\widehat s_h
\right].
$$

Under the stated pointwise CLT, variance consistency, undersmoothing, and regularity conditions, its coverage converges to $1-\alpha$ for the fixed interior point and fixed density. No simultaneous or function class uniform claim follows.

More generally, if $b_h(x)/s_h(x)\to\delta$, the feasible statistic converges to $\mathcal N(\delta,1)$ and the interval's limiting coverage is

$$
\Phi\lbrace z_{1-\alpha/2}-\delta\rbrace -
\Phi\lbrace-z_{1-\alpha/2}-\delta\rbrace.
$$

This reduces to $1-\alpha$ when $\delta=0$.

Let $L$ be a twice differentiable pilot kernel and define

$$
\widehat f_b^{(2)}(x) =
\frac1{nb^3}\sum_{i=1}^n
L^{(2)}\mkern-3mu\left(\frac{X_i-x}{b}\right).
$$

The bias corrected estimator is

$$
\widehat f_{\mathrm{bc}}(x) =
\widehat f_h(x)
-h^2\frac{\mu_2(K)}2\widehat f_b^{(2)}(x) =
\frac1n\sum_{i=1}^nU_{i,h,b}(x),
$$

where

$$
U_{i,h,b}(x) =
\frac1hK\mkern-3mu\left(\frac{X_i-x}{h}\right)
-h^2\frac{\mu_2(K)}{2b^3}
L^{(2)}\mkern-3mu\left(\frac{X_i-x}{b}\right).
$$

Therefore a direct robust variance estimator is

$$
\widehat s_{\mathrm{rbc}}^2(x) =
\frac1{n(n-1)}
\sum_{i=1}^n
\lbrace U_{i,h,b}(x)-\overline U_{h,b}(x)\rbrace^2,
\qquad
\overline U_{h,b}(x)=\widehat f_{\mathrm{bc}}(x).
$$

This empirical variance automatically includes the pilot derivative's variability and its covariance with the original KDE. Algebraically, the population variance equals

$$
\mathbb V(\widehat f_h)
+h^4\frac{\mu_2(K)^2}{4}\mathbb V(\widehat f_b^{(2)})
-h^2\mu_2(K)\mathrm{Cov}(\widehat f_h,\widehat f_b^{(2)}).
$$

If $h\asymp n^{-1/5}$, $b/h$ converges to a positive finite constant, and sufficient smoothness makes the post correction bias negligible relative to $\widehat s_{\mathrm{rbc}}(x)$, then the robust bias corrected statistic is asymptotically standard Normal. The corresponding interval is centered at $\widehat f_{\mathrm{bc}}(x)$ and uses $\widehat s_{\mathrm{rbc}}(x)$. Keeping the original standard error drops terms of the same first order magnitude when $b$ is comparable to $h$ and is generally invalid.

[Back to the solution map](#solution-map)

<a id="problem-3"></a>

## Problem 3. Local polynomial regression and boundary adaptation

[Problem statement](../topic-11-kernel-nonparametric-methods.md#problem-3)

Using the notation of Week 11,

$$
\widehat b(x)=(R_x'W_xR_x)^{-1}R_x'W_xY
$$

and hence

$$
\widehat\mu(x)
=e_0'(R_x'W_xR_x)^{-1}R_x'W_xY
=\sum_i\ell_i(x)Y_i,
$$

where

$$
\ell_i(x)=e_0'(R_x'W_xR_x)^{-1}r_q\mkern-3mu\left(\frac{X_i-x}{h}\right)K\mkern-3mu\left(\frac{X_i-x}{h}\right).
$$

If $\mu(t)$ is a polynomial of degree at most $q$, then $\mu(X_i)=r_q((X_i-x)/h)'b_x$ for a coefficient vector $b_x$ whose intercept is $\mu(x)$. In the absence of noise, the weighted residual sum of squares is zero at $b_x$. Nonsingularity makes this solution unique, so the fitted intercept equals $\mu(x)$.

Apply this to $\mu(t)=1$ and $\mu(t)=t-x$. Linearity of the smoother gives

$$
\sum_i\ell_i(x)=1
$$

and

$$
\sum_i\ell_i(x)(X_i-x)=0.
$$

At $x=0$, the population local constant smooth of $\mu(t)=a+bt$ is

$$
\frac{\int_0^\infty K(t/h)(a+bt)f_X(t)\mkern3mu dt}
{\int_0^\infty K(t/h)f_X(t)\mkern3mu dt}.
$$

After $t=hu$ and continuity of $f_X$ at zero, this equals

$$
a+bh
\frac{\int_0^\infty uK(u)\mkern3mu du}
{\int_0^\infty K(u)\mkern3mu du}
+o(h),
$$

which generally has first order boundary bias. Local linear regression reproduces $a+bt$ exactly, so its fitted value at zero is $a$ at the population smoothing level. Negative weights allow the weighted first moment to be zero despite the one sided design; positivity of weights would prevent that cancellation.

Let $\mathbf X=(X_1,\ldots,X_n)$ and $\sigma^2(t)=\mathbb V(Y_i\mid X_i=t)$. Conditional on $\mathbf X$, the smoother weights are fixed, so

$$
\mathrm{Bias}\lbrace\widehat\mu(x)\mid\mathbf X\rbrace =
\sum_{i=1}^n\ell_i(x)\mu(X_i)-\mu(x)
$$

and

$$
\mathbb V\lbrace\widehat\mu(x)\mid\mathbf X\rbrace =
\sum_{i=1}^n\ell_i(x)^2\sigma^2(X_i).
$$

The exact conditional MSE is therefore

$$
\left\lbrace\sum_{i=1}^n\ell_i(x)\mu(X_i)-\mu(x)\right\rbrace^2
+
\sum_{i=1}^n\ell_i(x)^2\sigma^2(X_i).
$$

For one dimensional local linear regression at an interior point, under the stated smoothness and design conditions,

$$
\mathrm{Bias}\lbrace\widehat\mu(x)\mid\mathbf X\rbrace =
\frac{h^2\mu_2(K)}2\mu''(x)+o_{\mathbb{P}}(h^2)
$$

and

$$
\mathbb V\lbrace\widehat\mu(x)\mid\mathbf X\rbrace =
\frac{\sigma^2(x)R(K)}{nhf_X(x)}\lbrace1+o_{\mathbb{P}}(1)\rbrace.
$$

Thus the conditional MSE has order $h^4+(nh)^{-1}$. Balancing its two terms gives $h_{\mathrm{MSE}}\asymp n^{-1/5}$ and optimal MSE order $n^{-4/5}$.

[Back to the solution map](#solution-map)

<a id="problem-4"></a>

## Problem 4. Density and local polynomial rates under dimension and derivatives

[Problem statement](../topic-11-kernel-nonparametric-methods.md#problem-4)

For a density derivative of total order $s$, the pointwise MSE has order

$$
h^{2p}+\frac1{nh^{d+2s}}.
$$

Balancing its terms gives

$$
h_{\mathrm{MSE},s}\asymp n^{-1/(2p+d+2s)},
\qquad
\mathrm{MSE}_{\mathrm{opt},s}
\asymp n^{-2p/(2p+d+2s)}.
$$

For a local polynomial estimate of a regression derivative of total order $\nu$, the conditional MSE has order

$$
h^{2r_{q,\nu}}+\frac1{nh^{d+2\nu}}.
$$

Therefore

$$
h_{\mathrm{MSE},\nu}\asymp n^{-1/(2r_{q,\nu}+d+2\nu)},
\qquad
\mathrm{MSE}_{\mathrm{opt},\nu}
\asymp n^{-2r_{q,\nu}/(2r_{q,\nu}+d+2\nu)}.
$$

For level estimation with $s=\nu=0$ and $p=r_{q,\nu}=2$, both estimators have

$$
h_{\mathrm{MSE}}\asymp n^{-1/(4+d)},
\qquad
\mathrm{MSE}_{\mathrm{opt}}\asymp n^{-4/(4+d)}.
$$

In one dimension these orders are $n^{-1/5}$ and $n^{-4/5}$. Equal rates do not imply equal constants: the KDE variance depends on $f(x)$ and the kernel, whereas local polynomial regression variance depends on $\sigma^2(x)$, $f_X(x)$, and its equivalent kernel.

Adding one covariate increases the denominator governing the bandwidth and rate by one; increasing either derivative order by one increases the variance exponent by two. For the KDE bias, kernel moment cancellations still eliminate powers below $p$ in the convolution expansion for $f^{(s)}$, provided $f$ has at least $p+s$ derivatives. Thus the leading derivative estimation bias is $h^p$, not $h^{p-s}$. For local polynomials, polynomial reproduction gives the generic exponent $r_{q,\nu}=q+1-\nu$, while symmetry and parity can eliminate that term at an interior point and increase the exponent by one; the corresponding derivatives of $\mu$ must exist.

These balances establish estimator specific upper bound rate heuristics. For one dimensional density estimation with $s=0$ and smoothness $p=\beta$, the resulting squared error order $n^{-2\beta/(2\beta+1)}$ matches the standard pointwise and squared $L_2$ minimax benchmark. That match does not prove minimaxity: after defining the function class and loss, one must still show that every estimator has worst case risk bounded below by a constant multiple of the same rate. Nor does the pointwise calculation establish the $L_\infty$ benchmark, whose stochastic term and optimal rate include the logarithmic penalty from controlling all evaluation points simultaneously.

[Back to the solution map](#solution-map)

<a id="problem-5"></a>

## Problem 5. Higher order kernels and positivity

[Problem statement](../topic-11-kernel-nonparametric-methods.md#problem-5)

Differentiating the polynomial shows that the interior critical points include $u=\pm\sqrt{5/7}$, where

$$
K(u)=-\frac{15}{56}.
$$

Thus the estimator may be negative. Since $f(x)\geq0$, projection of any real number $a$ onto $[0,\infty)$ cannot increase its distance to $f(x)$. Therefore

$$
|\max\lbrace\widehat f_h(x),0\rbrace-f(x)|
\leq
|\widehat f_h(x)-f(x)|
$$

pointwise, and the squared error inequality follows. Truncation does not preserve the integral: $\int\max\lbrace\widehat f_h,0\rbrace$ need not equal one. Renormalization changes the estimate at every point, so its bias and risk require a separate analysis.

[Back to the solution map](#solution-map)

<a id="problem-6"></a>

## Problem 6. Three different meanings of uniform

[Problem statement](../topic-11-kernel-nonparametric-methods.md#problem-6)

For fixed $x$ and $f$:

$$
\mathbb{P}_f\lbrace f(x)\in\mathcal C_n(x)\rbrace\to1-\alpha.
$$

Simultaneous coverage for fixed $f$ is

$$
\mathbb{P}_f\lbrace f(x)\in\mathcal C_n(x)\ \forall x\in\mathcal X\rbrace
\to1-\alpha.
$$

Pointwise validity at $x$, uniformly over $\mathcal F$, can be stated as

$$
\liminf_n\inf_{f\in\mathcal F}
\mathbb{P}_f\lbrace f(x)\in\mathcal C_n(x)\rbrace
\geq1-\alpha
$$

for the fixed $x$. Combining both uniformities gives

$$
\liminf_n\inf_{f\in\mathcal F}
\mathbb{P}_f\lbrace f(x)\in\mathcal C_n(x)\ \forall x\in\mathcal X\rbrace
\geq1-\alpha.
$$

A pointwise CLT alone proves none of the stronger statements. A supremum of coverage probabilities only asks whether the best covered member of the class has adequate coverage; it does not protect the worst case.

[Back to the solution map](#solution-map)
