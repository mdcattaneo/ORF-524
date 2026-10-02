# ORF 524 Practice Module 5: Solutions

**Status:** Released solutions  
**Last updated:** October 2, 2026  
**Practice module:** [Module 5](../topic-05-convergence-limit-theorems.md)  
**Primary chapter:** [Week 5](../../lectures/week-05.md)

Attempt each problem before consulting its solution. After studying a solution, close it and reconstruct the argument, including the assumptions used at each step.

<a id="solution-map"></a>

## Solution map

| Problem | Solution |
|---|---|
| [1. Convergence modes and stochastic orders](../topic-05-convergence-limit-theorems.md#problem-1) | [Solution](#problem-1) |
| [2. Laws of large numbers and Studentization](../topic-05-convergence-limit-theorems.md#problem-2) | [Solution](#problem-2) |
| [3. Triangular array and multivariate CLTs](../topic-05-convergence-limit-theorems.md#problem-3) | [Solution](#problem-3) |
| [4. A non Gaussian boundary limit](../topic-05-convergence-limit-theorems.md#problem-4) | [Solution](#problem-4) |
| [5. Continuous mapping and Slutsky under stress](../topic-05-convergence-limit-theorems.md#problem-5) | [Solution](#problem-5) |
| [6. Diagnose failed limit theorems](../topic-05-convergence-limit-theorems.md#problem-6) | [Solution](#problem-6) |

<a id="problem-1"></a>

## Problem 1. Convergence modes and stochastic orders

[Problem statement](../topic-05-convergence-limit-theorems.md#problem-1)

For every $\varepsilon>0$, Markov's inequality gives

$$
\mathbb{P}(|X_n-X|>\varepsilon)
\leq\frac{\mathbb E[|X_n-X|^p]}{\varepsilon^p}\to0.
$$

Thus $L^p$ convergence implies convergence in probability.

If $X_n\rightsquigarrow X$, choose continuity points $-M,M$ of the c.d.f. of $X$ with $\mathbb{P}(|X|>M)$ arbitrarily small. Weak convergence makes both tails of $X_n$ close to those of $X$ for all sufficiently large $n$, proving $X_n=O_{\mathbb{P}}(1)$. A larger $M$ handles the finitely many initial terms if the definition is required for every $n$.

For the first false converse, let $Z$ be Rademacher, set $X_n=(-1)^nZ$, and take $X=Z$. Every $X_n$ has the law of $X$, so $X_n\rightsquigarrow X$. Along odd $n$, $\mathbb{P}(|X_n-X|>1)=1$, so convergence in probability fails. For the second, let $U\sim\mathsf{Uniform}(0,1)$ and $X_n=n\mathbf 1\lbrace U\leq1/n\rbrace$. Then $X_n\to_{\mathbb{P}}0$ but $\mathbb E[X_n]=1$.

After division by the normalizations, the product claim reduces to $O_{\mathbb{P}}(1)o_{\mathbb{P}}(1)=o_{\mathbb{P}}(1)$. Bound the $O_{\mathbb{P}}(1)$ factor with high probability and then use convergence of the other factor to zero. For the sum,

$$
\frac{|A_n+C_n|}{a_n+b_n}
\leq\frac{a_n}{a_n+b_n}\frac{|A_n|}{a_n}
+\frac{b_n}{a_n+b_n}\frac{|C_n|}{b_n},
$$

which is bounded in probability.

[Back to the solution map](#solution-map)

<a id="problem-2"></a>

## Problem 2. Laws of large numbers and Studentization

[Problem statement](../topic-05-convergence-limit-theorems.md#problem-2)

Let $\overline\mu_n=n^{-1}\sum_i\mu_i$. Independence gives

$$
\mathbb V\mkern-3mu\left(\frac1n\sum_iX_i\right)
=\frac1{n^2}\sum_i\mathbb V[X_i]\to0.
$$

Chebyshev yields concentration around $\overline\mu_n$, and $\overline\mu_n\to\mu$ completes the proof.

In the dependent case,

$$
\mathbb V(\overline X_n)
=\frac{\sigma^2}{n}
+\frac2{n^2}\sum_{h=1}^{n-1}(n-h)\rho(h).
$$

Given $\eta>0$, choose $H$ so that $|\rho(h)|<\eta$ for $h\geq H$. The finitely many low lag terms are $O(1/n)$ after normalization, while the remaining absolute contribution is at most $2\eta$. Letting first $n\to\infty$ and then $\eta\downarrow0$ gives variance convergence to zero; Chebyshev again applies.

For the triangular array, independence within each row gives

$$
\mathbb V\mkern-3mu\left[
\frac1{v_n}\sum_i\lbrace X_{ni}-\mathbb E[X_{ni}]\rbrace
\right]
=\frac1{v_n^2}\sum_i\mathbb V[X_{ni}]
\to0.
$$

Chebyshev proves convergence in probability to zero. Independence eliminates off diagonal covariances in parts 1 and 3. In part 2 the covariance decay condition controls those terms instead. In every case the normalization is the scale relative to which the variance must vanish.

For the sample variance, expanding around the population mean gives

$$
\begin{aligned}
S_n^2
&=\frac1{n-1}\sum_{i=1}^n(X_i-\overline X_n)^2\\
&=\frac{n}{n-1}
\left\lbrace
\frac1n\sum_{i=1}^n(X_i-\mu)^2-(\overline X_n-\mu)^2
\right\rbrace.
\end{aligned}
$$

The weak LLN applied to $(X_i-\mu)^2$ gives $n^{-1}\sum_i(X_i-\mu)^2\to_{\mathbb{P}}\sigma^2$, while $\overline X_n\to_{\mathbb{P}}\mu$ and continuous mapping give $(\overline X_n-\mu)^2\to_{\mathbb{P}}0$. Since $n/(n-1)\to1$, Slutsky yields $S_n^2\to_{\mathbb{P}}\sigma^2$. Because $\sigma>0$, continuous mapping also gives $S_n\to_{\mathbb{P}}\sigma$ and $\sigma/S_n\to_{\mathbb{P}}1$. Therefore

$$
\frac{\sqrt n(\overline X_n-\mu)}{S_n}
=\frac{\sqrt n(\overline X_n-\mu)}{\sigma}
\frac{\sigma}{S_n}
\rightsquigarrow\mathcal N(0,1)
$$

by the classical CLT and Slutsky. Independence and identical distribution support the LLNs and CLT; the finite second moment makes the variance target finite and supplies the classical CLT; the chosen denominator and $\sqrt n$ normalization determine the respective probability and distributional limits.

[Back to the solution map](#solution-map)

<a id="problem-3"></a>

## Problem 3. Triangular array and multivariate CLTs

[Problem statement](../topic-05-convergence-limit-theorems.md#problem-3)

For the Bernoulli array,

$$
s_n^2=\sum_i p_{ni}(1-p_{ni})
\geq n\varepsilon(1-\varepsilon)\to\infty.
$$

Also $|X_{ni}|\leq1$. For every fixed $\delta>0$, eventually $\delta s_n>1$, so every indicator $\mathbf 1\lbrace|X_{ni}|>\delta s_n\rbrace$ is zero. The Lindeberg sum is then zero, and Lindeberg--Feller gives the asserted CLT. Boundedness works because the total standard deviation diverges; without that fact, the Lindeberg cutoff need not exceed the bound on individual summands.

For the vector result, fix $a\in\mathbb R^d$. The scalar variables $a'Y_i$ are i.i.d. with mean $a'\mu$ and variance $a'\Sigma a$, so

$$
a'\sqrt n(\overline Y_n-\mu)
\rightsquigarrow\mathcal N(0,a'\Sigma a).
$$

This remains a valid degenerate limit when $a'\Sigma a=0$. These are the projections of $\mathcal N_d(0,\Sigma)$, so Cramer--Wold gives the vector limit even when $\Sigma$ is singular.

An LLN uses a normalization under which the centered error vanishes in probability. A CLT uses a finer normalization under which the error has a nondegenerate distributional limit. Thus a CLT describes the first order shape and scale that an LLN alone does not supply.

[Back to the solution map](#solution-map)

<a id="problem-4"></a>

## Problem 4. A non Gaussian boundary limit

[Problem statement](../topic-05-convergence-limit-theorems.md#problem-4)

For $0\leq m\leq\theta$,

$$
\mathbb{P}_\theta(M_n\leq m)=\left(\frac m\theta\right)^n.
$$

For every $\delta\in(0,\theta)$,

$$
\mathbb{P}_\theta(|M_n-\theta|>\delta)
=\left(1-\frac\delta\theta\right)^n\to0,
$$

so $M_n\to_{\mathbb{P}}\theta$. For fixed $t\geq0$ and eventually $t<n$,

$$
\mathbb{P}_\theta\mkern-3mu\left\lbrace\frac{n(\theta-M_n)}{\theta}>t\right\rbrace
=\left(1-\frac tn\right)^n\to e^{-t}.
$$

This proves the exponential limit and implies $n(\theta-M_n)/\theta=O_{\mathbb{P}}(1)$, hence $\theta-M_n=O_{\mathbb{P}}(n^{-1})$. The error is one sided, has a faster rate than root $n$, and has a non-Gaussian limit.

Because $M_n\leq\theta$ almost surely,

$$
\begin{aligned}
\mathbb{P}_\theta\lbrace\theta\in C_n(X)\rbrace
&=\mathbb{P}_\theta\lbrace M_n\geq\theta\alpha^{1/n}\rbrace\\
&=1-\mathbb{P}_\theta\lbrace M_n<\theta\alpha^{1/n}\rbrace
=1-\alpha.
\end{aligned}
$$

The interval is exact for every $n$. The exponential result is instead a large sample approximation to the scaled endpoint error.

[Back to the solution map](#solution-map)

<a id="problem-5"></a>

## Problem 5. Continuous mapping and Slutsky under stress

[Problem statement](../topic-05-convergence-limit-theorems.md#problem-5)

Slutsky's theorem gives joint convergence $(X_n,Y_n)\rightsquigarrow(X,c)$. Continuous mapping then yields

$$
X_n+Y_n\rightsquigarrow X+c,\qquad
X_nY_n\rightsquigarrow cX,
$$

and, if $c\neq0$, $X_n/Y_n\rightsquigarrow X/c$.

Let $g(x)=\mathbf 1\lbrace x\geq0\rbrace$, so $g(0)=1$. Although $X_n=Z/n\to_{\mathbb{P}}0$, $\mathbb{P}\lbrace g(X_n)\neq g(0)\rbrace=1/2$ for every $n$. The limit hits the discontinuity with positive probability.

Finally, let $U_n=Z$ for every $n$, let $V_n=Z$ for even $n$, and let $V_n=-Z$ for odd $n$. Both marginals are standard Normal for every $n$, but $U_n+V_n$ alternates between $2Z$ and zero. It therefore has no limiting distribution. Joint information, not merely marginal convergence, is needed.

[Back to the solution map](#solution-map)

<a id="problem-6"></a>

## Problem 6. Diagnose failed limit theorems

[Problem statement](../topic-05-convergence-limit-theorems.md#problem-6)

The sample mean in the first construction is

$$
\overline X_n=Z+\overline\varepsilon_n\to_{\mathbb{P}}Z,
$$

not zero. The common shock prevents the dependence from averaging away.

Averages of i.i.d. standard Cauchy variables remain standard Cauchy, so they do not converge in probability to zero. The usual finite first moment LLN assumption fails.

In the triangular array, the standardized row sum is exactly the Rademacher variable $X_{n1}$ for every $n$, so its distribution is not asymptotically Normal. For any $\delta\in(0,1)$ the Lindeberg ratio equals one; a single term carries all row variance.

[Back to the solution map](#solution-map)
