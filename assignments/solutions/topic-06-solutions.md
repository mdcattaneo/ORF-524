# ORF 524 Practice Module 6: Solutions

**Status:** Released solutions  
**Last updated:** October 2, 2026  
**Practice module:** [Module 6](../topic-06-delta-asymptotic-inference.md)  
**Primary chapter:** [Week 6](../../lectures/week-06.md)

Attempt each problem before consulting its solution. After studying a solution, close it and reconstruct the argument, including the assumptions used at each step.

<a id="solution-map"></a>

## Solution map

| Problem | Solution |
|---|---|
| [1. First and second order delta methods](../topic-06-delta-asymptotic-inference.md#problem-1) | [Solution](#problem-1) |
| [2. A complete asymptotic chain for sample variance](../topic-06-delta-asymptotic-inference.md#problem-2) | [Solution](#problem-2) |
| [3. Pointwise level without uniform level](../topic-06-delta-asymptotic-inference.md#problem-3) | [Solution](#problem-3) |
| [4. A transformed, studentized confidence interval](../topic-06-delta-asymptotic-inference.md#problem-4) | [Solution](#problem-4) |
| [5. Variance stabilization for a Bernoulli proportion](../topic-06-delta-asymptotic-inference.md#problem-5) | [Solution](#problem-5) |
| [6. A general second order delta method](../topic-06-delta-asymptotic-inference.md#problem-6) | [Solution](#problem-6) |
| [7. Higher order bias and MSE calculations](../topic-06-delta-asymptotic-inference.md#problem-7) | [Solution](#problem-7) |

<a id="problem-1"></a>

## Problem 1. First and second order delta methods

[Problem statement](../topic-06-delta-asymptotic-inference.md#problem-1)

The Bernoulli CLT gives

$$
\sqrt n(\widehat p_n-p)\rightsquigarrow\mathcal N(0,p(1-p)).
$$

Since $g'(p)=1-2p$, for $p\neq1/2$,

$$
\sqrt n\lbrace g(\widehat p_n)-g(p)\rbrace
\rightsquigarrow
\mathcal N(0,(1-2p)^2p(1-p)).
$$

At $p=1/2$, the exact identity

$$
g(\widehat p_n)-g(1/2)=-(\widehat p_n-1/2)^2
$$

shows that the expression scaled by $\sqrt n$ converges in probability to zero, while

$$
n\lbrace g(\widehat p_n)-g(1/2)\rbrace
\rightsquigarrow-\frac14\chi_1^2.
$$

For fixed $p\neq1/2$, a consistent variance estimator is

$$
\widehat V_n=(1-2\widehat p_n)^2\widehat p_n(1-\widehat p_n).
$$

An asymptotic interval is

$$
g(\widehat p_n)
\pm z_{1-\alpha/2}\sqrt{\widehat V_n/n}.
$$

As $p\to1/2$, the first derivative and first order variance vanish while the correct rate changes to $n$ and the limit becomes nonnormal. The fixed $p$ approximation is not uniformly valid over neighborhoods containing $1/2$.

[Back to the solution map](#solution-map)

<a id="problem-2"></a>

## Problem 2. A complete asymptotic chain for sample variance

[Problem statement](../topic-06-delta-asymptotic-inference.md#problem-2)

Use

$$
\widehat\sigma_n^2
=\frac1n\sum_i(X_i-\mu)^2-(\overline X_n-\mu)^2.
$$

The LLN and continuous mapping give convergence to $\sigma^2$. Since $S_n^2=n\widehat\sigma_n^2/(n-1)$, it is consistent as well.

Subtracting $\sigma^2$ and multiplying by $\sqrt n$ gives

$$
\sqrt n(\widehat\sigma_n^2-\sigma^2)
=\frac1{\sqrt n}\sum_i\lbrace(X_i-\mu)^2-\sigma^2\rbrace
-\sqrt n(\overline X_n-\mu)^2.
$$

The last term is $o_{\mathbb{P}}(1)$. The first term obeys a CLT because its summand has finite variance under the fourth moment assumption. Hence

$$
\sqrt n(\widehat\sigma_n^2-\sigma^2)
\rightsquigarrow\mathcal N(0,V),
\qquad
V=\mathbb E[(X_1-\mu)^4]-\sigma^4.
$$

A consistent estimator is

$$
\widehat V_n
=\frac1n\sum_{i=1}^n
\lbrace(X_i-\overline X_n)^2-\widehat\sigma_n^2\rbrace^2.
$$

Laws of large numbers for the second and fourth centered sample moments give $\widehat V_n\to_{\mathbb{P}}V$. Moreover,

$$
\sqrt n(S_n^2-\widehat\sigma_n^2)
=\frac{\sqrt n}{n-1}\widehat\sigma_n^2=o_{\mathbb{P}}(1).
$$

If $V>0$, a feasible interval is

$$
\left[
\widehat\sigma_n^2-z_{1-\alpha/2}\sqrt{\widehat V_n/n},
\widehat\sigma_n^2+z_{1-\alpha/2}\sqrt{\widehat V_n/n}
\right].
$$

It has pointwise asymptotic coverage for each fixed distribution satisfying the conditions. It is neither exact nor uniformly valid over an unspecified class.

[Back to the solution map](#solution-map)

<a id="problem-3"></a>

## Problem 3. Pointwise level without uniform level

[Problem statement](../topic-06-delta-asymptotic-inference.md#problem-3)

The rejection probability is $e^{-n\theta}$. It converges to zero for every fixed $\theta>0$. However,

$$
\sup_{\theta\in(0,1]}e^{-n\theta}=1
$$

for every $n$, with the supremum approached as $\theta\downarrow0$. For example, $\theta_n=n^{-2}$ gives rejection probability $e^{-1/n}\to1$.

Pointwise asymptotic level fixes each $\theta$ before taking $n\to\infty$. Uniform asymptotic level takes a supremum over the null at every $n$ before the limit. This test has pointwise asymptotic level zero but not uniform asymptotic level zero.

[Back to the solution map](#solution-map)

<a id="problem-4"></a>

## Problem 4. A transformed, studentized confidence interval

[Problem statement](../topic-06-delta-asymptotic-inference.md#problem-4)

The CLT and the delta method with derivative $1/\mu$ give

$$
\sqrt n(\log\overline X_n-\log\mu)
\rightsquigarrow
\mathcal N\mkern-3mu\left(0,\frac{\sigma^2}{\mu^2}\right).
$$

Because $\overline X_n\to_{\mathbb{P}}\mu>0$, the arbitrary extension outside $\lbrace\overline X_n>0\rbrace$ is asymptotically irrelevant. Write $S_n=\sqrt{S_n^2}$. A consistent asymptotic variance estimator is

$$
\widehat V_n=\frac{S_n^2}{\overline X_n^2}.
$$

An interval for $\eta=\log\mu$ is

$$
\log\overline X_n
\pm z_{1-\alpha/2}\frac{S_n}{\sqrt n\mkern3mu\overline X_n}.
$$

Exponentiating the endpoints gives

$$
\left[
\overline X_n\exp\mkern-3mu\left\lbrace-z_{1-\alpha/2}\frac{S_n}{\sqrt n\mkern3mu\overline X_n}\right\rbrace,
\overline X_n\exp\mkern-3mu\left\lbrace z_{1-\alpha/2}\frac{S_n}{\sqrt n\mkern3mu\overline X_n}\right\rbrace
\right].
$$

For each fixed distribution satisfying the assumptions, coverage tends to $1-\alpha$. The conclusion is pointwise asymptotic, not exact or uniformly valid over an unspecified class. A symmetric interval on the log scale becomes multiplicatively, rather than additively, symmetric on the original scale.

[Back to the solution map](#solution-map)

<a id="problem-5"></a>

## Problem 5. Variance stabilization for a Bernoulli proportion

[Problem statement](../topic-06-delta-asymptotic-inference.md#problem-5)

The derivative is

$$
g'(p)=\frac1{\sqrt{p(1-p)}}.
$$

The delta method therefore gives limiting variance $g'(p)^2p(1-p)=1$. An interval on the transformed scale is

$$
\left[g(\widehat p_n)-\frac{z_{1-\alpha/2}}{\sqrt n},
g(\widehat p_n)+\frac{z_{1-\alpha/2}}{\sqrt n}\right].
$$

Since $g^{-1}(u)=\sin^2(u/2)$ is increasing on $[0,\pi]$, intersect the transformed interval with $[0,\pi]$ and apply this inverse to its endpoints. The derivative is unbounded at $0$ and $1$, and boundary sequences need different analysis, so the argument is pointwise on the interior.

[Back to the solution map](#solution-map)

<a id="problem-6"></a>

## Problem 6. A general second order delta method

[Problem statement](../topic-06-delta-asymptotic-inference.md#problem-6)

A second order Taylor expansion gives

$$
g(T_n)-g(\theta)
=\frac12g''(\theta)(T_n-\theta)^2+o_{\mathbb{P}}\lbrace(T_n-\theta)^2\rbrace.
$$

Because $\sqrt n(T_n-\theta)=O_{\mathbb{P}}(1)$, the right side is $O_{\mathbb{P}}(n^{-1})$. Multiplication by $\sqrt n$ therefore gives convergence in probability to zero. Multiplication by $n$ and Slutsky give

$$
n\lbrace g(T_n)-g(\theta)\rbrace
\rightsquigarrow
\frac12g''(\theta)V\chi_1^2.
$$

This limit is one sided according to the sign of $g''(\theta)$ and is generally neither centered nor Normal. In Problem 1 at $p=1/2$, $g''(p)=-2$ and $V=p(1-p)=1/4$, producing $-\chi_1^2/4$.

[Back to the solution map](#solution-map)

<a id="problem-7"></a>

## Problem 7. Higher order bias and MSE calculations

[Problem statement](../topic-06-delta-asymptotic-inference.md#problem-7)

Let $Y_i=X_i-\mu$, so $\delta_n=n^{-1}\sum_iY_i$ and $\mathbb E[Y_i]=0$. In an expansion of a product of centered independent variables, a term has zero expectation whenever an index appears exactly once. Therefore

$$
\mathbb E[\delta_n^2]=\frac{nm_2}{n^2}=\frac{m_2}{n},
\qquad
\mathbb E[\delta_n^3]=\frac{nm_3}{n^3}=\frac{m_3}{n^2}.
$$

For the fourth moment, the nonzero terms either use one index four times or use two distinct indices twice each. There are $n$ terms of the first type and $3n(n-1)$ ordered pairings of the second type, so

$$
\begin{aligned}
\mathbb E[\delta_n^4]
&=\frac{nm_4+3n(n-1)m_2^2}{n^4}\\
&=\frac{3m_2^2}{n^2}+\frac{m_4-3m_2^2}{n^3}.
\end{aligned}
$$

The moment assumption and a standard moment inequality for centered sums give $\mathbb E[|\delta_n|^5]=O(n^{-5/2})$ and $\mathbb E[|\delta_n|^8]=O(n^{-4})$. Taylor expansion through order four has a remainder bounded by a constant times $|\delta_n|^5$. Taking expectations and substituting the moment formulas yields

$$
\begin{aligned}
\mathbb E[g(\overline X_n)]
={}&g(\mu)+\frac{g^{(2)}(\mu)m_2}{2n}\\
&+\frac1{n^2}\left\lbrace
\frac{g^{(3)}(\mu)m_3}{6}
+\frac{g^{(4)}(\mu)m_2^2}{8}
\right\rbrace
+o(n^{-2}).
\end{aligned}
$$

In particular,

$$
\mathrm{Bias}\lbrace g(\overline X_n)\rbrace^2
=\frac{g^{(2)}(\mu)^2m_2^2}{4n^2}+o(n^{-2}).
$$

For the MSE, write $a=g^{(1)}(\mu)$, $b=g^{(2)}(\mu)/2$, and $c=g^{(3)}(\mu)/6$. Taylor expansion through order three gives

$$
g(\mu+\delta_n)-g(\mu)
=a\delta_n+b\delta_n^2+c\delta_n^3+O(\delta_n^4).
$$

Squaring gives

$$
a^2\delta_n^2+2ab\delta_n^3+(b^2+2ac)\delta_n^4+R_n,
$$

where $\mathbb E[|R_n|]=o(n^{-2})$ by the fifth and eighth absolute moment bounds above. Substitution produces

$$
\begin{aligned}
\mathbb E[\lbrace g(\overline X_n)-g(\mu)\rbrace^2]
={}&\frac{g^{(1)}(\mu)^2m_2}{n}\\
&+\frac1{n^2}\left\lbrace
g^{(1)}(\mu)g^{(2)}(\mu)m_3
+g^{(1)}(\mu)g^{(3)}(\mu)m_2^2
+\frac34g^{(2)}(\mu)^2m_2^2
\right\rbrace
+o(n^{-2}).
\end{aligned}
$$

When $g^{(1)}(\mu)\neq0$, the leading MSE is $g^{(1)}(\mu)^2m_2/n$. When $g^{(1)}(\mu)=0$ and $g^{(2)}(\mu)\neq0$, the leading MSE is $3g^{(2)}(\mu)^2m_2^2/(4n^2)$. The second order delta method gives

$$
n\lbrace g(\overline X_n)-g(\mu)\rbrace
\rightsquigarrow
\frac12g^{(2)}(\mu)m_2\chi_1^2,
$$

whose squared second moment is $3g^{(2)}(\mu)^2m_2^2/4$, matching the expansion.

For $g(x)=x^2$, the bias is exactly $m_2/n$, and

$$
\mathbb E[(\overline X_n^2-\mu^2)^2]
=\frac{4\mu^2m_2}{n}
+\frac{4\mu m_3+3m_2^2}{n^2}
+\frac{m_4-3m_2^2}{n^3}.
$$

At $\mu=0$, the rate changes: $n\overline X_n^2\rightsquigarrow m_2\chi_1^2$ and $n^2\mathbb E[\overline X_n^4]\to3m_2^2$.

[Back to the solution map](#solution-map)
