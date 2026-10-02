# ORF 524 Practice Module 2: Solutions

**Status:** Released solutions  
**Last updated:** October 2, 2026  
**Practice module:** [Module 2](../topic-02-models-likelihood-sufficiency.md)  
**Primary chapter:** [Week 2](../../lectures/week-02.md)

Attempt each problem before consulting its solution. After studying a solution, close it and reconstruct the argument, including the assumptions used at each step.

<a id="solution-map"></a>

## Solution map

| Problem | Solution |
|---|---|
| [1. Identification in a threshold model](../topic-02-models-likelihood-sufficiency.md#problem-1) | [Solution](#problem-1) |
| [2. Normalizations in fixed effects models](../topic-02-models-likelihood-sufficiency.md#problem-2) | [Solution](#problem-2) |
| [3. Uniform likelihood with moving support](../topic-02-models-likelihood-sufficiency.md#problem-3) | [Solution](#problem-3) |
| [4. How much does sufficiency reduce?](../topic-02-models-likelihood-sufficiency.md#problem-4) | [Solution](#problem-4) |
| [5. A censored observation needs a mixed measure](../topic-02-models-likelihood-sufficiency.md#problem-5) | [Solution](#problem-5) |
| [6. Must an MLE be a function of a sufficient statistic?](../topic-02-models-likelihood-sufficiency.md#problem-6) | [Solution](#problem-6) |
| [7. Exact Gaussian sampling calculations](../topic-02-models-likelihood-sufficiency.md#problem-7) | [Solution](#problem-7) |

<a id="problem-1"></a>

## Problem 1. Identification in a threshold model

[Problem statement](../topic-02-models-likelihood-sufficiency.md#problem-1)

Since $(Z_i-\mu)/\sigma\sim\mathcal N(0,1)$,

$$
p(\theta)=\mathbb{P}_\theta(Y_i=1)
=\Phi\mkern-3mu\left(\frac{\kappa-\mu}{\sigma}\right)=\Phi(\delta(\theta)).
$$

Thus $Y_i\sim\mathsf{Bernoulli}(\Phi(\delta))$. Because $\Phi$ is strictly increasing, two parameters are observationally equivalent exactly when

$$
\frac{\kappa-\mu}{\sigma}
=\frac{\widetilde\kappa-\widetilde\mu}{\widetilde\sigma}.
$$

The full parameter is not identified, while $\delta$ is identified. If $(\mu,\sigma)$ are known, $\kappa=\mu+\sigma\delta$ is identified. If only $\kappa$ is known, the data identify the ratio $(\kappa-\mu)/\sigma$, not $\mu$ and $\sigma$ separately. If only $\sigma$ is known, the difference $\kappa-\mu=\sigma\delta$ is identified, not its two components. Under $\mu=0$ and $\sigma=1$, $\kappa=\delta$ is identified.

For true success probability $p_0$,

$$
M(p)=p_0\log p+(1-p_0)\log(1-p).
$$

Its derivative is $(p_0-p)/[p(1-p)]$ and its second derivative at the stationary point is negative; equivalently, $M(p_0)-M(p)$ is the Bernoulli KL divergence. Hence $p_0$ is the unique maximizer. Strict monotonicity of $\Phi$ identifies $\delta_0=\Phi^{-1}(p_0)$, but all parameter triples in its equivalence class remain population maximizers unless a normalization is imposed.

[Back to the solution map](#solution-map)

<a id="problem-2"></a>

## Problem 2. Normalizations in fixed effects models

[Problem statement](../topic-02-models-likelihood-sufficiency.md#problem-2)

For arbitrary $a,b\in\mathbb R$, define

$$
\widetilde\mu=\mu+a+b,
\qquad
\widetilde\alpha_i=\alpha_i-a,
\qquad
\widetilde\beta_t=\beta_t-b.
$$

Then $\widetilde\mu+\widetilde\alpha_i+\widetilde\beta_t=\mu+\alpha_i+\beta_t$ for every cell, proving a two dimensional redundancy.

Let

$$
\overline m_{\cdot\cdot}=\frac1{NT}\sum_{i,t}m_{it},
\quad
\overline m_{i\cdot}=\frac1T\sum_tm_{it},
\quad
\overline m_{\cdot t}=\frac1N\sum_im_{it}.
$$

Under the two zero sum restrictions,

$$
\mu=\overline m_{\cdot\cdot},
\qquad
\alpha_i=\overline m_{i\cdot}-\overline m_{\cdot\cdot},
\qquad
\beta_t=\overline m_{\cdot t}-\overline m_{\cdot\cdot}.
$$

These formulas prove identification. Equality of the observable Normal laws also equates their variances, so $\sigma^2$ is identified regardless of location normalization. Contrasts such as $\alpha_i-\alpha_j$ and $\beta_t-\beta_s$ are invariant and identified. The unrestricted intercept $\mu$ is not.

[Back to the solution map](#solution-map)

<a id="problem-3"></a>

## Problem 3. Uniform likelihood with moving support

[Problem statement](../topic-02-models-likelihood-sufficiency.md#problem-3)

The likelihood is

$$
L_n(\theta;x)=\theta^{-n}
\mathbf 1\lbrace x_{(1)}\geq0,\mkern3mu M_n\leq\theta\rbrace.
$$

On $M_n>0$, feasible values satisfy $\theta\geq M_n$, so $\theta^{-n}$ decreases and $\widehat\theta_{\mathrm{ML}}=M_n$. The optimum is a support boundary, not an interior score root. If every observation equals zero, the likelihood is $\theta^{-n}$ for every $\theta>0$ and has no maximizer: its supremum is approached as $\theta\downarrow0$. This is a null event under every $\theta>0$, so defining the estimator arbitrarily there does not change its sampling properties.

The full c.d.f. is

$$
\mathbb{P}_\theta(M_n\leq m) =
\begin{cases}
0, & m<0,\\
(m/\theta)^n, & 0\leq m<\theta,\\
1, & m\geq\theta.
\end{cases}
$$

Therefore, up to irrelevant endpoint values,

$$
f_{M_n}(m;\theta)=\frac{nm^{n-1}}{\theta^n}\mathbf 1\lbrace0<m<\theta\rbrace.
$$

Factorization through $M_n$ proves sufficiency. Consider two sample points $x,y$ in the union of the model supports, $[0,\infty)^n$. If $M_n(x)=M_n(y)$, their likelihood ratio is one for every parameter value at which it is defined. If, say, $M_n(x)<M_n(y)$, then for $M_n(x)\leq\theta<M_n(y)$ only $x$ is feasible, whereas for $\theta\geq M_n(y)$ both samples are feasible. Thus the ratio depends on $\theta$. The likelihood ratio criterion proves that $M_n$ is minimal sufficient.

For one observation under the true $\theta_0$, a candidate $\theta<\theta_0$ assigns density zero to $(\theta,\theta_0)$, an event of positive true probability, so $Q(\theta)=-\infty$. For $\theta\geq\theta_0$,

$$
Q(\theta)=-\log\theta.
$$

This is uniquely maximized over that region at $\theta_0$. The moving support rules out candidates below truth, monotonicity selects the boundary, and injectivity of the Uniform parameterization turns the selected distribution into the unique parameter.

[Back to the solution map](#solution-map)

<a id="problem-4"></a>

## Problem 4. How much does sufficiency reduce?

[Problem statement](../topic-02-models-likelihood-sufficiency.md#problem-4)

1. This is the Exponential distribution parameterized by its mean. Its likelihood is

$$
\theta^{-n}\exp\mkern-3mu\left(-\frac{\sum_i x_i}{\theta}\right)
\prod_i\mathbf 1\lbrace x_i>0\rbrace.
$$

With natural parameter $\eta=-1/\theta\in(-\infty,0)$, it has the form

$$
f(x;\eta)=\mathbf 1\lbrace x>0\rbrace\exp\lbrace\eta x-A(\eta)\rbrace,
\qquad
A(\eta)=-\log(-\eta).
$$

The support is fixed and the natural parameter space is open, so this is a regular one parameter exponential family. For two samples in the common support, their likelihood ratio is independent of $\theta$ exactly when their sample sums agree. Hence $\sum_iX_i$ is minimal sufficient.

2. This is a shifted Exponential distribution. Its likelihood is

$$
\exp\mkern-3mu\left(-\sum_i x_i+n\theta\right)
\mathbf 1\lbrace\theta<X_{(1)}\rbrace.
$$

Its support depends on $\theta$, so it is not regular in the stated sense. Samples with the same minimum have a likelihood ratio independent of $\theta$; samples with different minima become feasible at different parameter values. The ratio criterion therefore shows that $X_{(1)}$ is minimal sufficient.

3. This is the Laplace location family. It is not a fixed dimensional regular exponential family in $\theta$. For samples $x,y$, the log likelihood ratio is constant in $\theta$ exactly when

$$
D(\theta)=\sum_i|x_i-\theta|-\sum_i|y_i-\theta|
$$

is constant. Away from sample points,

$$
D'(\theta)=2\left|\lbrace i:x_i<\theta\rbrace\right|-2\left|\lbrace i:y_i<\theta\rbrace\right|.
$$

This derivative vanishes almost everywhere exactly when the two empirical counting functions, and hence the two multisets, agree. Therefore the ordered samples agree. Conversely, equal order statistics plainly give ratio one. The vector $(X_{(1)},\ldots,X_{(n)})$ is minimal sufficient, and its dimension grows with $n$.

[Back to the solution map](#solution-map)

<a id="problem-5"></a>

## Problem 5. A censored observation needs a mixed measure

[Problem statement](../topic-02-models-likelihood-sufficiency.md#problem-5)

Let $z_c=(c-\mu)/\sigma$. Then

$$
\mathbb{P}(Y=c)=\Phi(z_c),
$$

while for $y>c$ the density is $\sigma^{-1}\phi((y-\mu)/\sigma)$. With respect to $\nu=\delta_c+\lambda|_{(c,\infty)}$, one density is

$$
p_{\mu,\sigma}(y)
=\Phi(z_c)\mathbf 1\lbrace y=c\rbrace
+\frac1\sigma\phi\mkern-3mu\left(\frac{y-\mu}{\sigma}\right)\mathbf 1\lbrace y>c\rbrace.
$$

Thus the sample likelihood is $\prod_i p_{\mu,\sigma}(y_i)$. A purely Lebesgue density integrates only the continuous component and discards the likelihood contribution of every censored observation.

[Back to the solution map](#solution-map)

<a id="problem-6"></a>

## Problem 6. Must an MLE be a function of a sufficient statistic?

[Problem statement](../topic-02-models-likelihood-sufficiency.md#problem-6)

Because $h(x)>0$ does not depend on $\theta$,

$$
\arg\max_\theta f(x;\theta)=\arg\max_\theta g(T(x),\theta).
$$

The argmax correspondence therefore depends only on $T(x)$. If it is a singleton, its unique value is a function of $T$. With ties, a measurable rule can select one maximizer for each value of $T$, when the needed measurable selection exists. But a rule allowed to inspect the full sample could choose different elements of the same argmax set at two samples sharing $T$; such a selection is an MLE but not a function of $T$. Existence and measurability must be established rather than inferred from factorization alone.

[Back to the solution map](#solution-map)

<a id="problem-7"></a>

## Problem 7. Exact Gaussian sampling calculations

[Problem statement](../topic-02-models-likelihood-sufficiency.md#problem-7)

Linearity and independence give

$$
\mathbb E[\overline X_n]=\mu,
\qquad
\mathbb V(\overline X_n)=\frac{\sigma^2}{n}.
$$

Because a linear combination of jointly Gaussian variables is Gaussian,

$$
\overline X_n\sim\mathcal N\mkern-3mu\left(\mu,\frac{\sigma^2}{n}\right).
$$

Expanding around $\mu$ and using $\sum_i(X_i-\mu)=n(\overline X_n-\mu)$ gives

$$
\begin{aligned}
\sum_i(X_i-\overline X_n)^2
&=\sum_i\lbrace(X_i-\mu)-(\overline X_n-\mu)\rbrace^2\\
&=\sum_i(X_i-\mu)^2-n(\overline X_n-\mu)^2.
\end{aligned}
$$

Taking expectations yields

$$
(n-1)\mathbb E[S_n^2]
=n\sigma^2-n\frac{\sigma^2}{n}
=(n-1)\sigma^2,
$$

so $S_n^2$ is unbiased.

Now $Z=(X-\mu\mathbf 1_n)/\sigma\sim\mathcal N(0,I_n)$. Orthogonality gives

$$
U=QZ\sim\mathcal N(0,I_n),
$$

so $U_1,\ldots,U_n$ are independent standard Normal variables. The specified first row gives

$$
U_1=\frac{\sqrt n(\overline X_n-\mu)}{\sigma}.
$$

The remaining rows form an orthonormal basis for the space orthogonal to $\mathbf 1_n$, and therefore

$$
\sum_{j=2}^nU_j^2
=Z'\left(I_n-\frac1n\mathbf 1_n\mathbf 1_n'\right)Z
=\frac{(n-1)S_n^2}{\sigma^2}.
$$

Consequently the last expression is $\chi^2_{n-1}$ and is independent of $U_1$, hence independent of $\overline X_n$. Because a chi squared variable with $\nu$ degrees of freedom has mean $\nu$ and variance $2\nu$,

$$
\mathbb E[S_n^2]=\sigma^2,
\qquad
\mathbb V(S_n^2)=\frac{2\sigma^4}{n-1}.
$$

Finally,

$$
\frac{\sqrt n(\overline X_n-\mu)}{S_n}
=\frac{U_1}{\sqrt{\left(\sum_{j=2}^nU_j^2\right)/(n-1)}}
\sim t_{n-1}.
$$

Writing $c=t_{n-1,1-\alpha/2}$, inversion gives the exact interval

$$
\left[
\overline X_n-c\frac{S_n}{\sqrt n},
\overline X_n+c\frac{S_n}{\sqrt n}
\right],
$$

whose coverage is $1-\alpha$ for every $\mu\in\mathbb R$ and $\sigma>0$.

[Back to the solution map](#solution-map)
