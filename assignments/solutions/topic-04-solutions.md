# ORF 524 Practice Module 4: Solutions

**Status:** Released solutions  
**Last updated:** October 2, 2026  
**Practice module:** [Module 4](../topic-04-testing-confidence-sets.md)  
**Primary chapter:** [Week 4](../../lectures/week-04.md)

Attempt each problem before consulting its solution. After studying a solution, close it and reconstruct the argument, including the assumptions used at each step.

<a id="solution-map"></a>

## Solution map

| Problem | Solution |
|---|---|
| [1. Exact Normal inference with unknown variance](../topic-04-testing-confidence-sets.md#problem-1) | [Solution](#problem-1) |
| [2. Neyman--Pearson, discreteness, and MLR](../topic-04-testing-confidence-sets.md#problem-2) | [Solution](#problem-2) |
| [3. P-values as random variables](../topic-04-testing-confidence-sets.md#problem-3) | [Solution](#problem-3) |
| [4. Test inversion with parameter dependent support](../topic-04-testing-confidence-sets.md#problem-4) | [Solution](#problem-4) |
| [5. Fisher's combination under independence](../topic-04-testing-confidence-sets.md#problem-5) | [Solution](#problem-5) |
| [6. A prediction interval is not a confidence interval](../topic-04-testing-confidence-sets.md#problem-6) | [Solution](#problem-6) |
| [7. Two sample Student inference](../topic-04-testing-confidence-sets.md#problem-7) | [Solution](#problem-7) |

<a id="problem-1"></a>

## Problem 1. Exact Normal inference with unknown variance

[Problem statement](../topic-04-testing-confidence-sets.md#problem-1)

Under Gaussian sampling,

$$
\frac{\sqrt n(\overline X_n-\mu)}\sigma\sim\mathcal N(0,1),
\qquad
\frac{(n-1)S_n^2}{\sigma^2}\sim\chi^2_{n-1},
$$

and the two quantities are independent. Their ratio gives $T_n(\mu)\sim t_{n-1}$.

For the one sided problem, reject when

$$
T_n(\mu_0)>t_{n-1,1-\alpha}.
$$

Under $(\mu,\sigma)$ this statistic has a noncentral $t$ distribution with noncentrality $\lambda=\sqrt n(\mu-\mu_0)/\sigma$. If $G_{\nu,\lambda}$ is its c.d.f., power is

$$
\beta(\mu,\sigma)
=1-G_{n-1,\lambda}(t_{n-1,1-\alpha}).
$$

This is increasing in $\lambda$. Over $\mu\leq\mu_0$, $\lambda\leq0$, so rejection probability is at most its boundary value $\alpha$; equality occurs at $\mu=\mu_0$ for every $\sigma$.

For fixed $\mu>\mu_0$, $\lambda\to\infty$, so power tends to one. Under $\mu_n=\mu_0+h/\sqrt n$ with fixed $\sigma$, the noncentrality is instead $h/\sigma$. Since the critical value tends to $z_{1-\alpha}$ and the noncentral $t$ statistic converges to $\mathcal N(h/\sigma,1)$, the limiting power is

$$
1-\Phi\mkern-3mu\left(z_{1-\alpha}-\frac h\sigma\right).
$$

For the point null, reject when

$$
|T_n(\mu_0)|>t_{n-1,1-\alpha/2}.
$$

Let $c=t_{n-1,1-\alpha/2}$ and $\lambda=\sqrt n(\mu-\mu_0)/\sigma$. Since $T_n(\mu_0)$ has a noncentral $t_{n-1}$ distribution under $(\mu,\sigma)$, the two sided power is

$$
\beta_2(\mu,\sigma)
=G_{n-1,\lambda}(-c)+1-G_{n-1,\lambda}(c).
$$

At $\mu=\mu_0$, $\lambda=0$ and $\beta_2(\mu_0,\sigma)=\alpha$. The power is symmetric about $\mu_0$, increases with $|\mu-\mu_0|$ for fixed $\sigma$, and tends to one as $|\mu-\mu_0|\to\infty$. One way to see the shape is to condition on the chi squared denominator: for every fixed denominator, the probability that a shifted standard Normal variable falls outside a symmetric interval increases with the absolute shift.

Collecting nonrejected values gives

$$
\left[
\overline X_n-t_{n-1,1-\alpha/2}\frac{S_n}{\sqrt n},
\overline X_n+t_{n-1,1-\alpha/2}\frac{S_n}{\sqrt n}
\right].
$$

Opposite alternatives favor opposite tails, so the two sided test is generally not unrestricted UMP. In the Gaussian model with unknown variance, the usual equal tail Student test is UMP among unbiased level $\alpha$ tests under the standard Normal family theorem: it is UMPU, a restricted conclusion.

[Back to the solution map](#solution-map)

<a id="problem-2"></a>

## Problem 2. Neyman--Pearson, discreteness, and MLR

[Problem statement](../topic-04-testing-confidence-sets.md#problem-2)

For $p_1>p_0$,

$$
\frac{f_{p_1}(x)}{f_{p_0}(x)}
=\left(\frac{1-p_1}{1-p_0}\right)^n
\left\lbrace\frac{p_1(1-p_0)}{p_0(1-p_1)}\right\rbrace^{T(x)}.
$$

The term in braces exceeds one, so the ratio increases strictly in $T$. Define

$$
\phi(T)=\mathbf 1\lbrace T>c\rbrace+\gamma\mathbf 1\lbrace T=c\rbrace,
$$

where

$$
\gamma
=\frac{\alpha-\mathbb{P}_{p_0}(T>c)}{\mathbb{P}_{p_0}(T=c)}\in[0,1].
$$

Then $\mathbb E_{p_0}[\phi]=\alpha$. Neyman--Pearson makes it most powerful against the specified $p_1$. The Bernoulli family has MLR in $T$, so the same upper tail rule is UMP level $\alpha$ for $p\leq p_0$ versus every $p>p_0$; the null rejection probability is largest at $p_0$. The first claim solves one constrained optimization problem. The second establishes simultaneous optimality over a collection of alternatives and requires the MLR theorem.

[Back to the solution map](#solution-map)

<a id="problem-3"></a>

## Problem 3. P-values as random variables

[Problem statement](../topic-04-testing-confidence-sets.md#problem-3)

For continuous $F_0$, the probability integral transform gives $F_0(T)\sim\mathsf{Uniform}(0,1)$, and therefore so does $1-F_0(T)$.

For the Normal composite null, let

$$
Z_0=\frac{\sqrt n(\overline X_n-\mu_0)}\sigma
\sim\mathcal N(\delta,1),
\qquad
\delta=\frac{\sqrt n(\mu-\mu_0)}\sigma\leq0.
$$

For $u\in[0,1]$,

$$
\mathbb{P}_\mu(p\leq u)
=\mathbb{P}_\mu(Z_0\geq z_{1-u})
=1-\Phi(z_{1-u}-\delta)
\leq u.
$$

Equality holds at $\delta=0$, proving both validity and boundary uniformity.

For the Binomial tail p-value, the event $\lbrace p(T)\leq u\rbrace$ is either empty or an upper tail $\lbrace T\geq k_u\rbrace$. By the definition of the smallest such attainable tail,

$$
\mathbb{P}_{p_0}(p(T)\leq u)
=\mathbb{P}_{p_0}(T\geq k_u)\leq u.
$$

Because the p-value has only finitely many attainable values, it cannot be continuously Uniform and is generally conservative. In all cases the probability is computed from the sampling law under the null; no probability distribution over the truth of the hypothesis has been introduced.

[Back to the solution map](#solution-map)

<a id="problem-4"></a>

## Problem 4. Test inversion with parameter dependent support

[Problem statement](../topic-04-testing-confidence-sets.md#problem-4)

For $0\leq u\leq1$,

$$
\mathbb{P}_\theta(M_n/\theta\leq u)=u^n.
$$

Thus

$$
\mathbb{P}_{\theta_0}(M_n/\theta_0<q_L)=\alpha/2,
\qquad
\mathbb{P}_{\theta_0}(M_n/\theta_0>q_U)=\alpha/2.
$$

The two rejection events are disjoint, so the test has exact size $\alpha$. Nonrejection is

$$
q_L\leq M_n/\theta_0\leq q_U,
$$

which is equivalent to

$$
M_n/q_U\leq\theta_0\leq M_n/q_L.
$$

Because the parameter space is $(0,\infty)$, the inverted confidence set is therefore

$$
C(X)
=\left[\frac{M_n}{q_U},\frac{M_n}{q_L}\right]\cap(0,\infty).
$$

Directly,

$$
\mathbb{P}_\theta\mkern-3mu\left(\frac{M_n}{q_U}\leq\theta\leq\frac{M_n}{q_L}\right)
=\mathbb{P}_\theta(q_L\leq M_n/\theta\leq q_U)=1-\alpha.
$$

Since $q_U<1$, the lower endpoint exceeds $M_n$ whenever $M_n>0$. Confidence sets are created by a chosen family of tests; no theorem requires them to contain the MLE. The likelihood maximizer and equal tail test optimize different criteria. If $M_n=0$, the displayed interval before intersection is $[0,0]$, so the actual inverted set in $(0,\infty)$ is empty. Equivalently, $q_L\leq M_n/\theta_0$ fails for every $\theta_0>0$. The likelihood has no maximizer in the open parameter space on this sample, and the sample has probability zero under every $\theta>0$, so the edge case does not affect coverage.

[Back to the solution map](#solution-map)

<a id="problem-5"></a>

## Problem 5. Fisher's combination under independence

[Problem statement](../topic-04-testing-confidence-sets.md#problem-5)

If $U\sim\mathsf{Uniform}(0,1)$, then

$$
\mathbb{P}(-2\log U\leq y)=1-e^{-y/2},\qquad y\geq0,
$$

so $-2\log U\sim\chi^2_2$. Independence makes the sum chi squared with $2K$ degrees of freedom. Dependence destroys the convolution argument. Discrete or merely valid p-values are not exactly Uniform, so the exact chi squared calibration generally fails even if they are independent.

[Back to the solution map](#solution-map)

<a id="problem-6"></a>

## Problem 6. A prediction interval is not a confidence interval

[Problem statement](../topic-04-testing-confidence-sets.md#problem-6)

With the mean parameterization, $X_i/\theta\sim\mathsf{Exponential}(1)$, so

$$
\frac{2n\overline X_n}{\theta}=\frac{2\sum_iX_i}{\theta}\sim\chi^2_{2n}.
$$

Independence of $X_{n+1}$ and the first $n$ observations gives

$$
\frac{X_{n+1}}{\overline X_n}
=\frac{(2X_{n+1}/\theta)/2}{(2n\overline X_n/\theta)/(2n)}
\sim F_{2,2n}.
$$

If $f_q$ denotes the $q$ quantile of this $F$ distribution, the equal tail prediction interval is

$$
\left[\overline X_n f_{\alpha/2},
\overline X_n f_{1-\alpha/2}\right].
$$

Its repeated sampling probability concerns containment of the random future observation $X_{n+1}$. A confidence interval instead concerns containment of the fixed unknown parameter $\theta$ by a random set.

[Back to the solution map](#solution-map)

<a id="problem-7"></a>

## Problem 7. Two sample Student inference

[Problem statement](../topic-04-testing-confidence-sets.md#problem-7)

Independence of the two sample means gives

$$
D\sim\mathcal N\mkern-3mu\left(\Delta,\sigma^2\left(\frac1m+\frac1n\right)\right).
$$

Within each sample, the Gaussian mean is independent of the sample variance and

$$
\frac{(m-1)S_X^2}{\sigma^2}\sim\chi^2_{m-1},
\qquad
\frac{(n-1)S_Y^2}{\sigma^2}\sim\chi^2_{n-1}.
$$

The two quantities are independent, so their sum is $\chi^2_{N-2}$. The pair of sample variances is also independent of the pair of sample means, proving both

$$
\frac{(N-2)S_p^2}{\sigma^2}\sim\chi^2_{N-2}
$$

and its independence from $D$.

Under a general $\Delta$,

$$
\frac{D-\Delta_0}{\sigma\sqrt{1/m+1/n}}
\sim\mathcal N(\lambda,1),
\qquad
\lambda=\frac{\Delta-\Delta_0}{\sigma\sqrt{1/m+1/n}}.
$$

Dividing by the independent square root of a $\chi^2_{N-2}/(N-2)$ variable proves that $T(\Delta_0)$ has the noncentral $t_{N-2}$ law with noncentrality $\lambda$. Under the null, $\lambda=0$ and the law is central.

For $H_0:\Delta=0$, the unrestricted maximum likelihood estimates of the means are $\overline X_m$ and $\overline Y_n$, while the restricted mean estimate is

$$
\widetilde\mu=\frac{m\overline X_m+n\overline Y_n}{N}.
$$

Let

$$
W=(m-1)S_X^2+(n-1)S_Y^2=(N-2)S_p^2.
$$

The unrestricted and restricted residual sums of squares are $W$ and

$$
W+\frac{mn}{N}D^2,
$$

respectively. Hence the corresponding variance MLEs are $\widehat\sigma_U^2=W/N$ and $\widehat\sigma_R^2=(W+mnD^2/N)/N$. At a Normal maximum the likelihood is proportional to the negative $N/2$ power of the residual variance, so

$$
\begin{aligned}
\Lambda(X,Y)
&=\left(\frac{\widehat\sigma_U^2}{\widehat\sigma_R^2}\right)^{N/2}\\
&=\left\lbrace1+\frac{mnD^2/N}{(N-2)S_p^2}\right\rbrace^{-N/2}\\
&=\left\lbrace1+\frac{T(0)^2}{N-2}\right\rbrace^{-N/2}.
\end{aligned}
$$

The likelihood ratio decreases with $|T(0)|$, so the level $\alpha$ test rejects when

$$
|T(0)|>c,
\qquad
c=t_{N-2,1-\alpha/2}.
$$

If $G_{\nu,\lambda}$ is the noncentral $t$ c.d.f., its exact power at $(\Delta,\sigma)$ is

$$
G_{N-2,\lambda}(-c)+1-G_{N-2,\lambda}(c),
\qquad
\lambda=\frac{\Delta}{\sigma\sqrt{1/m+1/n}}.
$$

If $m,n\to\infty$ with $m/N\to\rho\in(0,1)$ and fixed $\Delta\neq0$, then $|\lambda|\to\infty$ while $c\to z_{1-\alpha/2}$, so power tends to one. The equal tail Student test is UMP among unbiased level $\alpha$ tests under the standard two sample Normal family theorem. It is not unrestricted UMP because alternatives on opposite sides favor opposite tails.

Inverting the tests gives

$$
\left[
D-cS_p\sqrt{\frac1m+\frac1n},
D+cS_p\sqrt{\frac1m+\frac1n}
\right].
$$

The exact Student pivot proves coverage $1-\alpha$ for every $\mu_X,\mu_Y\in\mathbb R$ and $\sigma>0$.

[Back to the solution map](#solution-map)
