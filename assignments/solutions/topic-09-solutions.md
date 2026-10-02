# ORF 524 Practice Module 9: Solutions

**Status:** Released solutions  
**Last updated:** October 2, 2026  
**Practice module:** [Module 9](../topic-09-m-z-estimation.md)  
**Primary chapter:** [Week 9](../../lectures/week-09.md)

Attempt each problem before consulting its solution. After studying a solution, close it and reconstruct the argument, including the assumptions used at each step.

<a id="solution-map"></a>

## Solution map

| Problem | Solution |
|---|---|
| [1. The representation determines the theorem](../topic-09-m-z-estimation.md#problem-1) | [Solution](#problem-1) |
| [2. Uniform convergence and argmin consistency](../topic-09-m-z-estimation.md#problem-2) | [Solution](#problem-2) |
| [3. Nonlinear least squares and sandwich variance](../topic-09-m-z-estimation.md#problem-3) | [Solution](#problem-3) |
| [4. Wald, score, and likelihood ratio in a Poisson model](../topic-09-m-z-estimation.md#problem-4) | [Solution](#problem-4) |
| [5. One Newton step can recover first order efficiency](../topic-09-m-z-estimation.md#problem-5) | [Solution](#problem-5) |
| [6. Misspecified likelihood needs the sandwich](../topic-09-m-z-estimation.md#problem-6) | [Solution](#problem-6) |
| [7. A curved Gaussian likelihood laboratory](../topic-09-m-z-estimation.md#problem-7) | [Solution](#problem-7) |

<a id="problem-1"></a>

## Problem 1. The representation determines the theorem

[Problem statement](../topic-09-m-z-estimation.md#problem-1)

For the mean, take

$$
\rho(w,\theta)=(w-\theta)^2,
\qquad
\psi(w,\theta)=w-\theta.
$$

Both population problems identify $\theta_0=\mathbb E[W]$ when the required moments exist.

For a $\tau$ quantile, minimize $\mathbb E_n\rho_\tau(W-\theta)$. The exact sample optimality condition is the subgradient condition

$$
F_n(\widehat\theta_n^-)
\leq\tau\leq
F_n(\widehat\theta_n).
$$

A corresponding nonsmooth estimating function is

$$
\psi(w,\theta)=\mathbf 1\lbrace w\leq\theta\rbrace-\tau.
$$

Its population zero is $F_W(\theta)=\tau$ when the c.d.f. is continuous and the quantile is unique. With continuous sampling and the convention $\widehat\theta_n=\inf\lbrace t:F_n(t)\geq\tau\rbrace$, there are no ties almost surely and $0\leq F_n(\widehat\theta_n)-\tau\leq1/n$. Thus this choice is an approximate Z-root with residual $o(n^{-1/2})$, although the ordinary smooth Z theorem still does not apply. Atoms can produce larger jumps and set valued minimizers, for which the subgradient formulation is the correct one.

For least squares,

$$
\rho((y,x),\beta)=(y-x'\beta)^2,
\qquad
\psi((y,x),\beta)=x(y-x'\beta).
$$

The population target solves $\mathbb E[X(Y-X'\beta_0)]=0$ and equals $\mathbb E[XX']^{-1}\mathbb E[XY]$ when $\mathbb E[XX']$ is nonsingular. An interior differentiable minimum satisfies the Z equation.

For Poisson, the negative log likelihood contribution can be written as

$$
\rho(w,\theta)=\theta-w\log\theta
$$

up to data only constants for $\theta>0$. Extend it at $\theta=0$ by setting $\rho(0,0)=0$ and $\rho(w,0)=+\infty$ for $w>0$. Then the MLE exists throughout $[0,\infty)$ and equals $\overline W_n$, including the all zero sample. The equation $\psi(w,\theta)=w-\theta$ has the same root, while the likelihood score $w/\theta-1$ is only an interior representation. Under correct specification the population target is the Poisson mean. Regular likelihood theory requires an interior true value $\theta_0>0$.

On the probability one event that the sample maximum is positive, the Uniform MLE minimizes negative log likelihood with its support indicator and equals that maximum. On the null sample the open parameter likelihood has no maximizer, so the estimator is defined arbitrarily. The optimum otherwise lies at a moving support boundary, so differentiating the smooth part does not give a valid score characterization. Consistency is naturally handled by an extremum/argmax argument or directly through the distribution of the maximum, not a smooth Z theorem.

[Back to the solution map](#solution-map)

<a id="problem-2"></a>

## Problem 2. Uniform convergence and argmin consistency

[Problem statement](../topic-09-m-z-estimation.md#problem-2)

Fix $\varepsilon>0$ and let $\eta_\varepsilon>0$ be the separation gap. On the event

$$
\sup_\theta|Q_n(\theta)-Q(\theta)|<\eta_\varepsilon/3,
\qquad
r_n<\eta_\varepsilon/3,
$$

every $\theta$ outside the $\varepsilon$ ball satisfies

$$
Q_n(\theta)>Q(\theta_0)+2\eta_\varepsilon/3,
$$

while

$$
Q_n(\widehat\theta_n)
\leq Q_n(\theta_0)+r_n
<Q(\theta_0)+2\eta_\varepsilon/3.
$$

Thus the approximate minimizer lies inside the ball on an event whose probability tends to one.

For the moving spike sequence, every fixed $\theta$ is eventually outside the shrinking interval centered near one, so $Q_n(\theta)\to\theta^2$. Because the interval is closed, for all sufficiently large $n$ its left endpoint is an attained minimizer; the criterion there is approximately $-1$, below $Q_n(0)=0$. This minimizer converges to one rather than zero. Also $\sup_\theta|Q_n-Q|=2$, displaying the failed condition.

A usable uniform LLN is: for i.i.d. observations, compact $\Theta$, almost sure continuity of $\theta\mapsto\rho(W,\theta)$, and a measurable envelope $F$ satisfying $\sup_{\theta\in\Theta}|\rho(W,\theta)|\leq F(W)$ with $\mathbb E[F(W)]<\infty$, standard measurability conditions imply

$$
\sup_{\theta\in\Theta}|\mathbb{E}_n\rho(W,\theta)-\mathbb{P}_0\rho(W,\theta)|
\longrightarrow0
$$

almost surely. For squared loss on $\mathcal B$, continuity in $\beta$ is immediate and

$$
\sup_{\beta\in\mathcal B}(Y-X'\beta)^2
\leq2Y^2+2K^2\Vert X\Vert^2.
$$

The assumed second moments therefore provide an integrable envelope, so the theorem applies.

Pointwise convergence selects a possibly different high probability index for each fixed $\theta$; it gives no control at the random, moving $\widehat\theta_n$. Uniform convergence controls all candidates simultaneously. The optimization error is a single random scalar independent of the candidate index and is explicitly required to vanish.

For Z-estimation, sufficient parallel conditions are

$$
\sup_\theta\Vert\Psi_n(\theta)-\Psi(\theta)\Vert\to_{\mathbb{P}}0,
\quad
\inf_{\Vert\theta-\theta_0\Vert\geq\varepsilon}\Vert\Psi(\theta)\Vert>0,
\quad
\Vert\Psi_n(\widehat\theta_n)\Vert=o_{\mathbb{P}}(1).
$$

They imply $\widehat\theta_n\to_{\mathbb{P}}\theta_0$ by the same outside the ball argument.

[Back to the solution map](#solution-map)

<a id="problem-3"></a>

## Problem 3. Nonlinear least squares and sandwich variance

[Problem statement](../topic-09-m-z-estimation.md#problem-3)

Write $m_0(X)=\mu(X'\beta_0)$. Conditional mean zero of $Y-m_0(X)$ makes the cross term vanish:

$$
\mathbb E[\lbrace Y-\mu(X'\beta)\rbrace^2]
=\mathbb E[\lbrace Y-m_0(X)\rbrace^2]
+\mathbb E[\lbrace m_0(X)-\mu(X'\beta)\rbrace^2].
$$

Thus $\beta_0$ is uniquely identified if the second term is positive for every $\beta\neq\beta_0$.

Differentiating the sample criterion gives $-2\mathbb{E}_n\psi(W,\beta)$, so an interior solution satisfies $\mathbb{E}_n\psi(W,\widehat\beta_n)=0$. At $\beta_0$,

$$
\Psi_\beta(\beta_0)
=\mathbb E[\partial\psi(W,\beta_0)/\partial\beta']
=-H,
\qquad
H=\mathbb E[\dot\mu(X'\beta_0)^2XX'],
$$

because the term containing $\ddot\mu(X'\beta_0)\lbrace Y-m_0(X)\rbrace$ has conditional mean zero. Also

$$
B=\mathbb V[\psi(W,\beta_0)]
=\mathbb E[\sigma^2(X)\dot\mu(X'\beta_0)^2XX'].
$$

If $H$ is nonsingular and the consistency, local derivative ULLN, and CLT conditions hold, the Z theorem gives asymptotic variance $H^{-1}BH^{-1}$.

With residuals $\widehat\varepsilon_i=Y_i-\mu(X_i'\widehat\beta_n)$, define

$$
\widehat H_n=\mathbb{E}_n[\dot\mu(X'\widehat\beta_n)^2XX'],
$$

$$
\widehat B_n=\mathbb{E}_n[\widehat\varepsilon^2
\dot\mu(X'\widehat\beta_n)^2XX'],
\qquad
\widehat V_n=\widehat H_n^{-1}\widehat B_n\widehat H_n^{-1}.
$$

Appropriate plug in LLNs give consistency. Under conditional homoskedasticity, $B=\sigma^2H$ and $V=\sigma^2H^{-1}$.

[Back to the solution map](#solution-map)

<a id="problem-4"></a>

## Problem 4. Wald, score, and likelihood ratio in a Poisson model

[Problem statement](../topic-09-m-z-estimation.md#problem-4)

With the boundary convention $0\log0=0$, the extended log likelihood is

$$
\ell_n(\theta)=n\overline X_n\log\theta-n\theta+\text{constant},
$$

so $\widehat\theta_n=\overline X_n$. Exactly,

$$
\sqrt n(\widehat\theta_n-\theta)
=\frac1{\sqrt n}\sum_i(X_i-\theta)
\rightsquigarrow\mathcal N(0,\theta),
$$

and $\widehat V_n=\overline X_n$ is consistent.

In particular, the MLE equals zero on the all zero sample rather than failing to exist. The null value $\theta_0>0$ is interior, so this boundary convention does not change the local regular asymptotic calculation.

The event $\overline X_n=0$ has probability $e^{-n\theta_0}\to0$ under the null, so any fixed definition of $W_n$ there is asymptotically irrelevant. On $\overline X_n>0$, the Wald statistic divides the squared estimator displacement by estimated asymptotic variance. The sample score at the null is $n(\overline X_n-\theta_0)/\theta_0$ and sample information is $n/\theta_0$, giving $LM_n$. Comparing maximized log likelihoods gives the stated $LR_n$.

Under the null, $\overline X_n-\theta_0=O_{\mathbb{P}}(n^{-1/2})$. Then

$$
\frac1{\overline X_n}=\frac1{\theta_0}+o_{\mathbb{P}}(1),
$$

so

$$
W_n
=\frac{n(\overline X_n-\theta_0)^2}{\theta_0}
+o_{\mathbb{P}}(1).
$$

Also

$$
(\theta_0+h)\log(1+h/\theta_0)-h
=\frac{h^2}{2\theta_0}+O(h^3),
$$

which, evaluated at $h=\overline X_n-\theta_0$, gives the same expansion for $LR_n$ because $n(\overline X_n-\theta_0)^3=o_{\mathbb{P}}(1)$. The LM expression is already the common leading term. The standardized CLT then gives a $\chi_1^2$ limit.

Wald evaluates distance using the unrestricted estimate and its variance; LM evaluates slope and information at the restricted value $\theta_0$; LR uses both restricted and unrestricted optimized criterion values. Their denominators and nonlinear forms differ in finite samples.

[Back to the solution map](#solution-map)

<a id="problem-5"></a>

## Problem 5. One Newton step can recover first order efficiency

[Problem statement](../topic-09-m-z-estimation.md#problem-5)

The mean value theorem gives, for an intermediate $\theta_n^\star$,

$$
\Psi_n(\overline\theta_n)
=\Psi_n(\theta_0)
+\Psi_{\theta,n}(\theta_n^\star)(\overline\theta_n-\theta_0).
$$

Substitution into the update yields

$$
\begin{aligned}
\sqrt n(\widetilde\theta_n-\theta_0)
={}&\left\lbrace1-\Psi_{\theta,n}(\overline\theta_n)^{-1}
\Psi_{\theta,n}(\theta_n^\star)\right\rbrace
\sqrt n(\overline\theta_n-\theta_0)\\
&-\Psi_{\theta,n}(\overline\theta_n)^{-1}\sqrt n\Psi_n(\theta_0).
\end{aligned}
$$

Consistency at the $\sqrt n$ rate places both random evaluation points near $\theta_0$ and keeps the preliminary error $O_{\mathbb{P}}(1)$ after scaling. Uniform local derivative convergence makes the bracket $o_{\mathbb{P}}(1)$ and the inverse derivative converge to $[\Psi_\theta(\theta_0)]^{-1}$. This gives the claimed representation. Pointwise derivative convergence does not control either random evaluation point.

[Back to the solution map](#solution-map)

<a id="problem-6"></a>

## Problem 6. Misspecified likelihood needs the sandwich

[Problem statement](../topic-09-m-z-estimation.md#problem-6)

The expected Poisson log criterion, up to a constant, is

$$
M(\theta)=\mu_0\log\theta-\theta.
$$

It is uniquely maximized at $\theta=\mu_0$, and the sample maximizer is $\overline W_n$. For $\psi(W,\theta)=W/\theta-1$,

$$
\Psi_\theta(\mu_0)
=\mathbb E[-W/\mu_0^2]
=-\frac1{\mu_0},
\qquad
B=\mathbb V(W/\mu_0-1)=\frac{\tau^2}{\mu_0^2}.
$$

Therefore

$$
[\Psi_\theta(\mu_0)]^{-2}B=\tau^2,
$$

matching the ordinary CLT for the sample mean. The negative expected Hessian of the log likelihood is $1/\mu_0$, while the score variance is $\tau^2/\mu_0^2$; they agree exactly when $\tau^2=\mu_0$. Under overdispersion $\tau^2>\mu_0$, imposing inverse Poisson information reports asymptotic variance $\mu_0$ instead of $\tau^2$ and understates uncertainty.

[Back to the solution map](#solution-map)

<a id="problem-7"></a>

## Problem 7. A curved Gaussian likelihood laboratory

[Problem statement](../topic-09-m-z-estimation.md#problem-7)

Write $X=\theta(1+Z)$ with $Z\sim\mathcal N(0,1)$. Using $\mathbb E[Z]=\mathbb E[Z^3]=0$, $\mathbb E[Z^2]=1$, and $\mathbb E[Z^4]=3$ gives

$$
\mathbb E[X]=\theta,
\qquad
\mathbb E[X^2]=2\theta^2,
\qquad
\mathbb E[X^3]=4\theta^3,
\qquad
\mathbb E[X^4]=10\theta^4.
$$

Up to a constant independent of $\theta$, the sample log likelihood is

$$
\ell_n(\theta)
=-n\log\theta-\frac{nT_n}{2\theta^2}+\frac{nS_n}{\theta}-\frac n2.
$$

Factorization proves that $(S_n,T_n)$ is sufficient. For two samples $x$ and $y$, their likelihood ratio is independent of $\theta$ exactly when

$$
\sum_i x_i=\sum_i y_i
\qquad\text{and}\qquad
\sum_i x_i^2=\sum_i y_i^2.
$$

Indeed, the log ratio is a constant plus the first difference divided by $\theta$ minus half the second difference divided by $\theta^2$. A function $a/\theta+b/\theta^2$ is constant on $(0,\infty)$ only when $a=b=0$. The likelihood ratio criterion therefore proves minimal sufficiency.

Differentiation gives

$$
\ell_n'(\theta)
=\frac{n}{\theta^3}\lbrace T_n-S_n\theta-\theta^2\rbrace.
$$

On $T_n>0$, the quadratic equation has exactly one positive root,

$$
\widehat\theta_n
=\frac{\sqrt{S_n^2+4T_n}-S_n}{2}.
$$

The other root is negative. The score is positive near zero and negative for sufficiently large $\theta$, while the log likelihood tends to $-\infty$ at both ends, so the positive root is the unique global maximum. If $T_n=0$, every observation is zero. On that probability zero sample the likelihood is unbounded as $\theta\downarrow0$, so no maximizer exists in the open parameter space.

The law of large numbers gives $(S_n,T_n)\to_{\mathbb{P}}(\theta,2\theta^2)$. The estimator is a continuous function of this pair at the probability limit, and

$$
\frac{\sqrt{\theta^2+8\theta^2}-\theta}{2}=\theta,
$$

so $\widehat\theta_n\to_{\mathbb{P}}\theta$.

The raw moments imply

$$
\mathbb V(X)=\theta^2,
\qquad
\mathrm{Cov}(X,X^2)=4\theta^3-2\theta^3=2\theta^3,
$$

and

$$
\mathbb V(X^2)=10\theta^4-4\theta^4=6\theta^4.
$$

Thus the multivariate CLT gives

$$
\sqrt n\left\lbrace
\begin{pmatrix}
S_n\\
T_n
\end{pmatrix}
-\begin{pmatrix}
\theta\\
2\theta^2
\end{pmatrix}
\right\rbrace
\rightsquigarrow
\mathcal N\mkern-3mu\left(
0,
\begin{pmatrix}
\theta^2&2\theta^3\\
2\theta^3&6\theta^4
\end{pmatrix}
\right).
$$

For $g(s,t)=(\sqrt{s^2+4t}-s)/2$,

$$
\nabla g(\theta,2\theta^2)
=\begin{pmatrix}
-1/3\\
1/(3\theta)
\end{pmatrix}.
$$

Multiplying the covariance matrix on both sides by this gradient gives $\theta^2/3$, proving

$$
\sqrt n(\widehat\theta_n-\theta)
\rightsquigarrow
\mathcal N\mkern-3mu\left(0,\frac{\theta^2}{3}\right).
$$

For one observation, the score is

$$
s_\theta(X)
=-\frac1\theta-\frac{X}{\theta^2}+\frac{X^2}{\theta^3}
=\frac{Z^2+Z-1}{\theta}.
$$

Because $\mathbb E[(Z^2+Z-1)^2]=3$, the Fisher information per observation is $I(\theta)=3/\theta^2$. Its inverse is $\theta^2/3$, which agrees with the MLE asymptotic variance.

A feasible standard error for $\widehat\theta_n$ is

$$
\widehat{\mathrm{se}}(\widehat\theta_n)
=\frac{\widehat\theta_n}{\sqrt{3n}}.
$$

The direct Wald interval is

$$
\left[
\widehat\theta_n-z_{1-\alpha/2}\frac{\widehat\theta_n}{\sqrt{3n}},
\widehat\theta_n+z_{1-\alpha/2}\frac{\widehat\theta_n}{\sqrt{3n}}
\right].
$$

Since $\sqrt n(\log\widehat\theta_n-\log\theta)\rightsquigarrow\mathcal N(0,1/3)$, a positivity preserving interval is

$$
\left[
\widehat\theta_n\exp\mkern-3mu\left(-\frac{z_{1-\alpha/2}}{\sqrt{3n}}\right),
\widehat\theta_n\exp\mkern-3mu\left(\frac{z_{1-\alpha/2}}{\sqrt{3n}}\right)
\right].
$$

Both intervals have pointwise asymptotic coverage $1-\alpha$ for each fixed $\theta>0$. The second interval always respects the parameter space, while the first can have a negative lower endpoint in small samples.

[Back to the solution map](#solution-map)
