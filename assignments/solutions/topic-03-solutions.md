# ORF 524 Practice Module 3: Solutions

**Status:** Released solutions  
**Last updated:** October 2, 2026  
**Practice module:** [Module 3](../topic-03-decision-point-estimation.md)  
**Primary chapter:** [Week 3](../../lectures/week-03.md)

Attempt each problem before consulting its solution. After studying a solution, close it and reconstruct the argument, including the assumptions used at each step.

<a id="solution-map"></a>

## Solution map

| Problem | Solution |
|---|---|
| [1. Risk depends on the parameter space](../topic-03-decision-point-estimation.md#problem-1) | [Solution](#problem-1) |
| [2. Rao--Blackwellization in the Bernoulli model](../topic-03-decision-point-estimation.md#problem-2) | [Solution](#problem-2) |
| [3. Unbiasedness versus MSE in the Uniform model](../topic-03-decision-point-estimation.md#problem-3) | [Solution](#problem-3) |
| [4. Information, reparameterization, and a boundary failure](../topic-03-decision-point-estimation.md#problem-4) | [Solution](#problem-4) |
| [5. A covariance characterization of UMVU](../topic-03-decision-point-estimation.md#problem-5) | [Solution](#problem-5) |
| [6. Unbiased does not imply admissible in dimensions three and higher](../topic-03-decision-point-estimation.md#problem-6) | [Solution](#problem-6) |

<a id="problem-1"></a>

## Problem 1. Risk depends on the parameter space

[Problem statement](../topic-03-decision-point-estimation.md#problem-1)

Since $\overline X_n\sim\mathcal N(\theta,1/n)$,

$$
R(\theta,\delta_0)=\frac1n,
\qquad
R(\theta,\delta_c)=\frac{c^2}{n}+(1-c)^2\theta^2.
$$

For $0<c<1$, $R(0,\delta_c)=c^2/n<1/n$, but the second risk diverges as $|\theta|\to\infty$. The risks cross, so neither rule dominates on $\mathbb R$.

On $[-M,M]$, risk of $\delta_c$ is largest at $|\theta|=M$. With $c_M=nM^2/(1+nM^2)$,

$$
\sup_{|\theta|\leq M}R(\theta,\delta_{c_M})
=\frac{M^2}{1+nM^2}<\frac1n.
$$

Thus $\delta_{c_M}$ strictly dominates $\delta_0$ on the restricted space. Changing the parameter space changes the set of risk inequalities required for dominance. Admissibility is also relative to the class of competing rules.

For the prior $\theta\sim\mathcal N(0,v)$, Normal conjugacy gives

$$
\theta\mid X_1,\ldots,X_n
\sim
\mathcal N\left(\frac{nv}{1+nv}\overline X_n,\frac{v}{1+nv}\right).
$$

Thus the posterior mean is the Bayes rule under squared loss, and its Bayes risk is the expected posterior variance $v/(1+nv)$. For any rule $\delta$,

$$
\sup_{\theta\in\mathbb R}R(\theta,\delta)
\geq r(\Pi_v,\delta)
\geq \frac{v}{1+nv}.
$$

Letting $v\to\infty$ makes the lower bound tend to $1/n$. Since $R(\theta,\overline X_n)=1/n$ for every $\theta$, $\overline X_n$ attains this lower bound and is minimax. The argument uses only proper Normal priors.

[Back to the solution map](#solution-map)

<a id="problem-2"></a>

## Problem 2. Rao--Blackwellization in the Bernoulli model

[Problem statement](../topic-03-decision-point-estimation.md#problem-2)

The likelihood factors through $T$:

$$
f_\theta(x)=\theta^T(1-\theta)^{n-T}.
$$

If $\mathbb E_\theta[g(T)]=0$ for all $\theta\in(0,1)$, divide

$$
\sum_{t=0}^ng(t)\binom nt\theta^t(1-\theta)^{n-t}=0
$$

by $(1-\theta)^n$ and set $z=\theta/(1-\theta)$. A polynomial vanishing for all $z>0$ has all coefficients zero, so $T$ is complete.

Because $T$ is complete and sufficient, Bahadur's theorem also makes it minimal sufficient, up to the usual common null set qualification.

Conditional on $T=t$, all binary sequences with $t$ successes are equally likely. Hence

$$
\mathbb E[X_1\mid T=t]=\frac tn,
$$

and

$$
\mathbb V_\theta(T/n)=\frac{\theta(1-\theta)}n
\leq\theta(1-\theta)=\mathbb V_\theta(X_1).
$$

The inequality is strict for $n>1$ and is equality for $n=1$.

Similarly,

$$
\mathbb E[X_1X_2\mid T=t]
=\frac{t(t-1)}{n(n-1)},
$$

because the probability that two labeled positions are both among the $t$ successes is $t(t-1)/[n(n-1)]$. The two Rao--Blackwellized estimators are therefore

$$
\frac Tn
\quad\text{and}\quad
\frac{T(T-1)}{n(n-1)}.
$$

Sufficiency makes the conditional expectations parameter free functions of $T$ and ensures the Rao--Blackwell variance comparison. Completeness makes each unbiased function of $T$ unique, so Lehmann--Scheffe gives the UMVU conclusions. Basu's theorem further says that $T$ is independent of every ancillary statistic $A$ in this model.

[Back to the solution map](#solution-map)

<a id="problem-3"></a>

## Problem 3. Unbiasedness versus MSE in the Uniform model

[Problem statement](../topic-03-decision-point-estimation.md#problem-3)

The likelihood factors as

$$
L_n(\theta;x)
=\theta^{-n}\mathbf 1\lbrace M_n\leq\theta\rbrace
\prod_{i=1}^n\mathbf 1\lbrace x_i\geq0\rbrace,
$$

so $M_n$ is sufficient. If $\mathbb E_\theta[g(M_n)]=0$ for every $\theta>0$, then

$$
\frac n{\theta^n}\int_0^\theta g(m)m^{n-1}\mkern3mu dm=0.
$$

Thus the indefinite integral is zero for every upper limit, and differentiation gives $g(\theta)\theta^{n-1}=0$ almost everywhere. Hence $M_n$ is complete. The unbiased function

$$
\widehat\theta_U=\frac{n+1}{n}M_n
$$

is therefore UMVU.

Using $\mathbb E[M_n^2]=n\theta^2/(n+2)$,

$$
R(\theta,kM_n)
=\theta^2\left\lbrace
\frac{n}{n+2}k^2-\frac{2n}{n+1}k+1
\right\rbrace.
$$

The unique minimizer is

$$
k^\star=\frac{n+2}{n+1},
\qquad
R(\theta,k^\starM_n)=\frac{\theta^2}{(n+1)^2}.
$$

The first and second raw moment equations give

$$
\widehat\theta_{\mathrm{MM},1}=2\overline X_n,
\qquad
\widehat\theta_{\mathrm{MM},2}=\sqrt{\frac3n\sum_{i=1}^nX_i^2}.
$$

Both population maps, $\theta\mapsto\theta/2$ and $\theta\mapsto\theta^2/3$, are one-to-one on $(0,\infty)$, so each identifies $\theta$. Their empirical analogues are generally different, illustrating that method of moments depends on the selected moment condition.

The relevant comparisons are

$$
\begin{array}{c|c|c}
\text{estimator}&\text{bias}&\text{MSE}\\
\hline
M_n&-\theta/(n+1)&2\theta^2/\lbrace(n+1)(n+2)\rbrace\\
2\overline X_n&0&\theta^2/(3n)\\
\frac{n+1}{n}M_n&0&\theta^2/\lbrace n(n+2)\rbrace\\
\frac{n+2}{n+1}M_n&-\theta/(n+1)^2&\theta^2/(n+1)^2.
\end{array}
$$

Here $M_n$ is the likelihood maximizer on the probability one event $M_n>0$ and is an arbitrary estimator convention on the null sample. Also, $2\overline X_n$ follows from equating $\mathbb E[X]=\theta/2$ to the sample mean. The UMVU claim compares all unbiased estimators and improves on this method of moments estimator (strictly for $n>1$). The $k^\star$ estimator accepts bias to reduce MSE within $\mathcal D$, so the optimality statements concern different classes and criteria.

[Back to the solution map](#solution-map)

<a id="problem-4"></a>

## Problem 4. Information, reparameterization, and a boundary failure

[Problem statement](../topic-03-decision-point-estimation.md#problem-4)

The log likelihood is

$$
\ell_n(\theta)=-n\log\theta-\frac1\theta\sum_iX_i.
$$

Therefore

$$
s_\theta=-\frac n\theta+\frac{\sum_iX_i}{\theta^2},
\qquad
\mathbb E_\theta[s_\theta]=0,
\qquad
\mathcal I_n(\theta)=\frac n{\theta^2}.
$$

The CRLB for unbiased estimation of $\theta$ is $\theta^2/n$. Since $\mathbb V(\overline X_n)=\theta^2/n$, the sample mean attains it.

For $\eta=\log\theta$, the chain rule gives $s_\eta=\theta s_\theta$ and

$$
\mathcal I_n(\eta)=\theta^2\mathcal I_n(\theta)=n.
$$

The target $\tau(\eta)=e^\eta$ has derivative $e^\eta=\theta$, so its bound is again $\theta^2/n$. Although likelihood invariance gives $\widehat\eta=\log\overline X_n$, strict concavity of log and nondegeneracy give

$$
\mathbb E_\theta[\log\overline X_n]
<\log\mathbb E_\theta[\overline X_n]=\log\theta.
$$

For the Uniform likelihood, differentiating $-n\log\theta$ while ignoring the indicator produces $-n/\theta$, whose expectation is not zero. The support depends on $\theta$, so differentiation cannot be moved through the normalization integral as if the integration region were fixed.

For the generic two parameter information matrix,

$$
\mathcal I_n(\theta)^{-1}
=\frac1{ac-b^2}
\begin{pmatrix}
c&-b\\
-b&a
\end{pmatrix}.
$$

The nuisance adjusted bound for the first coordinate is therefore

$$
\frac{c}{ac-b^2}
=\frac1{a-b^2/c}
\geq \frac1a.
$$

It agrees with $1/a$ exactly when $b=0$, that is, when the two score coordinates are orthogonal.

[Back to the solution map](#solution-map)

<a id="problem-5"></a>

## Problem 5. A covariance characterization of UMVU

[Problem statement](../topic-03-decision-point-estimation.md#problem-5)

Let $D=U-T$. Every $T_c=T+cD$ is unbiased, and

$$
\mathbb V_\theta(T_c)
=\mathbb V_\theta(T)+c^2\mathbb V_\theta(D)
+2c\mathrm{Cov}_\theta(T,D).
$$

If the covariance is zero, setting $c=1$ shows $\mathbb V(U)=\mathbb V(T)+\mathbb V(D)\geq\mathbb V(T)$, so $T$ is UMVU. Conversely, if $T$ is UMVU, the displayed quadratic minus $\mathbb V(T)$ must be nonnegative for every $c$. If $\mathbb V(D)>0$, its minimum is $-\mathrm{Cov}(T,D)^2/\mathbb V(D)$, so the covariance must vanish. If $\mathbb V(D)=0$, then $D$ is constant almost surely and, because it has mean zero, is zero almost surely; its covariance with $T$ is again zero.

[Back to the solution map](#solution-map)

<a id="problem-6"></a>

## Problem 6. Unbiased does not imply admissible in dimensions three and higher

[Problem statement](../topic-03-decision-point-estimation.md#problem-6)

Write $r^2=\Vert X\Vert^2$. Expanding the risk gives

$$
R(\mu,\delta_a)-R(\mu,X)
=a^2\mathbb E[r^{-2}]
-2a\sum_{j=1}^m
\mathbb E[(X_j-\mu_j)X_j/r^2].
$$

For $g_j(x)=x_j/\Vert x\Vert^2$,

$$
\sum_{j=1}^m\partial_jg_j(x)
=\frac{m-2}{\Vert x\Vert^2}.
$$

The singularity at zero can be handled by applying Stein's identity first to $g_{j,\varepsilon}(x)=x_j/(\Vert x\Vert^2+\varepsilon)$ and then taking $\varepsilon\downarrow0$. The required inverse square expectation is finite when $m\geq3$. Stein's identity therefore yields

$$
R(\mu,\delta_a)-R(\mu,X)
=\lbrace a^2-2a(m-2)\rbrace\mathbb E_\mu[\Vert X\Vert^{-2}].
$$

The expectation is finite and positive for $m\geq3$. Strict improvement occurs for $0<a<2(m-2)$, and the quadratic coefficient is minimized at $a=m-2$.

The model has natural exponential family form

$$
f(x;\mu)=\phi_m(x)\exp\lbrace\mu'x-\Vert\mu\Vert^2/2\rbrace,
\qquad
\mu\in\mathbb R^m,
$$

so the natural parameter space is open and $X$ is complete and sufficient for $\mu$. Each $X_j$ is unbiased for $\mu_j$ and is therefore UMVU by Lehmann--Scheffe. This componentwise conclusion compares unbiased scalar estimators one coordinate at a time. The shrinkage rule is biased and is compared with $X$ as a vector valued rule under total squared error loss, a different class and criterion; componentwise UMVU therefore does not imply joint admissibility.

[Back to the solution map](#solution-map)
