# ORF 524 Practice Module 14: Solutions

**Status:** Released solutions  
**Last updated:** October 2, 2026  
**Practice module:** [Module 14](../topic-14-semiparametric-methods.md)  
**Primary chapter:** [Week 14](../../lectures/week-14.md)

Attempt each problem before consulting its solution. After studying a solution, close it and reconstruct the argument, including the assumptions used at each step.

<a id="solution-map"></a>

## Solution map

| Problem | Solution |
|---|---|
| [1. Influence functions and transformations](../topic-14-semiparametric-methods.md#problem-1) | [Solution](#problem-1) |
| [2. Integrated squared density and the Hoeffding decomposition](../topic-14-semiparametric-methods.md#problem-2) | [Solution](#problem-2) |
| [3. Orthogonal estimation in the partially linear model](../topic-14-semiparametric-methods.md#problem-3) | [Solution](#problem-3) |
| [4. Missing outcomes: IPW, first step adjustment, and augmentation](../topic-14-semiparametric-methods.md#problem-4) | [Solution](#problem-4) |
| [5. Weighted average derivatives](../topic-14-semiparametric-methods.md#problem-5) | [Solution](#problem-5) |
| [6. Cross fitting is not a rate assumption](../topic-14-semiparametric-methods.md#problem-6) | [Solution](#problem-6) |

<a id="problem-1"></a>

## Problem 1. Influence functions and transformations

[Problem statement](../topic-14-semiparametric-methods.md#problem-1)

Let $p_t$ be a smooth submodel with $p_0=p$ and score $s$, so that $\mathbb E[s(W)]=0$. Differentiating the mean gives

$$
\left.\frac{\partial}{\partial t}\mathbb E_t[W]\right|_{t=0} =
\mathbb E[Ws(W)] =
\mathbb E[(W-\mu)s(W)].
$$

Thus $\varphi_\mu(W)=W-\mu$ represents the derivative.

For the variance, use $\sigma^2=\mathbb E[W^2]-\mu^2$. Its derivative is

$$
\mathbb E[W^2s(W)]-2\mu\mathbb E[Ws(W)] =
\mathbb E[{(W-\mu)^2-\sigma^2}s(W)],
$$

where the term $-\sigma^2$ can be added because $\mathbb E[s]=0$. Hence

$$
\varphi_{\sigma^2}(W) =
(W-\mu)^2-\sigma^2.
$$

This derivation accounts for estimation of the mean. Treating the mean as fixed reaches the same first order expression only because the derivative of $\mathbb E[(W-a)^2]$ with respect to $a$ vanishes at $a=\mu$.

The chain rule gives

$$
\varphi_\tau(W) =
a'(\mu)(W-\mu).
$$

The sample mean has asymptotic variance $\mathbb V(W)$, the sample variance has asymptotic variance $\mathbb V(\lbrace W-\mu\rbrace^2)$, and the smooth transformation has asymptotic variance $a'(\mu)^2\mathbb V(W)$. These are precisely the corresponding delta method formulas.

A finite second moment suffices for the ordinary mean CLT. A finite fourth moment is a simple sufficient condition for the variance influence function to have finite variance. The smooth transformation again needs a finite second moment and differentiability of $a$. If $a'(\mu)=0$, the first order influence function is zero and a second order expansion may be needed for a nondegenerate limit.

[Back to the solution map](#solution-map)

<a id="problem-2"></a>

## Problem 2. Integrated squared density and the Hoeffding decomposition

[Problem statement](../topic-14-semiparametric-methods.md#problem-2)

Expanding the squared KDE and interchanging its finite sums with integration gives

$$
\int\widehat f_h(x)^2dx
=\frac1{n^2}\sum_{i=1}^n\sum_{j=1}^n
\int K_h(x-X_i)K_h(x-X_j)dx.
$$

With the change of variables $u=(x-X_i)/h$, the integral in each summand is

$$
h^{-d}\int K(u)K\mkern-3mu\left(u+\frac{X_i-X_j}{h}\right)du
=L_h(X_i-X_j).
$$

This proves the plug in V-statistic identity. Deleting the diagonal and rescaling gives the stated U-statistic.

Conditioning on $X_1=x$ now gives

$$
\mathbb E[L_h(x-X_2)] =
\int L_h(x-u)f(u)du
=(L_h*f)(x)=g_h(x).
$$

Therefore

$$
\theta_h =
\mathbb E[g_h(X_1)] =
\int f(x)(L_h*f)(x)dx.
$$

If $L$ integrates to one and its moments below order $p$ vanish, a Taylor expansion of $f(x-hu)$ gives

$$
(L_h*f)(x)-f(x)=O(h^p)
$$

in a norm strong enough to integrate against $f$. Thus

$$
\theta_h-\theta =
\int f(x)\lbrace(L_h*f)(x)-f(x)\rbrace dx
=O(h^p).
$$

Define

$$
q_h(x,y) =
L_h(x-y)-g_h(x)-g_h(y)+\theta_h.
$$

It is degenerate because $\mathbb E[q_h(x,X_2)]=0$. Adding and subtracting the two first projections inside the U-statistic gives

$$
\widehat\theta_h-\theta_h =
\frac2n\sum_{i=1}^n\lbrace g_h(X_i)-\theta_h\rbrace
+
\binom{n}{2}^{-1}\sum_{i<j}q_h(X_i,X_j).
$$

The second term is $R_{n,h}$.

The same convolution expansion that gives the bias also gives, under the stated smoothness in $L^2(\mathbb{P})$,

$$
\Vert g_h-f\Vert_{L^2(\mathbb{P})}=O(h^p),
\qquad
|\theta_h-\theta|=O(h^p).
$$

Consequently,

$$
\mathbb E\left[
\left\lbrace
2\lbrace g_h(X)-\theta_h\rbrace
-2\lbrace f(X)-\theta\rbrace
\right\rbrace^2
\right]
=O(h^{2p}).
$$

To bound the degenerate kernel, first note that boundedness of $f$ and square integrability of $L$ imply

$$
\begin{aligned}
\mathbb E[L_h(X_1-X_2)^2]
&=h^{-2d}\iint
L\mkern-3mu\left(\frac{x-y}{h}\right)^2
f(x)f(y)\mkern3mu dx\mkern3mu dy\\
&=h^{-d}\iint L(t)^2f(y+ht)f(y)\mkern3mu dt\mkern3mu dy
=O(h^{-d}).
\end{aligned}
$$

The remaining terms in $q_h$ have bounded second moments, so $\mathbb E[q_h(X_1,X_2)^2]=O(h^{-d})$. Degeneracy makes covariances between distinct pairs vanish, including pairs that share one index. Hence

$$
\mathbb E[R_{n,h}^2]
=\binom{n}{2}^{-1}\mathbb E[q_h(X_1,X_2)^2]
=O(n^{-2}h^{-d}).
$$

The deterministic bias, first projection, and degenerate projection are mutually orthogonal after centering. Therefore

$$
\begin{aligned}
\mathbb E[(\widehat\theta_h-\theta)^2]
&=(\theta_h-\theta)^2
+\frac4n\mathbb V\lbrace g_h(X)\rbrace
+\mathbb E[R_{n,h}^2]\\
&=O\mkern-3mu\left(
h^{2p}+\frac1n+\frac1{n^2h^d}
\right).
\end{aligned}
$$

For asymptotic linearity, the degenerate term satisfies

$$
\sqrt n R_{n,h} =
O_{\mathbb{P}}\mkern-3mu\left(\frac1{\sqrt{nh^d}}\right)
=o_{\mathbb{P}}(1)
$$

if $nh^d\to\infty$. The deterministic bias obeys

$$
\sqrt n(\theta_h-\theta)=O(\sqrt n h^p),
$$

so require $\sqrt n h^p\to0$.

For $h=n^{-a}$, the two requirements become

$$
\frac1{2p}<a<\frac1d.
$$

Such an exponent exists exactly when $p>d/2$. Thus the displayed $\sqrt n$ rate argument implicitly requires enough smoothness and a sufficiently high order kernel relative to dimension; stating the two limits separately does not guarantee that a bandwidth can satisfy both.

If $h\to0$ and $g_h\to f$ in $L^2(\mathbb{P})$, the triangular array first projection approaches

$$
\varphi(X)=2\lbrace f(X)-\theta\rbrace.
$$

Under a suitable Lindeberg condition,

$$
\sqrt n(\widehat\theta_h-\theta)
\rightsquigarrow
\mathcal N\mkern-3mu\left(0,4\mathbb V\lbrace f(X)\rbrace\right).
$$

For feasible inference, define

$$
\widehat g_{-i,h}(X_i) =
\frac1{n-1}\sum_{j\ne i}L_h(X_i-X_j),
\qquad
\widehat\varphi_i =
2\lbrace\widehat g_{-i,h}(X_i)-\widehat\theta_h\rbrace.
$$

Then

$$
\widehat V=\mathbb{E}_n[\widehat\varphi^2],
\qquad
\mathrm{se}(\widehat\theta_h)=\sqrt{\widehat V/n}.
$$

Consistency requires the same smoothing, moment, and leave one out approximation controls used for the linear representation.

For the V-statistic with the same kernel, symmetry gives

$$
\widehat\theta_h^V =
\frac1{n^2}\sum_{i\ne j}L_h(X_i-X_j)
+
\frac1{n^2}\sum_iL_h(0) =
\frac{n-1}{n}\widehat\theta_h
+
\frac{L(0)}{nh^d}.
$$

Making the diagonal term negligible after multiplication by $\sqrt n$ requires $\sqrt n h^d\to\infty$, stronger than $nh^d\to\infty$. For $h=n^{-a}$, the V-statistic's bias and diagonal requirements are $1/(2p)<a<1/(2d)$, which are compatible exactly when $p>d$.

[Back to the solution map](#solution-map)

<a id="problem-3"></a>

## Problem 3. Orthogonal estimation in the partially linear model

[Problem statement](../topic-14-semiparametric-methods.md#problem-3)

Conditional expectation yields

$$
\ell_0(X)=\theta_0m_0(X)+g_0(X).
$$

With $V=D-m_0(X)$, the original equation and conditional mean restrictions give

$$
\mathbb E[VY]=\theta_0\mathbb E[VD]
$$

and $\mathbb E[VD]=\mathbb E[V^2]$. Hence

$$
\theta_0=\frac{\mathbb E[VY]}{\mathbb E[VD]}.
$$

Residualizing $Y$ as well gives

$$
Y-\ell_0(X)=\theta_0V+\varepsilon.
$$

Since $\mathbb E[V\varepsilon]=0$,

$$
\theta_0 =
\frac{\mathbb E[V\lbrace Y-\ell_0(X)\rbrace]}{\mathbb E[V^2]}
$$

provided $J_0=\mathbb E[V^2]>0$.

For a common series matrix $P$, the residual maker $M_P=I-P(P'P)^{-1}P'$ gives

$$
\widehat\theta=\frac{D'M_PY}{D'M_PD}.
$$

This is the Frisch--Waugh--Lovell form of the joint series regression. Local polynomial nuisance regressions implement the same population residualization through local fits rather than a global basis. Either first step may be trained outside the evaluation fold.

At the truth the score is

$$
\psi(W;\theta_0,m_0,\ell_0)=V\varepsilon,
$$

which has mean zero. For square integrable perturbations $\delta\ell$ and $\delta m$, the two directional derivatives of the population score are

$$
-\mathbb E[V\delta\ell(X)]=0
$$

and

$$
\mathbb E[-\delta m(X)\varepsilon+\theta_0V\delta m(X)]=0.
$$

The first equality uses $\mathbb E[V\mid X]=0$; the second additionally uses $\mathbb E[\varepsilon\mid X]=0$.

The derivative with respect to $\theta$ is $-J_0$. Z-estimator linearization therefore gives

$$
\varphi_0(W)=J_0^{-1}V\varepsilon
$$

and

$$
\mathbb V(\varphi_0) =
J_0^{-2}\mathbb E[V^2\varepsilon^2].
$$

No homoskedastic factorization is justified without an additional assumption.

To find the nuisance remainder, put $\widehat m=m_0+\delta m$ and $\widehat\ell=\ell_0+\delta\ell$. Then at $\theta_0$,

$$
\psi(W;\theta_0,\widehat m,\widehat\ell) =
(V-\delta m)
(\varepsilon-\delta\ell+\theta_0\delta m).
$$

Taking expectations and eliminating all first order terms by conditional means gives

$$
\mathbb E[\psi(W;\theta_0,\widehat m,\widehat\ell)] =
\mathbb E[\delta m\mkern3mu\delta\ell]
-\theta_0\mathbb E[(\delta m)^2].
$$

Consequently,

$$
|\mathbb E[\psi]|
\le
\Vert\delta m\Vert_2\Vert\delta\ell\Vert_2
+|\theta_0|\Vert\delta m\Vert_2^2.
$$

A sufficient condition is

$$
\Vert\delta m\Vert_2
\lbrace\Vert\delta\ell\Vert_2+\Vert\delta m\Vert_2\rbrace
=o_{\mathbb{P}}(n^{-1/2}).
$$

For example, both errors being $o_{\mathbb{P}}(n^{-1/4})$ suffices.

For $K$ fold cross fitting, train $(\widehat m_{-k},\widehat\ell_{-k})$ outside fold $I_k$ and, for $i\in I_k$, form

$$
\widehat V_i=D_i-\widehat m_{-k}(X_i),
\qquad
\widehat U_i=Y_i-\widehat\ell_{-k}(X_i).
$$

Then

$$
\widehat\theta =
\frac{\sum_i\widehat V_i\widehat U_i}{\sum_i\widehat V_i^2}.
$$

Let

$$
\widehat\varepsilon_i =
\widehat U_i-\widehat\theta\widehat V_i,
\qquad
\widehat J=\mathbb{E}_n[\widehat V^2],
\qquad
\widehat\varphi_i=\widehat J^{-1}\widehat V_i\widehat\varepsilon_i.
$$

The estimated standard error is

$$
\sqrt{\mathbb{E}_n[\widehat\varphi^2]/n}.
$$

Orthogonality removes first order population sensitivity to nuisance error. Cross fitting separates nuisance training from score evaluation and helps control overfitting dependence. Nuisance accuracy makes the second order product negligible. The condition $J_0>0$ identifies the coefficient and keeps the score derivative invertible. None of these properties implies the others.

[Back to the solution map](#solution-map)

<a id="problem-4"></a>

## Problem 4. Missing outcomes: IPW, first step adjustment, and augmentation

[Problem statement](../topic-14-semiparametric-methods.md#problem-4)

Under independence of $G$ and $T$,

$$
\mathbb E[TG]=\mathbb E[T]\mathbb E[G]=p_0\theta_0.
$$

On the event $\sum_iT_i>0$,

$$
\widetilde\theta
=\frac{\sum_iT_iG_i}{\sum_iT_i},
$$

so it is exactly the average among observed outcomes. If $\mathbb E|G|<\infty$, the LLN gives $\mathbb E_n[TG]\to_{\mathbb{P}}p_0\theta_0$ and $\widehat p\to_{\mathbb{P}}p_0>0$, establishing consistency by Slutsky's theorem.

For the distributional result, expand the ratio around $(p_0\theta_0,p_0)$:

$$
\begin{aligned}
\sqrt n(\widetilde\theta-\theta_0)
&=\frac1{p_0\sqrt n}\sum_{i=1}^n
\left[
T_iG_i-p_0\theta_0-\theta_0(T_i-p_0)
\right]+o_{\mathbb{P}}(1)\\
&=\frac1{\sqrt n}\sum_{i=1}^n
\frac{T_i}{p_0}(G_i-\theta_0)+o_{\mathbb{P}}(1).
\end{aligned}
$$

Thus

$$
\varphi(W)=\frac{T}{p_0}(G-\theta_0),
\qquad
\Omega=\mathbb E[\varphi(W)^2]
=\frac{\mathbb V(G)}{p_0},
$$

where the last equality uses independence. A consistent estimator and standard error are

$$
\widehat\Omega
=\mathbb{E}_n\mkern-3mu\left[
\frac{T(G-\widetilde\theta)^2}{\widehat p^2}
\right],
\qquad
\widehat{\mathrm{se}}(\widetilde\theta)
=\sqrt{\widehat\Omega/n}.
$$

Consequently, a feasible level $1-\alpha$ interval is $\widetilde\theta\pm z_{1-\alpha/2}\sqrt{\widehat\Omega/n}$.

If $p_0$ were known, the oracle estimator $\mathbb{E}_n[TG]/p_0$ would have influence function

$$
\varphi_{\mathrm{or}}(W)
=\frac{TG}{p_0}-\theta_0
=\frac{T}{p_0}(G-\theta_0)
+\frac{\theta_0}{p_0}(T-p_0).
$$

Under independence, the two terms on the right are uncorrelated. Estimating $p_0$ through the ratio removes the second term, so the ratio estimator has variance smaller by $\theta_0^2(1-p_0)/p_0$. This is the same first step adjustment logic as a control variate calculation.

Under missing at random,

$$
\mathbb E\mkern-3mu\left[
\left.\frac{TG}{p_0(X)}\right|X
\right]
=\frac{p_0(X)}{p_0(X)}\mathbb E[G\mid T=1,X]
=m_0(X).
$$

Iterated expectation gives $\mathbb E[m_0(X)]=\theta_0$, proving the conditional IPW identity. To establish consistency with an estimated propensity, decompose

$$
\widehat\theta_{\mathrm{IPW}}-\theta_0
=\mathbb{E}_n\mkern-3mu\left[
TG\left\lbrace\frac1{\widehat p(X)}-\frac1{p_0(X)}\right\rbrace
\right]
+\mathbb{E}_n\mkern-3mu\left[\frac{TG}{p_0(X)}\right]-\theta_0.
$$

If $p_0\geq c>0$ and $\Vert\widehat p-p_0\Vert_\infty=o_{\mathbb{P}}(1)$, then $\inf_x\widehat p(x)\geq c/2$ with probability tending to one. On that event, the first term is bounded in absolute value by

$$
\frac{2\Vert\widehat p-p_0\Vert_\infty}{c^2}
\mathbb{E}_n[T|G|]
=o_{\mathbb{P}}(1),
$$

while an integrable envelope gives an LLN for the second term.

For the parametric propensity model, define the row vector

$$
A_0'
=-\mathbb E\mkern-3mu\left[
\frac{TG\dot p(X;\gamma_0)'}{p_0(X)^2}
\right].
$$

A Taylor expansion and the supplied first step representation give

$$
\sqrt n(\widehat\theta_{\mathrm{IPW}}-\theta_0)
=\frac1{\sqrt n}\sum_{i=1}^n
\left[
\frac{T_iG_i}{p_0(X_i)}-\theta_0
+A_0'\varphi_\gamma(W_i)
\right]+o_{\mathbb{P}}(1).
$$

Thus the bracketed term is the first step adjusted influence function and its variance $\Sigma$ is the asymptotic variance. If $\widehat\varphi_{\gamma,i}$ consistently estimates the first step influence value, define

$$
\widehat A'
=-\mathbb{E}_n\mkern-3mu\left[
\frac{TG\dot p(X;\widehat\gamma)'}{p(X;\widehat\gamma)^2}
\right],
$$

$$
\widehat\psi_i
=\frac{T_iG_i}{p(X_i;\widehat\gamma)}
-\widehat\theta_{\mathrm{IPW}}
+\widehat A'\widehat\varphi_{\gamma,i},
\qquad
\widehat\Sigma=\mathbb{E}_n[\widehat\psi^2].
$$

The feasible standard error is $\sqrt{\widehat\Sigma/n}$, and the corresponding Gaussian interval uses this standard error.

For augmentation, both identifying identities have already been established. For generic $m$ and $p$, conditioning on $X$ gives

$$
\mathbb E\mkern-3mu\left[
\left.\frac{T}{p(X)}\lbrace G-m(X)\rbrace\right|X
\right]
=\frac{p_0(X)}{p(X)}\lbrace m_0(X)-m(X)\rbrace.
$$

Therefore

$$
\mathbb E[\psi(W;\theta_0,m,p)]
=\mathbb E\mkern-3mu\left[
\lbrace m(X)-m_0(X)\rbrace
\left\lbrace1-\frac{p_0(X)}{p(X)}\right\rbrace
\right].
$$

The expression is zero if $m=m_0$ or if $p=p_0$. This is a population identity. If one component is correct but estimated very noisily, finite sample behavior and variance can still be poor.

At the truth, perturbations $\delta m$ and $\delta p$ produce derivatives

$$
\mathbb E\mkern-3mu\left[
\delta m(X)\left\lbrace1-\frac{T}{p_0(X)}\right\rbrace
\right]=0
$$

and

$$
-\mathbb E\mkern-3mu\left[
\frac{T\lbrace G-m_0(X)\rbrace}{p_0(X)^2}\delta p(X)
\right]=0.
$$

Thus the score is orthogonal. Since $1-p_0/p=(p-p_0)/p$, positivity bounds the population remainder by a constant times $\Vert m-m_0\Vert_2\Vert p-p_0\Vert_2$. With cross fitting, a representative sufficient rate condition is

$$
\Vert\widehat m-m_0\Vert_2
\Vert\widehat p-p_0\Vert_2
=o_{\mathbb{P}}(n^{-1/2}),
$$

together with consistency of each nuisance, propensities bounded away from zero with high probability, a finite variance influence function, a Lindeberg condition, and consistent variance estimation.

The augmented influence function is

$$
\varphi_0(W)
=m_0(X)-\theta_0
+\frac{T}{p_0(X)}\lbrace G-m_0(X)\rbrace.
$$

For cross fitting, train $(\widehat m_{-k},\widehat p_{-k})$ outside fold $I_k$ and set

$$
\widehat\theta_{\mathrm{AIPW}}
=\frac1n\sum_{k=1}^K\sum_{i\in I_k}
\left[
\widehat m_{-k}(X_i)
+\frac{T_iG_i-T_i\widehat m_{-k}(X_i)}{\widehat p_{-k}(X_i)}
\right].
$$

Estimated influence values equal the bracketed summand minus $\widehat\theta_{\mathrm{AIPW}}$, and the standard error is the square root of their sample second moment divided by $n$. If $p_0(x)=0$ on a set of positive probability, the outcome regression on that set is not identified from observed outcomes without more structure. Propensities near zero also create large inverse weights, possibly infinite influence function variance, and numerical instability.

[Back to the solution map](#solution-map)

<a id="problem-5"></a>

## Problem 5. Weighted average derivatives

[Problem statement](../topic-14-semiparametric-methods.md#problem-5)

The general target is

$$
\theta_w =
\int \mu'(x)w(x)f(x)\mkern3mu dx.
$$

Suppose $\mu$, $w$, and $f$ are sufficiently differentiable and integrable and the limits of $\mu(x)w(x)f(x)$ at the support endpoints are zero. Integration by parts yields

$$
\theta_w
=\left[\mu(x)w(x)f(x)\right]_{\partial\mathcal X}
-\int \mu(x)\lbrace w'(x)f(x)+w(x)f'(x)\rbrace\mkern3mu dx.
$$

Taking $w=f$ gives, with the boundary term displayed before it is set to zero,

$$
\theta_f =
\left[\mu(x)f(x)^2\right]_{\partial\mathcal X}
-\int \mu(x)\lbrace f(x)^2\rbrace'dx =
-2\int \mu(x)f(x)f'(x)dx.
$$

Since $\mu(X)=\mathbb E[Y\mid X]$,

$$
\theta_f=-2\mathbb E[Yf'(X)].
$$

For a differentiable kernel, define

$$
\widehat f_{-i,h}'(X_i) =
\frac1{(n-1)h^2}\sum_{j\ne i}
K'\mkern-3mu\left(\frac{X_i-X_j}{h}\right)
$$

and

$$
\widehat\theta_{f,h} =
-\frac2n\sum_iY_i\widehat f_{-i,h}'(X_i).
$$

If $K$ is symmetric, its derivative is odd. Pairing ordered terms gives

$$
\widehat\theta_{f,h} =
\binom{n}{2}^{-1}\sum_{i<j}u_h(W_i,W_j),
$$

where

$$
u_h(W_i,W_j) =
-\frac1{h^2}(Y_i-Y_j)
K'\mkern-3mu\left(\frac{X_i-X_j}{h}\right).
$$

Exchanging $i$ and $j$ changes the signs of both differences, so this kernel is symmetric. Deleting $i$ prevents a self interaction and exposes the pairwise U-statistic structure. A $\sqrt n$ rate proof first centers at the smoothed target, applies a Hoeffding projection to isolate an empirical average, bounds the degenerate pairwise remainder, and imposes a separate undersmoothing condition on deterministic bias.

For $X\in\mathbb R^d$, $\mu'(x)$ becomes $\nabla \mu(x)$ and $f'(x)$ becomes $\nabla f(x)$. A directional target $a'\mathbb E[f(X)\nabla \mu(X)]$ uses the directional derivatives $a'\nabla \mu$ and $a'\nabla f$, together with boundary conditions on the relevant faces or tails.

[Back to the solution map](#solution-map)

<a id="problem-6"></a>

## Problem 6. Cross fitting is not a rate assumption

[Problem statement](../topic-14-semiparametric-methods.md#problem-6)

The two terms are negligible relative to $n^{-1/2}$ when

$$
a+b>\frac12
\qquad\text{and}\qquad
2a>\frac12.
$$

Equivalently, require $a>1/4$ and $a+b>1/2$ for exact polynomial orders. If $a=b=1/4$, both terms are only $O(n^{-1/2})$, not $o(n^{-1/2})$. The equal rate statement becomes sufficient if both nuisance errors are $o_{\mathbb{P}}(n^{-1/4})$.

For example, $a=0.30$ and $b=0.21$ works: $2a=0.60$ and $a+b=0.51$. The pair $a=0.20$ and $b=0.40$ fails because the term involving squared $m$ error has exponent $0.40$ despite the product having exponent $0.60$.

Cross fitting controls the dependence created when a flexible nuisance learner is evaluated on its own training observations. It does not assert that either prediction error converges at any particular rate. Valid partially linear inference also needs, among other conditions, residual variation $J_0$ bounded away from zero, finite moments and a CLT for the influence function, consistent variance estimation, and suitable stability of the cross fitted denominator.

[Back to the solution map](#solution-map)
