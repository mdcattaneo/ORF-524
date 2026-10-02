# ORF 524 Practice Module 1: Solutions

**Status:** Released solutions  
**Last updated:** October 2, 2026  
**Practice module:** [Module 1](../topic-01-probability-prediction.md)  
**Primary chapter:** [Week 1](../../lectures/week-01.md)

Attempt each problem before consulting its solution. After studying a solution, close it and reconstruct the argument, including the assumptions used at each step.

<a id="solution-map"></a>

## Solution map

| Problem | Solution |
|---|---|
| [1. Loss determines the target](../topic-01-probability-prediction.md#problem-1) | [Solution](#problem-1) |
| [2. Inequalities and their assumptions](../topic-01-probability-prediction.md#problem-2) | [Solution](#problem-2) |
| [3. Conditional expectation as projection](../topic-01-probability-prediction.md#problem-3) | [Solution](#problem-3) |
| [4. The predictor class matters](../topic-01-probability-prediction.md#problem-4) | [Solution](#problem-4) |
| [5. An integrability boundary for iterated expectations](../topic-01-probability-prediction.md#problem-5) | [Solution](#problem-5) |
| [6. Continuity of probability](../topic-01-probability-prediction.md#problem-6) | [Solution](#problem-6) |

<a id="problem-1"></a>

## Problem 1. Loss determines the target

[Problem statement](../topic-01-probability-prediction.md#problem-1)

### 1. Squared loss

Let $\mu=\mathbb E[Y]$. Since $Y\in L^2$,

$$
\begin{aligned}
\mathbb E[(Y-a)^2]
&=\mathbb E[(Y-\mu+\mu-a)^2]\\
&=\mathbb E[(Y-\mu)^2]
+2(\mu-a)\mathbb E[Y-\mu]+(\mu-a)^2\\
&=\mathbb V[Y]+(\mu-a)^2.
\end{aligned}
$$

The second term is nonnegative and vanishes only at $a=\mu$. Hence the unique minimizer is $\mathbb E[Y]$.

### 2. Absolute loss

Write $R(a)=\mathbb E[|Y-a|]$. Its left and right derivatives are

$$
R_{-}'(a)=2F_Y(a^-)-1,
\qquad
R_{+}'(a)=2F_Y(a)-1.
$$

The convex function $R$ is minimized at $a$ exactly when $R_{-}'(a)\leq0\leq R_{+}'(a)$, or

$$
F_Y(a^-)\leq\frac12\leq F_Y(a).
$$

These are precisely the medians of $Y$. Integrability of $Y$ ensures that $R(a)$ is finite for every finite $a$.

### 3. Check loss

For fixed $y$ and $b>a$, direct consideration of the positions of $y,a,b$ gives

$$
\rho_\tau(y-b)-\rho_\tau(y-a)
=\int_a^b\lbrace\mathbf 1(y\leq t)-\tau\rbrace\mkern3mu dt.
$$

The integrand is bounded, so Fubini's theorem gives

$$
R_\tau(b)-R_\tau(a)
=\int_a^b\lbrace F_Y(t)-\tau\rbrace\mkern3mu dt,
\qquad
R_\tau(a)=\mathbb E[\rho_\tau(Y-a)].
$$

Equivalently, the left and right derivatives are

$$
R'_{\tau,-}(a)=F_Y(a^-)-\tau,
\qquad
R'_{\tau,+}(a)=F_Y(a)-\tau.
$$

Convexity therefore implies that $a$ minimizes $R_\tau$ exactly when

$$
F_Y(a^-)\leq\tau\leq F_Y(a).
$$

### 4. Two point distribution

Here $\mathbb E[Y]=1$, so the squared loss minimizer is $1$. For absolute loss,

$$
F_Y(a^-)\leq\frac12\leq F_Y(a)
$$

holds for every $a\in[0,2]$. For $\tau=1/4$, it holds only at $a=0$; for $\tau=3/4$, it holds only at $a=2$. Thus optimality has no unique meaning until the loss is specified, and a specified loss still need not produce a unique optimizer.

[Back to the solution map](#solution-map)

<a id="problem-2"></a>

## Problem 2. Inequalities and their assumptions

[Problem statement](../topic-01-probability-prediction.md#problem-2)

### 1. Monotonicity of $L^p$

Because $|x|^p\leq1+|x|^q$, the assumed $q$-th moment implies that $X\in L^p$. Let $W=|X|^p$ and $r=q/p>1$. Jensen's inequality for the convex function $u\mapsto u^r$ gives

$$
(\mathbb E[|X|^p])^{q/p}
\leq\mathbb E[|X|^q].
$$

Taking $q$-th roots yields $\Vert X\Vert_p\leq\Vert X\Vert_q$. Jensen uses expectation with respect to a probability measure: the total mass is one. For a general finite measure, a factor depending on the total mass appears.

### 2. Mean--median inequality

Because $m$ minimizes expected absolute deviation,

$$
\begin{aligned}
|\mu-m|
&=|\mathbb E[X-m]|\\
&\leq\mathbb E[|X-m|]\\
&\leq\mathbb E[|X-\mu|]\\
&\leq\sqrt{\mathbb E[(X-\mu)^2]}=\sigma.
\end{aligned}
$$

The steps use the triangle inequality for expectation, median optimality under absolute loss, and Cauchy--Schwarz (or Jensen) respectively. The finite second moment implies the needed first moments.

### 3. Markov and equality

Pointwise,

$$
W\geq t\mathbf 1\lbrace W\geq t\rbrace.
$$

Taking expectations and dividing by $t$ proves the bound. Equality holds exactly when the nonnegative difference $W-t\mathbf 1\lbrace W\geq t\rbrace$ is zero almost surely. Equivalently,

$$
\mathbb{P}(W\in\lbrace0,t\rbrace)=1.
$$

This includes the degenerate cases at $0$ or $t$.

### 4. Chebyshev bound for the sample mean

Set $\overline X_n=n^{-1}\sum_iX_i$. Independence and common finite variance give

$$
\mathbb V[\overline X_n]
=\frac1{n^2}\sum_{i=1}^n\mathbb V[X_i]
=\frac{\sigma^2}{n}.
$$

Apply Markov's inequality to the nonnegative random variable $(\overline X_n-\mu)^2$ at threshold $\varepsilon^2$:

$$
\mathbb{P}(|\overline X_n-\mu|\geq\varepsilon)
\leq\frac{\mathbb E[(\overline X_n-\mu)^2]}{\varepsilon^2}
=\frac{\sigma^2}{n\varepsilon^2}.
$$

Independence is used only to make cross covariances zero. It can be replaced by pairwise uncorrelatedness (together with the stated common means and variances).

[Back to the solution map](#solution-map)

<a id="problem-3"></a>

## Problem 3. Conditional expectation as projection

[Problem statement](../topic-01-probability-prediction.md#problem-3)

### 1. Conditional and total variance

Using that $m(X)$ is measurable with respect to $X$,

$$
\begin{aligned}
\mathbb V(Y\mid X)
&=\mathbb E[(Y-m(X))^2\mid X]\\
&=\mathbb E[Y^2\mid X]-2m(X)\mathbb E[Y\mid X]+m(X)^2\\
&=\mathbb E[Y^2\mid X]-m(X)^2.
\end{aligned}
$$

Taking expectations and using iterated expectations,

$$
\begin{aligned}
\mathbb E[\mathbb V(Y\mid X)]
&=\mathbb E[Y^2]-\mathbb E[m(X)^2],\\
\mathbb V[m(X)]
&=\mathbb E[m(X)^2]-\mathbb E[m(X)]^2.
\end{aligned}
$$

Since $\mathbb E[m(X)]=\mathbb E[Y]$, summing the two displays gives $\mathbb V[Y]$.

### 2. Projection identity

Write

$$
Y-g(X)=Y-m(X)+m(X)-g(X).
$$

The cross term vanishes because

$$
\begin{aligned}
\mathbb E[(Y-m(X))(m(X)-g(X))]
&=\mathbb E\mkern-3mu\left[(m(X)-g(X))
\mathbb E[Y-m(X)\mid X]\right]\\
&=0.
\end{aligned}
$$

Expanding the square proves the identity. The final term is nonnegative and equals zero exactly when $g(X)=m(X)$ almost surely, proving uniqueness in $L^2$ up to almost sure equality.

### 3. Perfect prediction

The assumptions imply $X,Y\in L^2$. Moreover,

$$
\mathbb E[XY]
=\mathbb E\mkern-3mu\left[X\mathbb E[Y\mid X]\right]
=\mathbb E[X^2].
$$

Therefore

$$
\mathbb E[(Y-X)^2]
=\mathbb E[Y^2]-2\mathbb E[XY]+\mathbb E[X^2]
=\mathbb E[Y^2]-\mathbb E[X^2]=0.
$$

A nonnegative random variable with expectation zero is zero almost surely, so $Y=X$ almost surely.

### 4. Necessity of both conditions

- Let $X=0$ and let $Y$ be a Rademacher random variable, taking $-1$ and $1$ with equal probabilities. Then $\mathbb E[Y\mid X]=0=X$, but $\mathbb E[Y^2]=1\neq0=\mathbb E[X^2]$ and $X\neq Y$ almost surely.
- Let $X$ and $Y$ be independent Rademacher random variables. Then $\mathbb E[X^2]=\mathbb E[Y^2]=1$, but $\mathbb E[Y\mid X]=0\neq X$ almost surely and $\mathbb{P}(X\neq Y)=1/2$.

[Back to the solution map](#solution-map)

<a id="problem-4"></a>

## Problem 4. The predictor class matters

[Problem statement](../topic-01-probability-prediction.md#problem-4)

### 1. Unrestricted predictor

The conditional mean is

$$
\mathbb E[Y\mid X]=X^2+\mathbb E[\varepsilon\mid X]=X^2.
$$

It is therefore the best predictor among square integrable functions of $X$. Its risk is

$$
\mathbb E[(Y-X^2)^2]
=\mathbb E[\varepsilon^2]
=\mathbb E[\mathbb E[\varepsilon^2\mid X]]
=\sigma^2.
$$

### 2. Best affine predictor

For $R(a,b)=\mathbb E[(Y-a-bX)^2]$, the normal equations are

$$
\mathbb E[Y-a-bX]=0,
\qquad
\mathbb E[X(Y-a-bX)]=0.
$$

For $X\sim\mathsf{Uniform}(-1,1)$,

$$
\mathbb E[X]=0,
\quad
\mathbb E[X^2]=\frac13,
\quad
\mathbb E[X^3]=0,
\quad
\mathbb E[X^4]=\frac15.
$$

Also,

$$
\mathbb E[\varepsilon]=0,
\qquad
\mathbb E[X\varepsilon]
=\mathbb E[X\mathbb E[\varepsilon\mid X]]=0.
$$

Hence $\mathbb E[Y]=1/3$ and $\mathbb E[XY]=0$. The normal equations give $a=1/3$ and $b=0$. The best affine predictor is the constant $1/3$.

### 3. Affine risk and excess risk

The cross term between $\varepsilon$ and $X^2-1/3$ vanishes conditionally on $X$. Thus

$$
\begin{aligned}
\mathbb E[(Y-1/3)^2]
&=\mathbb E[\varepsilon^2]
+\mathbb E[(X^2-1/3)^2]\\
&=\sigma^2+\left(\frac15-\frac19\right)\\
&=\sigma^2+\frac4{45}.
\end{aligned}
$$

The excess risk relative to unrestricted prediction is $4/45$.

### 4. Interpretation

The loss did not change, but the feasible set did. The unrestricted class contains the true conditional mean $X^2$; the affine class does not contain it (up to almost sure equality under this continuous distribution). Restricting the predictor class therefore changes the optimizer and raises the minimum risk.

[Back to the solution map](#solution-map)

<a id="problem-5"></a>

## Problem 5. An integrability boundary for iterated expectations

[Problem statement](../topic-01-probability-prediction.md#problem-5)

### 1. Density of the ratio

Use the transformation $(z,v)=(u/v,v)$, whose inverse is $(u,v)=(zv,v)$ and whose absolute Jacobian is $|v|$. Then

$$
\begin{aligned}
f_Z(z)
&=\int_{-\infty}^{\infty}
\frac{|v|}{2\pi}
\exp\mkern-3mu\left\lbrace-\frac{(1+z^2)v^2}{2}\right\rbrace\mkern3mu dv\\
&=\frac1\pi\int_0^\infty
v\exp\mkern-3mu\left\lbrace-\frac{(1+z^2)v^2}{2}\right\rbrace\mkern3mu dv\\
&=\frac1{\pi(1+z^2)}.
\end{aligned}
$$

### 2. Nonexistence of the mean

By symmetry,

$$
\mathbb E[|Z|]
=\frac2\pi\int_0^\infty\frac{z}{1+z^2}\mkern3mu dz
=\infty.
$$

The positive and negative parts of $Z$ both have infinite expectation. Their formal cancellation does not define an expectation, so $\mathbb E[Z]$ is undefined rather than zero.

### 3. Why the tower property is unavailable

For a fixed $v\neq0$, $Z\mid V=v$ has the distribution of $U/v$ and hence has mean zero. But

$$
\int \mathbb E[|Z|\mid V=v],d\mathbb{P}_V(v)=\mathbb E[|Z|]=\infty.
$$

Thus $Z$ is not integrable, and the ordinary $L^1$ conditional expectation $\mathbb E[Z\mid V]$ and the tower property for it are not available. A collection of finite pointwise signed means does not repair the failure of absolute integrability needed to interchange integrations.

[Back to the solution map](#solution-map)

<a id="problem-6"></a>

## Problem 6. Continuity of probability

[Problem statement](../topic-01-probability-prediction.md#problem-6)

### 1. Continuity from below

Let $B_1=A_1$ and $B_n=A_n\setminus A_{n-1}$ for $n\geq2$. The $B_n$ are disjoint, $A_n=\bigcup_{k=1}^nB_k$, and $A=\bigcup_{k=1}^\infty B_k$. Countable additivity gives

$$
\mathbb{P}(A_n)=\sum_{k=1}^n\mathbb{P}(B_k)
\longrightarrow
\sum_{k=1}^\infty\mathbb{P}(B_k)=\mathbb{P}(A).
$$

### 2. Continuity from above

Since $A_n\downarrow A$, the complements satisfy $A_n^c\uparrow A^c$. Therefore

$$
\mathbb{P}(A_n)
=1-\mathbb{P}(A_n^c)
\longrightarrow
1-\mathbb{P}(A^c)=\mathbb{P}(A).
$$

The subtraction uses $\mathbb{P}(\Omega)=1<\infty$. For a general measure, continuity from above requires that at least one set in the decreasing sequence have finite measure.

### 3. Countable union bound

For each $N$, the finite union bound gives

$$
\mathbb{P}\mkern-3mu\left(\bigcup_{n=1}^N A_n\right)
\leq\sum_{n=1}^N\mathbb{P}(A_n).
$$

The unions on the left increase to $\bigcup_{n=1}^\infty A_n$. By continuity from below and then monotonicity of partial sums,

$$
\mathbb{P}\mkern-3mu\left(\bigcup_{n=1}^\infty A_n\right)
\leq\sum_{n=1}^\infty\mathbb{P}(A_n),
$$

with the right side allowed to be infinite.

[Back to the solution map](#solution-map)
