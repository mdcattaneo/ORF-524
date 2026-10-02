# ORF 524 Practice Module 12: Solutions

**Status:** Released solutions  
**Last updated:** October 2, 2026  
**Practice module:** [Module 12](../topic-12-series-nonparametric-methods.md)  
**Primary chapter:** [Week 12](../../lectures/week-12.md)

Attempt each problem before consulting its solution. After studying a solution, close it and reconstruct the argument, including the assumptions used at each step.

<a id="solution-map"></a>

## Solution map

| Problem | Solution |
|---|---|
| [1. Projection, approximation, and integrated rates](../topic-12-series-nonparametric-methods.md#problem-1) | [Solution](#problem-1) |
| [2. Pointwise series inference](../topic-12-series-nonparametric-methods.md#problem-2) | [Solution](#problem-2) |
| [3. Leave one out cross validation for series least squares](../topic-12-series-nonparametric-methods.md#problem-3) | [Solution](#problem-3) |
| [4. Partitioning regression as a series estimator](../topic-12-series-nonparametric-methods.md#problem-4) | [Solution](#problem-4) |
| [5. Basis reparameterization and numerical conditioning](../topic-12-series-nonparametric-methods.md#problem-5) | [Solution](#problem-5) |
| [6. One regularization principle, two methods](../topic-12-series-nonparametric-methods.md#problem-6) | [Solution](#problem-6) |

<a id="problem-1"></a>

## Problem 1. Projection, approximation, and integrated rates

[Problem statement](../topic-12-series-nonparametric-methods.md#problem-1)

The projection first order condition is

$$
\mathbb E[p_k(X)\lbrace \mu(X)-p_k(X)'\beta_k\rbrace] =
\mathbb E[p_k(X)r_k(X)]
=0.
$$

Since $Y=p_k(X)'\beta_k+r_k(X)+\varepsilon$, the sample normal equation gives

$$
\widehat Q_k(\widehat\beta_k-\beta_k) =
\mathbb{E}_n[p_k(X)\lbrace\varepsilon+r_k(X)\rbrace],
$$

which yields the exact decomposition.

For any $g_k\in\mathcal G_k$, write $\mu-g_k=(\mu-\mu_k)+(\mu_k-g_k)$. Projection orthogonality gives

$$
\mathbb E[(\mu(X)-\mu_k(X))(\mu_k(X)-g_k(X))]=0,
$$

and therefore

$$
\Vert \mu-g_k\Vert_{L^2(\mathbb{P}_X)}^2 =
\Vert \mu-\mu_k\Vert_{L^2(\mathbb{P}_X)}^2
+
\Vert \mu_k-g_k\Vert_{L^2(\mathbb{P}_X)}^2.
$$

Taking $g_k=\widehat\mu_k$ is valid for each realized coefficient vector because $\widehat\mu_k-\mu_k$ remains in $\mathcal G_k$. Thus

$$
\Vert\widehat\mu_k-\mu\Vert_{L^2(\mathbb{P}_X)}^2 =
\Vert\widehat\mu_k-\mu_k\Vert_{L^2(\mathbb{P}_X)}^2
+
\Vert r_k\Vert_{L^2(\mathbb{P}_X)}^2.
$$

Because the population mean of the vector summand is zero and its total second moment is of order $k$, a vector mean square or concentration argument gives

$$
\left\Vert
\mathbb{E}_n[p_k(X)\lbrace\varepsilon+r_k(X)\rbrace]
\right\Vert =
O_{\mathbb{P}}\mkern-3mu\left(\sqrt{\frac{k}{n}}\right).
$$

If $\Vert\widehat Q_k^{-1}\Vert_{\mathrm{op}}=O_{\mathbb{P}}(1)$, the coefficient rate follows.

With $Q_k=I_k$,

$$
\Vert p_k'(\widehat\beta_k-\beta_k)\Vert_{L^2(\mathbb{P}_X)} =
\Vert\widehat\beta_k-\beta_k\Vert.
$$

The triangle inequality then gives the claimed $L^2$ rate. Squaring its order yields $k/n+k^{-2\alpha}$. Balancing terms gives

$$
k\asymp n^{1/(2\alpha+1)}
$$

and

$$
\Vert\widehat\mu_k-\mu\Vert_{L^2(\mathbb{P}_X)}^2 =
O_{\mathbb{P}}\lbrace n^{-2\alpha/(2\alpha+1)}\rbrace.
$$

If $\alpha=s/d$, the dimension and rate specialize to

$$
k\asymp n^{d/(2s+d)},
\qquad
\Vert\widehat\mu_k-\mu\Vert_{L^2(\mathbb{P}_X)}^2 =
O_{\mathbb{P}}\mkern-3mu\left(n^{-2s/(2s+d)}\right).
$$

A supremum norm result additionally needs uniform approximation control and a basis size factor such as $\zeta_k=\sup_x\Vert p_k(x)\Vert$; an integrated coefficient rate alone does not control the largest pointwise error.

[Back to the solution map](#solution-map)

<a id="problem-2"></a>

## Problem 2. Pointwise series inference

[Problem statement](../topic-12-series-nonparametric-methods.md#problem-2)

The exact decomposition gives

$$
\widehat\mu_k(x)-\mu_k(x) =
p_k(x)'\widehat Q_k^{-1}
\mathbb{E}_n[p_k(X)u_k].
$$

Replacing $\widehat Q_k$ by $Q_k$ at first order gives influence term

$$
\varphi_{k,x}(W) =
p_k(x)'Q_k^{-1}p_k(X)u_k.
$$

Its variance is

$$
V_k(x) =
p_k(x)'Q_k^{-1}\Omega_kQ_k^{-1}p_k(x).
$$

With residuals $\widehat u_{ik}=Y_i-p_k(X_i)'\widehat\beta_k$, define

$$
\widehat\Omega_k =
\mathbb{E}_n[p_k(X)p_k(X)'\widehat u_k^2]
$$

and

$$
\widehat V_k(x) =
p_k(x)'\widehat Q_k^{-1}
\widehat\Omega_k
\widehat Q_k^{-1}p_k(x).
$$

The standard error is $\sqrt{\widehat V_k(x)/n}$.

A representative primitive sufficient set is the following. The eigenvalues of $Q_k$ are uniformly bounded above and away from zero, $\mathbb E[u_k^4\mid X]$ is uniformly bounded, $V_k(x)\geq c\Vert p_k(x)\Vert^2$ at the evaluation point, and $\sup_t\Vert p_k(t)\Vert\leq\zeta_k$ with $\zeta_k^2k\log k/n\to0$. These deliberately strong conditions make the empirical Gram replacement negligible and imply a Lindeberg condition by keeping normalized leverage small. More generally, one may replace them with direct conditions that $\widehat Q_k^{-1}$ is close enough to $Q_k^{-1}$ and

$$
\frac{\max_{1\leq i\leq n}\left|p_k(x)'Q_k^{-1}p_k(X_i)u_{k,i}\right|}{\sqrt{nV_k(x)}}\to_{\mathbb{P}}0.
$$

Under either formulation,

$$
\frac{\sqrt n\lbrace\widehat\mu_k(x)-\mu_k(x)\rbrace}{\sqrt{V_k(x)}}
\rightsquigarrow\mathcal N(0,1).
$$

To recenter at $\mu(x)$, require

$$
\frac{|r_k(x)|}{\sqrt{V_k(x)/n}}\to0.
$$

For fixed $k$, $r_k(x)$ generally does not vanish, so the ordinary OLS interval targets the projection value.

[Back to the solution map](#solution-map)

<a id="problem-3"></a>

## Problem 3. Leave one out cross validation for series least squares

[Problem statement](../topic-12-series-nonparametric-methods.md#problem-3)

Put $A=P'P$, let $p_i=p_k(X_i)$ be a column vector, and let $b=P'Y$. The coefficient after deleting observation $i$ is

$$
\widehat\beta_{-i} =
(A-p_ip_i')^{-1}(b-p_iY_i).
$$

Sherman--Morrison gives

$$
(A-p_ip_i')^{-1} =
A^{-1}
+
\frac{A^{-1}p_ip_i'A^{-1}}{1-p_i'A^{-1}p_i}.
$$

Let $h_i=p_i'A^{-1}p_i=H_{ii}$ and $\widehat Y_i=p_i'A^{-1}b$. Algebra then gives

$$
p_i'\widehat\beta_{-i} =
\frac{\widehat Y_i-h_iY_i}{1-h_i}.
$$

Hence

$$
Y_i-p_i'\widehat\beta_{-i} =
\frac{Y_i-\widehat Y_i}{1-H_{ii}}
$$

and

$$
\mathrm{CV}(k) =
\frac1n\sum_i
\left(
\frac{Y_i-\widehat Y_i}{1-H_{ii}}
\right)^2.
$$

Training error is optimistic because the same outcome affects both the fit and its residual; high leverage observations can be fit especially closely. Training data fit candidate models, validation data or folds choose among them, and an untouched test set estimates final performance. Cross validation targets prediction loss. Pointwise inference additionally needs approximation bias small relative to the standard error, which prediction optimal $k$ need not ensure.

For density estimation,

$$
\int\lbrace\widehat f_h(x)-f(x)\rbrace^2\mkern3mu dx =
\int\widehat f_h(x)^2\mkern3mu dx
-2\mathbb E[\widehat f_h(X)]
+\int f(x)^2\mkern3mu dx.
$$

The final term is constant in $h$. Replacing the middle expectation by leave one out evaluation gives

$$
\int\widehat f_h(x)^2\mkern3mu dx
-\frac2n\sum_{i=1}^n\widehat f_{-i,h}(X_i).
$$

Unlike regression cross validation, there is no response residual; the criterion estimates integrated squared density error up to an additive constant.

[Back to the solution map](#solution-map)

<a id="problem-4"></a>

## Problem 4. Partitioning regression as a series estimator

[Problem statement](../topic-12-series-nonparametric-methods.md#problem-4)

Let $P_j$ be the $j$-th bin and take

$$
p_J(x) =
(\mathbf 1\lbrace x\in P_1\rbrace,\ldots,\mathbf 1\lbrace x\in P_J\rbrace)'.
$$

Least squares separates by bins, and each fitted coefficient is the sample mean in its bin.

Each count is marginally $\mathsf{Binomial}(n,1/J)$. A union bound gives

$$
\mathbb{P}\left(\min_jN_j=0\right)
\leq\sum_{j=1}^J\mathbb{P}(N_j=0)
=J(1-1/J)^n
\leq J e^{-n/J}.
$$

The logarithm of the last bound is $\log J-n/J$, which tends to negative infinity because $J\log J/n\to0$. Thus all bins are nonempty with probability tending to one.

For $Z_{ij}=\mathbf 1\lbrace X_i\in P_j\rbrace-1/J$, Bernstein's inequality gives, for any $t>0$,

$$
\mathbb{P}\mkern-3mu\left(
\left|\frac1n\sum_iZ_{ij}\right|>t
\right)
\leq
2\exp\mkern-3mu\left(
-\frac{nt^2}{2/J+2t/3}
\right).
$$

Take

$$
t=C\left(
\sqrt{\frac{\log J}{nJ}}+\frac{\log J}{n}
\right)
$$

and apply a union bound over $j$. For sufficiently large $C$, the resulting probability tends to zero. Moreover,

$$
Jt
=C\left(
\sqrt{\frac{J\log J}{n}}+\frac{J\log J}{n}
\right)
\longrightarrow0,
$$

which proves the stated $o_{\mathbb{P}}(J^{-1})$ conclusion and displays exactly where the growth condition enters.

On the nonempty bin event, conditional independence and homoskedasticity give

$$
\mathbb V\lbrace\widehat\mu_J(x)\mid\mathbf X\rbrace
=\sum_{j=1}^J
\mathbf 1\lbrace x\in P_j\rbrace\frac{\sigma^2}{N_j}.
$$

Therefore

$$
\int_0^1\mathbb V\lbrace\widehat\mu_J(x)\mid\mathbf X\rbrace\mkern3mu dx
=\frac{\sigma^2}{J}\sum_{j=1}^J\frac1{N_j}.
$$

The uniform count result implies $N_j=(n/J)\lbrace1+o_{\mathbb{P}}(1)\rbrace$ uniformly in $j$, and hence

$$
\int_0^1\mathbb V\lbrace\widehat\mu_J(x)\mid\mathbf X\rbrace\mkern3mu dx
=\frac{\sigma^2J}{n}\lbrace1+o_{\mathbb{P}}(1)\rbrace.
$$

Let

$$
\overline \mu_j =
J\int_{P_j}\mu(t)dt
$$

be the population projection coefficient in bin $j$, and let $c_j$ be the bin midpoint. A uniform first order expansion gives

$$
\mu(x)-\overline \mu_j =
\mu'(c_j)(x-c_j)+o(J^{-1})
$$

for $x\in P_j$ under the stated smoothness. The linear term has zero bin average, which is why replacing $\mu(c_j)$ by the actual bin projection does not change the leading term. Since

$$
\int_{P_j}(x-c_j)^2dx=\frac1{12J^3},
$$

summing yields leading integrated squared approximation error

$$
\frac1{12J^2}\int_0^1\mu'(x)^2dx
$$

under Riemann sum regularity.

To account for the random design within each bin, put

$$
\widetilde\mu_j=
\frac1{N_j}\sum_{i:X_i\in P_j}\mu(X_i).
$$

Because $\int_{P_j}\lbrace\overline\mu_j-\mu(x)\rbrace dx=0$, the cross term in

$$
\int_{P_j}\lbrace\widetilde\mu_j-\mu(x)\rbrace^2dx
$$

vanishes exactly. Summing over bins gives

$$
\int_0^1
\left[
\mathbb E\lbrace\widehat\mu_J(x)\mid\mathbf X\rbrace-\mu(x)
\right]^2dx =
\int_0^1\lbrace\mu_J(x)-\mu(x)\rbrace^2dx
+\frac1J\sum_{j=1}^J(\widetilde\mu_j-\overline\mu_j)^2.
$$

Conditional on $N_j$, the observations within bin $j$ are uniform on that bin. Uniform continuity of $\mu'$ implies $\mathbb V\lbrace\mu(X)\mid X\in P_j\rbrace=O(J^{-2})$ uniformly in $j$. Since $N_j\asymp n/J$ uniformly, averaging over the bins gives

$$
\frac1J\sum_{j=1}^J(\widetilde\mu_j-\overline\mu_j)^2
=O_{\mathbb{P}}\mkern-3mu\left(\frac1{nJ}\right)
=o_{\mathbb{P}}(J^{-2}),
$$

where the final relation uses $J/n\to0$.

Writing

$$
\mathcal B=\frac1{12}\int_0^1\mu'(x)^2dx,
$$

the conditional IMSE expansion is

$$
\mathrm{IMSE}(J\mid\mathbf X)
=\sigma^2\frac{J}{n}\lbrace1+o_{\mathbb{P}}(1)\rbrace
+\frac{\mathcal B}{J^2}\lbrace1+o_{\mathbb{P}}(1)\rbrace.
$$

Treating $J$ as continuous, differentiating the leading expression gives

$$
J_{\mathrm{IMSE}}
=\left(\frac{2\mathcal B}{\sigma^2}\right)^{1/3}n^{1/3},
\qquad
\mathrm{IMSE}_{\mathrm{opt}}
\asymp n^{-2/3}.
$$

The displayed optimal constant assumes $\mathcal B>0$. If $\mu$ is constant, the leading approximation term vanishes and this balance no longer determines the appropriate growth rate.

Data dependent boundaries make the basis random. Their selection error and dependence on outcomes or covariates must then be included; fixed partition calculations no longer apply automatically.

[Back to the solution map](#solution-map)

<a id="problem-5"></a>

## Problem 5. Basis reparameterization and numerical conditioning

[Problem statement](../topic-12-series-nonparametric-methods.md#problem-5)

Since $A_k$ is nonsingular,

$$
\lbrace\widetilde p_k(x)'\gamma:\gamma\in\mathbb R^k\rbrace =
\lbrace p_k(x)'\beta:\beta\in\mathbb R^k\rbrace.
$$

In sample, the two design matrices have the same column space, so their orthogonal projection matrices and fitted values coincide when both have full column rank. Coefficients differ by the inverse linear transformation.

Eigenvalues and condition numbers are coordinate dependent. A poor transformation can make the Gram matrix numerically ill conditioned even though the function space is unchanged. For theory, one often uses a population normalization satisfying $\mathbb E[\widetilde p_k(X)\widetilde p_k(X)']=I_k$ or uniformly bounded eigenvalues.

[Back to the solution map](#solution-map)

<a id="problem-6"></a>

## Problem 6. One regularization principle, two methods

[Problem statement](../topic-12-series-nonparametric-methods.md#problem-6)

Kernel balance gives

$$
h\asymp n^{-1/(2p+1)},
\qquad
\mathrm{MSE}\asymp n^{-2p/(2p+1)}.
$$

Series balance gives

$$
k\asymp n^{1/(2\alpha+1)},
\qquad
\mathrm{MSE}\asymp n^{-2\alpha/(2\alpha+1)}.
$$

Smaller bandwidth uses fewer nearby observations and permits more local variation. More series terms enlarge the approximation space. For kernels, leave observation $i$ out, predict $Y_i$ using each candidate $h$, and minimize held out squared residuals. For series, do the same over candidate $k$, using the hat matrix shortcut.

Prediction optimal tuning can leave bias of standard error order, and data driven tuning introduces selection effects not accounted for by a fixed tuning standard error. Either issue can invalidate a naive pointwise interval.

[Back to the solution map](#solution-map)
