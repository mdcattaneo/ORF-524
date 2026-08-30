# ORF 524 Practice Module 12: Series Based Nonparametric Methods

**Status:** Ready for Fall 2026  
**Last updated:** August 30, 2026  
**Primary chapter:** [Week 12: Series based nonparametric methods](../lectures/week-12.md)

## Purpose

This ungraded module develops series projection, estimation, inference, and tuning. Problems 1--4 form the core bank. Problems 5--6 examine basis invariance and the common regularization logic behind kernel and series methods.

<a id="problem-map"></a>

## Problem map

| Problem | Bank | Main task |
|---|---|---|
| [1. Projection, approximation, and integrated rates](#problem-1) | Core | Decompose projection, approximation, and stochastic errors and derive integrated rates. |
| [2. Pointwise series inference](#problem-2) | Core | Derive pointwise sandwich inference and identify the approximation bias condition. |
| [3. Leave one out cross validation for series least squares](#problem-3) | Core | Derive leave one out formulas for regression and density series estimators. |
| [4. Partitioning regression as a series estimator](#problem-4) | Core | Analyze occupancy, conditional bias and variance, and integrated MSE for partitions. |
| [5. Basis reparameterization and numerical conditioning](#problem-5) | Further | Separate mathematical basis invariance from numerical conditioning. |
| [6. One regularization principle, two methods](#problem-6) | Further | Compare kernel and series tuning, rates, cross validation, and inference. |

## Core practice

<a id="problem-1"></a>

### Problem 1. Projection, approximation, and integrated rates

**Chapter connections:** [Section 1.3: Population projection](../lectures/week-12.md#13-population-projection), [Section 2.1: Exact decomposition](../lectures/week-12.md#21-exact-decomposition), and [Section 3.1: Integrated rate and its source](../lectures/week-12.md#31-integrated-rate-and-its-source)

Let

$$
Y=\mu(X)+\varepsilon,
\qquad
\mathbb E[\varepsilon\mid X]=0,
$$

and suppose the basis is normalized so that $Q_k=\mathbb E[p_k(X)p_k(X)']=I_k$. Let $\beta_k$ be the population projection coefficient, $\mu_k(x)=p_k(x)'\beta_k$, and $r_k=\mu-\mu_k$.

1. Derive the projection normal equation and the exact decomposition

$$
\widehat\beta_k-\beta_k =
\widehat Q_k^{-1}\mathbb{E}_n[p_k(X)\lbrace\varepsilon+r_k(X)\rbrace].
$$

2. Prove that for every $g_k\in\mathcal G_k$,

$$
\Vert \mu-g_k\Vert_{L^2(\mathbb{P}_X)}^2 =
\Vert \mu-\mu_k\Vert_{L^2(\mathbb{P}_X)}^2
+
\Vert \mu_k-g_k\Vert_{L^2(\mathbb{P}_X)}^2.
$$

Apply the identity to $g_k=\widehat\mu_k$ and explain why it remains valid although $\widehat\mu_k$ is random.

3. Under suitable moment and Gram matrix conditions, show why

$$
\Vert\widehat\beta_k-\beta_k\Vert =
O_{\mathbb{P}}\mkern-3mu\left(\sqrt{\frac{k}{n}}\right).
$$

4. Deduce

$$
\Vert\widehat\mu_k-\mu\Vert_{L^2(\mathbb{P}_X)} =
O_{\mathbb{P}}\mkern-3mu\left(
\sqrt{\frac{k}{n}}+\Vert r_k\Vert_{L^2(\mathbb{P}_X)}
\right).
$$

5. If $\Vert r_k\Vert_{L^2(\mathbb{P}_X)}\lesssim k^{-\alpha}$, derive the integrated MSE optimal order of $k$ and the resulting squared error rate. Then specialize to smoothness $s$ in dimension $d$ by taking $\alpha=s/d$.

6. Explain why the same calculation does not by itself establish a supremum norm rate.

<!-- Source lineage: legacy Week 11 series slides; repairs the stochastic square-root scale and separates projection from approximation. -->

[Back to the problem map](#problem-map)

<a id="problem-2"></a>

### Problem 2. Pointwise series inference

**Chapter connections:** [Section 4.1: Fixed $k$ OLS inference for the projection](../lectures/week-12.md#41-fixed-k-ols-inference-for-the-projection), [Section 4.2: CLT plus approximation bias](../lectures/week-12.md#42-clt-plus-approximation-bias), and [Section 4.4: RBC inference in practice](../lectures/week-12.md#44-rbc-inference-in-practice)

Fix an evaluation point $x$. Let

$$
u_k=Y-p_k(X)'\beta_k,
\qquad
\Omega_k=\mathbb E[p_k(X)p_k(X)'u_k^2].
$$

1. Derive the leading linear representation for $\widehat\mu_k(x)-\mu_k(x)$.

2. Identify its asymptotic variance and construct a heteroskedasticity robust plug in estimator.

3. State a pointwise CLT for the projection $\mu_k(x)$, making explicit any leverage, moment, and growing $k$ conditions you invoke.

4. Give the additional approximation bias condition needed to use the same interval for $\mu(x)$.

5. Explain why a fixed $k$ OLS interval is generally an interval for $\mu_k(x)$ rather than for $\mu(x)$.

<!-- Source lineage: legacy sequential-asymptotics discussion; rewritten as a projection-versus-function inference audit. -->

[Back to the problem map](#problem-map)

<a id="problem-3"></a>

### Problem 3. Leave one out cross validation for series least squares

**Chapter connections:** [Section 5.1: Prediction risk](../lectures/week-12.md#51-prediction-risk) and [Section 5.2: Leave one out cross validation](../lectures/week-12.md#52-leave-one-out-cross-validation)

Let $P$ be the $n\times k$ series design matrix, assume $P'P$ is nonsingular, and define

$$
H=P(P'P)^{-1}P',
\qquad
\widehat Y=HY.
$$

Let $\widehat\mu_{-i,k}(X_i)$ be the prediction at $X_i$ obtained after deleting observation $i$.

1. Use the Sherman--Morrison identity to prove

$$
Y_i-\widehat\mu_{-i,k}(X_i) =
\frac{Y_i-\widehat Y_i}{1-H_{ii}},
$$

   whenever $H_{ii}\neq1$.

2. Derive the computational shortcut for $\mathrm{CV}(k)$.

3. Explain why training error is systematically optimistic relative to out of sample prediction error.

4. Distinguish the roles of training, validation, and test data.

5. Explain why a $k$ selected for prediction need not make series approximation bias negligible for pointwise inference.

6. Starting from integrated squared error for a density estimate, derive the leave one out density criterion

$$
\int\widehat f_h(x)^2\mkern3mu dx
-\frac2n\sum_{i=1}^n\widehat f_{-i,h}(X_i),
$$

up to a term that does not depend on $h$. Explain why this differs from regression cross validation.

<!-- Source lineage: Fall 2024 pset8, “Linear smoothers, cross-validation”; corrects the legacy slide's generic CV objective. -->

[Back to the problem map](#problem-map)

<a id="problem-4"></a>

### Problem 4. Partitioning regression as a series estimator

**Chapter connections:** [Section 1.2: Concrete bases](../lectures/week-12.md#12-concrete-bases), [Section 2.1: Exact decomposition](../lectures/week-12.md#21-exact-decomposition), and [Section 3.1: Integrated rate and its source](../lectures/week-12.md#31-integrated-rate-and-its-source)

Let $(Y_i,X_i)_{i=1}^n$ be i.i.d. copies of $(Y,X)$, where $X\sim\mathsf{Uniform}(0,1)$ and

$$
Y=\mu(X)+\varepsilon,
\qquad
\mathbb E[\varepsilon\mid X]=0,
\qquad
\mathbb V(\varepsilon\mid X)=\sigma^2.
$$

Partition $[0,1]$ into the $J$ equal width bins $P_1,\ldots,P_J$, let $N_j=\sum_i\mathbf 1\lbrace X_i\in P_j\rbrace$, and estimate $\mu(x)$ by the sample mean of the outcomes in the bin containing $x$, with an arbitrary value assigned to an empty bin. Assume $J\to\infty$, $J\log J/n\to0$, and $\mu$ is continuously differentiable with uniformly continuous derivative.

1. Express the estimator as series least squares using bin indicator basis functions.

2. Prove

$$
\mathbb{P}\left(\min_{1\leq j\leq J}N_j=0\right)
\leq J\left(1-\frac1J\right)^n
\leq J e^{-n/J}
\longrightarrow0.
$$

3. Use Bernstein's inequality and a union bound to show

$$
\max_{1\leq j\leq J}
\left|\frac{N_j}{n}-\frac1J\right|
=O_{\mathbb{P}}\mkern-3mu\left(
\sqrt{\frac{\log J}{nJ}}+\frac{\log J}{n}
\right)
=o_{\mathbb{P}}(J^{-1}).
$$

Explain exactly where $J\log J/n\to0$ is used.

4. Conditional on $\mathbf X=(X_1,\ldots,X_n)$ and on the event that every bin is nonempty, derive the exact integrated conditional variance and show

$$
\int_0^1\mathbb V\lbrace\widehat\mu_J(x)\mid\mathbf X\rbrace\mkern3mu dx
=\frac{\sigma^2J}{n}\lbrace1+o_{\mathbb{P}}(1)\rbrace.
$$

5. Let $\overline\mu_j=J\int_{P_j}\mu(t)dt$ and $\mu_J(x)=\sum_j\mathbf 1\lbrace x\in P_j\rbrace\overline\mu_j$. Decompose the integrated squared conditional bias into the population approximation error $\int_0^1\lbrace\mu_J(x)-\mu(x)\rbrace^2dx$ and the error from replacing each $\overline\mu_j$ by the sample average of $\mu(X_i)$ within its bin. Show that the cross term is zero and that the second term is $O_{\mathbb{P}}(1/(nJ))$.

6. Show that

$$
\int_0^1\lbrace\mu_J(x)-\mu(x)\rbrace^2dx
=\frac1{12J^2}\int_0^1\mu'(x)^2dx+o(J^{-2}).
$$

7. Combine the preceding results to obtain the conditional IMSE expansion, including its leading constants. Derive the IMSE optimal constant and order of $J$ and the resulting IMSE rate.

8. Explain which steps become more difficult with data dependent partition boundaries.

<!-- Source lineage: Fall 2024 pset9, “Nonparametric partitioning regression”; restores the detailed occupancy, conditional variance, conditional bias, and IMSE calculation while correcting the legacy bin occupancy growth condition. -->

[Back to the problem map](#problem-map)

## Further practice

<a id="problem-5"></a>

### Problem 5. Basis reparameterization and numerical conditioning

**Chapter connection:** [Section 2.2: Basis normalization](../lectures/week-12.md#22-basis-normalization)

Let $\widetilde p_k(x)=A_kp_k(x)$ for a nonsingular matrix $A_k$.

1. Show that the population approximation spaces generated by $p_k$ and $\widetilde p_k$ are the same.

2. Show that, when both sample design matrices have full column rank, their fitted values are identical.

3. Explain why Gram matrix eigenvalue conditions and numerical stability can nevertheless look very different under the two parameterizations.

4. What normalization would you choose for theoretical analysis?

<!-- New structural problem motivated by the normalized-basis construction in the legacy slides. -->

[Back to the problem map](#problem-map)

<a id="problem-6"></a>

### Problem 6. One regularization principle, two methods

**Chapter connection:** [Section 6: Kernel and series methods compared](../lectures/week-12.md#6-kernel-and-series-methods-compared)

Suppose a one dimensional kernel estimator has squared error order

$$
h^{2p}+\frac1{nh},
$$

and a series estimator has squared error order

$$
k^{-2\alpha}+\frac{k}{n}.
$$

1. Derive both optimal tuning orders and rates.

2. Explain why smaller $h$ and larger $k$ both increase flexibility.

3. Describe how leave one out cross validation can be used for each method.

4. Give two reasons why the tuning parameter selected for prediction may not be suitable for pointwise inference.

<!-- New synthesis problem connecting Weeks 11 and 12. -->

[Back to the problem map](#problem-map)

## Completion check

After the core, you should be able to:

1. distinguish $\mu_k$ from $\mu$;

2. derive the series estimation/approximation decomposition;

3. obtain integrated rates on the correct square root scale;

4. construct pointwise sandwich standard errors;

5. derive the linear smoother leave one out formula; and

6. recognize partitioning as a nonadaptive series estimator; and

7. derive its uniform occupancy, conditional bias, conditional variance, and IMSE expansion.
