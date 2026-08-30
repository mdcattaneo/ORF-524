# ORF 524 Practice Module 3: Decision Theory and Point Estimation

**Status:** Ready for Fall 2026  
**Last updated:** August 30, 2026  
**Primary chapter:** [Week 3: Decision theory and point estimation](../lectures/week-03.md)

## Purpose

This ungraded module compares estimators only after specifying the model, loss, parameter space, and competitor class. Problems 1--4 form the core bank. They also connect complete sufficiency to minimal sufficiency and ancillarity, and give a compact Bayes--minimax calculation. Problems 5--6 develop a general UMVU criterion and the multivariate shrinkage phenomenon.

<a id="problem-map"></a>

## Problem map

| Problem | Bank | Main task |
|---|---|---|
| [1. Risk depends on the parameter space](#problem-1) | Core | Compare risk, dominance, Bayes rules, and minimaxity under different parameter spaces. |
| [2. Rao--Blackwellization in the Bernoulli model](#problem-2) | Core | Construct UMVU estimators using conditioning, sufficiency, and completeness. |
| [3. Unbiasedness versus MSE in the Uniform model](#problem-3) | Core | Compare maximum likelihood, moment, UMVU, and minimum MSE estimators. |
| [4. Information, reparameterization, and a boundary failure](#problem-4) | Core | Apply information bounds and diagnose their failure under moving support. |
| [5. A covariance characterization of UMVU](#problem-5) | Further | Prove a general covariance criterion for UMVU estimation. |
| [6. Unbiased does not imply admissible in dimensions three and higher](#problem-6) | Further | Derive shrinkage risk and reconcile UMVU with inadmissibility. |

## How to work with this module

For each problem, state the parameter space, target, loss, and competitor class before making an optimality claim. Distinguish construction from comparison, identify the exact role of sufficiency and completeness, and verify every regularity condition before using score or information identities.

AI may be used on these ungraded problems under the syllabus policy and the staged protocol in [`AGENTS.md`](AGENTS.md). Begin with a genuine attempt, request the least revealing useful hint, and verify every generated claim line by line. The closed book midterms remain fully unaided.

For released exam practice, use the [question level first half guide](../exams/README.md#question-level-guide-for-the-first-half). Confirm that the listed material has been covered, attempt the selected question before opening its solution, and treat historical exam instructions as superseded by the current syllabus.

## Core practice

<a id="problem-1"></a>

### Problem 1. Risk depends on the parameter space

**Chapter connections:** [Section 1.2: Dominance and admissibility](../lectures/week-03.md#12-dominance-and-admissibility) and [Section 1.3: Bayes and minimax criteria](../lectures/week-03.md#13-bayes-and-minimax-criteria)

Let $X_1,\ldots,X_n$ be i.i.d. $\mathcal N(\theta,1)$ and use squared error loss to estimate $\theta$. Consider $\delta_0=\overline X_n$ and $\delta_c=c\overline X_n$.

1. Derive both risk functions. For fixed $c\in(0,1)$, show that neither estimator dominates the other when $\Theta=\mathbb R$.

2. Now let $\Theta=[-M,M]$ with $M>0$. Show that

$$
c_M=\frac{nM^2}{1+nM^2}
$$

   makes $\delta_{c_M}$ dominate $\delta_0$ on this restricted parameter space.

3. Explain why parts 1 and 2 do not contradict one another and why admissibility must name both a parameter space and a decision rule class.

4. Return to $\Theta=\mathbb R$ and put the prior $\theta\sim\mathcal N(0,v)$, $v>0$, on the parameter. Find the posterior mean and Bayes rule under squared loss, and show that its Bayes risk is $v/(1+nv)$. Use these proper priors as $v\to\infty$ to prove that $\overline X_n$ is minimax, with worst case risk $1/n$.

<!-- Source lineage: Fall 2024 pset3, “Admissibility”; rewritten to make parameter-space dependence and crossing risk functions explicit. -->

[Back to the problem map](#problem-map)

<a id="problem-2"></a>

### Problem 2. Rao--Blackwellization in the Bernoulli model

**Chapter connections:** [Section 3.1: Rao--Blackwellization](../lectures/week-03.md#31-rao--blackwellization), [Section 3.3: Bernoulli and exponential family completeness](../lectures/week-03.md#33-bernoulli-and-exponential-family-completeness), [Section 3.4: Lehmann--Scheffe theorem](../lectures/week-03.md#34-lehmann--scheffe-theorem), and [Section 3.6: Bahadur and Basu](../lectures/week-03.md#36-structural-consequences-bahadur-and-basu)

Let $X_1,\ldots,X_n$ be i.i.d. $\mathsf{Bernoulli}(\theta)$, $\theta\in(0,1)$, and let $T=\sum_iX_i$.

1. Prove that $T$ is sufficient and complete. What does Bahadur's theorem then imply?

2. Start from the unbiased estimator $U_1=X_1$ of $\theta$. Compute $\mathbb E_\theta[U_1\mid T]$ and compare its variance with that of $U_1$.

3. For $n\geq2$, start from the unbiased estimator $U_2=X_1X_2$ of $\theta^2$. Compute $\mathbb E_\theta[U_2\mid T]$.

4. Use Lehmann--Scheffe to identify the UMVU estimators of $\theta$ and $\theta^2$. State exactly where sufficiency and completeness enter. What does Basu's theorem imply about $T$ and any ancillary statistic $A$?

<!-- Source lineage: Bernoulli sufficiency/completeness and legacy “UMVU estimation”; adds direct Rao--Blackwell calculations absent from the recent core sets. -->

[Back to the problem map](#problem-map)

<a id="problem-3"></a>

### Problem 3. Unbiasedness versus MSE in the Uniform model

**Chapter connections:** [Section 2: Bias, variance, and restricted optimality](../lectures/week-03.md#2-bias-variance-and-restricted-optimality), [Sections 3.2--3.4: Completeness and Lehmann--Scheffe](../lectures/week-03.md#32-completeness), and [Section 5.1: Method of moments](../lectures/week-03.md#51-method-of-moments)

Let $X_1,\ldots,X_n$ be i.i.d. $\mathsf{Uniform}(0,\theta)$, $\theta>0$, and $M_n=X_{(n)}$.

1. Show that $M_n$ is complete as well as sufficient. Using

$$
\mathbb E_\theta[M_n]=\frac{n}{n+1}\theta,
$$

   find the UMVU estimator of $\theta$.

2. In the class $\mathcal D=\lbrace kM_n:k>0\rbrace$, derive the risk under squared loss and find the value of $k$ minimizing risk for every $\theta>0$.

3. Use $\mathbb E_\theta[X]=\theta/2$ and $\mathbb E_\theta[X^2]=\theta^2/3$ to derive the method of moments estimators based on the first and second raw moments. Explain why both population moments identify $\theta$ on $(0,\infty)$ but produce different estimators.

4. Compare the bias and MSE of the usual almost sure MLE version $M_n$, the first moment estimator $2\overline X_n$, the UMVU estimator, and the minimum MSE estimator in $\mathcal D$. Recall that the estimator is defined arbitrarily on the null sample where the likelihood has no maximizer in $(0,\infty)$.

5. Explain why the UMVU and minimum MSE answers differ without contradicting either optimality statement.

<!-- Source lineage: Fall 2024 pset3, “UMVU estimation, uniform distribution”; retains the bias--variance comparison, makes the competitor classes explicit, and adds alternative identifying moments. -->

[Back to the problem map](#problem-map)

<a id="problem-4"></a>

### Problem 4. Information, reparameterization, and a boundary failure

**Chapter connections:** [Section 4.1: Regular likelihood identities](../lectures/week-03.md#41-regular-likelihood-identities), [Sections 4.2--4.4: Cramer--Rao bounds and UMVU construction](../lectures/week-03.md#42-cramer--rao-lower-bound), and [Section 5.2: Maximum likelihood](../lectures/week-03.md#52-maximum-likelihood)

Let $X_1,\ldots,X_n$ be i.i.d. Exponential with mean parameter $\theta>0$:

$$
f(x;\theta)=\frac1\theta e^{-x/\theta}\mathbf 1\lbrace x>0\rbrace.
$$

1. Derive the score and Fisher information for $\theta$. Verify the mean zero score identity.

2. Derive the Cramer--Rao lower bound for unbiased estimators of $\theta$ and show that $\overline X_n$ attains it.

3. Reparameterize by $\eta=\log\theta$. Derive the information for $\eta$ and verify that the Cramer--Rao bound for the target $e^\eta=\theta$ is unchanged.

4. The invariance property of maximum likelihood gives the MLE $\widehat\eta=\log\overline X_n$. Use strict Jensen inequality to show that it is not unbiased for $\eta$.

5. Contrast this calculation with the $\mathsf{Uniform}(0,\theta)$ likelihood. Show that differentiating only its smooth part produces a “score” with nonzero expectation, and identify the violated regularity condition.

6. Suppose more generally that a two parameter regular model has information matrix

$$
\mathcal I_n(\theta)=
\begin{pmatrix}
a & b\\
b & c
\end{pmatrix},
\qquad a>0,\quad c>0,\quad ac-b^2>0.
$$

   Find the Cramer--Rao bound for unbiased estimation of the first coordinate when the second is a nuisance parameter. Compare it with $1/a$, and state when the two bounds agree.

<!-- Source lineage: Fall 2024 pset3, “CRLB, transformations,” combined with the Week 2 Uniform boundary example and a regular Exponential calculation. -->

[Back to the problem map](#problem-map)

## Further practice

<a id="problem-5"></a>

### Problem 5. A covariance characterization of UMVU

**Chapter connection:** [Section 2: Bias, variance, and restricted optimality](../lectures/week-03.md#2-bias-variance-and-restricted-optimality)

Let $\mathcal T_u$ be the class of unbiased estimators of a scalar target $\tau(\theta)$ with finite second moments. Fix $T\in\mathcal T_u$. Prove that

$$
T\text{ is UMVU}
\quad\Longleftrightarrow\quad
\mathrm{Cov}_\theta(T,U-T)=0
$$

for every $U\in\mathcal T_u$ and every $\theta$. Use the unbiased path $T_c=T+c(U-T)$ and handle separately the case $\mathbb V_\theta(U-T)=0$.

<!-- Source lineage: Fall 2024 pset3, “UMVU estimation.” -->

[Back to the problem map](#problem-map)

<a id="problem-6"></a>

### Problem 6. Unbiased does not imply admissible in dimensions three and higher

**Chapter connections:** [Section 1.2: Dominance and admissibility](../lectures/week-03.md#12-dominance-and-admissibility), [Section 2: Bias, variance, and restricted optimality](../lectures/week-03.md#2-bias-variance-and-restricted-optimality), and [Section 3.3: Full rank exponential family completeness](../lectures/week-03.md#33-bernoulli-and-exponential-family-completeness)

Let $X\sim\mathcal N_m(\mu,I_m)$, $m\geq3$, and estimate $\mu$ under loss $\Vert a-\mu\Vert^2$. For $a>0$, define

$$
\delta_a(X)=\left(1-\frac{a}{\Vert X\Vert^2}\right)X.
$$

You may use Stein's identity

$$
\mathbb E_\mu[(X_j-\mu_j)g(X)]
=\mathbb E_\mu[\partial_jg(X)].
$$

1. Show that

$$
R(\mu,\delta_a)-R(\mu,X)
=\lbrace a^2-2a(m-2)\rbrace\mathbb E_\mu[\Vert X\Vert^{-2}].
$$

2. Find the range of $a$ giving strict dominance and the best constant $a$ within this family.

3. Use the full rank exponential family result to explain why $X$ is complete and sufficient for $\mu$ and why each $X_j$ is UMVU for $\mu_j$. Then explain why this componentwise unbiased optimality does not prevent inadmissibility under the joint vector loss.

<!-- Source lineage: Fall 2024 pset3, optional “Shrinkage Estimator”; coefficient generalized to expose the full risk calculation and complete-sufficiency connection. -->

[Back to the problem map](#problem-map)

## Completion check

After the core, you should be able to explain:

1. why risk comparisons are functions rather than single numbers;

2. how Rao--Blackwell and completeness play different roles;

3. how complete sufficiency yields minimality and independence from ancillaries;

4. how Bayes risks can certify a minimax rule;

5. why minimum variance among unbiased estimators is not minimum MSE among all estimators; and

6. why score and information identities can fail with moving support;

7. why identification of a full model does not guarantee identification by a selected moment condition.
