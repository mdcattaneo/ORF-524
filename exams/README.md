# Past ORF 524 exams

This directory contains exams from previous offerings of ORF 524 and accompanying instructor
solutions. They are ungraded study resources for the redesigned Fall 2026 course.

Historical exam instructions—including duration, permitted notes, submission method, point totals,
and topic sequence—do not apply to Fall 2026. The current syllabus and instructor announcements
govern the current course and its closed book, unaided midterms.

- The PDF files at this directory's root are previously administered exam questions.
- [`solutions/`](solutions/): instructor solutions for study after a genuine attempt.
- [`AGENTS.md`](AGENTS.md): protocol for practicing with these materials and an AI assistant.

## Recommended workflow

1. Choose an exam using the alignment table below and confirm that the relevant weeks have been
   covered.
2. Attempt selected problems without looking at the corresponding solution. For a realistic exam
   simulation, work without AI, notes, or outside assistance.
3. After the attempt, use AI to diagnose the first gap, request staged hints, or audit the argument.
4. Consult the solution only after completing an attempt or deliberately switching to
   solution study mode.
5. Close the solution and reconstruct the argument unaided, then solve a nearby variation.

## Question level guide for the first half

The following entries are representative question level matches for Weeks 1--6. Use only the listed parts when a question also contains later material. These links support cumulative practice; they do not promise that the current midterms will reproduce a past format, question, or difficulty level.

| Current week | Representative released questions | Primary connection |
|---|---|---|
| Week 1 | [2020 Final, Question 2](ORF524-Fall2020-Final.pdf); [2023 Midterm, Question 3](ORF524-Fall2023-Midterm.pdf); [2024 Final, Question 1](ORF524-Fall2024-Final.pdf) | Absolute loss prediction, conditional expectation, MSE decomposition, and identification of predictors |
| Week 2 | [2022 Midterm, Question 1, parts 1--3](ORF524-Fall2022-Midterm.pdf); [2023 Midterm, Question 2](ORF524-Fall2023-Midterm.pdf) | Sufficiency under transformations and the likelihood ratio characterization of minimal sufficiency |
| Week 3 | [2022 Midterm, Question 1, parts 4--5, and Question 2, part 2](ORF524-Fall2022-Midterm.pdf); [2023 Midterm, Question 1, parts 1--6](ORF524-Fall2023-Midterm.pdf); [2024 Midterm, Question 1](ORF524-Fall2024-Midterm.pdf) | Completeness, UMVU estimation, risk, Rao--Blackwell improvement, method of moments, likelihood, and efficiency |
| Week 4 | [2022 Midterm, Question 2, part 3](ORF524-Fall2022-Midterm.pdf); [2023 Midterm, Question 1, part 7](ORF524-Fall2023-Midterm.pdf); [2025 Midterm, Question 1](ORF524-Fall2025-Midterm.pdf) | Test calibration, power, Neyman--Pearson, UMP testing, consistency, and confidence set inversion |
| Week 5 | [2019 Final, Question 1, parts 1 and 3](ORF524-Fall2019-Final.pdf); [2020 Final, Question 1, parts 3--4](ORF524-Fall2020-Final.pdf) | Laws of large numbers, stochastic control, multivariate CLTs, and the joint limit of the sample mean and variance |
| Week 6 | [2020 Final, Question 1, part 4, and Question 3, parts 3--8](ORF524-Fall2020-Final.pdf); [2021 Final, Question 2](ORF524-Fall2021-Final.pdf); [2024 Midterm, Question 2, parts 2--10](ORF524-Fall2024-Midterm.pdf) | Delta methods, Studentization, feasible intervals and tests, boundary limits, and pointwise versus uniform validity |

## Question level guide for Week 9

The following are the cleanest released matches after completing Week 9. Use only the listed parts when a question also develops later applications. Several released solutions use legacy $H/\Sigma$ notation or invoke a Taylor expansion; translate those arguments into the current $\Psi_\theta(\theta_0)/B$ notation and the coordinatewise mean value theorem supported by a local uniform LLN.

| Current week | Representative released questions | Primary connection |
|---|---|---|
| Week 9 | [2023 Final, Question 3, parts 1--2 and 5--7](ORF524-Fall2023-Final.pdf); [2021 Final, Question 3, parts 3--4](ORF524-Fall2021-Final.pdf); [2022 Final, Question 2, parts 2--5 and 7](ORF524-Fall2022-Final.pdf); [2025 Midterm, Question 2, parts 3--6](ORF524-Fall2025-Midterm.pdf) | Argmin consistency and uniform convergence; one step estimation; Z-estimator asymptotic linearity, sandwich variance, and feasible inference; misspecified likelihood influence functions and robust standard errors |

## Question level guide for Week 10

The following are the cleanest released matches after completing Week 10. The 2023 and 2025 Poisson questions use earlier GLS/FGLS or weighting language: interpret an efficient infeasible weight as fixed at the population truth, construct feasible weights from a preliminary estimator and freeze them in the second step, and do not silently differentiate a parameter dependent residual criterion as though it were a fixed weight criterion.

| Current week | Representative released questions | Primary connection |
|---|---|---|
| Week 10 | [2025 Final, Question 1](ORF524-Fall2025-Final.pdf) and [Question 2, parts 1--6](ORF524-Fall2025-Final.pdf); [2023 Final, Question 1, parts 1--5](ORF524-Fall2023-Final.pdf); [2022 Final, Question 2, part 10](ORF524-Fall2022-Final.pdf); [2024 Final, Question 3, part 2(c)--(e)](ORF524-Fall2024-Final.pdf) | Best linear approximation and OLS; regression score, Hessian, sandwich variance, and robust standard errors; infeasible and feasible variance weighting; comparison of least squares, weighted least squares, and likelihood efficiency; finite dimensional two step influence adjustment and inference |

The released solution to 2025 Final Question 1, part 2 contains a minor norm bound typo: the relevant bound is $\beta'\Sigma_X\beta\leq\Vert\beta\Vert^2\Vert\Sigma_X\Vert_{\mathrm{op}}$. In 2025 Final Question 2, part 6, read the efficiency answer as using the infeasible fixed weight $1/\mu_{\beta_0}(X)$, or its estimated and frozen feasible version. These corrections do not change the questions themselves.

For 2024 Final Question 3, part 2(c)--(e), take the identification and consistency setup in parts 2(a)--(b) as given and focus on the parametric first step adjustment, variance estimation, and inference. In part 2(c), the printed criterion is missing the observation indicator: use $n^{-1}\sum_{i=1}^n t_i\rho(y_i,x_i;\theta)/\widehat p_n(x_i)$, as required by part 2(b) and used in the released derivation. In part 2(e), center the interval at $\widehat\theta_n$. The full missing data interpretation is deliberately deferred to Week 14.

## Question level guide for Weeks 11, 12, and 14

The released archive has strong kernel and semiparametric calculations but no clean series regression problem aligned with Week 12. Use the current Week 12 assignment for that topic rather than forcing a misleading match to an older exam.

| Current week | Representative released questions | Primary connection |
|---|---|---|
| Week 11 | [2023 Final, Question 2](ORF524-Fall2023-Final.pdf); [2025 Final, Question 3](ORF524-Fall2025-Final.pdf) | Kernel density bias, variance, bandwidth choice, dependence conditions, and uniform consistency |
| Week 12 | No clean released match | Series projection, approximation error, growing dimension, and pointwise series inference are practiced in the current assignment module |
| Week 14 | [2021 Final, Question 1](ORF524-Fall2021-Final.pdf); [2022 Final, Question 1](ORF524-Fall2022-Final.pdf); [2024 Final, Question 3, part 2](ORF524-Fall2024-Final.pdf) | U-statistic and V-statistic decompositions, integrated squared density, semiparametric bandwidth conditions, and first step adjustment under missing outcomes |

## Alignment with the Fall 2026 course

The labels below identify the principal connections, not every prerequisite or subpart. Some exams
contain topics outside the redesigned core; the AI protocol requires those parts to be identified
before a student begins.

| Offering | Exam | Principal Fall 2026 alignment | Solution |
|---|---|---|---|
| 2019 | [Midterm](ORF524-Fall2019-Midterm.pdf) | Weeks 1--5: probability and finite sample inference | [Solution](solutions/ORF524-Fall2019-Midterm-solutions.pdf) |
| 2019 | [Final](ORF524-Fall2019-Final.pdf) | Weeks 1--6 and 9: asymptotics, point estimation, inference, misspecified likelihood | [Solution](solutions/ORF524-Fall2019-Final-solutions.pdf) |
| 2020 | [Midterm](ORF524-Fall2020-Midterm.pdf) | Weeks 1--5: probability, decision theory, and finite sample inference | [Solution](solutions/ORF524-Fall2020-Midterm-solutions.pdf) |
| 2020 | [Final](ORF524-Fall2020-Final.pdf) | Weeks 1--6 and 9: asymptotics, loss based estimation, and inference | [Solution](solutions/ORF524-Fall2020-Final-solutions.pdf) |
| 2021 | [Midterm](ORF524-Fall2021-Midterm.pdf) | Weeks 1--5: concentration, decision theory, and finite sample inference | [Solution](solutions/ORF524-Fall2021-Midterm-solutions.pdf) |
| 2021 | [Final](ORF524-Fall2021-Final.pdf) | Weeks 5--6, 9, and 11: large sample and uniform inference; also U- and V-statistics | [Solution](solutions/ORF524-Fall2021-Final-solutions.pdf) |
| 2022 | [Midterm](ORF524-Fall2022-Midterm.pdf) | Weeks 2--4: sufficiency, completeness, and finite sample inference | [Solution](solutions/ORF524-Fall2022-Midterm-solutions.pdf) |
| 2022 | [Final](ORF524-Fall2022-Final.pdf) | Weeks 9--11 and 14: conditional Z-estimation, optimal weighting, and semiparametric bandwidth selection | [Solution](solutions/ORF524-Fall2022-Final-solutions.pdf) |
| 2023 | [Midterm](ORF524-Fall2023-Midterm.pdf) | Weeks 1--4: prediction, identification, sufficiency, estimation, and testing | [Solution](solutions/ORF524-Fall2023-Midterm-solutions.pdf) |
| 2023 | [Final](ORF524-Fall2023-Final.pdf) | Weeks 9--11: Poisson regression, feasible GLS, M-estimation, and kernel density estimation | [Solution](solutions/ORF524-Fall2023-Final-solutions.pdf) |
| 2024 | [Midterm](ORF524-Fall2024-Midterm.pdf) | Weeks 2--3 and 6: Rao--Blackwell improvement, likelihood estimation, delta methods, and feasible asymptotic inference | [Solution](solutions/ORF524-Fall2024-Midterm-solutions.pdf) |
| 2024 | [Final](ORF524-Fall2024-Final.pdf) | Weeks 1, 4, 9--10, and 14: prediction, conservative inference, and M-estimation with missing data | [Solution](solutions/ORF524-Fall2024-Final-solutions.pdf) |
| 2025 | [Midterm](ORF524-Fall2025-Midterm.pdf) | Weeks 2, 4, and 9: uniformly most powerful testing and misspecified likelihood | [Solution](solutions/ORF524-Fall2025-Midterm-solutions.pdf) |
| 2025 | [Final](ORF524-Fall2025-Final.pdf) | Weeks 1 and 9--11: linear prediction, quasi-likelihood, weighted nonlinear least squares, and uniform kernel density consistency | [Solution](solutions/ORF524-Fall2025-Final-solutions.pdf) |
