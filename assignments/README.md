# ORF 524 Practice Assignments

**Status:** Complete practice sequence  
**Last updated:** October 2, 2026

## Purpose

The assignments are an ungraded practice system for ORF 524. They are built from the course's legacy problem sets, precepts, and selected released material, but reorganized around the redesigned weekly chapters.

The organizing principle is **topic based files with weekly pacing**: each file has a durable mathematical topic, and the course identifies which part of that file is recommended in a given week. For Weeks 1--6, the topic and week coincide naturally.

These modules are not discussion guides alone. They require proofs, likelihood and moment calculations, exact distribution theory, asymptotic expansions, variance calculations, standard errors, tests, and confidence intervals. Their concise problem maps describe the mathematical destination, not the amount of technical work required to reach it. The two closed book, unaided midterms require students to carry out such calculations independently, so serious written practice is an essential use of this bank.

## First half modules

| File | Primary chapter | Topic | Solutions |
|---|---|---|---|
| [`topic-01-probability-prediction.md`](topic-01-probability-prediction.md) | Week 1 | Probability, expectation, loss, and prediction | [Solutions](solutions/topic-01-solutions.md) |
| [`topic-02-models-likelihood-sufficiency.md`](topic-02-models-likelihood-sufficiency.md) | Week 2 | Models, identification, likelihood, and sufficiency | [Solutions](solutions/topic-02-solutions.md) |
| [`topic-03-decision-point-estimation.md`](topic-03-decision-point-estimation.md) | Week 3 | Decision theory and point estimation | [Solutions](solutions/topic-03-solutions.md) |
| [`topic-04-testing-confidence-sets.md`](topic-04-testing-confidence-sets.md) | Week 4 | Testing and confidence sets | [Solutions](solutions/topic-04-solutions.md) |
| [`topic-05-convergence-limit-theorems.md`](topic-05-convergence-limit-theorems.md) | Week 5 | Convergence and limit theorems | [Solutions](solutions/topic-05-solutions.md) |
| [`topic-06-delta-asymptotic-inference.md`](topic-06-delta-asymptotic-inference.md) | Week 6 | Delta method and asymptotic inference | [Solutions](solutions/topic-06-solutions.md) |

Modules 1--6 form the complete first half practice package. Every module has four core problems and a further bank whose size reflects the useful historical material for that topic. Recommended selections and pacing may be refined as the semester unfolds.

## Second half modules

Weeks 7 and 8 are Midterm 1 and fall break, respectively, so the second half numbering begins at Week 9. Week 13 has no matching module.

| File | Primary chapter | Topic | Solutions |
|---|---|---|---|
| [`topic-09-m-z-estimation.md`](topic-09-m-z-estimation.md) | Week 9 | M and Z estimation | [Solutions](solutions/topic-09-solutions.md) |
| [`topic-10-regression-two-step-estimation.md`](topic-10-regression-two-step-estimation.md) | Week 10 | Regression applications and two step estimation | [Solutions](solutions/topic-10-solutions.md) |
| [`topic-11-kernel-nonparametric-methods.md`](topic-11-kernel-nonparametric-methods.md) | Week 11 | Kernel based nonparametric methods | [Solutions](solutions/topic-11-solutions.md) |
| [`topic-12-series-nonparametric-methods.md`](topic-12-series-nonparametric-methods.md) | Week 12 | Series based nonparametric methods | [Solutions](solutions/topic-12-solutions.md) |
| [`topic-14-semiparametric-methods.md`](topic-14-semiparametric-methods.md) | Week 14 | Semiparametric methods | [Solutions](solutions/topic-14-solutions.md) |

## Structure within a module

Each module contains:

1. **Problem map:** a table near the top that summarizes every problem, distinguishes the core and further banks, and provides links to each problem and return links back to the map;

2. **Core practice:** four substantial problems that exercise the week's essential learning goals;

3. **Further practice:** optional problems that provide a boundary case, a methodological extension, a cumulative calculation laboratory, or additional technical depth; and

4. **Completion check:** a short list students can use to diagnose whether they can explain and execute the module's durable ideas without assistance.

Problem length matters more than problem count. A problem with several derivations may constitute most of the core practice by itself. The four core problems form a coverage bank, not an expectation that every student complete every part in one week. Before each module is used, the instructor or preceptor should identify a smaller recommended route, normally selected parts of two or three core problems, based on lecture progress. The remaining core and further parts stay available for independent, cumulative, and exam preparation practice.

## Hints and solutions

Solutions for all 11 modules are linked in the tables above and at the top of each practice module. The problem statements and solutions are kept in separate files so that students can attempt the problems before opening the corresponding solutions.

Begin with a genuine attempt. If needed, seek a strategic hint and then a more explicit intermediate hint before studying the full solution. Afterward, close the solution and reconstruct the argument, checking the assumptions and explaining why each step is valid.

AI may be used for these ungraded practice assignments under the syllabus policy. The [`AGENTS.md`](AGENTS.md) file gives AI systems a staged tutoring protocol: guided attempts and progressively stronger hints remain the default, while solution study is a deliberate choice. Released solutions are available for comparison and study; the learning protocol does not grant permission to use AI on an assessment.

## Use in precept

The weekly recommendation should identify a small number of problems for preparation. Precept can
then select parts that expose common misconceptions, compare solution strategies, or repair a
missing prerequisite. Precept should not attempt to work through the entire module.

## Provenance and mathematical review

Every adapted problem must carry a Markdown comment recording its legacy source file and original
problem title or number. Repeated versions across years count as one problem lineage. The new file
should use the notation and assumptions of the current chapter rather than mechanically preserving
legacy wording.

During maintenance and before any weekly release, every changed problem, hint, and solution must be
checked independently. In particular:

- check support and dominating measures;
- state integrability and differentiability assumptions;
- distinguish exact and asymptotic conclusions;
- use $\mathbb{P}$, $\to_{\mathbb{P}}$, $O_{\mathbb{P}}$, and $o_{\mathbb{P}}$ consistently; and
- verify that each part depends only on material already introduced or explicitly marked as an
  extension.

## Exam reserve policy during construction

Past exam questions are a separate source bank and are not copied verbatim into the routine assignment modules. The exam index provides a [question level first half guide](../exams/README.md#question-level-guide-for-the-first-half) for Weeks 1--6, guides for [Week 9](../exams/README.md#question-level-guide-for-week-9) and [Week 10](../exams/README.md#question-level-guide-for-week-10), and a combined guide for [Weeks 11, 12, and 14](../exams/README.md#question-level-guide-for-weeks-11-12-and-14), while keeping exam practice separate from routine weekly work. Past questions may later be:

- retained as protected models for constructing new midterms;
- released as clearly identified practice exam questions; or
- substantially adapted into synthesis problems after the instructor approves their release.

This separation prevents routine practice materials from accidentally exhausting the best
assessment questions.
