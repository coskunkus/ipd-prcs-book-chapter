# Statistical Inference for the Inverse Pareto Distribution under Progressive Type-II Random Censoring: Applications to Sustainable Reliability

R code for the book chapter by Indrajeet Kumar, Kapil Kumar, Ashish Raj and Coskun Kus.

## What the chapter does

The inverse Pareto distribution (IPD) is studied as a lifetime model under progressive Type-II random censoring, a design in which some surviving units are withdrawn according to a prefixed plan while others leave the experiment at random. The random censoring is described by the Koziol-Green model, which links the censoring distribution to the lifetime distribution through a proportional-hazards structure. Under that model the censoring parameter can be profiled out of the likelihood in closed form, so that maximum likelihood estimation reduces to a one-dimensional problem, and the full conditional posterior of the censoring parameter is a gamma distribution, so that the Bayesian computation needs a Metropolis-Hastings step for the shape parameter only.

## Files

| File | What it does |
| --- | --- |
| `code/core_ipd_prcs.R` | Shared core. Distribution functions, entropies, order-statistic means, the Koziol-Green mechanism, the PRCS sampler, the likelihood, the score, the observed information matrix, the estimation routines, the MCMC sampler, the HPD intervals, the figure device and the LaTeX table emitter. Every other script sources this file. |
| `code/check_formulas.R` | Verification of every closed form against numerical integration, finite differences, grid search and quadrature of the posterior, followed by Monte Carlo checks reported as z-scores. |
| `code/run_all.R` | Everything in the chapter except the Monte Carlo study: the figure, the engineering application, the parametric bootstrap of the Kolmogorov-Smirnov statistic and the LaTeX of the three tables. Under a minute. |
| `code/run_simulation.R` | The Monte Carlo study, in parallel and in parts. Two chains are run per replication, one under the informative prior of the chapter and one under a nearly flat one, and both interval definitions are computed from each, so the report carries six variants and a question about the choice of prior or interval never requires the study to be run again. This is the only long computation. |
| `code/publish_github.R` | This publisher. |

## Run order

```r
setwd("<the folder that contains code/>")
source("code/check_formulas.R")    # verification of every closed form
source("code/run_all.R")           # figure, application, tables
source("code/run_simulation.R")    # the long one, run in parts
```

The Monte Carlo study is sixty runs and takes about three and a half hours, which is longer than one sitting. It is therefore carried out in parts: set `PART` at the top of `code/run_simulation.R` to 1, run it, set it to 2 at the next opportunity, and so on to 4. Each run writes the result of a configuration to `results/sim_parts/` as soon as it finishes and skips whatever is already there, so an interruption costs at most the configuration in progress. Setting `PART` to 0 picks up everything that is still missing. When the sixty are complete the script assembles the report and the LaTeX of the four tables.

The scripts locate the project folder themselves, so they also run from the RStudio Source button and from `Rscript` without a `setwd()` call. They create `results/` and `figures/` if these are absent and write nothing anywhere else. Reports are plain text and figures are written both as `.eps` and as `.pdf`.

## Requirements

Base R only; no contributed package is needed. `parallel` and `stats` are part of the standard distribution.

## Reproducibility

Everything in the chapter except the Monte Carlo study is reproducible digit for digit. The application imposes a censoring design on a complete data set, so its numbers depend on a stream of pseudo-random numbers. That stream is not taken from the generator of R but from one written out in the core, the multiplicative congruential recursion `x <- 16807 x mod (2^31 - 1)` started from a seed that is stated in the report. The censoring times, the choice of the units withdrawn and the whole MCMC draw from it, in the same order, so `run_all.R` returns the same digits on any machine and under any version of R. The full conditional of the censoring parameter is a gamma with an integer shape, which is drawn as a sum of exponentials, and the normal proposal is drawn by the Box-Muller transformation, so that no step falls back on a library generator.

The Monte Carlo study is deliberately outside that guarantee, since two Monte Carlo runs of the same configuration agree within the Monte Carlo standard error rather than to the last digit. It uses the master seed `20260908`, set in `code/core_ipd_prcs.R`, and gives every replication a L'Ecuyer stream of its own, derived from that seed and from the position of the configuration in a fixed list. Replication j of configuration i is therefore the same experiment whatever part it was run in, however many workers were used, and on whatever machine.

Every iterative fit stops only when two conditions hold at once: the absolute change in the log-likelihood is below `1e-12` and the largest relative change in the parameters is below `1e-9`. The second condition is not optional, because under heavy censoring the log-likelihood is flat along a ridge and a rule based on the log-likelihood alone stops while the estimates are still moving. After convergence the score and the eigenvalues of the observed information matrix are checked and reported. Where a profile solution exists it is solved directly and the iteration is used as an independent check, with the largest relative difference between the two routes reported.

## Notation shared with the chapter

| Symbol | Meaning |
| --- | --- |
| `alpha` | shape parameter of the IPD |
| `beta` | Koziol-Green censoring parameter |
| `n` | number of units placed on test |
| `m` | number of observed event times |
| `k` | number of observed failures |
| `W` | progressive withdrawal scheme, `sum(W) + m = n` |
| `y`, `d` | observed event times and event indicators |
| `z_i` | `y_i / (1 + y_i)` |
| `J` | observed information matrix |
| `zeta`, `zeta0` | MCMC iterations and burn-in |

## License

MIT. See `LICENSE`.
