# =============================================================================
# publish_github.R
#
# Publishes the code of this chapter to GitHub.
#
# What is published: the scripts in code/, together with README.md, LICENSE,
# .gitignore and .gitattributes. The manuscript, the reports, the figures and
# the reference PDFs stay in the working folder. Anyone who clones the
# repository reproduces them by running the scripts.
#
# Two details that otherwise cost time.
#
#   1. The project path contains non-ASCII characters. Git is an external
#      program that receives its arguments in the system encoding, so a
#      Turkish folder name reaches it corrupted and every command fails with
#      "cannot change to ...". No path is therefore ever passed to git.
#      Everything is staged into an ASCII-only folder under TEMP and the
#      working directory is changed into it before each git call.
#
#   2. The repository is rebuilt from the current scripts on every run, so a
#      fresh history is created each time and the push has to replace what is
#      on GitHub. force_push is therefore TRUE. That is correct for a
#      repository that holds the code of one paper and nothing else. If a
#      push is ever refused, inspect the remote before overwriting it.
# =============================================================================

rm(list = ls())

local({
  f <- "code/core_ipd_prcs.R"
  if (!file.exists(f)) f <- "core_ipd_prcs.R"
  source(f)
})

PROJECT <- setup_project()


# -----------------------------------------------------------------------------
# Settings
# -----------------------------------------------------------------------------

GH_USER <- "coskunkus"
GH_REPO <- "ipd-prcs-book-chapter"

AUTHOR_NAME  <- "Coskun Kus"
AUTHOR_EMAIL <- "coskunkus@gmail.com"

force_push <- TRUE

BRANCH <- "main"

PAPER_TITLE <- paste0("Statistical Inference for the Inverse Pareto ",
                      "Distribution under Progressive Type-II Random ",
                      "Censoring: Applications to Sustainable Reliability")

PAPER_AUTHORS <- "Indrajeet Kumar, Kapil Kumar, Ashish Raj and Coskun Kus"

FILES <- c("core_ipd_prcs.R", "run_all.R", "run_simulation.R",
           "check_formulas.R", "publish_github.R")


# -----------------------------------------------------------------------------
# The staging folder, ASCII only
# -----------------------------------------------------------------------------

stage <- file.path(tempdir(), paste0("gh_", GH_REPO))

unlink(stage, recursive = TRUE, force = TRUE)
dir.create(file.path(stage, "code"), recursive = TRUE, showWarnings = FALSE)

for (f in FILES)
{
  src <- file.path(PROJECT, "code", f)

  if (!file.exists(src))
    stop("Missing file: code/", f)

  file.copy(src, file.path(stage, "code", f), overwrite = TRUE)
}

cat("Staged", length(FILES), "files in", stage, "\n")


# -----------------------------------------------------------------------------
# README.md
# -----------------------------------------------------------------------------

readme <- c(
  paste0("# ", PAPER_TITLE),
  "",
  paste0("R code for the book chapter by ", PAPER_AUTHORS, "."),
  "",
  "## What the chapter does",
  "",
  paste0("The inverse Pareto distribution (IPD) is studied as a lifetime ",
         "model under progressive Type-II random censoring, a design in ",
         "which some surviving units are withdrawn according to a prefixed ",
         "plan while others leave the experiment at random. The random ",
         "censoring is described by the Koziol-Green model, which links the ",
         "censoring distribution to the lifetime distribution through a ",
         "proportional-hazards structure. Under that model the censoring ",
         "parameter can be profiled out of the likelihood in closed form, so ",
         "that maximum likelihood estimation reduces to a one-dimensional ",
         "problem, and the full conditional posterior of the censoring ",
         "parameter is a gamma distribution, so that the Bayesian ",
         "computation needs a Metropolis-Hastings step for the shape ",
         "parameter only."),
  "",
  "## Files",
  "",
  "| File | What it does |",
  "| --- | --- |",
  paste0("| `code/core_ipd_prcs.R` | Shared core. Distribution functions, ",
         "entropies, order-statistic means, the Koziol-Green mechanism, the ",
         "PRCS sampler, the likelihood, the score, the observed information ",
         "matrix, the estimation routines, the MCMC sampler, the HPD ",
         "intervals, the figure device and the LaTeX table emitter. Every ",
         "other script sources this file. |"),
  paste0("| `code/check_formulas.R` | Verification of every closed form ",
         "against numerical integration, finite differences, grid search and ",
         "quadrature of the posterior, followed by Monte Carlo checks ",
         "reported as z-scores. |"),
  paste0("| `code/run_all.R` | Everything in the chapter except the Monte ",
         "Carlo study: the figure, the engineering application, the ",
         "parametric bootstrap of the Kolmogorov-Smirnov statistic and the ",
         "LaTeX of the three tables. Under a minute. |"),
  paste0("| `code/run_simulation.R` | The Monte Carlo study, in parallel and ",
         "in parts. Two chains are run per replication, one under the ",
         "informative prior of the chapter and one under a nearly flat one, ",
         "and both interval definitions are computed from each, so the ",
         "report carries six variants and a question about the choice of ",
         "prior or interval never requires the study to be run again. This ",
         "is the only long computation. |"),
  paste0("| `code/publish_github.R` | This publisher. |"),
  "",
  "## Run order",
  "",
  "```r",
  "setwd(\"<the folder that contains code/>\")",
  "source(\"code/check_formulas.R\")    # verification of every closed form",
  "source(\"code/run_all.R\")           # figure, application, tables",
  "source(\"code/run_simulation.R\")    # the long one, run in parts",
  "```",
  "",
  paste0("The Monte Carlo study is sixty runs and takes about three and a ",
         "half hours, which is longer than one sitting. It is therefore ",
         "carried out in parts: set `PART` at the top of ",
         "`code/run_simulation.R` to 1, run it, set it to 2 at the next ",
         "opportunity, and so on to 4. Each run writes the result of a ",
         "configuration to `results/sim_parts/` as soon as it finishes and ",
         "skips whatever is already there, so an interruption costs at most ",
         "the configuration in progress. Setting `PART` to 0 picks up ",
         "everything that is still missing. When the sixty are complete the ",
         "script assembles the report and the LaTeX of the four tables."),
  "",
  paste0("The scripts locate the project folder themselves, so they also ",
         "run from the RStudio Source button and from `Rscript` without a ",
         "`setwd()` call. They create `results/` and `figures/` if these are ",
         "absent and write nothing anywhere else. Reports are plain text and ",
         "figures are written both as `.eps` and as `.pdf`."),
  "",
  "## Requirements",
  "",
  paste0("Base R only; no contributed package is needed. `parallel` and ",
         "`stats` are part of the standard distribution."),
  "",
  "## Reproducibility",
  "",
  paste0("Everything in the chapter except the Monte Carlo study is ",
         "reproducible digit for digit. The application imposes a censoring ",
         "design on a complete data set, so its numbers depend on a stream ",
         "of pseudo-random numbers. That stream is not taken from the ",
         "generator of R but from one written out in the core, the ",
         "multiplicative congruential recursion ",
         "`x <- 16807 x mod (2^31 - 1)` started from a seed that is stated ",
         "in the report. The censoring times, the choice of the units ",
         "withdrawn and the whole MCMC draw from it, in the same order, so ",
         "`run_all.R` returns the same digits on any machine and under any ",
         "version of R. The full conditional of the censoring parameter is a ",
         "gamma with an integer shape, which is drawn as a sum of ",
         "exponentials, and the normal proposal is drawn by the Box-Muller ",
         "transformation, so that no step falls back on a library ",
         "generator."),
  "",
  paste0("The Monte Carlo study is deliberately outside that guarantee, ",
         "since two Monte Carlo runs of the same configuration agree within ",
         "the Monte Carlo standard error rather than to the last digit. It ",
         "uses the master seed `", MASTER_SEED, "`, set in ",
         "`code/core_ipd_prcs.R`, and gives every replication a L'Ecuyer ",
         "stream of its own, derived from that seed and from the position of ",
         "the configuration in a fixed list. Replication j of configuration ",
         "i is therefore the same experiment whatever part it was run in, ",
         "however many workers were used, and on whatever machine."),
  "",
  paste0("Every iterative fit stops only when two conditions hold at once: ",
         "the absolute change in the log-likelihood is below `1e-12` and the ",
         "largest relative change in the parameters is below `1e-9`. The ",
         "second condition is not optional, because under heavy censoring ",
         "the log-likelihood is flat along a ridge and a rule based on the ",
         "log-likelihood alone stops while the estimates are still moving. ",
         "After convergence the score and the eigenvalues of the observed ",
         "information matrix are checked and reported. Where a profile ",
         "solution exists it is solved directly and the iteration is used as ",
         "an independent check, with the largest relative difference between ",
         "the two routes reported."),
  "",
  "## Notation shared with the chapter",
  "",
  "| Symbol | Meaning |",
  "| --- | --- |",
  "| `alpha` | shape parameter of the IPD |",
  "| `beta` | Koziol-Green censoring parameter |",
  "| `n` | number of units placed on test |",
  "| `m` | number of observed event times |",
  "| `k` | number of observed failures |",
  "| `W` | progressive withdrawal scheme, `sum(W) + m = n` |",
  "| `y`, `d` | observed event times and event indicators |",
  "| `z_i` | `y_i / (1 + y_i)` |",
  "| `J` | observed information matrix |",
  "| `zeta`, `zeta0` | MCMC iterations and burn-in |",
  "",
  "## License",
  "",
  "MIT. See `LICENSE`."
)

writeLines(readme, file.path(stage, "README.md"), useBytes = TRUE)


# -----------------------------------------------------------------------------
# LICENSE, .gitignore and .gitattributes
# -----------------------------------------------------------------------------

license <- c(
  "MIT License",
  "",
  paste0("Copyright (c) ", format(Sys.Date(), "%Y"), " ", PAPER_AUTHORS),
  "",
  paste0("Permission is hereby granted, free of charge, to any person ",
         "obtaining a copy of this software and associated documentation ",
         "files (the \"Software\"), to deal in the Software without ",
         "restriction, including without limitation the rights to use, copy, ",
         "modify, merge, publish, distribute, sublicense, and/or sell copies ",
         "of the Software, and to permit persons to whom the Software is ",
         "furnished to do so, subject to the following conditions:"),
  "",
  paste0("The above copyright notice and this permission notice shall be ",
         "included in all copies or substantial portions of the Software."),
  "",
  paste0("THE SOFTWARE IS PROVIDED \"AS IS\", WITHOUT WARRANTY OF ANY KIND, ",
         "EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF ",
         "MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND ",
         "NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT ",
         "HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, ",
         "WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, ",
         "OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER ",
         "DEALINGS IN THE SOFTWARE.")
)

writeLines(license, file.path(stage, "LICENSE"), useBytes = TRUE)

writeLines(c(".Rhistory", ".RData", ".Rproj.user/", "*.Rproj",
             "results/", "figures/", "output/",
             "*.pdf", "*.eps", "*.aux", "*.log", "*.bbl", "*.blg",
             "*.out", "*.synctex.gz"),
           file.path(stage, ".gitignore"), useBytes = TRUE)

writeLines(c("* text=auto eol=lf",
             "*.R text",
             "*.md text"),
           file.path(stage, ".gitattributes"), useBytes = TRUE)


# -----------------------------------------------------------------------------
# Git
#
# The working directory is inside the ASCII staging folder for every call and
# no path is ever passed to git.
# -----------------------------------------------------------------------------

old_wd <- getwd()
on.exit(setwd(old_wd), add = TRUE)

# system2 passes its arguments to the shell without quoting them, so on
# Windows an argument that contains a space arrives at git as several
# arguments and git reads them as file names. That is what happens to the
# commit message and to the author name, both of which contain spaces.
# Anything that is not a plain token is therefore quoted here.
git_quote <- function(a)
{
  if (grepl("[[:space:]\"]", a))
    shQuote(a, type = if (.Platform$OS.type == "windows") "cmd" else "sh") else
      a
}

git <- function(...)
{
  setwd(stage)

  args <- vapply(c(...), git_quote, character(1), USE.NAMES = FALSE)
  cat("git", args, "\n")

  st <- system2("git", args, stdout = TRUE, stderr = TRUE)

  if (length(st))
    cat(paste0("  ", st, collapse = "\n"), "\n")

  status <- attr(st, "status")

  if (!is.null(status) && status != 0L)
    stop("git ", paste(args, collapse = " "), " failed with status ", status)

  invisible(st)
}

if (Sys.which("git") == "")
  stop("git was not found on the PATH.")

remote <- paste0("https://github.com/", GH_USER, "/", GH_REPO, ".git")

git("init", "-q")
git("symbolic-ref", "HEAD", paste0("refs/heads/", BRANCH))
git("config", "user.name", AUTHOR_NAME)
git("config", "user.email", AUTHOR_EMAIL)
git("config", "core.quotepath", "false")

git("add", "-A")
git("commit", "-q", "-m",
    paste0("R code for: ", PAPER_TITLE))

git("remote", "add", "origin", remote)

cat("\nPushing to ", remote, "\n", sep = "")

if (force_push)
  git("push", "-u", "--force", "origin", BRANCH) else
    git("push", "-u", "origin", BRANCH)

setwd(old_wd)

cat("\nDone. https://github.com/", GH_USER, "/", GH_REPO, "\n", sep = "")
