# =============================================================================
# run_simulation.R
#
# The Monte Carlo study of Section 7, in parts.
#
# The study is sixty configurations: ten withdrawal schemes at each of three
# values of alpha for the table of alpha, and the same ten schemes at each of
# three values of beta for the table of beta. Two chains are run per
# replication, one under the informative prior of the chapter and one under a
# nearly flat prior, because the chapter reports the two side by side, so that
# the contribution of the prior to the performance of the Bayes estimator can
# be read off.
#
# The whole study takes about three and a half hours, which is longer than one
# sitting. It is therefore carried out in parts. Set PART below to 1, run the
# script, set it to 2 at the next opportunity, and so on. Each configuration
# writes its own file under results/sim_parts as soon as it finishes, and a
# configuration whose file is already there is skipped, so an interrupted run
# loses at most the configuration it was working on and nothing else. When the
# sixty files are complete the script assembles the report and the LaTeX of
# the two tables. Setting PART to 0 runs everything that is still missing.
#
# The result of a configuration does not depend on the part it was run in, on
# the number of workers, or on the machine. Each replication carries a stream
# of its own, derived on the master from the master seed and from the position
# of the configuration in a fixed list, so replication j of configuration i is
# the same experiment however the work is divided.
#
# Writes results/simulation_report.txt, and results/sim_parts/*.rds, which are
# working files and may be deleted once the report is complete.
# =============================================================================

rm(list = ls())

local({
  f <- "code/core_ipd_prcs.R"
  if (!file.exists(f)) f <- "core_ipd_prcs.R"
  source(f)
})

setup_project()
source("code/core_ipd_prcs.R")

library(parallel)


# -----------------------------------------------------------------------------
# Settings
#
# Every quantity here is reported in the txt file and has to appear in the
# chapter, because the tables cannot be reproduced without it.
# -----------------------------------------------------------------------------

PART  <- 1L        # 1, 2, 3, 4 in turn; 0 runs everything still missing
PARTS <- 4L        # how many parts the sixty configurations are dealt into

QUICK <- FALSE     # TRUE exercises every code path in a few minutes

N     <- if (QUICK)   40L else 1000L    # replications per configuration
ZETA  <- if (QUICK) 2000L else 10000L   # MCMC iterations
ZETA0 <- if (QUICK)  400L else  2000L   # burn-in

CONF  <- 0.95

# The value at which the parameter that is not varying is held.
BETA_FIXED  <- 1.0     # for the table of alpha
ALPHA_FIXED <- 1.5     # for the table of beta

ALPHA_GRID <- c(1.5, 2.5, 0.5)
BETA_GRID  <- c(0.5, 1.0, 1.5)

N_CORES <- max(1L, detectCores() - 1L)

# A trial run under QUICK writes to a directory of its own, so that its
# small files can never be mistaken for results of the real study.
PARTS_DIR <- file.path("results",
                       if (QUICK) "sim_parts_quick" else "sim_parts")

dir.create(PARTS_DIR, showWarnings = FALSE, recursive = TRUE)


# The informative prior of the chapter: the mean of each gamma equals the true
# value of the parameter. A value outside the table falls back on the same
# rule with rate 2.
PRIOR_ALPHA <- list("1.5" = c(shape = 3.00, rate = 2.0),
                    "2.5" = c(shape = 5.00, rate = 2.0),
                    "0.5" = c(shape = 1.25, rate = 2.5))

PRIOR_BETA  <- list("0.5" = c(shape = 1.25, rate = 2.5),
                    "1"   = c(shape = 2.00, rate = 2.0),
                    "1.5" = c(shape = 3.00, rate = 2.0))

# The nearly flat prior. Small equal shape and rate leave the posterior
# governed by the likelihood alone.
FLAT <- c(shape = 0.001, rate = 0.001)

prior_alpha <- function(alpha)
{
  key <- format(alpha)

  if (!is.null(PRIOR_ALPHA[[key]])) PRIOR_ALPHA[[key]] else
    c(shape = alpha * 2, rate = 2)
}

prior_beta <- function(beta_par)
{
  key <- format(beta_par)

  if (!is.null(PRIOR_BETA[[key]])) PRIOR_BETA[[key]] else
    c(shape = beta_par * 2, rate = 2)
}


# -----------------------------------------------------------------------------
# The sixty configurations, in a fixed order
#
# The order fixes the seed of every configuration, so it must never be
# changed once a run has been started.
# -----------------------------------------------------------------------------

CONFIGS <- list()

for (a in ALPHA_GRID)
  for (nm in names(SCHEMES_SIM))
    CONFIGS[[length(CONFIGS) + 1L]] <-
      list(table = "alpha", scheme = nm, alpha = a, beta = BETA_FIXED)

for (b in BETA_GRID)
  for (nm in names(SCHEMES_SIM))
    CONFIGS[[length(CONFIGS) + 1L]] <-
      list(table = "beta", scheme = nm, alpha = ALPHA_FIXED, beta = b)

for (i in seq_along(CONFIGS))
{
  cf <- CONFIGS[[i]]

  CONFIGS[[i]]$index <- i
  CONFIGS[[i]]$part  <- ((i - 1L) %% PARTS) + 1L
  CONFIGS[[i]]$id    <- sprintf("%02d_%s_%s_a%s_b%s", i, cf$table, cf$scheme,
                                format(cf$alpha), format(cf$beta))
  CONFIGS[[i]]$file  <- file.path(PARTS_DIR,
                                  paste0(CONFIGS[[i]]$id, ".rds"))
}

done <- vapply(CONFIGS, function(cf) file.exists(cf$file), logical(1))

todo <- if (PART == 0L) which(!done) else
  which(!done & vapply(CONFIGS, function(cf) cf$part == PART, logical(1)))


# A rough cost model, so that the script can say how long it expects to run.
# The cost of a configuration grows with m, and the two chains together cost
# about two and a half times one chain. The constants come from the timings
# of earlier runs of this script.
cost_secs <- function(cf)
{
  m <- SCHEMES_SIM[[cf$scheme]]$m

  2.6 * (15.1 + 1.70 * m) * (N / 1000) * (ZETA / 10000)
}


say <- new_report(file.path("results", "simulation_report.txt"))

# tex_table returns a character vector, one line per element.
emit <- function(lines)
{
  for (l in lines) say(l)
  say("")
}

say(rule())
say("MONTE CARLO SIMULATION STUDY")
say("Inverse Pareto distribution under progressive Type-II random censoring")
say(rule())
say("Date            : ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"))
say("R version       : ", R.version.string)
say("Master seed     : ", MASTER_SEED)
say("Replications N  : ", N)
say("MCMC            : zeta = ", ZETA, ", burn-in = ", ZETA0,
    ", two chains per replication")
say("Nominal level   : ", CONF)
say("beta held fixed at   : ", BETA_FIXED, "  (table for alpha)")
say("alpha held fixed at  : ", ALPHA_FIXED, "  (table for beta)")
say("Workers         : ", N_CORES)
say("")
say("Priors, informative, the mean of each gamma equal to the true value:")

for (a in ALPHA_GRID)
  say(sprintf("  alpha = %.1f : Gamma(%.2f, %.2f)", a,
              prior_alpha(a)["shape"], prior_alpha(a)["rate"]))

for (b in BETA_GRID)
  say(sprintf("  beta  = %.1f : Gamma(%.2f, %.2f)", b,
              prior_beta(b)["shape"], prior_beta(b)["rate"]))

say(sprintf("Prior, nearly flat  : Gamma(%.3f, %.3f) for both parameters",
            FLAT["shape"], FLAT["rate"]))
say("")
say("Parts           : ", PARTS, ", this run is part ",
    if (PART == 0L) "0, everything still missing" else PART)
say("  A configuration is assigned to a part by its position, so PARTS must")
say("  not be changed once a run has been started. Should it be changed by")
say("  accident, a final run with PART set to 0 picks up whatever was left.")
say("Configurations  : ", length(CONFIGS), " in all, ", sum(done),
    " already finished, ", length(todo), " to be run now")

if (length(todo))
  say(sprintf("Expected run time of this part: about %.0f minutes on %d workers",
              sum(vapply(CONFIGS[todo], cost_secs, numeric(1))) / 60, N_CORES))

say("")


# -----------------------------------------------------------------------------
# One stream per replication
#
# The streams of a configuration are derived from the master seed and the
# position of the configuration, so they do not depend on the number of
# workers or on the part in which the configuration was run.
# -----------------------------------------------------------------------------

RNGkind("L'Ecuyer-CMRG")

make_streams <- function(seed, n)
{
  set.seed(seed)

  s <- vector("list", n)
  s[[1L]] <- .Random.seed

  for (i in seq_len(n - 1L))
    s[[i + 1L]] <- nextRNGStream(s[[i]])

  s
}


# -----------------------------------------------------------------------------
# One replication
#
# One data set, one fit and two chains. Every interval that has been in
# question is computed here from that same fit and those same chains, at no
# further cost, so that a question about the choice of interval never requires
# the study to be run again.
# -----------------------------------------------------------------------------

one_rep <- function(stream, n, m, W, alpha, beta_par, pa, pb)
{
  assign(".Random.seed", stream, envir = globalenv())

  tryCatch(
    {
      dd <- gen_prcs(n, m, W, alpha, beta_par)

      if (dd$m == dd$k)
        return(list(ok = FALSE, why = "no random censoring event"))

      ft <- fit_mle(dd, conf = CONF)

      if (!isTRUE(ft$ok))
        return(list(ok = FALSE, why = ft$msg))

      p <- (1 - CONF) / 2

      ci <- mcmc_prcs(dd, pa["shape"], pa["rate"], pb["shape"], pb["rate"],
                      zeta = ZETA, zeta0 = ZETA0,
                      start = c(ft$alpha, ft$beta))

      cf <- mcmc_prcs(dd, FLAT["shape"], FLAT["rate"],
                      FLAT["shape"], FLAT["rate"],
                      zeta = ZETA, zeta0 = ZETA0,
                      start = c(ft$alpha, ft$beta))

      qa <- unname(quantile(ci$alpha, c(p, 1 - p)))
      qb <- unname(quantile(ci$beta,  c(p, 1 - p)))
      fa <- unname(quantile(cf$alpha, c(p, 1 - p)))
      fb <- unname(quantile(cf$beta,  c(p, 1 - p)))

      ha <- hpd(ci$alpha, CONF); hb <- hpd(ci$beta, CONF)
      ga <- hpd(cf$alpha, CONF); gb <- hpd(cf$beta, CONF)

      list(ok = TRUE,

           # alpha
           a_mle_log  = c(ft$alpha, ft$ci_log[1, 1],   ft$ci_log[1, 2]),
           a_mle_wald = c(ft$alpha, ft$ci_plain[1, 1], ft$ci_plain[1, 2]),
           a_inf_hpd  = c(mean(ci$alpha), ha[1], ha[2]),
           a_inf_eqt  = c(mean(ci$alpha), qa[1], qa[2]),
           a_flt_hpd  = c(mean(cf$alpha), ga[1], ga[2]),
           a_flt_eqt  = c(mean(cf$alpha), fa[1], fa[2]),

           # beta
           b_mle_log  = c(ft$beta, ft$ci_log[2, 1],   ft$ci_log[2, 2]),
           b_mle_wald = c(ft$beta, ft$ci_plain[2, 1], ft$ci_plain[2, 2]),
           b_inf_hpd  = c(mean(ci$beta), hb[1], hb[2]),
           b_inf_eqt  = c(mean(ci$beta), qb[1], qb[2]),
           b_flt_hpd  = c(mean(cf$beta), gb[1], gb[2]),
           b_flt_eqt  = c(mean(cf$beta), fb[1], fb[2]),

           acc_inf = ci$acc_rate,
           acc_flt = cf$acc_rate,
           sig_inf = ci$sigma_alpha,
           score   = max(abs(ft$score)),
           pd      = isTRUE(ft$pd))
    },
    error = function(e) list(ok = FALSE, why = conditionMessage(e))
  )
}


VARIANTS <- c("mle_log", "mle_wald", "inf_hpd", "inf_eqt", "flt_hpd", "flt_eqt")


summarise_one <- function(est, lo, hi, truth)
{
  nn  <- length(est)
  err <- est - truth
  sq  <- err^2
  len <- hi - lo
  cov <- as.numeric(lo <= truth & truth <= hi)

  c(AE  = mean(est), AB = mean(err), MSE = mean(sq),
    AL  = mean(len), CP = mean(cov),
    seAE = sd(est) / sqrt(nn), seAB = sd(err) / sqrt(nn),
    seMSE = sd(sq) / sqrt(nn), seAL = sd(len) / sqrt(nn),
    seCP = sqrt(mean(cov) * (1 - mean(cov)) / nn))
}


# -----------------------------------------------------------------------------
# The cluster
#
# The core is not sourced on the workers. A worker cannot source a file whose
# path contains non-ASCII characters on Windows, which is exactly the case for
# this project, so the objects of the core are exported instead.
# -----------------------------------------------------------------------------

if (length(todo))
{
  # A trial replication on the master first. An error in the model or in the
  # settings then shows itself in a few seconds, rather than after the cluster
  # has been built and the first configuration has been dealt out.
  local({
    cf <- CONFIGS[[todo[1]]]
    sc <- SCHEMES_SIM[[cf$scheme]]
    st <- make_streams(MASTER_SEED + 1000L * todo[1], 1L)

    tr <- one_rep(st[[1L]], sc$n, sc$m, sc$W, cf$alpha, cf$beta,
                  prior_alpha(cf$alpha), prior_beta(cf$beta))

    if (!isTRUE(tr$ok))
    {
      say("Trial replication failed: ", tr$why)
      stop("The trial replication failed. The cluster was not started.")
    }

    say(sprintf("Trial replication succeeded on %s: MLE = %.4f, Bayes = %.4f informative, %.4f flat",
                cf$id, tr$a_mle_log[1], tr$a_inf_hpd[1], tr$a_flt_hpd[1]))
    say("")
  })

  cl <- makeCluster(N_CORES)

  clusterExport(cl,
                varlist = setdiff(ls(envir = globalenv()), c("say", "cl")),
                envir = globalenv())

  for (i in todo)
  {
    cf <- CONFIGS[[i]]
    sc <- SCHEMES_SIM[[cf$scheme]]

    pa <- prior_alpha(cf$alpha)
    pb <- prior_beta(cf$beta)

    st <- make_streams(MASTER_SEED + 1000L * i, N)

    t0  <- Sys.time()
    res <- parLapply(cl, st, one_rep, n = sc$n, m = sc$m, W = sc$W,
                     alpha = cf$alpha, beta_par = cf$beta, pa = pa, pb = pb)
    el  <- as.numeric(difftime(Sys.time(), t0, units = "secs"))

    ok   <- vapply(res, function(r) isTRUE(r$ok), logical(1))
    bad  <- res[!ok]
    res  <- res[ok]

    out <- list(index = i, id = cf$id, table = cf$table, scheme = cf$scheme,
                n = sc$n, m = sc$m, pattern = sc$pattern,
                alpha = cf$alpha, beta = cf$beta,
                prior_alpha = pa, prior_beta = pb,
                nused = length(res), nbad = length(bad), secs = el,
                acc_inf = mean(vapply(res, function(r) r$acc_inf, numeric(1))),
                acc_flt = mean(vapply(res, function(r) r$acc_flt, numeric(1))),
                sig_inf = mean(vapply(res, function(r) r$sig_inf, numeric(1))),
                score = max(vapply(res, function(r) r$score, numeric(1))),
                allpd = all(vapply(res, function(r) r$pd, logical(1))),
                reasons = if (length(bad))
                  table(vapply(bad, function(r) r$why, character(1))) else NULL,
                summ = list())

    for (par in c("a", "b"))
    {
      truth <- if (par == "a") cf$alpha else cf$beta

      for (v in VARIANTS)
      {
        key <- paste0(par, "_", v)
        g   <- function(j) vapply(res, function(r) r[[key]][j], numeric(1))

        out$summ[[key]] <- summarise_one(g(1), g(2), g(3), truth)
      }
    }

    saveRDS(out, cf$file)

    say(sprintf("finished %2d/%d  %-28s  n = %2d, m = %2d, %5.1f s, discarded %d",
                i, length(CONFIGS), cf$id, sc$n, sc$m, el, length(bad)))
  }

  stopCluster(cl)
  say("")
}


# -----------------------------------------------------------------------------
# The report
#
# Written from the files on disk, so it reflects everything finished so far
# whichever part produced it.
# -----------------------------------------------------------------------------

have <- vapply(CONFIGS, function(cf) file.exists(cf$file), logical(1))

say(rule())
say("RESULTS")
say(rule())
say(sum(have), " of ", length(CONFIGS), " configurations are on disk.")

if (!all(have))
{
  say("")
  say("Still missing, by part:")

  for (p in seq_len(PARTS))
  {
    miss <- vapply(CONFIGS, function(cf) cf$part == p, logical(1)) & !have

    if (any(miss))
      say(sprintf("  part %d : %d configurations", p, sum(miss)))
  }

  say("")
  say("Set PART to a part that is still missing and run this script again.")
}

say("")

R <- lapply(CONFIGS[have], function(cf) readRDS(cf$file))

LAB <- c(mle_log  = "MLE, log-transformed interval",
         mle_wald = "MLE, plain Wald interval",
         inf_hpd  = "Bayes, informative prior, HPD interval",
         inf_eqt  = "Bayes, informative prior, equal-tailed interval",
         flt_hpd  = "Bayes, nearly flat prior, HPD interval",
         flt_eqt  = "Bayes, nearly flat prior, equal-tailed interval")

for (tb in c("alpha", "beta"))
{
  sel <- Filter(function(r) r$table == tb, R)

  if (!length(sel))
    next

  par <- if (tb == "alpha") "a" else "b"

  say(rule("-"))
  say("Estimation of ", tb,
      if (tb == "alpha") "  (beta held fixed)" else "  (alpha held fixed)")
  say(rule("-"))

  for (r in sel)
  {
    say(sprintf("%s  n = %2d, m = %2d, %-11s  true alpha = %.1f, true beta = %.1f   used %d of %d, %.1f s",
                r$scheme, r$n, r$m, r$pattern, r$alpha, r$beta, r$nused, N,
                r$secs))
    say(sprintf("  acceptance %.4f informative, %.4f flat; proposal sd %.4f; max |score| %.2e; information positive definite %s",
                r$acc_inf, r$acc_flt, r$sig_inf, r$score, r$allpd))

    if (!is.null(r$reasons))
      for (j in seq_along(r$reasons))
        say("  discarded, ", names(r$reasons)[j], ": ", r$reasons[j])

    # Both parameters are reported for every configuration, and the parameter
    # of the table is marked. The estimation of alpha at a given value of beta
    # appears in the configurations of the other table as well, which is what
    # makes it possible to tell whether a block of the chapter was run at the
    # value of beta it claims.
    for (pp in c("a", "b"))
    {
      lab <- if (pp == "a") "alpha" else "beta"
      mk  <- if (pp == par) "  (the parameter of this table)" else ""

      say("    estimation of ", lab, mk)
      say("      variant                                              AE        AB       MSE        AL       CP")

      for (v in VARIANTS)
      {
        s  <- r$summ[[paste0(pp, "_", v)]]

        say(sprintf("      %-48s %9.4f %9.4f %9.4f %9.4f %8.3f",
                    LAB[v], s["AE"], s["AB"], s["MSE"], s["AL"], s["CP"]))

        say(sprintf("      %-48s %9.4f %9.4f %9.4f %9.4f %8.3f",
                    "    standard error", s["seAE"], s["seAB"], s["seMSE"],
                    s["seAL"], s["seCP"]))
      }

      say("")
    }
  }
}


# -----------------------------------------------------------------------------
# The LaTeX of the tables
#
# Four tables: for each parameter, one for the maximum likelihood estimator and
# one for the Bayes estimator under the two priors side by side.
#
# The two priors are reported together because the comparison between them is
# itself a finding. The informative prior has its mean at the true value, and
# under it the Bayes estimator attains a smaller mean squared error and a
# shorter interval than the MLE while covering more often than the nominal
# level. A procedure cannot ordinarily do all three at once, so the reader is
# entitled to see what the prior contributes, and the nearly flat prior, under
# which the Bayes estimator very nearly reproduces the MLE, is the reference
# against which that contribution is read.
#
# The confidence interval is the log-transformed one and the credible interval
# is the HPD interval. The remaining combinations are in the report above.
# -----------------------------------------------------------------------------

MLE_VAR <- "mle_log"
BAY_INF <- "inf_hpd"
BAY_FLT <- "flt_hpd"

if (all(have))
{
  say(rule())
  say("LATEX OF THE TABLES")
  say(rule())
  say("Intervals: ", LAB[MLE_VAR], ", and ", LAB[BAY_INF], " or ",
      LAB[BAY_FLT], ".")
  say("")

  for (tb in c("alpha", "beta"))
  {
    sel <- Filter(function(r) r$table == tb, R)
    par <- if (tb == "alpha") "a" else "b"
    sym <- if (tb == "alpha") "\\alpha" else "\\beta"
    tag <- if (tb == "alpha") "alpha" else "beta"

    # the first five columns, shared by both tables of this parameter
    stub <- character(0)
    prev <- NA_real_

    for (r in sel)
    {
      truth <- if (tb == "alpha") r$alpha else r$beta
      shown <- if (is.na(prev) || truth != prev) sprintf("%.1f", truth) else ""
      prev  <- truth

      stub <- c(stub, sprintf("%s & %d & %d & %s & %s", r$scheme, r$n, r$m,
                              r$pattern, shown))
    }

    # Table of the maximum likelihood estimator
    rows <- character(0)

    for (j in seq_along(sel))
    {
      s1 <- sel[[j]]$summ[[paste0(par, "_", MLE_VAR)]]

      rows <- c(rows, sprintf("%s & %.4f & %.4f & %.4f & %.4f & %.3f \\\\",
                              stub[j], s1["AE"], s1["AB"], s1["MSE"],
                              s1["AL"], s1["CP"]))
    }

    emit(tex_table(
      body = rows,
      header = paste0("Scheme & $n$ & $m$ & Pattern & $", sym,
                      "$ & AE & AB & MSE & AL & CP \\\\"),
      align = "ccclcccccc",
      caption = paste0("Simulation results for the maximum likelihood ",
                       "estimator under different PRCS schemes and values of $",
                       sym, "$. AL is the average length of the ",
                       "log-transformed confidence interval and CP is its ",
                       "coverage probability."),
      label = paste0("tab:sim_mle_", tag)))

    say("")

    # Table of the Bayes estimator under the two priors
    rows <- character(0)

    for (j in seq_along(sel))
    {
      s1 <- sel[[j]]$summ[[paste0(par, "_", BAY_INF)]]
      s2 <- sel[[j]]$summ[[paste0(par, "_", BAY_FLT)]]

      rows <- c(rows,
                sprintf(paste0("%s & %.4f & %.4f & %.4f & %.4f & %.3f",
                               " & %.4f & %.4f & %.4f & %.4f & %.3f \\\\"),
                        stub[j],
                        s1["AE"], s1["AB"], s1["MSE"], s1["AL"], s1["CP"],
                        s2["AE"], s2["AB"], s2["MSE"], s2["AL"], s2["CP"]))
    }

    emit(tex_table(
      body = rows,
      header = paste0("\\multirow{2}{*}{Scheme} & \\multirow{2}{*}{$n$} & ",
                      "\\multirow{2}{*}{$m$} & \\multirow{2}{*}{Pattern} & ",
                      "\\multirow{2}{*}{$", sym, "$} & ",
                      "\\multicolumn{5}{c}{Informative prior} & ",
                      "\\multicolumn{5}{c}{Nearly flat prior} \\\\\n",
                      "\\cmidrule(lr){6-10}\\cmidrule(lr){11-15}\n",
                      "& & & & & AE & AB & MSE & AL & CP & ",
                      "AE & AB & MSE & AL & CP \\\\"),
      align = "ccclccccccccccc",
      caption = paste0("Simulation results for the Bayes estimator under the ",
                       "two priors, for different PRCS schemes and values of $",
                       sym, "$. AL is the average length of the HPD credible ",
                       "interval and CP is its coverage probability."),
      label = paste0("tab:sim_bayes_", tag)))

    say("")
  }
}

say(rule())
say("run_simulation.R finished")
say(rule())
