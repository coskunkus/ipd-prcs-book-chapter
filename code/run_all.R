# =============================================================================
# run_all.R
#
# Everything in the chapter except the Monte Carlo study, which has its own
# runner because it is the only long computation.
#
# Produces
#   figures/IPD_CDF_convergence.eps and .pdf   Figure 1 of the chapter
#   results/run_all_report.txt                 one report, with the LaTeX of
#                                              every table it produces
#
# Reproducibility. Every number below comes from the portable generator
# defined in the core, so the results do not depend on R's random number
# generator. Running this script again, on any machine and under any version
# of R, gives the same digits.
#
# Run time: well under a minute.
# =============================================================================

rm(list = ls())
graphics.off()

local({
  f <- "code/core_ipd_prcs.R"
  if (!file.exists(f)) f <- "core_ipd_prcs.R"
  source(f)
})

setup_project()
source("code/core_ipd_prcs.R")

say <- new_report(file.path("results", "run_all_report.txt"))


# -----------------------------------------------------------------------------
# Settings for the application of Section 8
#
# The data of Nelson (1982) are complete, so a progressive Type-II random
# censoring design is imposed on them for illustration. Doing that needs a
# censoring parameter and a stream of pseudo-random numbers, and both are
# fixed here and reported, so that the tables can be reproduced.
# -----------------------------------------------------------------------------

BETA_IMPOSE <- 1            # equal chance of a failure and a censoring event
GEN_SEED    <- 20260909     # seed of the portable generator

PRIOR_SHAPE <- 2            # integer, so that the gamma full conditional of
                            # beta has an integer shape

ZETA  <- 20000L
ZETA0 <- 4000L

MISSION_TIMES <- c(5, 10, 30)   # minutes

# Parametric bootstrap of the Kolmogorov-Smirnov statistic. The shape
# parameter is estimated from the sample, so the null distribution of the
# statistic is not the tabulated one and a critical value read from a table
# would not be valid here. The bootstrap is run from a generator of its own,
# started from a seed of its own, so that adding it leaves every other number
# in this report unchanged.
BOOT_KS   <- 2000L
BOOT_SEED <- 20260910


say(rule())
say("APPLICATION, FIGURE AND TABLES OF THE CHAPTER")
say("Inverse Pareto distribution under progressive Type-II random censoring")
say(rule())
say("Date        : ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"))
say("R version   : ", R.version.string)
say("Generator   : Lehmer, x <- 16807 x mod 2147483647, seed ", GEN_SEED)
say("             (fully specified by the seed, so no result depends on R)")
say("beta used to impose the random censoring : ", BETA_IMPOSE)
say("MCMC        : zeta = ", ZETA, ", burn-in = ", ZETA0)
say("Bootstrap   : ", BOOT_KS, " parametric replications, seed ", BOOT_SEED)
say("")


# =============================================================================
# 1. Figure 1, the convergence of the extreme order statistics
# =============================================================================

say(rule("-"))
say("1. Figure 1, convergence of the normalized extreme order statistics")
say(rule("-"))

alpha_fig <- 3
i_upper   <- 3

n_outer  <- c(5, 10, 25, 100)
n_middle <- c(25, 100, 1000, 100000)

cols <- c("#E69F00", "#D55E00", "#0072B2", "#009E73")
ltys <- c(2, 3, 4, 5)

bU <- function(n) q_ipd(1 - 1 / n, alpha_fig)
bL <- function(n) q_ipd(1 / n, alpha_fig)

cdf_max_n <- function(x, n)
{
  z <- pmin(pmax(p_ipd(bU(n) * x, alpha_fig), 0), 1)

  out <- numeric(length(z))
  ok  <- z > 0

  out[ok]     <- exp(n * log(z[ok]))
  out[z >= 1] <- 1

  pmin(pmax(out, 0), 1)
}

cdf_min_n <- function(x, n)
{
  z <- pmin(pmax(p_ipd(bL(n) * x, alpha_fig), 0), 1)

  out <- -expm1(n * log1p(-z))

  out[z <= 0] <- 0
  out[z >= 1] <- 1

  pmin(pmax(out, 0), 1)
}

cdf_upper_n <- function(x, n, i)
{
  Fy <- p_ipd(bU(n) * x, alpha_fig)

  pbinom(i - 1, size = n, prob = pmin(pmax(1 - Fy, 0), 1))
}

cdf_upper_lim <- function(x, i)
{
  out <- numeric(length(x))
  ok  <- x > 0
  xx  <- x[ok]

  s <- rep(0, length(xx))

  for (j in 0:(i - 1))
    s <- s + xx^(-j) / factorial(j)

  out[ok] <- exp(-1 / xx) * s

  out
}

x_max   <- exp(seq(log(0.1), log(100), length.out = 3000))
x_min   <- seq(0, 2.2, length.out = 2000)
x_upper <- seq(0, 3, length.out = 3000)

draw_fig1 <- function()
{
  par(mfrow = c(1, 3), mar = c(4.4, 4.5, 2.8, 1.0),
      oma = c(0, 0, 0.5, 0), mgp = c(2.5, 0.8, 0))

  plot(x_max, exp(-1 / x_max), type = "l", lwd = 3, col = "black",
       ylim = c(0, 1), log = "x", xlab = "x (log scale)", ylab = "CDF",
       main = "(a) Sample maximum", xaxs = "i", yaxs = "i")

  for (k in seq_along(n_outer))
    lines(x_max, cdf_max_n(x_max, n_outer[k]), col = cols[k], lwd = 1.8,
          lty = ltys[k])

  legend("bottomright",
         legend = c("Limit: exp(-1/x)", paste0("n = ", n_outer)),
         col = c("black", cols), lty = c(1, ltys),
         lwd = c(3, rep(1.8, length(n_outer))), bty = "n", cex = 0.78)

  plot(x_min, -expm1(-x_min^alpha_fig), type = "l", lwd = 3, col = "black",
       ylim = c(0, 1), xlab = "x", ylab = "CDF",
       main = paste0("(b) Sample minimum, alpha = ", alpha_fig))

  for (k in seq_along(n_middle))
    lines(x_min, cdf_min_n(x_min, n_middle[k]), col = cols[k], lwd = 1.8,
          lty = ltys[k])

  legend("bottomright",
         legend = c(paste0("Limit: 1-exp(-x^", alpha_fig, ")"),
                    paste0("n = ", n_middle)),
         col = c("black", cols), lty = c(1, ltys),
         lwd = c(3, rep(1.8, length(n_middle))), bty = "n", cex = 0.78)

  plot(x_upper, cdf_upper_lim(x_upper, i_upper), type = "l", lwd = 3,
       col = "black", xlim = c(0, 3), ylim = c(0, 1), xlab = "x",
       ylab = "CDF",
       main = paste0("(c) Fixed upper extreme, i = ", i_upper),
       xaxs = "i", yaxs = "i")

  for (k in seq_along(n_outer))
    lines(x_upper, cdf_upper_n(x_upper, n_outer[k], i_upper), col = cols[k],
          lwd = 1.8, lty = ltys[k])

  legend("bottomright",
         legend = c("Limit: exp(-1/x)[1 + 1/x + 1/(2x^2)]",
                    paste0("n = ", n_outer)),
         col = c("black", cols), lty = c(1, ltys),
         lwd = c(3, rep(1.8, length(n_outer))), bty = "n", cex = 0.72)
}

save_figure("IPD_CDF_convergence", draw_fig1, width = 14.5, height = 4.8)

say("Written: figures/IPD_CDF_convergence.eps and .pdf")
say("")
say("Largest absolute difference between the finite-sample cdf and its limit")

for (n in n_outer)
  say(sprintf("  sample maximum,      n = %-6d  %.8f", n,
              max(abs(cdf_max_n(x_max, n) - exp(-1 / x_max)))))

for (n in n_middle)
  say(sprintf("  sample minimum,      n = %-6d  %.8f", n,
              max(abs(cdf_min_n(x_min, n) - (-expm1(-x_min^alpha_fig))))))

for (n in n_outer)
  say(sprintf("  fixed upper, i = %d,  n = %-6d  %.8f", i_upper, n,
              max(abs(cdf_upper_n(x_upper, n, i_upper) -
                        cdf_upper_lim(x_upper, i_upper)))))

say("")


# =============================================================================
# 2. The complete data set
# =============================================================================

say(rule("-"))
say("2. The data of Nelson (1982), Table 1.1 of Chapter 3, the 34 kV column")
say(rule("-"))

x <- BREAKDOWN34

alpha_c <- mle_complete_ipd(x)

say("Times to breakdown of an insulating fluid at 34 kV, in minutes")
say("  ", paste(formatC(x, format = "f", digits = 2), collapse = ", "))
say(sprintf("  n = %d, sum = %.2f, mean = %.4f (Nelson gives 14.3 on page 210)",
            length(x), sum(x), mean(x)))
say(sprintf("  median = %.2f, min = %.2f, max = %.2f",
            median(x), min(x), max(x)))
say("")
say(sprintf("Complete-sample MLE of alpha = %.4f", alpha_c))

# goodness of fit of the fitted IPD to the complete sample
Fh  <- p_ipd(sort(x), alpha_c)
nn  <- length(x)
Dp  <- max(seq_len(nn) / nn - Fh)
Dm  <- max(Fh - (seq_len(nn) - 1) / nn)
KS  <- max(Dp, Dm)

say(sprintf("Kolmogorov-Smirnov distance of the fitted IPD = %.4f", KS))

# Parametric bootstrap p-value. Each replication draws a sample of the same
# size from the fitted model, re-estimates the shape parameter from that
# sample and recomputes the statistic, so that the estimation of the parameter
# is reproduced inside the null distribution.
ks_stat <- function(v, a)
{
  Fv <- p_ipd(sort(v), a)
  nv <- length(v)

  max(max(seq_len(nv) / nv - Fv), max(Fv - (seq_len(nv) - 1) / nv))
}

gen_boot <- lcg_new(BOOT_SEED)
KSb      <- numeric(BOOT_KS)

for (b in seq_len(BOOT_KS))
{
  xb      <- q_ipd(lcg_unif(gen_boot, nn), alpha_c)
  KSb[b]  <- ks_stat(xb, mle_complete_ipd(xb))
}

KS_pval <- (1 + sum(KSb >= KS)) / (BOOT_KS + 1)

say(sprintf("Parametric bootstrap p-value, %d replications = %.4f",
            BOOT_KS, KS_pval))
say("")


# =============================================================================
# 3. The imposed censoring designs and the estimates
# =============================================================================

say(rule("-"))
say("3. Imposed progressive Type-II random censoring designs")
say(rule("-"))

a1 <- PRIOR_SHAPE
b1 <- PRIOR_SHAPE / alpha_c
a2 <- PRIOR_SHAPE
b2 <- PRIOR_SHAPE / BETA_IMPOSE

say(sprintf("Gamma priors: alpha ~ Gamma(%.4f, %.4f), beta ~ Gamma(%.4f, %.4f)",
            a1, b1, a2, b2))
say("  the prior means are the complete-sample estimate of alpha and the")
say("  value of beta used to impose the censoring")
say("")

rows_sch <- character(0)
rows_par <- character(0)
rows_rel <- character(0)

for (nm in names(SCHEMES_E))
{
  sc   <- SCHEMES_E[[nm]]

  # the extra space after a block is omitted on the last one, so that it
  # does not push against the bottom rule
  last <- identical(nm, names(SCHEMES_E)[length(SCHEMES_E)])

  # one generator per scheme, seeded from GEN_SEED and the scheme number,
  # so that the schemes do not share a stream
  gen <- lcg_new(GEN_SEED + as.integer(sub("E", "", nm)))

  dd <- impose_prcs(x, sc$m, sc$W, BETA_IMPOSE, alpha_impose = alpha_c,
                    gen = gen)

  ft <- fit_mle(dd)

  if (!isTRUE(ft$ok))
  {
    say("Scheme ", nm, ": ", ft$msg)
    next
  }

  ch <- mcmc_prcs(dd, a1, b1, a2, b2, zeta = ZETA, zeta0 = ZETA0, gen = gen)

  hpd_a <- hpd(ch$alpha)
  hpd_b <- hpd(ch$beta)

  # the same posterior means by one-dimensional quadrature, which is exact
  # because beta can be integrated out of the posterior
  qm <- post_moments_quad(dd, a1, b1, a2, b2, upper = max(200, 8 * alpha_c))

  say("Scheme ", nm, " (", sc$pattern, "), k = ", dd$k,
      " failures and ", dd$m - dd$k, " random censoring events")
  say("  y = ", paste(formatC(dd$y, format = "f", digits = 4),
                      collapse = ", "))
  say("  d = ", paste(dd$d, collapse = ", "))
  say("  W = ", paste(dd$W, collapse = ", "))
  say(sprintf("  MLE   alpha = %.4f  (%.4f, %.4f)   beta = %.4f  (%.4f, %.4f)",
              ft$alpha, ft$ci_log[1, 1], ft$ci_log[1, 2],
              ft$beta, ft$ci_log[2, 1], ft$ci_log[2, 2]))
  say(sprintf("  Bayes alpha = %.4f  (%.4f, %.4f)   beta = %.4f  (%.4f, %.4f)",
              mean(ch$alpha), hpd_a[1], hpd_a[2],
              mean(ch$beta), hpd_b[1], hpd_b[2]))
  say(sprintf("  log-likelihood = %.6f, max |score| = %.3e", ft$loglik,
              max(abs(ft$score))),
      ", information matrix positive definite: ", ft$pd)
  say(sprintf("  acceptance rate = %.4f, proposal sd = %.4f",
              ch$acc_rate, ch$sigma_alpha))
  say(sprintf("  check, posterior means by quadrature: alpha = %.4f, beta = %.4f",
              qm["mean_alpha"], qm["mean_beta"]))

  nw <- fit_mle_newton(dd, start = c(1, 1))

  if (isTRUE(nw$ok))
    say(sprintf("  check, Newton-Raphson agrees to %.3e in %d iterations",
                max(abs(nw$est - ft$est) / abs(ft$est)), nw$iter))
  else
    say("  Newton-Raphson: ", nw$msg)

  say("")

  # Only k is computed here; n, m and the withdrawal vector are the prefixed
  # design and are typeset in the chapter with braces over the repeated runs.
  # The row below therefore reports the design in plain form together with the
  # computed k, and the chapter's own table is checked against it.
  rows_sch <- c(rows_sch, paste0(
    nm, " & ", sc$n, " & ", sc$m, " & ", sc$n - sc$m, " & ",
    paste(sc$W, collapse = ","), " & ", sc$pattern, " & ", dd$k, " \\\\"))

  rows_par <- c(rows_par,
    paste0(nm, " & ", fmt4(ft$alpha), " & ", fmt4(mean(ch$alpha)),
           " & ", fmt4(ft$beta), " & ", fmt4(mean(ch$beta)), "\\\\"),
    paste0(" & (", fmt4(ft$ci_log[1, 1]), ",\\,", fmt4(ft$ci_log[1, 2]), ")",
           " & (", fmt4(hpd_a[1]), ",\\,", fmt4(hpd_a[2]), ")",
           " & (", fmt4(ft$ci_log[2, 1]), ",\\,", fmt4(ft$ci_log[2, 2]), ")",
           " & (", fmt4(hpd_b[1]), ",\\,", fmt4(hpd_b[2]), ")",
           if (last) "\\\\" else "\\\\[2pt]"),
    "")

  for (tm in MISSION_TIMES)
  {
    S_ml <- s_ipd(tm, ft$alpha)
    h_ml <- h_ipd(tm, ft$alpha)

    S_ps <- 1 - (tm / (1 + tm))^ch$alpha
    h_ps <- ch$alpha * tm^(ch$alpha - 1) /
      ((1 + tm) * ((1 + tm)^ch$alpha - tm^ch$alpha))

    hS <- hpd(S_ps)

    rows_rel <- c(rows_rel, paste0(
      nm, " & ", formatC(tm, width = 2, format = "d"),
      " & ", fmt4(S_ml), " & ", fmt4(h_ml),
      " & ", fmt4(mean(S_ps)),
      "\\; (", fmt4(hS[1]), ",\\,", fmt4(hS[2]), ")",
      " & ", fmt4(mean(h_ps)),
      if (tm != MISSION_TIMES[length(MISSION_TIMES)]) "\\\\" else
        if (last) "\\\\" else "\\\\[2pt]"))
  }
}


# =============================================================================
# 4. LaTeX for the tables of Section 8
# =============================================================================

say(rule("-"))
say("4. LaTeX for the tables of Section 8")
say(rule("-"))
say("Paste these blocks into the chapter. Do not retype the numbers from the")
say("printed output above, because rounding twice can move a cell by one digit.")
say("")

emit <- function(lines)
{
  for (l in lines) say(l)
  say("")
}


# The three tables of Section 8 are emitted in the exact form in which they
# stand in the chapter, wrapper included, so that a block can be pasted over
# the corresponding block of paper.tex without being reshaped by hand. The
# wrapper is the one the chapter uses throughout, a plain tabular with a
# stretched row height, and not the tabular* of the core emitter.
tex_plain <- function(body, header, align, caption, label, stretch)
{
  c("\\begin{table}[!h]",
    "\\centering",
    paste0("\\caption{", caption, "}"),
    paste0("\\label{", label, "}"),
    paste0("\\renewcommand{\\arraystretch}{", stretch, "}"),
    paste0("\\begin{tabular}{", align, "}"),
    "\\toprule",
    header,
    "\\midrule",
    "",
    body,
    "\\bottomrule",
    "\\end{tabular}",
    "\\end{table}")
}

say("The design of each scheme is prefixed, so only $k$ below is computed.")
say("The chapter typesets the withdrawal vector with braces over the repeated")
say("runs; check its k column against this one.")
say("")

emit(c("Scheme & n & m & n-m & W & Pattern & k", rows_sch))

emit(tex_plain(
  body    = rows_par,
  header  = paste0("\\multirow{2}{*}{} &\n",
                   "\\multicolumn{2}{c}{$\\alpha$} &\n",
                   "\\multicolumn{2}{c}{$\\beta$}\\\\\n",
                   "\\cmidrule(lr){2-3}\n",
                   "\\cmidrule(lr){4-5}\n",
                   "Scheme & MLE & Bayes & MLE & Bayes\\\\"),
  align   = "lcccc",
  caption = paste0("MLE and Bayes estimates of $\\alpha$ and $\\beta$ with ",
                   "95\\% intervals\nfor the insulating fluid data under ",
                   "different PRCS schemes."),
  label   = "tab:breakdown34_results",
  stretch = "1.20"))

emit(tex_plain(
  body    = rows_rel,
  header  = paste0("\\multirow{2}{*}{Scheme} & \\multirow{2}{*}{$t$} &\n",
                   "\\multicolumn{2}{c}{MLE} &\n",
                   "\\multicolumn{2}{c}{Bayes}\\\\\n",
                   "\\cmidrule(lr){3-4}\n",
                   "\\cmidrule(lr){5-6}\n",
                   "& & $\\widehat S_X(t)$ & $\\widehat h_X(t)$\n",
                   "& $\\widehat S_{X,B}(t)$ with HPD interval & ",
                   "$\\widehat h_{X,B}(t)$\\\\"),
  align   = "cccccc",
  caption = paste0("Estimates of the survival and hazard rate functions at ",
                   "three mission\ntimes, in minutes, for the insulating ",
                   "fluid data under different PRCS\nschemes."),
  label   = "tab:reliability_hazard",
  stretch = "1.15"))

say(rule())
say("run_all.R finished")
say(rule())
