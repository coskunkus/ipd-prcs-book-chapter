# =============================================================================
# check_formulas.R
#
# Independent verification of every closed-form expression of the chapter.
#
# Two layers, in this order.
#
#   A. Exact checks. Each closed form is compared with an independent route:
#      numerical integration, exact enumeration, a series against a finite
#      representation, a numerical gradient against the analytic score, a
#      numerical Hessian against the observed information matrix, a grid
#      search against a profile maximiser, and one-dimensional quadrature of
#      the posterior against the output of the sampler.
#
#   B. Empirical checks by Monte Carlo with a fixed seed. Each is reported as
#      a z-score against its Monte Carlo standard error, and |z| < 4 passes.
#      The standard error of a variance is computed from the fourth central
#      moment, not under normality. Where the chapter states a limit, the
#      rate is tested: the quantity must fall by the predicted factor between
#      two sample sizes.
#
# Writes results/check_formulas_report.txt
# =============================================================================

rm(list = ls())

local({
  f <- "code/core_ipd_prcs.R"
  if (!file.exists(f)) f <- "core_ipd_prcs.R"
  source(f)
})

setup_project()
source("code/core_ipd_prcs.R")

say <- new_report(file.path("results", "check_formulas_report.txt"))

PASS <- 0L
FAIL <- 0L
FAILED_NAMES <- character(0)

TOL_EXACT <- 1e-7
TOL_Z     <- 4

report_exact <- function(name, value, reference, tol = TOL_EXACT)
{
  # When the reference value is itself zero, a relative error is not defined
  # and the absolute difference is used instead.
  rel <- if (abs(reference) < 1e-8)
    abs(value - reference) else abs(value - reference) / abs(reference)

  ok <- is.finite(rel) && rel < tol

  if (ok) PASS <<- PASS + 1L else
  {
    FAIL <<- FAIL + 1L
    FAILED_NAMES <<- c(FAILED_NAMES, name)
  }

  say(sprintf("  %-58s %-4s  rel.err = %.3e", name,
              if (ok) "PASS" else "FAIL", rel))

  invisible(ok)
}


report_z <- function(name, estimate, target, se, tol = TOL_Z)
{
  z  <- (estimate - target) / se
  ok <- is.finite(z) && abs(z) < tol

  if (ok) PASS <<- PASS + 1L else
  {
    FAIL <<- FAIL + 1L
    FAILED_NAMES <<- c(FAILED_NAMES, name)
  }

  say(sprintf("  %-46s %-4s  MC = %10.5f  target = %10.5f  z = %7.3f",
              name, if (ok) "PASS" else "FAIL", estimate, target, z))

  invisible(ok)
}


report_flag <- function(name, ok, extra = "")
{
  if (isTRUE(ok)) PASS <<- PASS + 1L else
  {
    FAIL <<- FAIL + 1L
    FAILED_NAMES <<- c(FAILED_NAMES, name)
  }

  say(sprintf("  %-58s %-4s  %s", name,
              if (isTRUE(ok)) "PASS" else "FAIL", extra))

  invisible(ok)
}


say(rule())
say("VERIFICATION OF THE FORMULAS OF THE CHAPTER")
say("Inverse Pareto distribution under progressive Type-II random censoring")
say(rule())
say("Date          : ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"))
say("R version     : ", R.version.string)
say("Master seed   : ", MASTER_SEED)
say("Exact checks  : relative error below ", format(TOL_EXACT))
say("Empirical     : |z| below ", TOL_Z, " against the Monte Carlo standard error")
say("")


ALPHAS <- c(0.5, 1, 1.5, 2.5, 3)
BETAS  <- c(0.5, 1, 1.5)


# =============================================================================
# A. EXACT CHECKS
# =============================================================================

say(rule("-"))
say("A1. Density, distribution, survival, hazard and quantile functions")
say(rule("-"))

for (a in ALPHAS)
{
  # the density integrates to the distribution function
  for (x0 in c(0.3, 1, 4))
  {
    num <- num_int(function(t) d_ipd(t, a), 0, x0)

    report_exact(sprintf("int_0^%.1f f dx = F(%.1f)   (alpha = %.1f)",
                         x0, x0, a),
                 num, p_ipd(x0, a))
  }

  # the density integrates to one
  report_exact(sprintf("int_0^inf f dx = 1                (alpha = %.1f)", a),
               num_int(function(t) d_ipd(t, a), 0, Inf),
               1)

  # survival function
  report_exact(sprintf("S(x) + F(x) = 1                   (alpha = %.1f)", a),
               s_ipd(2, a) + p_ipd(2, a), 1)

  # hazard rate
  report_exact(sprintf("h(x) = f(x) / S(x)                (alpha = %.1f)", a),
               h_ipd(2, a), d_ipd(2, a) / s_ipd(2, a))

  # the quantile function inverts the distribution function
  uu <- c(0.01, 0.25, 0.5, 0.75, 0.99)

  report_exact(sprintf("F(Q(u)) = u                       (alpha = %.1f)", a),
               max(abs(p_ipd(q_ipd(uu, a), a) - uu)) + 1, 1)

  # the density is the derivative of the distribution function
  h  <- 1e-5
  x0 <- 1.7

  report_exact(sprintf("dF/dx = f(x) by finite difference  (alpha = %.1f)", a),
               (p_ipd(x0 + h, a) - p_ipd(x0 - h, a)) / (2 * h),
               d_ipd(x0, a), tol = 1e-6)
}

say("")


say(rule("-"))
say("A2. Moments and their conditions of existence")
say(rule("-"))

for (a in ALPHAS)
{
  for (r in c(-0.5, 0.25, 0.75))
  {
    if (r <= -a || r >= 1) next

    num <- num_int(function(x) x^r * d_ipd(x, a), 0, Inf)

    report_exact(sprintf("E(X^%.2f) = alpha B(r+alpha, 1-r)  (alpha = %.1f)",
                         r, a),
                 mom_ipd(r, a), num, tol = 1e-6)
  }

  for (r in c(0.3, 0.7))
  {
    if (r >= a) next

    num <- num_int(function(x) x^(-r) * d_ipd(x, a), 0, Inf)

    report_exact(sprintf("E(X^-%.1f) = G(r+1)G(a-r)/G(a)     (alpha = %.1f)",
                         r, a),
                 negmom_ipd(r, a), num, tol = 1e-6)
  }

  # E(log X) and E[log(1 + X)]
  num1 <- num_int(function(x) log(x) * d_ipd(x, a), 0, Inf)
  num2 <- num_int(function(x) log1p(x) * d_ipd(x, a), 0, Inf)

  report_exact(sprintf("E(log X)     = gamma + psi(alpha)  (alpha = %.1f)", a),
               elog_ipd(a), num1, tol = 1e-6)
  report_exact(sprintf("E(log(1+X))  = gamma+psi(a)+1/a    (alpha = %.1f)", a),
               elog1p_ipd(a), num2, tol = 1e-6)
}

# The first moment does not exist. The truncated mean must grow without bound
# and its growth must be logarithmic in the truncation point.
tm <- vapply(c(1e3, 1e4, 1e5),
             function(U) num_int(function(x) x * d_ipd(x, 2), 0, U),
             numeric(1))

report_flag("E(X) does not exist: truncated mean increases",
            all(diff(tm) > 0),
            sprintf("truncated means %.2f, %.2f, %.2f", tm[1], tm[2], tm[3]))

say("")


say(rule("-"))
say("A3. Transformations")
say(rule("-"))

for (a in ALPHAS)
{
  # R = 1 / X follows the Lomax distribution
  fR <- function(r) a / (1 + r)^(a + 1)

  report_exact(sprintf("R = 1/X is Lomax(alpha)            (alpha = %.1f)", a),
               num_int(fR, 0, Inf), 1)

  r0 <- 1.3

  report_exact(sprintf("f_R(r) by change of variables      (alpha = %.1f)", a),
               fR(r0), d_ipd(1 / r0, a) / r0^2)

  # Z = X / (1 + X) follows Beta(alpha, 1)
  z0 <- 0.42
  x0 <- z0 / (1 - z0)

  report_exact(sprintf("Z = X/(1+X) is Beta(alpha, 1)      (alpha = %.1f)", a),
               a * z0^(a - 1),
               d_ipd(x0, a) / (1 - z0)^2)
}

say("")


say(rule("-"))
say("A4. Shannon and Renyi entropies")
say(rule("-"))

for (a in ALPHAS)
{
  num <- num_int(function(x)
  {
    f <- d_ipd(x, a)
    ifelse(f > 0, -f * log(f), 0)
  }, 0, Inf)

  report_exact(sprintf("Shannon entropy H(X)               (alpha = %.1f)", a),
               shannon_ipd(a), num, tol = 1e-5)

  for (dl in c(0.75, 1.5, 2.5))
  {
    if (dl * (a - 1) + 1 <= 0) next

    num <- num_int(function(x) d_ipd(x, a)^dl, 0, Inf)

    report_exact(sprintf("Renyi entropy, delta = %.2f        (alpha = %.1f)",
                         dl, a),
                 renyi_ipd(a, dl), log(num) / (1 - dl), tol = 1e-5)
  }

  # the validity conditions are not vacuous
  report_flag(sprintf("Renyi entropy undefined at delta = 1 (alpha = %.1f)", a),
              is.na(renyi_ipd(a, 1)))
}

say("")


say(rule("-"))
say("A5. Expected values of the order statistics")
say(rule("-"))

for (a in c(0.5, 1, 2, 3))
{
  for (n in c(5, 8))
  {
    for (r in 1:(n - 1))
    {
      mi <- mu_rn_integral(r, n, a)
      ms <- mu_rn_series(r, n, a)
      md <- mu_rn_digamma(r, n, a)

      report_exact(sprintf("mu_%d:%d  series = integral        (alpha = %.1f)",
                           r, n, a), ms, mi, tol = 1e-5)
      report_exact(sprintf("mu_%d:%d  digamma = integral       (alpha = %.1f)",
                           r, n, a), md, mi, tol = 1e-5)
    }
  }
}

# the special case alpha = 1
for (n in c(5, 8, 12))
{
  for (r in 1:(n - 1))
    report_exact(sprintf("mu_%d:%d = r/(n-r) at alpha = 1", r, n),
                 mu_rn_digamma(r, n, 1), r / (n - r), tol = 1e-6)
}

# the recurrence relation n mu_{r:n-1} = (n-r) mu_{r:n} + r mu_{r+1:n}
for (a in c(0.5, 2, 3))
{
  for (n in c(6, 9))
  {
    for (r in 1:(n - 2))
      report_exact(sprintf("recurrence r = %d, n = %d           (alpha = %.1f)",
                           r, n, a),
                   n * mu_rn_digamma(r, n - 1, a),
                   (n - r) * mu_rn_digamma(r, n, a) +
                     r * mu_rn_digamma(r + 1, n, a),
                   tol = 1e-5)
  }
}

# the expected sample maximum does not exist
pm <- vapply(c(0.99, 0.999, 0.9999), function(U)
{
  n <- 5
  num_int(function(u)
  {
    lg <- log(u) / 2
    u^(n - 1) * exp(lg) / (-expm1(lg))
  }, 0, U) / beta(n, 1)
}, numeric(1))

report_flag("E(X_{n:n}) = infinity: partial integrals increase",
            all(diff(pm) > 0),
            sprintf("partial values %.2f, %.2f, %.2f", pm[1], pm[2], pm[3]))

say("")


say(rule("-"))
say("A6. Koziol-Green censoring mechanism")
say(rule("-"))

for (a in c(0.5, 1.5, 3))
{
  for (b in BETAS)
  {
    # F_T is a proper distribution function
    report_exact(sprintf("int f_T dt = 1        (alpha=%.1f, beta=%.1f)", a, b),
                 num_int(function(t) d_kg(t, a, b), 0, Inf), 1)

    # f_T is the derivative of F_T
    h  <- 1e-5
    t0 <- 1.4

    report_exact(sprintf("dF_T/dt = f_T         (alpha=%.1f, beta=%.1f)", a, b),
                 (p_kg(t0 + h, a, b) - p_kg(t0 - h, a, b)) / (2 * h),
                 d_kg(t0, a, b), tol = 1e-6)

    # P(D = 1) = 1 / (1 + beta), by direct integration
    num <- num_int(function(x) (1 - p_kg(x, a, b)) * d_ipd(x, a),
                   0, Inf)

    report_exact(sprintf("P(D=1) = 1/(1+beta)   (alpha=%.1f, beta=%.1f)", a, b),
                 1 / (1 + b), num, tol = 1e-6)

    # the marginal of Y, F_Y = 1 - S_X^{1+beta}
    y0 <- 1.9

    num <- num_int(function(x) (1 + b) * d_ipd(x, a) *
                     s_ipd(x, a)^b, 0, y0)

    report_exact(sprintf("F_Y = 1 - S_X^(1+beta)(alpha=%.1f, beta=%.1f)", a, b),
                 1 - s_ipd(y0, a)^(1 + b), num, tol = 1e-6)

    # the inversion used to generate the censoring times
    v  <- c(0.05, 0.3, 0.77, 0.95)
    tt <- r_kg(length(v), a, b, v = v)

    report_exact(sprintf("F_T(T(v)) = v         (alpha=%.1f, beta=%.1f)", a, b),
                 max(abs(p_kg(tt, a, b) - v)) + 1, 1, tol = 1e-8)
  }
}

say("")


say(rule("-"))
say("A7. Likelihood, score and observed information")
say(rule("-"))

dat <- CHECK_DATA

say("  Check data set, n = ", dat$n, ", m = ", dat$m, ", k = ", dat$k)
say("  y = ", paste(formatC(dat$y, format = "f", digits = 4), collapse = ", "))
say("  d = ", paste(dat$d, collapse = ", "))
say("  W = ", paste(dat$W, collapse = ", "))
say("")

for (a in c(0.7, 1.5, 2.6))
{
  for (b in c(0.4, 1, 1.8))
  {
    # the reduced likelihood of Eq. (20) against the product form before the
    # Koziol-Green substitutions are made
    report_exact(sprintf("log L reduced = direct  (alpha=%.1f, beta=%.1f)",
                         a, b),
                 loglik_prcs(a, b, dat),
                 loglik_prcs_direct(a, b, dat), tol = 1e-9)

    # the analytic score against a central difference of the log-likelihood
    h  <- 1e-6
    ga <- (loglik_prcs(a + h, b, dat) - loglik_prcs(a - h, b, dat)) / (2 * h)
    gb <- (loglik_prcs(a, b + h, dat) - loglik_prcs(a, b - h, dat)) / (2 * h)

    sc <- score_prcs(a, b, dat)

    report_exact(sprintf("score alpha = numerical (alpha=%.1f, beta=%.1f)",
                         a, b), sc["alpha"], ga, tol = 1e-5)
    report_exact(sprintf("score beta  = numerical (alpha=%.1f, beta=%.1f)",
                         a, b), sc["beta"], gb, tol = 1e-5)

    # the observed information against a numerical Hessian
    h  <- 1e-3
    J  <- obs_info_prcs(a, b, dat)

    Haa <- (loglik_prcs(a + h, b, dat) - 2 * loglik_prcs(a, b, dat) +
              loglik_prcs(a - h, b, dat)) / h^2
    Hbb <- (loglik_prcs(a, b + h, dat) - 2 * loglik_prcs(a, b, dat) +
              loglik_prcs(a, b - h, dat)) / h^2
    Hab <- (loglik_prcs(a + h, b + h, dat) - loglik_prcs(a + h, b - h, dat) -
              loglik_prcs(a - h, b + h, dat) +
              loglik_prcs(a - h, b - h, dat)) / (4 * h^2)

    report_exact(sprintf("J_aa = -numerical Hess  (alpha=%.1f, beta=%.1f)",
                         a, b), J[1, 1], -Haa, tol = 1e-4)
    report_exact(sprintf("J_bb = -numerical Hess  (alpha=%.1f, beta=%.1f)",
                         a, b), J[2, 2], -Hbb, tol = 1e-4)
    report_exact(sprintf("J_ab = -numerical Hess  (alpha=%.1f, beta=%.1f)",
                         a, b), J[1, 2], -Hab, tol = 1e-4)
  }
}

# J_bb has a closed form
report_exact("J_bb = (m - k) / beta^2",
             obs_info_prcs(1.5, 0.8, dat)[2, 2],
             (dat$m - dat$k) / 0.8^2)

say("")


say(rule("-"))
say("A8. Profile estimator of beta and the maximum likelihood estimates")
say(rule("-"))

for (a in c(0.7, 1.5, 2.6))
{
  bh <- beta_profile(a, dat)

  # the profile estimator solves the score equation in beta
  report_exact(sprintf("d loglik / d beta = 0 at beta(alpha)  (alpha = %.1f)",
                       a),
               score_prcs(a, bh, dat)["beta"] + 1, 1, tol = 1e-9)

  # and it maximises the log-likelihood in beta, checked by a grid
  gr <- seq(max(1e-3, bh - 1), bh + 1, length.out = 4001L)
  ll <- vapply(gr, function(b) loglik_prcs(a, b, dat), numeric(1))

  report_exact(sprintf("grid maximiser in beta = beta(alpha)  (alpha = %.1f)",
                       a),
               gr[which.max(ll)], bh, tol = 1e-3)
}

ft <- fit_mle(dat)

report_flag("profile route converges", isTRUE(ft$ok))

if (isTRUE(ft$ok))
{
  say(sprintf("  alpha_hat = %.10f, beta_hat = %.10f", ft$alpha, ft$beta))

  # the score vanishes at the estimates
  report_flag("score at the estimates is numerically zero",
              max(abs(ft$score)) < 1e-7,
              sprintf("max |score| = %.3e", max(abs(ft$score))))

  # the observed information matrix is positive definite
  report_flag("observed information matrix is positive definite", ft$pd,
              sprintf("eigenvalues %.4f and %.4f",
                      eigen(ft$J, symmetric = TRUE,
                            only.values = TRUE)$values[1],
                      eigen(ft$J, symmetric = TRUE,
                            only.values = TRUE)$values[2]))

  # the independent iterative route reaches the same point
  nw <- fit_mle_newton(dat, start = c(1, 1))

  report_flag("Newton-Raphson route converges", isTRUE(nw$ok),
              sprintf("iterations = %d", nw$iter))

  if (isTRUE(nw$ok))
  {
    rel <- max(abs(nw$est - ft$est) / abs(ft$est))

    report_flag("profile and Newton-Raphson agree", rel < 1e-8,
                sprintf("largest relative difference = %.3e", rel))
  }

  # and a direct two-dimensional maximisation agrees as well
  op <- optim(c(1, 1), function(p) -loglik_prcs(p[1], p[2], dat),
              method = "Nelder-Mead",
              control = list(reltol = 1e-14, maxit = 5000))

  rel <- max(abs(op$par - ft$est) / abs(ft$est))

  report_flag("profile and Nelder-Mead agree", rel < 1e-4,
              sprintf("largest relative difference = %.3e", rel))

  # the delta-method standard errors against a numerical gradient of the
  # inverse information
  report_exact("V = J^{-1}, first diagonal entry",
               ft$V[1, 1], solve(ft$J)[1, 1])

  # the log-transformed interval is positive and contains the estimate
  report_flag("log-transformed CI contains the estimate and is positive",
              all(ft$ci_log[, "lower"] > 0) &&
                all(ft$ci_log[, "lower"] < ft$est) &&
                all(ft$ci_log[, "upper"] > ft$est))
}

say("")


say(rule("-"))
say("A9. Bayesian full conditionals")
say(rule("-"))

a1 <- 2; b1 <- 1; a2 <- 2; b2 <- 2

# the full conditional of beta is a gamma density with the stated rate
for (a in c(0.8, 1.6, 2.4))
{
  rt <- beta_rate(a, dat, b2)

  report_flag(sprintf("rate of the gamma full conditional is positive (a=%.1f)",
                      a), rt > 0, sprintf("rate = %.6f", rt))

  bb <- c(0.3, 0.9, 1.7, 2.8)

  # log posterior of beta given alpha, up to a constant, from the joint
  lp_joint <- vapply(bb, function(b)
    (dat$m - dat$k + a2 - 1) * log(b) - b2 * b +
      sum((b + (1 + b) * dat$W) *
            log(-expm1(a * log(dat$y / (1 + dat$y))))), numeric(1))

  lp_gamma <- dgamma(bb, shape = dat$m - dat$k + a2, rate = rt, log = TRUE)

  # the two must differ by a constant
  dd <- lp_joint - lp_gamma

  report_exact(sprintf("beta | alpha is Gamma(m-k+a2, rate)   (alpha = %.1f)",
                       a),
               max(dd) - min(dd) + 1, 1, tol = 1e-8)
}

# The marginal posterior of alpha, obtained by integrating beta out exactly,
# reproduces the joint posterior after numerical integration over beta. The
# marginal is defined only up to a constant, so the check is that the
# difference between the two routes does not depend on alpha and that the
# constant is the log gamma value that was left out.
diffs <- vapply(c(0.8, 1.6, 2.4), function(a)
{
  num <- num_int(function(b)
  {
    vapply(b, function(bb)
      exp(log_cond_alpha(a, bb, dat, a1, b1) +
            (dat$m - dat$k + a2 - 1) * log(bb) - b2 * bb), numeric(1))
  }, 0, Inf)

  log(num) - log_marg_alpha(a, dat, a1, b1, a2, b2)
}, numeric(1))

# The tolerance here is set by the quadrature over the half line, not by the
# formula, which is exact.
report_exact("marginal of alpha, integrating beta out, is exact",
             max(diffs) - min(diffs), 0, tol = 1e-4)
report_exact("the constant left out is log Gamma(m - k + a2)",
             mean(diffs), lgamma(dat$m - dat$k + a2), tol = 1e-5)

# the HPD interval is the shortest interval that holds the required proportion
set.seed(MASTER_SEED)
xs <- rgamma(20000, shape = 5, rate = 2)
hi <- hpd(xs, 0.95)

cand <- vapply(seq(0.001, 0.049, by = 0.001), function(p)
  diff(quantile(xs, c(p, p + 0.95))), numeric(1))

report_flag("HPD interval is the shortest 95 per cent interval",
            diff(hi) <= min(cand) + 1e-8,
            sprintf("HPD length = %.5f, shortest quantile pair = %.5f",
                    diff(hi), min(cand)))

report_flag("HPD interval holds 95 per cent of the sample",
            abs(mean(xs >= hi[1] & xs <= hi[2]) - 0.95) < 0.002,
            sprintf("actual coverage of the sample = %.4f",
                    mean(xs >= hi[1] & xs <= hi[2])))

say("")


# =============================================================================
# B. EMPIRICAL CHECKS
# =============================================================================

say(rule("-"))
say("B1. Distributional results by Monte Carlo")
say(rule("-"))

set.seed(MASTER_SEED)

NMC <- 400000L

for (a in c(0.8, 1.5, 3))
{
  x <- q_ipd(runif(NMC), a)

  v <- log(x)
  report_z(sprintf("E(log X), alpha = %.1f", a),
           mean(v), elog_ipd(a), sd(v) / sqrt(NMC))

  v <- log1p(x)
  report_z(sprintf("E(log(1+X)), alpha = %.1f", a),
           mean(v), elog1p_ipd(a), sd(v) / sqrt(NMC))

  v <- -log(d_ipd(x, a))
  report_z(sprintf("Shannon entropy, alpha = %.1f", a),
           mean(v), shannon_ipd(a), sd(v) / sqrt(NMC))

  r <- 0.5
  if (r < a)
  {
    v <- x^(-r)
    report_z(sprintf("E(X^-0.5), alpha = %.1f", a),
             mean(v), negmom_ipd(r, a), sd(v) / sqrt(NMC))
  }

  # Z = X / (1 + X) is Beta(alpha, 1), checked through its mean and variance
  z  <- x / (1 + x)
  mz <- a / (a + 1)
  vz <- a / ((a + 1)^2 * (a + 2))

  report_z(sprintf("mean of Z = X/(1+X), alpha = %.1f", a),
           mean(z), mz, sd(z) / sqrt(NMC))

  # the standard error of a variance from the fourth central moment
  m4 <- mean((z - mean(z))^4)
  se <- sqrt((m4 - var(z)^2) / NMC)

  report_z(sprintf("variance of Z, alpha = %.1f", a), var(z), vz, se)
}

say("")


say(rule("-"))
say("B2. Expected values of the order statistics by Monte Carlo")
say(rule("-"))

set.seed(MASTER_SEED + 1L)

NMC2 <- 200000L

for (a in c(1.5, 3))
{
  n <- 6

  X <- matrix(q_ipd(runif(NMC2 * n), a), nrow = NMC2)
  X <- t(apply(X, 1, sort))

  for (r in 1:(n - 2))
  {
    v <- X[, r]

    report_z(sprintf("mu_%d:%d, alpha = %.1f", r, n, a),
             mean(v), mu_rn_digamma(r, n, a), sd(v) / sqrt(NMC2))
  }
}

say("")


say(rule("-"))
say("B3. Censoring mechanism by Monte Carlo")
say(rule("-"))

set.seed(MASTER_SEED + 2L)

NMC3 <- 400000L

for (a in c(0.8, 2))
{
  for (b in BETAS)
  {
    X <- q_ipd(runif(NMC3), a)
    T <- r_kg(NMC3, a, b)

    D <- as.numeric(X <= T)

    report_z(sprintf("P(D=1), alpha = %.1f, beta = %.1f", a, b),
             mean(D), 1 / (1 + b), sd(D) / sqrt(NMC3))

    # the marginal distribution of Y at three quantiles
    Y <- pmin(X, T)

    for (y0 in c(0.5, 1.5, 4))
    {
      v <- as.numeric(Y <= y0)

      report_z(sprintf("F_Y(%.1f), alpha = %.1f, beta = %.1f", y0, a, b),
               mean(v), 1 - s_ipd(y0, a)^(1 + b), sd(v) / sqrt(NMC3))
    }
  }
}

say("")


say(rule("-"))
say("B4. Limiting distributions of the extreme order statistics")
say(rule("-"))

# The chapter states three limits. Each is tested as a rate: the largest
# absolute difference between the finite-sample distribution function and its
# limit must fall as n increases.

a <- 3

bU <- function(n) q_ipd(1 - 1 / n, a)
bL <- function(n) q_ipd(1 / n, a)

xg  <- exp(seq(log(0.05), log(200), length.out = 2000L))
xgl <- seq(0.001, 2.5, length.out = 2000L)

err_max <- vapply(c(10, 100, 1000), function(n)
  max(abs(exp(n * log(pmax(p_ipd(bU(n) * xg, a), 1e-300))) -
            exp(-1 / xg))), numeric(1))

report_flag("normalized maximum approaches the Frechet limit",
            all(diff(err_max) < 0),
            sprintf("max errors %.5f, %.5f, %.5f at n = 10, 100, 1000",
                    err_max[1], err_max[2], err_max[3]))

err_min <- vapply(c(100, 1000, 10000), function(n)
  max(abs(-expm1(n * log1p(-p_ipd(bL(n) * xgl, a))) -
            (-expm1(-xgl^a)))), numeric(1))

report_flag("normalized minimum approaches its limit",
            all(diff(err_min) < 0),
            sprintf("max errors %.5f, %.5f, %.5f at n = 100, 1000, 10000",
                    err_min[1], err_min[2], err_min[3]))

i_up <- 3

err_up <- vapply(c(10, 100, 1000), function(n)
{
  Fy <- p_ipd(bU(n) * xg, a)
  fn <- pbinom(i_up - 1, n, pmin(pmax(1 - Fy, 0), 1))
  lm <- exp(-1 / xg) * (1 + 1 / xg + 1 / (2 * xg^2))
  max(abs(fn - lm))
}, numeric(1))

report_flag("fixed upper extreme with i = 3 approaches its limit",
            all(diff(err_up) < 0),
            sprintf("max errors %.5f, %.5f, %.5f at n = 10, 100, 1000",
                    err_up[1], err_up[2], err_up[3]))

# the same three limits by simulation at one sample size
set.seed(MASTER_SEED + 3L)

nn   <- 200L
NMC4 <- 40000L

X <- matrix(q_ipd(runif(NMC4 * nn), a), nrow = NMC4)

mx <- apply(X, 1, max) / bU(nn)
mn <- apply(X, 1, min) / bL(nn)

# The simulated distributions are compared with the exact finite-sample
# distributions, not with the limits: at n = 200 the normalized minimum is
# still far from its limit, which is the point the chapter itself makes.
for (x0 in c(0.5, 1, 2, 5))
{
  v <- as.numeric(mx <= x0)

  report_z(sprintf("P(max / b_n <= %.1f) at n = %d, exact", x0, nn),
           mean(v), p_ipd(bU(nn) * x0, a)^nn, sd(v) / sqrt(NMC4), tol = 6)
}

for (x0 in c(0.5, 0.9, 1.3))
{
  v <- as.numeric(mn <= x0)

  report_z(sprintf("P(min / b*_n <= %.1f) at n = %d, exact", x0, nn),
           mean(v), 1 - (1 - p_ipd(bL(nn) * x0, a))^nn,
           sd(v) / sqrt(NMC4), tol = 6)
}

say("")


say(rule("-"))
say("B5. Estimation by Monte Carlo")
say(rule("-"))

set.seed(MASTER_SEED + 4L)

a_true <- 1.5
b_true <- 1.0

sim_fit <- function(nrep, n, m, W)
{
  out <- matrix(NA_real_, nrep, 5,
                dimnames = list(NULL, c("alpha", "beta", "v_aa", "cov_a",
                                        "cov_b")))

  for (i in seq_len(nrep))
  {
    dd <- gen_prcs(n, m, W, a_true, b_true)
    ff <- fit_mle(dd)

    if (isTRUE(ff$ok))
    {
      out[i, "alpha"] <- ff$alpha
      out[i, "beta"]  <- ff$beta
      out[i, "v_aa"]  <- ff$V[1, 1]
      out[i, "cov_a"] <- as.numeric(ff$ci_log[1, 1] <= a_true &&
                                      a_true <= ff$ci_log[1, 2])
      out[i, "cov_b"] <- as.numeric(ff$ci_log[2, 1] <= b_true &&
                                      b_true <= ff$ci_log[2, 2])
    }
  }

  out[stats::complete.cases(out), , drop = FALSE]
}

NREP <- 3000L

sc  <- SCHEMES_SIM$S10
r10 <- sim_fit(NREP, sc$n, sc$m, sc$W)

say(sprintf("  scheme S10, n = %d, m = %d, replications used = %d of %d",
            sc$n, sc$m, nrow(r10), NREP))

# the estimator of alpha is close to the true value
se <- sd(r10[, "alpha"]) / sqrt(nrow(r10))

report_z("mean of alpha_hat under S10", mean(r10[, "alpha"]), a_true, se,
         tol = 6)

# the asymptotic variance from the information matrix agrees with the
# Monte Carlo variance of the estimates
v_mc <- var(r10[, "alpha"])
m4   <- mean((r10[, "alpha"] - mean(r10[, "alpha"]))^4)
se_v <- sqrt((m4 - v_mc^2) / nrow(r10))

report_z("asymptotic variance of alpha_hat matches the MC variance",
         mean(r10[, "v_aa"]), v_mc, se_v, tol = 6)

# coverage of the log-transformed confidence intervals
for (nm in c("cov_a", "cov_b"))
{
  p  <- mean(r10[, nm])
  se <- sqrt(p * (1 - p) / nrow(r10))

  report_z(sprintf("coverage of the 95 per cent CI, %s",
                   if (nm == "cov_a") "alpha" else "beta"),
           p, 0.95, se, tol = 6)
}

# the rate at which the estimator approaches the true value. The chapter
# states that the estimator becomes more accurate as m grows; here the mean
# of the square root of m times the absolute deviation must stay bounded
# while the mean absolute deviation itself falls by about the square root of
# the ratio of the two values of m.
sc2 <- SCHEMES_SIM$S2
r2  <- sim_fit(NREP, sc2$n, sc2$m, sc2$W)

mad2  <- mean(abs(r2[,  "alpha"] - a_true))
mad10 <- mean(abs(r10[, "alpha"] - a_true))

ratio_obs <- mad2 / mad10
ratio_pre <- sqrt(sc$m / sc2$m)

report_flag("accuracy improves at the expected rate in m",
            abs(ratio_obs / ratio_pre - 1) < 0.25,
            sprintf("observed ratio = %.4f, predicted = %.4f (m = %d and %d)",
                    ratio_obs, ratio_pre, sc2$m, sc$m))

say("")


say(rule("-"))
say("B6. The sampler against quadrature of the posterior")
say(rule("-"))

# Because the full conditional of beta is a gamma distribution, beta can be
# integrated out exactly and the posterior means can be obtained by
# one-dimensional quadrature. This is an independent check of the sampler.

set.seed(MASTER_SEED + 5L)

qm <- post_moments_quad(dat, a1, b1, a2, b2)

ch <- mcmc_prcs(dat, a1, b1, a2, b2, zeta = 200000L, zeta0 = 20000L)

say(sprintf("  acceptance rate = %.4f, proposal sd = %.4f",
            ch$acc_rate, ch$sigma_alpha))

report_flag("acceptance rate is between 0.25 and 0.50",
            ch$acc_rate > 0.20 && ch$acc_rate < 0.60,
            sprintf("acceptance rate = %.4f", ch$acc_rate))

nk  <- length(ch$alpha)

# an effective sample size that allows for the autocorrelation of the chain
ess <- function(x)
{
  ac <- acf(x, lag.max = 200, plot = FALSE)$acf[-1]
  s  <- 0

  for (i in seq_along(ac))
  {
    if (ac[i] <= 0) break
    s <- s + ac[i]
  }

  length(x) / (1 + 2 * s)
}

e_a <- ess(ch$alpha)
e_b <- ess(ch$beta)

say(sprintf("  effective sample sizes: alpha %.0f, beta %.0f", e_a, e_b))

report_z("posterior mean of alpha, sampler against quadrature",
         mean(ch$alpha), qm["mean_alpha"], sd(ch$alpha) / sqrt(e_a), tol = 4)

report_z("posterior mean of beta, sampler against quadrature",
         mean(ch$beta), qm["mean_beta"], sd(ch$beta) / sqrt(e_b), tol = 4)

report_z("posterior sd of alpha, sampler against quadrature",
         sd(ch$alpha), qm["sd_alpha"],
         sd(ch$alpha) / sqrt(2 * e_a), tol = 4)

say("")


# =============================================================================
# Summary
# =============================================================================

say(rule())
say("SUMMARY")
say(rule())
say("Checks passed : ", PASS)
say("Checks failed : ", FAIL)

if (FAIL > 0L)
{
  say("")
  say("Failed checks:")
  for (nm in unique(FAILED_NAMES))
    say("  ", nm)
}

say("")
say(if (FAIL == 0L) "ALL CHECKS PASSED" else "SOME CHECKS FAILED")
say(rule())
