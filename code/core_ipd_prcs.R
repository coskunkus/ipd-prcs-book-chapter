# =============================================================================
# core_ipd_prcs.R
#
# Shared core for the book chapter
#
#   "Statistical Inference for the Inverse Pareto Distribution under
#    Progressive Type-II Random Censoring: Applications to Sustainable
#    Reliability"
#   Indrajeet Kumar, Kapil Kumar, Ashish Raj, Coskun Kus
#
# Every other script sources this file, so the applications, the simulation
# and the verification scripts all use exactly one implementation.
#
# Base R only. No contributed package is required.
#
# Notation follows the chapter:
#   alpha  IPD shape parameter
#   beta   Koziol-Green censoring parameter
#   n      number of units placed on test
#   m      number of observed event times
#   k      number of observed failures, k = sum(d_i)
#   W      progressive withdrawal scheme, sum(W) + m = n
#   y, d   observed event times and event indicators
#   z_i    y_i / (1 + y_i)
# =============================================================================


# -----------------------------------------------------------------------------
# Project root, so that a script runs when it is sourced, when the Source
# button is pressed in RStudio, and when it is started with Rscript.
# -----------------------------------------------------------------------------

find_project_root <- function(core_name = "core_ipd_prcs.R")
{
  cand <- character(0)

  # 1. the file being sourced
  nf <- sys.nframe()
  if (nf > 0)
  {
    for (i in seq_len(nf))
    {
      of <- sys.frame(i)$ofile
      if (!is.null(of) && nzchar(of))
        cand <- c(cand, dirname(normalizePath(of, mustWork = FALSE)))
    }
  }

  # 2. the file open in the RStudio editor
  if (requireNamespace("rstudioapi", quietly = TRUE))
  {
    p <- try(rstudioapi::getSourceEditorContext()$path, silent = TRUE)

    if (!inherits(p, "try-error") && length(p) == 1L && nzchar(p))
      cand <- c(cand, dirname(normalizePath(p, mustWork = FALSE)))
  }

  # 3. the --file= argument of Rscript
  a <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)

  if (length(a))
    cand <- c(cand,
              dirname(normalizePath(sub("^--file=", "", a[1L]),
                                    mustWork = FALSE)))

  # 4. the working directory
  cand <- c(cand, normalizePath(getwd(), mustWork = FALSE))

  for (d in unique(cand))
  {
    dd <- d

    for (up in 1:3)
    {
      if (file.exists(file.path(dd, "code", core_name)))
        return(dd)

      if (file.exists(file.path(dd, core_name)))
        return(dirname(dd))

      dd <- dirname(dd)
    }
  }

  stop("Project root not found. Set the working directory to the folder ",
       "that contains code/", core_name, ".")
}


setup_project <- function()
{
  root <- find_project_root()
  setwd(root)

  dir.create("results", showWarnings = FALSE)
  dir.create("figures", showWarnings = FALSE)

  root
}


# -----------------------------------------------------------------------------
# Global settings used by every script
# -----------------------------------------------------------------------------

MASTER_SEED <- 20260908

# Bias and MSE are printed both on the raw scale, which is the scale used in
# the chapter, and multiplied by var_scale in a second table.
var_scale     <- 1e3
var_scale_tex <- "10^{3}"

# Convergence rule for every iterative fit. Both conditions must hold.
EPS_LOGLIK <- 1e-12
EPS_PAR    <- 1e-9
MAX_ITER   <- 500


# -----------------------------------------------------------------------------
# Robust numerical integration
#
# QUADPACK, which is what integrate() calls, reports
#   "roundoff error is detected in the extrapolation table"
# when a tight tolerance is asked for on a half line whose integrand is
# singular at the lower endpoint. Several integrals used to verify the
# formulas of this chapter are of that kind, for example E(X^-r) with a small
# shape parameter, whose integrand behaves as x^{alpha-r-1} near zero.
#
# The half line is therefore split at one, so that the endpoint singularity
# and the tail are handled separately, and the tolerance is relaxed in steps
# until a finite value is obtained. A failure returns NA rather than an
# error, so that one difficult integral is reported as a failed check instead
# of stopping the whole script.
# -----------------------------------------------------------------------------

num_int <- function(f, lower, upper, rel.tol = 1e-10, subdivisions = 5000L)
{
  one <- function(a, b, tol)
    tryCatch(integrate(f, a, b, rel.tol = tol,
                       subdivisions = subdivisions)$value,
             error = function(e) NA_real_)

  for (tol in unique(c(rel.tol, 1e-9, 1e-8, 1e-7, 1e-6, 1e-5)))
  {
    v <- if (is.finite(upper))
      one(lower, upper, tol)
    else
    {
      mid <- if (lower < 1) 1 else 2 * lower

      one(lower, mid, tol) + one(mid, upper, tol)
    }

    if (is.finite(v))
      return(v)
  }

  NA_real_
}


# -----------------------------------------------------------------------------
# A portable pseudo-random generator
#
# The applications of the chapter impose a censoring design on complete data,
# so their numbers depend on the stream of pseudo-random numbers used. A
# result that depends on the generator built into one version of one language
# is not reproducible elsewhere.
#
# The Lehmer generator below is therefore used for everything in run_all.R.
# It is fully specified by three integers, so the censored samples and the
# estimates reproduce to the last printed digit on any machine and under any
# version of R.
#
#   x_{k+1} = 16807 x_k  mod  2147483647,     u_k = x_k / 2147483647
#
# The multiplication stays below 2^53, so double precision is exact.
# The Monte Carlo study keeps R's own generator, because it does not need to
# be reproduced in a second language.
# -----------------------------------------------------------------------------

LCG_A <- 16807
LCG_M <- 2147483647

lcg_new <- function(seed)
{
  e <- new.env(parent = emptyenv())

  s <- seed %% LCG_M
  e$state <- if (s <= 0) s + LCG_M - 1 else s

  e
}


lcg_unif <- function(e, n = 1L)
{
  u <- numeric(n)

  for (i in seq_len(n))
  {
    e$state <- (LCG_A * e$state) %% LCG_M
    u[i] <- e$state / LCG_M
  }

  u
}


# Box-Muller, using the first of each pair
lcg_norm <- function(e, n = 1L)
{
  z <- numeric(n)

  for (i in seq_len(n))
  {
    u <- lcg_unif(e, 2L)
    z[i] <- sqrt(-2 * log(u[1])) * cos(2 * pi * u[2])
  }

  z
}


# Gamma with an integer shape, as a sum of exponentials. The shape needed
# here is m - k + a2, which is an integer whenever a2 is.
lcg_gamma_int <- function(e, shape, rate)
{
  if (abs(shape - round(shape)) > 1e-12)
    stop("lcg_gamma_int requires an integer shape.")

  -sum(log(lcg_unif(e, as.integer(round(shape))))) / rate
}


# -----------------------------------------------------------------------------
# Inverse Pareto distribution, Eqs. (1), (2) and (3) of the chapter
# -----------------------------------------------------------------------------

d_ipd <- function(x, alpha)
{
  out <- numeric(length(x))
  ok  <- x > 0

  out[ok] <- alpha * x[ok]^(alpha - 1) / (1 + x[ok])^(alpha + 1)

  out
}


p_ipd <- function(x, alpha)
{
  out <- numeric(length(x))
  ok  <- x > 0

  out[ok] <- (x[ok] / (1 + x[ok]))^alpha

  out
}


s_ipd <- function(x, alpha)
{
  1 - p_ipd(x, alpha)
}


h_ipd <- function(x, alpha)
{
  out <- numeric(length(x))
  ok  <- x > 0

  out[ok] <- alpha * x[ok]^(alpha - 1) /
    ((1 + x[ok]) * ((1 + x[ok])^alpha - x[ok]^alpha))

  out
}


# Quantile function, Eq. (3). Written with expm1 so that it stays accurate
# for u close to one.
q_ipd <- function(u, alpha)
{
  if (any(u <= 0 | u >= 1))
    stop("u must lie strictly between 0 and 1.")

  logq <- log(u) / alpha

  exp(logq) / (-expm1(logq))
}


# -----------------------------------------------------------------------------
# Moments and entropies
# -----------------------------------------------------------------------------

# E(X^r) = alpha * B(r + alpha, 1 - r),  valid for -alpha < r < 1
mom_ipd <- function(r, alpha)
{
  if (r <= -alpha || r >= 1)
    return(NA_real_)

  alpha * beta(r + alpha, 1 - r)
}


# E(X^{-r}) = Gamma(r + 1) Gamma(alpha - r) / Gamma(alpha),  0 < r < alpha
negmom_ipd <- function(r, alpha)
{
  if (r <= 0 || r >= alpha)
    return(NA_real_)

  exp(lgamma(r + 1) + lgamma(alpha - r) - lgamma(alpha))
}


# E(log X) = gamma + psi(alpha)
elog_ipd <- function(alpha)
{
  -digamma(1) + digamma(alpha)
}


# E[log(1 + X)] = gamma + psi(alpha) + 1 / alpha
elog1p_ipd <- function(alpha)
{
  -digamma(1) + digamma(alpha) + 1 / alpha
}


# Shannon entropy, H(X) = -log alpha + 2 [gamma + psi(alpha)] + (alpha + 1)/alpha
shannon_ipd <- function(alpha)
{
  -log(alpha) + 2 * (-digamma(1) + digamma(alpha)) + (alpha + 1) / alpha
}


# Renyi entropy of order delta
renyi_ipd <- function(alpha, delta)
{
  a <- delta * (alpha - 1) + 1
  b <- 2 * delta - 1

  if (delta <= 0.5 || a <= 0 || isTRUE(all.equal(delta, 1)))
    return(NA_real_)

  (delta * log(alpha) + lbeta(a, b)) / (1 - delta)
}


# -----------------------------------------------------------------------------
# Expected values of the order statistics, three representations
# -----------------------------------------------------------------------------

# Eq. (5), the beta integral representation
mu_rn_integral <- function(r, n, alpha)
{
  if (r >= n)
    return(Inf)

  f <- function(u)
  {
    lg <- log(u) / alpha
    u^(r - 1) * (1 - u)^(n - r) * exp(lg) / (-expm1(lg))
  }

  val <- num_int(f, 0, 1, rel.tol = 1e-12)

  val / beta(r, n - r + 1)
}


# Eq. (6), the positive beta series.
#
# The summand behaves as k^{-(n-r+1)}, so for r = n - 1 the series converges
# only as k^{-2} and the truncated sum is not accurate enough on its own. The
# tail beyond kmax is therefore added from the asymptotic form of the
# summand, B(a, b) ~ Gamma(b) a^{-b} for large a.
mu_rn_series <- function(r, n, alpha, kmax = 2e5)
{
  if (r >= n)
    return(Inf)

  k     <- 1:kmax
  lterm <- lbeta(r + k / alpha, n - r + 1) - lbeta(r, n - r + 1)

  s <- sum(exp(lterm))

  b    <- n - r + 1
  tail <- exp(lgamma(b) - lbeta(r, b)) * alpha / (b - 1) *
    (r + (kmax + 0.5) / alpha)^(-(b - 1))

  s + tail
}


# Eq. (7), the finite digamma representation
mu_rn_digamma <- function(r, n, alpha)
{
  if (r >= n)
    return(Inf)

  j <- 0:(n - r)

  s <- sum((-1)^j * choose(n - r, j) * digamma(1 + alpha * (r + j)))

  -alpha * exp(lgamma(n + 1) - lgamma(r) - lgamma(n - r + 1)) * s
}


# -----------------------------------------------------------------------------
# Koziol-Green random censoring
# -----------------------------------------------------------------------------

p_kg <- function(t, alpha, beta_par)
{
  1 - s_ipd(t, alpha)^beta_par
}


d_kg <- function(t, alpha, beta_par)
{
  beta_par * s_ipd(t, alpha)^(beta_par - 1) * d_ipd(t, alpha)
}


# Censoring times by inversion of the Koziol-Green cdf.
r_kg <- function(nn, alpha, beta_par, v = NULL)
{
  if (is.null(v))
    v <- runif(nn)

  q <- (1 - (1 - v)^(1 / beta_par))^(1 / alpha)

  q / (1 - q)
}


# -----------------------------------------------------------------------------
# Progressive Type-II random censoring
# -----------------------------------------------------------------------------

check_scheme <- function(n, m, W)
{
  if (length(W) != m)
    stop("W must have exactly m entries.")

  if (sum(W) + m != n)
    stop("sum(W) + m must equal n.")

  invisible(TRUE)
}


# Draw a PRCS sample. The unit with the smallest observed time among those
# still on test is recorded, after which W_i of the surviving units are
# withdrawn at random.
prcs_from_YD <- function(Y, D, m, W, gen = NULL)
{
  n     <- length(Y)
  alive <- seq_len(n)

  y <- numeric(m)
  d <- integer(m)

  for (i in seq_len(m))
  {
    j <- alive[which.min(Y[alive])]

    y[i] <- Y[j]
    d[i] <- D[j]

    alive <- alive[alive != j]

    if (W[i] > 0)
    {
      if (length(alive) < W[i])
        stop("Not enough surviving units to withdraw at stage ", i, ".")

      # With gen supplied the withdrawals are drawn from the portable
      # generator, one unit at a time, so that the units removed are fixed by
      # the seed alone. Without it R's own generator is used, which is what
      # the Monte Carlo study does.
      if (is.null(gen))
      {
        idx   <- sample.int(length(alive), W[i])
        alive <- alive[-idx]
      }
      else
      {
        for (jj in seq_len(W[i]))
        {
          u   <- lcg_unif(gen, 1L)
          idx <- floor(u * length(alive)) + 1

          alive <- alive[-idx]
        }
      }
    }
  }

  list(y = y, d = d, W = W, n = n, m = m, k = sum(d))
}


gen_prcs <- function(n, m, W, alpha, beta_par)
{
  check_scheme(n, m, W)

  X <- q_ipd(runif(n), alpha)
  T <- r_kg(n, alpha, beta_par)

  Y <- pmin(X, T)
  D <- as.integer(X <= T)

  prcs_from_YD(Y, D, m, W)
}


# Impose a PRCS design on a complete failure-time data set. The lifetimes are
# the observed data; only the censoring times are generated, from the
# Koziol-Green model with the shape parameter estimated from the complete
# sample and with the censoring parameter beta_impose.
impose_prcs <- function(x_complete, m, W, beta_impose, alpha_impose = NULL,
                        gen = NULL)
{
  n <- length(x_complete)

  check_scheme(n, m, W)

  if (is.null(alpha_impose))
    alpha_impose <- mle_complete_ipd(x_complete)

  # The censoring times are drawn first, n of them, and the withdrawals
  # afterwards. The order in which the generator is consumed is part of the
  # definition of the design.
  v <- if (is.null(gen)) runif(n) else lcg_unif(gen, n)

  T <- r_kg(n, alpha_impose, beta_impose, v = v)

  Y <- pmin(x_complete, T)
  D <- as.integer(x_complete <= T)

  out <- prcs_from_YD(Y, D, m, W, gen = gen)

  out$alpha_impose <- alpha_impose
  out$beta_impose  <- beta_impose

  out
}


# Closed-form MLE of alpha from a complete sample, obtained from
# d/d alpha [ n log alpha + alpha sum log z_i - sum log(1 + x_i) ] = 0.
mle_complete_ipd <- function(x)
{
  z <- x / (1 + x)

  -length(x) / sum(log(z))
}


# -----------------------------------------------------------------------------
# Likelihood, Eqs. (20) and (21)
# -----------------------------------------------------------------------------

zvec <- function(dat)
{
  dat$y / (1 + dat$y)
}


loglik_prcs <- function(alpha, beta_par, dat)
{
  if (alpha <= 0 || beta_par <= 0)
    return(-Inf)

  y <- dat$y
  W <- dat$W
  m <- dat$m
  k <- dat$k

  z  <- y / (1 + y)
  lz <- log(z)

  la <- alpha * lz          # log(z^alpha)

  # log(1 - z^alpha), accurate when z^alpha is close to one
  l1 <- log(-expm1(la))

  if (any(!is.finite(l1)))
    return(-Inf)

  m * log(alpha) + (m - k) * log(beta_par) +
    (alpha - 1) * sum(log(y)) - (alpha + 1) * sum(log(1 + y)) +
    sum((beta_par + (1 + beta_par) * W) * l1)
}


# The likelihood written as the product over the observed events, before the
# Koziol-Green substitutions are made. Used only by check_formulas.R.
loglik_prcs_direct <- function(alpha, beta_par, dat)
{
  y <- dat$y
  d <- dat$d
  W <- dat$W

  fX <- d_ipd(y, alpha)
  SX <- s_ipd(y, alpha)
  fT <- d_kg(y, alpha, beta_par)
  ST <- 1 - p_kg(y, alpha, beta_par)

  sum(d * log(fX * ST) + (1 - d) * log(fT * SX) + W * log(SX * ST))
}


score_prcs <- function(alpha, beta_par, dat)
{
  y <- dat$y
  W <- dat$W
  m <- dat$m
  k <- dat$k

  z  <- y / (1 + y)
  lz <- log(z)
  za <- exp(alpha * lz)
  om <- -expm1(alpha * lz)        # 1 - z^alpha

  s_alpha <- m / alpha + sum(lz) -
    sum((beta_par + (1 + beta_par) * W) * za * lz / om)

  s_beta <- (m - k) / beta_par + sum((1 + W) * log(om))

  c(alpha = s_alpha, beta = s_beta)
}


# Observed information matrix, the three elements given in the chapter.
obs_info_prcs <- function(alpha, beta_par, dat)
{
  y <- dat$y
  W <- dat$W
  m <- dat$m
  k <- dat$k

  z  <- y / (1 + y)
  lz <- log(z)
  za <- exp(alpha * lz)
  om <- -expm1(alpha * lz)

  J_aa <- m / alpha^2 +
    sum((beta_par + (1 + beta_par) * W) * za * lz^2 / om^2)

  J_bb <- (m - k) / beta_par^2

  J_ab <- sum((1 + W) * za * lz / om)

  matrix(c(J_aa, J_ab, J_ab, J_bb), 2, 2,
         dimnames = list(c("alpha", "beta"), c("alpha", "beta")))
}


# Profile estimator of beta, Eq. (24).
beta_profile <- function(alpha, dat)
{
  if (dat$m == dat$k)
    return(NA_real_)

  z  <- dat$y / (1 + dat$y)
  om <- -expm1(alpha * log(z))

  den <- sum((1 + dat$W) * log(om))

  # For a very large alpha every z_i^alpha underflows, the denominator
  # becomes zero and the profile estimator is not defined.
  if (!is.finite(den) || den >= 0)
    return(NA_real_)

  -(dat$m - dat$k) / den
}


# The one-dimensional equation in alpha obtained after beta is profiled out.
profile_score <- function(alpha, dat)
{
  b <- beta_profile(alpha, dat)

  if (!is.finite(b) || b <= 0)
    return(NA_real_)

  score_prcs(alpha, b, dat)["alpha"]
}


# -----------------------------------------------------------------------------
# Maximum likelihood estimation
#
# Primary route: solve the profiled one-dimensional equation.
# Independent route: a two-dimensional Newton-Raphson iteration that stops
# only when the change in the log-likelihood and the largest relative change
# in the parameters are both below their tolerances.
# -----------------------------------------------------------------------------

fit_mle <- function(dat, lower = 1e-3, upper = 200, conf = 0.95)
{
  if (dat$m == dat$k)
    return(list(ok = FALSE, msg = "m = k, beta lies on the boundary"))

  gr <- exp(seq(log(lower), log(upper), length.out = 400L))
  sc <- vapply(gr, function(a) profile_score(a, dat), numeric(1))

  ok <- is.finite(sc)

  if (sum(ok) < 2L)
    return(list(ok = FALSE, msg = "profile score not computable"))

  gr <- gr[ok]
  sc <- sc[ok]

  sgn <- which(sc[-length(sc)] * sc[-1] < 0)

  if (!length(sgn))
    return(list(ok = FALSE, msg = "no sign change in the profile score"))

  i  <- sgn[1L]
  rt <- uniroot(function(a) profile_score(a, dat),
                lower = gr[i], upper = gr[i + 1L], tol = 1e-12)

  alpha_hat <- rt$root
  beta_hat  <- beta_profile(alpha_hat, dat)

  if (!is.finite(beta_hat) || beta_hat <= 0)
    return(list(ok = FALSE, msg = "profile estimator of beta not positive"))

  J <- obs_info_prcs(alpha_hat, beta_hat, dat)

  ev  <- eigen(J, symmetric = TRUE, only.values = TRUE)$values
  pdJ <- all(ev > 0)

  V <- tryCatch(solve(J), error = function(e) NULL)

  if (is.null(V) || any(diag(V) <= 0))
    return(list(ok = FALSE, msg = "information matrix not invertible"))

  se <- sqrt(diag(V))
  zq <- qnorm(1 - (1 - conf) / 2)

  est <- c(alpha = alpha_hat, beta = beta_hat)

  ci_plain <- cbind(lower = pmax(0, est - zq * se),
                    upper = est + zq * se)

  ci_log <- cbind(lower = est * exp(-zq * se / est),
                  upper = est * exp( zq * se / est))

  list(ok        = TRUE,
       alpha     = alpha_hat,
       beta      = beta_hat,
       est       = est,
       se        = se,
       J         = J,
       V         = V,
       pd        = pdJ,
       score     = score_prcs(alpha_hat, beta_hat, dat),
       loglik    = loglik_prcs(alpha_hat, beta_hat, dat),
       ci_plain  = ci_plain,
       ci_log    = ci_log)
}


# Independent iterative fit, used as a check of the profiled solution and of
# the analytic information matrix.
fit_mle_newton <- function(dat, start = c(1, 1),
                           eps_loglik = EPS_LOGLIK,
                           eps_par    = EPS_PAR,
                           max_iter   = MAX_ITER)
{
  th <- start
  ll <- loglik_prcs(th[1L], th[2L], dat)

  for (it in seq_len(max_iter))
  {
    s <- score_prcs(th[1L], th[2L], dat)
    J <- obs_info_prcs(th[1L], th[2L], dat)

    step <- tryCatch(solve(J, s), error = function(e) NULL)

    if (is.null(step))
      return(list(ok = FALSE, msg = "singular information matrix", iter = it))

    lam <- 1
    repeat
    {
      new <- th + lam * step

      if (all(new > 0))
      {
        ll_new <- loglik_prcs(new[1L], new[2L], dat)

        if (is.finite(ll_new) && ll_new >= ll - 1e-8)
          break
      }

      lam <- lam / 2

      if (lam < 1e-12)
        return(list(ok = FALSE, msg = "step halving failed", iter = it))
    }

    d_ll  <- abs(ll_new - ll)
    d_par <- max(abs(new - th) / pmax(abs(th), 1e-12))

    th <- new
    ll <- ll_new

    # Both conditions must hold. The condition on the parameters is not
    # optional: under heavy censoring the log-likelihood is flat along a
    # ridge and a rule based on the log-likelihood alone stops while the
    # estimates are still moving.
    if (d_ll < eps_loglik && d_par < eps_par)
      return(list(ok = TRUE, alpha = th[1L], beta = th[2L],
                  est = c(alpha = th[1L], beta = th[2L]),
                  loglik = ll, iter = it,
                  score = score_prcs(th[1L], th[2L], dat)))
  }

  list(ok = FALSE, msg = "no convergence within max_iter", iter = max_iter)
}


# -----------------------------------------------------------------------------
# Reliability and hazard rate at a fixed time
# -----------------------------------------------------------------------------

surv_hazard <- function(t, alpha)
{
  c(S = s_ipd(t, alpha), h = h_ipd(t, alpha))
}


# -----------------------------------------------------------------------------
# Bayesian estimation, Metropolis-Hastings within Gibbs
# -----------------------------------------------------------------------------

# rate of the full conditional Gamma(m - k + a2, rate) of beta
beta_rate <- function(alpha, dat, b2)
{
  z  <- dat$y / (1 + dat$y)
  om <- -expm1(alpha * log(z))

  b2 - sum((1 + dat$W) * log(om))
}


# log of the conditional posterior of alpha, up to an additive constant
log_cond_alpha <- function(alpha, beta_par, dat, a1, b1)
{
  if (alpha <= 0)
    return(-Inf)

  y <- dat$y
  W <- dat$W

  z  <- y / (1 + y)
  om <- -expm1(alpha * log(z))

  if (any(om <= 0))
    return(-Inf)

  (dat$m + a1 - 1) * log(alpha) - b1 * alpha +
    (alpha - 1) * sum(log(y)) - (alpha + 1) * sum(log(1 + y)) +
    sum((beta_par + (1 + beta_par) * W) * log(om))
}


mcmc_prcs <- function(dat, a1, b1, a2, b2,
                      zeta = 10000L, zeta0 = 2000L,
                      sigma_alpha = NULL, start = NULL, tune = TRUE,
                      gen = NULL)
{
  # gen supplied: every draw comes from the portable generator, so that the
  # chain is fixed by the seed alone. gen NULL: R's own generator, which is
  # what the Monte Carlo study uses.
  draw_u <- function(n) if (is.null(gen)) runif(n) else lcg_unif(gen, n)
  draw_z <- function(n) if (is.null(gen)) rnorm(n) else lcg_norm(gen, n)

  draw_gamma <- function(shape, rate)
    if (is.null(gen)) rgamma(1, shape = shape, rate = rate) else
      lcg_gamma_int(gen, shape, rate)

  if (!is.null(gen) && abs(a2 - round(a2)) > 1e-12)
    stop("With the portable generator a2 must be an integer, so that the ",
         "shape m - k + a2 of the full conditional of beta is an integer.")

  if (is.null(start))
  {
    ft <- fit_mle(dat)

    start <- if (isTRUE(ft$ok)) c(ft$alpha, ft$beta) else c(1, 1)
  }

  if (is.null(sigma_alpha))
    sigma_alpha <- 0.15

  shape_beta <- dat$m - dat$k + a2

  # Short pilot run, so that the acceptance rate of the one-dimensional
  # random-walk update falls between about 0.25 and 0.50.
  if (tune)
  {
    for (round_i in 1:4)
    {
      acc <- 0
      a   <- start[1L]
      b   <- start[2L]

      for (j in 1:250)
      {
        as_ <- exp(log(a) + sigma_alpha * draw_z(1))

        lr <- log_cond_alpha(as_, b, dat, a1, b1) -
          log_cond_alpha(a, b, dat, a1, b1) + log(as_) - log(a)

        if (is.finite(lr) && log(draw_u(1)) < lr)
        {
          a   <- as_
          acc <- acc + 1
        }

        b <- draw_gamma(shape_beta, beta_rate(a, dat, b2))
      }

      rate <- acc / 250

      if (rate > 0.25 && rate < 0.50)
        break

      sigma_alpha <- sigma_alpha * exp((rate - 0.35) * 2)
    }
  }

  alpha_ch <- numeric(zeta)
  beta_ch  <- numeric(zeta)

  a <- start[1L]
  b <- start[2L]

  lc  <- log_cond_alpha(a, b, dat, a1, b1)
  acc <- 0L

  for (j in seq_len(zeta))
  {
    as_ <- exp(log(a) + sigma_alpha * draw_z(1))

    lc_s <- log_cond_alpha(as_, b, dat, a1, b1)

    # the ratio carries the Jacobian factor alpha* / alpha of the
    # log-scale proposal
    lr <- lc_s - lc + log(as_) - log(a)

    if (is.finite(lr) && log(draw_u(1)) < lr)
    {
      a   <- as_
      lc  <- lc_s
      acc <- acc + 1L
    }

    b <- draw_gamma(shape_beta, beta_rate(a, dat, b2))

    lc <- log_cond_alpha(a, b, dat, a1, b1)

    alpha_ch[j] <- a
    beta_ch[j]  <- b
  }

  keep <- (zeta0 + 1L):zeta

  list(alpha       = alpha_ch[keep],
       beta        = beta_ch[keep],
       acc_rate    = acc / zeta,
       sigma_alpha = sigma_alpha,
       zeta        = zeta,
       zeta0       = zeta0)
}


# Highest posterior density interval, Chen and Shao (1999).
hpd <- function(x, level = 0.95)
{
  x <- sort(x)
  N <- length(x)
  r <- floor(level * N)

  if (r < 1L || r >= N)
    return(c(lower = NA_real_, upper = NA_real_))

  j <- which.min(x[(1L + r):N] - x[1L:(N - r)])

  c(lower = x[j], upper = x[j + r])
}


# -----------------------------------------------------------------------------
# Marginal posterior of alpha by quadrature
#
# The full conditional of beta is a gamma distribution, so beta can be
# integrated out exactly. This gives a one-dimensional posterior for alpha
# and hence exact posterior means, which check_formulas.R compares with the
# output of the sampler.
# -----------------------------------------------------------------------------

log_marg_alpha <- function(alpha, dat, a1, b1, a2, b2)
{
  if (alpha <= 0)
    return(-Inf)

  y <- dat$y
  W <- dat$W

  z  <- y / (1 + y)
  om <- -expm1(alpha * log(z))

  if (any(om <= 0))
    return(-Inf)

  rate <- b2 - sum((1 + W) * log(om))

  if (rate <= 0)
    return(-Inf)

  (dat$m + a1 - 1) * log(alpha) - b1 * alpha +
    (alpha - 1) * sum(log(y)) - (alpha + 1) * sum(log(1 + y)) +
    sum(W * log(om)) -
    (dat$m - dat$k + a2) * log(rate)
}


post_moments_quad <- function(dat, a1, b1, a2, b2,
                              lower = 1e-6, upper = 200, ngrid = 20001L)
{
  gr <- seq(lower, upper, length.out = ngrid)
  lm <- vapply(gr, log_marg_alpha, numeric(1), dat, a1, b1, a2, b2)

  mx <- max(lm[is.finite(lm)])
  w  <- exp(lm - mx)
  w[!is.finite(w)] <- 0

  trap <- function(f) sum((f[-1] + f[-ngrid]) / 2 * diff(gr))

  cst <- trap(w)

  e_alpha  <- trap(w * gr) / cst
  e_alpha2 <- trap(w * gr^2) / cst

  rate_gr <- vapply(gr, beta_rate, numeric(1), dat, b2)
  shape   <- dat$m - dat$k + a2

  e_beta  <- trap(w * (shape / rate_gr)) / cst
  e_beta2 <- trap(w * (shape * (shape + 1) / rate_gr^2)) / cst

  c(mean_alpha = e_alpha,
    sd_alpha   = sqrt(e_alpha2 - e_alpha^2),
    mean_beta  = e_beta,
    sd_beta    = sqrt(e_beta2 - e_beta^2))
}


# -----------------------------------------------------------------------------
# Figures
#
# Rendered to a temporary local file and copied into figures/, because
# postscript() and pdf() fail on some systems when they write directly to a
# synchronised drive or to a path with non-ASCII characters.
# -----------------------------------------------------------------------------

save_figure <- function(name, draw, width = 14.5, height = 4.8)
{
  dir.create("figures", showWarnings = FALSE)

  tmp_eps <- tempfile(fileext = ".eps")
  tmp_pdf <- tempfile(fileext = ".pdf")

  postscript(tmp_eps, width = width, height = height, horizontal = FALSE,
             onefile = FALSE, paper = "special", family = "serif")
  tryCatch(draw(), finally = if (dev.cur() > 1) dev.off())

  pdf(tmp_pdf, width = width, height = height, family = "serif")
  tryCatch(draw(), finally = if (dev.cur() > 1) dev.off())

  file.copy(tmp_eps, file.path("figures", paste0(name, ".eps")),
            overwrite = TRUE)
  file.copy(tmp_pdf, file.path("figures", paste0(name, ".pdf")),
            overwrite = TRUE)

  unlink(c(tmp_eps, tmp_pdf))

  invisible(file.path("figures", paste0(name, c(".eps", ".pdf"))))
}


# -----------------------------------------------------------------------------
# Reporting
# -----------------------------------------------------------------------------

fmt4 <- function(x)
{
  ifelse(is.na(x), "NA", formatC(x, format = "f", digits = 4))
}


new_report <- function(path)
{
  con <- file(path, open = "wt", encoding = "UTF-8")

  function(...)
  {
    txt <- paste0(...)
    cat(txt, "\n", sep = "", file = con, append = TRUE)
    cat(txt, "\n", sep = "")
    invisible(NULL)
  }
}


rule <- function(ch = "=", n = 78)
{
  paste(rep(ch, n), collapse = "")
}


# One emitter for every table in the chapter, so that all of them carry the
# same font size and the same column spacing and are flush to both margins.
tex_table <- function(body, header, align, caption, label,
                      note = NULL, midrules = NULL)
{
  out <- c("\\begin{table}[!htbp]",
           "\\centering",
           paste0("\\caption{", caption, "}"),
           paste0("\\label{", label, "}"),
           "\\footnotesize",
           "\\setlength{\\tabcolsep}{2pt}",
           paste0("\\begin{tabular*}{\\textwidth}{@{\\extracolsep{\\fill}}",
                  align, "}"),
           "\\toprule",
           header,
           "\\midrule")

  for (i in seq_along(body))
  {
    out <- c(out, body[i])

    if (!is.null(midrules) && i %in% midrules)
      out <- c(out, "\\midrule")
  }

  out <- c(out, "\\bottomrule", "\\end{tabular*}")

  if (!is.null(note))
    out <- c(out, paste0("\\\\[2pt]\\footnotesize ", note))

  c(out, "\\end{table}")
}


# -----------------------------------------------------------------------------
# The censoring schemes of the chapter
# -----------------------------------------------------------------------------

rep_scheme <- function(...)
{
  as.integer(c(...))
}


# Table 1 of the chapter
SCHEMES_SIM <- list(
  S1  = list(n = 25, m = 18, pattern = "Early",
             W = rep_scheme(rep(1, 7), rep(0, 11))),
  S2  = list(n = 30, m = 22, pattern = "Early",
             W = rep_scheme(rep(1, 8), rep(0, 14))),
  S3  = list(n = 35, m = 25, pattern = "Late",
             W = rep_scheme(rep(0, 15), rep(1, 10))),
  S4  = list(n = 40, m = 28, pattern = "Late",
             W = rep_scheme(rep(0, 16), rep(1, 12))),
  S5  = list(n = 45, m = 30, pattern = "Uniform",
             W = rep_scheme(rep(c(1, 0), 15))),
  S6  = list(n = 50, m = 35, pattern = "Middle",
             W = rep_scheme(rep(0, 10), rep(1, 15), rep(0, 10))),
  S7  = list(n = 55, m = 38, pattern = "Increasing",
             W = rep_scheme(rep(0, 26), rep(1, 7), rep(2, 5))),
  S8  = list(n = 60, m = 40, pattern = "Decreasing",
             W = rep_scheme(rep(2, 5), rep(1, 10), rep(0, 25))),
  S9  = list(n = 70, m = 50, pattern = "Alternating",
             W = rep_scheme(rep(c(1, 0), 20), rep(0, 10))),
  S10 = list(n = 80, m = 55, pattern = "Mixed",
             W = rep_scheme(rep(2, 5), rep(0, 10), rep(1, 15), rep(0, 25)))
)


# Table 6 of the chapter, the designs imposed on the insulating fluid data
SCHEMES_E <- list(
  E1 = list(n = 19, m = 14, pattern = "Early",
            W = rep_scheme(rep(1, 5), rep(0, 9))),
  E2 = list(n = 19, m = 14, pattern = "Late",
            W = rep_scheme(rep(0, 9), rep(1, 5))),
  E3 = list(n = 19, m = 14, pattern = "Uniform",
            W = rep_scheme(1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 0, 0, 0, 0)),
  E4 = list(n = 19, m = 14, pattern = "Middle",
            W = rep_scheme(rep(0, 4), rep(1, 5), rep(0, 5)))
)


# -----------------------------------------------------------------------------
# The engineering data set of Section 8
#
# Times to breakdown of an insulating fluid between electrodes at a voltage of
# 34 kV, in minutes. Nineteen specimens, all of them completely observed.
#
# Source: Nelson, W. (1982). Applied Life Data Analysis. Wiley, New York.
#         Table 1.1 of Chapter 3, page 105, the 34 kV column. The table gives
#         times to breakdown at seven voltages, from 26 kV to 38 kV.
#
# Check: Nelson gives the sample mean as 14.3 minutes on page 210. The values
#        below give 272.82 / 19 = 14.3589, which rounds to that.
# -----------------------------------------------------------------------------

BREAKDOWN34 <- c(0.19, 0.78, 0.96, 1.31, 2.78, 3.16, 4.15, 4.67, 4.85, 6.50,
                 7.35, 8.01, 8.27, 12.06, 31.75, 32.52, 33.91, 36.71, 72.89)


# -----------------------------------------------------------------------------
# A fixed PRCS data set, used by check_formulas.R, so that every check is
# carried out on the same likelihood and the same numbers.
# -----------------------------------------------------------------------------

CHECK_DATA <- list(
  y = c(0.4123, 0.6871, 0.9042, 1.1537, 1.4419,
        1.8265, 2.3391, 3.0142, 4.2077, 6.5310),
  d = c(1L, 1L, 0L, 1L, 0L, 1L, 1L, 0L, 1L, 0L),
  W = c(1L, 0L, 1L, 0L, 1L, 0L, 1L, 0L, 0L, 0L),
  n = 14L,
  m = 10L
)
CHECK_DATA$k <- sum(CHECK_DATA$d)
