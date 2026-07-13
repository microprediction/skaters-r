# Sticky / lattice projection: mean-preserving near-Dirac atoms at the
# exact values a series revisits. Port of skaters/sticky.py.

sticky <- function(base, k = 1L, propensity_alpha = 0.05, spike_frac = 0.005,
                   thresh_mult = 1.8, max_atoms = 6L, prune_eps = 1e-6) {
  function(y, state = NULL) {
    if (is.null(state)) state <- list(base = NULL, vals = numeric(0), wts = numeric(0))

    r <- base(y, state$base)
    state$base <- r$state
    dists <- r$dists

    # Recency-weighted frequency table of exact values (insertion-ordered,
    # like the Python dict): decay, drop below eps, then credit y.
    wts <- state$wts * (1.0 - propensity_alpha)
    keep <- wts >= prune_eps
    vals <- state$vals[keep]
    wts <- wts[keep]
    hit <- which(vals == y)
    if (length(hit) > 0) {
      wts[hit[1]] <- wts[hit[1]] + propensity_alpha
    } else {
      vals <- c(vals, y)
      wts <- c(wts, propensity_alpha)
    }
    state$vals <- vals
    state$wts <- wts

    # Lattice atoms: revisited values above the noise floor, top by weight
    # (stable sort, insertion order breaks ties, as in Python).
    thr <- thresh_mult * propensity_alpha
    sel <- wts > thr
    av <- vals[sel]
    aw <- wts[sel]
    if (length(aw) > 0) {
      o <- order(-aw)
      av <- av[o]; aw <- aw[o]
      if (length(aw) > max_atoms) {
        av <- av[seq_len(max_atoms)]
        aw <- aw[seq_len(max_atoms)]
      }
    }

    out <- vector("list", length(dists))
    for (i in seq_along(dists)) {
      d <- dists[[i]]
      if (length(aw) == 0) {
        out[[i]] <- d
        next
      }
      sw <- sum(aw)
      P <- min(sw, 0.999)
      pc <- 1.0 - P
      atom_mean <- sum(aw * av) / sw
      spike_std <- max(spike_frac * dist_std(d), 1e-9)
      if (pc <= 1e-9) {
        out[[i]] <- dist_new(aw / sw, av, rep(spike_std, length(aw)))
        next
      }
      mu <- dist_mean(d)
      delta <- P * (mu - atom_mean) / pc
      out[[i]] <- dist_new(c(P * (aw / sw), pc * d$w),
                           c(av, d$m + delta),
                           c(rep(spike_std, length(aw)), d$s))
    }
    list(dists = out, state = state)
  }
}
