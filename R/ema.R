# EMA skater: conjugation of a leaf with the EMA transform.
# Port of skaters/ema.py.

ema <- function(alpha = 0.05, k = 1L) {
  conjugate(leaf(k), ema_transform(alpha), k)
}
