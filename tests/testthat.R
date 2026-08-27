# testthat harness. The substantive assertions live in tests/testthat/ and are
# PORTED from the original hand-rolled scripts assertion-for-assertion; the
# parity vectors and the 1e-6 tolerance are the five-port contract and must
# never be weakened in a test-framework migration. (A migration that rewrites
# assertions instead of porting them is how parity goes falsely green.)
library(testthat)
library(skaters)
test_check("skaters")
