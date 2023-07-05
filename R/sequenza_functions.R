### functions from sequenza version 2.1.2
### copied over by Ariana Huebner

types.matrix <- function (CNt.min, CNt.max, CNn = 2) {
    cn.ratio.vect <- seq(from = CNt.min / CNn, to = CNt.max / CNn, by = 1 / CNn)
    CNt <- cn.ratio.vect * CNn
    mut.comb <- lapply(CNt, FUN = function(x) seq(from = 0, to = x))
    times.muts <- sapply(mut.comb, length)
    data.frame(CNn = CNn, CNt = rep(CNt, times = times.muts), Mt = unlist(mut.comb))
}

theoretical.mufreq <- function (Mt, CNt, CNn = 2, cellularity) {
    normal.alleles <- (CNt - Mt) * cellularity + CNn * (1 - cellularity)
    all.alleles <- (CNt * cellularity) + CNn * (1 - cellularity)
    1 - (normal.alleles / all.alleles)
}

mufreq.dpois <- function (mufreq, mufreq.model, depth.t, seq.errors = 0.01, ...) {
    mufreq.model[mufreq.model == 0] <- seq.errors
    n.success <- round(mufreq * depth.t, 0)
    dpois(x = n.success, lambda = mufreq.model * depth.t, ...)
}