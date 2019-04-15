library(ggplot2)
library(dplyr)

set.seed(12345);

thetas <- c(0.01, 0.02, 0.03, 0.05, 0.1);
m <- length(thetas);
ns <- rep(100000, m);

n.events <- rbinom(m, ns, thetas);
n.nonevents <- ns - n.events;

groups <- unlist(mapply(function(g, n) rep(g, n), 1:m, ns, SIMPLIFY=FALSE));
events <- unlist(mapply(function(n1, n0) c(rep(1, n1), rep(0, n0)), n.events, n.nonevents, SIMPLIFY=FALSE));

d <- data.frame(
	group = groups,
	event = events
);

# cohort analysis

ct <- with(d, table(group, event));
hs.cohort <- lapply(2:m,
	function(j) {
		fisher.test(ct[c(1, j), ])
	}
);

# case-control analysis with uniform sampling

n.case <- 1000;
n.control <- 2000;

d.cc <- rbind(
	sample_n(filter(d, event == 1), n.case),	
	sample_n(filter(d, event == 0), n.control)
);

ct.cc <- with(d.cc, table(group, event));
hs.cc <- lapply(2:m,
	function(j) {
		fisher.test(ct.cc[c(1, j), ])
	}
);

# case-control analysis with stratified sampling

prevalence.group1 <- ns[1] / sum(ns);

d.ccs <- rbind(
	sample_n(filter(d, event == 1, group == 1), round(n.case * prevalence.group1)),
	sample_n(filter(d, event == 1, group != 1), round(n.case * (1 - prevalence.group1))),
	sample_n(filter(d, event == 0, group == 1), round(n.control * prevalence.group1)),
	sample_n(filter(d, event == 0, group != 1), round(n.control * (1 - prevalence.group1)))
);

ct.ccs <- with(d.ccs, table(group, event));
hs.ccs <- lapply(2:m,
	function(j) {
		fisher.test(ct.ccs[c(1, j), ])
	}
);

or.h.cohort <- unlist(lapply(hs.cohort, function(h) h$estimate));
or.h.cc <- unlist(lapply(hs.cc, function(h) h$estimate));
or.h.ccs <- unlist(lapply(hs.ccs, function(h) h$estimate));

s <- data.frame(
	group = 1:m,
	theta = thetas,
	or_hat_cohort = c(1, or.h.cohort),
	or_hat_cc = c(1, or.h.cc),
	or_hat_ccs = c(1, or.h.ccs)
);

