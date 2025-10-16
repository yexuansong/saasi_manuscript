#!/bin/bash

# Simple ASR methods comparison
N_SIMS=${1:-10}
NTAXA=${2:-100}
MAX_TIME=${3:-100}
OUTPUT_DIR=${4:-"simple_asr"}

mkdir -p "$OUTPUT_DIR"

# Create R script
cat > "$OUTPUT_DIR/run_asr.R" << 'EOF'
library(saasi)
library(diversitree)
library(ape)
library(phytools)
library(tidytree)



#' Helper fn for `probability_density`.
#'
#' @noRd
c1 <- function(lambda, mu, psi) {
  return(abs(sqrt((lambda - mu - psi)^2 + (4 * lambda * psi))))
}

#' Helper fn for `probability_density`.
#'
#' @noRd
c2 <- function(lambda, mu, psi) {
  return(-(lambda - mu - psi) / (c1(lambda, mu, psi)))
}

#' Helper fn for `probability_density`.
#'
#' @noRd
q <- function(lambda, mu, psi, t) {
  d1 <- 2 * (1 - (c2(lambda, mu, psi))^2)
  d2 <- exp(-c1(lambda, mu, psi) * t) * (1 - c2(lambda, mu, psi))^2
  d3 <- exp(c1(lambda, mu, psi) * t) * (1 + c2(lambda, mu, psi))^2
  return(d1 + d2 + d3)
}

#' Helper fn for `mle_lm`.
#'
#' @noRd
probability_density <- function(lambda, mu, psi,
                                internal_node_times, leaf_node_times) {
  # Add checks for invalid inputs
  sorted_x <- sort(internal_node_times)
  sorted_y <- sort(leaf_node_times)
  d1 <- length(internal_node_times) * log(lambda)
  d2 <- sum(-log(q(lambda, mu, psi, sorted_x)))
  d3 <- sum(log(psi) + log(q(lambda, mu, psi, sorted_y)))
  ret <- d1 + d2 + d3
  return(ret)
}

#' Estimate speciation/extinction rates for a tree
#'
#' This is done with a maximum likelihood method,  implemented mainly by
#' \href{https://github.com/yexuansong}{@yexuansong}, from methods described in
#' \href{https://doi.org/10.1093/molbev/msr217}{Stadler et al. (2012)}.
#'
#' @param phy A `phylo` phylogenetic tree (`ape` format).
#' @param lambda An initial "guess" for speciation, used in subsequent formulae.
#' @param mu An initial "guess" for extinction, used in subsequent formulae.
#' @param psi Sampling rate.
#' @param method See `method` parameter in \link{mle}.
#' @param lower A two-element vector containing the lower bound of the speciation 
#' and extinction values (geater than 0), the default is c(0.001,0.001).
#' @param upper A two-element vector containing the upper bound of the speciation 
#' and extinction values, should be similar to the sampling rate, still need to find
#' a better justification for the upcoming version.
#' @return A two-element vector containing speciation and extinction values
#' estimated from `phy`.
#'@export
mle_lm <- function(phy, lambda, mu, psi, method = "L-BFGS-B", lower = c(0.001,0.001), upper) {
  node_depths <- ape::node.depth.edgelength(phy)
  node_times <- max(node_depths) - node_depths
  # Total number of nodes == number of non-leaf nodes * 2 + 1
  nnode <- phy[["Nnode"]] * 2 + 1
  nleaf_node <- nnode - phy[["Nnode"]]
  leaf_node_times <- node_times[1:nleaf_node]
  internal_node_times <- node_times[(nleaf_node + 1):nnode]
  
  negative_log_likelihood <- function(lambda, mu) {
    positive_ret <- probability_density(lambda, mu, psi,
                                        internal_node_times, leaf_node_times)
    return(-positive_ret)
  }
  
  fit <- stats4::mle(negative_log_likelihood,
                     start = list(lambda = lambda, mu = mu),
                     method = method,
                     lower = lower,
                     upper = upper)
  
  ret <- unname(c(stats4::coef(fit)[1], stats4::coef(fit)[2]))
  names(ret) = c("speciation","extinction")
  return(ret)
}


args <- commandArgs(trailingOnly = TRUE)
sim_num <- as.numeric(args[1])
ntaxa <- as.numeric(args[2])
max_time <- as.numeric(args[3])
output_dir <- args[4]

set.seed(sim_num)

# Generate tree
pars <- data.frame(state=c(1,2), prior=c(0.5,0.5), lambda=c(3,1.5), mu=c(0.05,0.1), psi=c(0.1,1))
q_matrix <- matrix(0.5, 2, 2)
q_matrix[2,1] <- 0.45
diag(q_matrix) <- NA

phy <- sim_bds_tree(pars, q_matrix, x0=1, max_taxa=ntaxa, max_t=max_time, include_extinct=FALSE)
h <- history.from.sim.discrete(phy, 1:2)

# Process tree
node_depths <- node.depth.edgelength(phy)
tmrca <- max(node_depths)
tips_to_drop <- phy$tip.label[abs(node_depths[1:length(phy$tip.label)] - tmrca) <= 0.01]
if(length(tips_to_drop) > 0) {
    phy <- drop.tip(phy, tips_to_drop)
    phy$tip.state <- phy$tip.state[setdiff(names(phy$tip.state), tips_to_drop)]
}

phy_info_before_ds <- as_tibble(phy)%>% mutate(State = c(factor(h$tip.state),factor(h$node.state))[label])
write.csv(phy_info_before_ds, file=paste0(output_dir,"/true_tree_before_ds_",sim_num,".csv"))

# Export tree and metadata
write.tree(phy, file=paste0(output_dir, "/tree_before_ds_", sim_num, ".nwk"))
metadata_df_before_ds <- data.frame(
    node = phy$tip.label,
    state = paste0("State", phy$tip.state),
    numeric_state = phy$tip.state
)
write.csv(metadata_df_before_ds, file=paste0(output_dir, "/metadata_before_ds_", sim_num, ".csv"), row.names=FALSE)
write.table(metadata_df_before_ds, file=paste0(output_dir, "/metadata_before_ds_", sim_num, ".txt"), sep=",", row.names=FALSE)

frac <- pars$psi[2]/pars$psi[1]
num_drop <- table(phy$tip.state)[2] - as.numeric(floor(table(phy$tip.state)[2]/frac)) # number of tips to drop
tip_to_drop <- sample( names(phy$tip.state[phy$tip.state == 2]),num_drop)
phy2 <- drop.tip(phy, tip_to_drop)

true_phy_info_new2 <- as_tibble(phy2) %>% mutate(State = c(factor(h$tip.state),factor(h$node.state))[label])
phy2$tip.state <- phy2$tip.state[setdiff(names(phy2$tip.state), tip_to_drop)]
phy <- phy2

phy_info <- as_tibble(phy)%>% mutate(State = c(factor(h$tip.state),factor(h$node.state))[label])
write.csv(phy_info, file=paste0(output_dir,"/true_tree_",sim_num,".csv"))

# Export tree and metadata
write.tree(phy, file=paste0(output_dir, "/tree_", sim_num, ".nwk"))
metadata_df <- data.frame(
    node = phy$tip.label,
    state = paste0("State", phy$tip.state),
    numeric_state = phy$tip.state
)
write.csv(metadata_df, file=paste0(output_dir, "/metadata_", sim_num, ".csv"), row.names=FALSE)
write.table(metadata_df, file=paste0(output_dir, "/metadata_", sim_num, ".txt"), sep=",", row.names=FALSE)

timing_results <- data.frame(
    simulation = integer(0),
    method = character(0),
    time_seconds = numeric(0),
    success = logical(0)
)

# ACE
ace_start <- Sys.time()
ace_success <- FALSE
ace_result <- ace(phy$tip.state,phy,type="discrete",model="ARD")
write.csv(ace_result$lik.anc, file=paste0(output_dir,"/ace_",sim_num,".csv"))
ace_success <- TRUE
ace_time <- as.numeric(Sys.time() - ace_start)
timing_results <- rbind(timing_results, data.frame(simulation=sim_num, method="ace", time_seconds=ace_time, success=ace_success))

# SAASI

pars <- data.frame(state=c(1,2), prior=c(0.5,0.5), lambda=c(3,1.5), mu=c(0.05,0.1), psi=c(0.1,0.1))

saasi_start <- Sys.time()
saasi_success <- FALSE
tryCatch({
    saasi_result <- saasi(phy, pars, q_matrix)
    write.csv(saasi_result, file=paste0(output_dir, "/saasi_", sim_num, ".csv"))
    saasi_success <- TRUE
}, error = function(e) {})
saasi_time <- as.numeric(Sys.time() - saasi_start)
timing_results <- rbind(timing_results, data.frame(simulation=sim_num, method="saasi", time_seconds=saasi_time, success=saasi_success))

# SAASI with pars estimations
saasi_par_start <- Sys.time()
saasi_par_success <- FALSE
tryCatch({
	estimates <- mle_lm(phy,lambda = 2, mu = 0.1, psi = 1, lower = c(0.001,0.001), upper = c(5,5))
	q_matrix_est <- extract_ace_q(ace_result)
	pars_est <- data.frame(state=c(1,2),
                       prior=c(0.5,0.5),
                       lambda=rep(estimates[1],2),
                       mu=rep(estimates[2],2),
                       psi=c(.1,.1))
	saasi_est_result <- saasi(phy,pars_est, q_matrix_est)
	write.csv(saasi_est_result, file=paste0(output_dir,"/saasi_pars_est_",sim_num,".csv"))
	saasi_par_success <- TRUE
}, error = function(e){})
saasi_par_est_time <- as.numeric(Sys.time() - saasi_par_start)
timing_results <- rbind(timing_results,data.frame(simulation=sim_num,method = "saasi_par_est",time_seconds=saasi_par_est_time,success=saasi_par_success))

# SIMMAP
simmap_start <- Sys.time()
simmap_success <- FALSE
tryCatch({
    smap_model <- fitMk(phy, phy$tip.state, model="ARD", pi="fitzjohn")
    phy_smap <- simmap(smap_model)
    sm <- summary(phy_smap)
    simmap_lik <- as.data.frame(sm$ace)
    write.csv(simmap_lik, file=paste0(output_dir, "/simmap_", sim_num, ".csv"))
    simmap_success <- TRUE
}, error = function(e) {})
simmap_time <- as.numeric(Sys.time() - simmap_start)
timing_results <- rbind(timing_results, data.frame(simulation=sim_num, method="simmap", time_seconds=simmap_time, success=simmap_success))

# Save timing
write.csv(timing_results, file=paste0(output_dir, "/timing_", sim_num, ".csv"), row.names=FALSE)
EOF

# Run simulations
for i in $(seq 1 $N_SIMS); do
    # Create subfolder for this simulation
    SIM_DIR="$OUTPUT_DIR/sim_$i"
    mkdir -p "$SIM_DIR"
    
    # Run R methods
    Rscript "$OUTPUT_DIR/run_asr.R" $i $NTAXA $MAX_TIME "$SIM_DIR"
    
    # TreeTime
    cd "$SIM_DIR"
    mkdir -p treetime_results
    treetime_start=$(date +%s)
    treetime mugration --tree "tree_${i}.nwk" --states "metadata_${i}.csv" --attribute numeric_state --confidence --outdir "treetime_results" --verbose 0 > /dev/null 2>&1
    treetime_success=$?
    treetime_time=$(($(date +%s) - treetime_start))
    
    # PastML
    mkdir -p pastml_results
    pastml_start=$(date +%s)
    pastml --tree "tree_${i}.nwk" --data "metadata_${i}.txt" --columns numeric_state --work_dir "pastml_results" --data_sep "," > /dev/null 2>&1
    pastml_success=$?
    pastml_time=$(($(date +%s) - pastml_start))
    
    cd - > /dev/null
    
    # Save external method timing in simulation folder
    echo "simulation,method,time_seconds,success" > "$SIM_DIR/external_timing.csv"
    echo "$i,treetime,$treetime_time,$([[ $treetime_success -eq 0 ]] && echo true || echo false)" >> "$SIM_DIR/external_timing.csv"
    echo "$i,pastml,$pastml_time,$([[ $pastml_success -eq 0 ]] && echo true || echo false)" >> "$SIM_DIR/external_timing.csv"
done
