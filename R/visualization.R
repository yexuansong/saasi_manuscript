library(ape)
library(ggplot2)
library(dplyr)
library(tidyr)

base_dir <- "result_baseline/"  

# Accuracy function
saasi_accuracy <- function(tree, asr, true.states) {
  reconstructed_states <- apply(asr, 1, which.max)
  comparison_data <- data.frame(
    true = tail(true.states$State, length(tree$tip.label) - 1),
    reconstructed = reconstructed_states
  )
  
  accuracy <- mean(comparison_data$true == comparison_data$reconstructed)
  
  prob_accuracy <- mean(sapply(1:nrow(asr), function(i) {
    asr[i, true.states$State[length(tree$tip.label) + i]]
  }))
  
  return(list(
    acc = accuracy,
    prob_acc = prob_accuracy
  ))
}

# Find all simulation directories
sim_dirs <- list.dirs(base_dir, recursive = FALSE)
sim_dirs <- sim_dirs[grep("sim_\\d+", sim_dirs)]
sim_nums <- as.numeric(gsub(".*sim_(\\d+)", "\\1", basename(sim_dirs)))

all_results <- data.frame()
all_timing <- data.frame()

treetime_corrections <- c("1_0", "2_0", "5_0", "10_0")
treetime_labels <- c("1.0", "2.0", "5.0", "10.0")

for(i in seq_along(sim_dirs)) {
  sim_dir <- sim_dirs[i]
  sim_num <- sim_nums[i]
  
  cat("Processing simulation", sim_num, "\n")
  
  tryCatch({
    true_tree <- read.csv(file.path(sim_dir, paste0("true_tree_", sim_num, ".csv")))
    tree <- read.tree(file.path(sim_dir, paste0("tree_", sim_num, ".nwk")))
    
    tree_size <- tree$Nnode * 2 + 1
    
    timing_file_1 <- file.path(sim_dir, paste0("timing_", sim_num, ".csv"))
    if(file.exists(timing_file_1)) {
      timing_data_1 <- read.csv(timing_file_1)
      timing_data_1$simulation <- sim_num
      timing_data_1$tree_size <- tree_size
      all_timing <- rbind(all_timing, timing_data_1)
    }
    
    timing_file_2 <- file.path(sim_dir, "external_timing.csv")
    if(file.exists(timing_file_2)) {
      timing_data_2 <- read.csv(timing_file_2)
      timing_data_2$tree_size <- tree_size
      all_timing <- rbind(all_timing, timing_data_2)
    }
    
    # ACE
    tryCatch({
      ace_result <- read.csv(file.path(sim_dir, paste0("ace_", sim_num, ".csv")))
      ace_acc <- saasi_accuracy(tree, ace_result[2:3], true_tree)
      all_results <- rbind(all_results, data.frame(
        simulation = sim_num, method = "ACE", 
        accuracy = ace_acc$acc, prob_accuracy = ace_acc$prob_acc,
        tree_size = tree_size
      ))
    }, error = function(e) cat("ACE failed for sim", sim_num, "\n"))
    
    # SIMMAP
    tryCatch({
      simmap_result <- read.csv(file.path(sim_dir, paste0("simmap_", sim_num, ".csv")))
      simmap_subset <- head(simmap_result[2:3], length(true_tree$State) - length(tree$tip.label))
      simmap_acc <- saasi_accuracy(tree, simmap_subset, true_tree)
      all_results <- rbind(all_results, data.frame(
        simulation = sim_num, method = "SIMMAP", 
        accuracy = simmap_acc$acc, prob_accuracy = simmap_acc$prob_acc,
        tree_size = tree_size
      ))
    }, error = function(e) cat("SIMMAP failed for sim", sim_num, "\n"))
    
    # SAASI
    tryCatch({
      saasi_result <- read.csv(file.path(sim_dir, paste0("saasi_", sim_num, ".csv")))
      saasi_acc <- saasi_accuracy(tree, saasi_result[2:3], true_tree)
      all_results <- rbind(all_results, data.frame(
        simulation = sim_num, method = "SAASI", 
        accuracy = saasi_acc$acc, prob_accuracy = saasi_acc$prob_acc,
        tree_size = tree_size
      ))
    }, error = function(e) cat("SAASI failed for sim", sim_num, "\n"))
    
    # SAASI with pars est
    tryCatch({
      saasi_pars_result <- read.csv(file.path(sim_dir, paste0("saasi_pars_est_", sim_num, ".csv")))
      saasi_pars_acc <- saasi_accuracy(tree, saasi_pars_result[2:3], true_tree)
      all_results <- rbind(all_results, data.frame(
        simulation = sim_num, method = "SAASI_pars", 
        accuracy = saasi_pars_acc$acc, prob_accuracy = saasi_pars_acc$prob_acc,
        tree_size = tree_size
      ))
    }, error = function(e) cat("SAASI_pars failed for sim", sim_num, "\n"))
    
    # PastML
    tryCatch({
      pastml_file <- file.path(sim_dir, "pastml_results", 
                               "marginal_probabilities.character_numeric_state.model_F81.tab")
      if(file.exists(pastml_file)) {
        pastml_result <- read.table(pastml_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
        pastml_subset <- pastml_result[grep("^(root|nd\\d+)$", pastml_result$node), ]
        ace_result <- read.csv(file.path(sim_dir, paste0("ace_", sim_num, ".csv")))
        pastml_matched <- pastml_subset[match(ace_result$X, pastml_subset$node), ]
        pastml_acc <- saasi_accuracy(tree, pastml_matched[2:3], true_tree)
        all_results <- rbind(all_results, data.frame(
          simulation = sim_num, method = "PastML", 
          accuracy = pastml_acc$acc, prob_accuracy = pastml_acc$prob_acc,
          tree_size = tree_size
        ))
      }
    }, error = function(e) cat("PastML failed for sim", sim_num, "\n"))
    
    # TreeTime with different sampling bias corrections
    for(j in seq_along(treetime_corrections)) {
      correction_suffix <- treetime_corrections[j]
      correction_label <- treetime_labels[j]
      
      tryCatch({
        treetime_file <- file.path(sim_dir, paste0("treetime_results_", correction_suffix), "confidence.csv")
        if(file.exists(treetime_file)) {
          treetime_result <- read.csv(treetime_file)
          treetime_subset <- treetime_result[grep("^(nd\\d+)$", treetime_result$X.name), ]
          ace_result <- read.csv(file.path(sim_dir, paste0("ace_", sim_num, ".csv")))
          treetime_matched <- treetime_subset[match(ace_result$X, treetime_subset$X.name), ]
          treetime_acc <- saasi_accuracy(tree, treetime_matched[2:3], true_tree)
          all_results <- rbind(all_results, data.frame(
            simulation = sim_num, 
            method = paste0("TreeTime_", correction_label), 
            accuracy = treetime_acc$acc, 
            prob_accuracy = treetime_acc$prob_acc,
            tree_size = tree_size
          ))
        }
      }, error = function(e) cat("TreeTime", correction_label, "failed for sim", sim_num, "\n"))
    }
    
  }, error = function(e) {
    cat("Failed to process simulation", sim_num, ":", e$message, "\n")
  })
}

# Filter by tree size
all_results <- all_results[all_results$tree_size >= 50, ]
all_timing <- all_timing[all_timing$tree_size >= 50, ]
all_timing <- all_timing[all_timing$time_seconds <= 50, ]

# Process timing data - standardize method names
all_timing <- all_timing %>%
  mutate(method = case_when(
    grepl("treetime_1\\.0", method, ignore.case = TRUE) ~ "TreeTime_1.0",
    grepl("treetime_2\\.0", method, ignore.case = TRUE) ~ "TreeTime_2.0",
    grepl("treetime_5\\.0", method, ignore.case = TRUE) ~ "TreeTime_5.0",
    grepl("treetime_10\\.0", method, ignore.case = TRUE) ~ "TreeTime_10.0",
    grepl("^treetime$", method, ignore.case = TRUE) ~ "TreeTime_1.0",  # fallback for old format
    grepl("pastml", method, ignore.case = TRUE) ~ "PastML",
    grepl("ace", method, ignore.case = TRUE) ~ "ACE",
    grepl("simmap", method, ignore.case = TRUE) ~ "SIMMAP",
    grepl("est", method, ignore.case = TRUE) ~ "SAASI_pars",
    grepl("saasi", method, ignore.case = TRUE) ~ "SAASI",
    TRUE ~ method
  ))

# Save results
write.csv(all_results, file.path(base_dir, "accuracy_results.csv"), row.names = FALSE)
write.csv(all_timing, file.path(base_dir, "timing_results.csv"), row.names = FALSE)

p5 <- ggplot(all_results, aes(x = method, y = accuracy)) +
  geom_boxplot(alpha = 0.3, outlier.shape = NA) +
  geom_jitter(aes(color = tree_size), width = 0.2, alpha = 0.7, size = 2) +
  scale_color_gradientn(colors = rainbow(length(unique(all_results$tree_size)))) +
  scale_x_discrete(limits = c("ACE", "PastML", "TreeTime_1.0", "TreeTime_2.0", "TreeTime_5.0", "TreeTime_10.0", 
                              "SIMMAP", "SAASI", "SAASI_pars"),
                   labels = c("ACE", "PastML", "TreeTime\n(1.0)", "TreeTime\n(2.0)", "TreeTime\n(5.0)", "TreeTime\n(10.0)",
                              "SIMMAP", "SAASI", "SAASI*")) +
  labs(title = "ASR Method Accuracy Comparison",
       x = "Method", y = "Accuracy") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ylim(0, 1)

ggsave(file.path(base_dir, "accuracy_boxplot_tree_size.png"), p5, width = 12, height = 6, dpi = 300)

p6 <- ggplot(all_results, aes(x = method, y = prob_accuracy)) +
  geom_boxplot(alpha = 0.3, outlier.shape = NA) +
  geom_jitter(aes(color = tree_size), width = 0.2, alpha = 0.7, size = 2) +
  scale_color_gradientn(colors = rainbow(length(unique(all_results$tree_size)))) +
  scale_x_discrete(limits = c("ACE", "PastML", "TreeTime_1.0", "TreeTime_2.0", "TreeTime_5.0", "TreeTime_10.0",
                              "SIMMAP", "SAASI", "SAASI_pars"),
                   labels = c("ACE", "PastML", "TreeTime\n(1.0)", "TreeTime\n(2.0)", "TreeTime\n(5.0)", "TreeTime\n(10.0)",
                              "SIMMAP", "SAASI", "SAASI*")) +
  labs(title = "ASR Method Probability Accuracy Comparison",
       x = "Method", y = "Probability Accuracy") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ylim(0, 1)

ggsave(file.path(base_dir, "prob_accuracy_boxplot_tree_size.png"), p6, width = 12, height = 6, dpi = 300)

p7 <- ggplot(all_timing, aes(x = method, y = time_seconds)) +
  geom_boxplot(alpha = 0.3, outlier.shape = NA) +
  geom_jitter(aes(color = tree_size), width = 0.2, alpha = 0.7, size = 2) +
  scale_color_gradientn(colors = rainbow(length(unique(all_timing$tree_size)))) +
  scale_x_discrete(limits = c("ACE", "PastML", "TreeTime_1.0", "TreeTime_2.0", "TreeTime_5.0", "TreeTime_10.0",
                              "SIMMAP", "SAASI", "SAASI_pars"),
                   labels = c("ACE", "PastML", "TreeTime\n(1.0)", "TreeTime\n(2.0)", "TreeTime\n(5.0)", "TreeTime\n(10.0)",
                              "SIMMAP", "SAASI", "SAASI*")) +
  labs(title = "ASR Method Running Time Comparison",
       x = "Method", y = "Time (seconds)") +
  theme_minimal() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(file.path(base_dir, "timing_boxplot.png"), p7, width = 12, height = 6, dpi = 300)
