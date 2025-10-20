# Complete ASR Accuracy Analysis
library(ape)
library(ggplot2)
library(dplyr)

# Set your base directory
base_dir <- "result_baseline_small_tree/"  

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

# Initialize results dataframe
all_results <- data.frame()

# Process each simulation
for(i in seq_along(sim_dirs)) {
  sim_dir <- sim_dirs[i]
  sim_num <- sim_nums[i]
  
  cat("Processing simulation", sim_num, "\n")
  
  tryCatch({
    # Read true states and tree
    true_tree <- read.csv(file.path(sim_dir, paste0("true_tree_", sim_num, ".csv")))
    tree <- read.tree(file.path(sim_dir, paste0("tree_", sim_num, ".nwk")))
    
    # Calculate tree size
    tree_size <- tree$Nnode * 2 + 1
    
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
    }, error = function(e) cat("SAASI failed for sim", sim_num, "\n"))
    
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
    
    # TreeTime
    tryCatch({
      treetime_file <- file.path(sim_dir, "treetime_results", "confidence.csv")
      if(file.exists(treetime_file)) {
        treetime_result <- read.csv(treetime_file)
        treetime_subset <- treetime_result[grep("^(nd\\d+)$", treetime_result$X.name), ]
        ace_result <- read.csv(file.path(sim_dir, paste0("ace_", sim_num, ".csv")))
        treetime_matched <- treetime_subset[match(ace_result$X, treetime_subset$X.name), ]
        treetime_acc <- saasi_accuracy(tree, treetime_matched[2:3], true_tree)
        all_results <- rbind(all_results, data.frame(
          simulation = sim_num, method = "TreeTime", 
          accuracy = treetime_acc$acc, prob_accuracy = treetime_acc$prob_acc,
          tree_size = tree_size
        ))
      }
    }, error = function(e) cat("TreeTime failed for sim", sim_num, "\n"))
    
  }, error = function(e) {
    cat("Failed to process simulation", sim_num, ":", e$message, "\n")
  })
}

all_results <- all_results[all_results$tree_size >= 50, ]


# Save results
write.csv(all_results, file.path(base_dir, "accuracy_results.csv"), row.names = FALSE)

# Create box plots with gradient coloring based on tree size
p5 <- ggplot(all_results, aes(x = method, y = accuracy)) +
  geom_boxplot(alpha = 0.3, outlier.shape = NA) +  # Make boxplots more transparent and remove outliers
  geom_jitter(aes(color = tree_size), width = 0.2, alpha = 0.7, size = 2) +
  scale_color_gradientn(colors = rainbow(all_results$tree_size))+
  scale_x_discrete(limits = c("ACE", "PastML", "TreeTime", "SIMMAP", "SAASI", "SAASI_pars"),
                   labels = c("ACE", "PastML", "TreeTime", "SIMMAP", "SAASI", "SAASI*")) +
  labs(title = "ASR Method Accuracy Comparison",
       x = "Method", y = "Accuracy") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ylim(0, 1)

ggsave(file.path(base_dir, "accuracy_boxplot_tree_size.png"), p5, width = 10, height = 6, dpi = 300)

p6 <- ggplot(all_results, aes(x = method, y = prob_accuracy)) +
  geom_boxplot(alpha = 0.3, outlier.shape = NA) +  # Make boxplots more transparent and remove outliers
  geom_jitter(aes(color = tree_size), width = 0.2, alpha = 0.7, size = 2) +
  #scale_color_gradient(low = "lightblue", high = "darkred", name = "Tree Size") +
  scale_color_gradientn(colors = rainbow(all_results$tree_size))+
  scale_x_discrete(limits = c("ACE", "PastML", "TreeTime", "SIMMAP", "SAASI", "SAASI_pars"),
                   labels = c("ACE", "PastML", "TreeTime", "SIMMAP", "SAASI", "SAASI*")) +
  labs(title = "ASR Method Probability Accuracy Comparison",
       x = "Method", y = "Probability Accuracy") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ylim(0, 1)

ggsave(file.path(base_dir, "prob_accuracy_boxplot_tree_size.png"), p6, width = 10, height = 6, dpi = 300)

# Summary statistics including tree size
summary_stats <- all_results %>%
  group_by(method) %>%
  summarise(
    n = n(),
    mean_accuracy = mean(accuracy, na.rm = TRUE),
    sd_accuracy = sd(accuracy, na.rm = TRUE),
    mean_prob_accuracy = mean(prob_accuracy, na.rm = TRUE),
    sd_prob_accuracy = sd(prob_accuracy, na.rm = TRUE),
    mean_tree_size = mean(tree_size, na.rm = TRUE),
    sd_tree_size = sd(tree_size, na.rm = TRUE),
    min_tree_size = min(tree_size, na.rm = TRUE),
    max_tree_size = max(tree_size, na.rm = TRUE),
    .groups = 'drop'
  )

print(summary_stats)
write.csv(summary_stats, file.path(base_dir, "accuracy_summary.csv"), row.names = FALSE)

# Print tree size range for reference
cat("Tree size range:", min(all_results$tree_size, na.rm = TRUE), "-", max(all_results$tree_size, na.rm = TRUE), "\n")

cat("Analysis complete! Check files:\n")
cat("- accuracy_results.csv\n")
cat("- accuracy_boxplot_tree_size.png\n") 
cat("- prob_accuracy_boxplot_tree_size.png\n")
cat("- accuracy_summary.csv\n")