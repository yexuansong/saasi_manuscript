# Complete ASR Accuracy Analysis - Multiple Directories Combined
library(ape)
library(ggplot2)
library(dplyr)

# Define your base directories
base_dirs <- c(
  "result_baseline_small_tree_subsample/"
)


saasi_accuracy_downsample <- function(tree,true_tree, true_tree_before_downsample,
                                      true_tree_before_downsample_csv,asr, true.states) {
  diff <- get_diff_hist(tree, true_tree,true_tree_before_downsample
                        ,true_tree_before_downsample_csv)
  reconstructed_states <- apply(asr, 1, which.max)
  comparison_data <- data.frame(
    true = tail(true.states$State, length(tree$tip.label) - 1),
    reconstructed = reconstructed_states
  )
  
  accuracy <- mean(comparison_data$true == comparison_data$reconstructed)
  
  prob_accuracy <- mean(sapply(1:nrow(asr), function(i) {
    asr[i, true.states$State[length(tree$tip.label) + i]]
  }))
  correct_count <- sum(comparison_data$true == comparison_data$reconstructed)
  accuracy_downsample <- correct_count / (length(tree$tip.label) - 1 + diff)
  
  prob_sum <- sum(sapply(1:nrow(asr), function(i) {
    asr[i, true.states$State[length(tree$tip.label) + i]]
  }))
  prob_accuracy_downsample <- prob_sum / (nrow(asr) + diff)
  
  return(list(
    acc = accuracy,
    prob_acc = prob_accuracy,
    acc_down = accuracy_downsample,
    prob_acc_down = prob_accuracy_downsample
  ))
}

# Initialize combined results dataframe
all_results <- data.frame()

# Loop through each base directory
for(base_idx in seq_along(base_dirs)) {
  base_dir <- base_dirs[base_idx]
  scenario <- scenario_names[base_idx]
  
  cat("\n=== Processing scenario:", scenario, "===\n")
  
  # Find all simulation directories
  sim_dirs <- list.dirs(base_dir, recursive = FALSE)
  sim_dirs <- sim_dirs[grep("sim_\\d+", sim_dirs)]
  sim_nums <- as.numeric(gsub(".*sim_(\\d+)", "\\1", basename(sim_dirs)))
  
  # Process each simulation
  for(i in seq_along(sim_dirs)) {
    sim_dir <- sim_dirs[i]
    sim_num <- sim_nums[i]
    
    cat("Processing simulation", sim_num, "\n")
    
    ####
    
    tryCatch({
      # Read true states and tree
      tree <- read.tree(file.path(sim_dir, paste0("tree_", sim_num, ".nwk")))
      true_tree <- read.csv(file.path(sim_dir, paste0("true_tree_", sim_num, ".csv")))
      true_tree_before_downsample <-  read.tree(file.path(sim_dir, paste0("tree_before_ds_", sim_num, ".nwk")))
      true_tree_before_downsample_csv <-  read.csv(file.path(sim_dir, paste0("true_tree_before_ds_", sim_num, ".csv")))

      # Calculate tree size
      tree_size <- tree$Nnode * 2 + 1
      
      # ACE
      tryCatch({
        ace_result <- read.csv(file.path(sim_dir, paste0("ace_", sim_num, ".csv")))
        ace_acc <- saasi_accuracy_downsample(tree,
                                             true_tree,
                                             true_tree_before_downsample,
                                             true_tree_before_downsample_csv,
                                             ace_result[2:3], true_tree)
        all_results <- rbind(all_results, data.frame(
          scenario = scenario,
          simulation = sim_num, 
          method = "ACE", 
          accuracy = ace_acc$acc, 
          prob_accuracy = ace_acc$prob_acc,
          accuracy_down = ace_acc$acc_down,
          prob_accuracy_down = ace_acc$prob_acc_down,
          tree_size = tree_size
        ))
      }, error = function(e) cat("ACE failed for sim", sim_num, "\n"))
      
      # SIMMAP
      tryCatch({
        simmap_result <- read.csv(file.path(sim_dir, paste0("simmap_", sim_num, ".csv")))
        simmap_subset <- head(simmap_result[2:3], length(true_tree$State) - length(tree$tip.label))
        simmap_acc <- saasi_accuracy_downsample(tree,
                                                true_tree,
                                                true_tree_before_downsample,
                                                true_tree_before_downsample_csv,
                                                simmap_subset, true_tree)
        all_results <- rbind(all_results, data.frame(
          scenario = scenario,
          simulation = sim_num, 
          method = "SIMMAP", 
          accuracy = simmap_acc$acc, 
          prob_accuracy = simmap_acc$prob_acc,
          accuracy_down = simmap_acc$acc_down,
          prob_accuracy_down = simmap_acc$prob_acc_down,
          tree_size = tree_size
        ))
      }, error = function(e) cat("SIMMAP failed for sim", sim_num, "\n"))
      
      # SAASI
      tryCatch({
        saasi_result <- read.csv(file.path(sim_dir, paste0("saasi_", sim_num, ".csv")))
        saasi_acc <- saasi_accuracy_downsample(tree,
                                               true_tree,
                                               true_tree_before_downsample,
                                               true_tree_before_downsample_csv,
                                               saasi_result[2:3], true_tree)
        all_results <- rbind(all_results, data.frame(
          scenario = scenario,
          simulation = sim_num, 
          method = "SAASI", 
          accuracy = saasi_acc$acc, 
          prob_accuracy = saasi_acc$prob_acc,
          accuracy_down = saasi_acc$acc_down,
          prob_accuracy_down = saasi_acc$prob_acc_down,
          tree_size = tree_size
        ))
      }, error = function(e) cat("SAASI failed for sim", sim_num, "\n"))
      
      # SAASI with pars est
      tryCatch({
        saasi_pars_result <- read.csv(file.path(sim_dir, paste0("saasi_pars_est_", sim_num, ".csv")))
        saasi_pars_acc <- saasi_accuracy_downsample(tree,
                                                    true_tree,
                                                    true_tree_before_downsample,
                                                    true_tree_before_downsample_csv,
                                                    saasi_pars_result[2:3], true_tree)
        all_results <- rbind(all_results, data.frame(
          scenario = scenario,
          simulation = sim_num, 
          method = "SAASI_pars", 
          accuracy = saasi_pars_acc$acc, 
          prob_accuracy = saasi_pars_acc$prob_acc,
          accuracy_down = saasi_pars_acc$acc_down,
          prob_accuracy_down = saasi_pars_acc$prob_acc_down,
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
          pastml_acc <- saasi_accuracy_downsample(tree,
                                                  true_tree,
                                                  true_tree_before_downsample,
                                                  true_tree_before_downsample_csv,
                                                  pastml_matched[2:3], true_tree)
          all_results <- rbind(all_results, data.frame(
            scenario = scenario,
            simulation = sim_num, 
            method = "PastML", 
            accuracy = pastml_acc$acc, 
            prob_accuracy = pastml_acc$prob_acc,
            accuracy_down = pastml_acc$acc_down,
            prob_accuracy_down = pastml_acc$prob_acc_down,
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
          treetime_acc <- saasi_accuracy_downsample(tree,
                                                    true_tree,
                                                    true_tree_before_downsample,
                                                    true_tree_before_downsample_csv,
                                                    treetime_matched[2:3], true_tree)
          all_results <- rbind(all_results, data.frame(
            scenario = scenario,
            simulation = sim_num, 
            method = "TreeTime", 
            accuracy = treetime_acc$acc, 
            prob_accuracy = treetime_acc$prob_acc,
            accuracy_down = treetime_acc$acc_down,
            prob_accuracy_down = treetime_acc$prob_acc_down,
            tree_size = tree_size
          ))
        }
      }, error = function(e) cat("TreeTime failed for sim", sim_num, "\n"))
      
    }, error = function(e) {
      cat("Failed to process simulation", sim_num, ":", e$message, "\n")
    })
  }
}

# Optional: Filter by tree size if needed
all_results <- all_results[all_results$tree_size >= 50, ]

# Save combined results
write.csv(all_results, "combined_accuracy_results.csv", row.names = FALSE)


# Create box plots with gradient coloring based on tree size
p5 <- ggplot(all_results, aes(x = method, y = accuracy)) +
  geom_boxplot(alpha = 0.3, outlier.shape = NA) +  # Make boxplots more transparent and remove outliers
  geom_jitter(aes(color = tree_size), width = 0.2, alpha = 0.7, size = 2) +
  scale_color_gradient(low = "blue", high = "red", name = "Tree Size") +
  #scale_color_gradientn(colors = rainbow(all_results$tree_size))+
  labs(title = "ASR Method Accuracy Comparison (Combined)",
       x = "Method", y = "Accuracy") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ylim(0, 1)


ggplot(all_results, aes(x = method)) +
  geom_boxplot(aes(y = accuracy), alpha = 0.3, outlier.shape = NA, fill = "lightblue") +
  geom_boxplot(aes(y = accuracy_down), alpha = 0.3, outlier.shape = NA, fill = "lightcoral") +
  geom_jitter(aes(y = accuracy, color = tree_size), width = 0.2, alpha = 0.7, size = 2) +
  geom_jitter(aes(y = accuracy_down, color = tree_size), width = 0.2, alpha = 0.7, size = 2, shape = 17) +
  scale_color_gradient(low = "blue", high = "red", name = "Tree Size") +
  labs(title = "ASR Method Accuracy Comparison (Combined)",
       x = "Method", y = "Accuracy") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ylim(0, 1)


ggplot(all_results, aes(x = method)) +
  geom_boxplot(aes(y = accuracy), alpha = 0.3, outlier.shape = NA, fill = "lightblue") +
  geom_boxplot(data = subset(all_results, !method %in% c("SAASI", "SAASI_pars")),
               aes(y = accuracy_down), alpha = 0.3, outlier.shape = NA, fill = "lightcoral") +
  geom_jitter(aes(y = accuracy, color = tree_size), width = 0.2, alpha = 0.7, size = 2) +
  geom_jitter(data = subset(all_results, !method %in% c("SAASI", "SAASI_pars")),
              aes(y = accuracy_down, color = tree_size), width = 0.2, alpha = 0.7, size = 2, shape = 17) +
  scale_color_gradient(low = "blue", high = "red", name = "Tree Size") +
  labs(title = "ASR Method Accuracy Comparison (Combined)",
       x = "Method", y = "Accuracy") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ylim(0, 1)


ggplot(all_results, aes(x = method)) +
  geom_violin(aes(y = accuracy), alpha = 0.3, fill = "lightblue") +
  geom_violin(data = subset(all_results, !method %in% c("SAASI", "SAASI_pars")),
              aes(y = accuracy_down), alpha = 0.3, fill = "lightcoral") +
  geom_jitter(aes(y = accuracy, color = tree_size), width = 0.2, alpha = 0.7, size = 2) +
  geom_jitter(data = subset(all_results, !method %in% c("SAASI", "SAASI_pars")),
              aes(y = accuracy_down, color = tree_size), width = 0.2, alpha = 0.7, size = 2, shape = 17) +
  scale_color_gradient(low = "blue", high = "red", name = "Tree Size") +
  labs(title = "ASR Method Accuracy Comparison",
       x = "Method", y = "Accuracy") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ylim(0, 1)

ggplot(all_results, aes(x = method)) +
  geom_violin(aes(y = accuracy, fill = "Without Downsample"), alpha = 0.3) +
  geom_violin(data = subset(all_results, !method %in% c("SAASI", "SAASI_pars")),
              aes(y = accuracy_down, fill = "With Downsample"), alpha = 0.3) +
  scale_fill_manual(values = c("Without Downsample" = "lightblue", 
                               "With Downsample" = "lightcoral"),
                    name = "Type") +
  scale_x_discrete(limits = c("ACE", "PastML", "TreeTime", "SIMMAP", "SAASI", "SAASI_pars"),
                   labels = c("ACE", "PastML", "TreeTime", "SIMMAP", "SAASI", "SAASI*")) +
  scale_y_continuous(limits = c(0, 1), sec.axis = dup_axis(name = NULL)) +
  labs(title = "ASR Method Accuracy Comparison",
       x = "Method", y = "Accuracy") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))




ggplot(all_results, aes(x = method)) +
  geom_boxplot(aes(y = accuracy, fill = "Without Downsample"), alpha = 0.3, outlier.shape = NA) +
  geom_boxplot(data = subset(all_results, !method %in% c("SAASI", "SAASI_pars")),
               aes(y = accuracy_down, fill = "With Downsample"), alpha = 0.3, outlier.shape = NA) +
  scale_fill_manual(values = c("Without Downsample" = "lightblue", 
                               "With Downsample" = "lightcoral"),
                    name = "Type") +
  scale_x_discrete(limits = c("ACE", "PastML", "TreeTime", "SIMMAP", "SAASI", "SAASI_pars"),
                   labels = c("ACE", "PastML", "TreeTime", "SIMMAP", "SAASI", "SAASI*")) +
  scale_y_continuous(limits = c(0, 1), sec.axis = dup_axis(name = NULL)) +
  labs(title = "ASR Method Accuracy Comparison (Combined)",
       x = "Method", y = "Accuracy") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


ggplot(all_results, aes(x = method)) +
  geom_boxplot(aes(y = accuracy, fill = "Without Downsample"), alpha = 0.3, outlier.shape = NA) +
  geom_boxplot(data = subset(all_results, !method %in% c("SAASI", "SAASI_pars")),
               aes(y = accuracy_down, fill = "With Downsample"), alpha = 0.3, outlier.shape = NA) +
  geom_jitter(aes(y = accuracy, color = tree_size), width = 0.2, alpha = 0.7, size = 2) +
  geom_jitter(data = subset(all_results, !method %in% c("SAASI", "SAASI_pars")),
              aes(y = accuracy_down, color = tree_size), width = 0.2, alpha = 0.7, size = 2, shape = 17) +
  scale_fill_manual(values = c("Without Downsample" = "lightblue", 
                               "With Downsample" = "lightcoral"),
                    name = "Type") +
  #scale_color_gradient(low = "blue", high = "red", name = "Tree Size") +
  scale_color_gradientn(colors = rainbow(all_results$tree_size))+
  scale_x_discrete(limits = c("ACE", "PastML", "TreeTime", "SIMMAP", "SAASI", "SAASI_pars"),
                   labels = c("ACE", "PastML", "TreeTime", "SIMMAP", "SAASI", "SAASI*")) +
  scale_y_continuous(limits = c(0, 1), sec.axis = dup_axis(name = NULL)) +
  labs(title = "ASR Method Accuracy Comparison",
       x = "Method", y = "Accuracy") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))



ggplot(all_results, aes(x = method)) +
  geom_violin(aes(y = accuracy, fill = "Without Downsample"), alpha = 0.3) +
  geom_violin(data = subset(all_results, !method %in% c("SAASI", "SAASI_pars")),
              aes(y = accuracy_down, fill = "With Downsample"), alpha = 0.3) +
  scale_fill_manual(values = c("Without Downsample" = "lightblue", 
                               "With Downsample" = "lightcoral"),
                    name = "") +
  labs(title = "ASR Method Accuracy Comparison",
       x = "Method", y = "Accuracy") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ylim(0, 1)

ggsave("combined_accuracy_boxplot_tree_size.png", p5, width = 10, height = 6, dpi = 300)

p6 <- ggplot(all_results, aes(x = method, y = prob_accuracy)) +
  geom_boxplot(alpha = 0.3, outlier.shape = NA) +  # Make boxplots more transparent and remove outliers
  geom_jitter(aes(color = tree_size), width = 0.2, alpha = 0.7, size = 2) +
  #scale_color_gradient(low = "lightblue", high = "darkred", name = "Tree Size") +
  scale_color_gradientn(colors = rainbow(all_results$tree_size))+
  labs(title = "ASR Method Probability Accuracy Comparison (Combined)",
       x = "Method", y = "Probability Accuracy") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ylim(0, 1)

ggsave("combined_prob_accuracy_boxplot_tree_size.png", p6, width = 10, height = 6, dpi = 300)

# Faceted by scenario with gradient
p7 <- ggplot(all_results, aes(x = method, y = accuracy)) +
  geom_boxplot(alpha = 0.3, outlier.shape = NA) +
  geom_jitter(aes(color = tree_size), width = 0.2, alpha = 0.7, size = 2) +
  #scale_color_gradient(low = "blue", high = "red", name = "Tree Size") +
  #scale_color_gradientn(colors = rainbow(all_results$tree_size))+
  scale_color_viridis_c(tree_size)+
  facet_wrap(~scenario) +
  labs(title = "ASR Method Accuracy by Scenario",
       x = "Method", y = "Accuracy") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ylim(0, 1)

ggsave("combined_accuracy_by_scenario_tree_size.png", p7, width = 12, height = 6, dpi = 300)

p8 <- ggplot(all_results, aes(x = method, y = prob_accuracy)) +
  geom_boxplot(alpha = 0.3, outlier.shape = NA) +
  geom_jitter(aes(color = tree_size), width = 0.2, alpha = 0.7, size = 2) +
  #scale_color_gradient(low = "lightblue", high = "darkred", name = "Tree Size") +
  scale_color_gradientn(colors = rainbow(all_results$tree_size))+
  facet_wrap(~scenario) +
  labs(title = "ASR Method Probability Accuracy by Scenario",
       x = "Method", y = "Probability Accuracy") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ylim(0, 1)
ggsave("combined_prob_accuracy_by_scenario_tree_size.png", p8, width = 12, height = 6, dpi = 300)



#new_df <- read.csv("combined_accuracy_results.csv")
# p9 <- ggplot(new_df,aes(x = tree_size)) + geom_histogram()
# new_df_filter <- new_df %>% filter(tree_size > 15)
# p10 <- ggplot(new_df_filter,aes(x = tree_size)) + geom_histogram()




# p7 <- ggplot(new_df_filter, aes(x = method, y = accuracy)) +
#   geom_boxplot(alpha = 0.3, outlier.shape = NA) +
#   geom_jitter(aes(color = tree_size), width = 0.2, alpha = 0.7, size = 2) +
#   scale_color_viridis_c(tree_size)+
#   facet_wrap(~scenario) +
#   labs(title = "ASR Method Accuracy by Scenario",
#        x = "Method", y = "Accuracy") +
#   theme_minimal() +
#   theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
#   ylim(0, 1)



# Summary statistics including tree size by scenario
summary_stats <- all_results %>%
  group_by(scenario, method) %>%
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
write.csv(summary_stats, "combined_accuracy_summary.csv", row.names = FALSE)

# Overall summary across all scenarios
overall_summary <- all_results %>%
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

print("\nOverall Summary:")
print(overall_summary)
write.csv(overall_summary, "overall_accuracy_summary.csv", row.names = FALSE)

# Print tree size range for reference
cat("\nTree size range:", min(all_results$tree_size, na.rm = TRUE), "-", max(all_results$tree_size, na.rm = TRUE), "\n")

cat("\nAnalysis complete! Check files:\n")
cat("- combined_accuracy_results.csv\n")
cat("- combined_accuracy_boxplot_tree_size.png\n") 
cat("- combined_prob_accuracy_boxplot_tree_size.png\n")
cat("- combined_accuracy_by_scenario_tree_size.png\n")
cat("- combined_prob_accuracy_by_scenario_tree_size.png\n")
cat("- combined_accuracy_summary.csv\n")
cat("- overall_accuracy_summary.csv\n")