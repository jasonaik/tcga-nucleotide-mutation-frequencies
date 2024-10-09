BiocManager::install("PoisonAlien/TCGAmutations")

library(TCGAmutations)
library(viridis)
library(tidytext)
library(tidyverse)
library(dplyr)
library(data.table)
library(parallel)
library(ggplot2)
library(parallel)
library(fastmatch)
library(Rcpp)
library(RColorBrewer)
library(cowplot)
library(grid)
library(gridExtra)
library(glue)
library(scales)

today_date <- as.character(Sys.Date())
plot_dir <- file.path("plots", today_date)

if (!file.exists(plot_dir)) {
  dir.create(plot_dir, recursive = TRUE)
}

library(dplyr)
library(tidyverse)

cancer_data <- as.data.frame(tcga_available())
cancer_dfs <- tibble(
  cancer_name = character(),
  data = list()
)

cancer_dfs_facet <- tibble(
  cancer_name = character(),
  data = list()
)

for (i in seq_along(cancer_data$Study_Abbreviation)) {
  # Load TCGA mutation data
  cancer_type <- cancer_data$Study_Abbreviation[i]
  if (cancer_type == "CCLE_2024Q2") {
    break
  }
#   cancer_name <- gsub("_", " ", cancer_data$Study_Name[i])
  mutations <- tcga_load(study = cancer_type)
  mutation_data <- as.data.frame(mutations@data)
  cancer_name <- cancer_type
  nucleotide_changes <- mutation_data[, c("Hugo_Symbol", "HGVSc")]

  tumor_samples <- mutation_data[, c("Tumor_Sample_Barcode")]

  num_tumors <- length(unique(tumor_samples))

  # group_by HGVSc and add frequency using mutate, then ungroup
  df <- nucleotide_changes %>%
    mutate(mutation = paste(Hugo_Symbol, HGVSc, sep = ", ")) %>%
    group_by(mutation) %>%
    mutate(freq = n()) %>%
    ungroup()

  # sort in descending order, distinct removes duplicate rows
  df_sorted <- df %>%
    arrange(desc(freq)) %>%
    distinct(mutation, .keep_all = TRUE)

  df_merged <- df_sorted %>%
    select(mutation, freq) %>%
    mutate(cancer_name_internal = cancer_name) %>%
    mutate(cancer_type = cancer_type)

  df_merged_trunc <- df_merged[1:50, ]
  df_merged_trunc_facet <- df_merged[1:10, ]

  df_merged_trunc <- df_merged_trunc %>% mutate(mutation = str_trunc(mutation, 20))
  df_merged_trunc_facet <- df_merged_trunc_facet %>% mutate(mutation = str_trunc(mutation, 20))

  # Convert the mutation column to a factor
  df_merged_trunc$freq <- df_merged_trunc$freq / num_tumors

  df_merged_trunc_facet$freq <- df_merged_trunc_facet$freq / num_tumors

  # append dfs to list
  cancer_dfs <- cancer_dfs %>% add_row(cancer_name, data = list(df_merged_trunc))
  cancer_dfs_facet <- cancer_dfs_facet %>% add_row(cancer_name, data = list(df_merged_trunc_facet))
}


library(ggplot2)
library(purrr)
library(gridExtra)

# Create the bar plot
create_plot <- function(df, name) {
  ggplot(data = df, aes(x = mutation, y = freq)) + # nolint: object_usage_linter.
    geom_bar(stat = "identity", fill = "blue") +
    geom_text(aes(label = round(freq, 3)), size = 3.5, hjust = 1.2, color = "white") +
    coord_flip() +
    labs(x = "Mutation", y = "Frequency", title = paste("Frequency of Mutations in ", name, sep = "")) +
    theme_minimal()
}

plots <- pmap(cancer_dfs, function(cancer_name, data) {
  create_plot(data, cancer_name)
})

plots


library(tidyverse)
library(tidytext)
library(glue)

# Split df into 4 because too big for one facet plot
cancer_types <- unique(cancer_dfs_facet$cancer_name)
plots_per_group <- 12
n_groups <- ceiling(length(cancer_types) / plots_per_group)
group_assignments <- rep(1:n_groups, each = plots_per_group)[seq_along(cancer_types)]
cancer_type_groups <- split(cancer_types, group_assignments)

create_facet_plot <- function(df_data, cancer_types_subset) {
  subset_data <- df_data %>% filter(cancer_name %in% cancer_types_subset) # nolint: object_usage_linter.

  list_of_tibbles <- subset_data$data
  combined_data <- bind_rows(list_of_tibbles)

  # Create plot
  ggplot(combined_data, aes(x = reorder_within(mutation, freq, cancer_type), y = freq)) + # nolint: object_usage_linter.
    geom_bar(stat = "identity", fill = "blue", width = 0.9) +
    geom_text(aes(label = round(freq, 3)), size = 1.5, hjust = 1.2, color = "white") +
    coord_flip() +
    labs(x = "Mutation", y = "Frequency") +
    facet_wrap(~cancer_type, scales = "free") +
    theme_minimal() +
    scale_x_reordered() +
    theme(
      axis.text.y = element_text(size = 7),
      axis.text.x = element_text(size = 5)
    )
}

# Create facet plots for each group of cancer types
facet_plots <- lapply(cancer_type_groups, function(cancer_types_subset) {
  create_facet_plot(cancer_dfs_facet, cancer_types_subset)
})


for (i in seq_along(facet_plots)) {
  pdf(glue("plots/facet_{i}.pdf"))
  print(facet_plots[[i]])
  dev.off()
}

print(facet_plots[[i]])


library(tidytext)

list_of_tibbles <- cancer_dfs_facet$data
combined_data <- bind_rows(list_of_tibbles)

# Create plot
one_plot <- ggplot(combined_data, aes(x = reorder_within(mutation, freq, cancer_type), y = freq)) + # nolint: object_usage_linter.
  geom_bar(stat = "identity", fill = "blue", width = 0.7) +
  coord_flip() +
  labs(x = "Mutation", y = "Frequency") +
  theme_minimal() +
  theme(
    axis.text.y = element_blank(),
    axis.ticks.x = element_blank(),
    axis.ticks.y = element_blank(),
    panel.grid.major.y = element_blank()
  ) +
  facet_wrap(~cancer_type, scales = "free_y")
pdf("plots/one_plot.pdf")

one_plot
dev.off()


library(tidytext)
library(data.table)

list_of_tibbles <- cancer_dfs_facet$data
combined_data <- bind_rows(list_of_tibbles)

# # Make dt with top ten most frequent genes
all_dt <- as.data.table(combined_data)
all_dt[, gene := tstrsplit(mutation, split = ",")[[1]]]

top_num <- 10

topten <- all_dt[, .(n = .N, mean_freq = mean(freq)), by = gene][order(-n)][1:top_num]

# Generate colors
color_pal <- rainbow(10)

topten <- topten[, bar_color := color_pal]

all_dt[, bar_color := ifelse(gene %in% topten$gene, gene, "No")]

color_pal <- append(color_pal, "#ffffffef")

# Create plot
top_ten_oneplot <- ggplot(all_dt, aes(x = reorder_within(mutation, freq, cancer_type), y = freq, fill = bar_color)) + # nolint: object_usage_linter.
  geom_bar(stat = "identity", width = 0.7) +
  coord_flip() +
  labs(x = "Mutation", y = "Frequency") +
  theme_minimal() +
  theme(
    axis.text.y = element_blank(),
    axis.ticks.x = element_blank(),
    axis.ticks.y = element_blank(),
    panel.grid.major.y = element_blank(),
    legend.title = element_blank(),
  ) +
  facet_wrap(~cancer_type, scales = "free_y") +
  scale_fill_manual(values = color_pal, breaks = c(topten$gene))


pdf(glue("plots/{Sys.Date()}/top_ten_oneplot.pdf"))

top_ten_oneplot
dev.off()


greedy_panel_curator <- function(patient_data, top_number) {
  panel <- data.table()

  patient_data[, covered := FALSE]

  mutation_coverage <- patient_data[, .(count = .N), by = mutation]

  num <- 0
  while (num != top_number) {
    # Get the best mutation
    # Make sure that all patients covered by mutation become TRUE
    best_mutations <- mutation_coverage[count == max(count), mutation]
    best_mutation <- sample(best_mutations, 1)

    # Note that the patient_data is not unique, i.e. there are multiple instances of the same patient bound to other mutations - find them
    patients_covered <- patient_data[mutation == best_mutation, patient]
    patient_data[patient %in% patients_covered, covered := TRUE]

    num <- num + 1
    new_panel_row <- data.table(number = num, panel_mutation = best_mutation, num_covered = length(patients_covered))
    panel <- rbind(panel, new_panel_row)

    # patient_data[mutation == best_mutation, covered := TRUE]
    mutation_coverage <- patient_data[covered == FALSE, .(count = .N), by = mutation]
    if (nrow(mutation_coverage) < 1) break
  }

  return(panel)
}


# Liver HCC 60% - LIHC
# Kidney RCC 5% - KICH, KIRC
# Skin-Melanoma 69% - SKCM
# Head-SCC 33% - HNSC
# Thy AdenoCA 23% - TGCT
# Lung SCC 6% - LUSC
# CNS-GBM 66% - GBM
# Lung AdenoCA 16% - LUAD
# Bladder 74% - BLCA
# Biliary-AdenoCA 9% - CHOL

# MAKE THIS A FUNCTION, RETURN DATA; GIVE THE FUNCTION A SELECTOR FOR NORMAL OR TERT

add_tert <- function(patient_mutation, add_yes) {
  if (add_yes) {
    tert_dt <- data.table(
      cancer_type = c("LIHC", "KICH", "KIRC", "SKCM", "HNSC", "TGCT", "LUSC", "GBM", "LUAD", "BLCA", "CHOL"),
      fraction = c(0.6, 0.05, 0.05, 0.69, 0.33, 0.23, 0.06, 0.66, 0.16, 0.74, 0.09)
    )

    # Function to perform the operation for each row
    process_row <- function(cancer_type_int, fraction, patient_mutation) {
      current_cancer_pm <- patient_mutation[cancer_type == cancer_type_int]
      current_unique_patients <- unique(current_cancer_pm[, patient])
      num_unique_patients <- length(current_unique_patients)
      spike_num <- round(num_unique_patients * fraction)

      set.seed(100)
      lucky_patients <- sample(current_unique_patients, size = spike_num)

      new_patient_mutation_row <- data.table(patient = lucky_patients, mutation = "TERT, c.1-124C>T", cancer_type = cancer_type_int)
      return(new_patient_mutation_row)
    }


    new_patient_mutation <- copy(patient_mutation)

    # Apply the function to each row of tert_dt
    new_patient_mutations <- rbindlist(lapply(seq_len(nrow(tert_dt)), function(i) {
      process_row(tert_dt$cancer_type[i], tert_dt$fraction[i], patient_mutation)
    }))

    # Combine the original patient_mutation with the new rows
    new_patient_mutation <- rbind(patient_mutation, new_patient_mutations)

    # for (i in seq_len(nrow(tert_dt))) {
    #   new_row <- process_row(tert_dt$cancer_type[i], tert_dt$fraction[i], patient_mutation)
    #   new_patient_mutation <- rbind(new_patient_mutation, new_row)
    # }

    return(new_patient_mutation)
  } else {
    return(patient_mutation)
  }
}


cbind_fill <- function(..., fill = NA) {
  # Get the list of data tables
  nm <- list(...)

  # Convert each to a matrix and find the maximum number of rows
  nm <- lapply(nm, function(x) if (nrow(x) == 0) matrix(fill, 0, ncol(x)) else as.matrix(x))
  n <- max(sapply(nm, nrow))

  # Fill each matrix with NA if it has fewer rows than the max
  filled <- lapply(nm, function(x) {
    if (nrow(x) < n) {
      fill_matrix <- matrix(fill, n - nrow(x), ncol(x))
      x <- rbind(x, fill_matrix)
    }
    return(x)
  })

  # Combine the matrices side by side
  result <- do.call(cbind, filled)

  # Convert back to a data.table
  return(as.data.table(result))
}

prepare_for_cow <- function(plot) {
  plot + theme(
    legend.position = "none",
    axis.title = element_blank(),
    plot.title = element_text(size = 10),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.y = element_blank()
  )
}

round_up_custom <- function(x) {
  if (x <= 10) {
    return(10)
  } else if (x <= 20) {
    return(20)
  } else if (x <= 50) {
    return(50)
  } else {
    return(100)
  }
}


# Cancer specific
test_and_plot_cs <- function(panel, patient_data, plot_title, filename, brewer_palette, pan_can_yes) {
  # Get coverages
  panel_coverage <- data.table()

  for (i in seq_along(cancer_data$Study_Abbreviation)) {
    current_cancer_type <- cancer_data$Study_Abbreviation[i]
    specific_patients <- patient_data[cancer_type == current_cancer_type]
    patient_mutation_list <- patient_data[, .(mutations = list(mutation)), by = .(patient, cancer_type)]
    unique_patient_data <- patient_mutation_list[cancer_type == current_cancer_type]

    unlisted_panel <- unlist(panel[cancer_type == current_cancer_type, panel])

    for (num in top_nums) {
      current_panel <- unlisted_panel[1:num]
      cancer_specific_hits <- unique_patient_data[, .(
        patient,
        mutations,
        top_number = num,
        hits = sapply(mutations, function(x) length(intersect(x, current_panel)))
      )]
      total_num_patients <- cancer_specific_hits[, .N]
      hit_fraction <- cancer_specific_hits[hits > 0, .N] / total_num_patients * 100

      new_row <- data.table(cancer_type = current_cancer_type, coverage = hit_fraction, top_number = num)
      panel_coverage <- rbind(panel_coverage, new_row)
    }
  }

  if (pan_can_yes) {
    current_cancer_type <- "PAN-CAN"
    specific_patients <- copy(patient_data)
    patient_mutation_list <- patient_data[, .(mutations = list(mutation)), by = .(patient, cancer_type)]
    unique_patient_data <- copy(patient_mutation_list)

    unlisted_panel <- unlist(panel[cancer_type == current_cancer_type, panel])

    for (num in top_nums) {
      current_panel <- unlisted_panel[1:num]

      cancer_specific_hits <- unique_patient_data[, .(
        patient,
        mutations,
        top_number = num,
        hits = sapply(mutations, function(x) length(intersect(x, current_panel)))
      )]

      total_num_patients <- cancer_specific_hits[, .N]
      hit_fraction <- cancer_specific_hits[hits > 0, .N] / total_num_patients * 100

      new_row <- data.table(cancer_type = "PANCAN", coverage = hit_fraction, top_number = num)
      panel_coverage <- rbind(panel_coverage, new_row)
    }
  }


  # Plot
  panel_coverage <- panel_coverage[order(-coverage, -top_number)]

  panel_coverage[, cancer_type := factor(cancer_type, levels = unique(cancer_type))]
  panel_coverage[, top_number := factor(top_number, levels = unique(top_number))]

  plot <- ggplot(panel_coverage, aes(x = cancer_type, y = coverage, fill = top_number)) +
    geom_bar(stat = "identity", position = "identity", width = 0.9) +
    labs(title = plot_title, x = "Cancer Type", y = "Percentage Coverage", fill = "Number of\nTop Mutations") +
    theme_minimal() +
    theme(
      axis.text.y = element_text(size = 10),
      axis.text.x = element_text(size = 10, angle = 90, vjust = 0.25, hjust = 1),
    ) +
    scale_y_continuous(breaks = seq(0, 100, by = 10), limits = c(0, 100)) +
    scale_fill_brewer(palette = brewer_palette)

  pdf(glue("plots/{Sys.Date()}/{filename}.pdf"))
  print(plot)
  dev.off()

  return(list(plot, panel_coverage))
}



# Pan-cancer
test_and_plot_pc <- function(panel, patient_data, plot_title, filename, brewer_palette) {
  # Get coverages
  panel_coverage <- data.table()

  unlisted_panel <- unlist(panel[, panel])

  for (i in seq_along(cancer_data$Study_Abbreviation)) {
    current_cancer_type <- cancer_data$Study_Abbreviation[i]
    specific_patients <- patient_data[cancer_type == current_cancer_type]
    patient_mutation_list <- patient_data[, .(mutations = list(mutation)), by = .(patient, cancer_type)]
    unique_patient_data <- patient_mutation_list[cancer_type == current_cancer_type]

    for (num in top_nums) {
      current_panel <- unlisted_panel[1:num]
      cancer_specific_hits <- unique_patient_data[, .(
        patient,
        mutations,
        top_number = num,
        hits = sapply(mutations, function(x) length(intersect(x, current_panel)))
      )]
      total_num_patients <- cancer_specific_hits[, .N]
      hit_fraction <- cancer_specific_hits[hits > 0, .N] / total_num_patients * 100

      new_row <- data.table(cancer_type = current_cancer_type, coverage = hit_fraction, top_number = num)
      panel_coverage <- rbind(panel_coverage, new_row)
    }
  }

  # For pan-can
  current_cancer_type <- "PAN-CAN"
  specific_patients <- copy(patient_data)
  patient_mutation_list <- patient_data[, .(mutations = list(mutation)), by = .(patient, cancer_type)]
  unique_patient_data <- copy(patient_mutation_list)

  for (num in top_nums) {
    current_panel <- unlisted_panel[1:num]

    cancer_specific_hits <- unique_patient_data[, .(
      patient,
      mutations,
      top_number = num,
      hits = sapply(mutations, function(x) length(intersect(x, current_panel)))
    )]

    total_num_patients <- cancer_specific_hits[, .N]
    hit_fraction <- cancer_specific_hits[hits > 0, .N] / total_num_patients * 100

    new_row <- data.table(cancer_type = "PANCAN", coverage = hit_fraction, top_number = num)
    panel_coverage <- rbind(panel_coverage, new_row)
  }

  # Plot
  panel_coverage <- panel_coverage[order(-coverage, -top_number)]

  panel_coverage[, cancer_type := factor(cancer_type, levels = unique(cancer_type))]
  panel_coverage[, top_number := factor(top_number, levels = unique(top_number))]

  plot <- ggplot(panel_coverage, aes(x = cancer_type, y = coverage, fill = top_number)) +
    geom_bar(stat = "identity", position = "identity", width = 0.9) +
    labs(title = plot_title, x = "Cancer Type", y = "Percentage Coverage", fill = "Number of\nTop Mutations") +
    theme_minimal() +
    theme(
      axis.text.y = element_text(size = 7),
      axis.text.x = element_text(size = 7, angle = 90, vjust = 0.25, hjust = 1),
      panel.grid.major.x = element_blank(),
    ) +
    scale_y_continuous(breaks = seq(0, 100, by = 10), limits = c(0, 100)) +
    scale_fill_brewer(palette = brewer_palette)

  pdf(glue("plots/{Sys.Date()}/{filename}.pdf"))
  print(plot)
  dev.off()

  return(list(plot, panel_coverage))
}


library(data.table)
library(TCGAmutations)

cancer_data <- as.data.frame(tcga_available())

patient_data <- data.table()
top_nums <- c(10, 20, 50, 100)
top_nums_minus_first <- top_nums[-1]

# Get all mutations and patients
for (i in seq_along(cancer_data$Study_Abbreviation)) {
  cancer_type <- cancer_data$Study_Abbreviation[i]
  test_mutations <- tcga_load(study = cancer_type)
  test_data <- as.data.table(test_mutations@data)
  current_data <- test_data[, c("Tumor_Sample_Barcode", "Hugo_Symbol", "HGVSc")]
  current_data <- current_data[, cancer_type := cancer_type]

  patient_data <- rbind(patient_data, current_data)
}

patient_mutation <- patient_data[, .(patient = tstrsplit(Tumor_Sample_Barcode, "-")[[3]], mutation = paste(Hugo_Symbol, HGVSc, sep = ", "), cancer_type)]
patient_mutation <- add_tert(patient_mutation, FALSE)

cancer_specific_coverage <- data.table(cancer_type = character(), coverage = numeric(), top_number = integer(), panel = list())

# Generate panel - cancer-specific
for (i in seq_along(cancer_data$Study_Abbreviation)) {
  current_cancer_type <- cancer_data$Study_Abbreviation[i]
  specific_patients <- patient_mutation[cancer_type == current_cancer_type]
  patient_data_mutation_list <- patient_mutation[, .(mutations = list(mutation)), by = .(patient, cancer_type)]
  unique_patient_data <- patient_data_mutation_list[cancer_type == current_cancer_type]

  # Generate max panel
  panel <- greedy_panel_curator(specific_patients, max(top_nums))

  for (num in top_nums) {
    current_panel <- panel[number <= num]$panel_mutation

    cancer_specific_hits <- unique_patient_data[, .(
      patient,
      mutations,
      top_number = num,
      hits = sapply(mutations, function(x) length(intersect(x, current_panel)))
    )]

    total_num_patients <- cancer_specific_hits[, .N]
    hit_fraction <- cancer_specific_hits[hits > 0, .N] / total_num_patients * 100

    new_row <- data.table(cancer_type = current_cancer_type, coverage = hit_fraction, top_number = num, panel = list(current_panel))
    cancer_specific_coverage <- rbind(cancer_specific_coverage, new_row)
  }
}

write.csv(cancer_specific_coverage, "csv-data/training_cancer_specific_panel.csv", row.names = FALSE)



cancer_specific_coverage <- cancer_specific_coverage[order(-coverage, -top_number)]

# Reorder the cancer_type factor based on the coverage
cancer_specific_coverage[, cancer_type := factor(cancer_type, levels = unique(cancer_type))]
cancer_specific_coverage[, top_number := factor(top_number, levels = unique(top_number))]

cancer_specific_panel_plot <- ggplot(cancer_specific_coverage, aes(x = cancer_type, y = coverage, fill = top_number)) +
  geom_bar(stat = "identity", position = "identity", width = 0.9) +
  labs(title = "Cancer-specific Panel Coverage (Greedy)", x = "Cancer Type", y = "Percentage Coverage", fill = "Number of\nTop Mutations") +
  theme_minimal() +
  theme(
    axis.text.y = element_text(size = 7),
    axis.text.x = element_text(size = 7, angle = 90, vjust = 0.5, hjust = 1)
  ) +
  scale_y_continuous(breaks = seq(0, 100, by = 10), limits = c(0, 100))


pdf("plots/cancer-specific-panel-coverage-greedy-new.pdf")
cancer_specific_panel_plot
dev.off()


top_nums <- c(10, 20, 50, 100)
pan_can_coverage <- data.table(cancer_type = character(), coverage = numeric(), top_number = integer())
pan_can_panel <- data.table(top_number = integer(), panel = list())
patient_data <- data.table()

for (i in seq_along(cancer_data$Study_Abbreviation)) {
  cancer_type <- cancer_data$Study_Abbreviation[i]
  test_mutations <- tcga_load(study = cancer_type)
  test_data <- as.data.table(test_mutations@data)
  current_data <- test_data[, c("Tumor_Sample_Barcode", "Hugo_Symbol", "HGVSc")]
  current_data <- current_data[, cancer_type := cancer_type]

  patient_data <- rbind(patient_data, current_data)
}

patient_mutation <- patient_data[, .(patient = tstrsplit(Tumor_Sample_Barcode, "-")[[3]], mutation = paste(Hugo_Symbol, HGVSc, sep = ", "), cancer_type)]
patient_mutation <- add_tert(patient_mutation, FALSE)

# Generate pan can panel
specific_patients <- copy(patient_mutation)
patient_data_mutation_list <- patient_mutation[, .(mutations = list(mutation)), by = .(patient, cancer_type)]
unique_patient_data <- copy(patient_data_mutation_list)

panel <- greedy_panel_curator(specific_patients, max(top_nums))

for (num in top_nums) {
  current_panel <- panel[number <= num]$panel_mutation

  pan_can_hits <- unique_patient_data[, .(
    patient,
    mutations,
    top_number = num,
    hits = sapply(mutations, function(x) length(intersect(x, current_panel)))
  )]

  total_num_patients <- pan_can_hits[, .N]
  hit_fraction <- pan_can_hits[hits > 0, .N] / total_num_patients * 100

  new_row <- data.table(cancer_type = "PAN-CAN", coverage = hit_fraction, top_number = num)
  pan_can_coverage <- rbind(pan_can_coverage, new_row)

  new_panel <- data.table(top_number = num, panel = list(current_panel))
  pan_can_panel <- rbind(pan_can_panel, new_panel)
}

fwrite(pan_can_panel, "csv-data/pan-cancer-panel.csv")


# Get fractions for all cancer types
for (i in seq_along(cancer_data$Study_Abbreviation)) {
  current_cancer_type <- cancer_data$Study_Abbreviation[i]
  specific_patients <- cleaned_patient_data[cancer_type == current_cancer_type]
  patient_data_mutation_list <- patient_mutation[, .(mutations = list(mutation)), by = .(patient, cancer_type)]
  unique_patient_data <- patient_data_mutation_list[cancer_type == current_cancer_type]

  # Generate panel
  for (num in top_nums) {
    panel <- unlist(pan_can_panel[top_number == num]$panel)

    cancer_specific_hits <- unique_patient_data[, .(
      patient,
      mutations,
      top_number = num,
      hits = sapply(mutations, function(x) length(intersect(x, panel)))
    )]

    total_num_patients <- cancer_specific_hits[, .N]
    hit_fraction <- cancer_specific_hits[hits > 0, .N] / total_num_patients * 100

    new_row <- data.table(cancer_type = current_cancer_type, coverage = hit_fraction, top_number = num)
    pan_can_coverage <- rbind(pan_can_coverage, new_row)
  }
}

fwrite(pan_can_coverage, "csv-data/pan_cancer_panel_coverage.csv")


pan_can_coverage <- pan_can_coverage[order(-coverage, -top_number)]

# Reorder the cancer_type factor based on the coverage
pan_can_coverage[, cancer_type := factor(cancer_type, levels = unique(cancer_type))]
pan_can_coverage[, top_number := factor(top_number, levels = unique(top_number))]

num_levels <- length(unique(pan_can_coverage$top_number))
reversed_viridis_colors <- rev(viridis(num_levels))

pan_can_panel_plot <- ggplot(pan_can_coverage, aes(x = cancer_type, y = coverage, fill = top_number)) +
  geom_bar(stat = "identity", position = "identity", width = 0.9) +
  labs(title = "Pan-Can Panel Coverage (Greedy)", x = "Cancer Type", y = "Percentage Coverage", fill = "Number of\nTop Mutations") +
  theme_minimal() +
  theme(
    axis.text.y = element_text(size = 7),
    axis.text.x = element_text(size = 7, angle = 90, vjust = 0.5, hjust = 1)
  ) +
  scale_y_continuous(breaks = seq(0, 100, by = 10), limits = c(0, 100)) +
  scale_fill_manual(values = reversed_viridis_colors)


pdf("plots/pan-can-panel-coverage-greedy-new.pdf")
pan_can_panel_plot
dev.off()


cv_fraction <- 0.6

# Split the dataset
cancer_data <- as.data.frame(tcga_available())

patient_data <- data.table()
patient_data_training <- data.table()
patient_data_test <- data.table()
top_nums <- c(10, 20, 50, 100)
top_nums_minus_first <- top_nums[-1]
set.seed(100)

# Get all mutations and patients
# for (i in 1) {
for (i in seq_along(cancer_data$Study_Abbreviation)) {
  cancer_type <- cancer_data$Study_Abbreviation[i]
  # cancer_type <- "BRCA"
  test_mutations <- tcga_load(study = cancer_type)
  test_data <- as.data.table(test_mutations@data)
  current_data <- test_data[, c("Tumor_Sample_Barcode", "Hugo_Symbol", "HGVSc")]
  current_data[, cancer_type := cancer_type]
  current_data <- current_data[, .(patient = tstrsplit(Tumor_Sample_Barcode, "-")[[3]], mutation = paste(Hugo_Symbol, HGVSc, sep = ", "), cancer_type)]

  patient_data <- rbind(patient_data, current_data)

  unique_patients <- current_data[, unique(patient)]
  training_patients <- sample(unique_patients, round(length(unique_patients) * 0.7))

  training_data <- current_data[patient %in% training_patients]
  test_data <- current_data[!patient %in% training_patients]

  patient_data_training <- rbind(patient_data_training, training_data)
  patient_data_test <- rbind(patient_data_test, test_data)
}

patient_mutation <- copy(patient_data)
patient_mutation <- add_tert(patient_mutation, FALSE)
patient_data_training <- unique(add_tert(patient_data_training, FALSE))
patient_data_test <- unique(add_tert(patient_data_test, FALSE))

# debug_dt <- data.table()
# debug_dt <- cbind_fill(debug_dt, patient_data_training, patient_data_test, fill = NA)

patient_data_mutation_list <- patient_data_training[, .(mutations = list(mutation)), by = .(patient, cancer_type)]

# Generate panel - cancer-specific
training_cancer_specific <- data.table()
panel_coverage <- data.table()
cv_coverage <- data.table()
# for (i in 1) {
for (i in seq_along(cancer_data$Study_Abbreviation)) {
  print(cancer_data$Study_Abbreviation[i])
  for (run_num in 1:10) {
    # CHANGE IF YOU WANT TRULY 'RANDOM', FOR REPRODUCIBILITY
    # set.seed(NULL)
    set.seed(run_num)

    current_cancer_type <- cancer_data$Study_Abbreviation[i]
    # current_cancer_type <- "BRCA"
    specific_patients <- patient_data_training[cancer_type == current_cancer_type]

    # sample 80% of training data
    unique_specific_patients <- specific_patients[, unique(patient)]
    eighty_training_patients <- sample(unique_specific_patients, round(length(unique_specific_patients) * cv_fraction))

    eighty_training_data <- specific_patients[patient %in% eighty_training_patients]

    # Generate max panel
    panel <- greedy_panel_curator(eighty_training_data, max(top_nums))

    panel[, cancer_type := current_cancer_type]
    panel_coverage <- rbind(panel_coverage, panel)

    new_row <- data.table(cancer_type = current_cancer_type, run = run_num, panel = list(panel$panel_mutation))

    training_cancer_specific <- rbind(training_cancer_specific, new_row)


    specific_pm_list <- patient_data_mutation_list[cancer_type == current_cancer_type]
    test_patients <- specific_pm_list[!patient %in% eighty_training_patients]

    cancer_specific_hits <- test_patients[, .(
      patient,
      mutations,
      top_number = num,
      hits = sapply(mutations, function(x) length(intersect(x, panel$panel_mutation)))
    )]

    total_num_patients <- cancer_specific_hits[, .N]
    hit_fraction <- cancer_specific_hits[hits > 0, .N] / total_num_patients * 100

    new_row <- data.table(cancer_type = current_cancer_type, coverage = hit_fraction)
    cv_coverage <- rbind(cv_coverage, new_row)
  }
}

dt_avg <- cv_coverage[, .(average_coverage = mean(coverage)), by = cancer_type]
dt_avg <- dt_avg[order(-average_coverage)]

# Reorder the cancer_type factor based on the coverage
dt_avg[, cancer_type := factor(cancer_type, levels = unique(cancer_type))]
# dt_avg[, top_number := factor(top_number, levels = unique(top_number))]

cancer_specific_panel_plot <- ggplot(dt_avg, aes(x = cancer_type, y = average_coverage)) +
  geom_bar(stat = "identity", position = "identity", width = 0.9, fill = "blue") +
  labs(title = "Cancer-specific Panel Coverage CV Average", x = "Cancer Type", y = "Percentage Coverage") +
  theme_minimal() +
  theme(
    axis.text.y = element_text(size = 7),
    axis.text.x = element_text(size = 7, angle = 90, vjust = 0.5, hjust = 1)
  ) +
  scale_y_continuous(breaks = seq(0, 100, by = 10), limits = c(0, 100))


pdf("plots/monte-carlo-cv-results.pdf")
cancer_specific_panel_plot
dev.off()


# debug_dt <- data.table()
# debug_dt <- cbind_fill(eighty_training_data, patient_data_mutation_list, unique_patient_data, training_cancer_specific, fill = NA)
# panels_minus_one <- panel_coverage[num_covered > 1, .(n = .N), by = .(panel_mutation, cancer_type, number)]

training_data_save <- copy(training_cancer_specific)
training_data_save$panel <- sapply(training_cancer_specific$panel, function(x) paste(x, collapse = ","))
write.csv(training_data_save, "csv-data/training_cancer_specific_panel.csv", row.names = FALSE)


all_training_dt <- copy(training_cancer_specific)

# Function to calculate unique mutation count for a given num_mutation
calculate_unique_mutations <- function(data, num_mutation) {
  only_n_mutations <- data[, .(mutation = unlist(lapply(panel, function(x) head(x, num_mutation)))), by = .(cancer_type, run)]
  unique_mutation_count <- only_n_mutations[, uniqueN(mutation)]
  return(unique_mutation_count)
}

final_panel <- data.table()

for (cancer in unique(all_training_dt$cancer_type)) {
  current_data <- all_training_dt[cancer_type == cancer]
  top_mutations <- current_data[, .(mutation = unlist(panel)), by = .(cancer_type, run)]
  unique_tops <- unique(top_mutations$mutation)
  if (length(unique_tops) < max(top_nums)) {
    # Unlist all mutations, merge into one list then make unique
    new_row <- data.table(cancer_type = cancer, panel = list(unique_tops))
    final_panel <- rbind(final_panel, new_row, fill = TRUE)
  } else {
    left <- 0
    right <- max(top_nums)
    target <- max(top_nums)
    current_num <- 0
    mid <- NA

    while (left <= right) {
      mid <- left + (right - left) %/% 2
      current_num <- calculate_unique_mutations(current_data, mid)

      if (abs(current_num - target) <= 5) {
        break
      }

      if (current_num < target) {
        left <- mid + 1
      } else {
        right <- mid - 1
      }
    }


    only_n_mutations <- current_data[, panel := lapply(panel, function(x) head(x, mid))]
    tops <- only_n_mutations[, .(mutation = unlist(panel)), by = .(cancer_type, run)]
    unique_tops <- unique(tops$mutation)

    new_row <- data.table(cancer_type = cancer, panel = list(unique_tops))
    final_panel <- rbind(final_panel, new_row, fill = TRUE)
  }
}

all_unlisted <- all_training_dt[, .(mutation = unlist(panel)), by = setdiff(names(all_training_dt), "panel")]
drop <- all_unlisted[, .(cancer_type, mutation)]
drop_count <- drop[, .(count = .N), by = .(mutation, cancer_type)][order(-count)]


final_panel_raw <- drop_count[, head(.SD, 100), by = cancer_type]
final_panel_raw <- final_panel_raw[, .(panel = list(mutation)), by = cancer_type]


# THIS IS FOR UNOPTIMIZED FINAL PANEL (RAW)
# plot on DISCOVERY SET
trained_cs_coverage_dis <- data.table()
for (i in seq_along(cancer_data$Study_Abbreviation)) {
  current_cancer_type <- cancer_data$Study_Abbreviation[i]
  specific_patients <- patient_data_training[cancer_type == current_cancer_type]
  patient_data_mutation_list <- patient_data_training[, .(mutations = list(mutation)), by = .(patient, cancer_type)]
  unique_patient_data <- patient_data_mutation_list[cancer_type == current_cancer_type]

  panel <- unlist(final_panel_raw[cancer_type == current_cancer_type, panel])

  for (num in top_nums) {
    current_panel <- panel[1:num]

    cancer_specific_hits <- unique_patient_data[, .(
      patient,
      mutations,
      top_number = num,
      hits = sapply(mutations, function(x) length(intersect(x, current_panel)))
    )]

    total_num_patients <- cancer_specific_hits[, .N]
    hit_fraction <- cancer_specific_hits[hits > 0, .N] / total_num_patients * 100

    new_row <- data.table(cancer_type = current_cancer_type, coverage = hit_fraction, top_number = num)
    trained_cs_coverage_dis <- rbind(trained_cs_coverage_dis, new_row)
  }
}

trained_cs_coverage_dis <- trained_cs_coverage_dis[order(-coverage, -top_number)]

# Reorder the cancer_type factor based on the coverage
trained_cs_coverage_dis[, cancer_type := factor(cancer_type, levels = unique(cancer_type))]
trained_cs_coverage_dis[, top_number := factor(top_number, levels = unique(top_number))]

unoptimized_cs_plot_dis <- ggplot(trained_cs_coverage_dis, aes(x = cancer_type, y = coverage, fill = top_number)) +
  geom_bar(stat = "identity", position = "identity", width = 0.9) +
  labs(title = "Cancer-specific Panel Coverage (Greedy, CV, Discovery)", x = "Cancer Type", y = "Percentage Coverage", fill = "Number of\nTop Mutations") +
  theme_minimal() +
  theme(
    axis.text.y = element_text(size = 7),
    axis.text.x = element_text(size = 7, angle = 90, vjust = 0.5, hjust = 1)
  ) +
  scale_y_continuous(breaks = seq(0, 100, by = 10), limits = c(0, 100)) +
  scale_fill_brewer(palette = "Set1")

pdf("plots/cancer-specific-panel-coverage-greedy-trained-discovery-unoptimized.pdf")
unoptimized_cs_plot_dis
dev.off()



# Get fractions for all cancer types when tested on TEST SET
trained_cs_coverage <- data.table()
for (i in seq_along(cancer_data$Study_Abbreviation)) {
  current_cancer_type <- cancer_data$Study_Abbreviation[i]
  # current_cancer_type <- "BRCA"
  specific_patients <- patient_data_test[cancer_type == current_cancer_type]
  patient_data_mutation_list <- patient_data_test[, .(mutations = list(mutation)), by = .(patient, cancer_type)]
  unique_patient_data <- patient_data_mutation_list[cancer_type == current_cancer_type]

  panel <- unlist(final_panel_raw[cancer_type == current_cancer_type, panel])

  for (num in top_nums) {
    current_panel <- panel[1:num]

    cancer_specific_hits <- unique_patient_data[, .(
      patient,
      mutations,
      top_number = num,
      hits = sapply(mutations, function(x) length(intersect(x, current_panel)))
    )]

    total_num_patients <- cancer_specific_hits[, .N]
    hit_fraction <- cancer_specific_hits[hits > 0, .N] / total_num_patients * 100

    new_row <- data.table(cancer_type = current_cancer_type, coverage = hit_fraction, top_number = num)
    trained_cs_coverage <- rbind(trained_cs_coverage, new_row)
  }
}

trained_cs_coverage <- trained_cs_coverage[order(-coverage, -top_number)]

# Reorder the cancer_type factor based on the coverage
trained_cs_coverage[, cancer_type := factor(cancer_type, levels = unique(cancer_type))]
trained_cs_coverage[, top_number := factor(top_number, levels = unique(top_number))]

unoptimized_cs_plot_test <- ggplot(trained_cs_coverage, aes(x = cancer_type, y = coverage, fill = top_number)) +
  geom_bar(stat = "identity", position = "identity", width = 0.9) +
  labs(title = "Cancer-specific Panel Coverage (Greedy, CV, Test)", x = "Cancer Type", y = "Percentage Coverage", fill = "Number of\nTop Mutations") +
  theme_minimal() +
  theme(
    axis.text.y = element_text(size = 5),
    axis.text.x = element_text(size = 5, angle = 90, vjust = 0.5, hjust = 1),
  ) +
  scale_y_continuous(breaks = seq(0, 100, by = 10), limits = c(0, 100)) +
  scale_fill_brewer(palette = "Set1")


pdf("plots/cancer-specific-panel-coverage-greedy-trained-unoptimized.pdf")
unoptimized_cs_plot_test
dev.off()

# THIS IS FOR OPTIMIZED PANEL
# plot on DISCOVERY SET
trained_cs_coverage <- data.table()
for (i in seq_along(cancer_data$Study_Abbreviation)) {
  current_cancer_type <- cancer_data$Study_Abbreviation[i]
  specific_patients <- patient_data_training[cancer_type == current_cancer_type]
  patient_data_mutation_list <- patient_data_training[, .(mutations = list(mutation)), by = .(patient, cancer_type)]
  unique_patient_data <- patient_data_mutation_list[cancer_type == current_cancer_type]

  panel <- unlist(final_panel[cancer_type == current_cancer_type, panel])

  for (num in top_nums) {
    current_panel <- panel[1:num]

    cancer_specific_hits <- unique_patient_data[, .(
      patient,
      mutations,
      top_number = num,
      hits = sapply(mutations, function(x) length(intersect(x, current_panel)))
    )]

    total_num_patients <- cancer_specific_hits[, .N]
    hit_fraction <- cancer_specific_hits[hits > 0, .N] / total_num_patients * 100

    new_row <- data.table(cancer_type = current_cancer_type, coverage = hit_fraction, top_number = num)
    trained_cs_coverage <- rbind(trained_cs_coverage, new_row)
  }
}

trained_cs_coverage <- trained_cs_coverage[order(-coverage, -top_number)]

# Reorder the cancer_type factor based on the coverage
trained_cs_coverage[, cancer_type := factor(cancer_type, levels = unique(cancer_type))]
trained_cs_coverage[, top_number := factor(top_number, levels = unique(top_number))]

cs_panel_dis_plot <- ggplot(trained_cs_coverage, aes(x = cancer_type, y = coverage, fill = top_number)) +
  geom_bar(stat = "identity", position = "identity", width = 0.9) +
  labs(title = "Cancer-specific Panel Coverage (Greedy, CV, Discovery, Optimized)", x = "Cancer Type", y = "Percentage Coverage", fill = "Number of\nTop Mutations") +
  theme_minimal() +
  theme(
    axis.text.y = element_text(size = 7),
    axis.text.x = element_text(size = 7, angle = 90, vjust = 0.5, hjust = 1)
  ) +
  scale_y_continuous(breaks = seq(0, 100, by = 10), limits = c(0, 100)) +
  scale_fill_brewer(palette = "Set2")




pdf("plots/cancer-specific-panel-coverage-greedy-trained-discovery-optimized.pdf")
cs_panel_dis_plot
dev.off()



# Get fractions for all cancer types when tested on TEST SET
trained_cs_coverage <- data.table()
for (i in seq_along(cancer_data$Study_Abbreviation)) {
  current_cancer_type <- cancer_data$Study_Abbreviation[i]
  # current_cancer_type <- "BRCA"
  specific_patients <- patient_data_test[cancer_type == current_cancer_type]
  patient_data_mutation_list <- patient_data_test[, .(mutations = list(mutation)), by = .(patient, cancer_type)]
  unique_patient_data <- patient_data_mutation_list[cancer_type == current_cancer_type]

  panel <- unlist(final_panel[cancer_type == current_cancer_type, panel])

  for (num in top_nums) {
    current_panel <- panel[1:num]

    cancer_specific_hits <- unique_patient_data[, .(
      patient,
      mutations,
      top_number = num,
      hits = sapply(mutations, function(x) length(intersect(x, current_panel)))
    )]

    total_num_patients <- cancer_specific_hits[, .N]
    hit_fraction <- cancer_specific_hits[hits > 0, .N] / total_num_patients * 100

    new_row <- data.table(cancer_type = current_cancer_type, coverage = hit_fraction, top_number = num)
    trained_cs_coverage <- rbind(trained_cs_coverage, new_row)
  }
}

trained_cs_coverage <- trained_cs_coverage[order(-coverage, -top_number)]

# Reorder the cancer_type factor based on the coverage
trained_cs_coverage[, cancer_type := factor(cancer_type, levels = unique(cancer_type))]
trained_cs_coverage[, top_number := factor(top_number, levels = unique(top_number))]

cs_panel_test_plot <- ggplot(trained_cs_coverage, aes(x = cancer_type, y = coverage, fill = top_number)) +
  geom_bar(stat = "identity", position = "identity", width = 0.9) +
  labs(title = "Cancer-specific Panel Coverage (Greedy, CV, Test, Optimized)", x = "Cancer Type", y = "Percentage Coverage", fill = "Number of\nTop Mutations") +
  theme_minimal() +
  theme(
    axis.text.y = element_text(size = 7),
    axis.text.x = element_text(size = 7, angle = 90, vjust = 0.5, hjust = 1)
  ) +
  scale_y_continuous(breaks = seq(0, 100, by = 10), limits = c(0, 100)) +
  scale_fill_brewer(palette = "Set2")


pdf("plots/cancer-specific-panel-coverage-greedy-trained-optimized.pdf")
cs_panel_test_plot
dev.off()


cs_legend <- get_legend(unoptimized_cs_plot_test)

# Remove legends from the plots
unoptimized_cs_plot_dis <- unoptimized_cs_plot_dis + theme(
  legend.position = "none",
  axis.title = element_blank(),
  title = element_blank(),
  axis.text.y = element_text(size = 5),
  axis.text.x = element_text(size = 5, angle = 90, vjust = 0.5, hjust = 1)
)
unoptimized_cs_plot_test <- unoptimized_cs_plot_test + theme(
  legend.position = "none",
  axis.title = element_blank(),
  title = element_blank(),
  axis.text.y = element_text(size = 5),
  axis.text.x = element_text(size = 5, angle = 90, vjust = 0.5, hjust = 1)
)
cs_panel_dis_plot <- cs_panel_dis_plot + theme(
  legend.position = "none",
  axis.title = element_blank(),
  title = element_blank(),
  axis.text.y = element_text(size = 5),
  axis.text.x = element_text(size = 5, angle = 90, vjust = 0.5, hjust = 1)
)
cs_panel_test_plot <- cs_panel_test_plot + theme(
  legend.position = "none",
  axis.title = element_blank(),
  title = element_blank(),
  axis.text.y = element_text(size = 5),
  axis.text.x = element_text(size = 5, angle = 90, vjust = 0.5, hjust = 1)
)

# Combine using cowplot
unopt_plots <- plot_grid(unoptimized_cs_plot_dis, cs_panel_dis_plot, unoptimized_cs_plot_test, cs_panel_test_plot, align = "vh", vjust = 1, scale = 1)

# Create common x and y labels
y.grob <- textGrob("Percentage Coverage", gp = gpar(fontsize = 10), rot = 90)
x.grob <- textGrob("Cancer Type", gp = gpar(fontsize = 10))

# Arrange the combined plot and labels
final_unopt_plots <- arrangeGrob(unopt_plots, left = y.grob, bottom = x.grob)

# Save to PDF
pdf("plots/unopt-vs-opt-cs-panel.pdf")
grid.draw(final_unopt_plots)
dev.off()


# NEED TO RUN THE CHUNK BELOW FIRST (pan-can one)

cv_fraction <- 0.6
TERT <- FALSE

cancer_data <- as.data.frame(tcga_available())

patient_data <- data.table()
patient_data_training <- data.table()
patient_data_test <- data.table()
colormap_genes <- data.table()
colormap_mutations <- data.table()

top_nums <- c(10, 20, 50, 100)
top_nums_minus_first <- top_nums[-1]
set.seed(100)

# Get all mutations and patients
for (i in seq_along(cancer_data$Study_Abbreviation)) {
  cancer_type <- cancer_data$Study_Abbreviation[i]
  test_mutations <- tcga_load(study = cancer_type)
  test_data <- as.data.table(test_mutations@data)
  raw_data <- test_data[, c("Tumor_Sample_Barcode", "Hugo_Symbol", "HGVSc")]
  raw_data[, cancer_type := cancer_type]

  new_genes <- raw_data[, .(gene = Hugo_Symbol)]
  colormap_genes <- rbind(colormap_genes, new_genes)

  current_data <- raw_data[, .(patient = tstrsplit(Tumor_Sample_Barcode, "-")[[3]], mutation = paste(Hugo_Symbol, HGVSc, sep = ", "), cancer_type)]

  new_mutations <- current_data[, .(mutation = mutation)]
  colormap_mutations <- rbind(colormap_mutations, new_mutations)

  patient_data <- rbind(patient_data, current_data)

  unique_patients <- current_data[, unique(patient)]
  training_patients <- sample(unique_patients, round(length(unique_patients) * 0.7))

  training_data <- current_data[patient %in% training_patients]
  test_data <- current_data[!patient %in% training_patients]

  patient_data_training <- rbind(patient_data_training, training_data)
  patient_data_test <- rbind(patient_data_test, test_data)
}

patient_mutation <- copy(patient_data)
patient_data_training <- unique(add_tert(patient_data_training, TERT))
patient_data_test <- unique(add_tert(patient_data_test, TERT))


# Generate pan can panel
specific_patients <- copy(patient_data_training)

# sample fraction of training data
unique_specific_patients <- specific_patients[, unique(patient)]
patient_data_mutation_list <- patient_data_training[, .(mutations = list(mutation)), by = .(patient, cancer_type)]

# Generate pan-can panel
top_cv_panel <- data.table()
for (run_num in 1:10) {
  set.seed(run_num)
  eighty_training_patients <- sample(unique_specific_patients, round(length(unique_specific_patients) * cv_fraction))

  eighty_training_data <- specific_patients[patient %in% eighty_training_patients]

  unique_patient_data <- patient_data_mutation_list[patient %in% eighty_training_patients]

  tops_panel <- eighty_training_data[, .(n = .N), by = mutation][order(-n)][1:max(top_nums)]

  new_row <- data.table(cancer_type = "PAN-CAN", run = run_num, panel = list(tops_panel$mutation))
  top_cv_panel <- rbind(top_cv_panel, new_row)
}


# Generate panel - cancer-specific
for (i in seq_along(cancer_data$Study_Abbreviation)) {
  print(cancer_data$Study_Abbreviation[i])
  for (run_num in 1:10) {
    set.seed(run_num)

    current_cancer_type <- cancer_data$Study_Abbreviation[i]
    specific_patients <- patient_data_training[cancer_type == current_cancer_type]

    # sample 80% of training data
    unique_specific_patients <- unique(specific_patients$patient)

    eighty_training_patients <- sample(unique_specific_patients, round(length(unique_specific_patients) * cv_fraction))

    eighty_training_data <- specific_patients[patient %in% eighty_training_patients]

    unique_patient_data <- patient_data_mutation_list[cancer_type == current_cancer_type]
    unique_patient_data <- unique_patient_data[patient %in% eighty_training_patients]

    # Generate top 100 mutations
    panel <- eighty_training_data[, .(n = .N), by = mutation][order(-n)][1:max(top_nums)]

    # panel <- greedy_panel_curator(eighty_training_data, max(top_nums))

    new_row <- data.table(cancer_type = current_cancer_type, run = run_num, panel = list(panel$mutation))

    top_cv_panel <- rbind(top_cv_panel, new_row)
  }
}

all_training_dt <- copy(top_cv_panel)

all_unlisted <- all_training_dt[, .(mutation = unlist(panel)), by = setdiff(names(all_training_dt), "panel")]
drop <- all_unlisted[, .(cancer_type, mutation)]
drop_count <- drop[, .(count = .N), by = .(mutation, cancer_type)][order(-count)]

final_panel <- drop_count[, .SD[1:min(100, .N)], by = cancer_type]
final_panel <- final_panel[, .(panel = list(mutation)), by = cancer_type]

cs_dis <- test_and_plot_cs(final_panel, patient_data_training, "Cancer Specific Panel Coverage (On Discovery)", "cs-tops-pluspc-ondis", "Set2", TRUE)
cs_test <- test_and_plot_cs(final_panel, patient_data_test, "Cancer Specific Panel Coverage", "cs-tops-pluspc-ontest", "Set2", TRUE)
cs_dis_plot <- cs_dis[[1]]
cs_test_plot <- cs_test[[1]]
pan_can_tops_dis_plot <- test_and_plot_pc(final_panel_tops, patient_data_training, "Pan-Cancer Panel (Tops, On Discovery)", "pan-can-tops-cved-ondis", "Set2")
pan_can_tops_test_plot <- test_and_plot_pc(final_panel_tops, patient_data_test, "Pan-Cancer Panel (Tops, On Test)", "pan-can-tops-cved-ontest", "Set2")




prepare_for_cow <- function(plot) {
  plot + theme(
    legend.position = "none",
    axis.title = element_blank(),
    plot.title = element_text(size = 10),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.y = element_blank(),
    axis.text.x = element_text(size = 8), 
  ) 
}

cs_dis_plot_copy <- copy(cs_dis_plot)

plots <- list(cs_dis_plot, pan_can_tops_dis_plot, cs_test_plot, pan_can_tops_test_plot)
prepared_plots <- lapply(plots, prepare_for_cow)

combined_plots <- plot_grid(
  plotlist = prepared_plots,
  align = "vh", axis = "tb", ncol = 2
)

# Create common x and y labels
y.grob <- textGrob("Percentage Coverage", gp = gpar(fontsize = 15), rot = 90)
x.grob <- textGrob("Cancer Type", gp = gpar(fontsize = 15))

# Arrange the combined plot and labels
final_plot_pc <- arrangeGrob(combined_plots, left = y.grob, bottom = x.grob)

pdf(glue("plots/{Sys.Date()}/cs-vs-pc-cow.pdf"), width = 17, height = 8.5)
grid.draw(final_plot_pc)
dev.off()


# # Prepare for colormap
# Process genes/mutations
unique_genes <- colormap_genes[, .(gene = unique(gene))]
unique_mutations <- colormap_mutations[, .(mutation = unique(mutation))]
uni_mut <- unique_mutations[, count := 0]
uni_genes <- unique_genes[, count := 0]

final_panel <- final_panel[cancer_type == "PAN-CAN", cancer_type := "PANCAN"]

for (i in seq_len(nrow(final_panel))) {
  current_panel <- final_panel$panel[[i]][1:20]
  uni_mut[mutation %in% current_panel, count := count + 1]
}

final_panel_genes <- final_panel[, genes := lapply(panel, function(x) unique(sapply(x, function(y) strsplit(y, ",")[[1]][1])))]

for (i in seq_len(nrow(final_panel_genes))) {
  current_panel <- final_panel_genes$genes[[i]][1:20]
  uni_genes[gene %in% current_panel, count := count + 1]
}

top_ten_genes <- uni_genes[order(-count)][1:20]
top_ten_mut <- uni_mut[order(-count)][1:20]

# Check in which top n panel do the genes/muts cover

mutation_panel_hits <- data.table()
mutation_panel_hits_complete <- data.table()
for (i in seq_len(nrow(final_panel))) {
  current_panel <- final_panel$panel[[i]]
  current_cancer_type <- final_panel$cancer_type[i]
  hit <- top_ten_mut[mutation %in% current_panel]

  new_row <- data.table(mutation = hit$mutation, cancer_type = current_cancer_type, top_number = as.factor(sapply(match(hit$mutation, current_panel), round_up_custom)))
  mutation_panel_hits <- rbind(mutation_panel_hits, new_row)
}


# FIND FREQUENCY OF MUTATIONS IN EACH CANCER TYPE + PANCAN
mutation_panel_hits <- mutation_panel_hits[, freq := 0]

for (i in seq_along(cancer_data$Study_Abbreviation)) {
  current_cancer_type <- cancer_data$Study_Abbreviation[i]
  current_mutations <- mutation_panel_hits[cancer_type == current_cancer_type]

  raw_data <- patient_data_test[cancer_type == current_cancer_type]
  total_patients <- length(raw_data[, unique(patient)])

  mutated_patients <- raw_data[mutation %in% current_mutations$mutation]

  mutation_fractions <- mutated_patients[, .(covered_patients = .N), by = mutation]
  mutation_fractions <- mutation_fractions[, cancer_type := current_cancer_type]
  mutation_fractions <- mutation_fractions[, freq := covered_patients / total_patients]

  mutation_panel_hits <- mutation_panel_hits[mutation_fractions, on = .(mutation, cancer_type), freq := i.freq]
}

# For pan cancer
current_cancer_type <- "PANCAN"
current_mutations <- mutation_panel_hits[cancer_type == current_cancer_type]

raw_data <- copy(patient_data_test)

total_patients <- length(raw_data[, unique(patient)])

mutated_patients <- raw_data[mutation %in% current_mutations$mutation]

mutation_fractions <- mutated_patients[, .(covered_patients = .N), by = mutation]
mutation_fractions <- mutation_fractions[, cancer_type := current_cancer_type]
mutation_fractions <- mutation_fractions[, freq := covered_patients / total_patients]

mutation_panel_hits <- mutation_panel_hits[mutation_fractions, on = .(mutation, cancer_type), freq := i.freq]

missing_ct_mut <- setdiff(final_panel$cancer_type, mutation_panel_hits$cancer_type)
placeholder_rows_mut <- data.table(mutation = mutation_panel_hits$mutation[1], cancer_type = missing_ct_mut, top_number = NA)
mutation_panel_hits_complete <- rbind(mutation_panel_hits, placeholder_rows_mut, fill = TRUE)

# This final panel is not unique by design to mimic the top 10 panel
final_panel_genes_nu <- final_panel[, genes := lapply(panel, function(x) sapply(x, function(y) strsplit(y, ",")[[1]][1]))]

gene_panel_hits <- data.table()
for (i in seq_len(nrow(final_panel_genes_nu))) {
  current_panel <- final_panel_genes_nu$genes[[i]]
  current_cancer_type <- final_panel_genes_nu$cancer_type[i]
  hit <- top_ten_genes[gene %in% current_panel]
  new_row <- data.table(gene = hit$gene, cancer_type = current_cancer_type, top_number = as.factor(sapply(match(hit$gene, current_panel), round_up_custom)))
  gene_panel_hits <- rbind(gene_panel_hits, new_row)
}

gene_panel_hits <- gene_panel_hits[, freq := 0]

patient_data_test_gene <- unique(patient_data_test[, .(gene = tstrsplit(mutation, ",")[[1]]), by = .(patient, cancer_type)])

for (i in seq_along(cancer_data$Study_Abbreviation)) {
  current_cancer_type <- cancer_data$Study_Abbreviation[i]
  current_genes <- gene_panel_hits[cancer_type == current_cancer_type]

  raw_data <- patient_data_test_gene[cancer_type == current_cancer_type]
  total_patients <- length(raw_data[, unique(patient)])

  mutated_patients <- raw_data[gene %in% current_genes$gene]

  gene_fractions <- mutated_patients[, .(covered_patients = .N), by = gene]
  gene_fractions <- gene_fractions[, cancer_type := current_cancer_type]
  gene_fractions <- gene_fractions[, freq := covered_patients / total_patients]

  gene_panel_hits <- gene_panel_hits[gene_fractions, on = .(gene, cancer_type), freq := i.freq]
}


current_cancer_type <- "PANCAN"
current_genes <- gene_panel_hits[cancer_type == current_cancer_type]

raw_data <- copy(patient_data_test_gene)
total_patients <- length(raw_data[, unique(patient)])

mutated_patients <- raw_data[gene %in% current_genes$gene]

gene_fractions <- mutated_patients[, .(covered_patients = .N), by = gene]
gene_fractions <- gene_fractions[, cancer_type := current_cancer_type]
gene_fractions <- gene_fractions[, freq := covered_patients / total_patients]

gene_panel_hits <- gene_panel_hits[gene_fractions, on = .(gene, cancer_type), freq := i.freq]


missing_ct_gene <- setdiff(final_panel$cancer_type, gene_panel_hits$cancer_type)
placeholder_rows_gene <- data.table(gene = gene_panel_hits$gene[1], cancer_type = missing_ct_gene, top_number = NA)
gene_panel_hits_complete <- rbind(gene_panel_hits, placeholder_rows_gene, fill = TRUE)



cs_test_panel_cov <- cs_test[[2]]$cancer_type
levels_cs_test <- levels(cs_test_panel_cov)

# Gene-cancer heatmap
gene_panel_hits_complete[freq > 0, count := .N, by = gene]

gene_panel_hits_complete[, cancer_type := factor(cancer_type, levels = levels_cs_test)]
gene_panel_hits_complete[, gene := factor(gene, levels = unique(gene[order(count)]))]

# Top 10 only
gene_panel_hits_complete <- gene_panel_hits_complete[top_number == 10 | is.na(top_number)]

missing_ct_gene <- setdiff(final_panel$cancer_type, gene_panel_hits_complete$cancer_type)
placeholder_rows_gene <- data.table(gene = gene_panel_hits_complete$gene[1], cancer_type = missing_ct_gene, top_number = NA)
gene_panel_hits_complete <- rbind(gene_panel_hits_complete, placeholder_rows_gene, fill = TRUE)

gene_colormap_plot <- ggplot(gene_panel_hits_complete, aes(cancer_type, gene, fill = freq)) +

  # # All tops
  # gene_colormap_plot <- ggplot(gene_panel_hits_complete, aes(cancer_type, gene, fill = top_number, alpha = freq)) +
  geom_tile(color = "white") +
  theme_minimal() +
  # scale_fill_brewer(palette = "Set2", direction = -1, na.value = "white", guide = "none") +
  # scale_alpha_continuous(range = c(0, 1), name = "Frequency", trans = "log10") +
  scale_fill_distiller(palette = "RdPu", name = "Frequency", na.value = "white", direction = 1) +
  labs(
    x = "Cancer Type",
    y = "Gene"
  ) +
  theme(
    legend.position = "right",
    axis.text.x = element_text(size = 10, , angle = 90, vjust = 0.35, hjust = 1),
    axis.text.y = element_text(size = 10),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.title.x = element_text(size = 15, margin = margin(t = 20)),
    axis.title.y = element_text(size = 15, margin = margin(r = 20)),
    title = element_text(size = 15),
  )

pdf(glue("plots/{Sys.Date()}/gene-heatmap.pdf"), width = 17, height = 8.5)
gene_colormap_plot
dev.off()



# Mutation-cancer heatmap
# mutation_panel_hits_complete[, mutation := gsub(" ", "\n", mutation)]

mutation_panel_hits_complete[freq > 0, count := .N, by = mutation]

mutation_panel_hits_complete[, cancer_type := factor(cancer_type, levels = levels_cs_test)]

mutation_panel_hits_complete[, mutation := factor(mutation, levels = unique(mutation[order(count)]))]

# Top 10 only
mutation_panel_hits_complete <- mutation_panel_hits_complete[top_number == 10 | is.na(top_number)]

missing_ct_mut <- setdiff(final_panel$cancer_type, mutation_panel_hits_complete$cancer_type)
placeholder_rows_mut <- data.table(mutation = mutation_panel_hits_complete$mutation[1], cancer_type = missing_ct_mut, top_number = NA)
mutation_panel_hits_complete <- rbind(mutation_panel_hits_complete, placeholder_rows_mut, fill = TRUE)

mutation_colormap_plot <- ggplot(mutation_panel_hits_complete, aes(cancer_type, mutation, fill = freq)) +

  # # All tops
  # mutation_colormap_plot <- ggplot(mutation_panel_hits_complete, aes(cancer_type, mutation, fill = top_number, alpha = freq)) +
  geom_tile(color = "white") +
  theme_minimal() +
  # scale_fill_brewer(palette = "Set2", direction = -1, na.value = "white", guide = "none") +
  # scale_alpha_continuous(range = c(0, 1), name = "Frequency", trans = "log10") +
  scale_fill_distiller(palette = "RdPu", name = "Frequency (Log)", na.value = "white", trans = "log10", direction = 1) +
  labs(
    x = "Cancer Type",
    y = "Mutations"
  ) +
  theme(
    legend.position = "right",
    axis.text.x = element_text(size = 10, , angle = 90, vjust = 0.35, hjust = 1),
    axis.text.y = element_text(size = 10),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.title.x = element_text(size = 15, margin = margin(t = 20)),
    axis.title.y = element_text(size = 15, margin = margin(r = 20)),
    title = element_text(size = 15),
  )

pdf(glue("plots/{Sys.Date()}/mutation-heatmap.pdf"), width = 17, height = 8.5)
mutation_colormap_plot
dev.off()

cs_dis_plot_copy <- cs_dis_plot_copy + theme(
  legend.text = element_text(size = 20),   
  legend.title = element_text(size = 23),  
  legend.key.size = unit(1.25, "cm")  
) 

plot_legend <- get_legend(cs_dis_plot_copy)

mutation_colormap_copy <- copy(mutation_colormap_plot)
gene_colormap_copy <- copy(gene_colormap_plot)

mutation_colormap_copy <- mutation_colormap_copy + theme(
  legend.text = element_text(size = 20),   
  legend.title = element_text(size = 23),  
  legend.key.size = unit(1.25, "cm")  
) 

gene_colormap_copy <- gene_colormap_copy + theme(
  legend.text = element_text(size = 20),   
  legend.title = element_text(size = 23),  
  legend.key.size = unit(1.25, "cm")  
) 

mutation_freq_legend <- get_legend(mutation_colormap_copy)
gene_freq_legend <- get_legend(gene_colormap_copy)

gene_colormap_plot <- gene_colormap_plot + theme(
  legend.position = "none", 
  axis.text.x = element_text(size = 18), 
  axis.text.y = element_text(size = 18),
  axis.title.y = element_text(size = 25),
  axis.title.x = element_text(size = 25)
  )
mutation_colormap_plot <- mutation_colormap_plot + theme(
  legend.position = "none", 
  axis.text.x = element_text(size = 18), 
  axis.text.y = element_text(size = 18),
  axis.title.y = element_text(size = 25),
  axis.title.x = element_text(size = 25)
  )


# Add the heatmap below the combined plots - GENE
added_gene_heatmap <- plot_grid(
  final_plot_pc, gene_colormap_plot,
  nrow = 2, scale = 0.95
)

pdf(glue("plots/{Sys.Date()}/bar-and-gene-heat.pdf"), width = 17, height = 17)
grid.draw(added_gene_heatmap)
dev.off()

gene_legend <- plot_grid(
  plot_legend, gene_freq_legend,
  nrow = 2
)

# Add the legend to the right side of the combined plots
added_legend_gene <- plot_grid(
  added_gene_heatmap, gene_legend,
  ncol = 2, rel_widths = c(7, 1), scale = 0.95
)

pdf(glue("plots/{Sys.Date()}/bar-gene-heat-legend.pdf"), width = 17, height = 17)
grid.draw(added_legend_gene)
dev.off()

# Add the heatmap below the combined plots - MUTATION
added_mutation_heatmap <- plot_grid(
  final_plot_pc, mutation_colormap_plot,
  nrow = 2, scale = 0.95
)

pdf(glue("plots/{Sys.Date()}/bar-and-mutation-heat.pdf"), width = 17, height = 17)
grid.draw(added_mutation_heatmap)
dev.off()


mutation_legend <- plot_grid(
  plot_legend, mutation_freq_legend,
  nrow = 2
)

# Add the legend to the right side of the combined plots
added_legend_mutation <- plot_grid(
  added_mutation_heatmap, mutation_legend,
  ncol = 2, rel_widths = c(7, 1), scale = 0.95
)

pdf(glue("plots/{Sys.Date()}/bar-mutation-heat-legend.pdf"), width = 17, height = 17)
grid.draw(added_legend_mutation)
dev.off()


prepared_test <- prepare_for_cow(cs_test_plot)

prepared_test <- prepared_test + labs(y = "Percentage Coverage") + theme(
    axis.title.y = element_text(size = 25), 
    plot.title = element_text(size = 30), 
    axis.text.x = element_text(size = 18),
    axis.text.y = element_text(size = 18)
)

# Final plot - test cs panel + heatmap
# Mutation
cs_panel_heat_mut <- plot_grid(
  prepared_test, mutation_colormap_plot,
  nrow = 2, scale = 0.97, align = "v"
)

final_plot <- plot_grid(
  cs_panel_heat_mut, mutation_legend,
  ncol = 2, rel_widths = c(7, 1), scale = 0.97
)

pdf(glue("plots/{Sys.Date()}/mutation-bar-heat-plot.pdf"), width = 18, height = 17)
grid.draw(final_plot)
dev.off()

png(glue("plots/{Sys.Date()}/mutation-bar-heat-plot.png"), width = 1333, height = 1250)
grid.draw(final_plot)
dev.off()



# Gene

cs_panel_heat_gene <- plot_grid(
  prepared_test, gene_colormap_plot,
  nrow = 2, scale = 0.97, align = "v"
)

final_plot <- plot_grid(
  cs_panel_heat_gene, gene_legend,
  ncol = 2, rel_widths = c(7, 1), scale = 0.97
)

pdf(glue("plots/{Sys.Date()}/gene-bar-heat-plot.pdf"), width = 18, height = 17)
grid.draw(final_plot)
dev.off()

png(glue("plots/{Sys.Date()}/gene-bar-heat-plot.png"), width = 1333, height = 1250)
grid.draw(final_plot)
dev.off()


cv_fraction <- 0.6
TERT <- FALSE

# Split the dataset
cancer_data <- as.data.frame(tcga_available())

patient_data <- data.table()
patient_data_training <- data.table()
patient_data_test <- data.table()
top_nums <- c(10, 20, 50, 100)
top_nums_minus_first <- top_nums[-1]
set.seed(100)

training_pan_can <- data.table()
training_pan_can_tops <- data.table()
pan_can_coverage <- data.table()
pan_can_tops_coverage <- data.table()
pan_can_panel <- data.table()
pan_can_tops_panel <- data.table()

# Get all mutations and patients
for (i in seq_along(cancer_data$Study_Abbreviation)) {
  cancer_type <- cancer_data$Study_Abbreviation[i]
  test_mutations <- tcga_load(study = cancer_type)
  test_data <- as.data.table(test_mutations@data)
  current_data <- test_data[, c("Tumor_Sample_Barcode", "Hugo_Symbol", "HGVSc")]
  current_data[, cancer_type := cancer_type]
  current_data <- current_data[, .(patient = tstrsplit(Tumor_Sample_Barcode, "-")[[3]], mutation = paste(Hugo_Symbol, HGVSc, sep = ", "), cancer_type)]

  patient_data <- rbind(patient_data, current_data)

  unique_patients <- current_data[, unique(patient)]
  # 70:30 Split
  training_patients <- sample(unique_patients, round(length(unique_patients) * 0.7))

  training_data <- current_data[patient %in% training_patients]
  test_data <- current_data[!patient %in% training_patients]

  patient_data_training <- rbind(patient_data_training, training_data)
  patient_data_test <- rbind(patient_data_test, test_data)
}

patient_mutation <- copy(patient_data)
patient_data_training <- unique(add_tert(patient_data_training, TERT))
patient_data_test <- unique(add_tert(patient_data_test, TERT))

# Generate pan can panel
specific_patients <- copy(patient_data_training)

# sample fraction of training data
unique_specific_patients <- specific_patients[, unique(patient)]
patient_data_mutation_list <- patient_data_training[, .(mutations = list(mutation)), by = .(patient, cancer_type)]


for (run_num in 1:10) {
  set.seed(run_num)
  eighty_training_patients <- sample(unique_specific_patients, round(length(unique_specific_patients) * cv_fraction))

  eighty_training_data <- specific_patients[patient %in% eighty_training_patients]

  unique_patient_data <- patient_data_mutation_list[patient %in% eighty_training_patients]

  # panel <- greedy_panel_curator(eighty_training_data, max(top_nums))
  # new_row <- data.table(run = run_num, panel = list(panel$panel_mutation))
  # training_pan_can <- rbind(training_pan_can, new_row)

  tops_panel <- eighty_training_data[, .(n = .N), by = mutation][order(-n)][1:max(top_nums)]
  new_row <- data.table(run = run_num, panel = list(tops_panel$mutation))
  training_pan_can_tops <- rbind(training_pan_can_tops, new_row)
}

# training_data_save_pc <- copy(training_pan_can)
# training_data_save_pc$panel <- sapply(training_pan_can$panel, function(x) paste(x, collapse = ","))
# write.csv(training_data_save_pc, "csv-data/training_pan_can_panel.csv", row.names = FALSE)

all_training_dt <- copy(training_pan_can_tops)

# # Function to calculate unique mutation count for a given num_mutation
# calculate_unique_mutations <- function(data, num_mutation) {
#   only_n_mutations <- data[, .(mutation = unlist(lapply(panel, function(x) head(x, num_mutation)))), by = run]
#   unique_mutation_count <- only_n_mutations[, uniqueN(mutation)]
#   return(unique_mutation_count)
# }

# final_panel <- data.table()

current_data <- copy(all_training_dt)
top_mutations <- all_training_dt[, .(mutation = unlist(panel)), by = run]
unique_tops <- unique(top_mutations$mutation)
# if (length(unique_tops) < max(top_nums)) {
#   # Unlist all mutations, merge into one list then make unique
#   new_row <- data.table(panel = list(unique_tops))
#   final_panel <- rbind(final_panel, new_row, fill = TRUE)
# } else {
#   left <- 0
#   right <- max(top_nums)
#   target <- max(top_nums)
#   current_num <- 0
#   mid <- NA

#   while (left <= right) {
#     mid <- left + (right - left) %/% 2
#     current_num <- calculate_unique_mutations(current_data, mid)

#     if (abs(current_num - target) <= 5) {
#       break
#     }

#     if (current_num < target) {
#       left <- mid + 1
#     } else {
#       right <- mid - 1
#     }
#   }


#   only_n_mutations <- current_data[, panel := lapply(panel, function(x) head(x, mid))]
#   tops <- only_n_mutations[, .(mutation = unlist(panel)), by = run]
#   unique_tops <- unique(tops$mutation)

#   new_row <- data.table(panel = list(unique_tops))
#   final_panel <- rbind(final_panel, new_row, fill = TRUE)
# }

# all_unlisted <- all_training_dt[, .(mutation = unlist(panel)), by = setdiff(names(all_training_dt), "panel")]
# drop <- all_unlisted[, .(mutation)]
# drop_count <- drop[, .(count = .N), by = mutation][order(-count)]

# final_panel_raw <- drop_count[1:100]
# final_panel_raw <- final_panel_raw[, .(panel = list(mutation))]

tops_unlisted <- training_pan_can_tops[, .(mutation = unlist(panel)), by = setdiff(names(all_training_dt), "panel")]
droptop <- tops_unlisted[, .(mutation)]
droptop_count <- droptop[, .(count = .N), by = mutation][order(-count)]
final_panel_tops <- droptop_count[1:100]
final_panel_tops <- final_panel_tops[, .(panel = list(mutation))]

# pan_can_unopt_dis_plot <- test_and_plot_pc(final_panel_raw, patient_data_training, "Pan-Cancer Panel (Greedy, On Discovery)", "pan-can-greedy-cved-ondis", "Set2")
# pan_can_unopt_test_plot <- test_and_plot_pc(final_panel_raw, patient_data_test, "Pan-Cancer Panel (Greedy, On Test)", "pan-can-greedy-cved-ontest", "Set2")
# pan_can_opt_dis_plot <- test_and_plot_pc(final_panel, patient_data_training, "Pan-Cancer Panel (Greedy, On Discovery, Optimized)", "pan-can-greedy-cved-opt-ondis", "Set2")
# pan_can_opt_test_plot <- test_and_plot_pc(final_panel, patient_data_test, "Pan-Cancer Panel (Greedy, On Test, Optimized)", "pan-can-greedy-cved-opt-ontest", "Set2")
pc_tops_dis <- test_and_plot_pc(final_panel_tops, patient_data_training, "Pan-Cancer Panel (On Discovery)", "pan-can-tops-cved-ondis", "Set2")
pc_tops_test <- test_and_plot_pc(final_panel_tops, patient_data_test, "Pan-Cancer Panel (On Test)", "pan-can-tops-cved-ontest", "Set2")
pan_can_tops_dis_plot <- pc_tops_dis[[1]]
pan_can_tops_test_plot <- pc_tops_test[[1]]


cancer_incidence_raw <- fread("csv-data/cancer_incidence_acs.csv")
cancer_incidence <- cancer_incidence_raw[!is.na(new_cases)]

pc_cov_base <- copy(pc_tops_test[[2]])

pc_cov_base[cancer_incidence, on = .(cancer_type), covered := round(i.new_cases * coverage / 100)]

total_us_patients <- cancer_incidence[, sum(new_cases)]
pc_cov_acs <- pc_cov_base[, sum_covered := sum(covered, na.rm = TRUE), by = top_number]

covered_us_patients <- pc_cov_acs[, .(top_number, covered = sum_covered)]
covered_us_patients <- unique(covered_us_patients)
covered_us_patients[, fraction_covered := 100 * covered / total_us_patients]

for (i in seq_along(covered_us_patients$top_number)) {
  current_top_num <- covered_us_patients$top_number[i]
  pc_cov_acs[cancer_type == "PANCAN" & top_number == current_top_num, coverage := covered_us_patients$fraction_covered[i]]
}

# Plot
pc_cov_acs <- pc_cov_acs[order(-coverage, top_number)]

pc_cov_acs[, cancer_type := factor(cancer_type, levels = unique(cancer_type))]
pc_cov_acs[, top_number := factor(top_number, levels = unique(top_number))]

pc_cov_acs_plot <- ggplot(pc_cov_acs, aes(x = cancer_type, y = coverage, fill = top_number)) +
  geom_bar(stat = "identity", position = "identity", width = 0.9) +
  labs(title = "Pan-Cancer Panel Coverage (Actual Pan-Can Freq)", x = "Cancer Type", y = "Percentage Coverage", fill = "Number of\nTop Mutations") +
  theme_minimal() +
  theme(
    axis.text.y = element_text(size = 10),
    axis.text.x = element_text(size = 10, angle = 90, vjust = 0.25, hjust = 1)
  ) +
  scale_y_continuous(breaks = seq(0, 100, by = 10), limits = c(0, 100)) +
  scale_fill_brewer(palette = "Set2")

pdf(glue("plots/{Sys.Date()}/pan-can-panel-actual-freq.pdf"))
print(pc_cov_acs_plot)
dev.off()

# CAreful with these guys
plot_legend <- get_legend(pc_cov_acs_plot)

mutation_freq_legend <- get_legend(mutation_colormap_plot)
gene_freq_legend <- get_legend(gene_colormap_plot)

gene_colormap_plot <- gene_colormap_plot + theme(legend.position = "none")
mutation_colormap_plot <- mutation_colormap_plot + theme(legend.position = "none")
pc_cov_acs_plot <- pc_cov_acs_plot + theme(legend.position = "none")


gene_legend <- plot_grid(
  plot_legend, gene_freq_legend,
  nrow = 2
)


mutation_legend <- plot_grid(
  plot_legend, mutation_freq_legend,
  nrow = 2
)

prepared_test <- prepare_for_cow(pc_cov_acs_plot)

prepared_test <- prepared_test + labs(y = "Percentage Coverage") + theme(axis.title.y = element_text(size = 15), plot.title = element_text(size = 20))

# Final plot - test pc panel + heatmap
# Mutation
cs_panel_heat_mut <- plot_grid(
  prepared_test, mutation_colormap_plot,
  nrow = 2, scale = 0.97, align = "v"
)

final_plot <- plot_grid(
  cs_panel_heat_mut, mutation_legend,
  ncol = 2, rel_widths = c(7, 1), scale = 0.97
)

pdf(glue("plots/{Sys.Date()}/mutation-bar-heat-pancan-topten.pdf"), width = 18, height = 17)
grid.draw(final_plot)
dev.off()

# Gene

cs_panel_heat_gene <- plot_grid(
  prepared_test, gene_colormap_plot,
  nrow = 2, scale = 0.97, align = "v"
)

final_plot <- plot_grid(
  cs_panel_heat_gene, gene_legend,
  ncol = 2, rel_widths = c(7, 1), scale = 0.97
)

pdf(glue("plots/{Sys.Date()}/gene-bar-heat-pancan-topten.pdf"), width = 18, height = 17)
grid.draw(final_plot)
dev.off()



# prepare_for_cow <- function(plot) {
#   plot + theme(
#     legend.position = "none",
#     axis.title = element_blank(),
#     plot.title = element_text(size = 10),
#     panel.grid.major.x = element_blank(),
#     panel.grid.minor.y = element_blank()
#   )
# }

# plot_legend <- get_legend(pan_can_tops_test_plot)

# plots <- list(pan_can_unopt_dis_plot, pan_can_opt_dis_plot, pan_can_tops_dis_plot, pan_can_unopt_test_plot, pan_can_opt_test_plot, pan_can_tops_test_plot)
# prepared_plots <- lapply(plots, prepare_for_cow)

# combined_plots <- plot_grid(
#   plotlist = prepared_plots,
#   align = "vh", axis = "tb", ncol = 3
# )

# # Add the legend to the right side of the combined plots
# pan_can_cow <- plot_grid(
#   combined_plots, plot_legend,
#   ncol = 2, rel_widths = c(6, 1)
# )

# # Create common x and y labels
# y.grob <- textGrob("Percentage Coverage", gp = gpar(fontsize = 15), rot = 90)
# x.grob <- textGrob("Cancer Type", gp = gpar(fontsize = 15))

# # Arrange the combined plot and labels
# final_plot_pc <- arrangeGrob(pan_can_cow, left = y.grob, bottom = x.grob)

# pdf(glue("plots/{Sys.Date()}/pan-can-cowplot.pdf"), width = 17, height = 8.5)
# grid.draw(final_plot_pc)
# dev.off()


library(viridis)
library(tidytext)
library(data.table)
library(tidyverse)

cancer_data <- as.data.frame(tcga_available())

top_nums <- c(10, 20, 50, 100)

top_nums_minus_first <- top_nums[-1]

all_tops <- data.table(mutation = character(), n = integer(), top_number = integer())
all_unique_patient_data <- data.table(patient = character(), mutations = list())

coverage_within_cancer_type <- data.table(cancer_type = character(), coverage = numeric(), top_number = integer())
hits_within_cancer_type <- data.table(cancer_type = character(), hits = numeric(), frequency_of_hit = numeric(), perc_freq_hit = numeric(), top_number = integer())
cancer_specific_total_hits <- data.table()

for (i in seq_along(cancer_data$Study_Abbreviation)) {
  cancer_type <- cancer_data$Study_Abbreviation[i]
  test_mutations <- tcga_load(study = cancer_type)
  test_data <- as.data.table(test_mutations@data)
  patient_data <- test_data[, c("Tumor_Sample_Barcode", "Hugo_Symbol", "HGVSc")]

  # Add patient and mutation columns
  patient_mutation <- patient_data[, .(patient = tstrsplit(Tumor_Sample_Barcode, "-")[[3]], mutation = paste(Hugo_Symbol, HGVSc, sep = ", "))]

  # Group mutations into list per unique patient + count number of mutations per patient
  unique_patient_data <- patient_mutation[, .(mutations = list(mutation)), by = patient]

  for (num in top_nums) {
    top_nucleotides <- patient_mutation[, .(n = .N, top_number = num), by = mutation][order(-n)][1:num]

    within_cancer_patient_hits <- unique_patient_data[, .(
      patient,
      mutations,
      top_number = num,
      hits = sapply(mutations, function(x) length(intersect(x, top_nucleotides$mutation)))
    )]

    total_num_patients <- within_cancer_patient_hits[, .N]

    more_than_zero <- within_cancer_patient_hits[hits > 0, .N] / total_num_patients * 100

    new_row <- data.table(cancer_type = cancer_type, coverage = more_than_zero, top_number = num)

    coverage_within_cancer_type <- rbind(coverage_within_cancer_type, new_row)

    cumulated_hits <- within_cancer_patient_hits[, .(frequency_of_hit = .N), by = hits]
    new_hits <- cumulated_hits[, .(hits, frequency_of_hit, perc_freq_hit = frequency_of_hit / total_num_patients * 100, top_number = num, cancer_type = cancer_type)]

    hits_within_cancer_type <- rbind(hits_within_cancer_type, new_hits)

    all_tops <- rbind(all_tops, top_nucleotides)

    cancer_specific_total_hits <- rbind(cancer_specific_total_hits, within_cancer_patient_hits[, cancer_type := cancer_type])
  }

  all_unique_patient_data <- rbind(all_unique_patient_data, unique_patient_data)
}

coverage_within_cancer_type <- coverage_within_cancer_type[order(-coverage, -top_number)]

# Reorder the cancer_type factor based on the coverage
coverage_within_cancer_type[, cancer_type := factor(cancer_type, levels = unique(cancer_type))]
coverage_within_cancer_type[, top_number := factor(top_number, levels = unique(top_number))]

coverage_within_cancer_plot <- ggplot(coverage_within_cancer_type, aes(x = cancer_type, y = coverage, fill = top_number)) +
  geom_bar(stat = "identity", position = "identity", width = 0.9) +
  labs(title = "Cancer-specific Panel Coverage", x = "Cancer Type", y = "Percentage Coverage", fill = "Number of\nTop Mutations") +
  theme_minimal() +
  theme(
    axis.text.y = element_text(size = 7),
    axis.text.x = element_text(size = 7, angle = 90, vjust = 0.5, hjust = 1)
  ) +
  scale_y_continuous(breaks = seq(0, 100, by = 10), limits = c(0, 100))


pdf("plots/coverage_within_cancer_types_plot.pdf")
coverage_within_cancer_plot
dev.off()

# Make mutations unique
all_tops_unique <- all_tops[, .(n = sum(n)), by = .(mutation, top_number)][order(-n)]

# Duplicate mutations for each top_number, yea i actually dk why i did this, what's this for??
patient_data_all <- all_unique_patient_data[, top_number := 10]

patient_data_expanded <- patient_data_all[rep(seq_len(.N), each = 4)]

topnum_assignment <- rep(c(patient_data_all$top_number[1], top_nums_minus_first), length.out = nrow(patient_data_expanded))

patient_data_expanded[, top_number := topnum_assignment][order(patient)]

# # Find hits - this is the next step to change
# patient_hits <- patient_data_expanded[, hits := sapply(mutations, function(x) length(intersect(x, all_topten$mutation)))][order(-hits)]

patient_hits <- data.table(patient = character(), mutations = list(), top_number = numeric(), hits = integer())
pan_can_fractions <- numeric(length(top_nums))

for (i in seq_along(top_nums)) {
  num <- top_nums[i]

  subset_tops <- all_tops_unique[top_number == num]
  # Calculate hits for the current top_number
  current_patient_hits <- patient_data_expanded[top_number == num, .(
    patient,
    mutations,
    top_number = num,
    hits = sapply(mutations, function(x) length(intersect(x, subset_tops$mutation)))
  )]

  # Calculate fraction_covered for current_patient_hits
  pan_can_fractions[i] <- current_patient_hits[hits > 0, .N] / current_patient_hits[, .N] * 100

  # # Make master list, why do i need this again?
  # patient_hits <- rbind(patient_hits, current_patient_hits)
}

# Coverage per cancer type
coverage_per_cancer <- data.table(cancer_type = character(), coverage = numeric(), top_number = integer())

# append pan_cancer
pan_can_row <- data.table(cancer_type = "PAN-CAN", coverage = pan_can_fractions, top_number = top_nums)
coverage_per_cancer <- rbind(coverage_per_cancer, pan_can_row)

# append all cancer types
for (i in seq_along(cancer_data$Study_Abbreviation)) {
  cancer_type <- cancer_data$Study_Abbreviation[i]
  test_mutations <- tcga_load(study = cancer_type)
  test_data <- as.data.table(test_mutations@data)
  patient_data <- test_data[, c("Tumor_Sample_Barcode", "Hugo_Symbol", "HGVSc")]

  # Add patient and mutation columns
  patient_mutation <- patient_data[, .(patient = tstrsplit(Tumor_Sample_Barcode, "-")[[3]], mutation = paste(Hugo_Symbol, HGVSc, sep = ", "))]

  # Group mutations into list per unique patient + count number of mutations per patient
  unique_patient_data <- patient_mutation[, .(mutations = list(mutation)), by = patient]

  for (num in top_nums) {
    subset_tops <- all_tops_unique[top_number == num]

    per_cancer_patient_hits <- unique_patient_data[, .(
      patient,
      mutations,
      top_number = num,
      hits = sapply(mutations, function(x) length(intersect(x, subset_tops$mutation)))
    )]

    more_than_zero <- per_cancer_patient_hits[hits > 0, .N] / per_cancer_patient_hits[, .N] * 100

    new_row <- data.table(cancer_type = cancer_type, coverage = more_than_zero, top_number = num)

    coverage_per_cancer <- rbind(coverage_per_cancer, new_row)
  }
}

coverage_per_cancer <- coverage_per_cancer[order(-coverage, -top_number)]

# Reorder the cancer_type factor based on the coverage
coverage_per_cancer[, cancer_type := factor(cancer_type, levels = unique(cancer_type))]
coverage_per_cancer[, top_number := factor(top_number, levels = unique(top_number))]

# Extract the viridis colors and reverse them
num_levels <- length(unique(coverage_per_cancer$top_number))
reversed_viridis_colors <- rev(viridis(num_levels))

coverage_per_cancer_plot <- ggplot(coverage_per_cancer, aes(x = cancer_type, y = coverage, fill = top_number)) +
  geom_bar(stat = "identity", position = "identity", width = 0.9) +
  labs(title = "Percentage Coverage Pan-Can Panel", x = "Cancer Type", y = "Percentage Coverage", fill = "Number of\nTop Mutations") +
  theme_minimal() +
  theme(
    axis.text.y = element_text(size = 7),
    axis.text.x = element_text(size = 7, angle = 90, vjust = 0.5, hjust = 1)
  ) +
  scale_fill_manual(values = reversed_viridis_colors) +
  scale_y_continuous(breaks = seq(0, 100, by = 10), limits = c(0, 100))


pdf("plots/coverage_per_cancer_plot.pdf")
coverage_per_cancer_plot
dev.off()


# # Create plot
# ggplot(combined_data, aes(x = reorder_within(mutation, freq, cancer_type), y = freq)) + # nolint: object_usage_linter.
#   geom_bar(stat = "identity", fill = "blue", width = 0.9) +
#   geom_text(aes(label = round(freq, 3)), size = 1.5, hjust = 1.2, color = "white") +
#   coord_flip() +
#   labs(x = "Mutation", y = "Frequency") +
#   facet_wrap(~cancer_type, scales = "free") +
#   theme_minimal() +
#   scale_x_reordered() +
#   theme(
#     axis.text.y = element_text(size = 7),
#     axis.text.x = element_text(size = 5)
#   )

# # This is for a pie chart and cumulative frequency plot, only suitable for one top number
# grouped_data <- patient_hits[, .(count = .N), by = .(hits > 0)]

# grouped_data[, group := ifelse(hits == TRUE, "Hit", "No Hit")]
# grouped_data[, prop := count / sum(count) * 100]

# grouped_data[, ypos := -((cumsum(count) - 0.5 * count) - sum(count))]
# grouped_data[order(-count)]

# fraction_pie <- ggplot(grouped_data, aes(x = "", y = count, fill = group)) +
#   geom_bar(stat = "identity", width = 1) +
#   coord_polar("y", start = 0) +
#   theme_void() +
#   theme(legend.position = "none") +
#   geom_text(aes(y = ypos, label = paste0(group, "\n", round(prop, 1), "%")), color = "black", size = 6) +
#   scale_fill_manual(values = c("Hit" = "#73daff", "No Hit" = "#ff6961"))

# pdf("plots/fraction_covered_pie.pdf")

# fraction_pie
# dev.off()

# cumulative_freq_plot <- ggplot(patient_hits, aes(hits, y = 1 - after_stat(y))) +
#   stat_ecdf(geom = "step", color = "purple", linewidth = 1) +
#   labs(x = "Number of Hits", y = "%") +
#   scale_y_continuous(breaks = seq(0, 1, 0.1), labels = seq(0, 100, 10)) +
#   scale_x_continuous(breaks = seq(0, max(patient_hits$hits), 1)) +
#   theme_minimal()

# pdf("plots/cumulative_freq_plot.pdf")
# cumulative_freq_plot
# dev.off()


# Get all cancer types
cancer_types <- unique(cancer_data$Study_Abbreviation)

# Function to load mutations individually because TCGAmutations does not support loading all
load_mutations <- function(cancer_type) {
  mutations <- tcgaLoad(study = cancer_type) # nolint: object_usage_linter.
  mutation_data <- as.data.frame(mutations@data)
  mutation_data$cancer_type <- cancer_type
  return(mutation_data)
}

pan_cancer_list <- lapply(cancer_types, load_mutations)

pan_cancer_df <- bind_rows(pan_cancer_list)


nucleotide_changes <- pan_cancer_df[, c("Hugo_Symbol", "HGVSc", "cancer_type")]

tumor_samples <- pan_cancer_df[, c("Tumor_Sample_Barcode")]

num_tumors <- length(unique(tumor_samples))

# group_by HGVSc and add frequency using mutate, then ungroup
df <- nucleotide_changes %>%
  mutate(mutation = paste(Hugo_Symbol, HGVSc, sep = ", ")) %>%
  group_by(mutation) %>%
  mutate(freq = n()) %>%
  ungroup()

# sort in descending order, distinct removes duplicate rows
df_sorted <- df %>%
  arrange(desc(freq)) %>%
  distinct(mutation, .keep_all = TRUE)

df_merged <- df_sorted %>%
  select(mutation, freq, cancer_type)

df_merged_trunc <- df_merged[1:50, ]

# Convert the mutation column to a factor
df_merged_trunc$mutation <- factor(df_merged_trunc$mutation, levels = df_merged$mutation)
df_merged_trunc$freq <- df_merged_trunc$freq / num_tumors




all_mutations <- as.data.table(select(df_merged, mutation, freq))

patient_hits <- data.table(patient = character(), mutations = list(), top_number = numeric(), hits = integer())
pan_can_fractions <- numeric(length(top_nums))

for (i in seq_along(top_nums)) {
  num <- top_nums[i]

  current_top <- all_mutations[order(-freq)][1:num]

  # Calculate hits for the current top_number
  current_patient_hits <- all_unique_patient_data[, .(
    patient,
    mutations,
    top_number = num,
    hits = sapply(mutations, function(x) length(intersect(x, current_top$mutation)))
  )]

  # Calculate fraction_covered for current_patient_hits
  pan_can_fractions[i] <- current_patient_hits[hits > 0, .N] / current_patient_hits[, .N] * 100
}


pan_can_tops <- data.table(cancer_type = character(), coverage = numeric(), top_number = integer())

# append pan_cancer
pan_can_row <- data.table(cancer_type = "PAN-CAN", coverage = pan_can_fractions, top_number = top_nums)
pan_can_tops <- rbind(pan_can_tops, pan_can_row)

hits_pan_can <- data.table(cancer_type = character(), hits = numeric(), frequency_of_hit = numeric(), perc_freq_hit = numeric(), top_number = integer())
pancan_total_hits <- data.table()

top_nums <- c(10, 20, 50, 100)

for (i in seq_along(cancer_data$Study_Abbreviation)) {
  cancer_type <- cancer_data$Study_Abbreviation[i]
  test_mutations <- tcga_load(study = cancer_type)
  test_data <- as.data.table(test_mutations@data)
  patient_data <- test_data[, c("Tumor_Sample_Barcode", "Hugo_Symbol", "HGVSc")]

  # Add patient and mutation columns
  patient_mutation <- patient_data[, .(patient = tstrsplit(Tumor_Sample_Barcode, "-")[[3]], mutation = paste(Hugo_Symbol, HGVSc, sep = ", "))]

  # Group mutations into list per unique patient + count number of mutations per patient
  unique_patient_data <- patient_mutation[, .(mutations = list(mutation)), by = patient]


  for (num in top_nums) {
    current_top <- all_mutations[order(-freq)][1:num]

    per_cancer_patient_hits <- unique_patient_data[, .(
      patient,
      mutations,
      top_number = num,
      hits = sapply(mutations, function(x) length(intersect(x, current_top$mutation)))
    )]

    total_num_patients <- per_cancer_patient_hits[, .N]
    more_than_zero <- per_cancer_patient_hits[hits > 0, .N] / total_num_patients * 100

    new_row <- data.table(cancer_type = cancer_type, coverage = more_than_zero, top_number = num)

    pan_can_tops <- rbind(pan_can_tops, new_row)

    cumulated_hits <- per_cancer_patient_hits[, .(frequency_of_hit = .N), by = hits]
    new_hits <- cumulated_hits[, .(hits, frequency_of_hit, perc_freq_hit = frequency_of_hit / total_num_patients * 100, top_number = num, cancer_type = cancer_type)]

    hits_pan_can <- rbind(hits_pan_can, new_hits)

    pancan_total_hits <- rbind(pancan_total_hits, per_cancer_patient_hits[, cancer_type := cancer_type])
  }
}

pan_can_tops <- pan_can_tops[order(-coverage, -top_number)]

# Reorder the cancer_type factor based on the coverage
pan_can_tops[, cancer_type := factor(cancer_type, levels = unique(cancer_type))]
pan_can_tops[, top_number := factor(top_number, levels = unique(top_number))]

num_levels <- length(unique(pan_can_tops$top_number))
reversed_viridis_colors <- rev(viridis(num_levels))

coverage_within_cancer_plot <- ggplot(pan_can_tops, aes(x = cancer_type, y = coverage, fill = top_number)) +
  geom_bar(stat = "identity", position = "identity", width = 0.9) +
  labs(title = "Pan-Can Panel Coverage", x = "Cancer Type", y = "Percentage Coverage", fill = "Number of\nTop Mutations") +
  theme_minimal() +
  theme(
    axis.text.y = element_text(size = 7),
    axis.text.x = element_text(size = 7, angle = 90, vjust = 0.5, hjust = 1)
  ) +
  scale_y_continuous(breaks = seq(0, 100, by = 10), limits = c(0, 100)) +
  scale_fill_manual(values = reversed_viridis_colors)


pdf("plots/pan-can_panel_coverage.pdf")
coverage_within_cancer_plot
dev.off()


pancan_plot <- ggplot(data = df_merged_trunc, aes(x = mutation, y = freq)) +
  geom_bar(stat = "identity", fill = "blue") +
  geom_text(aes(label = round(freq, 3)), size = 3.5, hjust = 1.2, color = "white") +
  coord_flip() +
  labs(x = "Mutation", y = "Frequency", title = "Pan-cancer Frequency of Mutations") +
  theme_minimal()

pdf("plots/pancan_plot.pdf")

pancan_plot
dev.off()



# Diagnostic plots - distribution of hits
library(glue)

for (num in top_nums) {
  cancer_specific_hits <- hits_within_cancer_type[top_number == num]
  cancer_specific_hits[, hit_type := "Cancer Specific"]
  pan_can_hits <- hits_pan_can[top_number == num]
  pan_can_hits[, hit_type := "Pan-Cancer"]

  current_all_hits <- rbind(cancer_specific_hits, pan_can_hits)

  distribution_plot <- ggplot(current_all_hits, aes(x = hits, y = perc_freq_hit, colour = hit_type, group = hit_type)) +
    geom_point() +
    geom_line() +
    labs(title = glue("Distribution of Hits of Top {num} Mutations"), y = "% Frequency", x = "Number of Hits") +
    theme_minimal() +
    theme(
      axis.ticks.x = element_blank(),
      axis.ticks.y = element_blank(),
      panel.grid.major.y = element_blank(),
      legend.title = element_blank(),
    ) +
    facet_wrap(~cancer_type, scales = "free_x")

  pdf(glue("plots/diagnostic-plots/hit-distribution-top-{num}.pdf"))
  print(distribution_plot)
  dev.off()
}

# Targeted plots
targeted_hits_within <- hits_within_cancer_type[cancer_type %in% c("UCEC", "STAD", "LUSC", "BRCA")]
targeted_hits_pancan <- hits_pan_can[cancer_type %in% c("UCEC", "STAD", "LUSC", "BRCA")]
for (num in top_nums) {
  cancer_specific_hits <- targeted_hits_within[top_number == num]
  cancer_specific_hits[, hit_type := "Cancer Specific"]
  pan_can_hits <- targeted_hits_pancan[top_number == num]
  pan_can_hits[, hit_type := "Pan-Cancer"]

  current_all_hits <- rbind(cancer_specific_hits, pan_can_hits)

  distribution_plot <- ggplot(current_all_hits, aes(x = hits, y = perc_freq_hit, colour = hit_type, group = hit_type)) +
    geom_point() +
    geom_line() +
    labs(title = glue("Distribution of Hits of Top {num} Mutations"), y = "% Frequency", x = "Number of Hits") +
    theme_minimal() +
    theme(
      axis.ticks.x = element_blank(),
      axis.ticks.y = element_blank(),
      panel.grid.major.y = element_blank(),
      legend.title = element_blank(),
    ) +
    facet_wrap(~cancer_type, scales = "free_x")

  pdf(glue("plots/diagnostic-plots/targeted-hit-distribution-top-{num}.pdf"))
  print(distribution_plot)
  dev.off()
}

# Density
total_targeted_hits_within <- cancer_specific_total_hits[cancer_type %in% c("UCEC", "STAD", "LUSC", "BRCA")]
total_targeted_hits_pancan <- pancan_total_hits[cancer_type %in% c("UCEC", "STAD", "LUSC", "BRCA")]

for (num in top_nums) {
  cancer_specific_hits <- total_targeted_hits_within[top_number == num]
  cancer_specific_hits[, hit_type := "Cancer Specific"]
  pan_can_hits <- total_targeted_hits_pancan[top_number == num]
  pan_can_hits[, hit_type := "Pan-Cancer"]

  current_all_hits <- rbind(cancer_specific_hits, pan_can_hits)

  distribution_plot <- ggplot(current_all_hits, aes(x = hits, fill = hit_type, group = hit_type)) +
    geom_density(alpha = 0.6) +
    labs(title = glue("Distribution of Hits of Top {num} Mutations"), x = "Number of Hits") +
    theme_minimal() +
    theme(
      axis.ticks.x = element_blank(),
      axis.ticks.y = element_blank(),
      panel.grid.major.y = element_blank(),
      legend.title = element_blank(),
    ) +
    facet_wrap(~cancer_type, scales = "free")

  pdf(glue("plots/diagnostic-plots/density-hit-distribution-top-{num}.pdf"))
  print(distribution_plot)
  dev.off()
}



# Coverage within cancer types
cancer_data <- as.data.frame(tcga_available())

top_nums <- c(10, 20, 50, 100)

top_nums_minus_first <- top_nums[-1]

all_tops <- data.table(mutation = character(), n = integer(), top_number = integer())
all_unique_patient_data <- data.table(patient = character(), mutations = list())

coverage_within_cancer_type <- data.table(cancer_type = character(), coverage = numeric(), top_number = integer())


cancer_type <- "UCEC"
test_mutations <- tcga_load(study = cancer_type)
test_data <- as.data.table(test_mutations@data)
patient_data <- test_data[, c("Tumor_Sample_Barcode", "Hugo_Symbol", "HGVSc")]

# Add patient and mutation columns
patient_mutation <- patient_data[, .(patient = tstrsplit(Tumor_Sample_Barcode, "-")[[3]], mutation = paste(Hugo_Symbol, HGVSc, sep = ", "))]

# Group mutations into list per unique patient + count number of mutations per patient
unique_patient_data <- patient_mutation[, .(mutations = list(mutation)), by = patient]

for (num in top_nums) {
  top_nucleotides <- patient_mutation[, .(n = .N, top_number = num), by = mutation][order(-n)][1:num]

  within_cancer_patient_hits <- unique_patient_data[, .(
    patient,
    mutations,
    top_number = num,
    hits = sapply(mutations, function(x) length(intersect(x, top_nucleotides$mutation)))
  )]

  more_than_zero <- within_cancer_patient_hits[, sum(hits)] / within_cancer_patient_hits[, .N] * 100

  new_row <- data.table(cancer_type = cancer_type, coverage = more_than_zero, top_number = num)

  coverage_within_cancer_type <- rbind(coverage_within_cancer_type, new_row)

  all_tops <- rbind(all_tops, top_nucleotides)
}

# pan-can panel
all_mutations <- data.table(mutation = character(), n = integer())

for (i in seq_along(cancer_data$Study_Abbreviation)) {
  # Load TCGA mutation data
  cancer_type <- cancer_data$Study_Abbreviation[i]
  cancer_name <- gsub("_", " ", cancer_data$Study_Name[i])
  mutations <- tcgaLoad(study = cancer_type)
  mutation_data <- as.data.table(mutations@data)
  nucleotide_changes <- mutation_data[, c("Hugo_Symbol", "HGVSc")]

  nucleotide_changes[, mutation := paste(Hugo_Symbol, HGVSc, sep = ", ")]
  mutation_counts <- nucleotide_changes[, .(n = .N), by = mutation]

  all_mutations <- rbind(all_mutations, mutation_counts)
}

all_mutations <- all_mutations[, .(n = sum(n)), by = mutation][order(-n)]

patient_hits <- data.table(patient = character(), mutations = list(), top_number = numeric(), hits = integer())
pan_can_fractions <- numeric(length(top_nums))

for (i in seq_along(top_nums)) {
  num <- top_nums[i]

  current_top <- all_mutations[order(-n)][1:num]

  # Calculate hits for the current top_number
  current_patient_hits <- unique_patient_data[, .(
    patient,
    mutations,
    top_number = num,
    hits = sapply(mutations, function(x) length(intersect(x, current_top$mutation)))
  )]

  # Calculate fraction_covered for current_patient_hits
  pan_can_fractions[i] <- current_patient_hits[, sum(hits)] / current_patient_hits[, .N]
}


pan_can_tops <- data.table(cancer_type = character(), coverage = numeric(), top_number = integer())

# append pan_cancer
pan_can_row <- data.table(cancer_type = "UCEC-pan", coverage = pan_can_fractions, top_number = top_nums)


coverage_within_cancer_type <- rbind(coverage_within_cancer_type, pan_can_row)

coverage_within_cancer_type <- coverage_within_cancer_type[order(-coverage, -top_number)]

# Reorder the cancer_type factor based on the coverage
coverage_within_cancer_type[, cancer_type := factor(cancer_type, levels = unique(cancer_type))]
coverage_within_cancer_type[, top_number := factor(top_number, levels = unique(top_number))]

coverage_within_cancer_plot <- ggplot(coverage_within_cancer_type, aes(x = cancer_type, y = coverage, fill = top_number)) +
  geom_bar(stat = "identity", position = "identity", width = 0.9) +
  labs(title = "weird", x = "Cancer Type", y = "Percentage Coverage", fill = "Number of\nTop Mutations") +
  theme_minimal() +
  theme(
    axis.text.y = element_text(size = 7),
    axis.text.x = element_text(size = 7, angle = 90, vjust = 0.5, hjust = 1)
  )


pdf("plots/coverage_within_cancer_types_ucec.pdf")
coverage_within_cancer_plot
dev.off()

patient_mutation_count <- patient_mutation[, .(mutation_count = .N), by = patient][order(mutation_count)]

# This is currently only for UCEC
mutation_histogram <- ggplot(patient_mutation_count, aes(x = mutation_count)) +
  geom_histogram(fill = "#69b3a2", color = "#e9ecef", alpha = 0.8) +
  scale_x_log10() +
  labs(y = "Number of Patients", x = "Mutation Count")
pdf("plots/mutation_count_histogram.pdf")
mutation_histogram
dev.off()

# pan-can panel
all_mutations <- data.table(mutation = character(), n = integer())

for (i in seq_along(cancer_data$Study_Abbreviation)) {
  # Load TCGA mutation data
  cancer_type <- cancer_data$Study_Abbreviation[i]
  cancer_name <- gsub("_", " ", cancer_data$Study_Name[i])
  mutations <- tcgaLoad(study = cancer_type)
  mutation_data <- as.data.table(mutations@data)
  nucleotide_changes <- mutation_data[, c("Hugo_Symbol", "HGVSc")]

  nucleotide_changes[, mutation := paste(Hugo_Symbol, HGVSc, sep = ", ")]
  mutation_counts <- nucleotide_changes[, .(n = .N), by = mutation]

  all_mutations <- rbind(all_mutations, mutation_counts)
}

all_mutations <- all_mutations[, .(n = sum(n)), by = mutation][order(-n)]

patient_hits <- data.table(patient = character(), mutations = list(), top_number = numeric(), hits = integer())
pan_can_fractions <- numeric(length(top_nums))

for (i in seq_along(top_nums)) {
  num <- top_nums[i]

  current_top <- all_mutations[order(-n)][1:num]

  # Calculate hits for the current top_number
  current_patient_hits <- unique_patient_data[, .(
    patient,
    mutations,
    top_number = num,
    hits = sapply(mutations, function(x) length(intersect(x, current_top$mutation)))
  )]

  # Calculate fraction_covered for current_patient_hits
  pan_can_fractions[i] <- current_patient_hits[hits > 0, .N] / current_patient_hits[, .N] * 100
}


pan_can_tops <- data.table(cancer_type = character(), coverage = numeric(), top_number = integer())

# append pan_cancer
pan_can_row <- data.table(cancer_type = "UCEC-pan", coverage = pan_can_fractions, top_number = top_nums)
pan_can_tops <- rbind(pan_can_tops, pan_can_row)

top_nums <- c(10, 20, 50, 100)


cancer_type <- "UCEC"
test_mutations <- tcga_load(study = cancer_type)
test_data <- as.data.table(test_mutations@data)
patient_data <- test_data[, c("Tumor_Sample_Barcode", "Hugo_Symbol", "HGVSc")]

# Add patient and mutation columns
patient_mutation <- patient_data[, .(patient = tstrsplit(Tumor_Sample_Barcode, "-")[[3]], mutation = paste(Hugo_Symbol, HGVSc, sep = ", "))]

# Group mutations into list per unique patient + count number of mutations per patient
unique_patient_data <- patient_mutation[, .(mutations = list(mutation)), by = patient]


for (num in top_nums) {
  current_top <- all_mutations[order(-n)][1:num]

  per_cancer_patient_hits <- unique_patient_data[, .(
    patient,
    mutations,
    top_number = num,
    hits = sapply(mutations, function(x) length(intersect(x, current_top$mutation)))
  )]
  more_than_zero <- per_cancer_patient_hits[hits > 0, .N] / per_cancer_patient_hits[, .N] * 100

  new_row <- data.table(cancer_type = cancer_type, coverage = more_than_zero, top_number = num)

  pan_can_tops <- rbind(pan_can_tops, new_row)
}

pan_can_tops <- pan_can_tops[order(-coverage, -top_number)]

# Reorder the cancer_type factor based on the coverage
pan_can_tops[, cancer_type := factor(cancer_type, levels = unique(cancer_type))]
pan_can_tops[, top_number := factor(top_number, levels = unique(top_number))]

num_levels <- length(unique(pan_can_tops$top_number))
reversed_viridis_colors <- rev(viridis(num_levels))

coverage_within_cancer_plot <- ggplot(pan_can_tops, aes(x = cancer_type, y = coverage, fill = top_number)) +
  geom_bar(stat = "identity", position = "identity", width = 0.9) +
  labs(title = "Pan-Can Panel Coverage", x = "Cancer Type", y = "Percentage Coverage", fill = "Number of\nTop Mutations") +
  theme_minimal() +
  theme(
    axis.text.y = element_text(size = 7),
    axis.text.x = element_text(size = 7, angle = 90, vjust = 0.5, hjust = 1)
  ) +
  scale_y_continuous(breaks = seq(0, 100, by = 10), limits = c(0, 100)) +
  scale_fill_manual(values = reversed_viridis_colors)


pdf("plots/pan-can_panel_ucec.pdf")
coverage_within_cancer_plot
dev.off()

