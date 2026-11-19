# Run this script from the IA-week project folder.

invisible(list2env(local({
  # Load the official G2F 2023 phenotypic data.
  phenotype_file = file.path(
    "G2F", "GenomesToFields_G2F_data_2023",
    "a._2023_phenotypic_data", "g2f_2023_phenotypic_clean_data.csv"
  )
  phenotypes = read.csv(phenotype_file, check.names = FALSE)

  count_distinct = function(x) {
    valid = !is.na(x) & !grepl("^[[:space:]]*$", x, useBytes = TRUE)
    length(unique(x[valid]))
  }

  # Summarize experiments described in the phenotypic guide.
  descriptions = c(
    G2F_Main = "Main field trial with check hybrids linking trial locations.",
    YS = "Yellow Stripe trial with reference hybrids evaluated across years.",
    HIP_Hybrid = "High-intensity phenotyping trial with 22 hybrids.",
    HIP_Inbred = "High-intensity phenotyping trial with 22 inbred lines."
  )

  experiments = intersect(names(descriptions), phenotypes[["Experiment"]])
  experiment_summary = do.call(rbind, lapply(experiments, function(experiment) {
    data = phenotypes[which(phenotypes[["Experiment"]] == experiment), ]
    data.frame(
      Experiment = experiment,
      Plots = nrow(data),
      Pedigrees = count_distinct(data[["Pedigree"]]),
      Field_Locations = count_distinct(data[["Field-Location"]]),
      Description = descriptions[[experiment]]
    )
  }))

  # Use G2F_Main because it has the most distinct pedigrees
  # and includes field evaluations across multiple locations.
  main = subset(phenotypes, Experiment == "G2F_Main")

  # Identify check pedigrees from the original Family labels.
  check_pedigrees = unique(main[["Pedigree"]][main[["Family"]] %in% "CHECK"])

  # Remove plots with missing grain yield or a discard flag.
  discard = tolower(trimws(main[["Plot Discarded [enter 'yes' or blank]"]]))
  missing_yield = is.na(main[["Grain Yield (bu/A)"]])
  discarded = discard %in% "yes"
  keep = !missing_yield & !discarded

  plot_filter = data.frame(
    Original_Plots = nrow(main),
    Missing_Yield = sum(missing_yield),
    Flagged_For_Discard = sum(discarded),
    Both = sum(missing_yield & discarded),
    Removed = sum(!keep),
    Retained = sum(keep)
  )

  # Keep the plot-level fields needed for grain yield analysis.
  columns = c(
    "Field-Location", "Replicate", "Block", "Plot_ID",
    "Range", "Pass", "Pedigree", "Family"
  )
  yield_data = main[keep, columns, drop = FALSE]

  # Express grain yield in kg/ha, a more familiar unit for our team.
  corn_bu_ac_to_kg_ha = (56 * 0.45359237) / 0.40468564224
  yield_data$GY = main[["Grain Yield (bu/A)"]][keep] * corn_bu_ac_to_kg_ha

  # Verify that each plot is unique within its field location.
  if (anyDuplicated(yield_data[, c("Field-Location", "Plot_ID")])) {
    stop("Duplicate plot IDs found within a field location.")
  }

  # Read parental IDs from pedigrees written as Parent1/Parent2.
  pedigree = yield_data[["Pedigree"]]
  valid = !is.na(pedigree) & grepl("^[^/]+/[^/]+$", pedigree, useBytes = TRUE)
  parents = strsplit(pedigree[valid], "/", fixed = TRUE)

  yield_data$DH1 = NA_character_
  yield_data$DH2 = NA_character_
  yield_data$DH1[valid] = trimws(vapply(parents, `[`, character(1L), 1L))
  yield_data$DH2[valid] = trimws(vapply(parents, `[`, character(1L), 2L))

  # Keep unresolved pedigrees available for plot-level analysis.
  unresolved = is.na(yield_data$DH1) | is.na(yield_data$DH2) |
    yield_data$DH1 %in% "" | yield_data$DH2 %in% ""

  yield_data[unresolved, c("DH1", "DH2")] = NA_character_
  parental_review = as.data.frame(
    table(Pedigree = pedigree[unresolved]),
    responseName = "Plots"
  )

  yield_data = yield_data[, c(
    "Field-Location", "Replicate", "Block", "Plot_ID", "Range", "Pass",
    "Pedigree", "DH1", "DH2", "Family", "GY"
  )]
  rownames(yield_data) = NULL

  # Count blocks and plots within each location and replicate.
  design_summary = aggregate(
    list(Blocks = yield_data$Block, Plots = yield_data$Plot_ID),
    by = list(
      Field_Location = yield_data[["Field-Location"]],
      Replicate = yield_data$Replicate
    ),
    FUN = count_distinct
  )

  # Count distinct replicates once per observed pedigree-location pair.
  replication = aggregate(
    list(Replicates = yield_data$Replicate),
    by = list(
      Field_Location = yield_data[["Field-Location"]],
      Pedigree = yield_data$Pedigree
    ),
    FUN = count_distinct
  )
  replication$Entry_Type = ifelse(
    replication$Pedigree %in% check_pedigrees, "Check", "Other"
  )
  replication_summary = with(replication, as.data.frame(
    table(Entry_Type, Replicates),
    responseName = "Pedigree_Location_Pairs"
  ))

  # Use the final analysis data to illustrate the field layout at DEH1.
  deh1 = subset(yield_data, `Field-Location` == "DEH1")
  deh1$Entry_Type = ifelse(deh1$Pedigree %in% check_pedigrees, "Check", "Other")

  # Count plots and distinct pedigrees separately.
  summarize_entries = function(data) {
    is_check = data$Entry_Type == "Check"
    data.frame(
      Plots = nrow(data),
      Check_Plots = sum(is_check),
      Other_Plots = sum(!is_check),
      Check_Pedigrees = count_distinct(data$Pedigree[is_check]),
      Other_Pedigrees = count_distinct(data$Pedigree[!is_check])
    )
  }

  # A replicate contains several smaller blocks.
  by_replicate = do.call(rbind, lapply(split(deh1, deh1$Replicate), function(x) {
    block_sizes = table(x$Block)
    data.frame(
      Replicate = x$Replicate[1L],
      Blocks = length(block_sizes),
      summarize_entries(x),
      Min_Plots_Per_Block = min(block_sizes),
      Max_Plots_Per_Block = max(block_sizes)
    )
  }))

  # Block numbers are reused across replicates.
  by_block = do.call(rbind, lapply(
    split(deh1, list(deh1$Replicate, deh1$Block), drop = TRUE),
    function(x) {
      data.frame(
        Replicate = x$Replicate[1L],
        Block = x$Block[1L],
        summarize_entries(x)
      )
    }
  ))
  by_block = by_block[order(by_block$Replicate, by_block$Block), ]

  # Range and Pass show where each plot sits in the field.
  example_block = subset(
    deh1, Replicate == 1 & Block == 1,
    select = c(Range, Pass, Plot_ID, Pedigree, Entry_Type, GY)
  )
  example_block = example_block[order(
    example_block$Range, suppressWarnings(as.integer(example_block$Pass))
  ), ]
  rownames(example_block) = NULL

  # Guide students from the whole trial to an individual plot.
  layout_guide = data.frame(
    Level = c("Field location", "Replicate", "Block", "Plot", "Pedigree"),
    Example = c("DEH1", "Replicate 1", "Block 1", "One row in the block table",
                "The entry planted in that plot"),
    How_to_Read = c(
      "One trial location in the 2023 dataset.",
      "A large section of the trial containing several blocks.",
      "A smaller group of plots within a replicate.",
      "The unit where grain yield was measured.",
      "A genetic entry that can appear in several plots."
    )
  )

  list(
    yield_data = yield_data,
    report_tables = list(
      experiments = experiment_summary,
      plot_filter = plot_filter,
      parental_review = parental_review,
      design = design_summary,
      replication = replication_summary,
      design_example = list(
        guide = layout_guide,
        by_replicate = by_replicate,
        by_block = by_block,
        example_block = example_block
      )
    )
  )
}), envir = environment()))

# Review the prepared dataset.
print(dim(yield_data))
print(report_tables$plot_filter, row.names = FALSE)

# Read the DEH1 layout from replicate to block to individual plot.
cat("\nDEH1 | How to read the field layout\n")
print(report_tables$design_example$guide, row.names = FALSE)

cat("\nDEH1 | All replicates\n")
print(report_tables$design_example$by_replicate, row.names = FALSE)

cat("\nDEH1 | Blocks within replicate 1\n")
print(subset(report_tables$design_example$by_block, Replicate == 1),
      row.names = FALSE)

cat("\nDEH1 | Plots in replicate 1, block 1 | GY in kg/ha\n")
print(report_tables$design_example$example_block, row.names = FALSE)

library(asreml)

# Keep a pedigree-to-parent lookup for the genomic analysis.
parent_lookup = unique(yield_data[, c("Pedigree", "DH1", "DH2")])

if (anyDuplicated(parent_lookup$Pedigree)) {
  stop("Some pedigrees have conflicting parental IDs.")
}

# Fit the first-stage model separately at each field location.
local_results = lapply(
  split(yield_data, yield_data[["Field-Location"]]),
  function(plots) {
    location = plots[["Field-Location"]][1L]
    n_pedigrees = length(unique(plots$Pedigree))
    message("Fitting ", location, " | ", nrow(plots), " plots")

    plots$Pedigree = factor(plots$Pedigree)
    plots$Replicate = factor(plots$Replicate)
    plots$Block = factor(plots$Block)

    diagnostic = data.frame(
      Field_Location = location,
      Plots = nrow(plots),
      Pedigrees = n_pedigrees,
      Converged = NA,
      Estimable_Means = NA_integer_,
      Status = "Fit error",
      Message = ""
    )

    tryCatch({
      fit = asreml(
        fixed = GY ~ Pedigree + Replicate,
        random = ~ Replicate:Block,
        data = plots
      )

      diagnostic$Converged = isTRUE(fit$converge)

      if (!diagnostic$Converged) {
        diagnostic$Status = "Not converged"
        return(list(blue = NULL, diagnostic = diagnostic))
      }

      # Request one adjusted mean per pedigree.
      diagnostic$Status = "Prediction error"
      pvals = predict(fit, classify = "Pedigree")$pvals
      valid = pvals$status %in% "Estimable" &
        is.finite(pvals$predicted.value) &
        is.finite(pvals$std.error)

      diagnostic$Estimable_Means = sum(valid)
      diagnostic$Status = if (all(valid) && sum(valid) == n_pedigrees) {
        "OK"
      } else {
        "Check prediction coverage"
      }

      blue = data.frame(
        Field_Location = rep(location, sum(valid)),
        Pedigree = as.character(pvals$Pedigree[valid]),
        BLUE = pvals$predicted.value[valid],
        SE = pvals$std.error[valid]
      )

      list(blue = blue, diagnostic = diagnostic)
    }, error = function(e) {
      diagnostic$Message = conditionMessage(e)
      list(blue = NULL, diagnostic = diagnostic)
    })
  }
)

# Combine diagnostics and estimable means across locations.
fit_summary = do.call(rbind, lapply(local_results, `[[`, "diagnostic"))
rownames(fit_summary) = NULL
print(fit_summary, row.names = FALSE)

blue_data = do.call(rbind, lapply(local_results, `[[`, "blue"))

if (is.null(blue_data) || !nrow(blue_data)) {
  stop("No adjusted means were obtained. Check fit_summary.")
}

# Attach the parental IDs and keep one row per pedigree and location.
parent_row = match(blue_data$Pedigree, parent_lookup$Pedigree)
blue_data$DH1 = parent_lookup$DH1[parent_row]
blue_data$DH2 = parent_lookup$DH2[parent_row]

blue_data = blue_data[, c(
  "Field_Location", "Pedigree", "DH1", "DH2", "BLUE", "SE"
)]
rownames(blue_data) = NULL

stopifnot(!anyDuplicated(blue_data[, c("Field_Location", "Pedigree")]))
print(dim(blue_data))
print(head(blue_data), row.names = FALSE)