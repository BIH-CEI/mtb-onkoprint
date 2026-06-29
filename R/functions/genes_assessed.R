get_gene_list_names_from_methdos <- function(Methods, mutation_type) {
  list_of_methods <- str_split_1(Methods, pattern = ", ")

  # Filters the list of platforms and returns the names of all gene lists needed
  panel_list <- filter(
    panel_per_platform_lookup,
    .data[[mutation_type]] == "yes",
    Platform %in% list_of_methods
  )
  if(mutation_type == "Fusions detected"){
    Gene_list_to_use <- panel_list[["Gene list Fusions"]]
  } else{
    Gene_list_to_use <- panel_list[["Gene list SNVs"]]
  }
    
  return(Gene_list_to_use)
}

get_fusion_gene_list_names <- function(Methods) {
  get_gene_list_names_from_methdos(Methods, mutation_type = "Fusions detected")
}

get_SNV_gene_list_names <- function(Methods) {
  get_gene_list_names_from_methdos(Methods, mutation_type = "SNVs detected")
}

get_genes_from_gene_lists <- function(gene_lists_to_use) {
  assessed_genes <- select(gene_panels, gene_symbol, all_of(gene_lists_to_use)) |>
    pivot_longer(!gene_symbol, names_to = "panel", values_to = "genes") |>
    filter(genes != "not assessed") |>
    pull(genes) |>
    unique()
  return(assessed_genes)
}

get_list_of_assessed_genes_per_patient <- function(all_methods, mutation_type = "SNV") {
  if (mutation_type == "SNV") {
    gene_lists_to_use <- get_SNV_gene_list_names(Methods = all_methods)
  }

  if (mutation_type == "Fusion") {
    gene_lists_to_use <- get_fusion_gene_list_names(Methods = all_methods)
  }
  if ("all_assessed" %in% gene_lists_to_use) {
    return(c("all_assessed"))
  }
  # If no genes were assessed, return an empty vector
  if (length(gene_lists_to_use) == 0) {
    return(NULL)
  }
  assessed_genes <- get_genes_from_gene_lists(gene_lists_to_use)
  return(assessed_genes)
}

check_if_gene_was_assessed <- function(patient_index, gene_symbol) {
  assessed_genes <- genes_assessed_per_patient[[patient_index]]

  if (is.null(assessed_genes)) {
    return(FALSE)
  } else if ("all_assessed" %in% assessed_genes) {
    return(TRUE)
  } else {
    gene_assessed_bool <- gene_symbol %in% assessed_genes
    return(gene_assessed_bool)
  }
}

set_assessment_value <- function(SNV_value, patient_index, col_name) {
  gene_symbol <- str_remove(col_name, "SNV_")
  if (!is.na(SNV_value)) {
    return(SNV_value)
  } else {
    if (patient_index == 4) {
      gene_assessed_bool <- check_if_gene_was_assessed(patient_index, gene_symbol)
    }
    gene_assessed_bool <- check_if_gene_was_assessed(patient_index, gene_symbol)
  }
  if (gene_assessed_bool) {
    return("assessed")
  } else {
    return("not assessed")
  }
}


get_patients_assessed_per_gene <- function(list_of_genes_assessed_per_patient){
  # Adapted based on code suggestion from ChatGPT:
  
  # All genes assessed
  possible_values <- list_of_genes_assessed_per_patient |>
    unlist() |>
    unique() |>
    setdiff("all_assessed") |>
    sort()
  
  
  df <- imap_dfr(list_of_genes_assessed_per_patient, \(vals, id) {
    tibble(ID = id) |>
      bind_cols(
        setNames(
          lapply(
            possible_values,
            \(v) as.integer(v %in% vals || "all_assessed" %in% vals)
          ),
          possible_values
        ) |>
          as_tibble()
      ) |>
      mutate(all_assessed = as.integer("all_assessed" %in% vals))
  })
  return(df)
}

# Testing if code works. Unit-test like

# get_all_SNV_platforms("CGP, OFA")
# get_all_fusion_platforms("CGP, OFA")

# # TODO Test with incorrect/missing methods
#
# # These should give ABCD
# get_list_of_assessed_genes_per_patient("OFA")
# get_list_of_assessed_genes_per_patient("Archer, Basic-NGS")
# get_list_of_assessed_genes_per_patient("CGP, Basic-NGS")
#
# # These should crash
# get_list_of_assessed_genes_per_patient(all_methods = "Archer", mutation_type = "SNV")
# # And this one BC
# get_list_of_assessed_genes_per_patient(all_methods = "Archer", mutation_type = "Fusion")
#
# # Not sure
# get_list_of_assessed_genes_per_patient(all_methods = "Qiagen DHS-3501Z (keine Fusionsanalyse), Archer", mutation_type = "SNV")

# Should work for some but not all patients
# patient_index = 4
# gene_symbol = "A"
# check_if_gene_was_assessed(patient_index, gene_symbol)
