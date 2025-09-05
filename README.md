# Educational Resilience Profiles in Greece (PISA data)
Code for replicating the analyses in "Discovering Profiles of Resilient Students in PISA: A Distance-Based Approach to Clustering Mixed-Type Educational Data"

---
title: "Educational Resilience Clustering Analysis using Akhanli & Hennig Framework"
subtitle: "Mixed-Type Distance Clustering for PISA Disadvantaged Students with Variable Weighting"
author: "Your Name"
date: "`r Sys.Date()`"
output:
  html_document:  
    toc: true
    toc_depth: 3
    toc_float: true
    number_sections: true
    theme: united
    highlight: tango
    code_folding: show
  pdf_document:
    toc: true
    number_sections: true
params:
  data_dir: "/Users/amarkos/PISA_Data/"
  target_countries: "GRC"
  disadvantaged_threshold: 0.25
  k_range_min: 2
  k_range_max: 8
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(
  echo = TRUE,
  warning = FALSE,
  message = FALSE,
  fig.width = 12,
  fig.height = 8,
  cache = TRUE,
  comment = NA
)
options(scipen = 999)
```


# Executive Summary

This analysis applies the Akhanli & Hennig (2020) composite clustering validity framework to identify educational resilience profiles among disadvantaged PISA students. The analysis uses mixed-type distance methods with variable weighting to emphasize key educational factors for intervention purposes.


# Library Setup and Configuration

```{r libraries, message=FALSE, warning=FALSE}
cat("=== INITIALIZING CLUSTERING ENVIRONMENT ===\n")

library(haven)
library(data.table)
library(dtplyr)
library(labelled)
library(psych)
library(ggplot2)
library(tidyr)
library(dplyr)
library(gridExtra)
library(viridis)
library(scales)
library(fpc)
library(mclust)
library(fastDummies)
library(mice)
library(corrplot)
library(factoextra)
library(cluster)
library(MASS)
library(moments)

if (!requireNamespace("manydist", quietly = TRUE)) {
  stop("manydist package required but not available. Install with: install.packages('manydist')")
}
library(manydist)

setDTthreads(0)
set.seed(42)
cat("✓ Analysis environment configured\n")
```

```{r configuration}
cat("\n=== CONFIGURING ANALYSIS PARAMETERS ===\n")

data_dir <- params$data_dir
target_countries <- params$target_countries
DISADVANTAGED_THRESHOLD <- params$disadvantaged_threshold

path <- function(fname) file.path(data_dir, fname)

CLUSTERING_PARAMS <- list(
  k_min = 2,
  k_max = 8,
  description = "Intervention-focused resilience profiles"
)

DISTANCE_METHODS <- list(
  gower_baseline = list(
    name = "Gower_Baseline",
    params = list(
      preset = "gower",
      commensurable = TRUE
    ),
    description = "Commensurable Gower distance for mixed-type data"
  )
)

cat("✓ Parameters configured\n")
```

# Data Loading and Processing

```{r data-loading}
cat("\n=== LOADING AND PROCESSING PISA DATA ===\n")

pisa_cycles <- list(
  "2015" = list(
    student_file = "CY6_MS_CMB_STU_QQQ.sav",
    school_file = "CY6_MS_CMB_SCH_QQQ.sav",
    year = 2015,
    var_mapping = list(
      student = c("PARED" = "PAREDINT"),
      school = c()
    )
  ),
  "2018" = list(
    student_file = "CY07_MSU_STU_QQQ.sav",
    school_file = "CY07_MSU_SCH_QQQ.sav",
    year = 2018,
    var_mapping = list(
      student = c(),
      school = c("STRATIO" = "STRATIO")
    )
  ),
  "2022" = list(
    student_file = "CY08MSP_STU_QQQ.sav",
    school_file = "CY08MSP_SCH_QQQ.sav",
    year = 2022,
    var_mapping = list(
      student = c(),
      school = c("SMRATIO" = "STRATIO")
    )
  )
)

core_student_vars <- list(
  ids = c("CNT", "CNTSCHID", "CNTSTUID"),
  achievement = c(paste0("PV", 1:10, "MATH"), paste0("PV", 1:10, "READ"), paste0("PV", 1:10, "SCIE")),
  ses_background = c("ESCS", "HISEI", "PAREDINT", "HOMEPOS", "ICTRES", "WEALTH", "HEDRES"),
  demographics = c("ST004D01T", "IMMIG", "REPEAT"),
  motivation = c("ANXTEST", "ANXMAT", "MATHMOT", "MATHEFF", "SCIEEFF", "JOYSCIE", "INTMAT",
                "GRWTHMND", "WORKMAST", "RESILIENCE", "MASTGOAL", "GFOFAIL", "COMPETE",
                "INTBRSCI", "INSTSCIE", "MOTIVAT", "JOYREAD", "INTSCIE", "READINTEREST", 
                "GLOBCOMP", "MATHINT"),
  psychological = c("PERSEV", "GRWTHMND", "WORKMAST", "COMPETE", "EPIST", "BELONG",
                   "GCSELFEFF", "METASUM", "ADAPTIVITY", "PERSEVAGR", "GROSAGR", "STRESAGR",
                   "EMOCOAGR", "CURIOAGR", "COOPAGR", "EMPATAGR", "ASSERAGR", "CREATBELIEF"),
  social_support = c("TEACHSUP", "BELONG", "FAMSUP", "RELATST", "EMOSUPS", "EMOSUPSCI", "CURSUPP",
                    "BEINGBULLIED", "BULLIED", "FEELSAFE", "SOCCON", "SOCONPA", "PARINVOL", "PARSUPMAT"),
  learning_behavior = c("HOMWRK", "PERFEED", "DISCLIMA", "DISCLISCI", "DIRINS", "TRUANCY",
                       "ADINST", "SMINS", "SADDINST", "STUDYHMW", "EXERPRAC", "WORKPAY", "WORKHOME", 
                       "SKIPPING", "TARDYSD", "COGACRCO", "COGACMCO", "EXPOFA", "EXPO21ST", "DIGLRN"),
  wellbeing = c("LIFESAT", "ATTSCHL", "EUDMO", "SWBP", "PSYCHSYM", "SOCCON", "BODYIMA", "EXPWB"),
  perceptions = c("PQSCHOOL", "PISADIFF", "PQMOTHER", "PQFATHER"),
  weights = "W_FSTUWT"
)

core_school_vars <- list(
  ids = c("CNT", "CNTSCHID"),
  characteristics = c("SCHSIZE", "SCHLTYPE", "STRATIO", "PRIVATESCH", "CLSIZE", "SMRATIO"),
  resources = c("STAFFSHORT", "EDUSHORT", "RATCMP1", "SCIERES", "CREACTIV", "RATCMP2", "TOTAT", "DIGRES"),
  teacher_quality = c("PROAT5AB", "PROAT5AM", "PROAT6", "PROATCE", "PROSTAT", "PROSTCE",
                     "PROSTMAS", "PROATCE", "SCIENTCH", "TCHENTHU", "TCHFAIR", "TCHSUPPORT"),
  school_climate = c("STUBEHA", "TEACHBEHA", "SCMCEG", "NEGSCLIM", "ENCOURPG"),
  leadership = c("LEAD", "LEADPD", "LEADTCH", "SCHAUT", "TEACHPART", "INSTLEAD", "TCHPART"),
  community = c("SCHCOM", "PARPART", "EQUITYPOL", "CREATCURR"),
  weights = "W_SCHGRNRABWT"
)

safe_read_pisa_with_mapping <- function(file_path, var_list, var_mapping,
                                       target_countries = NULL, cycle_year = NULL) {
  tryCatch({
    file_vars <- names(read_sav(file_path, n_max = 0))
    reverse_mapping <- setNames(names(var_mapping), var_mapping)
    vars_to_request <- sapply(var_list, function(std_var) {
      if (std_var %in% names(reverse_mapping)) {
        reverse_mapping[[std_var]]
      } else {
        std_var
      }
    })
    available_vars <- intersect(vars_to_request, file_vars)
    if (length(available_vars) == 0) {
      warning("No requested variables found in ", basename(file_path))
      return(NULL)
    }
    dt <- setDT(read_sav(file_path, col_select = any_of(available_vars)))
    for (old_name in names(var_mapping)) {
      new_name <- var_mapping[[old_name]]
      if (old_name %in% names(dt)) {
        setnames(dt, old_name, new_name)
      }
    }
    if (!is.null(target_countries) && "CNT" %in% names(dt)) {
      dt <- dt[CNT %in% target_countries]
    }
    if (!is.null(cycle_year)) {
      dt[, CYCLE := cycle_year]
    }
    return(dt)
  }, error = function(e) {
    cat("Error reading file:", basename(file_path), "\n")
    return(NULL)
  })
}

# Load all PISA cycles
all_student_vars <- unlist(core_student_vars, use.names = FALSE)
all_school_vars <- unlist(core_school_vars, use.names = FALSE)
pisa_student_list <- list()
pisa_school_list <- list()

for (cycle_name in names(pisa_cycles)) {
  cycle_info <- pisa_cycles[[cycle_name]]
  cat("\n--- Loading PISA", cycle_info$year, "---\n")
  
  stu_file <- path(cycle_info$student_file)
  stu_data <- safe_read_pisa_with_mapping(
    stu_file, all_student_vars, cycle_info$var_mapping$student,
    target_countries, cycle_info$year
  )
  if (!is.null(stu_data) && nrow(stu_data) > 0) {
    pisa_student_list[[cycle_name]] <- stu_data
    cat("Students loaded:", nrow(stu_data), "\n")
  }
  
  sch_file <- path(cycle_info$school_file)
  sch_data <- safe_read_pisa_with_mapping(
    sch_file, all_school_vars, cycle_info$var_mapping$school,
    target_countries, cycle_info$year
  )
  if (!is.null(sch_data) && nrow(sch_data) > 0) {
    pisa_school_list[[cycle_name]] <- sch_data
    cat("Schools loaded:", nrow(sch_data), "\n")
  }
}

# Process PISA cycles with pooled plausible values
process_pisa_cycle <- function(stu_dt, sch_dt, cycle_year) {
  if (is.null(stu_dt) || nrow(stu_dt) == 0) return(NULL)
  
  cat(sprintf("Missing MATH plausible values: %d, READ: %d, SCIE: %d\n",
              sum(rowSums(is.na(stu_dt[, .SD, .SDcols = patterns("MATH$")])) > 0),
              sum(rowSums(is.na(stu_dt[, .SD, .SDcols = patterns("READ$")])) > 0),
              sum(rowSums(is.na(stu_dt[, .SD, .SDcols = patterns("SCIE$")])) > 0)))
  
  escs_q25_national <- quantile(stu_dt$ESCS, probs = 0.25, na.rm = TRUE)
  
  orig_resilient_props <- numeric(10)
  oecd_resilient_props <- numeric(10)
  revised_resilient_props <- numeric(10)
  orig_resilient_vars <- numeric(10)
  oecd_resilient_vars <- numeric(10)
  revised_resilient_vars <- numeric(10)
  
  for (pv in 1:10) {
    stu_dt[, `:=`(
      MATH_PV = get(paste0("PV", pv, "MATH")),
      READ_PV = get(paste0("PV", pv, "READ")),
      SCIE_PV = get(paste0("PV", pv, "SCIE"))
    )]
    stu_dt[, ACHIEVEMENT_PV := (MATH_PV + READ_PV + SCIE_PV) / 3]
    
    math_q75_national <- quantile(stu_dt$MATH_PV, probs = 0.75, na.rm = TRUE)
    read_q75_national <- quantile(stu_dt$READ_PV, probs = 0.75, na.rm = TRUE)
    scie_q75_national <- quantile(stu_dt$SCIE_PV, probs = 0.75, na.rm = TRUE)
    achievement_q75_national <- quantile(stu_dt$ACHIEVEMENT_PV, probs = 0.75, na.rm = TRUE)
    
    stu_dt[, DISADVANTAGED := ifelse(ESCS <= escs_q25_national, 1, 0)]
    
    MATH_LEVEL3_THRESHOLD <- 482.38
    READ_LEVEL3_THRESHOLD <- 480.18
    SCIE_LEVEL3_THRESHOLD <- 484.14
    
    stu_dt[DISADVANTAGED == 1, paste0("RESILIENT_PV", pv) := ifelse(
      MATH_PV >= MATH_LEVEL3_THRESHOLD &
        READ_PV >= READ_LEVEL3_THRESHOLD &
        SCIE_PV >= SCIE_LEVEL3_THRESHOLD, 1, 0)]
    
    stu_dt[DISADVANTAGED == 1, paste0("RESILIENT_OECD_PV", pv) := ifelse(
      ESCS <= escs_q25_national & ACHIEVEMENT_PV >= achievement_q75_national, 1, 0)]
    
    stu_dt[DISADVANTAGED == 1, paste0("RESILIENT_REVISED_PV", pv) := ifelse(
      ESCS <= escs_q25_national &
        MATH_PV >= math_q75_national &
        READ_PV >= read_q75_national &
        SCIE_PV >= scie_q75_national, 1, 0)]
    
    disadv_dt <- stu_dt[DISADVANTAGED == 1]
    n_disadv <- nrow(disadv_dt)
    orig_resilient_props[pv] <- mean(disadv_dt[[paste0("RESILIENT_PV", pv)]], na.rm = TRUE)
    oecd_resilient_props[pv] <- mean(disadv_dt[[paste0("RESILIENT_OECD_PV", pv)]], na.rm = TRUE)
    revised_resilient_props[pv] <- mean(disadv_dt[[paste0("RESILIENT_REVISED_PV", pv)]], na.rm = TRUE)
    orig_resilient_vars[pv] <- var(disadv_dt[[paste0("RESILIENT_PV", pv)]], na.rm = TRUE) / n_disadv
    oecd_resilient_vars[pv] <- var(disadv_dt[[paste0("RESILIENT_OECD_PV", pv)]], na.rm = TRUE) / n_disadv
    revised_resilient_vars[pv] <- var(disadv_dt[[paste0("RESILIENT_REVISED_PV", pv)]], na.rm = TRUE) / n_disadv
  }
  
  pool_proportions <- function(props, vars, n) {
    mean_prop <- mean(props)
    within_var <- mean(vars)
    between_var <- var(props)
    total_var <- within_var + between_var * (1 + 1/length(props))
    se <- sqrt(total_var)
    return(list(mean = mean_prop, se = se, n = round(mean_prop * n)))
  }
  
  n_disadv <- nrow(stu_dt[DISADVANTAGED == 1])
  orig_pooled <- pool_proportions(orig_resilient_props, orig_resilient_vars, n_disadv)
  oecd_pooled <- pool_proportions(oecd_resilient_props, oecd_resilient_vars, n_disadv)
  revised_pooled <- pool_proportions(revised_resilient_props, revised_resilient_vars, n_disadv)
  
  stu_dt[, RESILIENT := orig_pooled$mean]
  stu_dt[, RESILIENT_OECD := oecd_pooled$mean]
  stu_dt[, RESILIENT_REVISED := revised_pooled$mean]
  stu_dt[, RESILIENT_SE := orig_pooled$se]
  stu_dt[, RESILIENT_OECD_SE := oecd_pooled$se]
  stu_dt[, RESILIENT_REVISED_SE := revised_pooled$se]
  
  disadv_dt <- stu_dt[DISADVANTAGED == 1]
  
  if (!is.null(sch_dt) && nrow(sch_dt) > 0) {
    disadv_dt <- merge(disadv_dt, sch_dt, by = c("CNT", "CNTSCHID", "CYCLE"), all.x = TRUE)
  }
  
  disadv_dt[, grep("PV[0-9]+", names(disadv_dt)) := NULL]
  disadv_dt[, c("MATH_PV", "READ_PV", "SCIE_PV", "ACHIEVEMENT_PV") := NULL]
  
  return(list(
    data = disadv_dt,
    pooled_results = list(
      original = orig_pooled,
      oecd = oecd_pooled,
      revised = revised_pooled
    )
  ))
}

processed_list <- mapply(process_pisa_cycle, pisa_student_list, pisa_school_list,
                         names(pisa_cycles), SIMPLIFY = FALSE)

cat("\n=== DATA PROCESSING SUMMARY ===\n")
for (cycle in names(processed_list)) {
  if (!is.null(processed_list[[cycle]])) {
    n_students <- nrow(processed_list[[cycle]]$data)
    n_resilient_orig <- processed_list[[cycle]]$pooled_results$original$n
    resilience_rate_orig <- round(processed_list[[cycle]]$pooled_results$original$mean * 100, 1)
    resilience_se_orig <- round(processed_list[[cycle]]$pooled_results$original$se * 100, 2)
    
    cat("- Cycle", cycle, ":", n_students, "disadvantaged students\n")
    cat(" Original resilient (Level 3 all domains):", n_resilient_orig, "(", resilience_rate_orig, "%, SE:", resilience_se_orig, "%)\n")
  }
}
```

# Variable Selection

```{r variable_selection}
cat("\n=== CYCLE-SPECIFIC VARIABLE SELECTION ===\n")

identify_cycle_variables <- function(cycle_data, cycle_name) {
  cat(sprintf("→ Variable analysis for PISA %s...\n", cycle_name))
  
  all_available_vars <- names(cycle_data)
  cat(sprintf("  Total variables in dataset: %d\n", length(all_available_vars)))
  
  core_vars <- c("ESCS")
  
  extended_vars <- list(
    ses_background = c("HISEI", "PAREDINT", "HOMEPOS", "ICTRES", "CULTPOSS", "HEDRES", "WEALTH"),
    demographics = c("ST004D01T", "IMMIG", "REPEAT"),
    motivation = c("ANXMAT", "MATHEFF", "MATHMOT", "MATHPERS", "SCIEEFF", "JOYSCIE", "INTMAT", 
                   "JOYREAD", "WORKMAST", "MASTGOAL", "GFOFAIL", "COMPETE", "MOTIVAT"),
    psychological = c("PERSEVAGR", "GROSAGR", "STRESAGR", "EMOCOAGR", "CURIOAGR", "COOPAGR", 
                      "EMPATAGR", "ASSERAGR", "METASUM", "ADAPTIVITY", "GCSELFEFF"),
    learning_behavior = c("STUDYHMW", "EXERPRAC", "WORKPAY", "WORKHOME", "SKIPPING", "TARDYSD",
                          "PERFEED", "DISCLIM", "COGACRCO", "COGACMCO", "EXPOFA", "EXPO21ST"),
    social_emotional = c("BELONG", "TEACHSUP", "RELATST", "BULLIED", "FEELSAFE", "FAMSUP", 
                         "CURSUPP", "EMOSUPS", "SOCCON", "SOCONPA"),
    wellbeing = c("LIFESAT", "PSYCHSYM", "BODYIMA", "ATTSCHL", "EXPWB"),
    school_context = c("CLSIZE", "MCLSIZE", "SCHSIZE", "SCHLTYPE", "STRATIO", "SMRATIO", 
                       "STAFFSHORT", "EDUSHORT", "RATCMP1", "RATCMP2", "STUBEHA", "TEACHBEHA",
                       "NEGSCLIM", "ENCOURPG", "EDULEAD", "INSTLEAD", "TCHPART"),
    ict_usage = c("ICTRES", "ICTSCH", "ICTHOME", "ICTQUAL", "ICTSUBJ", "ICTENQ", 
                  "ICTFEED", "ICTOUT", "ICTEFFIC"),
    creativity = c("CREATEFF", "CREATSCH", "CREATFAM", "CREATAS", "CREATOOS", 
                   "CREATOP", "OPENART", "IMAGINE")
  )
  
  # Pattern-based discovery
  educational_patterns <- c("TEACH", "ANXI", "EFFI", "JOY", "MOTI", "BELON", "STRESS", "EMO")
  pattern_vars <- character(0)
  for (pattern in educational_patterns) {
    pattern_matches <- all_available_vars[grepl(pattern, all_available_vars, ignore.case = TRUE)]
    pattern_vars <- c(pattern_vars, pattern_matches)
  }
  
  exclude_patterns <- c("^PV[0-9]", "^W_", "^CNT", "ID$", "^ST00[0-9]Q[0-9][0-9]", "FLAG", "ADMIN", "^PA", "^TC")
  for (exclude_pattern in exclude_patterns) {
    pattern_vars <- pattern_vars[!grepl(exclude_pattern, pattern_vars)]
  }
  extended_vars$pattern_based <- unique(pattern_vars)
  
  available_vars <- all_available_vars
  cycle_clustering_vars <- core_vars[core_vars %in% available_vars]
  
  cat("\n  Category analysis:\n")
  for (category in names(extended_vars)) {
    category_vars <- extended_vars[[category]]
    available_category_vars <- intersect(category_vars, available_vars)
    
    coverage <- if(length(category_vars) > 0) length(available_category_vars) / length(category_vars) else 0
    include_category <- if(category == "pattern_based") {
      length(available_category_vars) > 0
    } else {
      coverage >= 0.2
    }
    
    if (include_category) {
      cycle_clustering_vars <- c(cycle_clustering_vars, available_category_vars)
      status_symbol <- "✓"
    } else {
      status_symbol <- "✗"
    }
    
    cat(sprintf("    %s %s: %d/%d variables (%.1f%% coverage)\n", 
                status_symbol, category, 
                length(available_category_vars), length(category_vars), coverage * 100))
  }
  
  cycle_clustering_vars <- unique(cycle_clustering_vars)
  
  # Remove achievement variables
  achievement_vars <- c("MATH_AVG", "READ_AVG", "SCIE_AVG", "ACHIEVEMENT_OVERALL")
  cycle_clustering_vars <- setdiff(cycle_clustering_vars, achievement_vars)
  
  suspected_categorical <- cycle_clustering_vars[grepl("^ST[0-9].*Q[0-9]|IMMIG|SCHLTYPE|REPEAT|SKIP|TARDY|GENDER|CNT$", cycle_clustering_vars)]
  suspected_continuous <- setdiff(cycle_clustering_vars, suspected_categorical)
  
  cat(sprintf("\n✓ Selected %d clustering variables for PISA %s\n", 
              length(cycle_clustering_vars), cycle_name))
  cat(sprintf("  → Estimated: %d continuous, %d categorical\n", 
              length(suspected_continuous), length(suspected_categorical)))
  
  return(list(
    clustering_vars = cycle_clustering_vars,
    suspected_categorical = suspected_categorical,
    suspected_continuous = suspected_continuous,
    total_available = length(all_available_vars)
  ))
}

cycle_variable_analysis <- list()
for (cycle_name in names(processed_list)) {
  if (!is.null(processed_list[[cycle_name]])) {
    cycle_variable_analysis[[cycle_name]] <- identify_cycle_variables(
      processed_list[[cycle_name]]$data, cycle_name
    )
  }
}

cat("\n✓ Variable analysis completed for all cycles\n")
```

# Data Preparation with MICE Imputation

```{r data-preparation}
cat("\n=== DATA PREPARATION WITH MICE IMPUTATION ===\n")

prepare_clustering_data <- function(cycle_data, cycle_vars, cycle_name,
                                   missing_threshold = 0.3, mice_iterations = 5) {
  
  cat(sprintf("→ Preparing clustering dataset for PISA %s...\n", cycle_name))
  
  essential_cols <- c("RESILIENT", "RESILIENT_OECD", "CYCLE", "W_FSTUWT", "CNTSTUID", "CNTSCHID")
  available_essential <- intersect(essential_cols, names(cycle_data))
  
  all_vars <- c(cycle_vars, available_essential)
  clustering_data <- cycle_data[, all_vars, with = FALSE]
  
  for (col_name in names(clustering_data)) {
    if (inherits(clustering_data[[col_name]], "haven_labelled")) {
      clustering_data[[col_name]] <- as.numeric(clustering_data[[col_name]])
    }
  }
  
  cat("→ Analyzing missing value patterns...\n")
  missing_summary <- clustering_data[, lapply(.SD, function(x) sum(is.na(x))/length(x)),
                                     .SDcols = cycle_vars]
  
  good_vars <- names(missing_summary)[missing_summary <= missing_threshold]
  dropped_vars <- names(missing_summary)[missing_summary > missing_threshold]
  
  good_vars_for_clustering <- setdiff(good_vars, c("RESILIENT", "RESILIENT_OECD"))
  
  if (length(dropped_vars) > 0) {
    cat(sprintf("⚠ Dropped %d variables with >%.0f%% missingness\n",
                length(dropped_vars), missing_threshold * 100))
  }
  
  final_vars <- c(good_vars_for_clustering, available_essential)
  clustering_data <- clustering_data[, final_vars, with = FALSE]
  
  categorical_vars <- c("ST004D01T", "IMMIG", "SCHLTYPE", "TESTLANG", "PARINVOL", "EQUITYPOL")
  pattern_categorical <- good_vars_for_clustering[grepl("^ST[0-9].*Q[0-9]|IMMIG|SCHLTYPE|TESTLANG|PARINVOL|EQUITYPOL|SKIP|TARDY|GENDER|CNT$", good_vars_for_clustering)]
  categorical_vars <- unique(c(categorical_vars, pattern_categorical))
  available_categorical <- intersect(categorical_vars, good_vars_for_clustering)
  
  continuous_vars <- setdiff(good_vars_for_clustering, available_categorical)
  
  for (var in available_categorical) {
    if (var %in% names(clustering_data)) {
      clustering_data[[var]] <- as.factor(clustering_data[[var]])
    }
  }
  
  remaining_missing <- clustering_data[, lapply(.SD, function(x) sum(is.na(x))),
                                       .SDcols = good_vars_for_clustering]
  total_missing <- sum(unlist(remaining_missing))
  
  if (total_missing > 0) {
    cat(sprintf("→ Performing MICE imputation for %d remaining missing values...\n", total_missing))
    
    mice_data <- clustering_data[, good_vars_for_clustering, with = FALSE]
    mice_df <- as.data.frame(mice_data)
    
    for (col_name in names(mice_df)) {
      if (inherits(mice_df[[col_name]], "haven_labelled")) {
        mice_df[[col_name]] <- as.numeric(mice_df[[col_name]])
      }
    }
    
    mice_methods <- mice::make.method(mice_df)
    for (var in names(mice_df)) {
      if (is.factor(mice_df[[var]])) {
        n_levels <- length(levels(mice_df[[var]]))
        if (n_levels == 2) {
          mice_methods[var] <- "logreg"
        } else if (n_levels > 2) {
          mice_methods[var] <- "polyreg"
        }
      } else if (is.numeric(mice_df[[var]])) {
        mice_methods[var] <- "pmm"
      }
    }
    
    tryCatch({
      mice_result <- mice::mice(mice_df, m = mice_iterations, method = mice_methods,
                                printFlag = FALSE, seed = 42)
      completed_mice_df <- mice::complete(mice_result, 1)
      completed_mice_dt <- setDT(completed_mice_df)
      
      for (var in good_vars_for_clustering) {
        clustering_data[[var]] <- completed_mice_dt[[var]]
      }
      
      cat("✓ MICE imputation completed successfully\n")
      
    }, error = function(e) {
      cat(sprintf("⚠ MICE imputation failed: %s\n", e$message))
    })
    
  } else {
    cat("✓ No missing values detected\n")
  }
  
  initial_n <- nrow(clustering_data)
  essential_missing <- rowSums(is.na(clustering_data[, available_essential, with = FALSE]))
  clustering_data <- clustering_data[essential_missing == 0]
  final_n <- nrow(clustering_data)
  
  if (final_n < initial_n) {
    cat(sprintf("⚠ Removed %d cases with missing essential variables\n", initial_n - final_n))
  }
  
  cat(sprintf("✓ Final dataset: %d students (%.1f%% retention)\n",
              final_n, 100 * final_n / initial_n))
  
  if ("W_FSTUWT" %in% names(clustering_data)) {
    clustering_data[, W_FSTUWT_NORM := W_FSTUWT / mean(W_FSTUWT, na.rm = TRUE)]
  }
  
  return(list(
    data = clustering_data,
    clustering_vars = good_vars_for_clustering,
    categorical_vars = available_categorical,
    continuous_vars = continuous_vars,
    cycle = cycle_name,
    essential_vars = available_essential
  ))
}

cycle_prepared_data <- list()
for (cycle_name in names(processed_list)) {
  if (!is.null(processed_list[[cycle_name]])) {
    cycle_vars <- cycle_variable_analysis[[cycle_name]]$clustering_vars
    cycle_prepared_data[[cycle_name]] <- prepare_clustering_data(
      processed_list[[cycle_name]]$data, cycle_vars, cycle_name
    )
  }
}

cat("\n✓ Data preparation completed for all cycles\n")
```

# Enhanced Clustering with Variable Weighting

```{r clustering-setup}
cat("\n=== CLUSTERING WITH VARIABLE WEIGHTING ===\n")

clean_clustering_data <- function(data) {
  clean_data <- as.data.frame(data)
  for (col_name in names(clean_data)) {
    if (inherits(clean_data[[col_name]], "haven_labelled")) {
      clean_data[[col_name]] <- as.numeric(clean_data[[col_name]])
    }
  }
  list_cols <- sapply(clean_data, is.list)
  if (any(list_cols)) {
    clean_data <- clean_data[, !list_cols, drop = FALSE]
  }
  return(clean_data)
}

calculate_weighted_distances <- function(cycle_data_info, distance_methods) {
  cycle_name <- cycle_data_info$cycle
  cat(sprintf("\n→ Calculating weighted distance matrices for PISA %s...\n", cycle_name))
  
  distance_data <- cycle_data_info$data[, cycle_data_info$clustering_vars, with = FALSE]
  distance_data <- na.omit(distance_data)
  distance_data <- clean_clustering_data(distance_data)
  
  for (var in names(distance_data)) {
    if (any(is.na(distance_data[[var]]))) {
      if (is.factor(distance_data[[var]])) {
        mode_val <- names(sort(table(distance_data[[var]]), decreasing = TRUE))[1]
        distance_data[is.na(distance_data[[var]]), var] <- mode_val
      } else {
        median_val <- median(distance_data[[var]], na.rm = TRUE)
        distance_data[is.na(distance_data[[var]]), var] <- median_val
      }
    }
  }
  
  continuous_vars <- cycle_data_info$continuous_vars
  
  # Apply intervention-focused variable weighting
  var_weights <- rep(1.0, ncol(distance_data))
  names(var_weights) <- names(distance_data)
  
  # SES variables (2.0x - important context)
  ses_vars <- intersect(names(distance_data), c("ESCS", "HISEI", "HOMEPOS", "ICTRES", "WEALTH"))
  var_weights[ses_vars] <- 2.0
  
  # Psychological/Motivation variables (1.8x - key for intervention)
  psych_vars <- intersect(names(distance_data),
                          c("ANXMAT", "ANXTEST", "MATHEFF", "SCIEEFF", "JOYSCIE",
                            "WORKMAST", "MASTGOAL", "GFOFAIL", "COMPETE", "METASUM", 
                            "ADAPTIVITY", "GCSELFEFF", "PERSEVAGR", "GROSAGR", "STRESAGR"))
  var_weights[psych_vars] <- 1.8
  
  # Social Support variables (1.5x)
  social_vars <- intersect(names(distance_data), c("BELONG", "EMOSUPS", "FAMSUP", "RELATST", "TEACHSUP"))
  var_weights[social_vars] <- 1.5
  
  # Learning Behavior variables (1.3x)
  learning_vars <- intersect(names(distance_data),
                             c("STUDYHMW", "EXERPRAC", "WORKPAY", "WORKHOME", "SKIPPING", "PERFEED", "DISCLIM"))
  var_weights[learning_vars] <- 1.3
  
  cat(sprintf("  Variable weights applied:\n"))
  cat(sprintf("    SES (2.0x): %d variables\n", length(ses_vars)))
  cat(sprintf("    Psychology/Motivation (1.8x): %d variables\n", length(psych_vars)))
  cat(sprintf("    Social Support (1.5x): %d variables\n", length(social_vars)))
  cat(sprintf("    Learning Behavior (1.3x): %d variables\n", length(learning_vars)))
  
  # Apply weights to continuous variables
  weighted_distance_data <- distance_data
  for (var in names(var_weights)) {
    if (var %in% continuous_vars) {
      weighted_distance_data[[var]] <- distance_data[[var]] * var_weights[var]
    }
  }
  
  cycle_distance_matrices <- list()
  
  for (method_name in names(distance_methods)) {
    cat(sprintf("  → %s with variable weighting...", method_name))
    method_info <- distance_methods[[method_name]]
    
    tryCatch({
      start_time <- Sys.time()
      dist_matrix <- do.call(manydist::mdist, c(list(x = weighted_distance_data), method_info$params))
      end_time <- Sys.time()
      calc_time <- as.numeric(difftime(end_time, start_time, units = "secs"))
      
      cycle_distance_matrices[[method_name]] <- list(
        matrix = as.dist(dist_matrix),
        method = method_info$name,
        description = method_info$description,
        calculation_time = calc_time,
        variable_weights = var_weights,
        weighted_data_applied = TRUE
      )
      
      cat(sprintf(" %.1fs ✓\n", calc_time))
      
    }, error = function(e) {
      cat(sprintf(" ERROR: %s\n", e$message))
    })
  }
  
  return(cycle_distance_matrices)
}

extract_clusterbenchstats_matrix <- function(cbs_result) {
  if (is.null(cbs_result$sstat)) {
    if (!is.null(cbs_result$qstat)) {
      stat_component <- cbs_result$qstat
    } else {
      return(NULL)
    }
  } else {
    stat_component <- cbs_result$sstat
  }
  
  methods <- stat_component$name
  if (is.null(methods)) methods <- stat_component$method
  G_values <- stat_component$minG:stat_component$maxG
  stat_names <- stat_component$statistics
  
  n_rows <- length(methods) * length(G_values)
  stats_matrix <- matrix(NA, nrow = n_rows, ncol = length(stat_names))
  colnames(stats_matrix) <- stat_names
  
  row_names <- character(n_rows)
  row_idx <- 1
  for (method in methods) {
    for (g in G_values) {
      row_names[row_idx] <- paste(method, g, sep = ".")
      row_idx <- row_idx + 1
    }
  }
  rownames(stats_matrix) <- row_names
  
  row_idx <- 1
  for (method_idx in 1:length(methods)) {
    method_data <- stat_component[[method_idx]]
    
    if (is.list(method_data)) {
      for (g_idx in seq_along(G_values)) {
        pos <- g_idx + 1
        if (pos <= length(method_data)) {
          g_stats <- method_data[[pos]]
          
          if (is.list(g_stats) && length(g_stats) > 0) {
            for (col_idx in 1:length(stat_names)) {
              stat_name <- stat_names[col_idx]
              if (stat_name %in% names(g_stats)) {
                value <- g_stats[[stat_name]]
                if (is.numeric(value) && length(value) == 1 && is.finite(value)) {
                  stats_matrix[row_idx, col_idx] <- value
                }
              }
            }
          }
        }
        row_idx <- row_idx + 1
      }
    } else {
      row_idx <- row_idx + length(G_values)
    }
  }
  
  valid_rows <- rowSums(!is.na(stats_matrix)) > 0
  if (sum(valid_rows) == 0) {
    return(NULL)
  }
  
  return(stats_matrix[valid_rows, , drop = FALSE])
}

develop_intervention_weights <- function(available_stats) {
  weights <- rep(0.1, length(available_stats))
  names(weights) <- available_stats
  
  for (stat in available_stats) {
    if (grepl("boot|stability", stat, ignore.case = TRUE)) {
      weights[stat] <- 1.8
    } else if (grepl("silhouette|silwidth|asw", stat, ignore.case = TRUE)) {
      weights[stat] <- 1.5
    } else if (grepl("within|ave\\.within|homogeneity", stat, ignore.case = TRUE)) {
      weights[stat] <- 1.3
    } else if (grepl("entropy", stat, ignore.case = TRUE)) {
      weights[stat] <- 1.4
    } else if (grepl("separation|sep|minsep", stat, ignore.case = TRUE)) {
      weights[stat] <- 0.4
    }
  }
  
  return(weights)
}

find_best_intervention_solution <- function(distance_results, cycle_name) {
  best_solution <- NULL
  best_score <- -Inf
  
  for (dist_method in names(distance_results)) {
    dist_result <- distance_results[[dist_method]]
    
    if (!is.null(dist_result$analysis)) {
      stats_matrix <- extract_clusterbenchstats_matrix(dist_result$analysis)
      
      if (!is.null(stats_matrix)) {
        available_cols <- colnames(stats_matrix)
        weights <- develop_intervention_weights(available_cols)
        
        scores <- apply(stats_matrix, 1, function(row) {
          if (all(is.na(row))) return(NA)
          row[!is.finite(row)] <- NA
          valid_indices <- !is.na(row)
          if (sum(valid_indices) == 0) return(NA)
          weighted.mean(row[valid_indices], weights[valid_indices], na.rm = TRUE)
        })
        
        if (length(scores) > 0 && any(!is.na(scores))) {
          best_idx <- which.max(scores)
          best_name <- names(scores)[best_idx]
          
          if (scores[best_idx] > best_score) {
            best_score <- scores[best_idx]
            best_solution <- list(
              solution = best_name,
              score = best_score,
              distance_method = dist_method,
              statistics = stats_matrix[best_name, ]
            )
          }
        }
      }
    }
  }
  
  return(best_solution)
}

perform_intervention_clustering <- function(cycle_data, cycle_name) {
  cat(sprintf("\n=== INTERVENTION CLUSTERING: PISA %s ===\n", cycle_name))
  
  cycle_distances <- calculate_weighted_distances(cycle_data, DISTANCE_METHODS)
  
  clustermethod <- c("disthclustCBI", "pamkCBI")
  clustermethodpars <- list()
  clustermethodpars[[1]] <- list(method = "ward.D2")
  clustermethodpars[[2]] <- list(diss = TRUE, usepam = TRUE)
  methodnames <- c("Ward", "PAM")
  distmethod <- rep(TRUE, length(clustermethod))
  ncinput <- rep(TRUE, length(clustermethod))
  
  G_range <- CLUSTERING_PARAMS$k_min:CLUSTERING_PARAMS$k_max
  
  distance_results <- list()
  
  for (dist_method_name in names(cycle_distances)) {
    if (!is.null(cycle_distances[[dist_method_name]])) {
      cat(sprintf("\n→ Running clustering with %s distance...\n", dist_method_name))
      
      dist_matrix <- cycle_distances[[dist_method_name]]$matrix
      
      tryCatch({
        cbs_result <- clusterbenchstats(
          data = dist_matrix,
          G = G_range,
          diss = TRUE,
          scaling = FALSE,
          clustermethod = clustermethod,
          methodnames = methodnames,
          distmethod = distmethod,
          ncinput = ncinput,
          clustermethodpars = clustermethodpars,
          npstats = FALSE,
          useboot = FALSE,
          trace = FALSE,
          useallmethods = FALSE,
          useallg = FALSE,
          nnruns = 15,
          kmruns = 0,
          fnruns = 15,
          avenruns = 15
        )
        
        distance_results[[dist_method_name]] <- list(
          analysis = cbs_result,
          distance_info = cycle_distances[[dist_method_name]]
        )
        
        cat("    ✓ Clustering completed\n")
        
      }, error = function(e) {
        cat(sprintf("    ✗ Clustering failed: %s\n", e$message))
      })
    }
  }
  
  best_solution <- find_best_intervention_solution(distance_results, cycle_name)
  
  return(list(
    cycle = cycle_name,
    distance_results = distance_results,
    optimal_solution = best_solution,
    distance_matrices = cycle_distances
  ))
}

enhanced_clustering_results <- list()
for (cycle_name in names(cycle_prepared_data)) {
  if (!is.null(cycle_prepared_data[[cycle_name]])) {
    enhanced_clustering_results[[cycle_name]] <- perform_intervention_clustering(
      cycle_prepared_data[[cycle_name]], cycle_name
    )
  }
}

cat("\n✓ Enhanced clustering analysis completed for all cycles\n")
```

# Stability-First Clustering (Akhanli & Hennig Framework)

```{r stability-first-clustering}
cat("\n=== AKHANLI & HENNIG STABILITY-FIRST CLUSTERING ===\n")

filter_solutions_by_size <- function(clusters, min_size_pct = 0.05) {
  cluster_sizes <- table(clusters)
  min_n <- ceiling(min_size_pct * length(clusters))
  return(!any(cluster_sizes < min_n))
}

perform_stability_first_clustering <- function(enhanced_clustering_results, stability_threshold = 0.4) {
  stability_calibrated_results <- list()
  
  for (cycle_name in names(enhanced_clustering_results)) {
    result <- enhanced_clustering_results[[cycle_name]]
    
    cat(sprintf("\n→ Applying stability-first approach for PISA %s...\n", cycle_name))
    
    stable_solution <- NULL
    candidate_solutions <- list()
    
    for (dist_method in names(result$distance_results)) {
      if (!is.null(result$distance_results[[dist_method]])) {
        
        analysis <- result$distance_results[[dist_method]]$analysis
        
        if (!is.null(analysis)) {
          stats_matrix <- extract_clusterbenchstats_matrix(analysis)
          
          if (!is.null(stats_matrix)) {
            dist_matrix <- result$distance_matrices[[dist_method]]$matrix
            
            for (solution_name in rownames(stats_matrix)) {
              parts <- strsplit(solution_name, "\\.")[[1]]
              method_name <- parts[1]
              k_value <- as.numeric(parts[2])
              
              # Generate clusters
              if (method_name == "Ward") {
                hc <- hclust(dist_matrix, method = "ward.D2")
                test_clusters <- cutree(hc, k = k_value)
              } else if (method_name == "PAM") {
                test_clusters <- pam(dist_matrix, k = k_value, diss = TRUE)$clustering
              }
              
              if (!filter_solutions_by_size(test_clusters, min_size_pct = 0.05)) {
                next
              }
              
              cat(sprintf("    Testing stability: %s %s k=%d...", dist_method, method_name, k_value))
              
              tryCatch({
                if (method_name == "Ward") {
                  boot_result <- clusterboot(
                    data = dist_matrix, B = 25, distances = TRUE, bootmethod = "boot",
                    clustermethod = disthclustCBI, k = k_value, method = "ward.D2",
                    cut = "number", seed = 42, count = FALSE
                  )
                } else if (method_name == "PAM") {
                  boot_result <- clusterboot(
                    data = dist_matrix, B = 25, distances = TRUE, bootmethod = "boot",
                    clustermethod = pamkCBI, k = k_value, usepam = TRUE, seed = 42, count = FALSE
                  )
                }
                
                overall_stability <- mean(boot_result$bootmean, na.rm = TRUE)
                cat(sprintf(" Stability: %.3f\n", overall_stability))
                
                if (overall_stability >= stability_threshold) {
                  original_scores <- stats_matrix[solution_name, ]
                  
                  candidate_solutions[[paste(dist_method, solution_name, sep = "_")]] <- list(
                    solution = solution_name,
                    method = method_name,
                    k = k_value,
                    distance_method = dist_method,
                    stability_score = overall_stability,
                    original_scores = original_scores,
                    passed_stability = TRUE
                  )
                  
                  cat(sprintf("      ✓ PASSED stability filter (%.3f >= %.3f)\n", 
                              overall_stability, stability_threshold))
                } else {
                  cat(sprintf("      ✗ FAILED stability filter (%.3f < %.3f)\n", 
                              overall_stability, stability_threshold))
                }
                
              }, error = function(e) {
                cat(sprintf(" ERROR: %s\n", e$message))
              })
            }
          }
        }
      }
    }
    
    if (length(candidate_solutions) > 0) {
      cat(sprintf("  → Found %d stable solutions\n", length(candidate_solutions)))
      
      available_stats <- names(candidate_solutions[[1]]$original_scores)
      weights <- develop_intervention_weights(available_stats)
      
      best_score <- -Inf
      best_solution <- NULL
      
      for (sol_name in names(candidate_solutions)) {
        sol <- candidate_solutions[[sol_name]]
        
        valid_indices <- !is.na(sol$original_scores) & is.finite(sol$original_scores)
        
        if (sum(valid_indices) > 0) {
          composite_score <- weighted.mean(
            sol$original_scores[valid_indices], 
            weights[valid_indices], 
            na.rm = TRUE
          )
          
          stability_weighted_score <- composite_score * sol$stability_score
          
          if (stability_weighted_score > best_score) {
            best_score <- stability_weighted_score
            best_solution <- sol
            best_solution$composite_score <- composite_score
            best_solution$final_score <- stability_weighted_score
          }
        }
      }
      
      stable_solution <- best_solution
      
      if (!is.null(best_solution)) {
        cat(sprintf("  ✓ Best stable solution: %s %s k=%d (Stability: %.3f, Final: %.3f)\n",
                    best_solution$distance_method, best_solution$method, 
                    best_solution$k, best_solution$stability_score, best_solution$final_score))
      }
      
    } else {
      cat(sprintf("  ✗ No stable solutions found (threshold %.3f)\n", stability_threshold))
    }
    
    stability_calibrated_results[[cycle_name]] <- list(
      original_results = result,
      stable_solution = stable_solution,
      stability_threshold = stability_threshold
    )
  }
  
  return(stability_calibrated_results)
}

extract_final_clusters <- function(stability_calibrated_results) {
  all_cluster_assignments <- list()
  
  for (cycle_name in names(stability_calibrated_results)) {
    result <- stability_calibrated_results[[cycle_name]]
    
    cat(sprintf("\n→ Extracting clusters for PISA %s...\n", cycle_name))
    
    stable_sol <- result$stable_solution
    
    if (!is.null(stable_sol)) {
      original_result <- result$original_results
      dist_method <- stable_sol$distance_method
      dist_matrix <- original_result$distance_matrices[[dist_method]]$matrix
      
      method_name <- stable_sol$method
      k_value <- stable_sol$k
      
      if (method_name == "Ward") {
        hc <- hclust(dist_matrix, method = "ward.D2")
        clusters <- cutree(hc, k = k_value)
      } else if (method_name == "PAM") {
        clusters <- pam(dist_matrix, k = k_value, diss = TRUE)$clustering
      }
      
      all_cluster_assignments[[cycle_name]] <- list(
        clusters = clusters,
        method = method_name,
        k = k_value,
        distance_method = dist_method,
        cycle = cycle_name,
        solution_name = stable_sol$solution,
        stability_score = stable_sol$stability_score,
        final_score = stable_sol$final_score
      )
      
      cluster_dist <- table(clusters)
      cat(sprintf("  ✓ %s k=%d (Stability: %.3f, Final: %.3f)\n",
                  stable_sol$solution, k_value, stable_sol$stability_score, stable_sol$final_score))
      cat(sprintf("    Distribution: %s\n", paste(names(cluster_dist), "=", cluster_dist, collapse = ", ")))
      
    } else {
      cat(sprintf("  ✗ No stable solution found\n"))
    }
  }
  
  return(all_cluster_assignments)
}

add_clusters_to_data <- function(cycle_prepared_data, cluster_assignments) {
  enhanced_data_with_clusters <- list()
  
  for (cycle_name in names(cycle_prepared_data)) {
    if (!is.null(cycle_prepared_data[[cycle_name]])) {
      enhanced_data_with_clusters[[cycle_name]] <- cycle_prepared_data[[cycle_name]]
      
      if (cycle_name %in% names(cluster_assignments)) {
        cluster_info <- cluster_assignments[[cycle_name]]
        enhanced_data_with_clusters[[cycle_name]]$data$CLUSTER_INTERVENTION <- cluster_info$clusters
        enhanced_data_with_clusters[[cycle_name]]$cluster_info_intervention <- cluster_info
      }
    }
  }
  
  return(enhanced_data_with_clusters)
}

# Execute the Akhanli & Hennig approach
akhanli_hennig_results <- perform_stability_first_clustering(enhanced_clustering_results, stability_threshold = 0.4)
final_clusters <- extract_final_clusters(akhanli_hennig_results)
enhanced_data_with_clusters <- add_clusters_to_data(cycle_prepared_data, final_clusters)

cat("\n=== AKHANLI & HENNIG COMPLIANCE SUMMARY ===\n")
for (cycle_name in names(akhanli_hennig_results)) {
  result <- akhanli_hennig_results[[cycle_name]]
  sol <- result$stable_solution
  
  cat(sprintf("\nPISA %s:\n", cycle_name))
  if (!is.null(sol)) {
    cat(sprintf("  INTERVENTION: ✓ STABLE (%.3f) %s k=%d\n", 
                sol$stability_score, sol$solution, sol$k))
  } else {
    cat(sprintf("  INTERVENTION: ✗ NO STABLE SOLUTION FOUND\n"))
  }
}

cat("\n✓ Akhanli & Hennig (2020) framework implemented\n")
cat("✓ Stability filtering applied before composite scoring\n")
cat("✓ Variable weighting maintained within stable solutions\n")
```

# Cluster Analysis and Profiling

```{r cluster-profiling}
cat("\n=== COMPREHENSIVE CLUSTER ANALYSIS ===\n")

analyze_cluster_characteristics <- function(enhanced_data_with_clusters, cycle_name) {
  cat(sprintf("\n=== DETAILED CLUSTER ANALYSIS: PISA %s ===\n", cycle_name))
  
  cycle_data <- enhanced_data_with_clusters[[cycle_name]]
  
  if (!"CLUSTER_INTERVENTION" %in% names(cycle_data$data)) {
    cat("No intervention clusters found\n")
    return(NULL)
  }
  
  data_with_clusters <- cycle_data$data
  clusters <- data_with_clusters$CLUSTER_INTERVENTION
  n_clusters <- length(unique(clusters))
  
  key_continuous_vars <- intersect(c("ESCS", "HOMEPOS", "ICTRES", "WEALTH", "HEDRES", 
                                     "ANXMAT", "ANXTEST", "MATHEFF", "SCIEEFF", "JOYSCIE", 
                                     "TEACHSUP", "BELONG", "FAMSUP", "EMOSUPS"), 
                                   names(data_with_clusters))
  
  key_categorical_vars <- intersect(c("ST004D01T", "IMMIG", "SCHLTYPE", "REPEAT"), 
                                    names(data_with_clusters))
  
  resilience_vars <- intersect(c("RESILIENT", "RESILIENT_OECD"), names(data_with_clusters))
  
  cat(sprintf("Analyzing %d clusters with %d key continuous and %d categorical variables\n", 
              n_clusters, length(key_continuous_vars), length(key_categorical_vars)))
  
  # Overall statistics
  cat("\n=== OVERALL SAMPLE STATISTICS ===\n")
  for (var in key_continuous_vars[1:5]) {
    overall_mean <- mean(data_with_clusters[[var]], na.rm = TRUE)
    overall_sd <- sd(data_with_clusters[[var]], na.rm = TRUE)
    cat(sprintf("%s: M=%.2f, SD=%.2f\n", var, overall_mean, overall_sd))
  }
  
  # Analyze each cluster
  for (cluster_id in sort(unique(clusters))) {
    cluster_mask <- clusters == cluster_id
    cluster_data <- data_with_clusters[cluster_mask, ]
    n_students <- nrow(cluster_data)
    
    cat(sprintf("\n=== CLUSTER %d (n=%d, %.1f%%) ===\n", 
                cluster_id, n_students, 100*n_students/nrow(data_with_clusters)))
    
    # Resilience outcomes
    cat("RESILIENCE OUTCOMES:\n")
    for (var in resilience_vars) {
      count <- sum(cluster_data[[var]], na.rm = TRUE)
      rate <- mean(cluster_data[[var]], na.rm = TRUE) * 100
      overall_rate <- mean(data_with_clusters[[var]], na.rm = TRUE) * 100
      diff <- rate - overall_rate
      comparison <- if (abs(diff) < 1) "similar" else if (diff > 0) "higher" else "lower"
      
      cat(sprintf("  %s: %.0f/%d (%.1f%% vs %.1f%% overall - %s)\n", 
                  var, count, n_students, rate, overall_rate, comparison))
    }
    
    # Key characteristics
    cat("\nKEY CHARACTERISTICS:\n")
    for (var in key_continuous_vars) {
      cluster_mean <- mean(cluster_data[[var]], na.rm = TRUE)
      overall_mean <- mean(data_with_clusters[[var]], na.rm = TRUE)
      pooled_sd <- sd(data_with_clusters[[var]], na.rm = TRUE)
      effect_size <- (cluster_mean - overall_mean) / pooled_sd
      
      magnitude <- if (abs(effect_size) < 0.2) "small" 
                   else if (abs(effect_size) < 0.5) "medium" 
                   else "large"
      
      direction <- if (effect_size > 0.1) "above average" 
                   else if (effect_size < -0.1) "below average" 
                   else "average"
      
      cat(sprintf("  %s: M=%.2f (%.2f vs %.2f overall, %s %s difference)\n", 
                  var, cluster_mean, cluster_mean, overall_mean, magnitude, direction))
    }
    
    # Demographic profile
    if (length(key_categorical_vars) > 0) {
      cat("\nDEMOGRAPHIC PROFILE:\n")
      for (var in key_categorical_vars) {
        if (var %in% names(cluster_data)) {
          var_table <- table(cluster_data[[var]], useNA = "ifany")
          var_props <- prop.table(var_table) * 100
          
          mode_category <- names(var_table)[which.max(var_table)]
          mode_pct <- max(var_props)
          
        overall_table <- table(data_with_clusters[[var]], useNA = "ifany")
          overall_props <- prop.table(overall_table) * 100
          overall_mode_pct <- overall_props[mode_category]
          
          cat(sprintf("  %s: %s (%.1f%% vs %.1f%% overall)\n", 
                      var, mode_category, mode_pct, overall_mode_pct))
        }
      }
    }
  }
  
  # Cross-cluster comparison
  cat("\n=== VARIABLES THAT DIFFERENTIATE CLUSTERS MOST ===\n")
  
  effect_sizes <- sapply(key_continuous_vars, function(var) {
    cluster_means <- tapply(data_with_clusters[[var]], clusters, mean, na.rm = TRUE)
    overall_sd <- sd(data_with_clusters[[var]], na.rm = TRUE)
    max_diff <- max(cluster_means, na.rm = TRUE) - min(cluster_means, na.rm = TRUE)
    return(max_diff / overall_sd)
  })
  
  top_differentiating <- head(sort(effect_sizes, decreasing = TRUE), 5)
  
  for (var in names(top_differentiating)) {
    cluster_means <- tapply(data_with_clusters[[var]], clusters, mean, na.rm = TRUE)
    cat(sprintf("%s (effect size: %.2f):\n", var, top_differentiating[[var]]))
    for (i in 1:length(cluster_means)) {
      cat(sprintf("  Cluster %d: %.2f\n", i, cluster_means[i]))
    }
    cat("\n")
  }
  
  return(invisible(list(
    n_clusters = n_clusters,
    cluster_sizes = table(clusters),
    top_differentiating_vars = names(top_differentiating),
    effect_sizes = effect_sizes
  )))
}

# Run analysis for all cycles
cluster_analysis_results <- list()
for (cycle_name in names(enhanced_data_with_clusters)) {
  if ("CLUSTER_INTERVENTION" %in% names(enhanced_data_with_clusters[[cycle_name]]$data)) {
    cluster_analysis_results[[cycle_name]] <- analyze_cluster_characteristics(enhanced_data_with_clusters, cycle_name)
  }
}

# Additional summary function for intervention clusters
summarize_intervention_clusters <- function(enhanced_data_with_clusters) {
  cat("\n=== INTERVENTION CLUSTER SUMMARY ACROSS CYCLES ===\n")
  
  for (cycle_name in names(enhanced_data_with_clusters)) {
    if ("CLUSTER_INTERVENTION" %in% names(enhanced_data_with_clusters[[cycle_name]]$data)) {
      cycle_data <- enhanced_data_with_clusters[[cycle_name]]
      data_with_clusters <- cycle_data$data
      clusters <- data_with_clusters$CLUSTER_INTERVENTION
      
      cat(sprintf("\n--- PISA %s INTERVENTION SUMMARY ---\n", cycle_name))
      
      # Cluster sizes and resilience rates
      cluster_summary <- data.table(
        Cluster = sort(unique(clusters))
      )
      
      for (cluster_id in cluster_summary$Cluster) {
        cluster_mask <- clusters == cluster_id
        cluster_data <- data_with_clusters[cluster_mask, ]
        
        cluster_summary[Cluster == cluster_id, `:=`(
          N = nrow(cluster_data),
          Percentage = round(100 * nrow(cluster_data) / nrow(data_with_clusters), 1),
          Resilience_Rate = round(mean(cluster_data$RESILIENT, na.rm = TRUE) * 100, 1),
          ESCS_Mean = round(mean(cluster_data$ESCS, na.rm = TRUE), 2)
        )]
      }
      
      print(cluster_summary)
      
      # Key differentiating variables for this cycle
      if (cycle_name %in% names(cluster_analysis_results)) {
        top_vars <- cluster_analysis_results[[cycle_name]]$top_differentiating_vars[1:3]
        cat(sprintf("Top differentiating variables: %s\n", paste(top_vars, collapse = ", ")))
      }
    }
  }
}

summarize_intervention_clusters(enhanced_data_with_clusters)

cat("\n✓ Comprehensive cluster analysis completed\n")
```
          
# MDS Visualization

```{r mds-visualization, fig.width=10, fig.height=8}
# MDS Visualization
cat("\n=== MDS VISUALIZATION FOR CLUSTERING RESULTS ===\n")

create_mds_cluster_plot <- function(enhanced_data_with_clusters, cycle_name) {
  cat(sprintf("\n→ Creating MDS plot for PISA %s intervention clusters...\n", cycle_name))
  
  cycle_data <- enhanced_data_with_clusters[[cycle_name]]
  
  if (!"CLUSTER_INTERVENTION" %in% names(cycle_data$data)) {
    cat("No intervention clusters found\n")
    return(NULL)
  }
  
  data_with_clusters <- cycle_data$data
  clusters <- data_with_clusters$CLUSTER_INTERVENTION
  n_clusters <- length(unique(clusters))
  
  clustering_vars <- cycle_data$clustering_vars
  if (length(clustering_vars) == 0) {
    cat("No clustering variables available for MDS\n")
    return(NULL)
  }
  
  # Prepare data for MDS
  cluster_matrix <- data_with_clusters[, clustering_vars, with = FALSE]
  complete_cases <- complete.cases(cluster_matrix)
  
  if (sum(complete_cases) < 50) {
    cat("Too few complete cases for reliable MDS\n")
    return(NULL)
  }
  
  cluster_matrix_complete <- cluster_matrix[complete_cases, ]
  clusters_complete <- clusters[complete_cases]
  data_complete <- data_with_clusters[complete_cases, ]
  
  tryCatch({
    # Calculate Gower distance matrix
    dist_matrix <- cluster::daisy(cluster_matrix_complete, metric = "gower")
    
    cat(sprintf("Calculated distance matrix for %d observations\n", nrow(cluster_matrix_complete)))
    
    # Perform ordinal MDS
    classical_mds <- cmdscale(dist_matrix, k = 2)
    ordinal_mds <- isoMDS(dist_matrix, y = classical_mds, k = 2, 
                          maxit = 100, trace = FALSE)
    
    stress_value <- ordinal_mds$stress
    stress_interpretation <- if (stress_value < 0.05) "Excellent" 
                            else if (stress_value < 0.10) "Good" 
                            else if (stress_value < 0.20) "Fair" 
                            else "Poor"
    
    cat(sprintf("Ordinal MDS Stress: %.4f (%s)\n", stress_value, stress_interpretation))
    
    # Create plotting data
    plot_data <- data.frame(
      MDS1 = ordinal_mds$points[, 1],
      MDS2 = ordinal_mds$points[, 2],
      Cluster = as.factor(clusters_complete),
      RESILIENT = as.factor(data_complete$RESILIENT),
      ESCS = data_complete$ESCS
    )
    
    # Calculate cluster-specific resilience rates for labels
    cluster_labels <- c()
    for (cluster_id in sort(unique(clusters_complete))) {
      cluster_mask <- clusters_complete == cluster_id
      resilience_rate <- mean(data_complete$RESILIENT[cluster_mask], na.rm = TRUE) * 100
      n_students <- sum(cluster_mask)
      
      label <- sprintf("Cluster %d\n(n=%d, %.1f%% resilient)", 
                       cluster_id, n_students, resilience_rate)
      cluster_labels <- c(cluster_labels, label)
    }
    
    # Main MDS plot with clusters and resilience
    p1 <- ggplot(plot_data, aes(x = MDS1, y = MDS2)) +
      geom_point(aes(color = Cluster, shape = RESILIENT), 
                 alpha = 0.7, size = 2.5) +
      stat_ellipse(aes(color = Cluster), type = "norm", level = 0.68, alpha = 0.3) +
      scale_color_viridis_d(name = "Cluster", labels = cluster_labels) +
      scale_shape_manual(name = "Resilient", 
                         values = c("0" = 16, "1" = 17),
                         labels = c("0" = "Not Resilient", "1" = "Resilient")) +
      labs(
        title = sprintf("MDS Clustering Visualization: PISA %s (Intervention)", cycle_name),
        subtitle = sprintf("Gower Distance | Stress: %.4f (%s) | %d variables", 
                           stress_value, stress_interpretation, length(clustering_vars)),
        x = "MDS Dimension 1",
        y = "MDS Dimension 2",
        caption = sprintf("Akhanli & Hennig Framework | %d complete observations", nrow(plot_data))
      ) +
      theme_minimal() +
      theme(
        plot.title = element_text(size = 14, face = "bold"),
        plot.subtitle = element_text(size = 12),
        legend.position = "right",
        panel.grid.minor = element_blank()
      )
    
    # ESCS overlay plot
    p2 <- ggplot(plot_data, aes(x = MDS1, y = MDS2)) +
      geom_point(aes(color = ESCS), size = 2.5, alpha = 0.8) +
      scale_color_gradient2(name = "ESCS\n(SES Index)", 
                            low = "red", mid = "white", high = "blue", 
                            midpoint = median(plot_data$ESCS, na.rm = TRUE)) +
      stat_ellipse(aes(group = Cluster), color = "gray50", linetype = "dashed", alpha = 0.5) +
      labs(
        title = sprintf("Socioeconomic Status Distribution: PISA %s", cycle_name),
        subtitle = sprintf("Overlaid on MDS coordinates | Stress: %.4f", stress_value),
        x = "MDS Dimension 1", y = "MDS Dimension 2"
      ) +
      theme_minimal() +
      theme(panel.grid.minor = element_blank())
    
    # Print plots
    print(p1)
    print(p2)
    
    return(list(
      main = p1, 
      escs = p2, 
      stress = stress_value, 
      n_obs = nrow(plot_data),
      cluster_labels = cluster_labels
    ))
    
  }, error = function(e) {
    cat(sprintf("Error in MDS calculation: %s\n", e$message))
    return(NULL)
  })
}

# Generate MDS plots for all cycles
mds_plots <- list()
for (cycle_name in names(enhanced_data_with_clusters)) {
  if ("CLUSTER_INTERVENTION" %in% names(enhanced_data_with_clusters[[cycle_name]]$data)) {
    plots <- create_mds_cluster_plot(enhanced_data_with_clusters, cycle_name)
    if (!is.null(plots)) {
      mds_plots[[cycle_name]] <- plots
    }
  }
}

# MDS Summary
cat("\n=== MDS VISUALIZATION SUMMARY ===\n")
for (cycle_name in names(mds_plots)) {
  plots <- mds_plots[[cycle_name]]
  cat(sprintf("✓ %s: Stress=%.4f (%s), N=%d observations\n", 
              cycle_name, plots$stress, 
              if (plots$stress < 0.05) "Excellent" else if (plots$stress < 0.10) "Good" else if (plots$stress < 0.20) "Fair" else "Poor",
              plots$n_obs))
}

cat("\n✓ MDS visualizations completed\n")
cat("✓ Using Gower distance for mixed-type variables\n")
cat("✓ Resilience patterns overlaid on MDS coordinates\n")
cat("✓ Confidence ellipses show cluster boundaries\n")
```

