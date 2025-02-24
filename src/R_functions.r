library(edgeR)
library(limma)
library(dplyr)


#' Use edgeR's filterByExpr to remove low count genes from count matrix
#'
#' @param inpath String specifying path to csv file of gene expression data
#' @param outpath String specifying path to store filtered csv file
#' @param design Design matrix
edgeR_filterByExpression <- function(inpath, outpath, design) {
  x <- read.csv(inpath, row.names = 1)
  x <- data.matrix(x, rownames.force = TRUE)
  y <- DGEList(counts = x)

  if (design == "paired") {
    if (ncol(x) %% 2 != 0) {
      stop("Paired-design matrix must have even number of columns")
    }
    N <- ncol(x) / 2
    patient <- factor(c(seq(N), seq(N)))
    condition <- factor(c(rep("N", N), rep("T", N))) # normal vs tumor (control vs treatment)
    design_mat <- model.matrix(~ patient + condition)
  } else if (design == "unpaired") {
    if (ncol(x) %% 2 != 0) {
      stop("Design matrix must have even number of columns")
    }
    N <- ncol(x) / 2
    condition <- factor(c(rep("N", N), rep("T", N))) # normal vs tumor (control vs treatment)
    design_mat <- model.matrix(~condition)
  }

  if (design == "none") {
    keep <- filterByExpr(y)
  } else {
    keep <- filterByExpr(y, design = design_mat)
  }

  y <- y[keep, , keep.lib.sizes = FALSE]
  write.csv(y$counts, outpath, row.names = TRUE)
}

#' Run edgeR
#'
#' @param x: dataframe of counts
#' @param design: design matrix, if "paired" constructs design matrix from data assuming x is of the form: k control cols followed by k treatment cols
#' @param overwrite: logical, whether to overwrite existing results table if present
#' @param filter_expr: logical, whether to remove low counts using edgeR's filterByExpr function
#' @param top_tags: int or "Inf", store the results of the most significant genes only
#' @param lfc: float, logFC threshold when testing for DE
#' @param cols_to_keep: list of output table columns to save
run_edgeR <- function(x, outfile, design, overwrite = FALSE, filter_expr = FALSE, top_tags = "Inf", verbose = FALSE,
                      lfc = 0, cols_to_keep = "all", test = "qlf", meta_only = FALSE, check_gof = FALSE, N_control = 0, N_treat = 0) {
  suppressPackageStartupMessages(require("edgeR"))
  suppressPackageStartupMessages(require("limma"))

  # Check if files already exist
  if (!overwrite && file.exists(outfile)) {
    print("Existing table not overwritten")
    return()
  }

  if (design == "paired") {
    if (ncol(x) %% 2 != 0) {
      stop("Paired-design matrix must have even number of columns")
    }
    N <- ncol(x) / 2
    patient <- factor(c(seq(N), seq(N)))
    condition <- factor(c(rep("N", N), rep("T", N))) # normal vs tumor (control vs treatment)
    design <- model.matrix(~ patient + condition)
  } else if (design == "unpaired") {
    if (ncol(x) %% 2 != 0) {
      stop("Design matrix must have even number of columns")
    }
    N <- ncol(x) / 2
    condition <- factor(c(rep("N", N), rep("T", N)))
    design <- model.matrix(~condition)
  } else if (design == "unpaired_asymmetric") {
    if (N_control == 0) {
      stop("Design matrix has no control columns")
    }
    if (N_treat == 0) {
      stop("Design matrix has no treatment/condition columns")
    }
    condition <- factor(c(rep("N", N_control), rep("T", N_treat)))
    design <- model.matrix(~condition)
  } else if (grepl("\\.csv$", design, ignore.case = TRUE)) {
    print("Constructing design matrix from df")
    covariate_df <- read.csv(design)
    if (!("Condition" %in% colnames(covariate_df))) {
      stop("Error: 'Condition' column not found in dataframe")
    }

    covariate_df <- covariate_df %>%
      mutate_if(is.character, as.factor)

    other_vars <- setdiff(names(covariate_df), c("Condition", "X", "Sample"))
    print("Warning: hard-coded col names in design matrix")
    formula <- as.formula(paste("~", paste(c(other_vars, "Condition"), collapse = " + ")))
    print(paste("Formula:", formula))
    design <- model.matrix(formula, data = covariate_df)

    # Ensure your count matrix and metadata align
    if ("X" %in% names(covariate_df)) {
      rownames(covariate_df) <- covariate_df$X
    } else if ("Sample" %in% names(covariate_df)) {
      rownames(covariate_df) <- covariate_df$Sample
    } else {
      stop("Sample names not found in covariate df")
    }
  }

  if (verbose) {
    rank <- qr(design)$rank
    print(paste("Rank:", rank, "Cols:", ncol(design)))
    print(design)
  }

  y <- DGEList(counts = x)

  print(length(rownames(design)))
  print(length(colnames(y)))
  rownames(design) <- colnames(y)

  if (filter_expr) {
    keep <- filterByExpr(y, design = design)
    y <- y[keep, , keep.lib.sizes = FALSE]
  }

  y <- calcNormFactors(y)
  y <- estimateDisp(y, design, robust = TRUE)

  if (meta_only) {
    return(y)
  }

  if (test == "lrt") {
    fit <- glmFit(y, design)
  } else {
    fit <- glmQLFit(y, design)
  }

  # Goodness-of-fit
  if (check_gof) {
    res.gof <- gof(fit, plot = FALSE)
    file_name <- basename(outfile)
    new_file_name <- paste0("gof.", file_name)
    new_file_path <- file.path(dirname(outfile), new_file_name)
    write.csv(res.gof$gof.pvalues, new_file_path)
    print(paste("Saved gof in", new_file_path))
  }

  if (lfc > 0) {
    result <- glmTreat(fit, lfc = lfc)
  } else if (test == "lrt") {
    result <- glmLRT(fit)
  } else {
    result <- glmQLFTest(fit)
  } # omit coef (edgeR user's guide p. 39)

  table <- topTags(result, n = top_tags) # adjust.method="BH"

  if (any(cols_to_keep != "all")) {
    if (typeof(cols_to_keep) == "list") cols_to_keep <- unlist(cols_to_keep)
    table <- table[, cols_to_keep]
  }
  write.csv(table, outfile)
}


# DESeq2
run_deseq2 <- function(
    x, outfile, design = "paired", overwrite = FALSE, print_summary = FALSE, cols_to_keep = "all",
    size_factors_only = FALSE, lfc = 0, shrink_lfc = FALSE, shrink_method = "apeglm") {
  if (!overwrite && file.exists(outfile)) {
    print("Existing table not overwritten")
    return()
  }

  suppressPackageStartupMessages(require("DESeq2"))

  if (design == "paired") {
    if (ncol(x) %% 2 != 0) {
      stop("Paired-design matrix must have even number of columns")
    }
    N <- ncol(x) / 2
    patients <- factor(c(seq(N), seq(N)))
    Condition <- factor(c(rep("N", N), rep("T", N))) # normal vs tumor (control vs treatment)
    array <- array(factor(c(patients, Condition)), dim = c(length(patients), 2))
    coldata <- data.frame(array, row.names = colnames(x))
    colnames(coldata) <- c("patient", "Condition")

    dds <- DESeqDataSetFromMatrix(
      countData = x,
      colData = coldata,
      design = ~ patient + Condition
    )
  } else if (design == "unpaired") {
    if (ncol(x) %% 2 != 0) {
      stop("Design matrix must have even number of columns")
    }
    N <- ncol(x) / 2
    Condition <- factor(c(rep("N", N), rep("T", N))) # normal vs tumor (control vs treatment)
    array <- array(factor(Condition), dim = c(length(Condition)))
    coldata <- data.frame(array, row.names = colnames(x))
    colnames(coldata) <- c("Condition")
    dds <- DESeqDataSetFromMatrix(
      countData = x,
      colData = coldata,
      design = ~Condition
    )
  } else {
    print("Constructing design matrix from df")
    covariate_df <- read.csv(design)
    if (!("Condition" %in% colnames(covariate_df))) {
      stop("Error: 'Condition' column not found in dataframe")
    }

    covariate_df <- covariate_df %>%
      mutate_if(is.character, as.factor)

    other_vars <- setdiff(names(covariate_df), c("Condition", "X", "Sample"))
    print("Warning: hard-coded col names in design matrix")
    formula <- as.formula(paste("~", paste(c(other_vars, "Condition"), collapse = " + ")))
    print(paste("Formula:", formula))
    design <- model.matrix(formula, data = covariate_df)

    # Ensure your count matrix and metadata align
    if ("X" %in% names(covariate_df)) {
      rownames(covariate_df) <- covariate_df$X
    } else if ("Sample" %in% names(covariate_df)) {
      rownames(covariate_df) <- covariate_df$Sample
    } else {
      stop("Sample names not found in covariate df")
    }

    x <- x[, rownames(covariate_df)] # Align count data with metadata
    # Create DESeqDataSet object
    dds <- DESeqDataSetFromMatrix(
      countData = x,
      colData = covariate_df,
      design = formula
    )
  }



  if (size_factors_only) {
    return(sizeFactors(estimateSizeFactors(dds)))
  }

  dds <- DESeq(dds)
  contrastname <- resultsNames(dds)[grepl("Condition", resultsNames(dds))]
  res <- results(dds, name = contrastname, lfcThreshold = lfc, altHypothesis = "greaterAbs", test = "Wald")

  if (shrink_lfc) {
    print("shrinking lfc")
    res <- lfcShrink(dds, coef = contrastname, type = shrink_method)
  }

  if (print_summary) {
    print(summary(res))
  }

  names(res)[names(res) == "padj"] <- "FDR" # for consistency with edgeR
  names(res)[names(res) == "log2FoldChange"] <- "logFC"
  # names(res)[names(res) == "pvalue"] <- "PValue"

  if (any(cols_to_keep != "all")) {
    if (typeof(cols_to_keep) == "list") cols_to_keep <- unlist(cols_to_keep)
    res <- res[, cols_to_keep]
  }

  write.csv(res, outfile)
}
