#' Make unique names
#'
#' \code{make_unique} generates unique identifiers
#' for a proteomics dataset based on "name" and "id" columns.
#'

make_unique <- function(proteins, names, ids, delim = ";") {
  # Show error if inputs are not the required classes
  assertthat::assert_that(is.data.frame(proteins),
    is.character(names),
    length(names) == 1,
    is.character(ids),
    length(ids) == 1,
    is.character(delim),
    length(delim) == 1)

  col_names <- colnames(proteins)
  # Show error if inputs do not contain required columns
  if(!names %in% col_names) {
    stop("'", names, "' is not a column in '",
         deparse(substitute(proteins)), "'",
         call. = FALSE)
  }
  if(!ids %in% col_names) {
    stop("'", ids, "' is not a column in '",
         deparse(substitute(proteins)), "'",
         call. = FALSE)
  }

  # If input is a tibble, convert to data.frame
  if(tibble::is_tibble(proteins))
    proteins <- as.data.frame(proteins)

  # Select the name and id columns, and check for NAs
  double_NAs <- apply(proteins[,c(names, ids)], 1, function(x) all(is.na(x)))
  if(any(double_NAs)) {
    stop("NAs in both the 'names' and 'ids' columns")
  }

  # Take the first identifier per row and make unique names.
  # If there is no name, the ID will be taken.
  proteins_unique <- proteins %>%
    mutate(name = gsub(paste0(delim, ".*"), "", get(names)),
      ID = gsub(paste0(delim, ".*"), "", get(ids)),
      name = make.unique(ifelse(name == "" | is.na(name), ID, name)))
  return(proteins_unique)
}
#' Obtain the longest common prefix
#'
#' \code{get_prefix} returns the longest common prefix
#' of the supplied words.
get_prefix <- function(words) {
  # Show error if input is not the required class
  assertthat::assert_that(is.character(words))

  # Show error if 'words' contains 1 or less elements
  if(length(words) <= 1) {
    stop("'words' should contain more than one element")
  }
  # Show error if 'words' contains NA
  if(any(is.na(words))) {
    stop("'words' contains NAs")
  }

  # Truncate words to smallest name
  minlen <- min(nchar(words))
  truncated <- substr(words, 1, minlen)

  # Show error if one of the elements is shorter than one character
  if(minlen < 1) {
    stop("At least one of the elements is too short")
  }

  # Get identifical characters
  mat <- data.frame(strsplit(truncated, ""), stringsAsFactors = FALSE)
  identical <- apply(mat, 1, function(x) length(unique(x)) == 1)

  # Obtain the longest common prefix
  prefix <- as.logical(cumprod(identical))
  paste(mat[prefix, 1], collapse = "")
}

#' Obtain the longest common suffix
#'
#' \code{get_suffix} returns the longest common suffix
#' of the supplied words.

get_suffix <- function(words) {
  # Show error if input is not the required class
  assertthat::assert_that(is.character(words))

  # Show error if 'words' contains 1 or less elements
  if(length(words) <= 1) {
    stop("'words' should contain more than one element")
  }
  # Show error if 'words' contains NA
  if(any(is.na(words))) {
    stop("'words' contains NAs")
  }

  # Truncate words to smallest name
  minlen <- min(nchar(words))
  truncated <- substr(words, nchar(words) - minlen + 1, nchar(words))

  # Show error if one of the elements is shorter than one character
  if(minlen < 1) {
    stop("At least one of the elements is too short")
  }

  # Reverse characters wihtin word
  rev_string <- function(str) {
    paste(rev(strsplit(str, "")[[1]]), collapse = "")
  }
  rev_truncated <- vapply(truncated, rev_string, character(1))

  # Get identifical characters
  mat <- data.frame(strsplit(rev_truncated, ""), stringsAsFactors = FALSE)
  identical <- apply(mat, 1, function(x) length(unique(x)) == 1)

  # Obtain the longest common prefix
  prefix <- as.logical(cumprod(identical))
  rev_string(paste(mat[prefix, 1], collapse = ""))
}
# Short internal function to delete the longest common prefix
delete_prefix <- function(words) {
  # Get prefix
  prefix <- get_prefix(words)
  # Delete prefix from words
  gsub(paste0("^", prefix), "", words)
}

# Short internal function to delete the longest common suffix
delete_suffix <- function(words) {
  # Get prefix
  suffix <- get_suffix(words)
  # Delete prefix from words
  gsub(paste0(suffix, "$"), "", words)
}
#' Data.frame to SummarizedExperiment object
#' conversion using an experimental design
#'
#' \code{make_se} creates a SummarizedExperiment object
#' based on two data.frames: the protein table and experimental design.
#'

make_se <- function(proteins_unique, columns, expdesign) {
  # Show error if inputs are not the required classes
  assertthat::assert_that(is.data.frame(proteins_unique),
                          is.integer(columns),
                          is.data.frame(expdesign))

  # Show error if inputs do not contain required columns
  if(any(!c("name", "ID") %in% colnames(proteins_unique))) {
    stop("'name' and/or 'ID' columns are not present in '",
         deparse(substitute(proteins_unique)),
         "'.\nRun make_unique() to obtain the required columns",
         call. = FALSE)
  }
  if(any(!c("label", "condition", "replicate") %in% colnames(expdesign))) {
    stop("'label', 'condition' and/or 'replicate' columns",
         "are not present in the experimental design",
         call. = FALSE)
  }
  if(any(!apply(proteins_unique[, columns], 2, is.numeric))) {
    stop("specified 'columns' should be numeric",
         "\nRun make_se_parse() with the appropriate columns as argument",
         call. = FALSE)
  }

  # If input is a tibble, convert to data.frame
  if(tibble::is_tibble(proteins_unique))
    proteins_unique <- as.data.frame(proteins_unique)
  if(tibble::is_tibble(expdesign))
    expdesign <- as.data.frame(expdesign)

  # Select the assay data
  rownames(proteins_unique) <- proteins_unique$name
  raw <- proteins_unique[, columns]
  raw[raw == 0] <- NA
  raw <- log2(raw)
  

  # Generate the colData from the experimental design
  # and match these with the assay data
  expdesign <- mutate(expdesign, condition = make.names(condition)) %>%
    unite(ID, condition, replicate, remove = FALSE)
  rownames(expdesign) <- expdesign$ID

  matched <- match(make.names(delete_prefix(expdesign$label)),
                   make.names(delete_prefix(colnames(raw))))
  if(any(is.na(matched))) {
    stop("None of the labels in the experimental design match ",
         "with column names in 'proteins_unique'",
         "\nRun make_se() with the correct labels in the experimental design",
         "and/or correct columns specification")
  }

  colnames(raw)[matched] <- expdesign$ID
  raw <- raw[, !is.na(colnames(raw))][rownames(expdesign)]

  # Select the rowData
  row_data <- proteins_unique[, -columns]
  rownames(row_data) <- row_data$name

  # Generate the SummarizedExperiment object
  se <- SummarizedExperiment(assays = as.matrix(raw),
                                  colData = expdesign,
                                  rowData = row_data)
  return(se)
}

#' Data.frame to SummarizedExperiment object
#' conversion using parsing from column names
#'
#' \code{make_se_parse} creates a SummarizedExperiment object
#' based on a single data.frame.
#'

make_se_parse <- function(proteins_unique, columns,
                          mode = c("char", "delim"), chars = 1, sep = "_") {
  # Show error if inputs are not the required classes
  assertthat::assert_that(is.data.frame(proteins_unique),
                          is.integer(columns),
                          is.character(mode),
                          is.numeric(chars),
                          length(chars) == 1,
                          is.character(sep),
                          length(sep) == 1)

  # Show error if inputs do not contain required columns
  mode <- match.arg(mode)

  # Show error if inputs do not contain required columns
  if(any(!c("name", "ID") %in% colnames(proteins_unique))) {
    stop("'name' and/or 'ID' columns are not present in '",
         deparse(substitute(proteins_unique)),
         "'.\nRun make_unique() to obtain the required columns",
         call. = FALSE)
  }
  if(any(!apply(proteins_unique[, columns], 2, is.numeric))) {
    stop("specified 'columns' should be numeric",
         "\nRun make_se_parse() with the appropriate columns as argument",
         call. = FALSE)
  }

  # If input is a tibble, convert to data.frame
  if(tibble::is_tibble(proteins_unique))
    proteins_unique <- as.data.frame(proteins_unique)

  # Select the assay values
  rownames(proteins_unique) <- proteins_unique$name
  raw <- proteins_unique[, columns]
  raw[raw == 0] <- NA
  raw <- log2(raw)
  colnames(raw) <- delete_prefix(colnames(raw)) %>% make.names()

  # Select the rowData
  row_data <- proteins_unique[, -columns]
  rownames(row_data) <- row_data$name

  # Generate the colData
  if(mode == "char") {
    col_data <- data.frame(label = colnames(raw), stringsAsFactors = FALSE) %>%
      mutate(condition = substr(label, 1, nchar(label) - chars),
             replicate = substr(label, nchar(label) + 1 - chars, nchar(label))) %>%
      unite(ID, condition, replicate, remove = FALSE)
  }
  if(mode == "delim") {
    col_data <- data.frame(label = colnames(raw), stringsAsFactors = FALSE) %>%
      separate(label, c("condition", "replicate"), sep = sep,
               remove = FALSE, extra = "merge") %>%
      unite(ID, condition, replicate, remove = FALSE)
  }
  rownames(col_data) <- col_data$ID
  colnames(raw)[match(col_data$label, colnames(raw))] <- col_data$ID
  raw <- raw[, !is.na(colnames(raw))]


  # Generate the SummarizedExperiment object
  se <- SummarizedExperiment(assays = as.matrix(raw),
                             colData = col_data,
                             rowData = row_data)
  return(se)
}


#' Filter on missing values
#'
#' \code{filter_missval} filters a proteomics dataset based on missing values.
#' The dataset is filtered for proteins that have a maximum of
#' 'thr' missing values in at least one condition.
#'

filter_missval <- function(se, thr = 0) {
  # Show error if inputs are not the required classes
  if(is.integer(thr)) thr <- as.numeric(thr)
  assertthat::assert_that(inherits(se, "SummarizedExperiment"),
                          is.numeric(thr),
                          length(thr) == 1)

  # Show error if inputs do not contain required columns
  if(any(!c("name", "ID") %in% colnames(rowData(se, use.names = FALSE)))) {
    stop("'name' and/or 'ID' columns are not present in '",
         deparse(substitute(se)),
         "'\nRun make_unique() and make_se() to obtain the required columns",
         call. = FALSE)
  }
  if(any(!c("label", "condition", "replicate") %in% colnames(colData(se)))) {
    stop("'label', 'condition' and/or 'replicate' columns are not present in '",
         deparse(substitute(se)),
         "'\nRun make_se() or make_se_parse() to obtain the required columns",
         call. = FALSE)
  }
  max_repl <- max(colData(se)$replicate)
  if(thr < 0 | thr > max_repl) {
    stop("invalid filter threshold applied",
         "\nRun filter_missval() with a threshold ranging from 0 to ",
         max_repl)
  }

  # Make assay values binary (1 = valid value)
  bin_data <- assay(se)
  idx <- is.na(assay(se))
  bin_data[!idx] <- 1
  bin_data[idx] <- 0

  # Filter se on the maximum allowed number of
  # missing values per condition (defined by thr)
  keep <- bin_data %>%
    data.frame() %>%
    rownames_to_column() %>%
    gather(ID, value, -rowname) %>%
    left_join(., data.frame(colData(se)), by = "ID") %>%
    group_by(rowname, condition) %>%
    summarize(miss_val = n() - sum(value)) %>%
    filter(miss_val <= thr) %>%
    spread(condition, miss_val)
  se_fltrd <- se[keep$rowname, ]
  return(se_fltrd)
}

#' Filter proteins based on missing values
#'
#' \code{filter_proteins} filters a proteomic dataset based on missing values.
#' Different types of filtering can be applied, which range from only keeping
#' proteins without missing values to keeping proteins with a certain percent
#' valid values in all samples or keeping proteins that are complete
#' in at least one condition.
#'

filter_proteins <- function(se, type = c("complete", "condition", "fraction"),
                            thr = NULL, min = NULL) {
  # Show error if inputs are not the required classes
  assertthat::assert_that(inherits(se, "SummarizedExperiment"))
  type <- match.arg(type)

  # Show error if inputs do not contain required columns
  if(any(!c("name", "ID") %in% colnames(rowData(se, use.names = FALSE)))) {
    stop("'name' and/or 'ID' columns are not present in '",
         deparse(substitute(se)),
         "'\nRun make_unique() and make_se() to obtain the required columns",
         call. = FALSE)
  }
  if(any(!c("label", "condition", "replicate") %in% colnames(colData(se)))) {
    stop("'label', 'condition' and/or 'replicate' columns are not present in '",
         deparse(substitute(se)),
         "'\nRun make_se() or make_se_parse() to obtain the required columns",
         call. = FALSE)
  }

  if(type == "complete") {
    keep <- !apply(assay(se), 1, function(x) any(is.na(x)))
    filtered <- se[keep,]
  }
  if(type == "condition") {
    assertthat::assert_that(is.numeric(thr),
                            length(thr) == 1)
    max_repl <- max(colData(se)$replicate)
    if(thr < 0 | thr > max_repl) {
      stop("invalid filter threshold 'thr' applied",
           "\nRun filter() with a threshold ranging from 0 to ",
           max_repl)
    }

    filtered <- filter_missval(se, thr = thr)
  }
  if(type == "fraction") {
    assertthat::assert_that(is.numeric(min),
                            length(min) == 1)
    if(min < 0 | min > 1) {
      stop("invalid filter threshold 'min' applied",
           "\nRun filter() with a percent ranging from 0 to 1")
    }

    bin_data <- assay(se)
    idx <- is.na(assay(se))
    bin_data[!idx] <- 1
    bin_data[idx] <- 0

    # Filter se on the maximum allowed number of
    # missing values per condition (defined by thr)
    keep <- bin_data %>%
      as.data.frame() %>%
      rownames_to_column() %>%
      gather(ID, value, -rowname) %>%
      group_by(rowname) %>%
      summarize(n = n(),
                valid = sum(value),
                frac = valid / n) %>%
      filter(frac >= min)
    filtered <- se[keep$rowname, ]
  }
  return(filtered)
}

#' Normalization using vsn
#'
#' \code{normalize_vsn} performs variance stabilizing transformation
#' using the \code{\link[vsn]{vsn-package}}.
#'

normalize_vsn <- function(se) {
  # Show error if inputs are not the required classes
  assertthat::assert_that(inherits(se, "SummarizedExperiment"))

  # Variance stabilization transformation on assay data
  se_vsn <- se
  vsn.fit <- vsn::vsnMatrix(2 ^ assay(se_vsn))
  assay(se_vsn) <- vsn::predict(vsn.fit, 2 ^ assay(se_vsn))
  return(se_vsn)
}
#' Summarising phosphosites to proteins
#'
#' Summarising phosphosite-level information to proteins for performing
#' downstream gene-centric analyses.
#'
ptmCollapse <- function(mat, id, stat, by = "min") {
    listMat <- split(data.frame(mat, stat), id)

    matNew <- c()
    if (by == "min") {
        matNew <- as.matrix(do.call(rbind, lapply(listMat, function(x) {
            if (nrow(x) == 1) {
                as.numeric(x[1, seq_len(ncol(x) - 1)])
            } else {
                x[order(as.numeric(x[, ncol(x)]))[1], seq_len(ncol(x) - 1)]
            }
        })))
    } else if (by == "max") {
        matNew <- as.matrix(do.call(rbind, lapply(listMat, function(x) {
            if (nrow(x) == 1) {
                as.numeric(x[1, seq_len(ncol(x) - 1)])
            } else {
                x[order(as.numeric(x[, ncol(x)]), decreasing = TRUE)[1],
                    seq_len(ncol(x) - 1)]
            }
        })))
    } else if (by == "mid") {
        matNew <- as.matrix(do.call(rbind, lapply(listMat, function(x) {
            if (nrow(x) == 1) {
                as.numeric(x[1, seq_len(ncol(x) - 1)])
            } else {
                mid <- round(nrow(x)/2)
                x[order(as.numeric(x[, ncol(x)]))[mid], seq_len(ncol(x) - 1)]
            }
        })))
    } else {
        stop("Unrecognised way of collapsing the data!")
    }

    return(matNew)
}
