# Shared helpers for the PC-field association matrix and the field x field
# collinearity matrix. Both run smartcor::smart_cormat() with its defaults, so
# method selection is smartcor's own (see ?smart_cor for the type-pair table).
#
# Deviations from the defaults:
#   * types = <explicit>: our FIELD_TYPES are mapped to smartcor's vocabulary
#     instead of being re-detected, so the two classifications cannot drift.
#   * bootstrap = FALSE: only affects confidence intervals, which we don't use
#     (p-values stay analytic under every setting, per ?smart_cor).
#
# Caveat: continuous x categorical pairs use Cramer's V on a quantile-binned
# continuous variable (?smart_cor). With few samples the binned table is
# sparse, so those p-values are unreliable; the `method` column lets you
# filter them.

# FIELD_TYPES (Snakefile) use smartcor's vocabulary, plus two types of our own
# that carry no diagnostic value and are dropped before smartcor is called.
SMARTCOR_FIELD_TYPES <- c("binary", "categorical", "ordinal", "continuous")
EXCLUDED_FIELD_TYPES <- c("constant", "identifier")

# Fields that participate in diagnostics: exactly those in `field_types`
# (Snakefile FIELD_TYPES), minus constant/identifier. Returns a named
# character vector (name = field).
diagnostic_field_types <- function(coldata, field_types) {
  ft <- unlist(field_types)
  absent <- setdiff(names(ft), colnames(coldata))
  if (length(absent) > 0) {
    stop("ERROR: field_types lists field(s) not present in coldata: ",
         paste(absent, collapse = ", "))
  }
  ft <- ft[!(ft %in% EXCLUDED_FIELD_TYPES)]
  unknown <- setdiff(unique(ft), SMARTCOR_FIELD_TYPES)
  if (length(unknown) > 0) {
    stop("ERROR: unsupported field type(s): ",
         paste(unknown, collapse = ", "))
  }
  ft
}

# as.character() before as.numeric(): coldata fields are factors (relevel() in
# deseq2_fit.R) and as.numeric() on a factor returns level codes.
prepare_smartcor_vector <- function(x, type, name) {
  switch(
    type,
    continuous = {
      v <- suppressWarnings(as.numeric(as.character(x)))
      if (any(is.na(v) & !is.na(x))) {
        stop("ERROR: field '", name, "' is typed continuous but has ",
             "non-numeric values.")
      }
      v
    },
    binary = {
      f <- factor(as.character(x))
      if (nlevels(f) > 2) {
        stop("ERROR: field '", name, "' is typed binary but has ",
             nlevels(f), " levels.")
      }
      as.integer(f) - 1L
    },
    categorical = factor(as.character(x)),
    ordinal = {
      v <- suppressWarnings(as.numeric(as.character(x)))
      if (identical(is.na(v), is.na(x))) {
        v
      } else {
        message("Ordinal field '", name, "' is non-numeric; ordering levels ",
                "alphabetically. Recode numerically if that is wrong.")
        factor(as.character(x), levels = sort(unique(as.character(x))),
               ordered = TRUE)
      }
    },
    stop("ERROR: unsupported field type '", type, "' for field '", name, "'")
  )
}

# Fields -> data frame of smartcor-ready columns plus the matching `types` list. `field_types` is the named vector from diagnostic_field_types().
prepare_association_data <- function(coldata, field_types) {
  columns <- lapply(names(field_types), function(f) {
    prepare_smartcor_vector(coldata[[f]], field_types[[f]], f)
  })
  names(columns) <- names(field_types)
  data <- as.data.frame(columns, check.names = FALSE, stringsAsFactors = FALSE)
  rownames(data) <- rownames(coldata)
  list(
    data = data,
    types = as.list(field_types)
  )
}

# smart_cormat() over all column pairs -> one row per unique pair. Failed pairs
# come back NA; smartcor reports why as warnings (verbose = TRUE, the default),
# which end up in this rule's log.
compute_associations <- function(data, types) {
  matrix <- smartcor::smart_cormat(data, types = types, bootstrap = FALSE)
  tidy_matrix <- as.data.frame(smartcor::tidy(matrix))
  needed <- c("var_x", "var_y", "x_type", "y_type", "estimate", "method",
              "method_label", "p.value", "n_complete")
  missing <- setdiff(needed, colnames(tidy_matrix))
  if (length(missing) > 0) {
    stop("ERROR: tidy(smart_cormat()) lacks expected column(s): ",
         paste(missing, collapse = ", "), ". Found: ",
         paste(colnames(tidy_matrix), collapse = ", "))
  }
  data.frame(
    var_x = tidy_matrix$var_x, var_y = tidy_matrix$var_y,
    type_x = tidy_matrix$x_type, type_y = tidy_matrix$y_type,
    estimate = tidy_matrix$estimate, p_value = tidy_matrix$p.value,
    method = tidy_matrix$method, method_label = tidy_matrix$method_label,
    n_complete = tidy_matrix$n_complete, stringsAsFactors = FALSE
  )
}
