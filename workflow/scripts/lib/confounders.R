# Shared helper for scripts that put confounder fields into a design.
#
# Every coldata column is a string parsed from the sample label (and
# deseq2_fit.R turns each field into a factor via relevel()), while
# model.matrix() decides factor-vs-numeric from the R type. A confounder
# declared "continuous" must therefore be cast explicitly, or it gets one
# dummy column per distinct value.
#
# Returns `coldata` with the columns of all confounders declared
# "continuous" converted to numeric; every other column is left untouched.
# `coldata` can be a data.frame or a DataFrame (e.g. colData(dds)); the
# return value has the same class.
cast_continuous_confounders <- function(coldata, confounder_field_names,
                                        confounder_field_types) {
  for (f in confounder_field_names) {
    if (identical(confounder_field_types[[f]], "continuous")) {
      # as.character() first: as.numeric() on a factor returns the level
      # codes, not the values.
      values <- as.character(coldata[[f]])
      numeric_values <- suppressWarnings(as.numeric(values))
      if (anyNA(numeric_values)) {
        stop(
          "ERROR: confounder '", f, "' is declared continuous but has ",
          "non-numeric value(s): ",
          paste(unique(values[is.na(numeric_values)]), collapse = ", ")
        )
      }
      coldata[[f]] <- numeric_values
    }
  }
  coldata
}
