# ============================================================
# Reference preparation
# ============================================================

weighted_af <- function(data, af_columns, an_columns) {
  if (length(af_columns) != length(an_columns)) {
    stop("AF and AN column lists must have equal length.")
  }

  missing_columns <- setdiff(
    c(af_columns, an_columns),
    names(data)
  )

  if (length(missing_columns) > 0L) {
    stop(
      "Missing reference columns: ",
      paste(missing_columns, collapse = ", ")
    )
  }

  weighted_sum <- Reduce(
    `+`,
    Map(
      function(af, an) data[[af]] * data[[an]],
      af_columns,
      an_columns
    )
  )

  total_an <- Reduce(
    `+`,
    lapply(an_columns, function(an) data[[an]])
  )

  weighted_sum / total_an
}


prepare_reference_data <- function(input_file) {
  chr22 <- read.table(
    input_file,
    header = TRUE
  )

  chr22$AF_AFR <- weighted_af(
    chr22,
    af_columns = c(
      "AFesn", "AFgwd", "AFlwk", "AFmsl",
      "AFyri", "AFbantukenya", "AFmandenka", "AFyoruba"
    ),
    an_columns = c(
      "ANesn", "ANgwd", "ANlwk", "ANmsl",
      "ANyri", "ANbantukenya", "ANmandenka", "ANyoruba"
    )
  )

  chr22$AF_EUR <- weighted_af(
    chr22,
    af_columns = c(
      "AFceu", "AFgbr", "AFibs", "AFtsi", "AFbasque",
      "AFfrench", "AFitalian", "AForcadian",
      "AFsardinian", "AFtuscan"
    ),
    an_columns = c(
      "ANceu", "ANgbr", "ANibs", "ANtsi", "ANbasque",
      "ANfrench", "ANitalian", "ANorcadian",
      "ANsardinian", "ANtuscan"
    )
  )

  chr22$AF_EAS <- weighted_af(
    chr22,
    af_columns = c(
      "AFcdx", "AFchb", "AFchs", "AFjpt", "AFkhv",
      "AFcambodian", "AFdai", "AFdaur", "AFhan",
      "AFhezhen", "AFjapanese", "AFlahu", "AFmiaozu",
      "AFmongola", "AFnaxi", "AForoqen", "AFshe",
      "AFtu", "AFtujia", "AFxibo", "AFyizu"
    ),
    an_columns = c(
      "ANcdx", "ANchb", "ANchs", "ANjpt", "ANkhv",
      "ANcambodian", "ANdai", "ANdaur", "ANhan",
      "ANhezhen", "ANjapanese", "ANlahu", "ANmiaozu",
      "ANmongola", "ANnaxi", "ANoroqen", "ANshe",
      "ANtu", "ANtujia", "ANxibo", "ANyizu"
    )
  )

  chr22$AF_SAS <- weighted_af(
    chr22,
    af_columns = c("AFbeb", "AFgih", "AFitu", "AFpjl", "AFstu"),
    an_columns = c("ANbeb", "ANgih", "ANitu", "ANpjl", "ANstu")
  )

  chr22$AF_IAM <- weighted_af(
    chr22,
    af_columns = c("AFcolombian", "AFkaritiana", "AFmaya", "AFpima"),
    an_columns = c("ANcolombian", "ANkaritiana", "ANmaya", "ANpima")
  )

  chr22 <- chr22[
    (
      (chr22$AF_AFR > 0.01 & chr22$AF_AFR < 0.99) |
      (chr22$AF_EAS > 0.01 & chr22$AF_EAS < 0.99) |
      (chr22$AF_EUR > 0.01 & chr22$AF_EUR < 0.99) |
      (chr22$AF_IAM > 0.01 & chr22$AF_IAM < 0.99) |
      (chr22$AF_SAS > 0.01 & chr22$AF_SAS < 0.99)
    ),
  ]

  chr22 <- na.omit(chr22)

  chr22 %>%
    dplyr::select(
      POS,
      REF,
      ALT,
      AF_AFR,
      AF_EAS,
      AF_EUR,
      AF_IAM,
      AF_SAS
    ) %>%
    dplyr::rename(
      rsID = POS,
      AFR = AF_AFR,
      EAS = AF_EAS,
      EUR = AF_EUR,
      IAM = AF_IAM,
      SAS = AF_SAS
    )
}
