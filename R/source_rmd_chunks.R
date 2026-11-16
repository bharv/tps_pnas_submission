# =============================================================================
# source_rmd_chunks.R
# -----------------------------------------------------------------------------
# Evaluate selected, LABELLED code chunks of an R Markdown file, without knitting
# it. The appendix scripts use this so that model runners, priors, helpers and
# figure functions live in exactly one place (the main-text notebook) and the
# appendix re-uses them verbatim instead of carrying a drifting copy.
#
#   source_rmd_chunks("01_tps_catch_models/tps_catch_wq_models.Rmd",
#                     labels = c("libraries", "config", "load-data"))
#
# Chunks are matched on the label that follows ```{r  (e.g. ```{r config}).
# Chunks are evaluated in the order they appear in the file, not the order of
# `labels`. Missing labels are an error, so a renamed chunk fails loudly.
# =============================================================================

rmd_chunks <- function(rmd) {
  lines  <- readLines(rmd, warn = FALSE)
  starts <- grep("^```\\{r", lines)
  ends   <- grep("^```\\s*$", lines)
  out    <- list()
  for (s in starts) {
    e     <- ends[ends > s][1]
    label <- sub("^```\\{r[ ,]*([^,} ]*).*$", "\\1", lines[s])
    if (!nzchar(label)) label <- sprintf("unlabelled-line-%d", s)
    out[[label]] <- if (e > s + 1) lines[(s + 1):(e - 1)] else character(0)
  }
  out
}

source_rmd_chunks <- function(rmd, labels, envir = parent.frame(), echo = FALSE) {
  ch <- rmd_chunks(rmd)
  missing <- setdiff(labels, names(ch))
  if (length(missing)) stop("chunk label(s) not found in ", rmd, ": ",
                            paste(missing, collapse = ", "))
  for (lab in names(ch)[names(ch) %in% labels]) {
    if (echo) message("-- sourcing chunk '", lab, "' from ", basename(rmd))
    eval(parse(text = ch[[lab]], keep.source = FALSE), envir = envir)
  }
  invisible(TRUE)
}
