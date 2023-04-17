apsrtable_new <- function (..., se = c("robust", "vcov", "both", "pval"), model.names = NULL, 
          model.counter = 1, digits = 2, stars = 1, lev = 0.05, align = c("left", 
                                                                          "center", "right"), order = c("lr", "rl", "longest"), 
          notes = list(se.note, stars.note), omitcoef = NULL, coef.names = NULL, 
          coef.rows = 2, multicolumn.align = c("center", "left", "right"), 
          col.hspace = NULL, Sweave = FALSE, float = "table", Minionfig = FALSE, 
          label = NULL, caption = NULL, caption.position = c("above", 
                                                             "below")) 
{
  x <- list()
  myenv <- new.env()
  signif.stars <- TRUE
  order <- match.arg(order, c("lr", "rl", "longest"))
  opts <- match.call(expand.dots = FALSE)
  se <- match.arg(se, c("robust", "vcov", "both", "pval"))
  align <- substr(align, 1, 1)
  align <- match.arg(align, c("l", "c", "r"))
  multicolumn.align <- match.arg(substr(multicolumn.align, 
                                        1, 1), c("c", "l", "r"))
  adigits <- ifelse(align == "c", -1, digits)
  caption.position <- match.arg(substr(caption.position, 1, 
                                       1), c("a", "b"))
  models <- list(...)
  nmodels <- length(models)
  coef.cols <- ifelse(coef.rows == 2, 1, 2)
  if (is.null(col.hspace)) {
    col.hspace <- ifelse(coef.cols == 1, "", "2em")
  }
  if (col.hspace != "") {
    col.hspace <- paste("@{\\hspace{", col.hspace, "}}", 
                        sep = "")
  }
  colspec <- paste("{", align, paste(rep(paste("D{.}{.}{", 
                                               rep(adigits, coef.cols), "}", sep = "", collapse = ""), 
                                         nmodels), collapse = col.hspace), "}", collapse = "")
  if (float == "longtable") {
    long <- TRUE
    floatspec <- paste("\\begin{", float, "}", colspec, 
                       "\n", ifelse(caption.position == "a", paste("\\caption{", 
                                                                   caption, "}\n\\label{", label, "}", sep = ""), 
                                    ""), sep = "")
  }
  else {
    long <- FALSE
    floatspec <- paste(ifelse(!Sweave, paste("\\begin{", 
                                             float, "}[!ht]\n", ifelse(caption.position == "a", 
                                                                       paste("\\caption{", caption, "}\n\\label{", 
                                                                             label, "}", sep = ""), ""), sep = ""), ""), 
                       paste("\n\\begin{tabular}", colspec, sep = ""))
  }
  x <- paste(floatspec, ifelse(long, "\\\\", ""), sep = "")
  if (Minionfig) {
    x <- c(x, "%Uncomment the following line and the end one to change figure versions\n%if you are using a full-featured family such as Minion Pro.\n\\figureversion{tabular}\n")
  }
  model.summaries <- lapply(models, function(x) {
    s <- try(apsrtableSummary(x), silent = TRUE)
    if (inherits(s, "try-error")) {
      s <- summary(x)
    }
    if (!is.null(x$se) && se != "vcov") {
      est <- coef(x)
      if (inherits(x$se, "matrix")) {
        x$se <- sqrt(diag(x$se))
      }
      s$coefficients[, 3] <- tval <- est/x$se
      e <- try(s$coefficients[, 4] <- 2 * pt(abs(tval), 
                                             length(x$residuals) - x$rank, lower.tail = FALSE), 
               silent = TRUE)
      if (inherits(e, "try-error")) {
        s$coefficients[, 4] <- 2 * pnorm(abs(tval), 
                                         lower.tail = FALSE)
      }
      s$se <- x$se
    }
    if (se == "pval") {
      s$coefficients[, 2] <- s$coefficients[, 4]
    }
    return(s)
  })
  if (se == "pval") {
    se.note <- pval.note
  }
  if (is.null(model.names)) {
    m.first = model.counter
    m.last = m.first + (nmodels - 1)
    model.names = paste("Model", m.first:m.last)
  }
  else if (!is.null(model.names) && (length(model.names) < 
                                     nmodels)) {
    m.first = length(model.names) + 1
    model.names = c(model.names, paste("Model", m.first:nmodels))
  }
  coefnames <- orderCoef(model.summaries, order = order)
  incl <- rep(TRUE, length(coefnames))
  names(incl) <- coefnames
  if (!is.null(omitcoef)) {
    omitcoef <- unlist(sapply(omitcoef, eval, envir = myenv))
    incl[omitcoef] <- FALSE
  }
  model.summaries <- coefPosition(model.summaries, coefnames)
  if (!is.null(coef.names)) {
    if (length(coef.names) != sum(incl)) {
      warning("Supplied coef.names not the same length as output. Check automatic names before supplying 'pretty' names.\n")
    }
    coefnames[incl] <- coef.names
  }
  else {
    coefnames[incl] <- sanitize(coefnames[incl])
  }
  out.table <- lapply(model.summaries, function(x) {
    var.pos <- attr(x, "var.pos")
    model.out <- model.se.out <- star.out <- rep(NA, length(coefnames))
    model.out[var.pos] <- x$coefficients[, 1]
    star.out[var.pos] <- apsrStars(x$coefficients, stars = stars, 
                                   lev = lev, signif.stars = TRUE)
    model.out <- ifelse(!is.na(model.out), paste(formatC(model.out, 
                                                         digits = digits, format = "f"), star.out), "")
    model.se.out[var.pos] <- x$coefficients[, 2]
    if (!is.null(x$se) & se %in% c("robust", "both")) {
      model.se.out[var.pos] <- x$se
    }
    model.se.out <- ifelse(!is.na(model.se.out), paste("(", 
                                                       formatC(model.se.out, digits = digits, format = "f"), 
                                                       ")", sep = ""), "")
    if (se == "both" && !is.null(x$se)) {
      model.se.out[var.pos] <- ifelse(model.se.out != 
                                        "", paste(model.se.out, " [", formatC(x$coefficients[, 
                                                                                             2], digits = digits, format = "f"), "]", sep = ""), 
                                      "")
    }
    if (coef.rows == 2) {
      model.out <- rep(model.out[incl], each = 2)
      model.se.out <- rep(model.se.out[incl], each = 2)
      pos.se <- (1:length(model.out))[(1:length(model.out)%%2 == 
                                         0)]
      model.out[pos.se] <- model.se.out[pos.se]
    }
    else {
      model.out <- model.out[incl]
      model.out <- cbind(model.out, model.se.out[incl])
    }
    attr(model.out, "model.info") <- modelInfo(x)
    return(model.out)
  })
  out.matrix <- matrix(unlist(out.table), length(coefnames[incl]) * 
                         coef.rows, nmodels * coef.cols)
  out.matrix <- cbind(rep(coefnames[incl], each = coef.rows), 
                      out.matrix)
  if (coef.rows == 2) {
    out.matrix[(row(out.matrix)[, 1]%%2 == 0), 1] <- ""
  }
  out.info <- lapply(out.table, attr, "model.info")
  info.names <- orderCoef(out.info)
  out.info <- coefPosition(out.info, orderCoef(out.info))
  out.info <- lapply(out.info, function(x) {
    var.pos <- attr(x, "var.pos")
    model.out <- rep("", length(info.names))
    model.out[var.pos] <- coef(x)
    return(model.out)
  })
  out.info <- matrix(unlist(out.info), length(info.names), 
                     nmodels)
  out.info <- cbind(as.character(info.names), out.info)
  if (coef.rows == 2) {
    out.matrix <- rbind(c("%", model.names), out.matrix)
  }
  outrows <- nrow(out.matrix)
  if (coef.cols == 1) {
    out.matrix <- rbind(out.matrix, out.info)
    out.matrix[, -1] <- format(out.matrix[, -1])
    out.matrix[, 1] <- format(out.matrix)[, 1]
    out.matrix <- apply(out.matrix, 1, paste, collapse = " & ")
    out.info <- out.matrix[(1 + outrows):length(out.matrix)]
    out.matrix <- out.matrix[1:outrows]
  }
  else {
    out.matrix <- format(out.matrix)
    out.matrix <- apply(out.matrix, 1, paste, collapse = " & ")
    out.info[, -1] <- format(out.info[, -1])
    out.info[, -1] <- sapply(as.matrix(out.info[, -1]), 
                             function(x) {
                               paste("\\multicolumn{", coef.cols, "}{", multicolumn.align, 
                                     "}{", x, "}", sep = "")
                             })
    out.info[, 1] <- format(out.info[, 1])
    out.info <- apply(out.info, 1, paste, collapse = " & ")
  }
  headrow <- paste("\n\\hline \n", paste(" &", paste("\\multicolumn{", 
                                                     coef.cols, "}{", multicolumn.align, "}{", model.names, 
                                                     "}", collapse = " & ")), "\\\\ \\hline\n")
  if (long) {
    headrow <- paste(headrow, "\\endhead\n", sep = "")
  }
  x <- c(x, headrow)
  x <- c(x, paste(out.matrix, collapse = "\\\\ \n"))
  x <- c(x, "\\\\\n")
  x <- c(x, paste(out.info, collapse = "\\\\ \n"))
  se <- ifelse((se != "vcov" && sum(unlist(lapply(model.summaries, 
                                                  function(x) !is.null(x$se))) > 0)), "robust", "vcov")
  thenotes <- as.list(1:length(notes))
  thenotes[!sapply(notes, is.function)] <- notes[!sapply(notes, 
                                                         is.function)]
  thenotes[sapply(notes, is.function)] <- lapply(notes[sapply(notes, 
                                                              is.function)], do.call, args = list(env = myenv))
  x <- c(x, "\\\\ \\hline\n")
  notes <- lapply(thenotes, function(x) {
    paste("\\multicolumn{", (nmodels * coef.cols) + 1, "}{l}{\\footnotesize{", 
          x, "}}", sep = "")
  })
  x <- c(x, paste(notes, collapse = "\\\\\n"))
  if (!long) {
    x <- c(x, "\n\\end{tabular}")
  }
  if (long) {
    x <- c(x, "\n\\end{longtable}")
  }
  if (caption.position == "b") {
    x <- c(x, paste("\n\\caption{", caption, "}\n\\label{", 
                    label, "}", sep = ""))
  }
  x <- c(x, "\n")
  if (Minionfig) {
    x <- c(x, "\n\\figureversion{proportional}\n")
  }
  if (!Sweave & !long) {
    x <- c(x, paste("\\end{", float, "}\n", sep = ""))
  }
  class(x) <- "apsrtable"
  return(x)
}
