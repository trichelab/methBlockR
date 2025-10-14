#' get and plot weights for one or more factors (combine with enrichment)
#' 
#' @param x     RcppML::nmf object or SingleCellExperiment-like holder of one
#' @param j     which factor (1)
#' @param k     how many top weights to plot (30)
#' @param what  what to plot (if anything) ("logdensitymix")
#'
#' @return      a patchworked ggplot object and/or the data for the plot
#'
#' @import      patchwork
#' @import      ggplot2 
#' @import      RcppML
#' @import      mixR
#'
#' @export
#'
nmfWeights <- function(SCE, j=1, k=30, what=c("logdensitymix", "logdensity", "mix")) {

  if (!is(SCE, "nmf")) { 
    # could look for a LinearEmbeddingMatrix, or 
    stopifnot("nmf_fit" %in% names(metadata(SCE)))
    nmf <- metadata(SCE)$nmf_fit
  } else { 
    nmf <- SCE 
  }
  stopifnot(j <= ncol(nmf@w))

  what <- match.arg(what)
  if (what == "logdensitymix") { 
    ldm <- .logdensmix(nmf@w[, j])
    fit <- attr(ldm, "fit")
    toplot <- ldm[rev(order(ldm$weight))[seq_len(k)], ]
    ymin <- min(toplot$weight)
    ymax <- max(toplot$weight) 
    p1 <- ggplot(toplot) +  
            aes(x=reorder(name, raw, decreasing=TRUE), 
                y=raw, fill=weight) + 
            geom_col() + 
            xlab("Feature") + 
            ylab("Weight") + 
            theme_classic() + 
            theme(axis.text.x=element_text(angle=45, hjust=1),
                  legend.position="none") + 
            ggtitle(paste("Factor", j, "(top", k, "feature weights)"))
    p2 <- plot(fit, xlab="Normalized log-weight")
    p1 + p2
  } else {
    stop("Not supported yet") 
  }

}


# helper fn
.gt0 <- function(x) x[x > 0]


# helper fn
.logdensmix <- function(xx) {
  message("Fitting mixture model...")
  x <- log(.gt0(xx) + .Machine$double.eps)
  x <- x - min(x)
  x <- x / max(x)
  names(x) <- names(.gt0(xx))
  sel <- select(x=x, ncomp=1:3, tol=1e-3)
  best <- which(sel$best == "*")
  ncomp <- sel$ncomp[best]
  ev <- sel$equal.var[best] == "Y"
  fit <- mixfit(x=x, ncomp=ncomp, ev=ev, tol=1e-3)
  res <- data.frame(raw=.gt0(xx),
                    name=names(fit$data),
                    weight=fit$data,
                    highlight=(apply(fit$comp.prob, 1, which.max) == ncomp))
  attr(res, "fit") <- fit
  attr(res, "ncomp") <- ncomp 
  return(res)
}

