library(purrr)
library(splines2)
################################################################################
################################################################################
surv_mspline_rw <- function(times,X,betas, gammas,degree,knots,bknots) {
    covxs <- as.vector(X %*% betas)
    i_spline_basis_evals <- iSpline(times, knots=knots, degree=degree,Boundary.knots=bknots,
                                intercept=FALSE)
    prods <- as.vector(i_spline_basis_evals %*% as.matrix(gammas))
    do.call(rbind,map(covxs, ~exp(-prods*exp(.))))
}
################################################################################
surv_mspline_simplex <- function(times,X,betas,intercept, gammas,degree,knots,bknots) {
    covxs <- as.vector(X %*% betas + intercept)
    i_spline_basis_evals <- iSpline(times, knots=knots, degree=degree,Boundary.knots=bknots,
                                intercept=FALSE)
    prods <- as.vector(i_spline_basis_evals %*% as.matrix(gammas))
    do.call(rbind,map(covxs, ~exp(-prods*exp(.))))
}
################################################################################
################################################################################
get_briers_exp_ <- function(i,j,post) {
    survs <- surv_exp(unique_times, X, 
    post[i,j,sprintf("betas[%d]", ncol(X))],
    post[i,j,"intercept"])
    pec(object=list("model"=survs),
    formula=Surv(time, event)~1, 
    data=df, exact=TRUE, cens.model="marginal", splitMethod="none",B=0,verbose=FALSE)$AppErr$model
}
################################################################################
get_briers_weibull_ <- function(i,j,post) {
    survs <- surv_weibull(unique_times, X, 
    post[i,j,sprintf("betas[%d]", ncol(X))],
    post[i,j,"intercept"],
    post[i,j,"alpha"]
    )
    pec(object=list("model"=survs),
    formula=Surv(time, event)~1, 
    data=df, exact=TRUE, cens.model="marginal", splitMethod="none",B=0,verbose=FALSE)$AppErr$model
}
################################################################################
get_briers_gamma_ <- function(i,j,post) {
    survs <- surv_gamma(unique_times, X, 
    post[i,j,sprintf("betas[%d]", ncol(X))],
    post[i,j,"intercept"],
    post[i,j,"alpha"]
    )
    pec(object=list("model"=survs),
    formula=Surv(time, event)~1, 
    data=df, exact=TRUE, cens.model="marginal", splitMethod="none",B=0,verbose=FALSE)$AppErr$model
}
################################################################################
get_briers_mspline_rw_ <- function(i,j,post) {
    survs <- surv_mspline_rw(unique_times, X, 
    post[i,j,sprintf("betas[%d]", ncol(X))],
    post[i,j,sprintf("gammas[%d]", 1:nbasis)], 
    mspline_degree,
    knots,attr(i_spline_basis_evals, "Boundary.knots"))
    pec(object=list("model"=survs),
    formula=Surv(time, event)~1, 
    data=df, exact=TRUE, cens.model="marginal", splitMethod="none",B=0,verbose=FALSE)$AppErr$model
}
################################################################################
get_briers_mspline_simplex_ <- function(i,j,post) {
    survs <- surv_mspline_simplex(unique_times, X, 
    post[i,j,sprintf("betas[%d]", ncol(X))],
    post[i,j,"intercept"],
    post[i,j,sprintf("gammas[%d]", 1:nbasis)], 
    mspline_degree,
    knots,attr(i_spline_basis_evals, "Boundary.knots"))
    pec(object=list("model"=survs),
    formula=Surv(time, event)~1, 
    data=df, exact=TRUE, cens.model="marginal", splitMethod="none",B=0,verbose=FALSE)$AppErr$model
}
################################################################################
get_brier_score_df <- function(post, get_briers) {
    nsamps_per_chain <- dim(post)[1]
    nchains <- dim(post)[2]
    purrr::map2_dfr(rep(1:nsamps_per_chain, nchains),rep(1:nchains, each=nsamps_per_chain),
                    ~tibble(brierscore=get_briers(.x,.y,post),time=unique_times)) %>%
        group_by(time) %>%
        summarise(mean=mean(brierscore) , low=quantile(brierscore, probs = c(0.05)), up=quantile(brierscore, probs = c(0.975)))
}
################################################################################
################################################################################
target_func <- function(t, u, gammas, degree, knots, bknots) {
    #abs(u - as.vector(iSpline(x, knots = knots, degree = degree, intercept = FALSE, Boundary.knots=bknots)[1,]))
    u - t(as.matrix(iSpline(t, knots = knots, degree = degree, intercept = FALSE, Boundary.knots=bknots)[1,]))%*% as.matrix(gammas)
}
################################################################################
samp_mspline_simplex_ <- function(x, betas, intercept, gammas, degree, knots, bknots, tmax) {
    H0_max <- t(as.matrix(iSpline(tmax, knots = knots, degree = degree, intercept = FALSE, Boundary.knots=bknots)[1,]))%*% as.matrix(gammas) 
    covx <- t(as.matrix(x)) %*% as.matrix(betas) + intercept
    rrx <- exp(covx)
    U <- rexp(1,rate=rrx)
    iteration_count <- 1
    while(U > H0_max) { # H0 is non-decreasing in time and we do ppc on study-length (conditioning)
        if(iteration_count > 1000) print("WARNING MAX ITER")
        U <- rexp(1,rate=rrx)
        iteration_count <- iteration_count+1
    }
    uniroot(target_func, c(0,tmax), u=U, gammas=gammas, degree=degree, knots=knots, bknots=bknots)$root 
}
################################################################################
samp_mspline_simplex <- function(X, betas, intercept, gammas, degree, knots, bknots, tmax) {
    purrr::map_dbl(1:nrow(X), ~samp_mspline_simplex_(X[.,], betas, intercept, gammas, degree, knots, bknots, tmax))
}
################################################################################
pp_rep_samps_mspline_simplex <- function(nsamps, X, post, degree, knots, bknots,nbasis, tmax) {
    nsamps_per_chain <- dim(post)[1]
    nchains <- dim(post)[2]
    NC <- ncol(X)
    xindices <- rep(1:nsamps_per_chain, nchains)
    yindices <- rep(1:nchains, each=nsamps_per_chain)
    indices  <- sample(1:length(xindices),nsamps)
    xindices <- xindices[indices]
    yindices <- yindices[indices]
   
    samp_times_lst <- purrr::map2(xindices, yindices, ~samp_mspline_simplex(X, post[.x,.y,sprintf("betas[%d]", 1:NC)], post[.x,.y,"intercept"],
                                                post[.x,.y,sprintf("gammas[%d]", 1:nbasis)], 
                                                mspline_degree,
                                                knots,attr(i_spline_basis_evals, "Boundary.knots"), time_max))
    do.call(rbind,samp_times_lst)
}
################################################################################
################################################################################
surv_exp <- function(times,X,betas,intercept) {
    covxs <- exp(as.vector(X %*% betas + intercept))
    do.call(rbind,map(covxs, ~pexp(times, rate=., lower.tail=FALSE) ))
}
################################################################################
surv_weibull <- function(times,X,betas,intercept, shape) {
    covxs <- exp(as.vector(X %*% betas + intercept))
    do.call(rbind,map(covxs, ~pweibull(times, scale=.,shape=shape, lower.tail=FALSE) ))
}
################################################################################
surv_gamma <- function(times,X,betas,intercept, shape) {
    covxs <- exp(as.vector(X %*% betas + intercept))
    do.call(rbind,map(covxs, ~pgamma(times, rate=.,shape=shape, lower.tail=FALSE) ))
}
################################################################################
################################################################################
# code below taken from 
# https://raw.githubusercontent.com/cran/RcmdrPlugin.KMggplot2/f38dcfb4f5ea3137e3c1cc3c96b4a1961461765e/R/geom-stepribbon.r
################################################################################
#' Step ribbon plots.
#'
#' \code{geom_stepribbon} is an extension of the \code{geom_ribbon}, and
#' is optimized for Kaplan-Meier plots with pointwise confidence intervals
#' or a confidence band.
#'
#' @section Aesthetics:
#' \Sexpr[results=rd,stage=build]{ggplot2:::rd_aesthetics("geom", "ribbon")}
#'
#' @seealso
#'   \code{\link[ggplot2:geom_ribbon]{geom_ribbon}} \code{geom_stepribbon}
#'   inherits from \code{geom_ribbon}.
#' @inheritParams ggplot2:::geom_ribbon
#' @param kmplot If \code{TRUE}, missing values are replaced by the previous
#' values. This option is needed to make Kaplan-Meier plots if the last
#' observation has event, in which case the upper and lower values of the
#' last observation are missing. This processing is optimized for results
#' from the survfit function.
#' @examples
#' huron <- data.frame(year = 1875:1972, level = as.vector(LakeHuron))
#' h <- ggplot(huron, aes(year))
#' h + geom_stepribbon(aes(ymin = level - 1, ymax = level + 1), fill = "grey70") +
#'     geom_step(aes(y = level))
#' h + geom_ribbon(aes(ymin = level - 1, ymax = level + 1), fill = "grey70") +
#'     geom_line(aes(y = level))
#' @rdname geom_stepribbon
#' @importFrom ggplot2 layer GeomRibbon
#' @export
geom_stepribbon <- function(
  mapping = NULL, data = NULL, stat = "identity", position = "identity",
  na.rm = FALSE, show.legend = NA, inherit.aes = TRUE, kmplot = FALSE, ...) {
  layer(
    data = data,
    mapping = mapping,
    stat = stat,
    geom = GeomStepribbon,
    position = position,
    show.legend = show.legend,
    inherit.aes = inherit.aes,
    params = list(
      na.rm = na.rm,
      kmplot = kmplot,
      ...
    )
  )
}

#' @rdname geom_stepribbon
#' @format NULL
#' @usage NULL
#' @export
GeomStepribbon <- ggproto(
  "GeomStepribbon", GeomRibbon, 
  
  extra_params = c("na.rm", "kmplot"),
  
  draw_group = function(data, panel_scales, coord, na.rm = FALSE) {
    if (na.rm) data <- data[complete.cases(data[c("x", "ymin", "ymax")]), ]
    data <- rbind(data, data)
    data <- data[order(data$x), ]
    data$x <- c(data$x[2:nrow(data)], NA)
    data <- data[complete.cases(data["x"]), ]
    GeomRibbon$draw_group(data, panel_scales, coord, na.rm = FALSE)
  },
  
  setup_data = function(data, params) {
    if (params$kmplot) {
      data <- data[order(data$PANEL, data$group, data$x), ]
      tmpmin <- tmpmax <- NA
      for (i in 1:nrow(data)) {
        if (is.na(data$ymin[i])) {
          data$ymin[i] <- tmpmin
        }
        if (is.na(data$ymax[i])) {
          data$ymax[i] <- tmpmax
        }
        tmpmin <- data$ymin[i]
        tmpmax <- data$ymax[i]
      }
    }
    data
  }
  
)
