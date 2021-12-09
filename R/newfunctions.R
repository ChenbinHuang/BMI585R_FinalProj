#' unscale the object return by scale function
#' @export
#'
#' @param x scaled object
#' @return the un-scaled matrix
#' @examples
#' mtcars |> scale() |> unscale()
unscale = function(x){
  centers = attr(x,"scaled:center")
  scales = attr(x,"scaled:scale")
  t(apply(x, 1, function(xx){xx*scales+centers}))
}

#' estimate the data with first npc pricipal components
#' @export
#'
#' @param x original numeric data
#' @param npc number of pricipal components
#' @return the estimation matrix
#' @examples
#' pcApprox(mtcars, 7)
pcApprox = function(x, npc){
  pca = stats::prcomp(x,center=T,scale.=T)
  recons = pca$x[,1:npc] %*% t(pca$rotation[,1:npc])
  t(apply(recons, 1, function(xx){xx*pca$scale + pca$center}))
}

#' Create a lollipop plot of the principal component loadings of the (potentially unscaled/uncentered) data x
#' @export
#'
#' @param x original numeric data
#' @return the plot object
#' @examples
#' pcLollipop(iris[,1:4])
#' pcLollipop(mtcars)
pcLollipop=function(x){
  pca = stats::prcomp(x,center=T,scale.=T)
  pc_long = pca$rotation |> data.frame() |> dplyr::mutate(variable=rownames(pca$rotation, do.NULL = F)) |> tidyr::pivot_longer( cols=!variable, values_to="value", names_to="PC")
  ggplot2::ggplot(data =pc_long, mapping = ggplot2::aes(variable,value)) +
    ggplot2::geom_segment( ggplot2::aes(x=variable, xend=variable, y=0, yend=value, color=PC)) +
    ggplot2::geom_point(ggplot2::aes(color=PC)) +
    ggplot2::facet_wrap(~PC, ncol = 1) +
    ggplot2::theme_minimal() + 
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5, hjust=1))
}
