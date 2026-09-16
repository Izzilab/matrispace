#' @title matrispace: Identification and Visualization of Spatially-Resolved ECM Gene Expression Patterns
#'
#' @description
#' \pkg{matrispace} is the R package companion to MatriSpace. It brings the
#' same extracellular matrix-focused workflow used in the web application to
#' local, scriptable analysis of spatial transcriptomics data. The package
#' includes functions for:
#'
#' \itemize{
#'   \item Matrisome scoring using UCell signatures
#'   \item ECM niche classification using the ScType algorithm
#'   \item Spatial statistics (LISA, Moran's I)
#'   \item Spatial co-expression analysis of unique heterotypic matrisome pairs
#'   \item Static visualizations for spatial data
#' }
#'
#' @section Main Analysis Functions:
#' \describe{
#'   \item{\code{\link{prepare_object}}}{Master preprocessing pipeline}
#'   \item{\code{\link{score_matrisome}}}{Calculate matrisome category scores}
#'   \item{\code{\link{annotate_ecm_niches}}}{Classify ECM niches with ScType}
#'   \item{\code{\link{compute_lisa}}}{Local spatial association statistics}
#'   \item{\code{\link{compute_morans_i}}}{Global spatial autocorrelation}
#'   \item{\code{\link{score_matrisome_pair_coexpression}}}{Matrisome-pair co-expression scoring}
#' }
#'
#' @section Data Objects:
#' \describe{
#'   \item{\code{\link{matrisome_db}}}{Gene-to-category mapping database}
#'   \item{\code{\link{matrisome_signatures}}}{UCell gene signatures}
#'   \item{\code{\link{ecm_signatures}}}{ECM niche signatures}
#'   \item{\code{\link{matrisome_pairs}}}{Unique heterotypic MatriComDB pairs}
#' }
#'
#' @section Visualization Functions:
#' \describe{
#'   \item{\code{\link{plot_spatial_feature}}}{Single feature spatial plot}
#'   \item{\code{\link{plot_matrisome}}}{Matrisome category plots}
#'   \item{\code{\link{plot_lisa}}}{LISA clustering visualization}
#'   \item{\code{\link{plot_spatial_blend}}}{Two-feature blend plot}
#' }
#'
#' @author Ayomide Oshinjo, Daiqing Chen, Petar Petrov, Valerio Izzi, Alexandra Naba
#'
#' @references
#' \url{https://github.com/Theayomideo/matrispace}
#'
#' @seealso
#' The companion MatriSpace web application:
#' \url{https://izzilab-matrispace.share.connect.posit.cloud/}
#'
#' @docType package
#' @name matrispace-package
#' @aliases matrispace
#'
#' @importFrom methods as new
#' @importFrom stats dist IQR median na.omit reorder sd setNames var
#' @importFrom dplyr `%>%` bind_rows count filter group_by left_join mutate n summarise ungroup
#' @importFrom ggplot2 aes annotate coord_cartesian coord_fixed element_blank
#'   element_text facet_wrap geom_col geom_point ggplot guide_legend guides
#'   labs margin scale_color_gradient scale_color_gradientn scale_color_identity
#'   scale_color_manual scale_fill_gradient scale_fill_gradientn scale_fill_manual
#'   scale_x_continuous scale_y_discrete scale_y_reverse theme theme_light theme_void
#' @importFrom Matrix colSums rowMeans rowSums
#' @importFrom Seurat CreateAssayObject DefaultAssay FetchData GetAssayData
#'   GetTissueCoordinates SCTransform SetAssayData SpatialFeaturePlot
#' @importFrom SeuratObject Assays
#' @importFrom UCell AddModuleScore_UCell
#' @importFrom HGNChelper checkGeneSymbols
#' @importFrom scales hue_pal percent_format rescale
#' @importFrom patchwork wrap_plots
#' @importFrom RColorBrewer brewer.pal
"_PACKAGE"
