#' @title Package Data Documentation
#' @name data
#' @description Documentation for package data objects.
NULL

#' Matrisome Gene Database
#'
#' A data frame mapping genes to matrisome categories and subcategories.
#'
#' @format A data frame with columns:
#' \describe{
#'   \item{division}{Main matrisome category (Core matrisome or Matrisome-associated)}
#'   \item{notes}{Additional notes}
#'   \item{gene}{HGNC gene symbol}
#'   \item{Gene Name}{Full gene name}
#'   \item{Synonyms}{Alternative gene names (pipe-separated)}
#'   \item{HGNC_IDs}{HGNC identifier}
#'   \item{HGNC_IDs Links}{Links to HGNC database}
#'   \item{UniProt_IDs}{UniProt identifiers}
#'   \item{Refseq_IDs}{RefSeq identifiers}
#'   \item{Notes_ignore}{Internal notes (ignored)}
#'   \item{ecm_subcategory}{Subcategory (e.g., Collagens, ECM Glycoproteins)}
#' }
#'
#' @source The Matrisome Project \url{https://matrisomeproject.org}
#'
#' @examples
#' data(matrisome_db)
#' head(matrisome_db)
#'
"matrisome_db"

#' Matrisome UCell Signatures
#'
#' Named list of gene vectors for MatriSpace matrisome categories and subcategories,
#' formatted for use with UCell scoring.
#'
#' @format A named list where:
#' \describe{
#'   \item{names}{Signature names (e.g., "collagens", "ecm_glycoproteins")}
#'   \item{values}{Character vectors of gene symbols}
#' }
#'
#' @details
#' Includes signatures for:
#' \itemize{
#'   \item 6 main categories: Collagens, ECM Glycoproteins, Proteoglycans,
#'         ECM Regulators, Secreted Factors, ECM-affiliated Proteins
#'   \item curated subcategories and gene families such as Basement Membrane,
#'         Laminins, Matricellular proteins, Mucins, and Peri-vascular ECM
#' }
#'
#' @examples
#' data(matrisome_signatures)
#' names(matrisome_signatures)
#'
"matrisome_signatures"

#' ECM Niche Signatures
#'
#' Gene signatures for three ECM niche types used for UCell scoring.
#'
#' @format A named list with three elements:
#' \describe{
#'   \item{Vascular}{Genes associated with vascular ECM}
#'   \item{Basement}{Genes associated with basement membrane}
#'   \item{Interstitial}{Genes associated with interstitial matrix}
#' }
#'
#' @examples
#' data(ecm_signatures)
#' names(ecm_signatures)
#'
"ecm_signatures"

#' MatriComDB Matrisome Pair Database
#'
#' Database of unique heterotypic MatriComDB pairs for spatial co-expression
#' analysis. Every pair contains at least one matrisome gene. The source data
#' include ECM ligand-receptor relationships and structural ECM interactions;
#' a row should not be interpreted as signaling unless its annotation supports
#' that interpretation.
#'
#' @format A data frame with columns:
#' \describe{
#'   \item{Ligand}{First pair-member gene symbol (legacy column name)}
#'   \item{Receptor}{Second pair-member gene symbol (legacy column name)}
#'   \item{Source}{Database source for the interaction}
#'   \item{Interaction_Type}{Type of interaction}
#'   \item{Database_Score}{Confidence score from source database}
#' }
#'
#' @details
#' Reciprocal records are collapsed to one alphabetically ordered pair and
#' homomeric records are excluded.
#'
#' @examples
#' data(lr_database)
#' head(lr_database)
#' nrow(lr_database)
#'
"lr_database"

#' ECM Niche Markers
#'
#' Marker genes for ECM niche classification using ScType.
#'
#' @format A list with two elements:
#' \describe{
#'   \item{gs_positive}{Positive marker genes per niche}
#'   \item{gs_negative}{Negative marker genes per niche}
#' }
#'
#' @examples
#' data(ecm_markers)
#' names(ecm_markers$gs_positive)
#'
"ecm_markers"
