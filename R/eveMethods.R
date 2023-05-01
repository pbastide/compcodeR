#' @title Get the tree from a phyloCompData object
#'
#' @description
#' Return the tree of a \code{phyloCompData} object,
#' removing replicates for use in EVE.
#' If no tree, return a star tree with unit height, and throw a warning.
#' 
#' @param cdata a phyloCompData object.
#' 
#' @return A tree of class \code{phylo}
#' 
#' @keywords internal
#' 
getTreeEVE <- function(cdata) {
  if (is.null(phylo.tree(cdata)) || length(phylo.tree(cdata)) == 0) {
    message("There were no tree in the data object. Using a star tree of unit height in EVE")
    ntaxa <- info.parameters(cdata)$n.samples
    tree <- ape::stree(ntaxa, "star")
    tree$edge.length <- rep(1, nrow(tree$edge))
    tree$tip.label <- rownames(sample.annotations(cdata))
    return(tree)
  } else {
    tree <- phylo.tree(cdata)
    sample.annotations <- sample.annotations(cdata)
    sample.annotations$id <- colnames(count.matrix(cdata))
    tree <- removeReplicatesFromTree(tree, sample.annotations, species = "id.species", id = "id")
    if (!ape::is.ultrametric(tree)) stop("The tree must be ultrametric.")
    return(tree)
  }
}

#' @title Add replicates to a tree
#'
#' @description
#' Utility function to add replicates to a tree, as tips with zero length branches.
#'
#' @param tree A phylogenetic tree with n tips.
#' @param traits A data frame containing at least two columns,
#' one with sample ids, and on with species names for each samples.
#' @param species Name of the column containing species names. Default to "species".
#' @param id Name of the column containing samples ids. Default to "id".
#'
#' @return A phylogenetic tree with as many tips as the number of rows in \code{traits},
#'  and clusters of tips with zero branch lengths corresponding to replicates.
#'
#' @keywords internal
#'
removeReplicatesFromTree <- function(tree, traits,
                                     species = "species", id = "id") {
  if (!requireNamespace("phytools", quietly = TRUE)) {
    stop("Package 'phytools' is needed for function 'add_replicates'.", call. = FALSE)
  }
  if (!is.null(traits[[id]]) && (length(traits[[id]]) != length(unique(traits[[id]])))){
    stop("The samples ids in column named ", id, " should be unique identifiers of the samples.")
  }
  if (is.null(traits[[species]])) {
    stop("The `traits` data frame should contain a column named ", species, " with species names for each sample. Plead adjust argument `species` accordingly.")
  }
  if (is.null(traits[[id]]) && id != "auto") {
    stop("The `traits` data frame should contain a column named ", id, " with sample ids. Plead adjust argument `id` accordingly.")
  }
  data_tree_cor <- match(traits[[id]], tree$tip.label)
  if (anyNA(data_tree_cor)) {
    # Species in data NOT in the tree
    stop("Species '", paste(unique(traits[[id]][is.na(data_tree_cor)]), collapse = "', '"), "' are in the data but not in the tree. Please remove them from the data before proceding (or use a tree that include them)." )
  }
  ## Make tree
  tree_norep <- tree
  # Remove replicates
  tree_norep <- ape::drop.tip(tree_norep, tree$tip.label[duplicated(traits[[species]])])
  # result
  return(tree_norep)
}

#' @importFrom utils capture.output
NULL

#' @title Get the edges in theta2 with parsimony
#'
#' @description
#' Utility function to paint the edges of the tree from conditions at the tips.
#'
#' @param cdata a \code{phyloCompData} object
#' @param tree_norep a phylogenetic tree without replicates
#' 
#' @return a vector specifying the regime of each branch of the tree, to be fed
#' to the \code{isTheta2edge} of function \code{evemodel::fitTwoTheta}
#'
#' @keywords internal
#'
getIsTheta2edge <- function(cdata, tree_norep = getTreeEVE(cdata)) {
  if (!requireNamespace("mvSLOUCH", quietly = TRUE)) {
    stop("Package 'mvSLOUCH' is needed for function 'getIsTheta2edge'.", call. = FALSE)
  }
  sample_annotations <- sample.annotations(cdata)
  if (!all(sapply(split(sample_annotations, sample_annotations$id.species), function(x) length(unique(x$condition)) == 1))) {
    stop("For evemodel, all the samples from a species must be in the same condition.", call. = FALSE)
  }
  ## Parsimony reconstruction
  unique_species_ind <- !duplicated(sample_annotations$id.species)
  cond_no_rep <- sample_annotations$condition[unique_species_ind]
  names(cond_no_rep) <- rownames(sample_annotations)[unique_species_ind]
  cond_no_rep <- cond_no_rep[match(names(cond_no_rep), tree_norep$tip.label)]
  oo <- capture.output(pars <- suppressMessages(mvSLOUCH::fitch.mvsl(tree_norep, cond_no_rep, acctran ="TRUE")))
  return(pars$branch_regimes == "2")
}

#' @title Log fold change
#' 
#' @description
#' Get the log fold change from a two theta fit with evemodel.
#' 
#'
#' @param twoThetaRes result of function \code{evemodel::fitTwoTheta}
#' @param isTheta2edge see matching parameter of \code{evemodel::fitTwoTheta}
#' @param tree_norep tree used in \code{evemodel::fitTwoTheta}
#'
#' @return vector of logFC
#'
#' @keywords internal
#'
getlogFCEVE <- function(twoThetaRes, isTheta2edge, tree_norep) {
  
  getDiffMeanTips <- function(par_eve) {
    thetas <- isTheta2edge 
    thetas[isTheta2edge] <- par_eve["theta2"]
    thetas[!isTheta2edge] <- par_eve["theta1"]
    alphas <- rep(par_eve["alpha"], nrow(tree_norep$edge))
    ee <- evemodel::expectedMeanOU(tree_norep, thetas, alphas, par_eve["theta1"])
    dd <- unique(ee[1:length(tree_norep)]) - par_eve["theta1"]
    return(dd[which.max(abs(dd))])
  }
  
  return(apply(twoThetaRes$par, 1, getDiffMeanTips))
}

#' Generate a \code{.Rmd} file containing code to perform differential expression analysis with \code{\link[evemodel]{twoThetaTest}}.
#' 
#' A function to generate code that can be run to perform differential expression analysis of RNAseq data (comparing two conditions) using the evemodel package. The code is written to a \code{.Rmd} file. This function is generally not called by the user, the main interface for performing differential expression analysis is the \code{\link{runDiffExp}} function.
#' 
#' For more information about the methods and the interpretation of the parameters, see the \code{\link[evemodel]{twoThetaTest}} package and the corresponding publications. 
#' 
#' @param data.path The path to a .rds file containing the \code{phyloCompData} object that will be used for the differential expression analysis.
#' @param result.path The path to the file where the result object will be saved.
#' @param codefile The path to the file where the code will be written.
#' @param norm.method The between-sample normalization method used to compensate for varying library sizes and composition in the differential expression analysis. The normalization factors are calculated using the \code{calcNormFactors} of the \code{edgeR} package. Possible values are \code{"TMM"}, \code{"RLE"}, \code{"upperquartile"} and \code{"none"}
#' @param length.normalization one of "none" (no correction), "TPM" or "RPKM" (default). See details.
#' @param data.transformation one of "log2", "asin(sqrt)" or "sqrt". Data transformation to apply to the normalized data.
#' @param empirical.p.values Boolean (default to FALSE). If TRUE, then an empirical null distribution is generated, using the parameters estimated from the oneThetaFits, and the \code{simOneTheta} function. See \code{evemodel} package vignette for more details.
#' @param n.genes.null.dist if \code{empirical.p.values=TRUE}, the number of genes to simulate under the null distribution. Default to 1000.
#' @param ... Further arguments to be passed to function \code{\link[evemodel]{twoThetaTest}}.
#' 
#' @details 
#' The \code{length.matrix} field of the \code{phyloCompData} object 
#' is used to normalize the counts, using one of the following formulas:
#' * \code{length.normalization="none"} : \eqn{CPM_{gi} = \frac{N_{gi} + 0.5}{NF_i \times \sum_{g} N_{gi} + 1} \times 10^6}
#' * \code{length.normalization="TPM"} : \eqn{TPM_{gi} = \frac{(N_{gi} + 0.5) / L_{gi}}{NF_i \times \sum_{g} N_{gi}/L_{gi} + 1} \times 10^6}
#' * \code{length.normalization="RPKM"} : \eqn{RPKM_{gi} = \frac{(N_{gi} + 0.5) / L_{gi}}{NF_i \times \sum_{g} N_{gi} + 1} \times 10^9}
#' 
#' where \eqn{N_{gi}} is the count for gene g and sample i,
#' where \eqn{L_{gi}} is the length of gene g in sample i,
#' and \eqn{NF_i} is the normalization for sample i,
#' normalized using \code{calcNormFactors} of the \code{edgeR} package.
#' 
#' The function specified by the \code{data.transformation} is then applied
#' to the normalized count matrix.
#' 
#' The "\eqn{+0.5}" and "\eqn{+1}" are taken from Law et al 2014,
#' and dropped from the normalization 
#' when the transformation is something else than \code{log2}.
#' 
#' The "\eqn{\times 10^6}" and "\eqn{\times 10^9}" factors are omitted when
#' the \code{asin(sqrt)} transformation is taken, as \eqn{asin} can only
#' be applied to real numbers smaller than 1.
#' 
#' @export 
#' @author Charlotte Soneson, Paul Bastide, Mélina Gallopin
#' @return The function generates a \code{.Rmd} file containing the code for performing the differential expression analysis. This file can be executed using e.g. the \code{knitr} package.
#' @references 
#' ROHLFS, Rori V. & NIELSEN, Rasmus. Phylogenetic ANOVA: The Expression Variance and Evolution Model for Quantitative Trait Evolution. 2015, 695–708. (64). ISSN: 1063-5157.
#' 
#' Law, C.W., Chen, Y., Shi, W. et al. (2014) voom: precision weights unlock linear model analysis tools for RNA-seq read counts. Genome Biol 15, R29.
#'
#' Musser, JM, Wagner, GP. (2015): Character trees from transcriptome data: Origin and individuation of morphological characters and the so‐called “species signal”. J. Exp. Zool. (Mol. Dev. Evol.) 324B: 588– 604.
#' 
#' @examples
#' try(
#' if (require(ape) && require(evemodel)) {
#' tmpdir <- normalizePath(tempdir(), winslash = "/")
#' set.seed(20200317)
#' tree <- ape::read.tree(
#'   text = "(((A1:0,A2:0,A3:0):1,B1:1):1,((C1:0,C2:0):1.5,(D1:0,D2:0):1.5):0.5);")
#' mydata.obj <- generateSyntheticData(dataset = "mydata", n.vars = 100, 
#'                                     samples.per.cond = 4, n.diffexp = 10, 
#'                                     tree = tree,
#'                                     id.species = factor(c("A", "A", "A", "B", "C", "C", "D", "D")),
#'                                     lengths.relmeans = rpois(100, 1000),
#'                                     lengths.dispersions = rgamma(100, 1, 1),
#'                                     output.file = file.path(tmpdir, "mydata.rds"))
#' ## Diff Exp
#' runDiffExp(data.file = file.path(tmpdir, "mydata.rds"), result.extent = "twoThetaTest", 
#'            Rmdfunction = "evemodel.twoThetaTest.createRmd", 
#'            output.directory = tmpdir,
#'            norm.method = "TMM",
#'            length.normalization = "RPKM")
#' generateCodeHTMLs(file.path(tmpdir, "mydata_twoThetaTest.rds"), tmpdir)
#' \dontrun{
#' file.show(file.path(tmpdir, "mydata_twoThetaTest_Code.html"))
#' }
#' })
evemodel.twoThetaTest.createRmd <- function(data.path, result.path, codefile, 
                                            norm.method,
                                            length.normalization = "RPKM",
                                            data.transformation = "log2",
                                            empirical.p.values = FALSE,
                                            n.genes.null.dist = 1000,
                                            ...) {
  codefile <- file(codefile, open = 'w')
  writeLines("### evemodel twoThetaTest", codefile)
  writeLines(paste("Data file: ", data.path, sep = ''), codefile)
  writeLines(c("```{r, echo = TRUE, eval = TRUE, include = TRUE, message = TRUE, error = TRUE, warning = TRUE}", 
               "require(evemodel)", 
               "require(limma)", 
               "require(edgeR)",
               paste("cdata <- readRDS('", data.path, "')", sep = '')), codefile)
  if (is.list(readRDS(data.path))) {
    writeLines("cdata <- convertListTophyloCompData(cdata)", codefile)
  }
  
  writeLines(c("is.valid <- check_phyloCompData(cdata)",
               "if (!(is.valid == TRUE)) stop('Not a valid phyloCompData object.')"),
             codefile)
  ## Tree and parameters
  writeLines(c("tree_rep <- getTree(cdata)"), codefile)
  writeLines(c("tree_norep <- getTreeEVE(cdata)"), codefile)
  writeLines(c("theta_2_vec <- getIsTheta2edge(cdata, tree_norep)"), codefile)
  writeLines(c("col_species <- tree_norep$tip.label[sample.annotations(cdata)$id.species]"), codefile)
  ## Normalization
  writeNormalization(norm.method, length.normalization, data.transformation, codefile)
  ## Apply analysis
  writeLines(c("", "# Analysis with EVE"),codefile)
  extra_args <- eval(substitute(alist(...)))
  extra_args <- sapply(extra_args, function(x) paste(" = ", x))
  extra_args <- paste(names(extra_args), extra_args, collapse = ", ")
  if (length(extra_args) > 1) extra_args <- paste0(", ", extra_args)
  writeLines(
    paste0("evemodel.results_list <- evemodel::twoThetaTest(tree = tree_norep, gene.data = data.trans, isTheta2edge = theta_2_vec, colSpecies = col_species", extra_args, ")"),
    codefile)
  if (empirical.p.values) {
    writeLines(c(
      "parEstim <- apply(evemodel.results_list$oneThetaRes$par, 2, median)",
      "set.seed(1289)",
      paste0("nullData <- simOneTheta(n = ", n.genes.null.dist, ", tree = tree_norep, colSpecies = col_species, theta = parEstim[\"theta\"], sigma2 = parEstim[\"sigma2\"], alpha = parEstim[\"alpha\"], beta = parEstim[\"beta\"])"),
      "test.nullData_full <- twoThetaTest(tree = tree_norep, gene.data = nullData, isTheta2edge = theta_2_vec, colSpecies = col_species)",
      "emp_cff <- ecdf(test.nullData_full$LRT)",
      "result.table <- data.frame(pvalue = 1 - emp_cff(evemodel.results_list$LRT), logFC = getlogFCEVE(evemodel.results_list$twoThetaRes, theta_2_vec, tree_norep))"),
      codefile)
  } else {
    writeLines(c(
      "result.table <- data.frame(pvalue = pchisq(evemodel.results_list$LRT, df = 1, lower.tail = FALSE), logFC = getlogFCEVE(evemodel.results_list$twoThetaRes, theta_2_vec, tree_norep))"),
      codefile)
  }
  writeLines(c(
    "result.table$score <- 1 - result.table$pvalue",
    "result.table$adjpvalue <- p.adjust(result.table$pvalue, 'BH')"),
    codefile)
  writeLines(c("", "# Save the results"),codefile)
  writeLines(c(
    "rownames(result.table) <- rownames(count.matrix(cdata))",
    "result.table(cdata) <- result.table", 
    "package.version(cdata) <- paste('evemodel,', packageVersion('evemodel'))",
    "package.version(cdata) <- paste('edgeR,', packageVersion('edgeR'))",
    "analysis.date(cdata) <- date()",
    paste("method.names(cdata) <- list('short.name' = 'evetwotheta', 'full.name' = '",
          paste('evemodel', packageVersion('evemodel'), '.', norm.method, '.',
                "lengthNorm.", length.normalization, '.',
                "dataTrans.", data.transformation, '.',
                "empNull.", empirical.p.values,
                ifelse(empirical.p.values, paste0(".nGenesNull.", n.genes.null.dist), ""),
                sep = ''),
          "')", sep = ''),
    "is.valid <- check_compData_results(cdata)",
    "if (!(is.valid == TRUE)) stop('Not a valid phyloCompData result object.')",
    paste("saveRDS(cdata, '", result.path, "')", sep = "")),
    codefile)  
  writeLines("print(paste('Unique data set ID:', info.parameters(cdata)$uID))", codefile)
  writeLines("sessionInfo()", codefile)
  writeLines("```", codefile)
  close(codefile)
}
