test_that("EVE vs phylolm comparisons", {
  skip_if_not_installed("phytools")
  skip_if_not_installed("mvSLOUCH")
  set.seed(1289)
  tdir <- tempdir()
  
  ## Phylogentic trees with without replicates
  tree_rep <- ape::read.tree(text = "(((A1:0,A2:0,A3:0):1,B1:1):1,((C1:0,C2:0):1.5,(D1:0,D2:0):1.5):0.5);")
  tree <- ape::read.tree(text = "((A1:1,B1:1):1,(C1:1.5,D1:1.5):0.5);")
  
  ## Sample annotations
  idcondition = c(1, 1, 1, 1, 2, 2, 2, 2)
  idspecies = factor(c("A", "A", "A", "B", "C", "C", "D", "D"))
  names(idcondition) <- names(idspecies) <- tree_rep$tip.label
  species_conditions <- c(1, 1, 2, 2)
  
  ## Simulated data
  testdat <- generateSyntheticData(
    dataset = "test", n.vars = 500, 
    n.diffexp = 50, 
    repl.id = 1, seqdepth = 1e5, 
    fraction.upregulated = 0.5, 
    between.group.diffdisp = FALSE, 
    filter.threshold.total = 1, 
    filter.threshold.mediancpm = 0, 
    fraction.non.overdispersed = 0, 
    tree = tree_rep,
    model.process = "OU",
    selection.strength = 2,
    id.condition = idcondition,
    id.species =  idspecies,
    lengths.relmeans = "auto",
    lengths.dispersions = "auto",
    output.file = file.path(tdir, "test.rds"),
    effect.size = 10
  )
  
  ## Check trees
  tree_rep_dat <- getTree(testdat)
  expect_equal(tree_rep_dat, tree_rep)
  
  tree_dat <- getTreeEVE(testdat)
  expect_equal(tree_dat, tree)
  
  ## Check thetas
  theta2_edges <- getIsTheta2edge(testdat, tree_dat)
  theta2_edges_num <- theta2_edges + 1
  cond_tips <- species_conditions[tree_dat$edge[, 2]]
  expect_true(all(theta2_edges_num[!is.na(cond_tips)] == cond_tips[!is.na(cond_tips)]))
  
  ## Data trans
  nf <- edgeR::calcNormFactors(count.matrix(testdat) / length.matrix(testdat), method = "TMM")
  data.trans <- phylolimma::lengthNormalizeRNASeq(count.matrix(testdat),
                                                  length.matrix(testdat),
                                                  normalisationFactor = nf,
                                                  lengthNormalization = "TPM",
                                                  dataTransformation = "log2")
  
  ## phylolm
  min_sigma2_error <- (.Machine$double.eps)^0.5 * max(ape::node.depth.edgelength(tree_rep_dat))
  cond <- as.factor(sample.annotations(testdat)$condition)
  
  all_res_lm <- apply(data.trans[1:100, ], 1,
                      function(dd) suppressWarnings(phylolm::phylolm(resp ~ condition, phy = tree_rep_dat,
                                                                     data = data.frame(resp = dd, condition = cond),
                                                                     model = "OUfixedRoot", measurement_error = TRUE,
                                                                     lower.bound = list(sigma2_error = min_sigma2_error))))
  
  res_lm <- as.data.frame(t(sapply(all_res_lm, extract_results_phylolm)))

  ## eve
  all_res_eve <- evemodel::twoThetaTest(tree = tree_dat,
                                        gene.data = data.trans[1:100, ],
                                        isTheta2edge = theta2_edges,
                                        colSpecies = tree_dat$tip.label[sample.annotations(testdat)$id.species])
  
  res_eve <- data.frame(pvalue = pchisq(all_res_eve$LRT, df = 1, lower.tail = F))
  
  expect_equal(2 * (all_res_eve$twoThetaRes$ll - all_res_eve$oneThetaRes$ll),
               all_res_eve$LRT)
  
  ## likelihoods are equal
  expect_equal(all_res_eve$twoThetaRes$ll,
               unname(sapply(all_res_lm, function(x) x$logLik)),
               tolerance = 1e-3) 
  
  ## sigma2_error are equal
  get_sigma2_error <- function(par_eve) {
    par_eve["sigma2"] / 2 / par_eve["alpha"] * par_eve["beta"]
  }
  expect_equal(apply(all_res_eve$twoThetaRes$par, 1, get_sigma2_error),
               unname(sapply(all_res_lm, function(x) x$sigma2_error)),
               tolerance = 1e-2)
  
  ## expectations are equal
  get_exp_tips <- function(par_eve) {
    thetas <- theta2_edges 
    thetas[theta2_edges] <- par_eve["theta2"]
    thetas[!theta2_edges] <- par_eve["theta1"]
    alphas <- rep(par_eve["alpha"], nrow(tree_dat$edge))
    ee <- evemodel::expectedMeanOU(tree_dat, thetas, alphas, par_eve["theta1"])
    return(unique(ee[1:length(tree_dat$tip.label)]))
  }
  means_eve <- apply(all_res_eve$twoThetaRes$par, 1, get_exp_tips)
  means_lm <- unname(sapply(all_res_lm, function(x) unique(predict(x))))

  expect_equal(means_eve, means_lm, tolerance = 1e-2)
  
  expect_equal(all_res_eve$twoThetaRes$par[, "theta1"],
               unname(sapply(all_res_lm, function(x) x$coefficients[1])),
               tol = 1e-2)
  
  ## logFC
  res_eve$logFC <- getlogFCEVE(all_res_eve$twoThetaRes, theta2_edges, tree_dat)
  
  expect_equal(res_eve$logFC, unname(unlist(res_lm$logFC)), tolerance = 1e-3)
  
  ## p values
  res_eve$adjpvalue <- p.adjust(res_eve$pvalue, 'BH')
  FD_eve <- sum(res_eve$adjpvalue[-(1:50)] <= 0.05)
  TD_eve <- sum(res_eve$adjpvalue[1:50] <= 0.05)
  res_lm$adjpvalue <- p.adjust(res_lm$pvalue, 'BH')
  FD_lm <- sum(res_lm$adjpvalue[-(1:50)] <= 0.05)
  TD_lm <- sum(res_lm$adjpvalue[1:50] <= 0.05)
  
  expect_equal(FD_lm, 0)
  expect_equal(FD_eve, 2)
  
})

test_that("EVE vs phylolm comparisons - bigger tree", {
  skip_if_not_installed("phytools")
  skip_if_not_installed("mvSLOUCH")
  set.seed(1289)
  tdir <- tempdir()
  
  ## Phylogentic trees with without replicates
  tree <- ape::rphylo(20, 0.1, 0)
  tree_rep <- add_replicates(tree, 3)
  tree$tip.label <- paste0(tree$tip.label, "_3")
  
  ## Sample annotations
  idcondition = rep(1, length(tree_rep$tip.label))
  idcondition[grepl("t5|t3|t2_|t17", tree_rep$tip.label)] <- 2
  idspecies = as.factor(rep(tree$tip.label, each = 3))
  names(idcondition) <- names(idspecies) <- tree_rep$tip.label
  species_conditions <- idcondition[1 + (0:19)*3]
  
  ## Simulated data
  testdat <- generateSyntheticData(
    dataset = "test", n.vars = 500, 
    n.diffexp = 50, 
    repl.id = 1, seqdepth = 1e5, 
    fraction.upregulated = 0.5, 
    between.group.diffdisp = FALSE, 
    filter.threshold.total = 1, 
    filter.threshold.mediancpm = 0, 
    fraction.non.overdispersed = 0, 
    tree = tree_rep,
    model.process = "OU",
    selection.strength = 2,
    id.condition = idcondition,
    id.species =  idspecies,
    lengths.relmeans = "auto",
    lengths.dispersions = "auto",
    output.file = file.path(tdir, "test.rds"),
    effect.size = 1
  )
  
  ## Check trees
  tree_rep_dat <- getTree(testdat)
  expect_equal(tree_rep_dat, tree_rep)
  
  tree_dat <- getTreeEVE(testdat)
  expect_equal(tree_dat, tree)
  
  ## Check thetas
  theta2_edges <- getIsTheta2edge(testdat, tree_dat)
  theta2_edges_num <- theta2_edges + 1
  cond_tips <- species_conditions[tree_dat$edge[, 2]]
  expect_true(all(theta2_edges_num[!is.na(cond_tips)] == cond_tips[!is.na(cond_tips)]))
  # plot(tree_dat)
  # ape::edgelabels(theta2_edges_num)
  
  ## Data trans
  nf <- edgeR::calcNormFactors(count.matrix(testdat) / length.matrix(testdat), method = "TMM")
  data.trans <- phylolimma::lengthNormalizeRNASeq(count.matrix(testdat),
                                                  length.matrix(testdat),
                                                  normalisationFactor = nf,
                                                  lengthNormalization = "TPM",
                                                  dataTransformation = "log2")
  
  ## phylolm
  min_sigma2_error <- (.Machine$double.eps)^0.5 * max(ape::node.depth.edgelength(tree_rep_dat))
  cond <- as.factor(sample.annotations(testdat)$condition)
  
  all_res_lm <- apply(data.trans[1:100, ], 1,
                      function(dd) suppressWarnings(phylolm::phylolm(resp ~ condition, phy = tree_rep_dat,
                                                                     data = data.frame(resp = dd, condition = cond),
                                                                     model = "OUrandomRoot", measurement_error = TRUE,
                                                                     lower.bound = list(sigma2_error = min_sigma2_error))))
  
  res_lm <- as.data.frame(t(sapply(all_res_lm, extract_results_phylolm)))
  
  ## eve
  all_res_eve <- evemodel::twoThetaTest(tree = tree_dat,
                                        gene.data = data.trans[1:100, ],
                                        isTheta2edge = theta2_edges,
                                        colSpecies = tree_dat$tip.label[sample.annotations(testdat)$id.species])
  
  res_eve <- data.frame(pvalue = pchisq(all_res_eve$LRT, df = 1, lower.tail = F))
  
  expect_equal(2 * (all_res_eve$twoThetaRes$ll - all_res_eve$oneThetaRes$ll),
               all_res_eve$LRT)
  
  ## likelihoods are similar
  expect_equal(all_res_eve$twoThetaRes$ll,
               unname(sapply(all_res_lm, function(x) x$logLik)),
               tolerance = 1e-1) 
  
  ## sigma2_error are equal
  get_sigma2_error <- function(par_eve) {
    par_eve["sigma2"] / 2 / par_eve["alpha"] * par_eve["beta"]
  }
  expect_equal(apply(all_res_eve$twoThetaRes$par, 1, get_sigma2_error),
               unname(sapply(all_res_lm, function(x) x$sigma2_error)),
               tolerance = 1e-1)
  
  ## expectations are equal
  get_exp_tips <- function(par_eve) {
    thetas <- theta2_edges 
    thetas[theta2_edges] <- par_eve["theta2"]
    thetas[!theta2_edges] <- par_eve["theta1"]
    alphas <- rep(par_eve["alpha"], nrow(tree_dat$edge))
    ee <- evemodel::expectedMeanOU(tree_dat, thetas, alphas, par_eve["theta1"])
    return(unique(round(ee[1:length(tree_dat$tip.label)], 5)))
  }
  means_eve <- apply(all_res_eve$twoThetaRes$par, 1, get_exp_tips)
  means_lm <- unname(sapply(all_res_lm, function(x) unique(predict(x))))
  
  expect_equal(means_eve, means_lm, tolerance = 1e-1)
  
  expect_equal(all_res_eve$twoThetaRes$par[, "theta1"],
               unname(sapply(all_res_lm, function(x) x$coefficients[1])),
               tol = 1e-1)
  
  ## logFC are different
  res_eve$logFC <- getlogFCEVE(all_res_eve$twoThetaRes, theta2_edges, tree_dat)
  
  expect_equal(res_eve$logFC, unname(unlist(res_lm$logFC)), tolerance = 1e0)
  
  ## p values
  res_eve$adjpvalue <- p.adjust(res_eve$pvalue, 'BH')
  FD_eve <- sum(res_eve$adjpvalue[-(1:50)] <= 0.05)
  TD_eve <- sum(res_eve$adjpvalue[1:50] <= 0.05)
  res_lm$adjpvalue <- p.adjust(res_lm$pvalue, 'BH')
  FD_lm <- sum(res_lm$adjpvalue[-(1:50)] <= 0.05)
  TD_lm <- sum(res_lm$adjpvalue[1:50] <= 0.05)
  
  expect_equal(FD_lm, 1)
  expect_equal(TD_lm, 16)
  expect_equal(FD_eve, 0)
  expect_equal(TD_eve, 0)
  
})

test_that("Errors", {
  skip_if_not_installed("phytools")
  skip_if_not_installed("phangorn")
  set.seed(1289)
  tdir <- tempdir()
  
  ## Phylogentic trees with without replicates
  tree_rep <- ape::read.tree(text = "(((A1:0,A2:0,A3:0):1,B1:1):1,((C1:0,C2:0):1.5,(D1:0,D2:0):1.5):0.5);")
  ## Sample annotations
  idcondition = c(1, 1, 2, 1, 2, 2, 2, 2)
  idspecies = factor(c("A", "A", "A", "B", "C", "C", "D", "D"))
  names(idcondition) <- names(idspecies) <- tree_rep$tip.label
  ## Simulated data
  testdat <- generateSyntheticData(
    dataset = "test", n.vars = 500, 
    n.diffexp = 50, 
    repl.id = 1, seqdepth = 1e5, 
    fraction.upregulated = 0.5, 
    between.group.diffdisp = FALSE, 
    filter.threshold.total = 1, 
    filter.threshold.mediancpm = 0, 
    fraction.non.overdispersed = 0, 
    tree = tree_rep,
    id.condition = idcondition,
    id.species =  idspecies,
    lengths.relmeans = "auto",
    lengths.dispersions = "auto",
    output.file = file.path(tdir, "test.rds")
  )
  
  expect_error(getIsTheta2edge(testdat),
               "For evemodel, all the samples from a species must be in the same condition.")
  
})

test_that("evemodel runComparison", {
  skip_if_not_installed("phytools")
  skip_if_not_installed("mvSLOUCH")
  set.seed(1289)
  tdir <- tempdir()
  # library(here)
  # tdir <- here()
  
  ## Phylogentic trees with without replicates
  tree <- ape::rphylo(20, 0.1, 0)
  tree_rep <- add_replicates(tree, 3)
  tree$tip.label <- paste0(tree$tip.label, "_3")
  
  ## Sample annotations
  idcondition = rep(1, length(tree_rep$tip.label))
  idcondition[grepl("t5|t3|t2_|t17", tree_rep$tip.label)] <- 2
  idspecies = as.factor(rep(tree$tip.label, each = 3))
  names(idcondition) <- names(idspecies) <- tree_rep$tip.label
  species_conditions <- idcondition[1 + (0:19)*3]
  
  ## Simulated data
  testdat <- generateSyntheticData(
    dataset = "test", n.vars = 100, 
    n.diffexp = 50, 
    repl.id = 1, seqdepth = 1e5, 
    fraction.upregulated = 0.5, 
    between.group.diffdisp = FALSE, 
    filter.threshold.total = 1, 
    filter.threshold.mediancpm = 0, 
    fraction.non.overdispersed = 0, 
    tree = tree_rep,
    model.process = "OU",
    selection.strength = 2,
    id.condition = idcondition,
    id.species =  idspecies,
    lengths.relmeans = "auto",
    lengths.dispersions = "auto",
    output.file = file.path(tdir, "test.rds"),
    effect.size = 10
  )
  
  ## phylolm
  runDiffExp(data.file = file.path(tdir, "test.rds"),
             result.extent = "phylolm", 
             Rmdfunction = "phylolm.createRmd", 
             model = "OUfixedRoot",
             measurement_error = TRUE,
             output.directory = tdir,
             norm.method = "TMM",
             length.normalization = "TPM")
  generateCodeHTMLs(file.path(tdir, "test_phylolm.rds"), tdir)
  res_lm <- readRDS(file.path(tdir, "test_phylolm.rds"))
  expect_true(!anyNA(res_lm@result.table))

  ## eve
  runDiffExp(data.file = file.path(tdir, "test.rds"),
             result.extent = "evemodel", 
             Rmdfunction = "evemodel.twoThetaTest.createRmd", 
             output.directory = tdir,
             norm.method = "TMM",
             length.normalization = "TPM")
  generateCodeHTMLs(file.path(tdir, "test_evemodel.rds"), tdir)
  res_eve <- readRDS(file.path(tdir, "test_evemodel.rds"))
  expect_true(!anyNA(res_eve@result.table))
  pos_test <- res_eve@result.table$adjpvalue <= 0.05
  expect_equal(sum(pos_test), 0)
  # FP <- sum(pos_test[-(1:50)])
  # TP <- sum(pos_test[1:50])
  # TPR <- TP / (FP + TP)
})
