context("Lengths and length factors")

test_that("computeFactorLengths is equivalent to trueProba", {
  set.seed(20200430)
  
  n.sample <- 3
  n.vars <- 5
  
  ## Mock data
  count.matrix <- matrix(c(5,3,8,17,23,42,10,13,27,752,615,1203,1507,1225,2455), ncol = n.sample, byrow = T)
  truemeans <- rowMeans(count.matrix)
  leng <- matrix(c(rep(1000, 2 * n.vars - 1), 2000, 1:n.vars*1000), ncol = n.sample, byrow = F)
  
  ## True proba
  trueMeansLeng <- diag(truemeans) %*% leng
  trueProbas <- trueMeansLeng %*% diag(1 / colSums(trueMeansLeng))
  
  ## nfact
  lfact <- computeFactorLengths(leng, truemeans, sum(truemeans))
  
  ## Equal
  expect_equal(trueProbas, diag(truemeans / sum(truemeans)) %*% lfact)
  
  ## Sequencing depths
  seqdepth <- 2500
  nfacts <- runif(n.sample, min = 0.7, max = 1.4)
  seq.depths <- nfacts * seqdepth
  
  ## expectations
  exp_1 <- trueProbas %*% diag(seq.depths)
  expect_equal(colSums(exp_1), seq.depths)
  
  exp_2 <- diag(truemeans / sum(truemeans)) %*% lfact %*% diag(seq.depths)
  expect_equal(exp_1, exp_2)
})

test_that("computeFactorLengths is equivalent to trueProba", {
  skip_if_not_installed("phytools")
  set.seed(20200430)
  
  ntaxa <- 10
  nrep <- 3
  n.sample <- ntaxa * nrep
  n.vars <- 50
  
  ## Mock tree
  set.seed("17900714")
  tree <- ape::rphylo(ntaxa, birth = 0.1, death = 0)
  # rescale to unit height
  tree$edge.length <- tree$edge.length / max(ape::node.depth.edgelength(tree))
  # rename tims
  tree$tip.label <- paste0("t", 1:ntaxa)
  tree <- add_replicates(tree, nrep)
  
  id.species <- rep(1:ntaxa, each = nrep)
  names(id.species) <- tree$tip.label
  id.species <- as.factor(id.species)
  
  id.species <- checkSpecies(id.species, "id.species", tree, 1e-10, TRUE)
  
  ## Mock data
  count.matrix <- matrix(rpois(n.sample * n.vars, 1:10 * 10), ncol = n.sample, byrow = T)
  truemeans <- rowMeans(count.matrix)
  leng <- matrix(c(rep(1000, 2 * n.vars - 1), 2000, 1:n.vars*1000), ncol = n.sample, byrow = F)
  lengths.relmeans <- 1:n.vars * 1000
  lengths.dispersions <- 1:n.vars / 10
  expect_warning(leng <- generateLengths(id.species, lengths.relmeans, lengths.dispersions),
                 "to get a non zero value. Please check")
  
  expect_equal(leng[, 1], leng[, 2])
  expect_equal(leng[, 4], leng[, 5])
  
  ## True proba
  trueMeansLeng <- diag(truemeans) %*% leng
  trueProbas <- trueMeansLeng %*% diag(1 / colSums(trueMeansLeng))
  
  ## nfact
  lfact <- computeFactorLengths(leng, truemeans, sum(truemeans))
  
  ## Equal
  expect_equal(trueProbas, diag(truemeans / sum(truemeans)) %*% lfact)
  
  ## Sequencing depths
  seqdepth <- 2500
  nfacts <- runif(n.sample, min = 0.7, max = 1.4)
  seq.depths <- nfacts * seqdepth
  
  ## expectations
  exp_1 <- trueProbas %*% diag(seq.depths)
  expect_equal(colSums(exp_1), seq.depths)
  
  exp_2 <- diag(truemeans / sum(truemeans)) %*% lfact %*% diag(seq.depths)
  expect_equal(exp_1, exp_2)
})

test_that("Zero lengths errors", {
  set.seed(20200430)
  
  tdir <- tempdir()
  
  n.vars <- 500
  tree_rep <- ape::read.tree(text = "(((A1:0,A2:0,A3:0):1,B1:1):1,((C1:0,C2:0):1.5,(D1:0,D2:0):1.5):0.5);")
  idcondition = c(1, 1, 1, 1, 2, 2, 2, 2)
  idspecies = factor(c("A", "A", "A", "B", "C", "C", "D", "D"))
  names(idcondition) <- names(idspecies) <- tree_rep$tip.label
  species_conditions <- c(1, 1, 2, 2)
  
  ## Zero length means
  lengths.relmeans <- rpois(n.vars, 1000)
  lengths.relmeans[c(2, 34, 76)] <- 0
  
  expect_error(
    generateSyntheticData(
      dataset = "test", n.vars = n.vars, 
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
      lengths.relmeans = lengths.relmeans,
      lengths.dispersions = "auto",
      output.file = file.path(tdir, "test.rds"),
      effect.size = 10
    ),
    "zero length expectation."
  )
  
  ## Zero length generated
  set.seed(1289)
  lengths.relmeans <- rpois(n.vars, 1000)
  lengths.dispersions <- runif(n.vars, 0, 1)
  lengths.relmeans[225] <- 1
  lengths.dispersions[225] <- 5
  expect_error(generateLengths(idspecies, lengths.relmeans, lengths.dispersions),
               "could not generate non-zero lengths for gene")
  lengths.relmeans[225] <- 4
  lengths.dispersions[225] <- 5
  expect_warning(generateLengths(idspecies, lengths.relmeans, lengths.dispersions),
                 "to get a non zero value. Please check that the associated")
  
  ## Zero length phylo
  set.seed(1289)
  lengths.relmeans <- rpois(n.vars, 1000)
  lengths.dispersions <- runif(n.vars, 0, 1)
  lengths.relmeans[225] <- 1
  lengths.dispersions[225] <- 100
  expect_error(generateLengthsPhylo(tree_rep, idspecies, lengths.relmeans, lengths.dispersions),
               "could not generate non-zero lengths for gene")
  lengths.relmeans[225] <- 1
  lengths.dispersions[225] <- 5
  expect_warning(generateLengthsPhylo(tree_rep, idspecies, lengths.relmeans, lengths.dispersions),
                 "to get a non zero value. Please check that the associated")
  
})