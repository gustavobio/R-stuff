# Get a partial hclust object
subtree <- function(hc, nl) {
  if (class(hc) != "hclust") stop("Not a hclust object")
  if (!all(nl %in% hc$labels)) stop("Some labels do not match those
    in the hclust object")
  nl <- hc$labels %in% nl
  m <- as.matrix(cophenetic(hc))
  d <- as.dist(m[nl, nl])
  return(update(hc, d = d))
}

# Descending leaves for each branch
desc.leaves <- function(hc) {
  signs <- sign(hc$merge)
  branches <- 1:NROW(hc$merge)
  descending.leaves <- vector("list", length(branches))
  heights <- numeric(length(branches))
  counter <- 1
  for (i in branches) {
    j <- i
    r.leaves <- c()
    while(length(j) > 0) {
      r <- hc$merge[j, ]
      j <- r[r > 0]
      r.leaves <- c(r.leaves, abs(r[r < 0]))
    }
    asc.node <- which(hc$merge == i, arr.ind = TRUE)[1]
    if (!is.na(asc.node)) {
      heights[i] <- hc$height[asc.node] - hc$height[i]
    } else {
      heights[i] <- hc$height[i]
    }
    descending.leaves[[counter]] <- r.leaves
    counter <- counter + 1
  }
  names(descending.leaves) <- branches
  return(list(dl = descending.leaves, length = heights))
}

# Functional entropy using only tips
fe.partial.leaves <- function(hc, abundances) {
  leaves <- which(hc$merge < 0)
  leaves.name <- abs(hc$merge[leaves])
  abundances <- abundances[hc$labels]
  heights <- c(hc$height, hc$height)
  rows <- ifelse(leaves <= NROW(hc$merge), leaves, leaves -
    NROW(hc$merge))
  fe.leaves <- heights[leaves] * abundances *
  log(abundances)
  return(fe.leaves)
}

# Functional entropy using only branches
fe.partial.branches <- function(hc, abundances) {
  abundances <- abundances[hc$labels]
  abundances <- abundances/sum(abundances)
  descending.leaves <- desc.leaves(hc)
  prop.ab <- sapply(descending.leaves$dl, function(x)
    {sum(abundances[x])})
  log.prop.ab <- log(prop.ab)
  return(descending.leaves$length * prop.ab * log(prop.ab))
}

# Functional entropy for a single site
# Divided between branches and tips to make testing easier
fe <- function(hc, abundances) {
  abundances <- abundances[abundances > 0]
  if (sum(abundances != 1)) abundances <- abundances/sum(abundances)
  if (sum(abundances) == 0) stop("None of the abundances is higher than
    0")
  if (class(hc) != "hclust") stop("The tree you passed does not look
    like a hclust object")
  if (is.null(names(abundances)) | !all(names(abundances) %in% hc$labels)) stop("Tree labels do not
    match site labels")
  return(abs(sum(fe.partial.leaves(hc, abundances),
        fe.partial.branches(hc, abundances))))
}

# Functional entropy for community data
fe.community <- function(hc, comm) {
  fe.sites <- numeric(NROW(comm))
  names(fe.sites) <- rownames(comm)
  for (i in 1:NROW(comm)) {
    species <- unlist(comm[i, comm[i,] >0])
    cluster <- subtree(hc, names(species))
    fe.sites[i] <- fe(cluster, species)
  }
  fe.sites
}

# Example:
sp <- 20
traits <- sapply(1:7, function(x) rnorm(sp))
colnames(traits) <- LETTERS[1:7]
rownames(traits) <- letters[1:sp]
hc <- hclust(dist(traits))
comm <- do.call(cbind, lapply(1:sp, function(x) sample(0:10, 50, re = T)))
colnames(comm) <- rownames(traits)

fe(hc, colSums(comm))
fe.community(hc, comm)
