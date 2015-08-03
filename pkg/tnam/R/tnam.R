# This file contains all necessary functions for estimating temporal or 
# cross-sectional network autocorrelation models. Written by Philip Leifeld.

# display version number and date when the package is loaded
.onAttach <- function(libname, pkgname) {
  desc  <- packageDescription(pkgname, libname)
  packageStartupMessage(
    'Package:  tnam\n', 
    'Version:  ', desc$Version, '\n', 
    'Date:     ', desc$Date, '\n', 
    'Authors:  Philip Leifeld (Eawag (ETH) and University of Bern)\n',
    '          Skyler J. Cranmer (The Ohio State University)\n'
  )
}

# check whether 'y' and 'networks' are OK and convert data types if necessary
checkDataTypes <- function(y, networks = NULL, lag = 0) {
  
  # data types of 'y' and 'networks'
  if (is.null(y)) {
    if (is.null(networks)) {
      stop("No 'y' and 'networks' arguments were provided.")
    } else if (class(networks) == "list") {
      time.steps <- length(networks)
      for (i in 1:length(networks)) {
        # go through networks list and check compatibility of each item
        if (class(networks[[i]]) == "matrix") {
          # OK; do nothing
        } else if (class(networks[[i]]) == "network") {
          # convert this network to a matrix and replace inside the list
          networks[[i]] <- as.matrix(networks[[i]])
        } else {
          stop("The 'networks' argument must contain networks or matrices.")
        }
      }
    } else if (class(networks) == "network") {
      # convert to matrix and wrap in list
      time.steps <- 1
      networks <- list(as.matrix(networks))
    } else if (class(networks) == "matrix") {
      # wrap in list
      time.steps <- 1
      networks <- list(networks)
    }
    # create 'y' list with same number of time steps as 'networks' (NA-filled)
    y <- list()
    for (i in 1:time.steps) {
      y[[i]] <- rep(NA, nrow(networks[[i]]))
    }
  } else if (class(y) == "list") {
    for (i in 1:length(y)) {
      if (is.integer(y[[i]])) {
        current.names <- names(y[[i]])
        y[[i]] <- as.numeric(y[[i]])
        names(y[[i]]) <- current.names
      }
    }
    if (is.null(networks)) {
      # insert NA matrices into new list
      networks <- list()
      for (i in 1:length(y)) {
        networks[[i]] <- matrix(NA, nrow = length(y[[i]]), 
            ncol = length(y[[i]]))
        names(networks[[i]]) <- names(y[[i]])
      }
    } else if (class(networks) == "list") {
      # OK; do nothing
    } else if (class(networks) == "network") {
      # convert network to matrix and wrap in list and fill up with NA matrices
      networks <- list(as.matrix(networks))
      if (length(y) > 1) {
        for (i in 2:length(y)) {
          #networks[[i]] <- matrix(NA, nrow = length(y[[i]]), 
          #    ncol = length(y[[i]]))
          networks[[i]] <- as.matrix(networks)
        }
      }
    } else if (class(networks) == "matrix") {
      # wrap matrix in list and fill up with NA matrices
      networks <- list(networks)
      if (length(y) > 1) {
        for (i in 2:length(y)) {
          #networks[[i]] <- matrix(NA, nrow = length(y[[i]]), 
          #    ncol = length(y[[i]]))
          networks[[i]] <- as.matrix(networks)
        }
      }
    } else {
      stop("Data type of the 'networks' argument could not be recognized.")
    }
    # check consistency of lengths
    if (length(y) != length(networks)) {
      stop("Arguments 'y' and 'networks' should be lists of the same length.")
    } else {
      time.steps <- length(y)
    }
  } else if (class(y) == "numeric" || class(y) == "integer") {
    if (class(y) == "integer") {
      nam <- names(y)
      y <- as.numeric(y)
      names(y) <- nam
    }
    y <- list(y)  # wrap y in list
    time.steps <- 1
    if (is.null(networks)) {
      # insert NA objects into new list
      networks <- list()
      for (i in 1:length(y)) {
        networks[[i]] <- matrix(NA, nrow = length(y[[i]]), 
            ncol = length(y[[i]]))
        rownames(networks[[i]]) <- names(y[[i]])
      }
    } else if (class(networks) == "matrix") {
      # wrap matrix in list and fill up with additional NA matrices if necessary
      networks <- list(networks)
      if (length(y) > 1) {
        for (i in 2:length(y)) {
          networks[[i]] <- matrix(NA, nrow = length(y[[i]]), 
              ncol = length(y[[i]]))
          rownames(networks[[i]]) <- names(y[[i]])
        }
      }
    } else if (class(networks) == "network") {
      # convert network to matrix, wrap in list and fill up with NA matrices
      networks <- list(as.matrix(networks))
      if (length(y) > 1) {
        for (i in 2:length(y)) {
          networks[[i]] <- matrix(NA, nrow = length(y[[i]]), 
              ncol = length(y[[i]]))
          rownames(networks[[i]]) <- names(y[[i]])
        }
      }
    } else if (class(networks) == "list") {
      # check length of networks and type of object contained in list
      if (length(networks) != 1) {
        stop("There is only one 'y' vector but there are multiple 'networks'.")
      } else if (class(networks[[1]]) == "matrix") {
        # OK; do nothing
      } else if (class(networks[[1]]) == "network") {
        networks[[1]] <- as.matrix(networks[[1]])
      } else {
        stop("Objects contained in 'networks' list could not be identified.")
      }
    } else {
      stop("Data type of the 'networks' argument could not be recognized.")
    }
  } else if (class(y) == "character") {
    if (is.null(networks)) {
      stop("The response variable must be a numeric vector.")
    } else if (class(networks) == "list") {
      # extract vertex attributes
      y <- lapply(networks, function(x) get.vertex.attribute(x, y))
      time.steps <- length(y)
    } else if (class(networks) == "matrix") {
      # nodal attribute cannot be extracted from a matrix
      stop(paste("If 'networks' contains a 'matrix' object, 'y' must be a", 
          "numeric vector."))
    } else if (class(networks) == "network") {
      # extract nodal attribute and wrap both y and networks in lists
      y <- list(get.vertex.attribute(networks, y))
      networks <- list(as.matrix(networks))
      time.steps <- 1
    } else {
      stop("Data type of the 'networks' argument could not be recognized.")
    }
  } else if (class(y) == "data.frame") {
    time.steps <- ncol(y)
    if (is.null(networks)) {
      # insert NA matrices into new list
      networks <- list()
      for (i in 1:length(y)) {
        mat <- matrix(NA, nrow = nrow(y), ncol = nrow(y))
        rownames(mat) <- rownames(y)
        networks[[i]] <- mat
      }
      # convert y data frame to a list
      l <- list()
      for (i in 1:ncol(y)) {
        if (is.integer(y[, i])) {
          y[, i] <- as.numeric(y[, i])
        } else if (!is.numeric(y[, i])) {
          stop(paste("Column", i, "of argument 'y' is not numeric."))
        }
        l[[i]] <- y[, i]
        names(l[[i]]) <- rownames(y)
      }
      y <- l
    } else if (class(networks) == "list") {
      # check compatibility and convert data frame to list
      if (length(networks) != ncol(y)) {
        stop("'y' and 'networks' have different dimensions.")
      }
      l <- list()
      for (i in 1:ncol(y)) {
        if (!is.numeric(y[, i])) {
          stop(paste("Column", i, "of argument 'y' is not numeric."))
        }
        l[[i]] <- y[, i]
        names(l[[i]]) <- rownames(y)
      }
      y <- l
    } else if (class(networks) == "matrix") {
      # check length and wrap matrix and first y column in list
      networks <- list(networks)
      if (ncol(y) != 1) {
       stop("There is only one network but there are multiple columns in 'y'.")
      }
      l <- list(y[, 1])
      names(l[[1]]) <- rownames(y)
      y <- l
    } else if (class(networks) == "network") {
      # check length and wrap network as matrix and first y column in list
      networks <- list(as.matrix(networks))
      if (ncol(y) != 1) {
       stop("There is only one network but there are multiple columns in 'y'.")
      }
      l <- list(y[, 1])
      names(l[[1]]) <- rownames(y)
      y <- l
    } else {
      stop("Data type of the 'networks' argument could not be recognized.")
    }
  } else if (class(y) == "matrix") {
    stop("The data type of the 'y' argument is 'matrix'. This is not allowed.")
  } else {
    stop("Data type of the 'y' argument could not be recognized.")
  }
  
  # check compatibility of dimensions and presence of labels at each time step
  for (i in 1:time.steps) {
    if (length(y[[i]]) != nrow(as.matrix(networks[[i]]))) {
      stop(paste0("'y' and 'networks' do not have compatible dimensions at ", 
          "t = ", i, ". Please use the 'preprocess' function to make them ", 
          "compatible."))
    }
    if (is.null(names(y[[i]])) && is.null(rownames(as.matrix(networks[[i]])))) {
      stop(paste("The matrices or networks or the 'y' variable must have", 
          "row names, vertex names or names attached to them if multiple", 
          "time points are present."))
    } else if (is.null(names(y[[i]])) && 
        !is.null(rownames(as.matrix(networks[[i]])))) {
      names(y[[i]]) <- rownames(networks[[i]])
    } else if (!is.null(names(y[[i]])) && 
        is.null(rownames(as.matrix(networks[[i]])))) {
      rownames(networks[[i]]) <- names(y[[i]])
    }
    if (!identical(names(y[[i]]), rownames(as.matrix(networks[[i]])))) {
      stop(paste0("The names of 'y' and the row names of 'networks' do not ", 
          "match at t = ", i, "."))
    }
  }
  
  # check 'lag'
  if (!is.numeric(lag)) {
    stop("The 'lag' argument must be numeric.")
  } else if (length(lag) > 1) {
    stop("The 'lag' argument must be of length 1.")
  } else if (lag < 0) {
    stop("The 'lag' argument must be >= 0.")
  } else if (time.steps - lag < 1) {
    if (time.steps == 1) {
      stop(paste("A lag of", lag, "was specified, but there is only", 
          "one time step."))
    } else {
      stop(paste("A lag of", lag, "was specified, but there are only", 
          time.steps, "time steps."))
    }
  }
  n <- lapply(y, length)
  
  objects <- list()
  objects$y <- y
  objects$networks <- networks
  objects$time.steps <- time.steps
  objects$n <- n
  
  # embed node labels
  if (time.steps == 1 && is.null(rownames(objects$networks[[1]]))) {
    nodelabels <- 1:objects$n[[1]]
  } else {
    nodelabels <- character()
    for (i in 1:objects$time.steps) {
      nodelabels <- c(nodelabels, rownames(objects$networks[[i]]))
    }
  }
  objects$nodelabels <- nodelabels
  
  return(objects)
}


# spatial network lag term with optional temporal lag and path distance decay
netlag <- function(y, networks, lag = 0, pathdist = 1, 
    decay = pathdist^-1, normalization = c("no", "row", "column", "complete"), 
    reciprocal = FALSE, center = FALSE, coefname = NULL, ...) {
  
  # check validity of arguments
  if (!is.numeric(pathdist)) {
    stop("'pathdist' must be numeric.")
  }
  if (length(pathdist) != length(decay)) {
    stop("'decay' and 'pathdist' must have the same length.")
  }
  if (any(!is.finite(decay)) || any(!is.finite(pathdist))) {
    stop("'pathdist' and/or 'decay' contain(s) infinite values.")
  }
  if (is.null(reciprocal) || length(reciprocal) > 1) {
    stop("The 'reciprocal' argument must be TRUE or FALSE.")
  }
  
  # get data in the right shape
  objects <- checkDataTypes(y = y, networks = networks, lag = lag)
  
  # do the computations
  results <- list()  # will contain vectors of results for each time step
  for (k in 1:objects$time.steps) {
    pdistmat <- geodist(objects$networks[[k]], ...)$gdist  # path dist mat.
    results[[k]] <- suppressWarnings(netLagCppLoop(objects$networks[[k]], 
        pdistmat, pathdist, decay, objects$y[[k]], normalization[1], 
        reciprocal))  # netlag loop in C++
  }
  
  # convert list of results into data frame and add time and response columns
  response <- numeric()
  time <- numeric()
  spatlag <- numeric()
  for (i in 1:objects$time.steps) {
    time <- c(time, rep(i, objects$n[[i]]))
    response <- c(response, objects$y[[i]])
    if (center == TRUE) {
      results[[i]] <- results[[i]] - mean(results[[i]], na.rm = TRUE)
    }
    spatlag <- c(spatlag, results[[i]])
  }
  
  if (length(lag) == 1 && lag == 0) {
    laglabel <- ""
  } else {
    laglabel <- paste0(".lag", paste(lag, collapse = "."))
  }
  pathdistlabel <- paste0(".pathdist", paste(pathdist, collapse = "."))
  if (any(decay != 1)) {
    decaylabel <- paste0(".decay", paste(decay, collapse = "."))
  } else {
    decaylabel <- ""
  }
  if (normalization[1] == "row") {
    normlabel <- ".rownorm"
  } else if (normalization[1] == "column" || normalization[1] == "col") {
    normlabel <- ".colnorm"
  } else {
    normlabel <- ""
  }
  if (is.null(coefname) || !is.character(coefname) || length(coefname) > 1) {
    coeflabel <- ""
  } else {
    coeflabel <- paste0(".", coefname)
  }
  label <- paste0("netlag", coeflabel, laglabel, pathdistlabel, decaylabel, 
      normlabel)
  dat <- data.frame(spatlag, time = time, node = objects$nodelabels, 
      response = response)
  dat$node <- as.character(dat$node)
  colnames(dat)[1] <- label
  attributes(dat)$lag <- lag
  return(dat)
}


# spatial weighted lag with row and column normalization; for weighted matrices
weightlag <- function(y, networks, lag = 0, normalization = c("no", "row", 
    "column"), center = FALSE, coefname = NULL) {
  
  # check arguments
  if (!is.numeric(lag)) {
    stop("The 'lag' argument must be an integer value or vector of integers.")
  }
  
  # get data in the right shape
  objects <- checkDataTypes(y = y, networks = networks, lag = lag)
  
  # do the computations
  results <- list()  # will contain vectors of results for each time step
  for (k in 1:objects$time.steps) {
    objects$networks[[k]][is.na(objects$networks[[k]])] <- 0
    if (normalization[1] == "row") {
      for (i in 1:nrow(objects$networks[[k]])) {
        rs <- rowSums(objects$networks[[k]])[i]
        for (j in 1:ncol(objects$networks[[k]])) {
          normalized <- objects$networks[[k]][i, j] / rs
          if (is.nan(normalized)) {
            normalized <- 0
          }
          objects$networks[[k]][i, j] <- normalized
        }
      }
    } else if (normalization[1] == "column") {
      for (i in 1:nrow(objects$networks[[k]])) {
        for (j in 1:ncol(objects$networks[[k]])) {
          cs <- colSums(objects$networks[[k]])[j]
          normalized <- objects$networks[[k]][i, j] / cs
          if (is.nan(normalized)) {
            normalized <- 0
          }
          objects$networks[[k]][i, j] <- normalized
        }
      }
    }
    result <- objects$networks[[k]] %*% objects$y[[k]]
    results[[k]] <- result
  }
  
  # convert list of results into data frame and add time and response columns
  response <- numeric()
  time <- numeric()
  spatlag <- numeric()
  for (i in 1:objects$time.steps) {
    time <- c(time, rep(i, objects$n[[i]]))
    response <- c(response, objects$y[[i]])
    if (center == TRUE) {
      results[[i]] <- results[[i]] - mean(results[[i]], na.rm = TRUE)
    }
    spatlag <- c(spatlag, results[[i]])
  }
  
  # create the label
  if (length(lag) == 1 && lag == 0) {
    laglabel <- ""
  } else {
    laglabel <- paste0(".lag", paste(lag, collapse = "."))
  }
  if (normalization[1] == "row") {
    normlabel <- ".rownorm"
  } else if (normalization[1] == "column" || normalization[1] == "col") {
    normlabel <- ".colnorm"
  } else {
    normlabel <- ""
  }
  if (is.null(coefname) || !is.character(coefname) || length(coefname) > 1) {
    coeflabel <- ""
  } else {
    coeflabel <- paste0(".", coefname)
  }
  label <- paste0("weightlag", coeflabel, laglabel, normlabel)
  
  dat <- data.frame(spatlag, time = time, node = objects$nodelabels, 
      response = response)
  dat$node <- as.character(dat$node)
  colnames(dat)[1] <- label
  attributes(dat)$lag <- lag
  return(dat)
}


# outcome variable of other actors * their similarity on another attribute
attribsim <- function(y, attribute, match = FALSE, lag = 0, 
    normalization = c("no", "row", "column"), center = FALSE, coefname = NULL) {
  
  # get data in the right shape
  objects <- checkDataTypes(y = y, networks = NULL, lag = lag)
  
  # check validity of 'attribute' argument and convert
  attrib <- checkDataTypes(y = attribute, networks = NULL, lag = lag)
  if (objects$time.steps != attrib$time.steps) {
    stop("'y' and 'attribute' must have the same number of time steps.")
  }
  if (!identical(attrib$n, objects$n)) {
    stop("'attribute' must have the same number of observations as 'y'.")
  }
  attrib <- attrib$y
  
  # do the computations
  results <- list()
  for (i in 1:objects$time.steps) {
    mat <- matrix(NA, nrow = length(attrib[[i]]), ncol = length(attrib[[i]]))
    for (j in 1:length(attrib[[i]])) {
      for (k in 1:length(attrib[[i]])) {
        if (match == TRUE) {  # node match on the attribute: 1 if both the same
          if (!is.na(attrib[[i]][j]) && !is.na(attrib[[i]][k]) && 
              attrib[[i]][j] == attrib[[i]][k]) {
            mat[j, k] <- 1
          } else {
            mat[j, k] <- 0
          }
        } else if (is.na(attrib[[i]][j]) || is.na(attrib[[i]][k])) {
          mat[j, k] <- NA
        } else {  # absolute dissimilarity
          mat[j, k] <- abs(attrib[[i]][j] - attrib[[i]][k])
        }
      }  # create matrix with absolute differences on the 'attribute'
    }
    mat <- mat / max(mat, na.rm = TRUE)  # standardize
    mat <- 1 - mat  # convert to similarities
    
    # normalization of the similarity matrix
    if (normalization[1] == "row") {
      for (j in 1:nrow(mat)) {
        rs <- rowSums(mat)[j]
        for (k in 1:ncol(mat)) {
          normalized <- mat[j, k] / rs
          if (is.nan(normalized)) {
            normalized <- 0
          }
          mat[j, k] <- normalized
        }
      }
    } else if (normalization[1] == "column") {
      for (j in 1:nrow(mat)) {
        for (k in 1:ncol(mat)) {
          cs <- colSums(mat)[k]
          normalized <- mat[j, k] / cs
          if (is.nan(normalized)) {
            normalized <- 0
          }
          mat[j, k] <- normalized
        }
      }
    }
    
    results[[i]] <- mat %*% objects$y[[i]]  # apply the weights to 'y' vector
  }
  
  # convert list of results into data frame and add time and response columns
  response <- numeric()
  time <- numeric()
  attribsim <- numeric()
  for (i in 1:objects$time.steps) {
    time <- c(time, rep(i, objects$n[[i]]))
    response <- c(response, objects$y[[i]])
    if (center == TRUE) {
      results[[i]] <- results[[i]] - mean(results[[i]], na.rm = TRUE)
    }
    attribsim <- c(attribsim, results[[i]])
  }
  
  # aggregate label
  if (is.null(coefname) || !is.character(coefname) || length(coefname) > 1) {
    coeflabel <- ""
  } else {
    coeflabel <- paste0(".", coefname)
  }
  if (lag == 0) {
    laglabel <- ""
  } else {
    laglabel <- paste0(".lag", paste(lag, collapse = "."))
  }
  label <- paste0("attribsim", coeflabel, laglabel)
  
  # aggregate and return data
  dat <- data.frame(attribsim, time = time, node = objects$nodelabels, 
      response = response)
  dat$node <- as.character(dat$node)
  colnames(dat)[1] <- label
  attributes(dat)$lag <- lag
  return(dat)
}


# interaction effect between multiple terms
interact <- function(x, y, lag = 0, center = FALSE, coefname = NULL) {
  
  # create interaction between x and y
  if (class(x) == "data.frame") {
    first <- x[, 1]
  } else {
    first <- x
  }
  if (class(y) == "data.frame") {
    second <- y[, 1]
  } else {
    second <- y
  }
  z <- first * second
  
  # create time and node variables
  if (class(x) == "data.frame") {
    time <- x[, 2]
    node <- as.character(x[, 3])
  } else if (class(y) == "data.frame") {
    time <- y[, 2]
    node <- as.character(y[, 3])
  } else {
    time <- rep(1, length(z))
  }
  
  # centering
  if (center == TRUE) {
    u <- unique(time)
    for (i in 1:length(u)) {
      z[time == u] <- z[time == u] - mean(z[time == u], na.rm = TRUE)
    }
  }
  
  # create label
  if (is.null(coefname) || !is.character(coefname) || length(coefname) > 1) {
    coeflabel <- ""
  } else {
    coeflabel <- paste0(".", coefname)
  }
  label <- paste0("interaction", coeflabel)
  
  # aggregate results, add response variable if available, and return result
  if (class(x) == "data.frame" && ncol(x) == 4) {
    dat <- data.frame(z, time, node, x[, 4])
    colnames(dat) <- c(label, "time", "node", "response")
  } else if (class(y) == "data.frame" && ncol(y) == 4) {
    dat <- data.frame(z, time, node, y[, 4])
    colnames(dat) <- c(label, "time", "node", "response")
  } else {
    dat <- data.frame(z, time, node)
    colnames(dat) <- c(label, "time", "node")
  }
  dat$node <- as.character(dat$node)
  if (is.null(lag)) {
    lag <- union(attributes(x)$lag, attributes(y)$lag)
  } else if (!is.numeric(lag[1])) {
    stop("'lag' must be a single integer value.")
  } else if (attributes(x)$lag != 0 || attributes(x)$lag != 0) {
    message(paste("Using the 'lag' argument from the interaction term, not", 
        "from the 'x' or 'y' arguments handed over to the interaction term."))
  }
  if (length(lag) > 1) {
    warning(paste("Interaction term: several 'lag' arguments are provided.", 
        "Only the first one is retained."))
  }
  if (is.null(attributes(dat)$lag)) {
    attributes(dat)$lag <- 0
  }
  return(dat)
}


# simple covariate or temporally lagged covariate
covariate <- function(y, lag = 0, exponent = 1, center = FALSE, 
    coefname = NULL) {
  
  # check validity of arguments
  if (is.null(lag) || !is.numeric(lag) || floor(lag) != lag) {
    stop("'lag' should be an integer value or a vector of integers.")
  }
  if (is.null(exponent)) {
    exponent <- 1
  } else if (length(exponent) > 1 || !is.numeric(exponent)) {
    stop("'exponent' should contain a single numeric value.")
  }
  objects <- checkDataTypes(y = y, networks = NULL, lag = lag)
  
  response <- numeric()
  time <- numeric()
  y <- numeric()
  for (i in 1:objects$time.steps) {
    time <- c(time, rep(i, objects$n[[i]]))
    response <- c(response, objects$y[[i]])
    if (center == TRUE) {
      result <- objects$y[[i]]^exponent - mean(objects$y[[i]]^exponent, 
          na.rm = TRUE)
    } else {
      result <- objects$y[[i]]^exponent
    }
    y <- c(y, result)
  }
  
  # aggregate label
  if (is.null(coefname) || !is.character(coefname) || length(coefname) > 1) {
    coeflabel <- ""
  } else {
    coeflabel <- paste0(".", coefname)
  }
  if (length(lag) == 1 && lag == 0) {
    laglabel <- ""
  } else {
    laglabel <- paste0(".lag", paste(lag, collapse = "."))
  }
  label <- paste0("covariate", coeflabel, laglabel)
  dat <- data.frame(y, time = time, node = objects$nodelabels, 
      response = response)
  dat$node <- as.character(dat$node)
  colnames(dat)[1] <- label
  attributes(dat)$lag <- lag
  return(dat)
}


# structural similarity term with optional temporal lag
structsim <- function(y, networks, lag = 0, method = c("euclidean", 
    "minkowski", "jaccard", "binary", "hamming"), center = FALSE, 
    coefname = NULL, ...) {
  
  # get data in the right shape
  objects <- checkDataTypes(y = y, networks = networks, lag = lag)
  
  # do the computations
  results <- list()  # will contain vectors of results for each time step
  for (k in 1:objects$time.steps) {
    # compute structural similarity matrix
    if (method[1] == "euclidean" || method[1] == "minkowski") {
      d <- as.matrix(dist(objects$networks[[k]], method = method[1], ...))
      mx <- max(d, na.rm = TRUE)
      s <- (mx - d) / mx  # convert dist to similarity and standardize [0; 1]
    } else if (method[1] == "binary") {
      d <- as.matrix(dist(objects$networks[[k]], method = method[1], ...))
      s <- 1 - d
    } else if (method[1] == "jaccard") {
      d <- as.matrix(vegdist(objects$networks[[k]], method = "jaccard", 
          na.rm = TRUE, ...))
      s <- 1 - d
    } else if (method[1] == "hamming") {
      d <- sedist(objects$networks[[k]], method = "hamming", ...) / 
          nrow(objects$networks[[k]])  # standardize to [0; 1]
      s <- 1 - d
    } else {
      stop("'method' argument was not recognized.")
    }
    diag(s) <- 0
    rm(d)
    
    # apply structural similarities as weight matrix
    s[is.na(s)] <- 0
    result <- s %*% objects$y[[k]]
    results[[k]] <- result
  }
  
  # convert list of results into data frame and add time and response columns
  response <- numeric()
  time <- numeric()
  structsim <- numeric()
  for (i in 1:objects$time.steps) {
    time <- c(time, rep(i, objects$n[[i]]))
    response <- c(response, objects$y[[i]])
    if (center == TRUE) {
      results[[i]] <- results[[i]] - mean(results[[i]], na.rm = TRUE)
    }
    structsim <- c(structsim, results[[i]])
  }
  
  # aggregate label
  if (is.null(coefname) || !is.character(coefname) || length(coefname) > 1) {
    coeflabel <- ""
  } else {
    coeflabel <- paste0(".", coefname)
  }
  if (lag == 0) {
    laglabel <- ""
  } else {
    laglabel <- paste0(".lag", paste(lag, collapse = "."))
  }
  label <- paste0("structsim", coeflabel, laglabel, ".", method[1])
  
  # aggregate and return data
  dat <- data.frame(structsim, time = time, node = objects$nodelabels, 
      response = response)
  dat$node <- as.character(dat$node)
  colnames(dat)[1] <- label
  attributes(dat)$lag <- lag
  return(dat)
}


# model term which indicates whether an actor has a certain degree centrality
centrality <- function(networks, type = c("indegree", "outdegree", "freeman", 
    "betweenness", "flow", "closeness", "eigenvector", "information", "load", 
    "bonpow"), directed = TRUE, lag = 0, rescale = FALSE, center = FALSE, 
    coefname = NULL, ...) {
  
  # check validity of arguments and prepare data
  if (is.null(directed) || !is.logical(directed)) {
    stop("'directed' must be TRUE or FALSE.")
  } else if (length(directed) != 1) {
    stop("The 'directed' argument must contain a single logical value only.")
  } else if (directed == FALSE) {
    gmode <- "graph"
  } else {
    gmode <- "digraph"
  }
  objects <- checkDataTypes(y = NULL, networks = networks, lag = lag)
  
  # do the computations
  centlist <- list()
  for (i in 1:objects$time.steps) {
    if (type[1] == "indegree") {
      cent <- degree(objects$networks[[i]], gmode = gmode, cmode = "indegree", 
          rescale = rescale, ...)
    } else if (type[1] == "outdegree") {
      cent <- degree(objects$networks[[i]], gmode = gmode, cmode = "outdegree", 
          rescale = rescale, ...)
    } else if (type[1] == "freeman") {
      cent <- degree(objects$networks[[i]], gmode = gmode, cmode = "freeman", 
          rescale = rescale, ...)
    } else if (type[1] == "betweenness") {
      cent <- betweenness(objects$networks[[i]], gmode = gmode, 
          rescale = rescale, ...)
    } else if (type[1] == "flow") {
      cent <- flowbet(objects$networks[[i]], gmode = gmode, rescale = rescale, 
          ...)
    } else if (type[1] == "closeness") {
      cent <- closeness(objects$networks[[i]], gmode = gmode, 
          rescale = rescale, ...)
    } else if (type[1] == "eigenvector") {
      cent <- evcent(objects$networks[[i]], gmode = gmode, rescale = rescale, 
          ...)
    } else if (type[1] == "information") {
      cent <- infocent(objects$networks[[i]], gmode = gmode, 
          rescale = rescale, ...)
    } else if (type[1] == "load") {
      cent <- loadcent(objects$networks[[i]], gmode = gmode, 
          rescale = rescale, ...)
    } else if (type[1] == "bonpow") {
      cent <- bonpow(objects$networks[[i]], gmode = gmode, 
          rescale = rescale, tol = 1e-20, ...)
    } else {
      stop("'type' argument was not recognized.")
    }
    centlist[[i]] <- cent
  }
  
  # aggregate data frame
  time <- numeric()
  y <- numeric()
  for (i in 1:objects$time.steps) {
    time <- c(time, rep(i, objects$n[[i]]))
    if (is.null(centlist[[i]])) {
      y <- c(y, rep(NA, objects$n[[i]]))
    } else {
      if (center == TRUE) {
        centlist[[i]] <- centlist[[i]] - mean(centlist[[i]], na.rm = TRUE)
      }
      y <- c(y, centlist[[i]])
    }
  }
  
  # aggregate label and results and return them
  if (is.null(coefname) || !is.character(coefname) || length(coefname) > 1) {
    coeflabel <- ""
  } else {
    coeflabel <- paste0(".", coefname)
  }
  if (lag == 0) {
    laglabel <- ""
  } else {
    laglabel <- paste0(".lag", paste(lag, collapse = "."))
  }
  label <- paste0(type[1], coeflabel, laglabel)
  dat <- data.frame(y, time = time, node = objects$nodelabels)
  dat$node <- as.character(dat$node)
  colnames(dat)[1] <- label
  attributes(dat)$lag <- lag
  return(dat)
}


# model term which indicates whether an actor has a certain degree centrality
degreedummy <- function(networks, deg = 0, type = c("indegree", "outdegree", 
    "freeman"), reverse = FALSE, directed = TRUE, lag = 0, center = FALSE, 
    coefname = NULL, ...) {
  
  # check validity of arguments and prepare data
  if (is.null(directed) || !is.logical(directed)) {
    stop("'directed' must be TRUE or FALSE.")
  } else if (length(directed) != 1) {
    stop("The 'directed' argument must contain a single logical value only.")
  } else if (directed == FALSE) {
    gmode <- "graph"
  } else {
    gmode <- "digraph"
  }
  if (is.null(deg) || !is.numeric(deg) || any(deg < 0)) {
    stop("'deg' must be a numeric value or vector of positive integers.")
  }
  objects <- checkDataTypes(y = NULL, networks = networks, lag = lag)
  
  # do the computations
  dummylist <- list()
  for (i in 1:objects$time.steps) {
    d <- degree(objects$networks[[i]], gmode = gmode, cmode = type[1], ...)
    d <- 1 * (d %in% deg)
    if (reverse == TRUE) {
      d <- d * -1 + 1
    }
    dummylist[[i]] <- d
  }
  
  # aggregate data frame
  time <- numeric()
  y <- numeric()
  for (i in 1:objects$time.steps) {
    time <- c(time, rep(i, objects$n[[i]]))
    if (is.null(dummylist[[i]])) {
      y <- c(y, rep(NA, objects$n[[i]]))
    } else {
      if (center == TRUE) {
        dummylist[[i]] <- dummylist[[i]] - mean(dummylist[[i]], na.rm = TRUE)
      }
      y <- c(y, dummylist[[i]])
    }
  }
  
  # aggregate label and return data
  if (is.null(coefname) || !is.character(coefname) || length(coefname) > 1) {
    coeflabel <- ""
  } else {
    coeflabel <- paste0(".", coefname)
  }
  if (lag == 0) {
    laglabel <- ""
  } else {
    laglabel <- paste0(".lag", paste(lag, collapse = "."))
  }
  degreelabel <- paste0(".degree", paste(deg, collapse = "."))
  label <- paste0(type[1], coeflabel, degreelabel, laglabel)
  dat <- data.frame(y, time = time, node = objects$nodelabels)
  dat$node <- as.character(dat$node)
  colnames(dat)[1] <- label
  attributes(dat)$lag <- lag
  return(dat)
}


# local clustering coefficient
clustering <- function(networks, directed = TRUE, lag = 0, center = FALSE, 
    coefname = NULL, ...) {
  
  # prepare arguments and data
  if (directed == TRUE) {
    mode <- "directed"
  } else {
    mode <- "undirected"
  }
  objects <- checkDataTypes(y = NULL, networks = networks, lag = lag)
  
  # do the computations
  results <- list()
  for (i in 1:objects$time.steps) {
    # local clustering coefficient from the igraph package
    g <- graph.adjacency(objects$networks[[i]], mode = mode)
    trans <- transitivity(g, type = "local", isolates = "zero")
    results[[i]] <- trans
  }
  
  # convert list of results into data frame and add time column
  time <- numeric()
  lcc <- numeric()
  for (i in 1:objects$time.steps) {
    time <- c(time, rep(i, objects$n[[i]]))
    if (is.null(results[[i]])) {
      lcc <- c(lcc, rep(NA, objects$n[[i]]))
    } else {
      if (center == TRUE) {
        results[[i]] <- results[[i]] - mean(results[[i]], na.rm = TRUE)
      }
      lcc <- c(lcc, results[[i]])
    }
  }
  
  # aggregate label and return results
  if (is.null(coefname) || !is.character(coefname) || length(coefname) > 1) {
    coeflabel <- ""
  } else {
    coeflabel <- paste0(".", coefname)
  }
  if (lag == 0) {
    laglabel <- ""
  } else {
    laglabel <- paste0(".lag", paste(lag, collapse = "."))
  }
  label <- paste0("clustering", coeflabel, laglabel)
  dat <- data.frame(lcc, time = time, node = objects$nodelabels)
  dat$node <- as.character(dat$node)
  colnames(dat)[1] <- label
  attributes(dat)$lag <- lag
  return(dat)
}


# spatial lag for k-clique co-members
cliquelag <- function(y, networks, k.min = 2, k.max = Inf, directed = TRUE, 
    lag = 0, normalization = c("no", "row", "column"), center = FALSE, 
    coefname = NULL) {
  
  # check arguments
  if (!is.numeric(lag)) {
    stop("The 'lag' argument must be an integer value or vector of integers.")
  }
  if (is.null(k.min) || !is.numeric(k.min) || length(k.min) > 1) {
    stop("The 'k.min' argument must be a single numeric value.")
  }
  if (is.null(k.max) || !is.numeric(k.max) || length(k.max) > 1) {
    stop("The 'k.max' argument must be a single numeric value.")
  }
  if (k.max < k.min) {
    stop("'k.min' must be smaller than 'k.max'.")
  }
  if (floor(k.min) != k.min || floor(k.max) != k.max) {
    stop("The 'k.min' and 'k.max' arguments should be integer values.")
  }
  if (directed == TRUE) {
    mode <- "digraph"
  } else {
    mode <- "graph"
  }
  
  # get data in the right shape
  objects <- checkDataTypes(y = y, networks = networks, lag = lag)
  
  # do the computations
  results <- list()  # will contain vectors of results for each time step
  for (k in 1:objects$time.steps) {
    objects$network[[k]][is.na(objects$network[[k]])] <- 0  # correct NAs in nw
    objects$y[[k]][is.na(objects$y[[k]])] <- 0  # correct NAs in y
    
    # retrieve clique comembership matrices by clique size
    w3d <- clique.census(objects$networks[[k]], mode = mode, 
        clique.comembership = "bysize", tabulate.by.vertex = FALSE, 
        enumerate = FALSE)$clique.comemb
    k.max <- min(c(k.max, dim(w3d)[1]))
    w <- matrix(0, nrow = nrow(objects$networks[[k]]), 
        ncol = ncol(objects$networks[[k]]))
    for (i in k.min:k.max) {  # sum up three-cliques, four-cliques etc.
      w <- w + w3d[i, , ]
    }
    diag(w) <- 0  # diagonal is not valid (self-cliques...)
    if (normalization[1] == "row") {  # row norm. = average clique alter effect
      for (i in 1:nrow(w)) {
        rs <- rowSums(w)[i]
        for (j in 1:ncol(w)) {
          normalized <- w[i, j] / rs
          if (is.nan(normalized)) {
            normalized <- 0
          }
          w[i, j] <- normalized
        }
      }
    } else if (normalization[1] == "column") {
      for (i in 1:nrow(w)) {
        for (j in 1:ncol(w)) {
          cs <- colSums(w)[j]
          normalized <- w[i, j] / cs
          if (is.nan(normalized)) {
            normalized <- 0
          }
          w[i, j] <- normalized
        }
      }
    }
    results[[k]] <- w %*% objects$y[[k]]
  }
  
  # convert list of results into data frame and add time and response columns
  response <- numeric()
  time <- numeric()
  cliquelag <- numeric()
  for (i in 1:objects$time.steps) {
    time <- c(time, rep(i, objects$n[[i]]))
    response <- c(response, objects$y[[i]])
    if (center == TRUE) {
      results[[i]] <- results[[i]] - mean(results[[i]], na.rm = TRUE)
    }
    cliquelag <- c(cliquelag, results[[i]])
  }
  
  # create the label
  if (length(lag) == 1 && lag == 0) {
    laglabel <- ""
  } else {
    laglabel <- paste0(".lag", paste(lag, collapse = "."))
  }
  if (normalization[1] == "row") {
    normlabel <- ".rownorm"
  } else if (normalization[1] == "column" || normalization[1] == "col") {
    normlabel <- ".colnorm"
  } else {
    normlabel <- ""
  }
  klabel <- paste0(".k.", k.min, ".", k.max)
  if (is.null(coefname) || !is.character(coefname) || length(coefname) > 1) {
    coeflabel <- ""
  } else {
    coeflabel <- paste0(".", coefname)
  }
  label <- paste0("cliquelag", coeflabel, klabel, laglabel, normlabel)
  
  dat <- data.frame(cliquelag, time = time, node = objects$nodelabels, 
      response = response)
  dat$node <- as.character(dat$node)
  colnames(dat)[1] <- label
  attributes(dat)$lag <- lag
  return(dat)
}


# function which aggregates data for glm analysis
tnamdata <- function(formula, center.y = FALSE) {
  
  # parse the formula
  if (class(formula) != "formula") {
    stop("'formula' must be a formula object.")
  }
  lhs <- deparse(formula[[2]])  # name of the response variable
  lhs <- eval(parse(text = lhs))  # get the actual response data
  rhs <- paste0(deparse(formula[[3]]), collapse = "")  # rhs of formula
  rhs <- gsub("\\s+", " ", rhs)  # get rid of redundant spaces
  rhs <- strsplit(rhs, " \\+ ")[[1]]  # parse separate formula elements
  
  # create data frame with response variable, time, and nodes
  time <- numeric()
  node <- character()
  response <- numeric()
  if (class(lhs) == "list") {
    for (i in 1:length(lhs)) {
      if (!is.numeric(lhs[[i]])) {
        stop(paste("The response variable should be numeric or a list of", 
            "numerics or a data frame with one time point per column."))
      }
      if (is.null(names(lhs[[i]])) || length(names(lhs[[i]])) != 
          length(lhs[[i]])) {
        stop(paste("The outcome variable must have node labels if multiple", 
            "time points are present."))
      }
      node <- c(node, names(lhs[[i]]))
      time <- c(time, rep(i, length(lhs[[i]])))
      if (center.y == TRUE) {
        lhs[[i]] <- lhs[[i]] - mean(lhs[[i]], na.rm = TRUE)
      }
      response <- c(response, lhs[[i]])
    }
  } else if (class(lhs) == "data.frame") {
    for (i in 1:ncol(lhs)) {
      if (!is.numeric(lhs[, i])) {
        stop(paste("The response variable should be numeric or a list of", 
            "numerics or a data frame with one time point per column."))
      }
      if (is.null(rownames(lhs)) || length(rownames(lhs)) != nrow(lhs)) {
        stop(paste("The outcome variable must have node labels if multiple", 
            "time points are present."))
      }
      node <- c(node, rownames(lhs))
      time <- c(time, rep(i, nrow(lhs)))
      if (center.y == TRUE) {
        lhs[, i] <- lhs[, i] - mean(lhs[, i], na.rm = TRUE)
      }
      response <- c(response, lhs[, i])
    }
  } else if (!is.numeric(lhs)) {
    stop("Data type of the response variable could not be recognized.")
  } else {
    response <- lhs
    if (center.y == TRUE) {
      response <- response - mean(response, na.rm = TRUE)
    }
    time <- rep(1, length(lhs))
    node <- as.character(1:length(lhs))
  }
  dat <- data.frame(response = response, time = time, node = node)
  
  # compute results according to rhs
  resultlist <- list()
  for (i in 1:length(rhs)) {
    result <- eval(parse(text = rhs[i]))
    resultlist[[i]] <- result
  }
  
  # check compatibility of labels
  for (i in 1:length(resultlist)) {
    for (j in 1:length(resultlist)) {
      itime <- length(unique(resultlist[[i]]$time))
      jtime <- length(unique(resultlist[[j]]$time))
      if ((itime > 1 || jtime > 1) && i < j) {
        inters <- length(intersect(resultlist[[i]]$node, resultlist[[j]]$node))
        if (inters == 0) {
          stop(paste("Model terms", i, "and", j, "do not have any", 
              "intersecting node labels. Please attach names, row names, or", 
              "vertex names to the 'y' or 'networks' argument."))
        }
      }
    }
  }
  
  # take care of the lags
  for (i in 1:length(resultlist)) {
    lag.i <- attributes(resultlist[[i]])$lag
    if (is.null(lag.i) || length(lag.i) == 0) {
      lag.i <- 0
    }
    resultlist[[i]]$time <- resultlist[[i]]$time + lag.i
  }
  
  # merge results with response variable and take care of lags
  for (i in 1:length(resultlist)) {
    dat <- merge(dat, resultlist[[i]], by = c("time", "node"), all.x = TRUE, 
        all.y = FALSE)
    colnames(dat)[3] <- "response"
    dat$node <- as.character(dat$node)
    if (ncol(resultlist[[i]]) == 4) {
      dat <- dat[, -ncol(dat)]
    }
  }
  dat <- dat[, c(3, 1, 2, 4:ncol(dat))]
  
  return(dat)
}


# temporal network autocorrelation model
tnam <- function(formula, family = gaussian, re.node = FALSE, 
    re.time = FALSE, time.linear = FALSE, time.quadratic = FALSE, 
    center.y = FALSE, na.action = na.omit, ...) {
  
  # prepare the data frame
  dat <- tnamdata(formula, center.y = center.y)
  
  # check if GLM is appropriate
  if (re.node == FALSE && re.time == FALSE && length(unique(dat$time)) > 1) {
    warning(paste("Different time points are available. You might want to use", 
        "a mixed effects model using arguments 're.time' and/or 're.node'."))
  }
  
  # take care of the node variable: keep as random effect or remove
  if (re.node == TRUE && length(unique(dat$time)) > 1) {
    glmest.node <- FALSE
  } else {
    dat <- dat[, -3]
    glmest.node <- TRUE
  }
  
  # take care of the time variable
  if (time.linear == TRUE && time.quadratic == TRUE && re.time == TRUE) {
    # T-T-T
    if (length(unique(dat$time)) > 1) {
      glmest.time <- FALSE
      dat$re.time <- dat$time  # add RE
      dat <- dat[, -2]  # remove linear effect
      warning(paste("Arguments 're.time' and 'time.linear' cannot be used", 
          "together. Omitting the linear time effect."))
      warning(paste("Arguments 're.time' and 'time.quadratic' cannot be used", 
          "together. Omitting the quadratic time effect."))
    } else {
      glmest.time <- TRUE
      dat <- dat[, -2]  # remove linear effect
      message("Time effects are ignored because only one time step is present.")
    }
  } else if (time.linear == TRUE && time.quadratic == TRUE && 
      re.time == FALSE) {
    # T-T-F
    glmest.time <- TRUE
    if (length(unique(dat$time)) > 1) {
      dat$time.squared <- dat$time^2
    } else {
      dat <- dat[, -2]  # remove linear effect
      message("Time effects are ignored because only one time step is present.")
    }
  } else if (time.linear == TRUE && time.quadratic == FALSE && 
      re.time == FALSE) {
    # T-F-F
    glmest.time <- TRUE
    if (length(unique(dat$time)) > 1) {
      # OK; do not modify anything
    } else {
      dat <- dat[, -2]  # remove linear effect
      message("Time effects are ignored because only one time step is present.")
    }
  } else if (time.linear == FALSE && time.quadratic == FALSE && 
      re.time == FALSE) {
    # F-F-F
    dat <- dat[, -2]  # remove linear effect
    glmest.time <- TRUE
  } else if (time.linear == FALSE && time.quadratic == TRUE && 
      re.time == FALSE) {
    # F-T-F
    glmest.time <- TRUE
    if (length(unique(dat$time)) > 1) {
      dat$time.squared <- dat$time^2  # create quadratic effect
    } else {
      message("Time effects are ignored because only one time step is present.")
    }
    dat <- dat[, -2]  # remove linear effect
  } else if (time.linear == FALSE && time.quadratic == TRUE && 
      re.time == TRUE) {
    # F-T-T
    if (length(unique(dat$time)) > 1) {
      glmest.time <- FALSE
      dat$re.time <- dat$time  # add RE
      message(paste("Arguments 're.time' and 'time.quadratic' cannot be used", 
          "together. Omitting the quadratic time effect."))
    } else {
      glmest.time <- TRUE
      message("Time effects are ignored because only one time step is present.")
    }
    dat <- dat[, -2]  # remove linear effect
  } else if (time.linear == FALSE && time.quadratic == FALSE && 
      re.time == TRUE) {
    # F-F-T
    if (length(unique(dat$time)) > 1) {
      glmest.time <- FALSE
      dat$re.time <- dat$time  # add RE
    } else {
      glmest.time <- TRUE
      message("Time effects are ignored because only one time step is present.")
    }
    dat <- dat[, -2]  # remove linear effect
  } else if (time.linear == TRUE && time.quadratic == FALSE && 
      re.time == TRUE) {
    # T-F-T
    if (length(unique(dat$time)) > 1) {
      glmest.time <- FALSE
      dat$re.time <- dat$time  # add RE
      dat <- dat[, -2]  # remove linear effect
      warning(paste("Arguments 're.time' and 'time.linear' cannot be used", 
          "together. Omitting the linear time effect."))
    } else {
      glmest.time <- TRUE
      dat <- dat[, -2]  # remove linear effect
      message("Time effects are ignored because only one time step is present.")
    }
  }
  
  if (glmest.node == FALSE || glmest.time == FALSE) {
    glmest <- FALSE
  } else {
    glmest <- TRUE
  }
  
  # estimate!
  if (glmest == TRUE) {  # GLM is necessary; no random effects
    model <- glm(dat, family = family, na.action = na.action, ...)
  } else {  # mixed-effects model (lme4) is required; random effects present
    if (is.character(family)) {
      family <- get(family, mode = "function", envir = parent.frame(2))
    } else if (is.function(family)) {
      family <- family()
    }
    if (isTRUE(all.equal(family, gaussian()))) {  # gaussian link: use lmer
      if (re.node == TRUE && re.time == TRUE) {
        model <- lme4::lmer(response ~ . - re.time - node + (1|re.time) + 
            (1|node), data = dat, na.action = na.action, ...)
      } else if (re.node == TRUE && re.time == FALSE) {
        model <- lme4::lmer(response ~ . - node + (1|node), data = dat, 
            na.action = na.action, ...)
      } else if (re.node == FALSE && re.time == TRUE) {
        model <- lme4::lmer(response ~ . -re.time + (1|re.time), data = dat, 
            na.action = na.action, ...)
      }
    } else {
      if (re.node == TRUE && re.time == TRUE) {  # other link function: glmer
        model <- lme4::glmer(response ~ . - re.time - node + (1|re.time) + 
            (1|node), data = dat, family = family, na.action = na.action, ...)
      } else if (re.node == TRUE && re.time == FALSE) {
        model <- lme4::glmer(response ~ . - node + (1|node), data = dat, 
            family = family, na.action = na.action, ...)
      } else if (re.node == FALSE && re.time == TRUE) {
        model <- lme4::glmer(response ~ . - re.time + (1|re.time), data = dat, 
            family = family, na.action = na.action, ...)
      }
    }
  }
  
  return(model)
}

