# 11.1.21 

#modifying tnamdata_mod to incorporate time
# 11.2.31- X has been modified and is the correct output. 
# note: created separate function to test output. The existing code + examples will still work with tnamdata_mod 


## modifying tnamdata_mod (ha)
tnamdata_mod_t <- function(formula, time_steps, center.y = FALSE) {
  print(formula)
  # parse the formula
  if (class(formula) != "formula") { 
    stop("'formula' must be a formula object.")
  }
  lhs <- deparse(formula[[2]])  # name of the response variable
  lhs <- eval(parse(text = lhs))  # get the actual response data
  rhs <- paste0(deparse(formula[[3]]), collapse = "")  # rhs of formula
  rhs <- gsub("\\s+", " ", rhs)  # get rid of redundant spaces
  rhs <- strsplit(rhs, " \\+ ")[[1]]  # parse separate formula elements
  rhs_split<- strsplit(rhs, '[(]')
  exog_list<- c("centrality", "clustering", "degreedummy", "interact", "covariate")
  endog_list<- c("attribsim", "cliquelag", "netlag", "weightlag", "structsim")
  W_listr<- c("W", "W_t") #add in term for just adding in full W matrices 
  rhs_exog<- vector()
  rhs_endog<- vector() 
  rhs_w<- vector() 
  
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
  
  #now deal w/ RHS
  for (i in 1:length(rhs)) {
    #print(rhs_split[[i]][1])
    if (rhs_split[[i]][1] %in% exog_list) {
      rhs_exog <- c(rhs_exog, rhs[i])
    } else if (rhs_split[[i]][1] %in% endog_list) { 
      rhs_endog <- c(rhs_endog, rhs[i])
    } else if (rhs_split[[i]][1] %in% W_listr) {
      rhs_w <- c(rhs_w, rhs[i])
    } else { 
      #print(rhs_split[[i]][1])
      print("you inputted an unspecified model term") #could add optional exog, endog for alt model terms 
    }
  }
  
  #build X (from rhs_exog)
  exog_resultlist <- list()
  for (i in 1:length(rhs_exog)) {
    result <- eval(parse(text = rhs_exog[i])) #print(result) -- matches here 
    exog_resultlist[[i]] <- result
  }

  for (i in 1:length(exog_resultlist)) {
    dat <- merge(dat, exog_resultlist[[i]], by = c("time", "node"), all.x = TRUE, 
                 all.y = FALSE, sort= FALSE)
    colnames(dat)[3] <- "response"
    dat$node <- as.character(dat$node)
    if (ncol(exog_resultlist[[i]]) == 4) {
      dat <- dat[, -ncol(dat)]
    }
  }
  X <- dat[, c(1, 3, 2, 4:ncol(dat))]
  
  #change data format! 
  X<- pivot_wider(data = X, 
                       id_cols = node, 
                       names_from = time, 
                       values_from = colnames(X)[4:ncol(X)])
  X<- as.data.frame(X)
  X<- X[,2:ncol(X)] #just the covariate terms...(for now -- don't need node label)
  
  #build W.list (from rhs_endog)
  W.list<- vector()
  for (i in 0:length(rhs_endog)) { 
    W.list<- c(W.list, eval(parse(text = rhs_endog[i])))
  }
  
  ## build actual W.list from W 
  #W.list<- vector()
  temp<- vector()
  if (length(rhs_w)!=0){
    for (i in 1:length(rhs_w)) { 
      if ((rhs_w[i]) == "W"){
        result <- eval(parse(text = rhs_w[i]))
        temp<- append(temp, list(result))
      }else { 
        result <- eval(parse(text = rhs_w[i]))
        for (i in 1:length(result)) {
          temp<- append(temp, list(result[[i]]))
        }
      }
    }
  }
  
  dat1 <- data.frame(response = response, time = time, node = node)
  
  ##change data structure of y to wide 
  Y<- pivot_wider(data = dat, 
              id_cols = node, 
              names_from = time, 
              values_from = response)
  
  Y<- as.data.frame(Y[,2:ncol(Y)])
  for (i in 1:time_steps) {
    names(Y)[i] <- paste("t", i, sep="")
  }

  #change data structure on W.list to list of lists, label 
  W.list<- append(temp,W.list)
  
  lili<- list() ##list of lists
  first_index<- 0
  for (i in 1:(length(W.list)/time_steps)) {
    li <- list() #list 
    for (j in 1:3) {
      li[[j]]<- W.list[[first_index + 1]] 
      first_index = first_index + 1
      names(li)[j] <- paste("t", j, sep="")
    }
    lili[[i]]<- li
    names(lili)[i]<-paste("W", i, sep="") ##this is where you'd change to have it say friendship, or other W inputs rather than W labelling
  }
  
  return (list("y"=Y,"X"= X, "W.list"= lili))
}