checkAdjacency <- function(nodes,graph){
  # Function to check for adjacency in the directed graph, to be used with apply
  v1 <- nodes[1]
  v2 <- nodes[2]
  if(v1 %in% parents(graph,v2) | v2 %in% parents(graph,v1)){
    return(TRUE)
  }else{
    return(FALSE)
  }
}

createFirstTree <- function(dagX,ordering,alpha_1,alpha_high,alpha_low,beta,gamma){
  graphWeights <- matrix(0,nrow=length(ordering),ncol=length(ordering))
  colnames(graphWeights) <- ordering.igraph.top
  rownames(graphWeights) <- ordering.igraph.top
  for(i in ordering){
    # define W as parent set of i
    W <- parents(dagX,i)
    if(length(W)==0){
      # do nothing
    }else{
      if(length(W)==1){
        graphWeights[i,W] <- graphWeights[i,W] + alpha_1
        graphWeights[W,i] <- graphWeights[W,i] + alpha_1
      }else{
        combW <- combn(W,2)
        beta_sets <- which(apply(combW,2,checkAdjacency,graph=dagX)==TRUE)
        gamma_sets <- which(apply(combW,2,checkAdjacency,graph=dagX)==FALSE)
        for(l in beta_sets){ # zwischen verbundenen Eltern werden betas gesetzt
          graphWeights[combW[1,l],combW[2,l]] <- graphWeights[combW[1,l],combW[2,l]] + beta
          graphWeights[combW[2,l],combW[1,l]] <- graphWeights[combW[2,l],combW[1,l]] + beta
          if(abs(TauX[combW[1,l],i])>abs(TauX[combW[2,l],i])){ # bei verbundenen Eltern werden die Kendall's tau Werte berÃ¼cksichtigt da sonst Kreise begÃ¼nstigt werden
            graphWeights[combW[1,l],i] <- graphWeights[combW[1,l],i] + alpha_high
            graphWeights[i,combW[1,l]] <- graphWeights[i,combW[1,l]] + alpha_high
            graphWeights[combW[2,l],i] <- graphWeights[combW[2,l],i] + alpha_low
            graphWeights[i,combW[2,l]] <- graphWeights[i,combW[2,l]] + alpha_low
          }else{
            graphWeights[combW[1,l],i] <- graphWeights[combW[1,l],i] + alpha_low
            graphWeights[i,combW[1,l]] <- graphWeights[i,combW[1,l]] + alpha_low
            graphWeights[combW[2,l],i] <- graphWeights[combW[2,l],i] + alpha_high
            graphWeights[i,combW[2,l]] <- graphWeights[i,combW[2,l]] + alpha_high
          }
        }
        for(m in gamma_sets){ # bei nicht verbundenen Eltern werden zwischen i und den v-structures gammas gesetzt
          graphWeights[combW[1,m],i] <- graphWeights[combW[1,m],i] + gamma
          graphWeights[combW[2,m],i] <- graphWeights[combW[2,m],i] + gamma
          graphWeights[i,combW[1,m]] <- graphWeights[i,combW[1,m]] + gamma
          graphWeights[i,combW[2,m]] <- graphWeights[i,combW[2,m]] + gamma
        }
      }
    }
  }
  graphAdjMatrix <- graph_from_adjacency_matrix(graphWeights,mode="undirected",weighted=TRUE)
  while(!is.connected(graphAdjMatrix)){
    d <- length(ordering)
    dsample <- sample(d,2)
    graphWeights[dsample[1],dsample[2]] <- graphWeights[dsample[1],dsample[2]] + 0.01
    graphWeights[dsample[2],dsample[1]] <- graphWeights[dsample[2],dsample[1]] + 0.01
    graphAdjMatrix <- graph_from_adjacency_matrix(graphWeights,mode="undirected",weighted=TRUE)
  }
  E(graphAdjMatrix)$weight <- -E(graphAdjMatrix)$weight
  g <- minimum.spanning.tree(graphAdjMatrix)
  ### transform first tree to RVM tree object - code borrowed from RVineStructureSelect (credits...)
  E(g)$name <- paste(as_edgelist(g)[, 1], as_edgelist(g)[, 2], sep = ",")
  for (i in 1:gsize(g)) {
    E(g)$conditionedSet[[i]] <- ends(g, i, names = TRUE)
  }
  oldTreeMatrix <- as.matrix(get.adjacency(g,attr="weight"))
  newWeights <- graphWeights + oldTreeMatrix # funktioniert da hier die negativen Werte die positiven aufheben wo eine edge vorliegt, bei kleinerdimensionierten graphen funktioniert das nicht mehr
  return(list(firstTree=g,graph=graphAdjMatrix,newWeights=newWeights,firstWeights=graphWeights))
}

checkProximity <- function(y,RVS,k,x){
  ### Function evaluates if placing the value y in the level-th row of the RVine Structure Matrix RVS in column of entry x
  ### is admissible - note that k==1 is d-th row
  ### Code is not intended to check consistency of entire R-Vine matrices, however it may check the validity of recently
  ### added rows and colums starting from the d-th row
  codeT1 <- codeT2 <- FALSE
  col.x <- which(diag(RVS)==x)
  d <- length(diag(RVS))
  if(y < 1 || y > d){
    return(FALSE)
  }
  if(y %in% RVS[min(d-k+2,d):d,col.x]){
    return(FALSE)
  }
  if((y %in% diag(RVS)[((col.x+1):d)])==FALSE){
    return(FALSE) # check whether y is in a diagonal element right from the column it shall be placed in
  }
  else{
    if(y %in% diag(RVS)[((col.x+1):d)]){
      col.y <- which(diag(RVS) %in% y)
      if(length(setdiff(as.vector(RVS[((d-k+2):d),col.x]),as.vector(RVS[((d-k+2):d),col.y])))==0){
        codeT1 <- TRUE
      }else{
        codeT1 <- FALSE
      }
    }
    m <- col.x+1
    while(codeT2 == FALSE && m <= (d-1) && codeT1 == FALSE){
      ### Works as designed since both sets have same majority and each column in a R-Vine Matrix contains each of its values only once
      checkSet.left <- union(y,as.vector(RVS[min(d,d-k+2):d,col.x]))
      checkSet.right <- union(as.vector(RVS[min(d,d-k+2):d,m]),RVS[m,m])
      if(length(setdiff(checkSet.left,checkSet.right))==0){
        codeT2 <- TRUE
      }
      else{
        codeT2 <- FALSE
      }
      m <- m+1
    }
  }
  return(codeT1 || codeT2)
}
