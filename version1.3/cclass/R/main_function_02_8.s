## dataset: the input dataset which need analysing
## ccol: the features of dataset which will be considered as dimensions
## cluster: maximum clusters will be compared(from 2 to 19)
## iter: iterated times of k-means and fcm algorithm
## dist: the method will be used for fcm and pam, choice: manhattan and euclidean

# cclass <- function(x, ...) UseMethod("cclass")

cclass= function (x, ccol = c(1:dim(x)[2]),cluster=19, iter = 100, dist="euclidean")
{
  
  dataset <- as.matrix(x[,ccol])

  if(length(dataset) > length(na.omit(dataset)))
  {
    dataset <- na.omit(dataset)
    warning("There are some missing data in your choosing range, they have been removed from matrix")
  }
  
  isNum <- sapply(dataset, is.numeric)
  for(i in length(isNum))
    if(isNum[i] == FALSE)
    {
      stop("Sorry, the data you chosen must be numeric.")
    }
  
  
  wide     = dim(dataset)[2]
  length   = dim(dataset)[1]
  
  ### load several libraries (some of them maybe are not needed, but I don't know anymore which one...):
  library(cclust)
  library(e1071)
  library(cluster)
  
  ### plot function has been defined here
  ### which will return the plot it made
  
  bip = function (x, choices = 1:2, scale = 1, pc.biplot = FALSE, col=NULL, color,...) 
  {
    if (length(choices) != 2) 
      stop("length of choices must be 2")
    if (!length(scores <- x$scores)) 
      stop(gettextf("object '%s' has no scores", deparse(substitute(x))), 
           domain = NA)
    lam <- x$sdev[choices]
    if (is.null(n <- x$n.obs)) 
      n <- 1
    lam <- lam * sqrt(n)
    if (scale < 0 || scale > 1) 
      warning("'scale' is outside [0, 1]")
    if (scale != 0) 
      lam <- lam^scale
    else lam <- 1
    if (pc.biplot) 
      lam <- lam/sqrt(n)
    bips(t(t(scores[, choices])/lam), t(t(x$loadings[, choices]) * lam), col = color, ...)
    legend("bottomright", c(color$legendName),
           pch=c(color$legendPch), col=c(color$legendColor), bty="o",cex = 1,y.intersp = 0.8)
    x <- recordPlot()
    return(x)
  }
    
  bips = function (x, y, var.axes = TRUE, col, col2="black", cex = rep(par("cex"), 2), 
                   xlabs = NULL, ylabs = NULL, expand = 1, xlim = NULL, ylim = NULL, 
                   arrow.len = 0.1, main = NULL, sub = NULL, xlab = NULL, ylab = NULL, 
                   ...) 
  {
    n <- nrow(x)
    p <- nrow(y)
    if (missing(xlabs)) {
      xlabs <- dimnames(x)[[1]]
      if (is.null(xlabs)) 
        xlabs <- 1:n
    }
    xlabs <- as.character(xlabs)
    dimnames(x) <- list(xlabs, dimnames(x)[[2]])
    if (missing(ylabs)) {
      ylabs <- dimnames(y)[[1]]
      if (is.null(ylabs)) 
        ylabs <- paste("Var", 1:p)
    }
    if (length(cex) == 1) 
      cex <- c(cex, cex)
    unsigned.range <- function(x) c(-abs(min(x)), abs(max(x)))
    rangx1 <- unsigned.range(x[, 1])
    rangx2 <- unsigned.range(x[, 2])
    rangy1 <- unsigned.range(y[, 1])
    rangy2 <- unsigned.range(y[, 2])
    if (missing(xlim) && missing(ylim)) 
      xlim <- ylim <- rangx1 <- rangx2 <- range(rangx1, rangx2)
    else if (missing(xlim)) 
      xlim <- rangx1
    else if (missing(ylim)) 
      ylim <- rangx2
    ratio <- max(rangy1/rangx1, rangy2/rangx2)/expand
    on.exit(par(op))
    op <- par(pty = "s")
    if (!is.null(main)) 
      op <- c(op, par(mar = par("mar") + c(0, 0, 1, 0)))
    plot(x, type = "n", xlim = xlim, ylim = ylim, col = col$colors,   
         xlab = xlab, ylab = ylab, sub = sub, main = main,...)
    
    for(i in 1 : n)
    {
      if(col$colors[i] == col$legendColor[1])
      {
        points(x[i,1],x[i,2],cex = cex[1], col = col$colors[i],pch = 4,...)
      }
      else 
      {
        points(x[i,1],x[i,2],cex = cex[1], col = col$colors[i],pch = 1,...)
      }
    }
    par(new = TRUE)
    plot(y, axes = FALSE, type = "n", xlim = xlim * ratio, ylim = ylim * 
           ratio, xlab = "", ylab = "", col = col2, ...)
    box(col = col2)
    text(y, labels = ylabs, cex = cex[2], col = col2, ...)   ### general use
    if (var.axes)                                            ### general use 
      arrows(0, 0, y[,1] * 0.8, y[,2] * 0.8, col = "black",   ### general use
             length = arrow.len)  ### general use
  }
  
  ## function to align the labels
  align <- function(final,label)
  {
    
    #function to pick out rows
    pick<-function(data,row,value,choose){
      length = dim(data)[1]
      row_length = table(data[,row] == value)[2]
      col_length = length(choose)
      result = matrix(nrow = row_length, ncol = col_length)
      
      a = 1
      b = 1
      
      for(i in 1:length){
        if (data[b,row] == value ){
          result[a,] = data[b,choose]
          a = a+1
        } 
        b = b+1
      }
      
      head(result)
      return (result) 
    }
    
    #function to compare distant
    dis<-function(a,b){
      sqrt(sum((a-b)^2))
    }
  
    length = dim(final)[1]
    wide   = dim(final)[2]
    length_label = length(label)  
    l = length(table(final[,label[1]])) 
    length_attr = wide - length_label 
    attr = 1:length_attr  
    
    
    pivot = array(0,c(l,length_attr,length_label))
    
    
    #calculate the centers of the first label attribute
    for(w in 1:length_label){ 
      for(u in 1:l){
        pivot[u,,w] <-colMeans(pick(final,label[w],u,attr))
      }
    }
    
    #initialize the rank array
    coff = array(0,c(l,l,length_label-1))
    
    

#  get the rank table based on the distance matrix
for(a in 1:(length_label-1)){ 
  for(b in 1:l){
    for(c in 1:l){ 
      distant = dis(pivot[c,,1],pivot[b,,a+1])
      coff[c,b,a] = distant
    }
  }
}

#   the distance matrix
    per = permutations(n = l)
    sum_row = nrow(per)
    sum = rep(data = 0,sum_row)
    rank = matrix(data = NA, nrow = l, ncol = l)

    for(i in 1:2)
    {
      for(j in 1:sum_row)
      {
        for(a in 1:(length_label-1)) 
        {
          for(x in 1:(length_label-1))
          {
            sum[j] = sum[j] + coff[x,per[a,x],i]
          }
        }
      }
      rank[,i] <- t(per[which.min(sum),])    
    }
  

    
    #according to the rank table, assign the new value to the orginal data 
    for(ccc in 1:length) {  
      for(aaa in 1:(length_label-1)){  
        final[ccc,label[1]+aaa] = rank[final[ccc,label[1]+aaa],aaa]
      }  
    }
    
    return(final)
    
  }

  ## function to select best number of clusters
  cluster_num =function (dataset,cluster)
  {
    ### apply K-means from 2 to 20 clusters and build matrix of validity indeces   
    cluster_table1 <- data.frame(index=c('Cal','Scott','Marr','TracW','Fried'))
    cluster_table2 <- data.frame(index=c('Fuz','PartDen','Xie','PartCoe'))
    isCluster <- TRUE
    ## Index function for pam and k-means
    plot_indices <- function(indexes,cluster)
    {
      diff_cal <- rep(NA,cluster)
      diff_scott <- rep(NA,cluster)
      diff_marr <- rep(NA,cluster)
      diff_tracw <- rep(NA,cluster)
      diff_fried <- rep(NA,cluster)
      output <- list()
      
      for (j in 2:cluster) {
        diff_cal[j] <- (indexes[j+1,1]-indexes[j,1])-(indexes[j,1]-indexes[j-1,1])}
      for (j in 2:cluster) {
        diff_scott[j] <- (indexes[j,6]-indexes[j-1,6])}
      for (j in 2:cluster) {
        diff_marr[j] <- (indexes[j+1,7]-indexes[j,7])-(indexes[j,7]-indexes[j-1,7])}
      for (j in 2:cluster) {
        diff_tracw[j] <- (indexes[j+1,10]-indexes[j,10])-(indexes[j,10]-indexes[j-1,10])}
      for (j in 2:cluster) {
        diff_fried[j] <- (indexes[j,11]-indexes[j-1,11])}
      output$bestCluster=c(which.min(diff_cal),which.max(diff_scott),
                           which.max(diff_marr),which.max(diff_tracw),which.max(diff_fried))
      output$clusterRank <- matrix(data = 0, nrow = 5, ncol = cluster - 1)
      for(i in 2:cluster)
      {
        output$clusterRank[1,i-1] <- rank(diff_cal)[i]
        output$clusterRank[2,i-1] <- cluster-rank(diff_scott)[i]
        output$clusterRank[3,i-1] <- cluster-rank(diff_marr)[i]
        output$clusterRank[4,i-1] <- cluster-rank(diff_tracw)[i]
        output$clusterRank[5,i-1] <- cluster-rank(diff_fried)[i]
      }
      return(output)
    }
    
    ## Patch the indexes for one cluster
    addClusterOne <- function(indexes,inputDataset)
    {
      
      ## Calculate scott index, marriot, tracew and friedman index when cluster is 1.
      rows <- nrow(inputDataset)
      oneCluster <- rep(as.integer(1),time=rows)
      
      
      ttww <- function(x, clsize, cluster) {
        n <- sum(clsize)
        k <- length(clsize)
        w <- 0
        tt <- cov(x) * n
        for (l in 1:k) w <- w + cov(x[cluster == l, ]) * clsize[l]
        zttw <- list(tt = tt, w = w)
        return(zttw)
      }
      
      scott <- function(zttw, clsize) {
        n <- sum(clsize)
        dettt <- prod(eigen(zttw$tt)$values)
        detw <- prod(eigen(zttw$w)$values)
        scott <- n * log(dettt/detw)
        return(scott)
      }
      marriot <- function(zttw, clsize) {
        k <- length(clsize)
        detw <- prod(eigen(zttw$w)$values)
        mar <- (k^2) * detw
        return(mar)
      }
      tracew <- function(zttw) {
        tracew <- sum(diag(zttw$w))
        return(tracew)
      }
      friedman <- function(zttw) {
        b <- zttw$tt - zttw$w
        fried <- sum(diag(solve(zttw$w) %*% b))
        return(fried)
      }
      zttw<-ttww(inputDataset,rows,oneCluster)
      zscott<-scott(zttw,rows)
      zmarriot<-marriot(zttw,rows)
      ztracew<-tracew(zttw)
      zfriedman<-friedman(zttw)
      
      ## Calinski
      indexes[1,1] <- 0
      ## scott
      indexes[1,6] <- zscott
      ## Marriot
      indexes[1,7] <- zmarriot
      ## TraceW
      indexes[1,10] <- ztracew
      ## Fried
      indexes[1,11] <- zfriedman
      
      return (indexes)
    }
    
## function to plot the fuzzy val indices
    plot_indices_fcm <- function(indexes_fcm,cluster)
    {
      diff_fhv <- indexes_fcm[,1]
      diff_pd <- indexes_fcm[,3]
      diff_xb <- indexes_fcm[,4]
      diff_pc <- indexes_fcm[,6]
      
      output <- list() 
      output$bestCluster=c(which.min(diff_fhv),which.max(diff_pd),which.min(diff_xb),
                           which.max(diff_pc))
      output$clusterRank <- matrix(data = 0, nrow = 4, ncol = cluster - 1)
      for(i in 2:cluster)
      {
        output$clusterRank[1,i-1] <- rank(diff_fhv)[i]
        output$clusterRank[2,i-1] <- cluster-rank(diff_pd)[i]
        output$clusterRank[3,i-1] <- rank(diff_xb)[i]
        output$clusterRank[4,i-1] <- cluster-rank(diff_pc)[i]
      }
      return(output)
    }
    
## Calculate best Cluster number
    bestCluster <- function(kmRankTable, pamRankTable, fcmRankTable, cluster)
    {
      finalRankTable <- rep(NA, cluster - 1)
      a <- 0
      b <- 0
      c <- 0
      for(i in 1:(cluster - 1))
      {
        for(l in 1:5)
        {
          c <- c + kmRankTable[l,i] * (3/8)
        }
        for(j in 1:5)
        {
          a <- a + pamRankTable[j,i] * (3/8)
        }
        for(k in 1:4)
        {
          b <- b + fcmRankTable[k,i] * (1/4)
        }
        finalRankTable[i] <- a + b + c
        a <- 0
        b <- 0
        c <- 0
      }
      return(which.min(finalRankTable)+1)
    }
    
    while(isCluster){
    isCluster <- FALSE
    kmIndexes <- matrix(data = 0, nrow = cluster + 1, ncol = 15)
    pamIndexes <- matrix(data = 0, nrow = cluster + 1, ncol = 15)
    fcmIndexes <- matrix(data = NA, nrow = cluster, ncol = 9)
    kmRankTable <- matrix(data = 0, nrow = 5, ncol = (cluster - 1))
    pamRankTable <- matrix(data = 0, nrow = 5, ncol = (cluster -1))
    fcmRankTable <- matrix(data = 0, nrow = 4, ncol = (cluster - 1))
    
    ### HCA algorithm
    h <- hclust(dist(dataset), method="average")
    
    ### K-means part
    for (i in 2:(cluster+1)) {
      initial <- tapply(dataset, list(rep(cutree(h,i), ncol(dataset)), col(dataset)), mean)
      cl <- kmeans(dataset,initial,iter)
      if(suppressWarnings(min(cl$size) == 1))
      {
        cluster <- i - 2
        isCluster <- TRUE
        warning(paste("The cluster of current dataset cannot be set more than",cluster))
        break
      }
      
      if("try-error" %in% class(try(suppressWarnings(clustIndex(cl,dataset,index="all")))))
      {
        cluster <- i - 2
        if(cluster < 2)
        {
          stop("The indexes cannot be calculated")
        }
        isCluster <- TRUE
        warning(paste("The cluster of current dataset cannot be set more than",cluster))
        break
      }
      kmIndexes[i,] <- suppressWarnings(clustIndex(cl,dataset,index="all"))
    }
    
    ## If cluster is reset, start another circle of this loop    
    if(isCluster)
    {
      next
    }
 
    kmIndexes <- addClusterOne(kmIndexes,dataset)   
    KM<-list()
    KM<-plot_indices(kmIndexes,cluster)
    cluster_table1<-cbind(cluster_table1,KM$bestCluster)
    kmRankTable <- KM$clusterRank
    
   miowss <- function(centri, cluster, dati) {
      retval <- rep(0, nrow(centri))
      x <- (dati - centri[cluster, ])^2
      for (k in 1:nrow(centri)) {
        retval[k] <- sum(x[cluster == k, ])
      }
      retval
    }
    
    ## PAM algorithm
    for (i in 2:(cluster + 1)) {
      cl <- pam(dataset, i, diss = inherits(dataset, "dist"), metric = dist,
                                  medoids = NULL, stand = FALSE, cluster.only = FALSE)
      if(suppressWarnings(min(cl$size) == 1))
      {
        cluster <- i - 2
        isCluster <- TRUE
        warning(paste("The cluster of current dataset cannot be set more than",cluster))
        break
      }
      cl$centers <- cl$medoids
      cl$cluster <- cl$clustering
      cl$size <- table(cl$cluster)[]
      cl$withins <- miowss(cl$centers, cl$cluster, dataset)
      
      if("try-error" %in% class(try(suppressWarnings(clustIndex(cl, dataset, index="all")))))
      {
        cluster <- i - 2
        if(cluster < 2)
        {
        stop("The indexes cannot be calculated")
        }
        isCluster <- TRUE
        warning(paste("The cluster of current dataset cannot be set more than",cluster))
        break
      }
      pamIndexes[i,] <- suppressWarnings(clustIndex(cl, dataset, index="all"))
    }
   
   ## If cluster is reset, start another circle of this loop 
   if(isCluster)
   {
      next
   } 
    pamIndexes <- addClusterOne(pamIndexes,dataset)
    
    PAM <- plot_indices(pamIndexes,cluster)
    cluster_table1<-cbind(cluster_table1,PAM$bestCluster)
    pamRankTable <- PAM$clusterRank
    
    for (i in 2:cluster) {
      initial <- tapply(dataset, list(rep(cutree(h,i), ncol(dataset)), col(dataset)), mean)
      cl <- cmeans(dataset,initial, iter,verbose=FALSE,dist = dist, method="cmeans",m=2,initial)
      if(suppressWarnings(min(cl$size) == 1 ))
      {
        cluster <- i - 2
        isCluster <- TRUE
        break
      }
      
      
      if("try-error" %in% class(try(suppressWarnings(fclustIndex(cl,dataset,index="all")))))
      {
        cluster <- i - 2
        if(cluster < 2)
        {
          stop("The indexes cannot be calculated")
        }
        isCluster <- TRUE
        warning(paste("The cluster of current dataset cannot be set more than",cluster))
        break
      }
      fcmIndexes[i,] <- suppressWarnings(fclustIndex(cl,dataset,index="all"))
    }
    
   ## If cluster is reset, start another circle of this loop 
   if(isCluster)
   {
     next
   }
    FCM <- list()
    FCM<-plot_indices_fcm(fcmIndexes,cluster)
    cluster_table2<-cbind(cluster_table2,FCM$bestCluster)
    fcmRankTable <- FCM$clusterRank
   }
    clusterNum <- 0
    clusterNum <- bestCluster(kmRankTable, pamRankTable, fcmRankTable,cluster)    
    return (clusterNum)
  }


  ## function to generate the representing color according to the number of cluster
  cluster_color <- function(cluster, final_dataset)
  {
    color_cluster <- list()
    dataset_number <- dim(final_dataset)[1]
    colors_use <- rainbow(cluster + 1)
    cluster_label = rep(NA, cluster + 1)
    cluster_label[1] <- "NC"
    colors_output <- rep(NA, dataset_number)
    pch_sty <- rep(4, dataset_number)
    for(i in 1 : cluster){
      cluster_label[i+1] <- i 
    }
    for(j in 1 : dataset_number)
    {
      for(k in 1 : (cluster + 1))
      {
        if(final_dataset$class[j] == cluster_label[k])
        {
          colors_output[j] <- colors_use[k]
          break
        }
      }
    }
    for(l in 1 : (cluster + 1))
    {
      cluster_label[l] <- paste("Cluster",cluster_label[l],sep=" ")
    }
    for(m in 2 : (cluster + 1))
    {
      pch_sty[m] <- 1
    }
    color_cluster$colors <- colors_output
    color_cluster$legendName <- cluster_label
    color_cluster$legendColor <- colors_use
    color_cluster$legendPch <- pch_sty
    return(color_cluster)
  }
  
  
  
  ## find out the number of labels (clusters)
  number = cluster_num(dataset, cluster)
  
  #create the tree of the orginal tree
  h <- hclust(dist(dataset), method="average")
  
  ##Kmeans
  initial <- tapply(dataset, list(rep(cutree(h,number), ncol(dataset)), col(dataset)), mean)
  cl_km <- kmeans(dataset,initial, iter)
  km <- data.frame(dataset, clust=cl_km$cluster)    
  
  ##pam
  cl_pam <- pam(dataset, number, diss = inherits(dataset, "dist"), metric = dist,
                 medoids = NULL, stand = FALSE, cluster.only = FALSE)
  pam <- data.frame(dataset, clust=cl_pam$clustering)  
  
  ##FCM
  initial <- tapply(dataset, list(rep(cutree(h,number), ncol(dataset)), col(dataset)), mean)
  cl <- cmeans(dataset,initial,iter,verbose=FALSE,dist= dist, method="cmeans",m=2,initial)
  fcm <- data.frame(dataset, clust=cl$cluster)
  
  #construct all results together
  classif <- data.frame(cbind(km$clust, pam$clust, fcm$clust))
  names(classif) <- c("km", "pam", "fcm")
  result <- data.frame(dataset,classif)
  result = as.matrix(result)
  
  #align label, according the number of clutering methods, the paramater can be changed
  aligned_result = align(result,(wide+1):(wide+3))
  
  interss = c()
  for(i in 1:length){
    if(aligned_result[i,wide+1] == aligned_result[i,wide+2] & aligned_result[i,wide+2] == aligned_result[i,wide+3] ){
      interss = c(interss,as.integer(aligned_result[i,wide+1]))
    }else{
      interss = c(interss,"NC")
    }
  }

  common <- list()  
  
common$cluster <- data.frame(dataset, km = cl_km$cluster, pam = cl_pam$clustering, fcm = cl$cluster, class=interss) 
  
  colors <- cluster_color(number, common$cluster)
  common$plot <- bip(princomp(common$cluster[,1:wide], cor=T),color = colors)
  
  common$number = number
  
  class(common) <- "cclass"
  common
}

print.cclass <- function(x)
{
  wide     = dim(x$cluster)[2]
  table    =  table((x$cluster)[,wide])
  length_label = length(table)
  
  cat("The best number of clusters is ", x$number,"\n\n")
  
  for(i in 1:(length_label-1) ){
    cat("The number of instances assigned to label", i,"is",table[i],"\n")
  }
    
  cat("The number of instances assigned to NC is",table[length_label],"\n")
}