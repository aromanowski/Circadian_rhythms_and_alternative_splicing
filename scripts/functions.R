# Function to get annotation data
getAnnnot <- function(df, md){
  #' \dontrun{
  #' getAnnot(df = dataframe where first column is named "Gene ID" and contains genes,
  #'          md = metadata where first column is named "Gene ID" and contains genes and has extra columns with annotations  
  #'          to add to the dataframe)
  #' }
  if (colnames(md)[1] != "Gene_ID") {
    colnames(md)[1] <- "Gene_ID"
  }
  if (colnames(df)[1] != "Gene_ID") {
    df <- data.frame(cbind(rownames(df)),df)
    colnames(df)[1] <- "Gene_ID"
  }
  res <- merge(x = md, y = df, by = "Gene_ID", all.y = TRUE) # Left inner join of df and getAnnnot <- function(df, md){
  res[is.na(res)] <- "-"
  return(res)
}

# Function to plot a blue-yellow heatmap 
# (modified from the tightClust package by George C. Tseng and Wing H. Wong)
plot.clust <- function (x, standardize.gene = TRUE, gene_labels = TRUE, order.sample = FALSE, 
                        plot.noise = TRUE, ...) 
{
  mirror.matrix <- function(x) {
    xx <- as.data.frame(x)
    xx <- rev(xx)
    xx <- as.matrix(xx)
    xx
  }
  rotate180.matrix <- function(x) {
    xx <- rev(x)
    dim(xx) <- dim(x)
    xx
  }
  flip.matrix <- function(x) {
    mirror.matrix(rotate180.matrix(x))
  }
  data <- x$data
  cluster <- x$cluster
  size <- x$size
  if (standardize.gene) 
    data <- t(scale(t(data)))
  l <- 1:max(cluster)
  if (plot.noise) 
    l <- c(l, -1)
  order.genes <- unlist(lapply(l, function(i) which(cluster == 
                                                      i)))
  data.plot <- data[order.genes, ]
  if (order.sample) 
    data.plot <- data.plot[, hclust(dist(t(data.plot)))$order]
  cuts <- cumsum(size)
  if (!plot.noise) 
    cuts <- cuts[-length(size)]
  nr <- dim(data.plot)[1]
  nc <- dim(data.plot)[2]
  cexCol <- 1/log10(nc)
  cexRow <- 1/log(nr, 4)
  image(x = 1:nc, y = 1:nr, z = t(flip.matrix(data.plot)), 
        col = colorRampPalette(c("blue", "gray","yellow"))(64), 
        axes = FALSE, xlab = "", ylab = "", ...)
  
  if (gene_labels == TRUE){
    axis(4, nr:1, labels = (row.names(data.plot)), las = 2, 
         line = -0.5, tick = 0, cex.axis = cexRow)}
  
  
  for (i in cuts) {
    lines(c(-1, nc + 1), rep(nr - i + 0.5, 2), col = "white")
  }
  print("Limits: ")
  print(min(data.plot))
  print(max(data.plot))
  
}

# Function to filter genes by read density
filterByRd<-function(df0, targets, min = 0.05, type = "any") 
{
  a<-table(targets$condition) 
  start=ncol(df0)-sum(a)+1
  genes.rd=df0[,start:ncol(df0)]/df0$effective_length
  write.table(genes.rd, file="genes.rd.txt")
  colnames(genes.rd)<-targets$condition
  list<-matrix(unlist(lapply(unique(colnames(genes.rd)), 
                             function(x) rowMeans(genes.rd[,colnames(genes.rd) == x]) > min)
  ),  
  nrow=nrow(genes.rd), 
  byrow=FALSE)
  if (type=="any") { 
    ii=rowSums(list)>0  
    df=df0[ii,]
  }else { 
    ii=rowSums(list)==ncol(list)
    df=df0[ii,]
  }
  return(df)
}

# Function to filter bins by read density
filterByRdBin<-function(dfGen, dfBin, targets, min = 0.05) 
{
  a <- table(targets$condition) 
  
  startG = ncol(dfGen)-sum(a)+1
  startB = ncol(dfBin)-sum(a)+1
  
  genes.rd = dfGen[,startG:ncol(dfGen)]/dfGen$effective_length
  
  bins.rd = dfBin[,startB:ncol(dfBin)]/dfBin$length
  colnames(bins.rd) <- targets$condition
  avRdBin <- matrix(unlist(lapply(unique(colnames(bins.rd)), 
                                  function(x) rowMeans(bins.rd[,colnames(bins.rd) == x]))),  
                    nrow=nrow(bins.rd), 
                    byrow=FALSE)
  
  colnames(genes.rd) <- targets$condition
  
  avRdGen<-matrix(unlist(lapply( unique(colnames(genes.rd)), 
                                 function(x) rowMeans(genes.rd[,colnames(genes.rd) == x]))),  
                  nrow=nrow(genes.rd), 
                  byrow=FALSE)
  
  te <- match(dfBin$locus, rownames(dfGen))
  gen.rdb = avRdGen[te,]
  
  bin.gen.rd = avRdBin/gen.rdb
  ii = rowSums(bin.gen.rd > min)>0
  
  dfBin = dfBin[ii,]
  return(dfBin)
}

# Function to create a standalone scale for heatmaps
image.scale <- function(z, zlim, col = heat.colors(12),
                        breaks, horiz=TRUE, ylim=NULL, xlim=NULL, ...){
  if(!missing(breaks)){
    if(length(breaks) != (length(col)+1)){stop("must have one more break than colour")}
  }
  if(missing(breaks) & !missing(zlim)){
    breaks <- seq(zlim[1], zlim[2], length.out=(length(col)+1)) 
  }
  if(missing(breaks) & missing(zlim)){
    zlim <- range(z, na.rm=TRUE)
    zlim[2] <- zlim[2]+c(zlim[2]-zlim[1])*(1E-3)#adds a bit to the range in both directions
    zlim[1] <- zlim[1]-c(zlim[2]-zlim[1])*(1E-3)
    breaks <- seq(zlim[1], zlim[2], length.out=(length(col)+1))
  }
  poly <- vector(mode="list", length(col))
  for(i in seq(poly)){
    poly[[i]] <- c(breaks[i], breaks[i+1], breaks[i+1], breaks[i])
  }
  xaxt <- ifelse(horiz, "s", "n")
  yaxt <- ifelse(horiz, "n", "s")
  if(horiz){YLIM<-c(0,1); XLIM<-range(breaks)}
  if(!horiz){YLIM<-range(breaks); XLIM<-c(0,1)}
  if(missing(xlim)) xlim=XLIM
  if(missing(ylim)) ylim=YLIM
  plot(1,1,t="n",ylim=ylim, xlim=xlim, xaxt=xaxt, yaxt=yaxt, xaxs="i", yaxs="i", ...)  
  for(i in seq(poly)){
    if(horiz){
      polygon(poly[[i]], c(0,0,1,1), col=col[i], border=NA)
    }
    if(!horiz){
      polygon(c(0,0,1,1), poly[[i]], col=col[i], border=NA)
    }
  }
}

# Function to check the operating system
get_os <- function(){
  sysinf <- Sys.info()
  if (!is.null(sysinf)){
    os <- sysinf['sysname']
    if (os == 'Darwin')
      os <- "osx"
  } else { ## mystery machine
    os <- .Platform$OS.type
    if (grepl("^darwin", R.version$os))
      os <- "osx"
    if (grepl("linux-gnu", R.version$os))
      os <- "linux"
  }
  tolower(os)
}

getEventInfo <- function(df, md){
  #' \dontrun{
  #' getAnnot(df = dataframe where first column is named "Event_ID" and contains genes,
  #'          md = metadata where first column is named "Event_ID" and contains genes and has extra columns with annotations  
  #'          to add to the dataframe)
  #' }
  if (colnames(md)[1] != "Event_ID") {
    md <- data.frame(cbind(rownames(md)),md)
    colnames(md)[1] <- "Event_ID"
  }
  if (colnames(df)[1] != "Event_ID") {
    df <- data.frame(cbind(rownames(df)),df)
    colnames(df)[1] <- "Event_ID"
  }
  res <- merge(x = md, y = df, by = "Event_ID", all.y = TRUE) # Left inner join of df and getAnnnot <- function(df, md){
  res[is.na(res)] <- "-"
  return(res)
}
