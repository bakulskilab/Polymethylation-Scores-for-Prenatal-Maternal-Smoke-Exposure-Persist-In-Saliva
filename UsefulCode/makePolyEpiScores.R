## ---- makePolyEpiScores.R
#' Title
#'
#' @param m a matrix; cleaned matrix of methylation values (0-1) from your data where rows are cpg sites and columns are samples; rownames should be of format cg########
#' @param b a list; a list of dataframes containing a column CpG with cpg sites named in format cg######## and a column Coef containing relavent coeficient, should be a named list (names will become score names)
#' @param transformopt character, the method for transforming the m values, default none, options zscore or meancenter
#'
#' @return scores
#' @export
#'
#' @examples
makePolyEpiScore<-function(m, b, transformopt='none'){
  #check inputs
  try(if(is.matrix(m)==F) stop('m is not a matrix'))
  try(if(is.list(b)==F) stop('b is not a list'))
  if(is.null(names(b))){names(b)=paste0('option', 1:length(b))}
  #generate poly methylation scores
  scores<-lapply(b, function(x){
    commonp<-Reduce(intersect,list(x%>%filter(!is.na(Coef))%>%pull(`CpG`), rownames(m)))%>% unique # filter to non-missing shared cpg sites
    beta1<-m[match(commonp, rownames(m)),] #pull shared cpg sites from methylation matrix 
    #apply transform to cpgs of the cleaned methylation matrix 
    print(paste0('transforming the methylation matrix:', transformopt))
    if(transformopt=='zscore'){beta1=t(scale(t(beta1)))} #center and scale the rows (cpgs) of the (shared cpg sites) methylation matrix
    if(transformopt=='meancenter'){beta1=t(scale(t(beta1), center=T, scale=F))} #center the rows (cpgs) of (shared cpgs sites) methylation matrix only
    x1<-x %>% filter(`CpG` %in% commonp) #pull shared cpg sites from coefficient matrix
    (as.numeric(x1[['Coef']]) %*% beta1)%>%unlist%>%c #matrix multiplication of shared coefficients and shared methylation
  })%>% bind_cols()
  return(scores)
}
