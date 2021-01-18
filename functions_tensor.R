# Matricization (Hoff (2015))
mat<-function(A,j)
{
  Aj<-t(apply(A,j,"c"))
  if(nrow(Aj)!=dim(A)[j])  { Aj<-t(Aj) }
  Aj
}

# Kronecker product
kron<-function(...)
{
  M<-list(...)
  if(is.list(M[[1]])) { M<-M[[1]] }
  JM<-1 ; for(j in 1:length(M)){ JM<-kronecker(JM,M[[j]]) }
  JM
}