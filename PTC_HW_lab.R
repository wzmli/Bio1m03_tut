library(dplyr)

two_way_chisq <- function(pt, npt, pnt, npnt){
  p_num <- pt + pnt
  np_num <- npt + npnt
  t_num <- pt + npt
  nt_num <- pnt + npnt
  expected_pt <- p_num*(t_num/(t_num+nt_num))
  expected_npt <- np_num*(t_num/(t_num+nt_num))
  expected_pnt <- p_num*(nt_num/(t_num+nt_num))
  expected_npnt <- np_num*(nt_num/(t_num+nt_num))
  
  chisq <- sum( 
    ((pt - expected_pt)^2)/expected_pt,
    ((npt - expected_npt)^2)/expected_npt,
    ((pnt - expected_pnt)^2)/expected_pnt,
    ((npnt - expected_npnt)^2)/expected_npnt
  )
  return(chisq)
} 

genotype_df <- function(AA,Aa,aa){
  AAdat <- data.frame(A1=rep(1,AA),A2=rep(1,AA))
  Aadat <- data.frame(A1=rep(1,Aa),A2=rep(0,Aa))
  aadat <- data.frame(A1=rep(0,aa),A2=rep(0,aa))
  dd <- rbind(AAdat,Aadat,aadat)
  return(dd)
}

HWE_chisq <- function(x){
  x[,3] <- x[,1] + x[,2]
  AA <- sum(x[,3]==0)
  Aa <- sum(x[,3]==1)
  aa <- sum(x[,3]==2)
  n <- sum(AA,aa,Aa)


  p <- (AA*2 + Aa)/(2*n)
  q <- (aa*2 + Aa)/(2*n)

  eAA <- n*p^2
  eAa <- n*2*p*q
  eaa <- n*q^2

  g_chisq <- sum( ((AA-eAA)^2)/eAA 
                , ((Aa - eAa)^2)/eAa
                , ((aa - eaa)^2)/eaa
                )
  return(g_chisq)
}
  
new_generation <- function(d){
  mixrow <- function(x){
    dd_p <- x[sample(1:nrow(x),nrow(x)),]
    return(dd_p)
  }
  
  d <- mixrow(d)
  d1 <- head(d,(nrow(d)/2))
  d2 <- tail(d,(nrow(d)/2))
  colnames(d2) <- c("A1p","A2p")
  
  new_d <- (cbind(d1,d2)
              %>% rowwise()
              %>% transmute(A1_1 = sample(c(A1,A1p),size=1,prob=c(0.5,0.5))
                         , A2_1 = sample(c(A2,A2p),size=1,prob=c(0.5,0.5))
                         , A1_2 = sample(c(A1,A1p),size=1,prob=c(0.5,0.5))
                         , A2_2 = sample(c(A2,A2p),size=1,prob=c(0.5,0.5))
              )
  )
  new_gen <- data.frame(A1=c(new_d$A1_1,new_d$A1_2)
                        , A2=c(new_d$A2_1,new_d$A2_2))
  return(new_gen)
}

##simple example
set.seed(101)
A1 <- sample(c(0,1),50,replace=TRUE,prob=c(0.8,0.2))
A2 <- sample(c(0,1),50, replace=TRUE,prob=c(0.5,0.5))
dd <- data.frame(A1,A2)
HWE_chisq(dd)
dd2 <- dd %>% new_generation
HWE_chisq(dd2)
dd3 <- dd2 %>% new_generation
HWE_chisq(dd3)
dd4 <- dd3 %>% new_generation
HWE_chisq(dd4)
dd5 <- dd4 %>% new_generation
HWE_chisq(dd5)


dd <- genotype_df(11,43,0)
HWE_chisq(dd)
