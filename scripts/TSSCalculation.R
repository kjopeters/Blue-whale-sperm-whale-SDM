
library(PresenceAbsence)## call libraries

getTSS <- function(model,data,df) {
  
  ids = seq(1,nrow(df),1)
  
trh=seq(0,1,0.001);nt=length(trh);#set up treshold test
mat<-cbind(ids,data,model);strh<-0;sTSS<-0;#save treshold and TSS 
for (i in 1:nt) {
  conf<-cmx(mat, threshold = trh[i], which.model = 1, na.rm = FALSE);#creating confusion matrix
  out.sens<-sensitivity(conf, st.dev = FALSE);#calculation sensitivity
  out.spe<-specificity(conf, st.dev = FALSE);#calculationn specificity
  out.tss<-(out.sens+out.spe)-1;#calculation TSS
  if (out.tss>sTSS){
    sTSS<-out.tss;strh<-trh[i];
  }
  rm(out.spe,out.sens,out.tss,conf)
}

return(data.frame(c(sTSS,strh), c("sTSS", "treshold")))

}
