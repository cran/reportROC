reportROC=function(gold,predictor,important="se",plot=TRUE){

  error=FALSE
  if(length(unique(gold))!=2){
    message("Error! The 'gold' variable must be binomial!");error=TRUE
  }
  if(class(predictor)!="numeric"){
    message("Error! The 'predictor' variable must be numeric!");error=TRUE
  }

  if(!error){

    roc=roc(gold~predictor,auc=TRUE,ci=TRUE)

    sensitivities=roc$sensitivities
    specificities=roc$specificities
    cutpoints=roc$thresholds
    yuedeng=sensitivities+specificities-1
    cut.s=cutpoints[yuedeng==max(yuedeng)]
    se.s=sensitivities[yuedeng==max(yuedeng)]
    sp.s=specificities[yuedeng==max(yuedeng)]
    if(important=="se"){
      cut=cut.s[se.s==max(se.s)]
      se=se.s[se.s==max(se.s)]
      sp=sp.s[se.s==max(se.s)]
    }
    if(important=="sp"){
      cut=cut.s[sp.s==max(sp.s)]
      se=se.s[sp.s==max(sp.s)]
      sp=sp.s[sp.s==max(sp.s)]
    }
    AUC=roc$auc
    AUC.SE=(roc$ci[3]-roc$ci[2])/1.96
    PLR=se/(1-sp)
    NLR=(1-se)/sp
    pre=rep(0,length(roc$predictor))
    pre[roc$predictor>=cut]=1
    pre=table(pre,roc$response)
    pre=as.vector(pre)
    PPV=pre[1]/(pre[1]+pre[3])
    NPV=pre[4]/(pre[4]+pre[2])
    rst=round(data.frame(Cutoff=cut,SEN=se,SPE=sp,AUC,AUC.SE,PLR,NLR,PPV,NPV),3)

    if(plot){
      par(mai=c(1,1,0.3,0.3))
      plot(1-specificities,sensitivities,col="black",xlab="1-Specificity",
           ylab="Sensetivity",type='l',lwd=1,lty=1,cex.lab=1,cex.axis=1)
      points(1-sp,se,col="grey",cex=1,pch=16)
      text(1-sp,se,pos=4,paste("(",round(1-sp,2),", ",round(se,2),")",sep=""))
      legend("bottomright",paste("AUC =",round(AUC,2)),col=1,bty="n",cex=1)
    }
    return(rst)
  }
}
