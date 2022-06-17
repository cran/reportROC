\name{reportROC}
\alias{reportROC}
\title{An Easy Way to Report ROC Analysis}
\description{
Provides an easy way to report the results of ROC analysis, including:
  1. an ROC curve. 2. the value of Cutoff,
  AUC (Area Under Curve), ACC (accuracy),
  SEN (sensitivity), SPE (specificity),
  PLR (positive likelihood ratio), NLR (negative likelihood ratio),
  PPV (positive predictive value), NPV (negative predictive value).
}
\keyword{ROC analysis}
\usage{
reportROC(gold,predictor,predictor.binary,important,positive,exact,plot,xlab,ylab)
}
\arguments{
  \item{gold}{numeric(0/1) or binary, the 'gold standard'; typically encoded with 0 (controls) and 1 (cases)}
  \item{predictor}{numeric, the predictor variable}
  \item{predictor.binary}{numeric(0/1) or binary, if this argument was used, other arguments including 'predictor' and 'important' would be disabled}
  \item{important}{'se' or 'sp', some Youden index maybe the same, and 'important' is to indicate which is more important between sensitivity and specificity}
  \item{positive}{logic, 'l': the larger predictor or predictor.binary indicates the 'cases', 's': the smaller one indicates the 'cases'}
  \item{exact}{logic, whether to calculate the exact p value for AUC}
  \item{plot}{logic, whether to plot the ROC curve with specific style}
  \item{xlab}{character, the name of X axis}
  \item{ylab}{character, the name of Y axis}
}
\value{
  \item{Curoff}{cutoff, only for numeric predictor}
  \item{AUC}{Area Under Curve, AUC}
  \item{AUC.SE}{the standard error of AUC}
  \item{AUC.low,AUC.up}{the 95 percent CI of AUC}
  \item{P}{the p value for AUC using one-sided test, which is different from SPSS using two-sided test}
  \item{ACC}{accuracy}
  \item{ACC.low,ACC.up}{the 95 percent CT of accuracy}
  \item{SEN,SEN.low,SEN.up}{the esitmate and the 95 percent CI of sensitivity}
  \item{SPE,SPE.low,SPE.up}{the esitmate and the 95 percent CI of specificity}
  \item{PLR,PLR.low,PLR.up}{the esitmate and the 95 percent CI of positive likelihood ratio}
  \item{NLR,NLR.low,NLR.up}{the esitmate and the 95 percent CI of negative likelihood ratio}
  \item{PPV,PPV.low,PPV.up}{the esitmate and the 95 percent CI of positive predictive value}
  \item{NPV,NPV.low,NPV.up}{the esitmate and the 95 percent CI of negative predictive value}
  \item{PPA,PPA.low,PPA.up}{the esitmate and the 95 percent CI of percentage of positive accordance}
  \item{NPA,NPA.low,NPA.up}{the esitmate and the 95 percent CI of percentage of negative accordance}
  \item{TPA,TPA.low,TPA.up}{the esitmate and the 95 percent CI of percentage of total accordance}
  \item{KAPPA,KAPPA.low,KAPPA.up}{the esitmate and the 95 percent CI of Kappa}
}
\author{
Zhicheng Du<dgdzc@hotmail.com>, Yuantao Hao<haoyt@mail.sysu.edu.cn>
}
\note{
Please feel free to contact us, if you have any advice and find any bug!

Update description:

version 2.0: 1. 95 percent CIs for AUC, SEN, SPE, PLR, NLR, PPV and NPV are available.

version 3.0: 1. binary predictor is available with the new argument 'predictor.binary'. 2. positive indicator is available with the new argument 'positive'.

version 3.1: 1. accuracy is available. 2. roc curve is available for binary predictor.

version 3.2: 1. data with missing values can be handled.

version 3.3: 1. fixed the bug of the same value in 'AUC' and 'AUC.low'.

version 3.4: 1. fixed the axis names of the ROC plot. Thank Cesar S. Rabak (csrabak@hotmial.com) for the useful feedback.

version 3.5: 1. add p value to the outputs.

version 3.6: 1. add PPA, NPA, TPA, and KAPPA to the outputs. 2. limit sensitivity and specificity to a maximum of one. Thank Tengfei Song (songtf@mail2.sysu.edu.cn) for the useful feedback. 3. add the argument 'exact' to avoid the wanring message of "cannot compute exact p-value" in wilcox.test. 4. add the 95 percent CI to AUC in the ROC cruve.

more functions will be included in 'reportROC' package!
}
\seealso{
nothing
}
\examples{
data(aSAH)
#for continuous variables
reportROC(gold=aSAH$outcome,predictor=aSAH$s100b,important="se",plot=TRUE)
#for binary variables
binary=rep(0,nrow(aSAH))
binary[aSAH$s100b>=0.205]=1
reportROC(gold=aSAH$outcome,predictor.binary=binary)
reportROC(gold=aSAH$outcome[1:50],predictor.binary=binary[1:50],exact=FALSE)
}

