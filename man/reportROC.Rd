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
reportROC(gold,predictor,predictor.binary,important,plot,xlab,ylab,positive)
}
\arguments{
  \item{gold}{numeric(0/1) or binary, the 'gold standard'; typically encoded with 0 (controls) and 1 (cases)}
  \item{predictor}{numeric, the predictor variable}
  \item{predictor.binary}{numeric(0/1) or binary, if this argument was used, other arguments including 'predictor' and 'important' would be disabled}
  \item{important}{'se' or 'sp', some Youden index maybe the same, and 'important' is to indicate which is more important between sensitivity and specificity}
  \item{plot}{logic, whether to plot the ROC curve with specific style}
  \item{xlab}{character, the name of X axis}
  \item{ylab}{character, the name of Y axis}
  \item{positive}{logic, 'l': the larger predictor or predictor.binary indicates the 'cases', 's': the smaller one indicates the 'cases'}
}
\value{
  \item{Curoff}{cutoff, only for numeric predictor}
  \item{AUC}{Area Under Curve, AUC}
  \item{AUC.SE}{the standard error of AUC}
  \item{AUC.low,AUC.up}{the 95 percent CI of AUC}
  \item{P}{the p value for AUC using one-sided test, which is different from SPSS using two-sided test}
  \item{ACC}{accuracy}
  \item{ACC.low,ACC.up}{the 95 percent CT of accuracy}
  \item{SEN}{sensitivity}
  \item{SEN.low,SEN.up}{the 95 percent CI of sensitivity}
  \item{SPE}{specificity}
  \item{SPE.low,SPE.up}{the 95 percent CI of specificity}
  \item{PLR}{positive likelihood ratio}
  \item{PLR.low,PLR.up}{the 95 percent CI of PLR}
  \item{NLR}{negative likelihood ratio}
  \item{NLR.low,NLR.up}{the 95 percent CI of NLR}
  \item{PPV}{positive predictive value}
  \item{PPV.low,PPV.up}{the 95 percent CI of PPV}
  \item{NPV}{negative predictive value}
  \item{NPV.low,NPV.up}{the 95 percent CI of NPV}
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

version 3.4: 1. fixed the axis names of the ROC plot. Thank Cesar S. Rabak (csrabak@hotmial.com) for the useful comments.

version 3.5: 1. add p value to the outputs.

more functions will be included in 'reportROC' package!
}
\seealso{
nothing
}
\examples{
data(aSAH)
reportROC(gold=aSAH$outcome,predictor=aSAH$s100b,important="se",plot=TRUE)

binary=rep(0,nrow(aSAH))
binary[aSAH$s100b>=0.205]=1
reportROC(gold=aSAH$outcome,predictor.binary=binary)
}

