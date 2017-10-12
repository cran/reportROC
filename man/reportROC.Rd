\name{reportROC}
\alias{reportROC}
\title{An Easy Way to Report ROC Analysis}
\description{
Provides an easy way to report the results of ROC analysis, including:
  1. an ROC curve. 2. the value of Cutoff, SEN (sensitivity), SPE (specificity),
  AUC (Area Under Curve), AUC.SE (the standard error of AUC),
  PLR (positive likelihood ratio), NLR (negative likelihood ratio),
  PPV (positive predictive value), NPV (negative predictive value).
}
\keyword{ROC analysis}
\usage{
reportROC(gold,predictor,important,plot)
}
\arguments{
  \item{gold}{0/1 or binomial factor, the 'gold standard'}
  \item{predictor}{numeric, the predictor variable}
  \item{important}{'se' or 'sp', some Youden index maybe the same, and 'important' is to indicate which is more important between sensitivity and specificity}
  \item{plot}{logic, whether to plot the ROC curve with specific style}
}
\value{
  \item{Curoff}{cutoff}
  \item{AUC}{Area Under Curve, AUC}
  \item{AUC.SE}{the standard error of AUC}
  \item{AUC.low,AUC.up}{the 95 percent CI of AUC}
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

more functions will be included in 'reportROC' package!
}
\seealso{
nothing
}
\examples{
data(aSAH)
reportROC(gold=aSAH$outcome,predictor=aSAH$s100b,important="se",plot=TRUE)
}

