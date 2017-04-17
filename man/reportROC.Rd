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
  \item{SEN}{sensitivity}
  \item{SPE}{specificity}
  \item{AUC}{Area Under Curve, AUC}
  \item{AUC.SE}{the standard error of AUC}
  \item{PLR}{positive likelihood ratio}
  \item{NLR}{negative likelihood ratio}
  \item{PPV}{positive predictive value}
  \item{NPV}{negative predictive value}
}
\author{
Zhicheng Du<dgdzc@hotmail.com>, Yuantao Hao<haoyt@mail.sysu.edu.cn>
}
\note{
Please feel free to contact us, if you have any advice and find any bug!

Update description:

more functions will be included in 'reportROC' package!
}
\seealso{
nothing
}
\examples{
data(aSAH)
reportROC(gold=aSAH$outcome,predictor=aSAH$s100b,important="se",plot=TRUE)
}

