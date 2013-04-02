\name{plot.trioGxE}
\alias{plot.trioGxE}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Graphical display of gene-environment interaction between a SNP and a continuous non-genetic 
attribute in case-parent trio data
}
\description{
The function \code{plot.trioGxE} uses the calculations made in \code{\link{trioGxE}} and 
plots the point- and interval-estimates of gene-environment interaction 
between a single nucleotide polymorphism (SNP) and a continuously varying environmental or non-genetic 
covariate in case-parent trio data.}
\usage{
\method{plot}{trioGxE}(x, se = TRUE, seWithGxE.only = TRUE, ylim = NULL, yscale = TRUE, 
     xlab = NULL, ylab = NULL, rugplot = TRUE, \ldots)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
A returned object produced by \code{\link{trioGxE}} function.
}
  \item{se}{
A logical or a positive number. 
When \code{TRUE} (default), upper and lower lines are added to the plots 
at 2 standard errors above and below the fitted values of the interaction functions.
When it is a positive number, lines are added at \code{se} standard errors above
and below the fitted interaction values.
When \code{FALSE}, no standard error lines are plotted.
}
  \item{seWithGxE.only}{
If \code{TRUE}, the associated standard errors reflect the uncertainty in the estimates of the 
gene-environment interaction functions only. 
If \code{FALSE}, the standard errors include the uncertainty in the genetic main effect estimates. 
}
  \item{ylim}{
Either a list holding two length-2 numeric vectors that give different y-coordinate ranges for 
the two plots, or a single length-2 vector that gives equal y-coordinate ranges for 
both plots.
}
  \item{yscale}{
If \code{TRUE} (default), the same y-axis scale is chosen for each plot. 
Ignored if \code{ylim} is supplied.
}
  \item{xlab}{
An optional string setting the title for the x-axis. 
}
  \item{ylab}{
An optional string setting the title for the y-axis. 
}
  \item{rugplot}{
Logical indicating whether to add rug representation of the data to the plots. 
Default (\code{TRUE}) adds rugs. 
}
  \item{\dots}{
Further graphical parameters passed to \code{\link{plot}}, such as \code{col}, \code{lwd}, etc.
}
}
\details{
The function produces two plots in a 2 x 1 layout. 
The first plot in the left panel displays the estimated gene-environment interaction (GxE) curve 
related to \eqn{{\rm GRR}_1}, the genotype relative risks (GRRs) among the individuals with one copy 
of the putative risk allele compared to those with zero copies. 
The right panel displays the estimated GxE curve related to \eqn{{\rm GRR}_2}, 
the GRRs among the individuals with two copies of the risk allele compared to those with one copy.

When \code{object$penmod="codominant"} (with \code{se=TRUE}), 
confidence intervals are plotted for both interaction curves that are related to \eqn{{\rm GRR}_1} and \eqn{{\rm GRR}_2}. 
When \code{object$penmod="dominant"}, the confidence intervals are plotted only in the left panel, 
but not in the right panel because \eqn{{\rm GRR}_2} is not estimated but set to be 1 under this penetrance mode.
Similarly, when \code{object$penmod="recessive"}, the confidence intervals are plotted only in the right panel, 
but not in the right panel because \eqn{{\rm GRR}_1} is not estimated but set to be 1 under this penetrance mode.
When \code{object$penmod="additive"}, equivalent confidence intervals 
are plotted in both panels, which display equivalent fitted curves. 
This is because \eqn{{\rm GRR}_1} and \eqn{{\rm GRR}_2} are set to be equivalent 
under the log-additive or multiplicative penetrance mode.

When \code{se} is \code{TRUE} or a positive number, standard error lines are plotted 
based on the calculations of the Bayesian posterior variance estimates of the generalized 
additive model parameters for GRRs (Wood, 2006). 
}
\references{
Shin, J.-H. (2012): \emph{Inferring gene-environment interaction from case-parent trio data: evaluation of and adjustment for spurious \eqn{G\times E} and development of a data-smoothing method to uncover true \eqn{G\times E}}, Ph.D. thesis, Statistics and Actuarial Science, Simon Fraser University: URL https://theses.lib.sfu.ca/sites/all/files/public_copies/etd7214-j-shin-etd7214jshin.pdf.\cr

Wood, S. (2006): \emph{Generalized Additive Models: An Introduction with R}, Boca Raton, FL: Chapman & Hall/CRC.
}
\author{
Ji-Hyung Shin <shin@sfu.ca>, 
Brad McNeney <mcneney@sfu.ca>, 
Jinko Graham <jgraham@sfu.ca>
}
%% ~Make other sections like Warning with \section{Warning }{....} ~
\seealso{
\code{\link{trioGxE}}, \code{\link{test.trioGxE}}, \code{\link{trioSim}}

}
\examples{
data(hypoTrioDat)

## fitting a co-dominant model to the hypothetical data
simfit <- trioGxE(data=hypoTrioDat,pgenos=c("parent1","parent2"),cgeno="child",cenv="attr",
                  k=c(5,5),knots=NULL,sp=NULL)

## produce the graphical display of the point- and interval-estimates of GxE curve 
plot.trioGxE(simfit) # or just plot(simfit)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{methods} \keyword{models} \keyword{smooth} \keyword{regression} \keyword{hplot}

