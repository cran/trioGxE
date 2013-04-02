\name{test.trioGxE}
\alias{test.trioGxE}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Test of gene-environment interaction between a SNP and a continuous non-genetic 
covariate from case-parent trio data.
}
\description{
Performs permutation test of gene-environment interaction based on 
the associated penalized maximum likelihood estimates
obtained by fitting a generalized additive model to case-parent trio data.
}
\usage{
test.trioGxE(object, data = NULL, nreps, level = 0.05, early.stop = FALSE, 
             fix.sp = FALSE, output = NULL, return.data = FALSE, 
             return.object = FALSE, ...)
}
\arguments{
  \item{object}{
A returned object from \code{\link{trioGxE}} function. 
When \code{NULL}, a data set of case-parent trios must be provided 
(through `\code{data}' argument).
}
  \item{data}{
Trio data to be passed into \code{\link{trioGxE}} when `\code{object}' is not provided.
}
  \item{nreps}{
Desired number of permutation replicates.
}
  \item{fix.sp}{
When \code{TRUE}, the approximated null distribution of the test statistic is obtained by 
computing the test statistic by fitting each simulated data set under 
fixed values of the smoothing parameters. 
When \code{FALSE} (Default), the null distribution is obtained by fitting each simulated data set 
while simultaneously estimating the smoothing parameter values. 
}
  \item{level}{
Desired significance level for the test.
}
  \item{early.stop}{
When \code{TRUE}, sampling is terminated early when the number of test statistics 
that are more extreme than or as extreme as the observed test statistic equals 
\code{nreps*level}. 
\eqn{\code{nreps} \times \code{level}} values larger 
than the observed test statistic are obtained.
}
  \item{output}{
A character string specifying the name of the output file that writes the values of 
the test statistics calculated for the actual and simulated data set. 
When \code{NULL} (Default), no written output file is produced.
}
  \item{return.data}{
When \code{TRUE}, the original data set is returned.
}
  \item{return.object}{
When \code{TRUE}, the fitting object for the original data set is returned.
}
  \item{\dots}{
Arguments passed to \code{\link{trioGxE}}: when \code{data} is provided, 
instead of \code{trioGxE} class object, parameters of \code{\link{trioGxE}} 
must be provided through \code{...}.
}
}

\details{
Suppose \eqn{k_1} and \eqn{k_2} are the numbers of knots used to represent 
the interaction functions \eqn{f_1} and \eqn{f_2}, respectively, via 
cubic regression spline functions.
Let \eqn{\bm{c}_1 = (c_{11},...,c_{1K_1-1})^\prime} and 
\eqn{\bm{c}_2 = (c_{21},...,c_{2K_2-1})^\prime} 
are the spline coefficient vectors for \eqn{f_1} and \eqn{f_2}
that satisfy model identifiability constraints.

The function \code{test.trioGxE} calculates test statistic \eqn{T}, 
%under the co-dominant penetrance mode, 
\deqn{ T = t(\hat{\bm{c}}){\rm V}^{-1}(\bm{c})\hat{\bm{c}}, }
where \eqn{\bm{c}=(\bm{c}_1^{\prime},\bm{c}_2^{\prime})^{\prime}}
and \eqn{V_c} is a square matrix of size \eqn{(k_1+k_2-2)}, 
formed by extracting the rows and columns, corresponding to the spline
coefficients from the Bayesian posterior variance-covariance matrix 
calculated in \code{\link{trioGxE}}.

If the actual data were fitted under the co-dominant penetrance mode 
(i.e., \code{object$penmod="codominant"}), 
the test statistic \eqn{T} represents an \emph{overall} test of GxE, where 
\deqn{{\rm H}_0: \bm{c}=0.}
Depending on the context, an investigator may also want to perform individual tests:
\eqn{{\rm H}_{01}: \bm{c}_1 = \bm{0}} and \eqn{{\rm H}_{02}: \bm{c}_2 = \bm{0}}.
For example, when the null hypothesis is rejected, the user may want to know which of 
the two interaction function is not zero (i.e., which curve is not flat). 
For the \emph{individual} tests, \code{test.trioGxE} calculates the 
permutation p-values based on the Monte-Carlo distributions of the individual 
test statistics \eqn{T_1} and \eqn{T_2}, where 
\deqn{T_h = t(\hat{\bm{c}}_h){\rm V}^{-1}(\bm{c}_h)\hat{\bm{c}}_h, h=1,2.}

Under the dominant, log-additive (multiplicative) or recessive penetrance model,  
\eqn{T} can be viewed as an individual test 
since \eqn{\bm{c}_2=\bm{0}}, \eqn{\bm{c}_1=\bm{c}_2} and \eqn{\bm{c}_1=\bm{0}}, 
respectively, under the dominant, log-additive and recessive models.
For example, under the dominant penetrance model, \eqn{T\equiv{T_1}} because 
\eqn{\bm{c}_2=\bm{0}}, and \eqn{T_2=0}. 

As the analysis is conditional on parental genotypes, 
the distribution of the test statistic under \eqn{{\rm H}_0} is calculated by 
shuffling the column that holds the values of the non-genetic covariate within mating types. 
This can be justifiable based on the fact that under no interaction, 
the SNP and the non-genetic covariate are independent of each other within 
a random affected trio when they are independent within a trio from the general population 
(Umbach and Weinberg, 2000).

The distribution of the test statistics can be obtained in two ways: 
either under fixed smoothing parameters (\code{fixed.sp=TRUE}) 
or under varying smoothing parameters (\code{fixed.sp=FALSE}).
Under the fixed smoothing parameters, the penalized iteratively re-weighted least squares
procedure is performed for each simulated data set under the same smoothing parameter values. 
Under varying smoothing parameters, smoothing parameters are estimated for each 
simulated data set. 
Therefore, the test under \code{fixed.sp=FALSE} accounts for the extra uncertainty introduced by 
the smoothing parameter estimation.

To save computation time, the user can use `early-termination' option (Besag and Clifford, 1991).
Under this option, sampling is terminated when the number of the simulated data sets reaches 
\code{nreps*{level}} \eqn{<} \code{nreps} 
when the evidence is not strong enough to reject the null hypothesis at the given significance level (\code{level}). 
For example, if the user specifies \code{nreps=1000} and \code{level=0.05}, the test terminates when the number of 
data sets that have test statistic values that are more extreme than or as extreme as the observed test statistic value reaches 50.}

\value{
\item{GxE.test}{
Either a 3- or 1-column matrix. 
When the actual data was fitted under a co-dominant penetrance mode 
(i.e., \code{object$penmod} = \code{"codominant"}), 
a 3-column matrix is returned, where the first column holds the values for \eqn{T} 
for the original and the generated data sets, 
and the second and third columns hold the values of 
\eqn{T_1} and \eqn{T_2}, respectively for the same data sets.
When the actual data was fitted under a non-co-dominant penetrance mode (e.g., dominant), 
\code{GxE.test} is retuned as a matrix with a sigle column holding \eqn{T}.
}
\item{p.value}{ If \code{object$penmod} = \code{"codominant"}, 
it is returned as a vector holding three values, 
where the first value indicates the overall p-value obtained from the distribution of \eqn{T}, 
and the other two values indicate the individual p-values obtained from the 
distributions of \eqn{T_1} and \eqn{T_2}. 
Under \code{object$penmod} is \code{dominant}, \code{additive} or \code{recessive}, 
it is returned as a single p-value calculated based on \eqn{T}. }
}
\references{
Umbach, D. and Weinberg, C. (2000). The use of case-parent triads to study joint effects of genotype 
and exposure. 
\emph{Am J Hum Gen}, \bold{66}:251-61. \cr

Besag, J. and P. Clifford (1991). Sequential Monte Carlo p-values. \emph{Biometrika}, \bold{78}:301-304.
}
\author{
Ji-Hyung Shin <shin@sfu.ca>, 
Brad McNeney <mcneney@sfu.ca>, 
Jinko Graham <jgraham@sfu.ca>
}
\seealso{
\code{\link{trioGxE}}, \code{\link{plot.trioGxE}}, \code{\link{trioSim}}
}
\examples{
data(hypoTrioDat)
example.fit <- trioGxE(hypoTrioDat, pgenos = c("parent1","parent2"), cgeno = "child",
                       cenv = "attr",penmod="codominant", k=c(5,5))
# A toy example with 'few' permutation replicates
example.test <- test.trioGxE(example.fit, nreps=10, early.stop = FALSE, 
                              output=NULL)
 

\dontrun{

## More proper examples of permutation tests with 5000 replicates

## Example1: does not generate an output file containing test statistic values
example.test1 <- test.trioGxE(example.fit, nreps=5000, early.stop = TRUE, 
                              output=NULL)
## Example 2: generates an output file 'myoutput.out' containing test statistic values 
example.test2 <- test.trioGxE(example.fit, nreps=5000, early.stop = TRUE, 
                              output="myoutput.out")
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{methods}
