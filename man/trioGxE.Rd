\name{trioGxE}
\alias{trioGxE}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Generalized additive model estimation of gene-environment interaction using data from case-parent trios
}
\description{
\code{trioGxE} estimates statistical interaction (GxE) between a single nucleotide polymorphism (SNP) 
and a continuous environmental or non-genetic attributes  
in case-parent trio data by fitting a generalized additive 
model (GAM) using a penalized iteratively re-weighted least squares algorithm.
}
\usage{
trioGxE(data, pgenos, cgeno, cenv, 
        penmod = c("codominant","dominant","additive","recessive"), 
        k = NULL, knots = NULL, sp = NULL, lsp0 = NULL, lsp.grid = NULL, 
        control = list(maxit = 100, tol = 1e-07, trace = FALSE), 
        testGxE = FALSE, return.data = TRUE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{data}{
a data frame with columns for parental genotypes, child genotypes and  
child environmental/non-genetic attribute. 
See `Details' section for the required format.}
\item{pgenos}{
a length-2 vector of character strings specifying the names of the two columns 
in \code{data} that hold parental genotypes.}
\item{cgeno}{
a character string specifying the name of the column in \code{data} that holds 
the child genotypes.
}
\item{cenv}{
a character string specifying the name of the column in \code{data} that holds 
the non-genetic attribute being examined for interaction with genotype.
}
\item{penmod}{
the penetrance mode of the genetic and interaction effects: \code{"codominant"} (default), 
\code{"dominant"}, \code{"additive"}, or \code{"recessive"}.}
\item{k}{
an optional vector or single value specifying the 
desired number(s) of knots to be used for the cubic spline basis construction of smooth function(s) representing GxE.
When \code{penmod="codominant"}, a length-2 vector with positive integers must be provided to specify the 
numbers of knots (or basis dimensions) for the two interaction functions. 
Otherwise, a single positive integer must be provided. 
The minimum value for each integer is 3.
The default basis dimension is either \code{k=c(5,5)} or \code{k=5}. 
See `Details' section for more information.}
\item{knots}{
knot positions for the cubic spline basis construction. 
When \code{penmod="codominant"}, a list of two vectors must be provided. 
For the other penetrance modes, a single vector must be provided.
When \code{NULL} (default), knots will be placed at equally-spaced quantiles of the distribution 
of \eqn{E} within trios from appropriate parental mating types. 
If both \code{knots} and \code{k} are provided, the argument \code{k} is ignored.
See `Details' section for more information. 
}
\item{sp}{
smoothing parameters for the interaction functions.
When \code{penmod="codominant"}, a vector with two non-negative 
numbers must be provided.
Otherwise, a single non-negative number must be provided.
When \code{NULL} (default), a double (under the co-dominant mode) 
or a single (under a non-co-dominant mode) 
1-dimensional grid search finds the optimal smoothing parameter values.
}
\item{lsp0}{
an optional length-2 numeric vector or a single numeric value used for 
choosing trial values of log smoothing parameters in the grid search for 
the optimal smoothing parameters.
When {\code{NULL}} (default), \code{trioGxE} takes the log of smoothing parameter estimates 
obtained by applying a likelihood approach that makes inference of GxE conditional 
on the parental genotypes, non-genetic attribute and partial information on child genotypes.}
\item{lsp.grid}{
trial values of log smoothing parameters used in the grid search for smoothing parameters. 
When \code{penmod=} \code{"codominant"}, a list of two vectors of lengths \eqn{\geq 2}
must be provided. As the vector is longer, the grid becomes more refined. 
When the penetrance mode is not co-dominant, a single vector must be provided.
When \code{lsp.grid=NULL} (default), the function take the vectors of length 6 obtained 
by using the truncated normal distributions constructed based on \code{lsp0}.}
\item{control}{
a list of convergence parameters for the penalized iteratively re-weighted 
least squares (PIRLS) procedure: 
	\tabular{ccrl}{
      	\tab\tab\code{maxit}:\tab positive integer giving the maximal number of PIRLS iterations \cr
      	\tab\tab\code{tol}:\tab positive convergence tolerance in terms the relative difference in penalized \cr 
	\tab\tab\tab deviances (\code{pdev}) between iterations: 
		   \eqn{{|\code{pdev} - \code{pdev}_{old}|}/{(|\code{pdev}| + 0.1)} < \code{tol}} \cr
      	\tab\tab\code{trace}:\tab logical indicating if output should be produced 
	     	for each PIRLS iteration.
	}
}
\item{testGxE}{
a logical specifying whether the fitting is for testing interaction. 
Default is \code{FALSE}. User should not modify this argument.
}
\item{return.data}{
a logical specifying whether the original data should be returned. 
If \code{TRUE} (default), the data is returned. 
}
\item{\dots}{ 
sets the arguments for \code{control}, which includes \code{maxit}, \code{tol} or \code{trace}.
}
}
\details{
\code{trioGxE} fits data from case-parent trios to a GAM with smooth functions 
representing gene-environment interaction (\eqn{G \times E}). 

The \code{data} input must be a data frame containing the following 4 columns (of any order):
\tabular{rll}{
[,1] \tab number (0, 1 or 2) of copies of a putative risk allele carried by the mother \cr
[,2] \tab number of copies of a putative risk allele carried by the father \cr
[,3] \tab number of copies of a putative risk allele carried by the affected child \eqn{(G)} \cr
[,4] \tab value of a continuous environmental/non-genetic attribute measured on the child \eqn{(E)}
} 
The function \code{trioGxE} does basic error checking to ensure that  
only the trios that are consistent with Mendelian segregation law 
with complete genotype, environment and parental genotype information. 
The function determines which trios are from \emph{informative} parental mating types.
An informative parental pair has at least one heterozygote; such parental pair can 
have offspring that are genetically different.
Under the assumption that the parents transmit the alleles to their child under 
Mendel's law, with no mutation, there are three types of informative mating types 
\eqn{G_p=1,2,3}: 
\tabular{rl}{
\eqn{G_p=1}:\tab if one parent is heterozygous, and the other parent is homozygous for 
the non-risk allele \cr
\eqn{G_p=2}:\tab if one parent is heterozygous, and the other parent is homozygous for 
the risk allele \cr
\eqn{G_p=3}:\tab if both parents are heterozygous
} 

Since GxE occurs when genotype relative risks vary with non-genetic attribute values \eqn{E=e},
GxE inference is based on the attribute-specific genotype relative risks, GRR\eqn{_h(e)}, expressed as 
\deqn{
{\rm GRR}_h(e) = 
\frac{P(D=1|G=h,E=e)}{P(D=1|G=h-1,E=e)} = \exp{(\gamma_h + f_h(e))}, ~~h=1,2,
}
where \eqn{D=1} indicates the affected status, 
\eqn{\gamma_1} and \eqn{\gamma_2} represent genetic main effect, 
and \eqn{f_1(e)} and \eqn{f_2(e)} are unknown smooth functions. 
The functions \eqn{f_h(e)} represent GxE since GRRs vary with \eqn{E=e} only when 
\eqn{f_1(e)\neq 0} or \eqn{f_2(e)\neq 0} vary with \eqn{E}. 
The expressions are followed by assuming a log-linear model of disease risk in Shin et al. (2010).

Under the co-dominant penetrance mode (i.e., \code{penetrance="codominant"}),
GRR\eqn{_1(e)} and GRR\eqn{_2(e)} are estimated using the information 
in the trios from the informative mating types \eqn{G_p=1}, 3 and those from \eqn{G_p=2}, 3, respectively.
Under a non-co-dominant penetrance mode, only one GRR function, GRR\eqn{(e)=\gamma+f(e)}, 
is estimated from an appropriate set of informative trios.
Under the dominant penetrance mode (i.e., \code{penetrance="dominant"}), 
because GRR\eqn{_2(e)} is 1 (i.e., \eqn{\gamma_2=0} and \eqn{f_2(e)=0}), 
GRR\eqn{(e)\equiv} GRR\eqn{_1(e)} is estimated based on the trios from \eqn{G_p=1} and 3.
Under the recessive penetrance mode (i.e., \code{penetrance="recessive"}), 
because GRR\eqn{_1(e)} is 1 (i.e., \eqn{\gamma_1=0} and \eqn{f_1(e)=0}),
GRR\eqn{(e)\equiv} GRR\eqn{_2(e)} is estimated based on the trios from \eqn{G_p=2} and 3.
Under the multiplicative or log-additive penetrance model (\code{penetrance="additive"}),
since GRR\eqn{_1(e)} = GRR\eqn{_2(e)} (i.e., \eqn{\gamma_1=\gamma_2} and \eqn{f_1=f_2}), 
GRR\eqn{(e)\equiv} GRR\eqn{_1(e) \equiv} GRR\eqn{_2(e)} is estimated based on all informative trios.

The interaction functions are approximated by cubic regression spline functions defined by 
the knots specified through the arguments \code{k} and \code{knots}. 
For each interaction function, at least three knots are chosen within the range of 
the observed non-genetic attribute values.
Under the co-dominant mode, \code{k[1]} and \code{k[2]} knots, respectively, located at 
\code{knots[[1]]} and \code{knots[[2]]} are used to construct the basis for 
\eqn{f_1(e)} and \eqn{f_2(e)}, respectively. 
By default, a total of 5 knots are placed for each interaction function:
three interior knots at the 25th, 50th and 75th quantiles and two boundary knots 
at the endpoints of the data in trios from \eqn{G_p=1} or \eqn{3}, for \eqn{f_1(e)}, 
and in trios from \eqn{G_p=2} or \eqn{3}, for \eqn{f_2(e)}. 
Similarly, under a non-co-dominant penetrance mode, when \code{knots=NULL}, 
\code{k} knots are chosen based on the data in trios from 
\eqn{G_p=1} and \eqn{3} (dominant mode); 
in trios from all informative mating types (log-additive mode); 
and in trios from \eqn{G_p=2} or \eqn{3} (recessive mode). 
A standard model identifiability constraint is imposed on each interaction function, 
which involves the sum of the interaction function over all observed attribute values of cases 
in the appropriate set of informative trios.

For smoothing parameter estimation, \code{trioGxE} finds the optimal values using either a double (if co-dominant) or 
a single 1-dimensional grid search constructed based on the arguments \code{lsp0} and \code{lsp.grid}.
When \code{lsp0 = NULL}, \code{trioGxE} takes the log smoothing parameter estimates obtained 
from a likelihood approach that makes inference of GxE conditional 
on parental mating types, non-genetic attribute and partial information on child genotypes (Duke, 2007).
When \code{lsp.grid = NULL}, \code{trioGxE} takes the following 6 numbers to be the trial values of 
the log-smoothing parameter for each interaction function:  
-20 and 20, \code{lsp0} and the quartiles of the truncated normal distributions 
constructed based on \code{lsp0}. 
At each trial value in \code{lsp.grid}, the prediction error criterion function, 
UBRE (un-biased risk estimator, is minimized to find the optimal smoothing parameter.
For more details on how to estimate the smoothing parameters, see Appendix B.3 in Shin (2012).

For variance estimation, \code{trioGxE} uses a Bayesian approach (Wood, 2006); 
the resulting Bayesian credible intervals have been shown to have good frequentist coverage probabilities 
as long as the bias is a relatively small fraction of the mean squared error (Marra and Wood, 2012) 
}
\value{
   \item{coefficients}{a vector holding the spline coefficient estimates for 
\eqn{\hat{f}_{1}} and/or \eqn{\hat{f}_{2}}. 
The length of the vector is equal to the total number of knots used for 
constructing the bases of the interaction curves. 
For example, under the default basis dimension with co-dominant penetrance mode, 
the vector has size 10 (i.e., 5 for \eqn{f_1} and the other 5 for \eqn{f_2}).}
   \item{control}{a list of convergence parameters used for the penalized iteratively re-weighted 
least squares (PIRLS) procedure}
   \item{data}{original data passed to \code{trioGxE()} as an argument: returned when \code{return.data=TRUE}}
   \item{edf}{a vector of effective degrees of freedom for the model parameters 
(see page 171 in Wood, 2006 for details).}
   \item{Gp}{a vector containing the values of parental mating types \eqn{G_p} (see `Details' for 
the definition)}
   \item{lsp0}{log smoothing parameter(s) used in the grid search. 
Not returned (i.e., \code{NULL}) when smoothing parameters were not estimated but provided by the user.}
   \item{lsp.grid}{trial values of the log smoothing parameter(s) used in the grid search.
Not returned when smoothing parameters were not estimated but provided by the user.}
   \item{penmod}{the penetrance mode under which the data were fitted.}
   \item{qrc}{a list containing the QR-decomposition results used for imposing the identifiability constraints. 
See \code{\link{qr}} for the list of values.}
   \item{smooth}{a list with components:
\tabular{ll}{
\code{model.mat}:\tab The design matrix of the problem. The number of rows is equal to \cr
\tab \eqn{n_1+n_2+2n_3}, where \eqn{n_m} is the number of \eqn{m}th informative \cr
\tab mating types. The number of columns is equal to the size of \code{coefficient}. \cr
\code{pen.mat}:\tab penalty matrix with size equal to the size of \code{coefficient}.\cr
\code{bs.dim}:\tab number of knots used for basis construction. \cr
\code{knots}:\tab knot positions used for basis construction.
} 
}
   \item{sp}{optimal smoothing parameter values calculated from UBRE optimization or smoothing parameter values provided by the user.}
   \item{sp.user}{logical whether or not the smoothing parmeter values were provided by the user. 
If \code{FALSE}, \code{sp} contains the smoothing parameter values estimated by the UBRE optimization.} 
   \item{terms}{list of character strings of column names in \code{data} corresponding to 
the child genotypes, parental genotypes and child's non-genetic attributes.}
   \item{triodata}{Formatted data passed into the internal fitting functions. To be used in \code{\link{test.trioGxE}} function.}
   \item{ubre}{the minimum value of the un-biased risk estimator (UBRE), measure of predictability for the interaction function 
estimators \eqn{\hat{f}_{1}} or \eqn{\hat{f}_{2}}.
Not returned when smoothing parameters were not estimated but provided by the user.}
   \item{ubre.val}{a list or a vector of ubre values corresponding to the trial values 
of smoothing parameter(s) in \code{lsp.grid}.}
   \item{Vp}{Bayesian posterior variance-covariance matrix for the coefficients. 
The size the matrix is the same as that of \code{coefficient}.}
}
\references{
Duke, L. (2007): 
\emph{A graphical tool for exploring SNP-by-environment interaction in
	case-parent trios}, 
M.Sc. Thesis, Statistics and Actuarial Science, Simon Fraser University: \cr
URL http://www.stat.sfu.ca/content/dam/sfu/stat/alumnitheses/Duke-2007.pdf.\cr

Marra, G., Wood, S.N. (2012).
Coverage properties of confidence intervals for generalized additive model components. 
\emph{Scand J Stat}, \bold{39}: 53-74. \cr

Shin, J.-H., McNeney, B. and Graham, J. (2010). 
On the use of allelic transmission rates for assessing gene-by-environment interaction 
in case-parent trios. 
\emph{Ann Hum Gen}, \bold{74}: 439-51. \cr

Shin, J.-H. (2012): 
\emph{Inferring gene-environment interaction from case-parent trio data: evaluation of and 
adjustment for spurious \eqn{G\times E} and development of a data-smoothing method to uncover 
true \eqn{G\times E}}, 
Ph.D. Thesis, Statistics and Actuarial Science, Simon Fraser University: 
URL https://theses.lib.sfu.ca/sites/all/files/public_copies/etd7214-j-shin-etd7214jshin.pdf.\cr

Wood, S. (2006): \emph{Generalized Additive Models: An Introduction with R}, Boca Raton, FL: 
Chapman & Hall/CRC.}

\author{
Ji-Hyung Shin <shin@sfu.ca>, 
Brad McNeney <mcneney@sfu.ca>, 
Jinko Graham <jgraham@sfu.ca>
}
\seealso{
\code{\link{trioSim}}, \code{\link{plot.trioGxE}}, \code{\link{test.trioGxE}}
}
\examples{

## fitting a co-dominant model
data(hypoTrioDat)
simfit <- trioGxE(data=hypoTrioDat,pgenos=c("parent1","parent2"),cgeno="child",cenv="attr",
                  k=c(5,5),knots=NULL,sp=NULL)

## fitting a dominant model to the hypothetical data
simfit.dom <- trioGxE(data=hypoTrioDat,pgenos=c("parent1","parent2"),cgeno="child",cenv="attr",
                      penmod="dom",k=5,knots=NULL,sp=NULL)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{methods} \keyword{models} \keyword{smooth} \keyword{regression}
