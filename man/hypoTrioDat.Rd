\name{hypoTrioDat}
\docType{data}
\alias{hypoTrioDat}
\title{Simulated data for a hypothetical example}
\description{
A simulated data set with 1000 informative case-parent trios to illustrate the \code{trioGxE} package.
}
\usage{data(hypoTrioDat)}
\format{data frame with columns:
\tabular{rlll}{
     \tab  column name \tab type \tab details \cr
[,1] \tab parent 1 \tab numeric \tab copy number of the putative risk allele carried by the first parent (0, 1, 2) \cr
[,2] \tab parent 2 \tab numeric \tab copy number of the putative risk allele carried by the other parent (0, 1, 2) \cr
[,3] \tab child \tab numeric \tab copy number of the putative allele carried by the child (0, 1, 2) \cr
[,4] \tab subpop \tab numeric \tab membership for two subpopulations 0, 1 \cr
[,5] \tab attr \tab numeric \tab non-genetic attribute value of the child}
}
\note{ 
The data contains 1000 informative case-parent trios, 
each of which has at least one heterozygous parent.
The trios were generated, using \code{\link{trioSim}} 
under recessive quadratic gene-environment interaction, 
from a stratified population composed of two equal-sized subpopulations. 
The two subpopulations have different distributions for the SNP and the non-genetic attribute. 
For the SNP, the risk allele frequencies for the first (\code{subpop=0}) and 
second (\code{subpop=1}) subpopulations are 0.1 and 0.9, respectively.
For the non-genetic attribute, the distributions in subpopulations 0 and 1 are Normal 
with subpopulation-specific means -0.8 and 0.8, respectively, and common variance 0.36.}
\seealso{\code{\link{trioSim}}}
\keyword{dataset}
