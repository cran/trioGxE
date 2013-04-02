\name{trioSim}
\alias{trioSim}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Simulate informative case-parent trios
}
\description{
\code{trioSim()} simulates parental genotypes, child genotypes, 
environmental attribute and sub-population membership on affected trios with 
\emph{informative} mating types
from a stratified population. All genotypes are at a test locus that is linked to a causal locus.
}
\usage{
trioSim(n, popfs, hapfs, edists, recomb = 0, riskmod, batchsize = 1000)
}
\arguments{
  \item{n}{
A desired number of informative case-parent trios to simulate.
}
  \item{popfs}{
A vector of sub-population frequencies whose length is equal to the number of sub-populations.
}
  \item{hapfs}{
A list comprised of vectors of haplotype frequencies. 
One haplotype frequency vector for each sub-population. 
See Details for the assumed order of haplotypes. 
}
  \item{edists}{
A list comprised of functions to simulate the environment attribute. 
One simulation function for each sub-population.
}
  \item{recomb}{
Recombination frequency between the test and causal locus. 
Currently not implemented. 
The function will stop execution if a non-zero value is specified.
}
  \item{riskmod}{
A function to evaluate the risk (probability) of disease. 
The function should take two arguments. The first is the child's genotype, 
and the second is the environmental attribute.
}
 \item{batchsize}{
Size of the batches of trios to simulate. See Details for more information.
}
}
\details{
The function simulates trios from a stratified population.
Population stratification is controlled by the user's choice of 
sub-population sizes, haplotype frequencies in each sub-population and the 
distribution of the environmental attribute in each sub-population.
Given sub-population sizes, the degree of population stratification increases with 
greater differences in the distributions of the haplotype frequency and the 
environmental attribute among sub-populations.

The function first simulates sub-population membership for each trio using 
the sub-population frequencies supplied by the user in the argument \code{popfs}.
Conditional on sub-population, parental haplotypes Hp are simulated assuming 
Hardy-Weinberg proportions using the subpopulation-specific haplotype frequencies 
in the argument \code{hapfs}. 
Haplotype frequencies should be in the order N0, N1, R0, R1, where N and R denote 
the non-risk and risk alleles at the causal locus, and 0 and 1 denote the non-index 
and index alleles at the test locus. 

To save computation time, we only considered informative parental mating types 
by simulating one parent from the conditional distribution given that the parent 
is heterozygous and simulating the other parent without any restrictions. 
Conditional on parental haplotypes, child haplotypes are sampled according to Mendel's laws.
From the sampled haplotypes of the parents and children, their genotypes for the causal 
and test loci are extracted. 

Assuming conditional independence between the gene and the environmental attribute given 
sub-population, the environmental attribute for each trio is simulated conditional on 
sub-population using the subpopulation-specific simulation functions in the argument \code{edists}.
Finally, disease status is simulated according to the risk model in the argument \code{riskmod};
only those trios with affected children are retained.

To speed up computation, the rejection sampling of trios is done in batches of size \code{'batchsize'} 
until a desired number of affected trios is obtained. 
In simulation studies we have performed, choosing \code{batchsize} on the order of 1/3 the 
desired number of trios appeared to be the fastest.
}
%
\value{
A data frame with columns
\item{parent1}{Test locus genotypes for one parent (heterozygous) coded as 0, 1, 2 
representing the number of copies of the index allele.}
\item{parent2}{Test locus genotype for the other parent.}
\item{child}{Test locus genotypes for the child.}
\item{subpop}{Sub-population membership for the trio. 
Sub-populations are numbered \eqn{0,1,...,k-1}, where \eqn{k} is the number of sub-populations.}
\item{attr}{The environmental attribute.}
}
%\references{
%% ~put references to the literature/web site here ~
%}
\author{
Ji-Hyung Shin <shin@sfu.ca>, 
Brad McNeney <mcneney@sfu.ca>, 
Jinko Graham <jgraham@sfu.ca>
}
%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{trioGxE}}, \code{\link{plot.trioGxE}}, \code{\link{test.trioGxE}}
}

\examples{
# Generate case-parent trio from a population composed of 
# two equal sized subpopulations.

# Set up list of functions to sample from each E distribution

e1<-function(n) {
return(rnorm(n,mean=(-0.8),sd=sqrt(1-.8^2)))
}
e2<-function(n) {
return(rnorm(n,mean=(0.8),sd=sqrt(1-.8^2)))
}

# Set up haplotype frequency distributions in the two subpopulations: 
# The first subpopulation has the risk allele frequency of 0.1, where as
# the second subpopulation's frequency is 0.9. 

# Set up risk model function.
## Simulate informative case-parent trios under additive linear GxE with a negative slope

riskmod<-function(G,E) {
n<-length(G)
# Baseline risk. Affects disease prevalence.
# The higher the prevalence, the less time wasted
# rejecting unaffected trios.
k<-(-2) 

betaG<-log(3)/2

# Interaction
betaGE<-(-0.1)

# quadratic GxE
rr<-exp(k+betaG*G + betaGE*G*E)
rr[rr>1]<-1 # It is up to the user to make sure there are
# no probabilities greater than one.

D<-rbinom(n=n,size=1,prob=rr)
return(D)
}


# Simulate trio data under haplotype-environment dependence 
# when marker locus is causal locus.
# allele frequency in subpop 0 is 0.1, allele frequency in subpop 1 is 0.9.
hapf1=c(0.9, 0, 0, 0.1)
hapf2=c(0.1, 0, 0, 0.9)
simdat.HEdep<-trioSim(n=3000,popfs=c(0.5,0.5),riskmod=riskmod,
                      edists=list(e1,e2),hapfs=list(hapf1,hapf2), 
                      recomb=0,batchsize=1000)

# Simulate trio data under haplotype-environment independence
# when marker locus is causal locus.
# allele frequency in subpop 0 and subpop 1 is 0.1.
hapf1=hapf2=c(0.9, 0, 0, 0.1)
simdat.HEindep<-trioSim(n=3000,popfs=c(0.5,0.5),riskmod=riskmod,
                        edists=list(e1,e2),hapfs=list(hapf1,hapf2), 
                        recomb=0,batchsize=1000)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ datagen }
% __ONLY ONE__ keyword per line
