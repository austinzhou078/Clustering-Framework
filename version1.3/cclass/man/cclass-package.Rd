\name{core_class}
\alias{core_class}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Core class selection framework
}
\description{
A framework to elucidate a set of core groups of a dataset
}
\usage{
cclass = function (x, ccol = c(1:dim(x)[2]),cluster=19, iter = 100, dist="euclidean")
}
\arguments{
  \item{x}{
the input dataset which need analysing
}
  \item{ccol}{
the features of dataset which will be considered as dimensions
}
  \item{cluster}{
 maximum clusters will be compared(from 2 to 19)
}
  \item{iter}{
iterated times of k-means and fcm algorithm(default 100 times)
}
  \item{dist}{
the method will be used for fcm and pam, choice: manhattan and euclidean
}
}
\details{
This is a framework to elucidate a set of core groups of a dataset. It specifies the application of different clustering algorithms, many statistical and visualization techniques. At the beginning of this framework, different clustering algorithms are applied. Then a set of common classes is defined in order to determine the fundamental characteristics of data expressed by different groups.
}
\value{
item{cluster}{the new dataset contain the core class}
item{plot }{the biplot of classes}
item{number }{the selected number of clusters }
}
\references{
P. Kellam, X. Liu, N. Martin, C. Orengo, S. Swift, and A. Tucker,
“Comparing, contrasting and combining clusters in viral gene expression data,” inProceedings of 6th Workshop on Intelligent Data
Analysis in Medicine, 2001.

Soria, D. Garibaldi, J.M.
"A novel framework to elucidate core classes in a dataset"
Evolutionary Computation (CEC), 2010 IEEE Congress on. Date of Conference: 18-23 July 2010 Page(s): 1 - 8 Print ISBN: 978-1-4244-6909-3 INSPEC Accession Number: 11567584 Conference Location : Barcelona DOI: 10.1109/CEC.2010.5586331 Publisher: IEEE

}
\author{
Xiaotian Zhou, Chunfan Gao, Daniele.Soria
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
result = cclass(iris, ccol=c(1:4), cluster = 15, iter = 50)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{cluster}
\keyword{core}
\keyword{framework}
