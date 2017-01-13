\name{wdbc}
\alias{wdbc}
\docType{data}
\title{
Wisconsin Diagnostic Breast Cancer (WDBC)
}
\description{
The data contain measurements on cells in suspicious lumps in a women's breast.
Features are computed from a digitized image of a fine needle aspirate (FNA) of a breast mass. They describe characteristics of the cell nuclei present in the image. All samples are classsified as either \emph{benign} or \emph{malignant}.
}
\usage{data(wdbc)}
\format{
\code{wdbc} is a \code{data.frame} with 31 columns. The first column indicates wether the sample is classified as benign (\code{B}) or malignant (\code{M}). The remaining columns contain measurements for 30 features.
}
\details{
Ten real-valued features are computed for each cell nucleus: \cr

    a) radius (mean of distances from center to points on the perimeter)  \cr
	b) texture (standard deviation of gray-scale values) \cr
	c) perimeter \cr
	d) area \cr
	e) smoothness (local variation in radius lengths) \cr
	f) compactness (perimeter^2 / area - 1.0) \cr
	g) concavity (severity of concave portions of the contour) \cr
	h) concave points (number of concave portions of the contour) \cr
	i) symmetry  \cr
	j) fractal dimension ("coastline approximation" - 1) \cr

The references listed below contain detailed descriptions of
how these features are computed. 

The mean, standard error, and "worst" or largest (mean of the three
largest values) of these features were computed for each image,
resulting in 30 features.
}
\note{
This breast cancer database was obtained from the University of Wisconsin  Hospitals, Madison from Dr. William H. Wolberg.
}
\source{
\url{https://archive.ics.uci.edu/ml/datasets/Breast+Cancer+Wisconsin+(Diagnostic)} \cr \cr
Bache, K. & Lichman, M. (2013). UCI Machine Learning Repository. Irvine, CA: University of California, School of Information and Computer Science. 
}
\references{
O. L. Mangasarian and W. H. Wolberg: "Cancer diagnosis via linear 
      programming",\cr SIAM News, Volume 23, Number 5, September 1990, pp 1 & 18.

William H. Wolberg and O.L. Mangasarian: "Multisurface method of 
      pattern separation for medical diagnosis applied to breast cytology", \cr
      Proceedings of the National Academy of Sciences, U.S.A., Volume 87, 
      December 1990, pp 9193-9196.


K. P. Bennett & O. L. Mangasarian: "Robust linear programming 
      discrimination of two linearly inseparable sets",\cr Optimization Methods
      and Software 1, 1992, 23-34 (Gordon & Breach Science Publishers).

}
\examples{
data(wdbc)
str(wdbc)
}
\keyword{datasets}
