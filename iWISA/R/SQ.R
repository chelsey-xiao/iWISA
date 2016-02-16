#'Estimation of Solar Quiet (SQ) Variation
#'
#'Function for calculating estimation of a non-constant solar quiet daily variation by using multiresolution analysis (MRA) and functional principal component analysis. The procedure removes the global signature of a ring current to eliminate the storm effect and then apply a wavelet filter to the data and use just the levels that contain the periodic component. In order to get better Sq estimate it is necessary to use records from at least two stations.
#'
#'@param data records of magnetic storm activitys from stations which are not included stations used in \code{"SAIndex"} to estimate index of SI.
#'@param si.v estimation of WISA.
#'@param wf Filter applied to each level j to remove the components which do not reflect storm activities. LA(8) filter is the default one.
#'@param n.levels  the maximum level J is typically chosen to be smaller than its largest admissible value log2(N), where N is the number of observations.
#'@param boundary boundary MODWT coefficients. The default set for boundary is reflection.
#'@param quantile the empirical probability, i.e., the proportion of the MODWT coefficients which exceed level-independent thresholding. Default value is 0.98.
#'
#'@details  \code{SQ} is a function for estimating solar quiet daily variation. Stations used in estimation of sq need be different from stations used in estimating SI. The resulting daily variation is non-constant, and its day-to-day variability is quantified by functional principal component scores. It uses \code{mra} to obtain MRA from the data after storm activity was removed. We compute the first PC rather than the first
#'PC of the raw magnetometer data with some seasonal adjustments which do not remove the storm activity from the Sq.
#'@return SQ the estimation of a non-constant solar quiet daily variation.
#'@references Maslova, P. Kokoszka, J. Sojka, L. Zhu (2010), Estimation of Sq variation by means of multiresolution and principal component analyses.
#'@seealso  \code{\link{SAIndex}}
#'@examples
#'## example data
#'data(record)
#'coord=matrix(c(124.43, 19.23, 53.77,140.18,68.68,202.00,71.89,293.85),nrow=2,ncol=4)
#'
#'## Estimation of SI index
#'index.sample<- SAIndex(record, coord, wf="la8")
#'si.v<-index.sample$SI
#'
#'## example data of stations which are different from the ones that are used to estimate SI
#'## index.
#'
#'## estimation of sq
#'SQ (datasq, si.v=si.v, wf = "la8", n.levels = 10, boundary = "reflection", quantile=0.98)
#'@import waveslim
#'@export



SQ<-
function(data, si.v=si.v, wf = "la8", n.levels = 10, boundary = "reflection",
            quantile=0.98){
J0<-7
J1<-10
if(is.vector(data)==T) {n.station<-1
                        N<-length(data)} else {n.station<-dim(data)[2] #each column of the data matrix are records for each individual station
N<-dim(data)[1]}         # number of observations for each station

data.ds<-matrix(NA, ncol=n.station, nrow = N)
# remove the storm activity from all stations
for(i in 1:n.station){data.ds[,i]<-data[,i] - si.v}
mra.sq <- matrix(data = 0, ncol = n.station, nrow = N)
for(i in 1:n.station){
# obtain MRA from the data after storm activity was removed
    data.mra<-mra(data.ds[,i], wf=wf, J = n.levels, boundary=boundary)
# MRA leveles used to extract periodic component
for(f in (J0+1):J1) mra.sq[,i] <- mra.sq[,i] + data.mra[[f]]
}
# extract periodic component
deco<-pca.SQ.new(mra.sq)
# Sq component
SQ <- deco$sq
return(SQ)
}
