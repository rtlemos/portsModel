% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/05ports.R
\docType{class}
\name{portsClass-class}
\alias{portsClass-class}
\alias{portsClass}
\title{Ports Class. Create an object that characterizes world ports, defined as polygons enveloping AIS-based anchorages}
\description{
Ports Class. Create an object that characterizes world ports, defined as polygons enveloping AIS-based anchorages
}
\section{Fields}{

\describe{
\item{\code{polygons}}{Object of class polygonsClass}

\item{\code{clusters}}{Object of class clustersClass}

\item{\code{anchorages}}{Object of class points2MatrixClass}

\item{\code{earthGeo}}{Object of class earthGeoClass}

\item{\code{wpi}}{World Ports Index object, of class points2MatrixClass}

\item{\code{cities1000}}{Geonames cities1000 object, of class points2MatrixClass}

\item{\code{uniqueTypes}}{Array of vessel types, found in anchorageTypes object}

\item{\code{types}}{Decomposition of anchorageTypes object into a list; each list element contains only records of a given type of vessel}

\item{\code{country}}{ISO-3166 2-letter country code}

\item{\code{counts}}{N*M matrix of counts, where N = number of ports, M = number of types of vessels}

\item{\code{nearestWPI}}{data.frame with N rows and 2 columns; 1st - index of nearest WPI port; 2nd - distance between polygon centroid and nearest WPI port}

\item{\code{nearestCity}}{data.frame with N rows and 2 columns; 1st - index of nearest city; 2nd - distance between polygon centroid and nearest city}

\item{\code{verbose}}{Output progress to console}
}}


