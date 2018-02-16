#' 
#' Convert clusters of anchorages into polygons
#'
#' @field allPolygons list List of all polygons
#' @field earthGeo ANY Object of class earthGeoClass
#' @field clusters ANY Object of class ClustersClass
#' @field centroids data.frame Longitude and latitude of polygon centroids
#' @field verbose Output progress to console
#' 
polygonsClass <- setRefClass(
  Class = "polygonsClass",
  fields = list(
    allPolygons = 'list',
    earthGeo = 'ANY',
    clusters = 'ANY',
    centroids = 'data.frame',
    verbose = 'logical'
  ),

  methods = list(

    initialize = function(clusters, autobuild = TRUE, verbose = TRUE) {
      .self$earthGeo = clusters$earthGeo
      .self$clusters = clusters
      .self$verbose = verbose
      if (autobuild) {
        .self$polyAll()
        .self$getCentroids()
      }
    },

    #
    # Algorithm that constructs polygons (outer loop)
    #
    polyAll = function() {
      .self$allPolygons = lapply(1:.self$clusters$nClust, FUN = function(i) {
        if (i %% 1000 == 0 && .self$verbose) print(paste0(i, " out of ", .self$clusters$nClust, " done."))
        .self$polyOne(i)
      })
      .self$allPolygons[sapply(.self$allPolygons, is.null)] <- NULL
    },

    #
    # Algorithm that constructs polygons (inner loop, one cluster)
    #
    polyOne = function(i) {
      crit = .self$clusters$groupID@x == i
      cc = .self$clusters$cIdx[crit]
      rr = .self$clusters$rIdx[crit]
      if (length(cc) > 2) {
        lon = mapply(cc, FUN = function(ccc) .self$earthGeo$lon[ccc] +
                       rnorm(1, mean = 0, sd = 0.0001))
        lat = mapply(rr, FUN = function(rrr) .self$earthGeo$lat[rrr] +
                       rnorm(1, mean = 0, sd = 0.0001))
        h = convex.hull(tri.mesh(lon, lat))
        data.frame(x = c(h$x, h$x[1]), y = c(h$y, h$y[1]))
      } else {
        NULL #excluding clusters with only 1 or 2 anchorages
      }
    },

    #
    # Compute the coordinates of the polygon centroids (outer loop)
    #
    getCentroids = function() {
      tmp = t(mapply(.self$allPolygons, FUN = .self$getOneCentroid))
      .self$centroids = data.frame(lon = tmp[,1], lat = tmp[,2])
    },

    #
    # Compute the coordinates of the polygon centroids (inner loop, one centroid)
    #
    getOneCentroid = function(myPoly) {
      c(mean(myPoly$x[2:length(myPoly$x)]), mean(myPoly$y[2:length(myPoly$y)]))
    },

    #
    # Write centroids to file
    #
    exportCentroids = function(outfile) {
      write.table(.self$centroids, file = outfile, row.names = FALSE, col.names = TRUE, sep = ",", quote = FALSE)
    }
  )
)
