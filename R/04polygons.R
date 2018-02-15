polygonsClass <- setRefClass(
  Class = "polygonsClass",
  fields = list(
    allPolygons = 'list',
    earthGeo = 'ANY',
    clusters = 'ANY',
    centroids = 'data.frame'
  ),

  methods = list(

    initialize = function(clusters, autobuild = TRUE) {
      .self$earthGeo = clusters$earthGeo
      .self$clusters = clusters
      if (autobuild) {
        .self$polyAll()
        .self$getCentroids()
      }
    },

    polyAll = function() {
      .self$allPolygons = lapply(1:.self$clusters$id, FUN = function(i) {
        if (i %% 1000 == 0) print(paste0(i, " out of ", .self$clusters$id, " done."))
        .self$polyOne(i)
      })
      .self$allPolygons[sapply(.self$allPolygons, is.null)] <- NULL
    },

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
        NULL
        #data.frame(x = .self$clusters$lon[cc], y = .self$clusters$lat[rr])
      }
    },

    getCentroids = function() {
      tmp = t(mapply(.self$allPolygons, FUN = .self$getOneCentroid))
      .self$centroids = data.frame(lon = tmp[,1], lat = tmp[,2])
    },

    getOneCentroid = function(myPoly) {
      c(mean(myPoly$x[2:length(myPoly$x)]), mean(myPoly$y[2:length(myPoly$y)]))
    },

    exportCentroids = function(outfile) {
      write.table(.self$centroids, file = outfile, row.names = FALSE, col.names = TRUE, sep = ",", quote = FALSE)
    }
  )
)
