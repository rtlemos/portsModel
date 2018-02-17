#'
#' Convert a dataset of points into a sparse matrix, for faster downstream computations
#'
#' @field id character ID of point
#' @field country character Country of point, if available
#' @field valid logical Is point valid or redundant?
#' @field mat Matrix Sparse matrix with 100 m resolution at the equator
#' @field nlat Number of latitude cells that go from 90N to 90S
#' @field nlon Number of longitude cells that go from 180W to 180E
#' @field idx numeric Total number of valid IDs
#' @field delta numeric Meridional spatial resolution, in degrees
#' @field verbose Output progress to console
#'
points2MatrixClass <- setRefClass(

  Class = 'points2MatrixClass',

  fields = list(id = 'character',
                country = 'character',
                valid = 'logical',
                mat = 'Matrix',
                nlon = "numeric",
                nlat = "numeric",
                idx = 'numeric',
                delta = 'numeric',
                verbose = 'logical'),

  methods = list(

    initialize = function(earthGeo, points, latName = 'lat', lonName = 'lon', id = 's2id', country = NULL, verbose = TRUE) {

      .self$nlon = earthGeo$nlon
      .self$nlat = earthGeo$nlat
      .self$verbose = verbose
      
      .self$id = points[[id]]
      if (!is.null(country)) {
        .self$country = points[[country]]
      }
      .self$valid = rep(TRUE, length(.self$id))
      .self$delta = earthGeo$lat[1] - earthGeo$lat[2]
      .self$mat = Matrix(nrow = earthGeo$nlat, ncol = earthGeo$nlon, 0)
      n = nrow(points)
      lat = points[[latName]]
      lon = points[[lonName]]
      .self$idx = 1
      
      mapply(1:n, FUN = function(i) {
        if (i %% 10000 == 0 && .self$verbose) print(paste0(i, " out of ", n, " done."))
        mylat = lat[i]
        mylon = lon[i]
        crd = .self$getLatLonIdx(mylat, mylon)
        .self$valid[i] = (.self$mat[crd[1], crd[2]] == 0)
        if (.self$valid[i]) {
          .self$mat[crd[1], crd[2]] = i
          .self$idx = .self$idx + 1
        }
      })
      .self$id = .self$id[.self$valid]
    },

    #
    # Obtain row and column indices, for a given (lat, lon) pair
    #
    getLatLonIdx = function(mylat, mylon) {
      rr = min(.self$nlat, max(1, round((90 - mylat) / .self$delta)))
      cc = min(.self$nlon, max(1, round((mylon + 180) / .self$delta)))
      c(rr, cc)
    }
  )

)