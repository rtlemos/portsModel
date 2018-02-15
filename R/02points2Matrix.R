points2MatrixClass <- setRefClass(

  Class = 'points2MatrixClass',

  fields = list(id = 'character',
                country = 'character',
                valid = 'logical',
                mat = 'Matrix',
                nlon = "numeric",
                nlat = "numeric",
                idx = 'numeric',
                delta = 'numeric'),

  methods = list(

    initialize = function(earthGeo, points, latName = 'lat', lonName = 'lon', id = 's2id', country = NULL) {

      .self$nlon = earthGeo$nlon
      .self$nlat = earthGeo$nlat
      
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
        if (i %% 10000 == 0) print(paste0(i, " out of ", n, " done."))
        mylat = lat[i]
        mylon = lon[i]
        crd = .self$getLatLonIdx(mylat, mylon)
        .self$valid[i] = (.self$mat[crd[1], crd[2]] == 0)
        if (.self$valid[i]) {
          .self$mat[crd[1], crd[2]] = .self$idx
          .self$idx = .self$idx + 1
        }
      })
      .self$id = .self$id[.self$valid]
    },

    getLatLonIdx = function(mylat, mylon) {
      rr = min(.self$nlat, max(1, round((90 - mylat) / .self$delta)))
      cc = min(.self$nlon, max(1, round((mylon + 180) / .self$delta)))
      c(rr, cc)
    }
  )

)

#eg = earthGeoClass()
#an = anchorages2MatrixClass("~/Documents/fishackathon/unnamed_anchorages_csv/unnamed_anchorages_20171120.csv", eg$lon, eg$lat)
