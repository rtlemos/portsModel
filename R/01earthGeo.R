#'
#' Calculate how many cells in longitude should be used, so that the search square is roughly 0.1 * 0.1 km^2
#'
#' @field earthCircumfKm Circumference of the Earth, in km
#' @field length1DegLatKm Length of one degree of latitude, in km
#' @field resKm Desired spatial resolution, in km
#' @field resDeg Meridional spatial resolution, in degrees
#' @field nlat Number of latitude cells that go from 90N to 90S
#' @field nlon Number of longitude cells that go from 180W to 180E
#' @field lat Array of latitudes (in decreasing order)
#' @field lon Array of longitudes (in increasing order)
#'
earthGeoClass <- setRefClass(
  Class = 'earthGeoClass',

  fields = list(earthCircumfKm = 'numeric',
                length1DegLatKm = 'numeric',
                resKm = 'numeric',
                resDeg = 'numeric',
                nlat = 'numeric',
                nlon = 'numeric',
                lat = 'numeric',
                lon = 'numeric'),

  methods = list(

    initialize = function() {
      .self$earthCircumfKm = 40075
      .self$length1DegLatKm = .self$earthCircumfKm / 360
      .self$resKm = 0.1
      .self$resDeg = .self$resKm / .self$length1DegLatKm
      .self$nlat = round(.self$earthCircumfKm / (.self$resKm * 2))
      .self$nlon = round(.self$earthCircumfKm / .self$resKm)
      .self$lat = seq(90 - .self$resDeg / 2, -90 + .self$resDeg / 2, by = -.self$resDeg)
      .self$lon = seq(-180 + .self$resDeg / 2, 180 - .self$resDeg / 2, by = .self$resDeg)
    },

    #
    # Length of one degree of longitude, in km, as a function of latitude
    #
    length1DegLonKm = function(latDeg) {
      .self$earthCircumfKm * cos(latDeg * pi / 180) / 360
    },

    #
    # Number of cells that must be traversed, zonally, in order to cover 4km. Capped at 200 by default.
    #
    deltaLon = function(latDeg, maxCells = 200) {
      min(maxCells, round(4 / ((.self$lon[2] - .self$lon[1]) * .self$length1DegLonKm(latDeg))))
    }
  )
)
