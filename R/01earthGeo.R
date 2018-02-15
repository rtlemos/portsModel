#
# Calculate how many cells in longitude should be used,
# so that the search square is roughly 0.5 * 0.5 km^2
#
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

    length1DegLonKm = function(latDeg) {
      .self$earthCircumfKm * cos(latDeg * pi / 180) / 360
    },

    deltaLon = function(latDeg, maxCells = 200) {
      min(maxCells, round(4 / ((.self$lon[2] - .self$lon[1]) * .self$length1DegLonKm(latDeg))))
    }
  )
)
