#'
#' Ports Class. Create an object that characterizes world ports, defined as polygons enveloping AIS-based anchorages
#' 
#' @field polygons Object of class polygonsClass
#' @field clusters Object of class clustersClass
#' @field anchorages Object of class points2MatrixClass
#' @field earthGeo Object of class earthGeoClass
#' @field wpi World Ports Index object, of class points2MatrixClass
#' @field cities1000 Geonames cities1000 object, of class points2MatrixClass
#' @field uniqueTypes Array of vessel types, found in anchorageTypes object
#' @field types Decomposition of anchorageTypes object into a list; each list element contains only records of a given type of vessel
#' @field country ISO-3166 2-letter country code
#' @field counts N*M matrix of counts, where N = number of ports, M = number of types of vessels
#' @field nearestWPI data.frame with N rows and 2 columns; 1st - index of nearest WPI port; 2nd - distance between polygon centroid and nearest WPI port
#' @field nearestCity data.frame with N rows and 2 columns; 1st - index of nearest city; 2nd - distance between polygon centroid and nearest city
#' @field verbose Output progress to console
#'
portsClass <- setRefClass(
  
  Class = "portsClass",
  
  fields = list(
    polygons = 'ANY',
    clusters = 'ANY',
    anchorages = 'ANY',
    earthGeo = 'ANY',
    wpi = 'ANY',
    cities1000 = 'ANY',
    uniqueTypes = 'character',
    types = 'list',
    country = 'character',
    counts = 'matrix',
    nearestWPI = 'data.frame',
    nearestCity = 'data.frame',
    verbose = 'logical'
  ),
  
  methods = list(
    
    initialize = function(polygons, wpi, cities1000, anchorageTypes, autobuild = TRUE, verbose = TRUE) {
      .self$polygons = polygons
      .self$clusters = polygons$clusters
      .self$earthGeo = .self$clusters$earthGeo
      .self$anchorages = .self$clusters$anchorages
      .self$wpi = wpi # to read original WPI shapefile, use readOGR(dsn = "/.../WPI_Shapefile/WPI.shp")
      .self$cities1000 = cities1000
      .self$verbose = verbose
        
      if (autobuild) {
        .self$nearestWPI = .self$getAllNearest(wpi, 1)
        .self$nearestCity = .self$getAllNearest(cities1000, 2)
        .self$country = .self$getCountry()
        
        .self$uniqueTypes = unique(anchorageTypes$gear_type)
        .self$types = lapply(.self$uniqueTypes, FUN = function(myType) anchorageTypes[anchorageTypes$gear_type == myType, ])
        .self$setTables()
      }
    },
    
    #
    # Retrieve port data
    #
    getData = function() {
      lapply(1:length(polygons$allPolygons), FUN = function(i) {
        myPoly = .self$polygons$allPolygons[[i]]
        myCentroid = .self$polygons$centroids[i,]
        myCounts = .self$counts[i,]
        myWPI = as.numeric(.self$nearestWPI[i, ])
        myCity = as.numeric(.self$nearestCity[i, ])
        list(lon = myPoly$x, lat = myPoly$y, centroidLon = myCentroid$lon, centroidLat = myCentroid$lat, counts = myCounts, 
             wpi = .self$wpi[myWPI[1]], wpiDistance = myWPI[2], city = .self$cities100[myCity[1]], cityDistance = myCity[2])
      })
    },
    
    #
    # Populate the country field
    # 
    getCountry = function() {
      mapply(.self$nearestCity[,1], FUN = function(id) .self$cities1000$country[id])
    },
    
    #
    # Given an object sparseObj of class points2MatrixClass, obtain, for each port, the index of the nearest point in sparseObj and the distance
    #
    getAllNearest = function(sparseObj, part) {
      n = nrow(.self$polygons$centroids)
      if (.self$verbose) print(paste0('Mapping, part ', part))
      raw = t(mapply(1:n, FUN = function(i) {
        if (i %% 1000 == 0 && .self$verbose) print(paste0(i, " out of ", n, " done."))
        .self$getOneNearest(sparseObj, i)
      }))
      data.frame(index = raw[,1], distance = raw[,2])
    },
    
    #
    # Inner loop of function getAllNearest (computations for a single port)
    # 
    getOneNearest = function(sparseObj, i) {
      crds = as.numeric(.self$polygons$centroids[i, ])
      rc = sparseObj$getLatLonIdx(mylat = crds[2], mylon = crds[1])
      searching = TRUE
      range = 10
      while (searching) {
        rb = c(max(1, rc[1] - range), min(.self$earthGeo$nlat, rc[1] + range))
        cb = c(max(1, rc[2] - range), min(.self$earthGeo$nlon, rc[2] + range))
        u = sparseObj$mat[rb[1]:rb[2], cb[1]:cb[2]]
        if (length(u@x) > 0) {
          rIdx = rb[1] + u@i
          cIdx = cb[1] + .self$columnIndexes(u@p) - 1
          d = .self$getDistanceKm(lat1 = crds[2], lon1 = crds[1], lat2 = .self$earthGeo$lat[rIdx], lon2 = .self$earthGeo$lon[cIdx])
          j = which.min(d)
          nearest = c(sparseObj$mat[rIdx[j], cIdx[j]], d[j])
          searching = FALSE
        } else {
          range = range * 10
        }
      }
      nearest
    },
    
    #
    # Given the p array of a sparse matrix, get the column indexes
    #
    columnIndexes = function(v) {
      unlist(lapply(1:(length(v) - 1), FUN = function(i) rep(i, v[i + 1] - v[i])))
    },
    
    #
    # Haversine distance calculator
    #
    getDistanceKm = function(lat1, lon1, lat2, lon2) {
      earthGeo$earthCircumfKm / pi * asin(sqrt(sin((lat1 - lat2) * pi / 360) ^ 2 + cos(lat1 * pi / 180) * cos(lat2 * pi / 180) * sin((lon1 - lon2) * pi / 360) ^ 2))
    },
    
    #
    # Populate the counts field
    #
    setTables = function() {
      if (.self$verbose) print('Constructing vessel type table')
      .self$counts = t(mapply(1:.self$clusters$nClust, FUN = function(i) {
        if (i %% 1000 == 0 && .self$verbose) print(paste0(i, " out of ", .self$clusters$nClust, " done."))
        setOneTable(i)
      }))
      colnames(.self$counts) = .self$uniqueTypes
    },
    
    #
    # Inner loop of function setTables (computations for a single port)
    # 
    setOneTable = function(i) {
      idx = which(.self$clusters$groupID@x == i)
      myS2id = .self$anchorages$id[idx]
      cnt = mapply(.self$types, FUN = function(mygear) {
        sum(unlist(mapply(myS2id, FUN = function(s) mygear$counts[mygear$s2id == s])))
      })
      cnt
    }
  )
)
