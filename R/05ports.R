portsClass <- setRefClass(
  
  Class = "portsClass",
  
  fields = list(
    lon = 'numeric',
    lat = 'numeric',
    polygons = 'ANY',
    clusters = 'ANY',
    anchorages = 'ANY',
    earthGeo = 'ANY',
    wpi = 'ANY',
    cities1000 = 'ANY',
    types = 'list',
    uniqueTypes = 'character',
    country = 'character',
    tbl = 'matrix',
    nearestWPI = 'data.frame',
    nearestCity = 'data.frame'
  ),
  
  methods = list(
    initialize = function(polygons, wpi, cities1000, anchorageTypes, autobuild = FALSE) {
      
      .self$polygons = polygons
      .self$clusters = polygons$clusters
      .self$earthGeo = .self$clusters$earthGeo
      .self$anchorages = .self$clusters$anchorages
      .self$wpi = wpi
      .self$cities1000 = cities1000
        
      if (autobuild) {
        .self$nearestWPI = .self$getAllNearest(wpi)
        .self$nearestCity = .self$getAllNearest(cities1000)
        .self$country = .self$getCountry()
        
        .self$uniqueTypes = unique(anchorageTypes$gear_type)
        .self$types = lapply(.self$uniqueTypes, FUN = function(myType) anchorageTypes[anchorageTypes$gear_type == myType, ])
        .self$setTables()
      }
    },
    
    getCountry = function() {
      .self$country = mapply(.self$nearestCity[,1], FUN = function(id) .self$cities1000$country[id])
    },
    
    getAllNearest = function(sparseObj) {
      n = nrow(.self$polygons$centroids)
      raw = t(mapply(1:n, FUN = function(i) {
        if (i %% 1000 == 0) print(paste0(i, " out of ", n, " done."))
        .self$getOneNearest(sparseObj, i)
      }))
      data.frame(index = raw[,1], distance = raw[,2])
    },
    
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
    
    columnIndexes = function(v) {
      unlist(lapply(1:(length(v) - 1), FUN = function(i) rep(i, v[i + 1] - v[i])))
    },
    
    getDistanceKm = function(lat1, lon1, lat2, lon2) {
      earthGeo$earthCircumfKm / pi * asin(sqrt(sin((lat1 - lat2) * pi / 360) ^ 2 + cos(lat1 * pi / 180) * cos(lat2 * pi / 180) * sin((lon1 - lon2) * pi / 360) ^ 2))
    },
    
    setTables = function() {
      .self$tbl = t(mapply(1:.self$clusters$id, FUN = function(i) {
        if (i %% 1000 == 0) print(paste0(i, " out of ", .self$clusters$id, " done."))
        setOneTable(i)
      }))
      colnames(.self$tbl) = .self$uniqueTypes
    },
    
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
