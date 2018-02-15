#' 
#' Aggregate anchorages into clusters
#'
#' @field anchorages ANY Object of class poins2Matrix
#' @field deltaR numeric Number of cells that must be traversed to cover 4 km in latitude
#' @field nn numeric Number of valid anchorages
#' @field groupID Matrix Sparse matrix indicating the cluster ID for each anchorage
#' @field rIdx numeric Row indices (in groupID) of valid anchorages
#' @field cIdx numeric Column indices (in groupID) of valid anchorages
#' @field nClust numeric Number of clusters found
#' @field earthGeo ANY Object of class earthGeoClass
#' @field verbose Output progress to console
#'
clustersClass <- setRefClass(

  Class = 'clustersClass',

  fields = list(anchorages = 'ANY',
                deltaR = 'numeric',
                nn = 'numeric',
                groupID = 'Matrix',
                rIdx = 'numeric',
                cIdx = 'numeric',
                nClust = 'numeric',
                earthGeo = 'ANY',
                verbose = 'logical'),

  methods = list(
    initialize = function(earthGeo, anchorages, autobuild = TRUE, verbose = TRUE) {
      .self$earthGeo = earthGeo
      .self$anchorages = anchorages
      .self$deltaR = earthGeo$deltaLon(0)
      .self$nn = length(anchorages$mat@x)
      .self$groupID = anchorages$mat
      .self$groupID@x[] = -1
      .self$rIdx = anchorages$mat@i + 1
      .self$cIdx = columnIndexes(anchorages$mat@p)
      .self$nClust = 0
      .self$verbose = verbose
      
      if (autobuild) .self$buildGroups()
    },

    #
    # Clustering algorithm (outer loop)
    #
    buildGroups = function() {
      for(i in 1:.self$nn) {
        if (i %% 10000 == 0 && .self$verbose) print(paste0(i, " out of ", .self$nn, " done."))
        buildOneGroup(i)
      }
    },

    #
    # Clustering algorithm (inner loop)
    #
    buildOneGroup = function(i) {
      if (.self$groupID[.self$rIdx[i], .self$cIdx[i]] == -1) {
        .self$nClust = .self$nClust + 1
        r = .self$rIdx[i]
        c = .self$cIdx[i]
        recursion(r, c, getDeltaR(r), getDeltaC(r))
      }
    },
    
    #
    # Extract a chunk of the big sparse matrix and check for unassigned anchorages in it; return their row and column indices
    #
    getSubset = function(rb, cb) {
      u = .self$groupID[rb[1]:rb[2], cb[1]:cb[2]]
      crit = (u@x == -1)
      rr = rb[1] + u@i
      cc = cb[1] + columnIndexes(u@p) - 1
      cbind(rr[crit], cc[crit])
    },

    #
    # Recursive clustering algorithm
    #
    recursion = function(r, c, deltaR, deltaC) {
      rb = c(max(1, r - deltaR), min(.self$earthGeo$nlat, r + deltaR))
      cb = c(max(1, c - deltaC), min(.self$earthGeo$nlon, c + deltaC))
      subset = getSubset(rb, cb)
      ns = nrow(subset)
      if (ns > 0) {
        for (k in 1:ns) .self$groupID[subset[k, 1], subset[k, 2]] = .self$nClust
        for (k in 1:ns) .self$recursion(subset[k, 1], subset[k, 2], deltaR, deltaC)
      }
    },

    #
    # Given the p array of a sparse matrix, get the column indexes
    #
    columnIndexes = function(v) {
      n = length(v)
      d = v[2:n] - v[1:(n - 1)]
      unlist(lapply(1:(n - 1), FUN = function(i) rep(i, d[i])))
    },

    #
    # Number of cells that must be traversed zonally to cover 4 km
    #
    getDeltaC = function(r) {
      .self$earthGeo$deltaLon(earthGeo$lat[r])
    },

    #
    # Number of cells that must be traversed meridionally to cover 4 km
    #
    getDeltaR = function(r) {
      .self$deltaR
    }
  )
)



