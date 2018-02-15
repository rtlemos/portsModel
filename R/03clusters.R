clustersClass <- setRefClass(

  Class = 'clustersClass',

  fields = list(anchorages = 'ANY',
                deltaR = 'numeric',
                nn = 'numeric',
                groupID = 'Matrix',
                rIdx = 'numeric',
                cIdx = 'numeric',
                iIdx = 'numeric',
                id = 'numeric',
                earthGeo = 'ANY'),

  methods = list(
    initialize = function(earthGeo, anchorages, autobuild = TRUE) {
      .self$earthGeo = earthGeo
      .self$anchorages = anchorages
      .self$deltaR = earthGeo$deltaLon(0)
      .self$nn = length(anchorages$mat@x)
      .self$groupID = anchorages$mat
      .self$groupID@x[] = -1
      .self$rIdx = anchorages$mat@i + 1
      .self$cIdx = columnIndexes(anchorages$mat@p)
      .self$id = 0
      
      if (autobuild) .self$buildGroups()
    },

    buildGroups = function() {
      for(i in 1:.self$nn) {
        if (i %% 10000 == 0) print(paste0(i, " out of ", .self$nn, " done."))
        buildOneGroup(i)
      }
    },

    buildOneGroup = function(i) {
      if (.self$groupID[.self$rIdx[i], .self$cIdx[i]] == -1) {
        .self$id = .self$id + 1
        r = .self$rIdx[i]
        c = .self$cIdx[i]
        recursion(r, c, getDeltaR(r), getDeltaC(r))
      }
    },
    
    getSubset = function(rb, cb) {
      u = .self$groupID[rb[1]:rb[2], cb[1]:cb[2]]
      crit = (u@x == -1)
      rr = rb[1] + u@i
      cc = cb[1] + columnIndexes(u@p) - 1
      cbind(rr[crit], cc[crit])
    },

    recursion = function(r, c, deltaR, deltaC) {
      rb = c(max(1, r - deltaR), min(.self$earthGeo$nlat, r + deltaR))
      cb = c(max(1, c - deltaC), min(.self$earthGeo$nlon, c + deltaC))
      subset = getSubset(rb, cb)
      ns = nrow(subset)
      if (ns > 0) {
        for (k in 1:ns) .self$groupID[subset[k, 1], subset[k, 2]] = .self$id
        for (k in 1:ns) .self$recursion(subset[k, 1], subset[k, 2], deltaR, deltaC)
      }
    },

    columnIndexes = function(v) {
      n = length(v)
      d = v[2:n] - v[1:(n - 1)]
      unlist(lapply(1:(n - 1), FUN = function(i) rep(i, d[i])))
    },

    getDeltaC = function(r) {
      .self$earthGeo$deltaLon(earthGeo$lat[r])
    },

    getDeltaR = function(r) {
      .self$deltaR
    }
  )
)



