find.common.peaks = function(peaks1, peaks2, tolerance) {
  if(class(peaks1) == 'data.frame')
    peaks1 = as.matrix(peaks1)
  if(class(peaks2) == 'data.frame')
    peaks2 = as.matrix(peaks2)
  
  indexs.1 = order(peaks1[, 2], decreasing = TRUE)
  
  is.used = rep(FALSE, nrow(peaks2))
  indexs.2 = sapply(indexs.1, function(i) {
    mz1 = peaks1[i, 1]
    candidates = which(!is.used & peaks2[, 1] >= mz1 * (1 - tolerance * 1e-6) & peaks2[, 1] <= mz1 * (1 + tolerance * 1e-6))
    j = candidates[which.min(abs(peaks2[candidates, 1] - mz1))]
    if(length(j) == 0)
      return(NA)
    else {
      is.used[j] <<- TRUE
      return(j)
    }
  })
  
  common.indexs = which(!is.na(indexs.2))
  if(length(common.indexs) == 0)
    list(
      common.peaks1 = peaks1[c(), ],
      common.peaks2 = peaks2[c(), ],
      peaks1 = peaks1,
      peaks2 = peaks2,
      tolerance = tolerance
    )
  else if(length(common.indexs) == 1)
    list(
      common.peaks1 = matrix(peaks1[indexs.1[common.indexs], ], nrow = 1),
      common.peaks2 = matrix(peaks2[indexs.2[common.indexs], ], nrow = 1),
      peaks1 = peaks1,
      peaks2 = peaks2,
      tolerance = tolerance
    )
  else
    list(
      common.peaks1 = peaks1[indexs.1[common.indexs], ],
      common.peaks2 = peaks2[indexs.2[common.indexs], ],
      peaks1 = peaks1,
      peaks2 = peaks2,
      tolerance = tolerance
    )
}

align.peaks = function(common.peaks) {
  peaks1.common = cbind(common.peaks[[1]][, 1], common.peaks[[1]][, 2])
  peaks2.common = cbind(common.peaks[[2]][, 1], common.peaks[[2]][, 2])
  peaks1 = cbind(common.peaks[[3]][, 1], common.peaks[[3]][, 2])
  peaks2 = cbind(common.peaks[[4]][, 1], common.peaks[[4]][, 2])
  if(nrow(peaks1.common) > 0) {
    get.noncommom.peaks.index = function(peaks, peaks.common) {
      sapply(1:nrow(peaks), function(i) {
        idx = which(peaks.common[, 1] == peaks[i, 1] & peaks.common[, 2] == peaks[i, 2])
        if(length(idx > 0)) {
          peaks.common[idx[1], 1] <<- NA
          FALSE
        }
        else {
          TRUE
        }
      })
    }
    
    mz2.noncommon = peaks2[get.noncommom.peaks.index(peaks2, peaks2.common), 1]
    mz1.noncommon = peaks1[get.noncommom.peaks.index(peaks1, peaks1.common), 1]
    peaks1.new = if(length(mz2.noncommon) > 0) rbind(peaks1, cbind(mz2.noncommon, 0)) else peaks1
    peaks2.new = if(length(mz1.noncommon) > 0) rbind(peaks2, cbind(mz1.noncommon, 0)) else peaks2
  }
  else {
    peaks1.new = rbind(peaks1, cbind(peaks2[, 1], 0))
    peaks2.new = rbind(peaks2, cbind(peaks1[, 1], 0))
  }
  peaks1.new = peaks1.new[order(peaks1.new[, 1]), ]
  peaks2.new = peaks2.new[order(peaks2.new[, 1]), ]
  peaks1.new[which(peaks1.new[, 2] == 0), 1] = NA
  peaks2.new[which(peaks2.new[, 2] == 0), 1] = NA
  colnames(peaks1.new) = c('mz1', 'int1')
  colnames(peaks2.new) = c('mz2', 'int2')
  return(cbind(peaks1.new, peaks2.new))
}

cosine = function(common.peaks) {
  peaks1.common = common.peaks[[1]]
  peaks2.common = common.peaks[[2]]
  peaks1 = common.peaks[[3]]
  peaks2 = common.peaks[[4]]
  peaks1.n = nrow(peaks1)
  peaks2.n = nrow(peaks2)
  peaks1.common.n = nrow(peaks1.common)
  if(peaks1.common.n == 0)
    return(0)
  norm1 = (sum((peaks1[, 2]) ^ 2)) ^ (0.5)
  norm2 = (sum((peaks2[, 2]) ^ 2)) ^ (0.5)
  cos.score = sum(peaks1.common[, 2] * peaks2.common[, 2]) / (norm1 * norm2)
  return(cos.score)
}

jaccard = function(common.peaks) {
  peaks1.common = common.peaks[[1]]
  peaks2.common = common.peaks[[2]]
  peaks1 = common.peaks[[3]]
  peaks2 = common.peaks[[4]]
  peaks1.n = nrow(peaks1)
  peaks2.n = nrow(peaks2)
  peaks1.common.n = nrow(peaks1.common)
  return(peaks1.common.n / (peaks1.n + peaks2.n - peaks1.common.n))
}

score = function(peaks1, peaks2, tolerance, method = c('cosine')) {
  peaks1 = normalize.intensity(peaks1) 
  peaks2 = normalize.intensity(peaks2)
  common.peaks = find.common.peaks(peaks1, peaks2, tolerance)
  scores = sapply(method, function(m) {
    do.call(m, args = list(common.peaks = common.peaks))
  })
  names(scores) = method
  scores
}