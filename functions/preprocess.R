normalize.intensity = function(data, max.intensity = 100) {
  normalized.intensity = data[, 2] / max(data[, 2]) * max.intensity
  normalized.intensity[normalized.intensity < 0] = 0
  new.data = data
  new.data[, 2] = normalized.intensity
  new.data
}

crop.mz = function(data, mz.lower = 4000, mz.upper = 12000) {
  data[data[, 1] >= mz.lower & data[, 1] <= mz.upper, ]
}

get.intensity = function(data, mz) {
  linear.interpolate.y = function(x1, y1, x2, y2, x) {
    if (y1 == y2)
      y1
    else if (x == x1)
      y1
    else if (x == x2)
      y2
    else
      y1 + (y2 - y1) / (x2 - x1) * (x - x1)
  }

  index = which(data[, 1] >= mz)
  index = if(length(index) == 0)
    nrow(data) + 1
  else
    index[1]
  if(index > nrow(data))
    0
  else if(data[index, 1] == mz)
    data[index, 2]
  else if(index == 1)
    0 
  else
    linear.interpolate.y(x1 = data[index - 1, 1], x2 = data[index, 1],
                         y1 = data[index - 1, 2], y2 = data[index, 2], mz)
}

get.baseline = function(data, window = 0.1, offset = 0) {
  window.mids = (function(data, window) {
    window.mids = c()
    window.mid = data[nrow(data), 1]
    while (window.mid > data[1, 1]) {
      window.mids = append(window.mids, window.mid)
      window.mid = window.mid - max(50.0, window.mid * window)
    }
    window.mids = append(window.mids, data[1, 1])
    window.mids = window.mids[order(window.mids)]
  })(data, window)
  
  noise = t(sapply(window.mids, function(window.mid) {
    window.indexs = which(data[, 1] >= window.mid * (1 - window) & data[, 1] <= window.mid * (1 + window))
    if(length(window.indexs) == 1)
      c(level = data[window.indexs, 2], width = 0)
    else
      (function(data) {
        level = median(data[, 2])
        width = median(abs(data[, 2] - level)) * 2
        c(level, width)
      })(data[window.indexs, ])
  }))
  
  gaussian.smooth = function(x, window = 5) {
    if (window > length(x))
      window = length(x)  
    if (window %% 2 != 0)
      window = window - 1
    
    ksize = window + 1
    kernel = sapply(0:(ksize - 1), function(i) {
      r = i - (ksize - 1) / 2
      k = exp(-(r * r / (ksize * ksize / 16)))
    })
    kernel = kernel / sum(kernel)
    
    sapply(0:(length(x) - 1), function(i) {
      sum(sapply(0:window, function(j) {
        index = as.integer(abs(i + j - window / 2))
        if (index >= length(x))
          index = index - 2 * (index - length(x) + 1)
        kernel[j + 1] * x[index + 1]
      }))
    })
  }
  
  level = gaussian.smooth(noise[, 1])
  width = gaussian.smooth(noise[, 2])
  width = abs(width)
  level = level - width * offset
  level[level < 0] = 0
  cbind(mz = window.mids, level = level, width = width)
}

subtract.baseline = function(data, baseline) {
  new.data = data
  if(nrow(baseline) == 1)
    new.data[, 2] = data[, 2] - baseline[1, 2]
  else
    new.data[, 2] = data[, 2] - sapply(data[, 1], function(mz) get.intensity(baseline, mz))
  new.data
}

find.peaks = function(data, baselinefun = get.baseline, 
                      absolute.threshold = 0,
                      relative.threshold = 0,
                      sn.threshold = 3, ...) {
  baseline = baselinefun(data, ...)
  baseline.level = baseline[, c(1, 2)]
  baseline.width = baseline[, c(1, 3)]
  
  local.maxima.indexs = (function(x) {
    pre = x - append(x[-length(x)], Inf, 0)
    nxt = x - append(x[-1], Inf)
    which(pre > 0 & nxt >= 0 | pre >= 0 & nxt > 0)
  })(data[, 2])
  
  get.width = function(data, mz, height) {
    linear.interpolate.x = function(x1, y1, x2, y2, y) {
      if (x1 == x2)
        x1
      if (y == y1)
        x1
      else if (y == y2)
        x2
      else
        x1 + (x2 - x1) / (y2 - y1) * (y - y1)
    }
    
    left.index = which(data[, 1] <= mz & data[, 2] <= height)
    left.index = if(length(left.index) == 0)
      1
    else
      left.index[length(left.index)]
    left.mz = if (data[left.index, 2] == height)
      data[left.index, 1]
    else
      linear.interpolate.x(
        x1 = data[left.index, 1],
        x2 = data[left.index + 1, 1],
        y1 = data[left.index, 2],
        y2 = data[left.index + 1, 2],
        y = height
      )
    
    right.index = which(data[, 1] >= mz & data[, 2] <= height)
    right.index = if(length(right.index) == 0)
      nrow(data)
    else
      right.index[1]
    right.mz = if (data[right.index, 2] == height)
      data[right.index, 1]
    else
      linear.interpolate.x(
        x1 = data[right.index - 1, 1],
        x2 = data[right.index, 1],
        y1 = data[right.index - 1, 2],
        y2 = data[right.index, 2],
        y = height
      )
    c(width = right.mz - left.mz, left = left.mz, right = right.mz)
  }
  
  candidates = (function() {
    candidates = list()
    prev.right = 0
    lapply(local.maxima.indexs, function(i) {
      mz = data[i, 1]
      intensity = data[i, 2]
      base = get.intensity(baseline.level, mz)
      noise = get.intensity(baseline.width, mz)
      sn = (intensity - base) / noise
      if(sn < sn.threshold)
        return(NULL)
      width = get.width(data, mz, (intensity - base) * 0.5 + base)
      peak = c(
        mz = mz,
        intensity = intensity,
        sn = sn,
        fwhm = width[1]
      )
      names(peak) = c('mz', 'intensity', 'sn', 'fwhm')
      if (length(candidates) == 0) {
        candidates <<- append(candidates, list(peak))
        prev.right <<- max(width[3], mz + width[1] / 2)
      }
      else {
        if (min(width[2], mz - width[1] / 2) < prev.right) {
          if (intensity > candidates[[length(candidates)]][2]) {
            candidates[[length(candidates)]] <<- peak
            prev.right <<- max(width[3], mz + width[1] / 2)
          }
        }
        else {
          candidates <<- append(candidates, list(peak))
          prev.right <<- max(width[3], mz + width[1] / 2)
        }
      }
    })
    candidates
  })()
  
  peaks = if(length(candidates) == 0) {
    m = matrix(nrow = 0, ncol = 4)
    colnames(m) = c('mz', 'intensity', 'sn', 'fwhm')
    m
  }
  else {
    intensity.threshold = max(max(sapply(candidates, function(peak) peak[2])) * relative.threshold, absolute.threshold)
    do.call(rbind, candidates[sapply(candidates, function(peak) peak[2] >= intensity.threshold)])
  }
}