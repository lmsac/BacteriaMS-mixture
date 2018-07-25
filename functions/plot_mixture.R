plot.compare.mixture.peaklist = function(sample.peaklist, peaklists, tolerance, coefficients = rep(1, length(peaklists)), col = c('cadetblue4', 'springgreen4', 'tomato3', 'orange4'), include.individuals = F, ...) {
  aligned.peaks = align.peaks.multi(find.common.peaks.multi(peaklists, tolerance))
  
  mzs = aligned.peaks[, 1:length(peaklists)]
  intensities = aligned.peaks[, (length(peaklists) + 1):(length(peaklists) * 2)]
  
  intensity = Reduce('+', lapply(1:length(peaklists), function(i) {
    intensities[, i] * coefficients[i]
  }))
  mz = apply(cbind(mzs, intensities), 1, function(row) {
    mz = row[1:(length(row) / 2)]
    inten = row[(length(row) / 2 + 1):length(row)]
    sum(mz * inten * coefficients, na.rm = TRUE) / sum(inten * coefficients)
  })
  
  indexs = which(!is.na(mz) & intensity != 0)
  peaks = cbind(mz[indexs], intensity[indexs])
  colnames(peaks) = c('mz', 'intensity')
  
  points = lapply(1:nrow(peaks), function(i) {
    mz = peaks[i, 1]
    inten = peaks[i, 2]
    x = intensities[indexs[i], ] * coefficients
    x = x / sum(x)
    lower = c(0, cumsum(x)[-length(x)]) * inten
    upper = cumsum(x) * inten
    list(mz = mz,
         lower = lower,
         upper = upper)
  })
  
  cp = find.common.peaks(sample.peaklist, peaks, tolerance)
  
  if(include.individuals) {
    spectra = c(peaklists, list(peaks, sample.peaklist))
    xmin = min(unlist(lapply(spectra, function(spectrum) { min(unlist(spectrum[, 1])) })))
    xmax = max(unlist(lapply(spectra, function(spectrum) { max(unlist(spectrum[, 1])) })))
    ymin = min(unlist(lapply(spectra, function(spectrum) { min(unlist(spectrum[, 2])) })))
    ymax = max(unlist(lapply(spectra, function(spectrum) { max(unlist(spectrum[, 2])) })))
    
    plot(
      c(0, 0),
      xlim = c(xmin, xmax),
      ylim = c(0, (ymax - ymin) * length(spectra)),
      xlab = "m/z",
      ylab = "intensity",
      yaxt="n"
    )
    
    lapply(1:length(peaklists), function(i) {
      x = apply(spectra[[i]], 1, function(point) { point[[1]] })
      y = apply(spectra[[i]], 1, function(point) { point[[2]] - ymin + (ymax - ymin) * (i - 1) })
      color.index = i %% length(col)
      if(color.index == 0) { color.index = length(col) }
      color = col[color.index]
      mz.presence = find.common.peaks(spectra[[i]], sample.peaklist, tolerance)[[1]][, 1]
      lapply(1:length(x), function(j) {
        if (x[j] %in% mz.presence)
          col. = color
        else
          col. = 'gray'
        lines(c(x[j], x[j]), c((ymax - ymin) * (i - 1), y[j]), col = col., ...)
      })
      lines(c(xmin * 0.5, xmax * 2), c((ymax - ymin) * (i - 1), (ymax - ymin) * (i - 1)))
    })
    
    lapply(points, function(p) {
      lapply(1:length(p$lower), function(i) {
        color.index = i %% length(col)
        if (color.index == 0) {
          color.index = length(col)
        }
        color = col[color.index]
        lines(c(p$mz, p$mz), c(p$lower[i] + (ymax - ymin) * length(peaklists), p$upper[i] + (ymax - ymin) * length(peaklists)), col = color, ...)
      })
    })
    lines(c(xmin * 0.5, xmax * 2), c((ymax - ymin) * length(peaklists), (ymax - ymin) * length(peaklists)))
    
    lapply(1:nrow(sample.peaklist), function(i) {
      lines(c(sample.peaklist[i, 1], sample.peaklist[i, 1]), c((ymax - ymin) * (length(peaklists) + 1), (ymax - ymin) * (length(peaklists) + 1) + sample.peaklist[i, 2]), col = 'gray', ...)
    })
    lapply(1:nrow(cp[[1]]), function(i) {
      lines(c(cp[[1]][i, 1], cp[[1]][i, 1]), c((ymax - ymin) * (length(peaklists) + 1), (ymax - ymin) * (length(peaklists) + 1) + cp[[1]][i, 2]), col = col[length(col)], ...)
    })
    lines(c(xmin * 0.5, xmax * 2), c((ymax - ymin) * (length(peaklists) + 1), (ymax - ymin) * (length(peaklists) + 1)))
  }
  else {
    xmin = min(peaks[, 1], sample.peaklist[, 1])
    xmax = max(peaks[, 1], sample.peaklist[, 1])
    ymax = max(peaks[, 2], sample.peaklist[, 2])
    
    plot(
      c(0, 0),
      xlim = c(xmin, xmax),
      ylim = c(-ymax, ymax),
      xlab = "m/z",
      ylab = "intensity"
    )
    
    lapply(points, function(p) {
      lapply(1:length(p$lower), function(i) {
        color.index = i %% length(col)
        if (color.index == 0) {
          color.index = length(col)
        }
        color = col[color.index]
        lines(c(p$mz, p$mz), c(-p$lower[i], -p$upper[i]), col = color, ...)
      })
    })
    
    lapply(1:nrow(sample.peaklist), function(i) {
      lines(c(sample.peaklist[i, 1], sample.peaklist[i, 1]), c(0, sample.peaklist[i, 2]), col = 'gray', ...)
    })
    lapply(1:nrow(cp[[1]]), function(i) {
      lines(c(cp[[1]][i, 1], cp[[1]][i, 1]), c(0, cp[[1]][i, 2]), col = col[length(col)], ...)
    })
    lines(c(xmin * 0.5, xmax * 2), c(0, 0))
  }
  return()
}
