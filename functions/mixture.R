mix.peaklists = function(peaklists, tolerance, coefficients = rep(1, length(peaklists))) {
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
  return(peaks)
}

find.optimal.mixture = function(sample.peaklist, reference.peaklists, tolerance) {
  peaklists = c(reference.peaklists, list(sample.peaklist))
  aligned.peaks = align.peaks.multi(find.common.peaks.multi(peaklists, tolerance))
  
  mzs = aligned.peaks[, 1:length(peaklists)]
  intensities = aligned.peaks[, (length(peaklists) + 1):(length(peaklists) * 2)]

  #' lm.result = lm(
  #'   formula = as.formula(paste0(
  #'     colnames(intensities)[length(colnames(intensities))],
  #'     '~',
  #'     paste(colnames(intensities)[-length(colnames(intensities))], collapse = '+'),
  #'     '+0'
  #'   )),
  #'   data = as.data.frame(intensities.1)
  #' )
  #' coefficients = lm.result$coefficients
  #' coefficients[which(coefficients < 0)] = 0
  
  nnls.result = nnls::nnls(intensities[, -ncol(intensities)], intensities[, ncol(intensities)])
  coefficients = nnls.result$x
  names(coefficients) = colnames(intensities)[-length(colnames(intensities))]
  
  intensity = Reduce('+', lapply(colnames(intensities)[-length(colnames(intensities))], function(col) {
    intensities[, col] * coefficients[col]
  }))
  mz = apply(cbind(mzs[, -ncol(mzs)], intensities[, -ncol(intensities)]), 1, function(row) {
    mz = row[1:(length(row) / 2)]
    inten = row[(length(row) / 2 + 1):length(row)]
    sum(mz * inten * coefficients, na.rm = TRUE) / sum(inten * coefficients)
  })
  
  indexs = which(!is.na(mz) & intensity != 0)
  peaks = cbind(mz[indexs], intensity[indexs])
  colnames(peaks) = c('mz', 'intensity')
  return(list(peaks = peaks, coefficients = coefficients))
}