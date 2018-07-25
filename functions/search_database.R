search.datebase.mixture = function(sample.peaklist, reference.spectra, tolerance, mix.number = 2, method = 'cosine') {
  cp.ratio = do.call(rbind, lapply(reference.spectra, function(ref) {
    score(sample.peaklist, ref, tolerance, method = 'jaccard')
  }))
  
  orders = local({
    orders = order(cp.ratio, decreasing = TRUE)
    strains = names(reference.spectra)[orders]
    species = sapply(strsplit(strains, ' '), function(x) paste(x[1], x[2]))
    orders[!duplicated(species)]
  })
  
  n = max(5, mix.number + 1)
  indexs = t(apply(gtools::combinations(n, mix.number), 1, function(indexs) {
    orders[indexs]
  }))
  
  mixtures = lapply(1:nrow(indexs), function(i) {
    mix = find.optimal.mixture(sample.peaklist, reference.spectra[indexs[i, ]], tolerance)
    coe = mix$coefficients
    coe.order = order(coe, decreasing = TRUE)
    list(reference = indexs[i, coe.order],
         coefficients = coe[coe.order],
         peaks = mix$peaks)
  })
  
  scores = do.call(rbind, lapply(mixtures, function(mixture) {
    score(mixture$peaks, sample.peaklist, tolerance, method = method)
  }))
  
  result = lapply(method, function(m) {
    ref.order = order(scores[, m], decreasing = TRUE)
    reference = lapply(mixtures[ref.order], function(mixture) names(reference.spectra)[mixture$reference])
    coefficients = lapply(mixtures[ref.order], function(mixture) mixture$coefficients)
    r = cbind(
      reference = do.call(rbind, reference),
      coefficients = do.call(rbind, coefficients),
      score = scores[, m][ref.order]
    )
    colnames(r) = c(paste0('reference', 1:mix.number), paste0('coefficient', 1:mix.number), 'score')
    r
  })
  names(result) = method
  return(result)
}