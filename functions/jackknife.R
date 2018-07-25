jackknife.peaklist = function(peaklist, ratio = 1/3, n = 100) {
  peak.number = nrow(peaklist)
  lapply(1:n, function(i) {
    indexs = sample.int(peak.number, peak.number * (1 - ratio))
    if (length(indexs) == 1)
      rbind(peaklist[indexs, ])
    else 
      peaklist[indexs, ]
  })
}

search.datebase.mixture.jackknife = function(sample.peaklist, reference.spectra, tolerance, mix.number = 2, method = 'cosine', jackknife.ratio = 1/3, jackknife.number = 100) {
  peaklists = jackknife.peaklist(sample.peaklist, ratio = jackknife.ratio, n = jackknife.number)
  do.call(rbind, lapply(peaklists, function(peaklist) {
    search.datebase.mixture(peaklist, reference.spectra, tolerance, mix.number, method)[[method]][1, ]
  }))
}

get.confidence.score = function(result) {
  jackknife.number = nrow(result)
  mix.number = length(grep('reference[0-9]+', colnames(result)))
  strains = unlist(result[, paste0('reference', 1:mix.number)])
  species = sapply(strsplit(strains, ' '), function(x) paste(x[1], x[2]))
  count = table(species)
  count = count[order(count, decreasing = TRUE)]
  count / jackknife.number
}
