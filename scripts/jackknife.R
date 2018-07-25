tolerance = 2000
mix.number = 6
jackknife.ratio = 1/3
jackknife.number = 100

setwd(REFERENCE_DIR)
reference.spectra = local({
  reference.spectra = lapply(list.files(pattern = '.txt', full.names = T), function(file) {
    normalize.intensity(read.table(file))
  })
  reference.names = sapply(list.files(pattern = '.txt'), function(file) {
    strsplit(file, '\\.')[[1]][1]
  })
  names(reference.spectra) = reference.names
  reference.spectra
})

sample.table = c(
  Ko = 'Klebsiella oxytoca ATCC 13182',
  El = 'Enterobacter cloacae ATCC 23373',
  Ec = 'Escherichia coli ATCC 25922',
  Kp = 'Klebsiella pneumoniae ATCC 700603',
  Pa = 'Pseudomonas aeruginosa ATCC 27853',
  Sa = 'Staphylococcus aureus ATCC 25923',
  Vp = 'Vibrio parahaemolyticus ATCC 17802',
  Va = 'Vibrio alginolyticus ATCC 17749',
  Vm = 'Vibrio mimicus ATCC 33653',
  Vf = 'Vibrio fluvialis ATCC 33809'
)

setwd(SAMPLE_DIR)
sample.strains = lapply(list.dirs(full.names = F, recursive = F), function(x) {
  lapply(list.dirs(x, full.names = F, recursive = F), function(y) {
    label = strsplit(gsub('([A-Z])', '-\\1', y), '-')[[1]][-1]
    indexs.1 = grep('1$', label)
    if (length(indexs.1) > 0)
      label = c(label[-indexs.1], label[indexs.1])
    as.character(unlist(sample.table[sub('1', '', label)]))
  })
})

results = lapply(list.dirs(full.names = F, recursive = F), function(x) {
  lapply(list.dirs(x, full.names = F, recursive = F), function(y) {
    setwd(paste0(x, '/', y))
    
    sample.peaklists = lapply(list.files(pattern = '.txt', full.names = T), read.table)
    
    combined.peaklist = combine.peaklists(sample.peaklists, tolerance)
    
    result = search.datebase.mixture.jackknife(combined.peaklist[, 1:2], reference.spectra, tolerance, mix.number = mix.number, jackknife.ratio = jackknife.ratio, jackknife.number = jackknife.number)
    
    setwd('../..')
    message(paste0(x, '/', y))
    result
  })
})

confidence.scores = lapply(results, function(r) {
  lapply(r, function(x) {
    get.confidence.score(x)
  })
})

local({
  dirs = list.dirs(full.names = F, recursive = F)
  lapply(1:length(sample.strains), function(i) {
    result.table = cbind(do.call(rbind, sample.strains[[i]]), 
                         do.call(rbind, lapply(confidence.scores[[i]], function(x) {
                           c(names(x[1:mix.number]), as.numeric(x[1:mix.number]))
                         })))
    
    colnames(result.table) = c(
      paste0('sample', 1:length(sample.strains[[i]][[1]])),
      paste0('reference', 1:mix.number),
      paste0('confidence', 1:mix.number)
    )
    write.csv(result.table, paste0(dirs[i], '/results.csv'), row.names = F)
  })
})
