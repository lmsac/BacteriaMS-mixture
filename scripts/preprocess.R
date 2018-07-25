#' setwd(...)

local({
  mz.lower = 4000
  mz.upper = 12000
  window = 0.015
  offset = 0.5
  
  dir.create('raw')
  # dir.create('normalized')
  dir.create('peaklists')
  
  lapply(list.files(pattern = '.txt'), function(file) {
    raw = read.table(file)
    
    raw = crop.mz(raw, mz.lower = mz.lower, mz.upper = mz.upper)
    
    baseline = get.baseline(raw, window = window, offset = offset)
    subbase = subtract.baseline(raw, baseline)
    normalized = normalize.intensity(subbase)
    
    peaklist = find.peaks(normalized, window = window, offset = offset)
    
    # write.table(
    #   normalized,
    #   file = paste0('normalized/normalized ', file),
    #   col.names = F,
    #   row.names = F,
    #   quote = F
    # )
    
    write.table(
      peaklist,
      file = paste0('peaklists/peaklist ', file),
      col.names = F,
      row.names = F,
      quote = F
    )
    
    file.rename(file, paste0('raw/', file))
    file
  })
})