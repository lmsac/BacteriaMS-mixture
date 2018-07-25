plot.sensitivity.errorrate = function(x, ...) {
  plot(c(0, 0), xlim = c(0, 1), ylim = c(0, 1), ...)
  lapply(seq(0, 1, 0.2), function(x) {
    lines(c(x, x), c(-0.2, 1.2), col = 'gray', lty = 2)
  })
  lapply(seq(0, 1, 0.2), function(x) {
    lines(c(-0.2, 1.2), c(x, x), col = 'gray', lty = 2)
  })
  lines(x = x[, 1], y = x[, 2], lwd = 2, col = 'springgreen4')
  lines(x = x[, 1], y = 1 - x[, 3], lwd = 2, col = 'tomato3')
}

plot.roc = function(x, add = F, col = 'springgreen4', ...) {
  if (!add) 
    plot(c(0, 0), xlim = c(0, 1), ylim = c(0, 1),
         xlab = 'Error Rate', ylab = 'Sensitivity', ...)
  data = x[!is.na(x[, 3]), 3:2]
  data[, 1] = 1 - data[, 1]
  data = rbind(c(1, 1), c(1, data[1, 2]), data, c(data[nrow(data), 1], 0), c(0, 0))
  lines(c(0,1), c(0,1), col = 'gray')
  lines(x = data[, 1], y = data[, 2], lwd = 2, col = col)
  sum(sapply(1:(nrow(data) - 1), function(i) (data[i + 1, 2] + data[i, 2]) * abs(data[i + 1, 1] - data[i, 1]) / 2))
}


results = lapply(list.files(pattern = '^(2|3).csv'), function(file) {
  read.csv(file)
})

results = lapply(list.files(pattern = '^(2|3)_\\w+.csv'), function(file) {
  read.csv(file)
})

results = lapply(list.files(pattern = '^(4|5|6).csv'), function(file) {
  read.csv(file)
})

rates = do.call(rbind, lapply(seq(0, 1, by = 0.01), function(threshold) {
  count = lapply(results, function(r) {
    r = as.matrix(r)
    sample.index = grep('sample', colnames(r))
    reference.index = grep('reference', colnames(r))
    
    score.index = grep('pvalue', colnames(r))
    if (length(score.index) > 0) {
      upper = FALSE
    }
    else {
      score.index = grep('confidence', colnames(r))
      if (length(score.index) > 0) {
        upper = TRUE
      }
      else {
        score.index = grep('coefficient', colnames(r))
        upper = TRUE
      }
    }
    
    count = lapply(1:nrow(r), function(i) {
      sample.species = sapply(strsplit(r[i, sample.index], ' '), function(x) paste(x[1], x[2]))
      reference.species = sapply(strsplit(r[i, reference.index], ' '), function(x) paste(x[1], x[2]))
      scores = as.numeric(r[i, score.index])
      reference.species = reference.species[if (upper) scores >= threshold else scores <= threshold]
      # decoy.index = grep('!DECOY', reference.species)
      # if (length(decoy.index) > 0)
      #   reference.species = reference.species[1:(decoy.index[1] - 1)]
      sample.count = length(sample.species)
      detected.count = length(reference.species)
      correct.count = sum(sample.species %in% reference.species)
      c(sample = sample.count, detected = detected.count, correct = correct.count)
    })
    
    c(sample = sum(sapply(count, function(x) x['sample'])),
      detected = sum(sapply(count, function(x) x['detected'])),
      correct = sum(sapply(count, function(x) x['correct'])))
  })
  sample.count = sum(sapply(count, function(x) x['sample']))
  detected.count = sum(sapply(count, function(x) x['detected']))
  correct.count = sum(sapply(count, function(x) x['correct']))
  
  c(threshold = threshold,
    detection = correct.count / sample.count,
    correct = correct.count / detected.count)
}))

plot.sensitivity.errorrate(rates, xlab = 'Threshold', ylab = 'Value')
local({
  idx = which.min(abs(rates[, 3] - 0.95))
  lines(c(rates[idx, 1], rates[idx, 1]), c(-0.2, 1.2), lwd = 2, lty = 2, col = 'navy')
  rates[(idx - 1):(idx + 1), ]
})

plot.roc(rates)
