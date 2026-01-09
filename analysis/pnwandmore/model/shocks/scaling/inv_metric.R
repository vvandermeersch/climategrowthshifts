
par(mfcol = c(4,3))
winfo <- extract_warmup_info(fit8stands_ov)
for(i in 1:4){
  
  main = ifelse(i == 1, '8 stands\n(initial version)', '')
  hist(winfo$inv_metric[[1]], main = main, xlab = 'inv_metric elements',
       col = '#9EB9F3', border = 'white')
  
}

winfo <- extract_warmup_info(fit8stands)
for(i in 1:4){
  
  main = ifelse(i == 1, '8 stands\n(4 cores)', '')
  hist(winfo$inv_metric[[1]], main = main, xlab = 'inv_metric elements',
       col = '#9EB9F3', border = 'white')
  
}

winfo <- extract_warmup_info(fit16stands)
for(i in 1:4){
  
  main = ifelse(i == 1, '16 stands\n(4 cores)', '')
  hist(winfo$inv_metric[[1]], main = main, xlab = 'inv_metric elements',
       col = '#9EB9F3', border = 'white')
  
}
