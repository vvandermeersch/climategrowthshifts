(fit8stands |> 
  rstan::get_adaptation_info() |> 
  (\(x) gsub('\\# ','',x))() |> 
  strsplit('\n'))[[1]][4] |> 
  (\(x) sprintf('c(%s)',x))() |> 
  (\(x) parse(text = x))() |> 
  eval() |>
  as.vector()


winfo <- extract_warmup_info(fit8stands)

extract_warmup_info <- function(fit) {
  adapt  <- lapply(rstan::get_adaptation_info(fit), strsplit, split="\\n")
  step_size  <- lapply(adapt, function(a) as.numeric(strsplit(a[[1]], " = ")[[2]][2]))
  inv_metric <- lapply(adapt, function(a) as.numeric(strsplit(sub("^# ", "", a[[1]])[4], ", ")[[1]]))
  list(step_size=step_size, inv_metric=inv_metric)
}
