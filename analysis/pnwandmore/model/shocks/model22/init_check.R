fit <- model$sample(data = data, threads_per_chain = 1,
                    chains = 4, parallel_chains = 4, seed = 1234567,
                    refresh = 10,
                    iter_warmup = 1000, iter_sampling = 1000,
                    save_warmup=FALSE)


init_check <- model$sample(
  data = data,
  chains = 4,
  parallel_chains = 4,
  threads_per_chain = 1,
  seed = 1234567,
  fixed_param = TRUE,
  iter_warmup = 0,
  iter_sampling = 1
)

init_check$draws()