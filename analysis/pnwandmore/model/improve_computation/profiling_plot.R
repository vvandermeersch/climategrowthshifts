time_5threads <- readRDS(file.path(wd, 'model/improve_computation/output/runtime', 'time_5threads.rds'))
time_5threads[5,2]
sum(as.numeric(time_5threads[1:4,2]))
time_5threads$modality <- '5 threads'

time_5threads_v2 <- readRDS(file.path(wd, 'model/improve_computation/output/runtime', 'time_5threads_v2.rds'))
time_5threads_v2[5,2]
sum(as.numeric(time_5threads_v2[1:4,2]))
time_5threads_v2$modality <- '5 threads, v2'

time_10threads_v2 <- readRDS(file.path(wd, 'model/improve_computation/output/runtime', 'time_10threads_v2.rds'))
time_10threads_v2[5,2]
sum(as.numeric(time_10threads_v2[1:4,2]))
time_10threads_v2$modality <- '10 threads, v2'


