


# Compute the sums of the matrices 
# using foreach with 4 cores
## cl <- makeCluster(4)
## registerDoParallel(cl)
## start_time <- Sys.time()
## sums <- foreach(mat = matrices) %dopar% sum_matrix(mat)
## end_time <- Sys.time()
## stopCluster(cl)