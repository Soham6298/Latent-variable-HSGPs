## This piece of code is to transfer the simulated data generated in R to pyro for fitting VIGPs

simdata <- readRDS('hsgp_simdata_se.rds')
for ( i in 1:n_sim) {
  for(j in 1:length(N)) {
    if (j == 1) {
      for (k in 1:length(dims)){
        if (k == 1) {
          filename <- sprintf("~/GPLVMpyro/Simdata20N/Simdata5D/Simdata[20ND5]%04d.csv", i)  # Pads numbers with 2 leading zeros
          write.csv(simdata[[i]][[j]][[k]], filename)
        }
        else if (k == 2) {
          filename <- sprintf("~/GPLVMpyro/Simdata20N/Simdata10D/Simdata[20ND10]%04d.csv", i)  # Pads numbers with 2 leading zeros
          write.csv(simdata[[i]][[j]][[k]], filename)
        }
        else {
          filename <- sprintf("~/GPLVMpyro/Simdata20N/Simdata20D/Simdata[20ND20]%04d.csv", i)  # Pads numbers with 2 leading zeros
          write.csv(simdata[[i]][[j]][[k]], filename)
        }
      }
    }
    else if (j == 2) {
      for (k in 1:length(dims)){
        if (k == 1) {
          filename <- sprintf("~/GPLVMpyro/Simdata50N/Simdata5D/Simdata[50ND5]%04d.csv", i)  # Pads numbers with 2 leading zeros
          write.csv(simdata[[i]][[j]][[k]], filename)
        }
        else if (k == 2) {
          filename <- sprintf("~/GPLVMpyro/Simdata50N/Simdata10D/Simdata[50ND10]%04d.csv", i)  # Pads numbers with 2 leading zeros
          write.csv(simdata[[i]][[j]][[k]], filename)
        }
        else {
          filename <- sprintf("~/GPLVMpyro/Simdata50N/Simdata20D/Simdata[50ND20]%04d.csv", i)  # Pads numbers with 2 leading zeros
          write.csv(simdata[[i]][[j]][[k]], filename)
        }
      }
    }
    else {
      for (k in 1:length(dims)){
        if (k == 1) {
          filename <- sprintf("~/GPLVMpyro/Simdata200N/Simdata5D/Simdata[200ND5]%04d.csv", i)  # Pads numbers with 2 leading zeros
          write.csv(simdata[[i]][[j]][[k]], filename)
        }
        else if (k == 2) {
          filename <- sprintf("~/GPLVMpyro/Simdata200N/Simdata10D/Simdata[200ND10]%04d.csv", i)  # Pads numbers with 2 leading zeros
          write.csv(simdata[[i]][[j]][[k]], filename)
        }
        else {
          filename <- sprintf("~/GPLVMpyro/Simdata200N/Simdata20D/Simdata[200ND20]%04d.csv", i)  # Pads numbers with 2 leading zeros
          write.csv(simdata[[i]][[j]][[k]], filename)
        }
      }
    }
  }
}
