
# The simulation studies presented below are based on the article 
# http://dergipark.gov.tr/download/article-file/343265

#-----------------------------------------------------------------------------------------#


# Study 1:  Item Parameter Recovery in IRT

install.packages("mirt")
library("mirt")

itemrecovery <- function(nitem, sample.size, seed) {
  #Set the seed and generate the parameters
  set.seed(seed)
  a <- as.matrix(round(rlnorm(nitem, meanlog = 0.3, sdlog = 0.2),3), ncol=1)
  b <- as.matrix(round(rnorm(nitem, mean = 0, sd = 1),3), ncol=1)
  c <- as.matrix(round(rbeta(nitem, shape1 = 20, shape2 = 90),3), ncol=1)
  ability <- as.matrix(round(rnorm(sample.size, mean = 0, sd = 1),3), ncol=1)
  
  #Simulate response data and estimate item parameters
  dat <- simdata(a = a, d = b, N = sample.size, itemtype = 'dich',
                 guess = c, Theta = ability)
  model3PL <- mirt(data=dat, 1, itemtype='3PL', SE=TRUE, verbose=FALSE)
  #Extract estimated item parameters and compute bias and RMSE
  parameters <- as.data.frame(coef(model3PL, simplify=TRUE)$items)
  bias.a <- round(mean(parameters[,1]-a), 3)
  bias.b <- round(mean(parameters[,2]-b), 3)
  bias.c <- round(mean(parameters[,3]-c), 3)
  rmse.a <- round(sqrt(mean((parameters[,1]-a)^2)), 3)
  rmse.b <- round(sqrt(mean((parameters[,2]-b)^2)), 3)
  rmse.c <- round(sqrt(mean((parameters[,3]-c)^2)), 3)
  
  #Combine the results in a single data set
  result <- data.frame(sample.size=sample.size, nitem=nitem,
                       bias.a=bias.a, bias.b=bias.b, bias.c=bias.c,
                       rmse.a=rmse.a, rmse.b=rmse.b, rmse.c=rmse.c)
  return(result)
}


#Generate 100 random integers as seeds
myseed <- sample.int(n = 1000000, size = 100)
write.csv(myseed, "simulation seeds.txt", row.names = FALSE)

#Define an empty data frame to store the simulation results
result <- data.frame(sample.size=0, nitem=0, bias.a=0, bias.b=0, bias.c=0,
                     rmse.a=0, rmse.b=0, rmse.c=0)
#Run the loop and return the results across 100 iterations
for (i in 1:length(myseed)) {
  result[i,] <- itemrecovery(nitem = 20, sample.size = 1000, seed = myseed[i])
}
round(colMeans(result), 3)

#-------------------------------------------------------------------------------------#

# Study 2: Detecting Differential Item Functioning in Multidimensional IRT

install.packages("doParallel")
library("doParallel")
library("mirt")
library("MASS")

detectDIF <- function(sample.size, DIF.size, cor, seed) {
  
  require("mirt")
  require("MASS")
  set.seed(seed)
  
  #Define multidimensional abilities for reference and focal groups
  theta.ref <- mvrnorm(n = sample.size[1], rep(0, 2), matrix(c(1,cor,cor,1),2,2))
  theta.foc <- mvrnorm(n = sample.size[2], rep(0, 2), matrix(c(1,cor,cor,1),2,2))
  
  #Generate item parameters for reference and focal groups
  a1 <- c(runif(n = 15, min = 1.1, max = 2.8), rep(0,15))
  a2 <- c(rep(0,15), runif(n = 15, min = 1.1, max = 2.8))
  a.ref <- as.matrix(cbind(a1, a2), ncol = 2)
  b1 <- runif(n = 30, min = 0.67, max = 2)
  b2 <- b1 - runif(n = 30, min = 0.67, max = 1.34)
  b3 <- b2 - runif(n = 30, min = 0.67, max = 1.34)
  b.ref <- as.matrix(cbind(b1, b2, b3), ncol = 3)
  
  #Uniform and nonuniform DIF for items 1, 7, 15, 16, 23, and 30
  b.foc <- b.ref
  b.foc[c(1,7,15,16,23,30),] <- b.foc[c(1,7,15,16,23,30),]+DIF.size[1]
  a.foc <- a.ref
  a.foc[c(1,7,15),1] <- a.foc[c(1,7,15),1]+DIF.size[2]
  a.foc[c(16,23,30),2] <- a.foc[c(16,23,30),2]+DIF.size[2]
  #Generate item responses according to MGRM
  ref <- simdata(a = a.ref, d = b.ref, itemtype = 'graded', Theta = theta.ref)
  foc <- simdata(a = a.foc, d = b.foc, itemtype = 'graded', Theta = theta.foc)
  dat <- rbind(ref, foc)
  #Define the group variable (0=reference; 1=focal) and test DIF using mirt
  group <- c(rep("0", sample.size[1]), rep("1", sample.size[2]))
  itemnames <- colnames(dat)
  model <- 'f1 = 1-15
            f2 = 16-30
            COV = F1*F2'
  model.mgrm <- mirt.model(model)
  #Test uniform DIF
  if(DIF.size[1]>0 & DIF.size[2]==0) {
    mod0 <- multipleGroup(data = dat, model = model.mgrm, group = group,
                          invariance = c(itemnames[-c(1,7,15,16,23,30)],
                                         'free_means', 'free_var'), verbose = FALSE)
    mod1 <- DIF(mod0, c('d1','d2','d3'), items2test = c(1,7,15,16,23,30))
    result <- data.frame(items=c(1,7,15,16,23,30),
                         DIF=c(mod1[[1]][2,8], mod1[[2]][2,8], mod1[[3]][2,8],
                               mod1[[4]][2,8], mod1[[5]][2,8], mod1[[6]][2,8]))
  } else
    
    #Test nonuniform DIF
    if(DIF.size[1]==0 & DIF.size[2]>0) {
      mod0 <- multipleGroup(data = dat, model = model.mgrm, group = group,
                            invariance = c(itemnames[-c(1,7,15,16,23,30)],
                                           'free_means', 'free_var'), verbose = FALSE)
      mod1 <- DIF(mod0, c('a1'), items2test = c(1,7,15))
      mod2 <- DIF(mod0, c('a2'), items2test = c(16,23,30))
      result <- data.frame(items=c(1,7,15,16,23,30),
                           DIF=c(mod1[[1]][2,8], mod1[[2]][2,8], mod1[[3]][2,8],
                                 mod2[[1]][2,8], mod2[[2]][2,8], mod2[[3]][2,8]))
    } else
      
      #Test type I error
      if(DIF.size[1]==0 & DIF.size[2]==0) {
        mod0 <- multipleGroup(data = dat, model = model.mgrm, group = group,
                              invariance = c(itemnames[-c(1,7,15,16,23,30)],
                                             'free_means', 'free_var'), verbose = FALSE)
        mod1 <- DIF(mod0, c('a1','d1','d2','d3'), items2test = c(1,7,15))
        mod2 <- DIF(mod0, c('a2','d1','d2','d3'), items2test = c(16,23,30))
        result <- data.frame(items=c(1,7,15,16,23,30),
                             DIF=c(mod1[[1]][2,8], mod1[[2]][2,8], mod1[[3]][2,8],
                                   mod2[[1]][2,8], mod2[[2]][2,8], mod2[[3]][2,8]))
      }
  return(result)
}

# Set the seed
myseed <- sample.int(n=10000, size = 30)
detectCores()
registerDoParallel(4)

# Run the simulation using parallel computing
result <- foreach (i = 1:30, .combine=rbind) %dopar% {
  detectDIF(sample.size=c(1000, 200), DIF.size=c(0.3, 0), cor=0.3, seed=myseed[i])
}
round(mean(ifelse(result$DIF < 0.05, 1, 0)),2)

#-------------------------------------------------------------------------------------#

# Study 3: Investigating Unidimensionality

install.packages("sirt")
library("sirt")
library("mirt")
library("MASS")

detectDIM <- function(sample.size, testLength1, testLength2, cor, seed) {
  
  require("mirt")
  require("MASS")
  require("sirt")
  set.seed(seed)
  
  # Generate item discrimination and difficulty parameters
  a1 <- c(runif(n = testLength1, min = 1.1, max = 2.8), rep(0, testLength2))
  a2 <- c(rep(0, testLength1), runif(n = testLength2, min = 1.1, max = 2.8))
  disc.matrix <- as.matrix(cbind(a1, a2), ncol = 2)
  difficulty <- runif(n = (testLength1 + testLength2), min = 0.67, max = 2)
  
  # Specify inter-dimensional correlations
  sigma <- matrix(c(1, cor, cor, 1), 2, 2)
  
  # Simulate data
  dataset <- simdata(disc.matrix, difficulty, sample.size, itemtype = 'dich',
                     sigma = sigma)
  
  # Analyze the simulated data by fitting a unidimensional model
  noharmOneFactorSolution <- noharm.sirt(dat = dataset , dimensions = 1)
  
  # Summarize the results
  result <- data.frame(sample.size = sample.size , testLength1 = testLength1,
                       testLength2 = testLength2, cor = cor,
                       RMSEA = noharmOneFactorSolution$rmsea)
  
  return(result)
}

# Set the seed
myseed <- sample.int(n = 1000000, size = 100)
write.csv(myseed, "simulation seeds.csv", row.names = FALSE)
write.table(matrix(c("Sample Size", "Test Length 1", "Test Length 2",
                     "Correlation","Mean RMSEA"), 1, 5), "results.csv", sep = ",",
            col.names = FALSE, row.names = FALSE)

# Create an empty data frame to save the results
result <- data.frame(sample.size = 0, testLength1 = 0, testLength2 = 0,
                     cor = 0, RMSEA = 0)

# Run all the input values through loops
sample.size <- c(500, 1000, 3000)
test.length2 <- c(10, 20, 30)
correlation <- c(0.3, 0.6, 0.9)
for (ss in 1:length(sample.size)) {
  for (tl in 1:length(test.length2)) {
    for (k in 1:length(correlation)) {
      for (i in 1:length(myseed)) {
        result[i, ] <-detectDIM(sample.size = sample.size[ss],testLength1 = 30,
                                testLength2 = test.length2[tl], cor = correlation[k],
                                seed = myseed[i])
      }
      meanRMSEAs <- round(colMeans(result), 3)
      write.table(matrix(meanRMSEAs, 1, 5), "results.csv", sep = ",",
                  col.names = FALSE, row.names = FALSE, append = TRUE)
    }
  }
}


# Read in the summary results in
summary <- read.csv("results.csv", header = TRUE)
summary$Sample.Size <- factor(paste0("Sample Size=",summary$Sample.Size),
                              levels = c("Sample Size=500", "Sample Size=1000",
                                         "Sample Size=3000"))
summary$Test.Length.2 <- factor(paste0("Test Length 2=",summary$Test.Length.2))
summary$Correlation <- factor(summary$Correlation)

library(lattice)
dotplot(Correlation ~ Mean.RMSEA | Sample.Size*Test.Length.2 , data = summary,
        pch=c(2), cex=1.5, xlab = "Average RMSEA", aspect=0.5, layout = c(3,3),
        ylab="Correlation", xlim=c(0,0.3))



