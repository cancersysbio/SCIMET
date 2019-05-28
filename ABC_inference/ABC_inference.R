freq_summary_u <- function(post.u){
      counts= as.matrix(table(post.u))
      rowname =  rownames(counts)
      final.u = data.frame(matrix(rep(0,10),ncol = 1, nrow = 10))
      colnames(final.u) = c("freq")
      rownames(final.u) = c("0.003","0.006","0.015","0.03","0.06","0.15","0.3","0.6","1.5","3")

      for(name in rowname){
            final.u[name,] = as.numeric(counts[name,])/sum(counts[,1])
      }
      return(final.u[,"freq"])
}


freq_summary_time <- function(post.time){
      counts= as.matrix(table(post.time))
      rowname =  rownames(counts)
      final.time = data.frame(matrix(rep(0,7),ncol = 1, nrow = 7))
      colnames(final.time) = c("freq")
      rownames(final.time) = c("3","4","5","6","7","8","9")

      for(name in rowname){
            final.time[name,] = as.numeric(counts[name,])/sum(counts[,1])
      }
      return(final.time[,"freq"])
}



########### model selection - N/N, N/S, S/N or S/S #########
#setwd("DIR/abc_summary")

setwd("/Users/huzheng/Desktop/Research/CurtisLab/mCRC_submission/mCRC_final/Codes/SCIMET/ABC_summary")
library(abc)

dat <- read.table("mCRC_9_summary_statistics.txt",header=T)

primet <- read.table("SimParamStats_4models_b0.55_280000tumors_complete.txt",header=T)
sumstat.sim = primet[,c(5:13)]
param.sim =   primet[,c(2,4)]
models <- as.vector(primet$model)

sink("mCRC_model_selection_b0.55_c1.txt")
name = c("Met","Model")
cat(name)
cat("\n")

for(k in 1:nrow(dat)){
      obs_ss <- dat[k,c(2:10)]
      #print(c("###########",as.vector(dat[k,1]),"##############"))
      model_sel <- postpr(as.numeric(obs_ss), models, sumstat.sim, tol=0.005, method="neuralnet", trace=FALSE)
      model_prob <- model_sel$pred
      optimal_model <- names(model_prob)[which(model_prob==max(model_prob))]
      cat(c(as.character(dat[k,1]),optimal_model))
      cat("\n")
}


############# Run ABC for each of the four models###############
primet <- read.table("SimParamStats_4models_b0.55_280000tumors_complete.txt",header=T)
dat <- read.table("mCRC_9_summary_statistics.txt",header=T)
met_names <- as.vector(dat$Met)
output_name = c("Met","Model", "u0.003","u0.006","u0.015","u0.03","u0.06","u0.15","u0.3","u0.6","u1.5","u3","gap","Nd_1e3","Nd_1e4","Nd_1e5","Nd_1e6","Nd_1e7","Nd_1e8","Nd_1e9", "u_median","log10.Nd_median","log10.Nd_1st_qu", "log10.Nd_3rd_qu")
	
output_file <- paste("mCRC_ABCposterior_4models_b0.55_c1.txt", sep="")
sink(output_file)
cat(output_name)
cat("\n")

models = c("NN", "NS", "SN", "SS")
for (m in models){
	primet1 <- primet[which(primet$model==m),]
	sumstat.sim <- primet1[,c(5:13)]
	param.sim <-   primet1[,c(2,4)]

	for(k in 1:nrow(dat)){
      	obs_ss <- as.numeric(dat[k,c(2:10)])
      	rej <- abc(target=obs_ss, param=param.sim, sumstat=sumstat.sim, tol=0.01, method="rejection", trace=FALSE)
      	rej <- data.frame(rej$unadj.values)
      	median_mu <- median(rej$mu)
      	median_Nd <- median(rej$met.time)
      	met_time_summary <- as.numeric(summary(rej$met.time))
      	Nd_1st_qu <- met_time_summary[2]
      	Nd_3rd_qu <- met_time_summary[5]
      
      	cat(met_names[k])
      	cat("\t")
      	cat(m)
      	cat("\t")
      	x = c()
      	x = c(x,freq_summary_u(rej$mu))
      	x = c(x,0)
      	x = c(x,freq_summary_time(rej$met.time))
      	x = c(x, median_mu, median_Nd, Nd_1st_qu, Nd_3rd_qu)
      	cat(x)
      	cat("\n")
	}

}

sink()


###########Select the posterior inference based on model selection##########

ABCposterior <- read.table("mCRC_ABCposterior_4models_b0.55_c1.txt",header=T)
ModelSel <- read.table("mCRC_model_selection_b0.55_c1.txt",header=T)

ABCposterior_optimal <- merge(ABCposterior, ModelSel, by.y=c("Met","Model"), sort = FALSE)
write.table(ABCposterior_optimal,"mCRC_ABCposterior_optimal_model_b0.55_c1.txt", sep="\t", row.names=F, quote=FALSE)


########## ploting the heatmap ##########
library(dplyr)
library(NMF)
library(RColorBrewer)


dat <- read.table("mCRC_ABCposterior_optimal_model_b0.55_c1_sorted.txt",header=T)
abcdat <- as.matrix(dat[,3:20])

colname <- c("0.003","0.006","0.015","0.03","0.06","0.15","0.3","0.6","1.5","3.0","","1e+3","1e+4","1e+5","1e+6","1e+7","1e+8","1e+9")

rowname <- dat$Met
rownames(abcdat) = rowname
colnames(abcdat) = colname

pdf("mCRC_ABCposterior_optimal_model_b0.55_c1_sorted.pdf",height=6, width=8, onefile=FALSE)
par(mar=c(5,5,4,4))

cl = colors()
colfunc <- colorRampPalette(c("white","#807dba","#54278f"))
aheatmap(abcdat,scale = "none", Rowv=NA, Colv=NA,col = colfunc(100),cexRow=1.2,cexCol=1.2,border_color="gray70")

dev.off()


