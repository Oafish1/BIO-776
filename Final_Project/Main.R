# Libraries
library(corrplot)
library(GEOquery)
library(Biobase)
library(stringr)
library(tidyr)
library(BisqueRNA)
library(MuSiC)
library(xbioc)
library(immunedeconv)
library(FARDEEP)
library(CDSeq)
#library(bseqsc)
library(SCDC)
set_cibersort_binary("CIBERSORT.R")
set_cibersort_mat("LM22.txt")

# 'Construct in Desired Format' = Change for differing datasets
### Testing Dataset
test_file = "lb_coarse_compiled_DS.csv"; tag="DS"
#test_file = "lb_coarse_compiled_CIV.csv"; tag="CIV"
#test_file = "lb_coarse_compiled_CIAS5.csv"; tag="CIAS5"
#test_file = "lb_coarse_compiled_CIAS6.csv"; tag="CIAS6"
test_dir = "GEO Data/Leaderboard-Compiled/"

### SC Dataset
Sgse = getGEO("GSE28490")
Sgse = Sgse[[1]]
celCol="cell type:ch1"
subCol="source_name_ch1"
# Construct in Desired Format
Sex = exprs(Sgse)[,-c(31:39,43:47)]
Spd = phenoData(Sgse)[-c(31:39,43:47)]
# Get Subject/Pool ids
Sids = paste("pool",as.integer(str_extract(Spd[[subCol]], '[0-9]+$')))
Spd[[subCol]] = Sids

# Sgse = getGEO("GSE72642")
# Sgse = Sgse[[1]]
# celCol="cell type:ch1"
# subCol="source_name_ch1"
# # Construct in Desired Format
# Sex = exprs(Sgse)[,-c(6,12,18)]
# Spd = phenoData(Sgse)[-c(6,12,18)]
# # Get Subject/Pool ids
# Sids = rep(1,5)*1
# Sids = append(Sids,rep(2,5))
# Sids = append(Sids,rep(3,5))
# Spd[[subCol]] = Sids

# Sgse = getGEO("GSE110085")
# Sgse = Sgse[[1]]
# celCol="cell type:ch1"
# subCol="source_name_ch1"
# set = "DS468"

# Modify Dataset
Sgse = ExpressionSet(Sex, phenoData=Spd)

### Bulk Dataset
Bpd = read.csv(file=paste(test_dir,test_file,sep=""), 
							 header=TRUE)
for (set in unique(Bpd["dataset.name"])[,1])
{
	f = paste(test_dir, set, "-native-gene-expr.csv",
						sep="")
	print(f)
	Bex_temp = read.csv(file=f, 
								 header=TRUE)
	rownames(Bex_temp) = Bex_temp[,1]
	colnames(Bex_temp) = paste(set,colnames(Bex_temp),sep="-")
	Bex_temp = Bex_temp[,-1]
	if (set != unique(Bpd["dataset.name"])[1,1])
		Bex = cbind(Bex, Bex_temp)
	else
		Bex = Bex_temp
}
Bex = Bex[rowSums(Bex)!=0,]

frac=1
# Random Selection
set.seed(42)
Bex = Bex[sample(1:nrow(Bex), round(frac*nrow(Bex)), replace=FALSE),]

Bgse = ExpressionSet(as.matrix(Bex))

### Hugo Bulk Dataset (CIBERSORT)
# Bpd = read.csv(file=paste(test_dir,test_file,sep=""), 
# 							 header=TRUE)
# for (set in unique(Bpd["dataset.name"])[,1])
# {
# 	f = paste(test_dir, set, "-hugo-gene-expr.csv",
# 						sep="")
# 	print(f)
# 	Bex_temp = read.csv(file=f, 
# 											header=TRUE)
# 	rownames(Bex_temp) = Bex_temp[,1]
# 	colnames(Bex_temp) = paste(set,colnames(Bex_temp),sep="-")
# 	Bex_temp = Bex_temp[,-1]
# 	if (set != unique(Bpd["dataset.name"])[1,1])
# 		Bex = cbind(Bex, Bex_temp)
# 	else
# 		Bex = Bex_temp
# }
# Bex = Bex[rowSums(Bex)!=0,]
# 
# Bgse_hugo = ExpressionSet(as.matrix(Bex))

### True Dataset
trueData = function(set)
{
	Bpd = read.csv(file=paste(test_dir,test_file,sep=""), 
								header=TRUE)
	#Bpd = Bpd[Bpd[,1]==set,]
	sample.id = paste(as.matrix(Bpd["dataset.name"]),as.matrix(Bpd["sample.id"]),sep="-")
	Bpd["sample.id"] = as.data.frame(sample.id)
	Bpd = Bpd[,-1]
	Bpd = spread(Bpd, "sample.id", "measured")
	rownames(Bpd) = Bpd[,1]
	Bpd = Bpd[,-1]
	Bpd[is.na(Bpd)] = 0
	t(t(Bpd)/colSums(Bpd))[order(row.names(Bpd)),]
}
true = trueData(set)

### Run Bisque
rbd = ReferenceBasedDecomposition(Bgse, Sgse, cell.types=celCol, 
														subject.names=subCol, 
														use.overlap=FALSE)
pred_BI = rbd$bulk.props
pred_BI = pred_BI[order(row.names(pred_BI)),]

### Run MuSiC
rbd = music_prop(Bgse, Sgse, clusters=celCol, samples=subCol)
pred_MU = t(rbd$Est.prop.weighted)
pred_MU = pred_MU[order(row.names(pred_MU)),]

### Run CIBERSORT
# pred_CI = deconvolute(exprs(Bgse_hugo), "cibersort")

### Run SCDC
Spd[[subCol]] = "1"
Sgse_SC = ExpressionSet(Sex, phenoData=Spd)
pred_SC = SCDC_prop_ONE(Bgse,Sgse_SC,ct.varname=celCol,sample=subCol,ct.sub=unique(Sgse[[celCol]]))
pred_SC = t(pred_SC$prop.est.mvw)
pred_SC = pred_SC[order(row.names(pred_SC)),]

### Display Sample Results
to_test = c("DS468-S1","DS468-S10")
pred_BI[,to_test]
pred_MU[,to_test]
pred_SC[,to_test]
true[,to_test]

### Modify Format
# Construct in Desired Format
#GEO2
true = true[-c(4:5),]
#GEO7
# true = true[-c(4:5),]
# true = true[-5,]
# true = t(t(true)/colSums(true))
# pred_BI = pred_BI[c(2,3,5,1,4),]
# pred_MU = pred_MU[c(2,3,5,1,4),]
# pred_SC = pred_SC[c(2,3,5,1,4),]

### Assess Error
err = function(predicted, actual)
{
	sqdiff = (actual - predicted)^2
	MSE = colSums(sqdiff,na.rm=TRUE)/dim(sqdiff)[1]
	return(list("MSE"=MSE, "sqdiff"=sqdiff))
}
ERR_BI = err(pred_BI, true)
ERR_MU = err(pred_MU, true)
ERR_SC = err(pred_SC, true)
ERR_BI$MSE[to_test]
ERR_MU$MSE[to_test]
ERR_SC$MSE[to_test]

### Generate Plots
ymax=.3
# Error by Sample
plot_mean = (ERR_BI$MSE+ERR_MU$MSE+ERR_SC$MSE)/3
plot(ERR_BI$MSE[order(plot_mean)], type="p", col="red", pch=24, bg="red",
		 xlab="Samples", ylab="MSE",
		 main=paste(tag,": ","Error by Sample",sep=""), ylim=c(0,ymax), xaxt = "n",
		 cex.axis=2,cex.lab=2)
abline(median(ERR_BI$MSE),0, col="red")

points(ERR_MU$MSE[order(plot_mean)], type="p", col="green", pch=22, bg="green")
abline(median(ERR_MU$MSE),0, col="green")

points(ERR_SC$MSE[order(plot_mean)], type="p", col="blue", pch=21, bg="blue")
abline(median(ERR_MU$MSE),0, col="blue")

legend(1,ymax,c("Bisque", "MuSiC", "SCDC", "Median"),
			 col=c("red", "green", "blue", "black"), cex=2,
			 pt.bg=c("red", "green", "blue", "black"), pch=c(24,22,21,-1), lty=c(0,0,0,1))

# Error by Dataset
error = integer(length(unique(Bpd["dataset.name"])[,1]))
to_err = ERR_BI$MSE
for (i in 1:length(unique(Bpd["dataset.name"])[,1]))
{
	set = unique(Bpd["dataset.name"])[,1][i]
	matches = grepl(set,names(to_err),fixed=TRUE)
	error[i] = sum(matches * to_err)/sum(matches)
}
plot(error, type="p", col="red", pch=24, bg="red",
		 xlab="Dataset", ylab="MSE", main=paste(tag,": ","Error by Sample",sep=""),
		 ylim=c(0,ymax), xaxt = "n")
abline(median(error),0, col="red")

error = integer(length(unique(Bpd["dataset.name"])[,1]))
to_err = ERR_MU$MSE
for (i in 1:length(unique(Bpd["dataset.name"])[,1]))
{
	set = unique(Bpd["dataset.name"])[,1][i]
	matches = grepl(set,names(to_err),fixed=TRUE)
	error[i] = sum(matches * to_err)/sum(matches)
}
points(error, type="p", col="green", pch=22, bg="green")
abline(median(error),0, col="green")

error = integer(length(unique(Bpd["dataset.name"])[,1]))
to_err = ERR_SC$MSE
for (i in 1:length(unique(Bpd["dataset.name"])[,1]))
{
	set = unique(Bpd["dataset.name"])[,1][i]
	matches = grepl(set,names(to_err),fixed=TRUE)
	error[i] = sum(matches * to_err)/sum(matches)
}
points(error, type="p", col="blue", pch=21, bg="blue")
abline(median(error),0, col="blue")

axis(1, at=1:length(unique(Bpd["dataset.name"])[,1]), labels=unique(Bpd["dataset.name"])[,1], las=2)
legend(1,ymax,c("Bisque", "MuSiC", "SCDC", "Median"),
			 col=c("red", "green", "blue", "black"), cex=0.8,
			 pt.bg=c("red", "green", "blue", "black"), pch=c(24,22,21,-1), lty=c(0,0,0,1))

# Error by Cell Type
error = rowSums(ERR_BI$sqdiff,na.rm=TRUE)/dim(ERR_BI$sqdiff)[2]
par(mar = c(10, 4, 4, 4) + 1)
plot(error, type="p", col="red", pch=24, bg="red",
		 xlab="Cell", ylab="MSE", main=paste(tag,": ","Error by Cell Type",sep=""), ylim=c(0,ymax),
		 xaxt = "n",cex.lab=1.5, cex.axis=1.5)
abline(median(error),0, col="red")

error = rowSums(ERR_MU$sqdiff,na.rm=TRUE)/dim(ERR_MU$sqdiff)[2]
points(error, type="p", col="green", pch=22, bg="green")
abline(median(error),0, col="green")

error = rowSums(ERR_SC$sqdiff,na.rm=TRUE)/dim(ERR_SC$sqdiff)[2]
points(error, type="p", col="blue", pch=22, bg="blue")
abline(median(error),0, col="blue")
axis(1, at=1:dim(ERR_BI$sqdiff)[1], labels=names(error), las=2, cex.axis=1.5)

legend(1,ymax,c("Bisque", "MuSiC", "SCDC", "Median", "Average Cell Pct"),
			 col=c("red", "green", "blue", "black", "black"), cex=1,
			 pt.bg=c("red", "green", "blue", "black", "black"), pch=c(24,22,21,-1,4),
			 lty=c(0,0,0,1,0))

pct = rowSums(true,na.rm=TRUE)/dim(true)[2]
par(new = TRUE)
plot(pct, type="p", col="black", pch=4, bg="black",
		 xaxt="n", yaxt="n", xlab="", ylab="", ylim=c(0,1), cex=1.5)
axis(side=4, at=c(0,.25,.5,.75,1), cex.axis=1.5)
mtext("Average Cell Percentage", side = 4, line = 3, cex=1.5)

### Corrplot
temp = exprs(Sgse)
colnames(temp) = Sgse[[celCol]]
corrplot(cor(temp[,c(1,11,16,21,26,31)]),tl.col="black")
# rownames(temp)=NULL
# res = as.data.frame(sapply(unique(colnames(temp)), function(x) rowMeans(temp[colnames(temp) == x],na.rm=TRUE)))
