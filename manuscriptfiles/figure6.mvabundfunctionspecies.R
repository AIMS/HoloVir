# Reading the data:
allspongedna=read.csv("revisedfunctions.swissprotkeywords.4mvabund.species.csv");
dim(allspongedna)
# That told me 10 rows and 4442 columns, which makes sense cause I have the 4439 columns of data and the additional columns of bubbleornot and species. 
allspongedna[,2:806]
# That's bubbleornot and species
library(vegan)
# Rarefying the data:
allspongednarare=rrarefy(allspongedna[,3:806], min(rowSums(allspongedna[3:806])))
rowSums(allspongednarare)
# Data was rarefied to 97068 cogs
library(mvabund)
# checking the relationship between the mean and the variance:
meanvar.plot(allspongednarare)
abline(0,1)
# and saving this plot:
dev.copy2eps(file="allspongedna.refseqgenes.meanvarplot.eps")
dev.copy2pdf(file="allspongedna.refseqgenes.meanvarplot.eps")
# Defining what "x" will be in the manyglm model:
xallspongedna=allspongedna[,2:3]
xallspongedna
# Making an mvabund object, which will be used to make glm models:
yallspongednarare=mvabund(allspongednarare)
# Making the gml models:
modelallspongedna=manyglm(yallspongednarare~species,data=xallspongedna)
save(modelallspongedna, file="modelallspongedna")
save(allspongednarare, file="allspongednarare")
save(xallspongedna,file="xallspongedna")
save(yallspongednarare, file="yallspongednarare")
# Plotting the manyglm models check if the distribution is "normal" (i.e you shouldn't see any pattern, the residues should be distributed equally in this space). 
# Also, when you do this, don't use all the residues, just use 50 chosen randomly. The formula is like this:
plot.manyglm(modelallspongedna, var.subset = sample(1:ncol(yallspongednarare), 50))
# Saving the plot:
dev.copy2eps(file="plotmanyglm_yallspongednafunctrare_a.eps")
dev.copy2pdf(file="plotmanyglm_yallspongednafunctrare_a.pdf")
# ...and making 3 more plots, just for replication purposes:
plot.manyglm(modelallspongedna, var.subset = sample(1:ncol(yallspongednarare), 50))
dev.copy2pdf(file="plotmanyglm_yallspongednafunctrare_b.pdf")
dev.copy2eps(file="plotmanyglm_yallspongednafunctrare_b.eps")
#plot.manyglm(modelcoelo4439, var.subset = sample(1:ncol(ycoelo4439rare), 50))
plot.manyglm(modelallspongedna, var.subset = sample(1:ncol(yallspongednarare), 50))
dev.copy2pdf(file="plotmanyglm_yallspongednafunctrare_c.pdf")
dev.copy2eps(file="plotmanyglm_yallspongednafunctrare_c.eps")
plot.manyglm(modelallspongedna, var.subset = sample(1:ncol(yallspongednarare), 50))
dev.copy2pdf(file="plotmanyglm_yallspongednafunctrare_d.pdf")
dev.copy2eps(file="plotmanyglm_yallspongednafunctrare_d.eps")
# and then i went to katana to run the anova
allspongednafunct_anova_unadjust_999bts=anova(modelallspongedna, p.uni = "unadjusted", nBoot=999, show.time = "all")
save(allspongednafunct_anova_unadjust_999bts, file="allspongednafunct_anova_unadjust_999bts")
allspongednafunct_anova_unadjust_999bts$table

# check the p-value for each of the cogs:
anovaunipunadj=allspongednafunct_anova_unadjust_999bts$uni.p
save(anovaunipunadj, file="anovaunipunadj")
# get the names of the cogs only, without the p-value:
anovaunipunadj_pvalue=anovaunipunadj["species",]
anovaunipunadj_pvalue
# get the names of the cogs for which the p-value is less than 0.02 (I did this with 0.05, 0.02 and 0.01):
namesof_functions_pvalueless002=names(which(anovaunipunadj_pvalue<0.02))
save(namesof_functions_pvalueless002, file="namesof_functions_pvalueless002")
# extract  from your original rarefied data file, the cogs that match these names:
allspongednarare_pvalueless002=allspongednarare[,colnames(allspongednarare) %in% namesof_functions_pvalueless002]
save(allspongednarare_pvalueless002, file="allspongednafunctrare_pvalueless002")
write.csv(allspongednarare_pvalueless002, file="allspongednafunctrare_pvalueless002.csv")
# make this a dataframe:
allspongednarare_pvalueless002_df=data.frame(allspongednarare_pvalueless002)
is.data.frame(allspongednarare_pvalueless002_df)
# add rownames to the dataframe (Patrick, this is some R stuff that you might have to adjust to your dataset. Maybe it's just my R skills that are not good enough yet).
rownames(allspongednarare_pvalueless002_df)=allspongedna$X
namesof_functions_pvalueless002=names(which(anovaunipunadj_pvalue<0.02))
save(namesof_functions_pvalueless002, file="namesof_functions_pvalueless002")
allspongednarare_pvalueless002=allspongednarare[,colnames(allspongednarare) %in% namesof_functions_pvalueless002]
save(allspongednarare_pvalueless002, file="allspongednafunctrare_pvalueless002")
write.csv(allspongednarare_pvalueless002, file="allspongednafunctrare_pvalueless002.csv")
allspongednarare_pvalueless002_df=data.frame(allspongednarare_pvalueless002)
is.data.frame(allspongednarare_pvalueless002_df)
rownames(allspongednarare_pvalueless002_df)=allspongedna$X
########################################################################
########################################################################
library(pheatmap)
pheatmap(log10(allspongednarare_pvalueless002_df+1), cluster_rows=F, cluster_cols=T, fontsize=5)
dev.copy2eps(file="pheatmap_allspongednafunctrare_pvalueless002_log.eps")
dev.copy2pdf(file="pheatmap_allspongednafunctrare_pvalueless002_log.pdf")
# make a heatmap with the log10 transformed data to visualise it better:
pheatmap(log10(allspongednarare_pvalueless002_df+1), cluster_rows=F, cluster_cols=T, fontsize=7)
dev.copy2eps(file="pheatmap_allspongednafunctrare_pvalueless002_log.eps")
dev.copy2pdf(file="pheatmap_allspongednafunctrare_pvalueless002_log.pdf")
# and then this is Shaun's trick to see the standard deviation between "seep" and "control" instead of the rarefied counts:
pheatmap(log10(allspongednarare_pvalueless002_df+1), cluster_rows=F, cluster_cols=T, fontsize=6)
dev.copy2pdf(file="pheatmap_allspongednafunctrare_pvalueless002_scalecol.pdf")
dev.copy2eps(file="pheatmap_allspongednafunctrare_pvalueless002_scalecol.eps")
#########################################################################
#########################################################################
# I now have to get the order of the cogs on the heatmap. This is actually not done directly from pheatmap. You have to go back to clustering in order to do that:
function_distance_matrix_less002=dist(t(allspongednarare_pvalueless002_df), method="euclidean")
# this last script stands for: "make a distance matrix of the TRANSPOSED ("t) dataframe called "logtenof_coeloblahblah" using the "euclidean" method."
function_hcluster_less002=hclust(function_distance_matrix_less002, "complete")
# and this last one stands for "do a COMPLETE clustering of the distance matrix called "cogs_distance_matrixblahblah".
function_hcluster_less002
function_hcluster_less002$labels
# this last one showed you the labels of the cluster, but they're not in the same order as what they are in the actual cluster. To get that you should do:
function_hcluster_less002$labels[function_hcluster_less002$order]
# and then you can name this:
hclust002_names_ordered_allspongedna=function_hcluster_less002$labels[function_hcluster_less002$order]
save(hclust002_names_ordered_allspongedna, file="hclust002_names_ordered_allspongedna")
write.table(hclust002_names_ordered_allspongedna, file="hclust002_names_ordered_allspongedna.txt")
# !!!! WARNING!!!! I noticed that the order of the cogs on the heatmap and the order of the cogs from the cluster are not always exactly the same. 
# It makes sense, there are so many cogs, that the clustering is probably not EXACTLY the same each time. 
# It is a problem though and so far I don't know how to fix it. 
save(function_distance_matrix_less002, file="function_distance_matrix_less002allspongedna")
save(function_hcluster_less002, file="function_hcluster_less002")
########################################################################
namesof_functions_pvalueless005=names(which(anovaunipunadj_pvalue<0.05))
save(namesof_functions_pvalueless005, file="namesof_functions_pvalueless005")
# extract  from your original rarefied data file, the cogs that match these names:
allspongednarare_pvalueless005=allspongednarare[,colnames(allspongednarare) %in% namesof_functions_pvalueless005]
save(allspongednarare_pvalueless005, file="allspongednafunctrare_pvalueless005")
write.csv(allspongednarare_pvalueless005, file="allspongednafunctrare_pvalueless005.csv")
# make this a dataframe:
allspongednarare_pvalueless005_df=data.frame(allspongednarare_pvalueless005)
is.data.frame(allspongednarare_pvalueless005_df)
# add rownames to the dataframe (Patrick, this is some R stuff that you might have to adjust to your dataset. Maybe it's just my R skills that are not good enough yet).
rownames(allspongednarare_pvalueless005_df)=allspongedna$X
########################################################################
########################################################################
library(pheatmap)
# make a heatmap with the log10 transformed data to visualise it better:
pheatmap(log10(allspongednarare_pvalueless005_df+1), cluster_rows=F, cluster_cols=T, fontsize=5)
dev.copy2eps(file="pheatmap_allspongednafunctrare_pvalueless005_log.eps")
dev.copy2pdf(file="pheatmap_allspongednafunctrare_pvalueless005_log.pdf")
# and then this is Shaun's trick to see the standard deviation between "seep" and "control" instead of the rarefied counts:
pheatmap((allspongednarare_pvalueless005_df+1), cluster_rows=T, cluster_cols=T, fontsize=5)
dev.copy2pdf(file="pheatmap_allspongednafunctrare_pvalueless005_scalecol.pdf")
dev.copy2eps(file="pheatmap_allspongednafunctrare_pvalueless005_scalecol.eps")
# I now have to get the order of the cogs on the heatmap. This is actually not done directly from pheatmap. You have to go back to clustering in order to do that:
function_distance_matrix_less005=dist(t(allspongednarare_pvalueless005_df), method="euclidean")
# this last script stands for: "make a distance matrix of the TRANSPOSED ("t) dataframe called "logtenof_coeloblahblah" using the "euclidean" method."
function_hcluster_less005=hclust(function_distance_matrix_less005, "complete")
# and this last one stands for "do a COMPLETE clustering of the distance matrix called "cogs_distance_matrixblahblah".
function_hcluster_less005
function_hcluster_less005$labels
# this last one showed you the labels of the cluster, but they're not in the same order as what they are in the actual cluster. To get that you should do:
function_hcluster_less005$labels[function_hcluster_less005$order]
# and then you can name this:
hclust005_names_ordered_allspongedna=function_hcluster_less005$labels[function_hcluster_less005$order]
save(hclust005_names_ordered_allspongedna, file="hclust005_names_ordered_allspongedna")
write.table(hclust005_names_ordered_allspongedna, file="hclust005_names_ordered_allspongedna.txt")
