#From count data to de analysis and annotation
library(edgeR)
library(dplyr)
library(plyr)
library(pheatmap)
library(biomaRt)
library(tidyverse)
library(forcats)
#Data load, ID mapping and duplicate removal
c50=read.csv('../fifty_samples.csv',row.names = 1)
ensembl = useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl")
ids <- getBM(attributes=c('ensembl_gene_id','ensembl_transcript_id','hgnc_symbol','entrezgene', 'refseq_mrna','interpro','interpro_description'),mart = ensembl)
idc=ids[complete.cases(ids),]
c50$ensembl_gene_id=rownames(c50)
idc1=idc[match(c50$ensembl_gene_id,idc$ensembl_gene_id),]
idc1=idc1[complete.cases(idc1),]
c50i=merge(c50,idc1,by='ensembl_gene_id')
c50i=c50i[complete.cases(c50i),]
o <- order(rowSums(c50i[,2:51]), decreasing=T)
c50o=c50i[o,]
d=duplicated(c50o$hgnc_symbol)
c50f=c50o[!d,]
rownames(c50f)=c50f$hgnc_symbol
cf1=c50f
#Remove low count genes in different groups and normalization
cls1=DGEList(counts = cf1[2:51],genes = cf1[,52:57])
gls1=cls1 #Preserve the original count data to use multiple time for different experiment
cls2=cls1
cls2$samples$group=c(rep('Control',15),rep('Case',35)) #Annotate case and controls
dim(cls2$samples)
head(cls2$samples)
q1=read.csv('../DRC_quart.csv')
q1=q1[complete.cases(q1),]
q1$ID=gsub('RCA','',q1$ID) #Convert 'BRCA' to onle 'B' from sample id
cls2$samples$quartile=q1$Quart_DRC[match(rownames(cls2$samples),q1$ID)] #Add DRC info into the count data
head(cls2$samples)
keep=rowSums(cpm(cls2)>1) >=2 #Remove low count genes in case and control
#We have to do it multiple times for different group consideration
cls3=cls2[keep, , keep.lib.sizes=F]
dim(cls3$counts)
dim(cls2$counts)

cls4=calcNormFactors(cls3) #Normalization
head(cls4$samples)
cls4=estimateDisp(cls4) #Dispersion
cls4$common.dispersion
etcc=exactTest(cls4,pair = c('Control','Case')) #DE classic approach
summary(decideTests(etcc))

dcc=model.matrix(~cls3$samples$group) #Experimental group for glm approach
dcc
dcc=model.matrix(~0+cls3$samples$group)
dcc
colnames(dcc)
colnames(dcc)=c('Case','Control')
dcc
clsg=estimateDisp(cls3,design = dcc)
clsg$common.dispersion
ftcc=glmQLFit(clsg,design = dcc)
qlcc=glmQLFTest(ftcc,contrast = c(1,-1)) #DE glm approach
summary(decideTests(qlcc))
topTags(qlcc)
topTags(etcc)
tgcc=topTags(etcc,n=14841,adjust.method = 'BH',sort.by = 'PValue',p.value = 0.05)
tgcc_glm=topTags(qlcc,n=14841,adjust.method = 'BH',sort.by = 'PValue',p.value = 0.05)
png('dispersion_and_de_case_control.png')
par(mfrow=c(2,2))
plotBCV(cls4,main='Case vs Control')
plotMD(etcc)
abline(h=c(-1,1),lwd=2)
plotBCV(clsg,main='Case vs Control')
plotMD(qlcc)
abline(h=c(-1,1),lwd=2)
dev.off()
png('distribution_of_counts.png',width = 400)
plot(density(colSums(cf1[,(2:51)])), col='red', xlab = 'Observed values', main = 'Distribution of counts',lwd=2)
lines(density(colSums(all_sample)),col='blue',lwd=2)
lines(density(colSums(cls4$counts)*cls4$samples$norm.factors),col='green',lwd=2)
legend('topleft',legend = c('all samples','without outliers','normalized'),col = c('blue','red','green'),lty = 1,lwd=2)
dev.off()

#DE in different quartiles
head(cls2$samples)
dim(cls2$counts)
cls2$samples$group=cls2$samples$quartile
dim(cls2$counts)
head(cls2$samples)
keep=rowSums(cpm(cls2)>1) >=4 #Remove low count genes from all quartiles
cls5=cls2[keep, , keep.lib.sizes=F]
dim(cls5$counts)
dim(cls2$counts)
dim(cls3$counts)
head(cls5$samples)
cls5=calcNormFactors(cls5)
head(cls5$samples)

cls6=estimateDisp(cls5)
cls6$common.dispersion
dall=model.matrix(~cls5$samples$group)
dall
dall=dall[,-1]
dall
colnames(dall)=levels(cls5$samples$group)
colnames(dall)=c('Q1','Q2','Q3','Q4')
dall
cls6$common.dispersion
cls7=estimateDisp(cls5,design = dall)
cls7$common.dispersion
clq2q1=exactTest(cls6,pair = c('Q1','Q2'))
summary(decideTests(clq2q1))
ftq2q1=glmQLFit(cls7,design = dall)
head(dall)
qlq2q1=glmQLFTest(ftq2q1, contrast = c(-1,1,0,0))
summary(decideTests(qlq2q1))
clq3q1=exactTest(cls6,pair = c('Q1','Q3'))
summary(decideTests(clq3q1))
qlq3q1=glmQLFTest(ftq2q1, contrast = c(-1,0,1,0))
summary(decideTests(qlq3q1))
clq4q1=exactTest(cls6,pair = c('Q1','Q4'))
summary(decideTests(clq4q1))
qlq4q1=glmQLFTest(ftq2q1, contrast = c(-1,0,0,1))
summary(decideTests(qlq4q1))
clq3q2=exactTest(cls6,pair = c('Q2','Q3'))
summary(decideTests(clq3q2))
qlq3q2=glmQLFTest(ftq2q1, contrast = c(0,-1,1,0))
summary(decideTests(qlq3q2))
clq4q2=exactTest(cls6,pair = c('Q2','Q4'))
summary(decideTests(clq4q2))
qlq4q2=glmQLFTest(ftq2q1, contrast = c(0,-1,0,1))
summary(decideTests(qlq4q2))
clq4q3=exactTest(cls6,pair = c('Q3','Q4'))
summary(decideTests(clq4q3))
qlq4q3=glmQLFTest(ftq2q1, contrast = c(0,0,-1,1))
summary(decideTests(qlq4q3))
qlq34q12=glmQLFTest(ftq2q1, contrast = c(-1,-1,1,1))
summary(decideTests(qlq34q12))
cls7=cls6
cls7$samples$group=fct_collapse(cls7$samples$group, Q1=c('Q1','Q2'),Q3=c('Q3','Q4'))
cls7$samples$group
clq34q12=exactTest(cls7,pair = c('Q1','Q3'))

#write summary
write.csv(as.data.frame(topTags(clq3q1,n=19, adjust.method = "BH", sort.by = "PValue", p.value = 0.05)),'domain_classic_q3q1.csv')
write.csv(as.data.frame(topTags(qlq3q1,n=19, adjust.method = "BH", sort.by = "PValue", p.value = 0.5)),'domain_glm_q3q1.csv')
write.csv(as.data.frame(topTags(qlq4q1,n=54, adjust.method = "BH", sort.by = "PValue", p.value = 0.5)),'domain_glm_q4q1.csv')
write.csv(as.data.frame(topTags(clq4q1,n=54, adjust.method = "BH", sort.by = "PValue", p.value = 0.05)),'domain_classic_q4q1.csv')
write.csv(as.data.frame(topTags(clq3q2, adjust.method = "BH", sort.by = "PValue", p.value = 0.05)),'domain_classic_q3q2.csv')
write.csv(as.data.frame(topTags(clq4q1, adjust.method = "BH", sort.by = "PValue", p.value = 0.5)),'domain_glm_q3q2.csv')
write.csv(as.data.frame(topTags(clq4q3, adjust.method = "BH", sort.by = "PValue")),'domain_classic_q4q3.csv')
write.csv(as.data.frame(topTags(qlq4q3, adjust.method = "BH", sort.by = "PValue")),'domain_glm_q4q3.csv')
write.csv(as.data.frame(topTags(qlq34q12, n=30, adjust.method = "BH", sort.by = "PValue")),'domain_glm_q34q12.csv')
write.csv(as.data.frame(topTags(clq34q12,n=30, adjust.method = "BH", sort.by = "PValue")),'domain_classic_q34q12.csv')
write.csv(as.data.frame(topTags(qlq34q12c, n=30, adjust.method = "BH", sort.by = "PValue")),'domain_glm_q34q12c.csv')
write.csv(as.data.frame(topTags(clq34q12, n=30, adjust.method = "BH", sort.by = "PValue")),'domain_classic_q34q12c.csv')

#Case only
cs1=DGEList(counts = cf1[,17:51],genes = cf1[,52:57])
cs1$samples
cs1$samples$group=cls6$samples$group[match(rownames(cs1$samples),rownames(cls6$samples))]
cs1$samples
keep=rowSums(cpm(cs1)>1) >=4
cs2=cs1[keep, , keep.lib.sizes=F]
cs1$samples
cs2$samples
cs2=calcNormFactors(cs1)
cs2$samples
cs3=estimateDisp(cs2)
cs3$common.dispersion
dcc
dcs=model.matrix(~cs3$samples$group)
dcs
dcs[,-1]
dcs=dcs[,-1]
colnames(dcs)
colnames(dcs)=c('Q1','Q2','Q3','Q4')
dcs
cs4=estimateDisp(cs2,design = dcs)
cs4$common.dispersion
cs3$common.dispersion
cs3$samples$group=fct_collapse(cs3$samples$group, Q1=c('Q1','Q2'),Q3=c('Q3','Q4'))
cs3$samples
cs2$samples
cs3$samples
clq34q12c=exactTest(cs3,pair = c('Q1','Q3'))
summary(decideTests(clq34q12c))
dim(cs2$counts)
dim(cs1$counts)
dim(cs3$counts)
keep=rowSums(cpm(cs1)>1) >= 4
cs2=cs1[keep, , keep.lib.sizes=F]
dim(cs1$counts)
dim(cs2$counts)
cs2=calcNormFactors(cs2)
cs2$samples$group=fct_collapse(cs2$samples$group, Q1=c('Q1','Q2'),Q3=c('Q3','Q4')) #Convert 4 groups into 2 for classic approach
cs3=estimateDisp(cs2)
cs3$common.dispersion
cs2$samples
cs3$common.dispersion
clq34q12c=exactTest(cs3,pair = c('Q1','Q3'))
summary(decideTests(clq34q12c))
dcs
cs1$samples
keep=rowSums(cpm(cs1)>1) >= 4
cs5=cs1[keep, , keep.lib.sizes=F]
dim(cs5$counts)
cs5=calcNormFactors(cs5)
cs5$samples
cs5=estimateDisp(cs5,design = dcs)
cs5$common.dispersion
ftq34q12c=glmQLFit(cs5,design = dcs)
qlq34q12c=glmQLFTest(ftq34q12c,contrast = c(-1,-1,1,1))
summary(decideTests(qlq34q12c))



#Functional annotaions for case vs control
?goana
goglm=goana(qlcc,geneid = qlcc$genes$entrezgene)
dim(goglm)
goglm1=topGO(goglm,ont='BP',sort = 'Down',n=21647)
dim(goglm1)
head(goglm1)
goglm2=topGO(goglm,ont='BP',sort = 'Up',n=21647)
head(goglm2)
write.csv(as.data.frame(goglm1),'go_glm_down_regulated.csv')
write.csv(as.data.frame(goglm2),'go_glm_up_regulated.csv')
gocls=goana(etcc, geneid = etcc$genes$entrezgene)
gocls2=topGO(gocls,ont='BP',sort = 'Up',n=21647)
dim(gocls2)
gocls3=topGO(gocls,ont='BP',sort = 'Down',n=21647)
write.csv(as.data.frame(gocls2),'go_classic_up_regulated.csv')
write.csv(as.data.frame(gocls3),'go_classic_Down_regulated.csv')
head(gocls3)
?kegga
kglm=kegga(qlcc,geneid = qlcc$genes$entrezgene)
?topKEGG
kglm2=topKEGG(kglm,ontology=c("BP", "CC", "MF"),sort = 'Up',number = Inf)
kglm2=topKEGG(kglm,sort = 'Up',number = Inf)
dim(kglm2)
head(kglm2)
kglm3=topKEGG(kglm,sort = 'Down',number = Inf)
head(kglm3)
kcls=kegga(etcc,geneid = etcc$genes$entrezgene)
kcls2=topKEGG(kcls,sort = 'Up',number = Inf)
kcls3=topKEGG(kcls,sort = 'Down',number = Inf)
head(kcls2)
write.csv(as.data.frame(kglm2),'pathway_glm_up_regulated.csv')
write.csv(as.data.frame(kglm3),'pathway_glm_down_regulated.csv')
write.csv(as.data.frame(kcls2),'pathway_classic_up_regulated.csv')
write.csv(as.data.frame(kcls3),'pathway_classic_down_regulated.csv')

#ER, PR and HER
pa1=read.csv('/Users/mzillur/rna_analysis/Pathologies.csv')
head(pa1)
dim(gls1$samples)
head(gls1$samples)
pa2=pa1[match(rownames(gls1$samples),pa1$ID),]
pa2
dim(pa2)
head(pa2)
pa3=pa2[,c('ER','PR','HER2')]
pa3
rownames(pa2)=pa2$ID
pa2
pa2=tail(pa2,15)
pa2
pa2=pa1[match(rownames(gls1$samples),pa1$ID),]
pa3=tail(pa2,35)
pa3
rownames(pa3)=pa3$ID
pa3
pa4=pa3[,c('ER','PR','HER2')]
pa4
head(cf1)
cf2=cf1[,(2:57)]
head(cf2)
cf3=cf2[,-2]
head(cf3)
dim(cf3)
er1=DGEList(counts = cf3[,1:34],genes = cf3[,35:37])
head(er1$counts)
head(er1$samples)
head(er1$genes)
dim(er1$counts)
er2=er1
er2$samples$group=pa4$ER[match(rownames(er1$samples),rownames(pa4)),]
er2$samples$er=pa4$ER
dim(pa4)
head(pa4)
pa5=pa4[-2,]
head(pa5)
dim(pa5)
er2$samples$er=pa5$ER
er2$samples
count(er2$samples$er)
pa5
er2$samples$group=er2$samples$ER
er2$samples
keep=rowSums(cpm(er2)>1) >=3 #Remove low count genes in case and control
er3=er2[keep, , keep.lib.sizes=F]
dim(er3$counts)
dim(er2$counts)
er3=calcNormFactors(er3)
er3$samples
er3=estimateDisp(er3)
er3$common.dispersion
plotBCV(er3)
erpn=exactTest(er3,pair = c('neg','pos'))
erpc=exactTest(er3,pair = c('Control','pos'))
ernc=exactTest(er3,pair = c('Control','neg'))
summary(decideTests(erpn))
summary(decideTests(erpc))
summary(decideTests(ernc))
write.csv(as.data.frame(topTags(erpn)),'domain_er_pos_neg_classic.csv')
write.csv(as.data.frame(topTags(erpc,n=2000,sort.by = 'PValue')),'domain_er_pos_con_classic.csv')
write.csv(as.data.frame(topTags(ernc,n=200,sort.by = 'PValue')),'domain_er_neg_con_classic.csv')
write.csv(as.data.frame(topTags(ergl_pn)),'domain_er_pos_neg_glm.csv')
write.csv(as.data.frame(topTags(ergl_pc,n=2000,sort.by = 'PValue')),'domain_er_pos_con_glm.csv')
write.csv(as.data.frame(topTags(ergl_nc,n=200,sort.by = 'PValue')),'domain_er_neg_con_glm.csv')
goerpc=goana(erpc,geneid = erpc$genes$entrezgene)
goernc=goana(ernc,geneid = ernc$genes$entrezgene)
goerglpc=goana(ergl_pc,geneid = ergl_pc$genes$entrezgene)
goerglnc=goana(ergl_nc,geneid = ergl_nc$genes$entrezgene)
write.csv(as.data.frame(topGO(goerpc, ont='BP',sort = 'Up',n=21392)),'go_classic_er_pos_con_up.csv')
write.csv(as.data.frame(topGO(goerpc, ont='BP',sort = 'Down',n=21392)),'go_classic_er_pos_con_down.csv')
write.csv(as.data.frame(topGO(goerglpc, ont='BP',sort = 'Up',n=21392)),'go_glm_er_pos_con_up.csv')
write.csv(as.data.frame(topGO(goerglpc, ont='BP',sort = 'Down',n=21392)),'go_glm_er_pos_con_down.csv')
write.csv(as.data.frame(topGO(goerglnc, ont='BP',sort = 'Down',n=21392)),'go_glm_er_neg_con_down.csv')
write.csv(as.data.frame(topGO(goerglnc, ont='BP',sort = 'Up',n=21392)),'go_glm_er_neg_con_up.csv')
write.csv(as.data.frame(topGO(goernc, ont='BP',sort = 'Down',n=21392)),'go_classic_er_neg_con_down.csv')
write.csv(as.data.frame(topGO(goernc, ont='BP',sort = 'Up',n=21392)),'go_classic_er_neg_con_up.csv')
kerpc=kegga(erpc,geneid = erpc$genes$entrezgene)
kernc=kegga(ernc,geneid = ernc$genes$entrezgene)
kerglpc=kegga(ergl_pc,geneid = ergl_pc$genes$entrezgene)
kerglnc=kegga(ergl_nc,geneid = ergl_nc$genes$entrezgene)
write.csv(as.data.frame(topKEGG(kerpc,sort='Up',n = Inf)),'kegg_er_pos_con_classic_up.csv')
write.csv(as.data.frame(topKEGG(kerpc,sort='Down',n = Inf)),'kegg_er_pos_con_classic_down.csv')
write.csv(as.data.frame(topKEGG(kerglpc,sort='Up',n = Inf)),'kegg_er_pos_con_glm_up.csv')
write.csv(as.data.frame(topKEGG(kerglpc,sort='Down',n = Inf)),'kegg_er_pos_con_glm_down.csv')
write.csv(as.data.frame(topKEGG(kernc,sort='Up',n = Inf)),'kegg_er_neg_con_classic_up.csv')
write.csv(as.data.frame(topKEGG(kernc,sort='Down',n = Inf)),'kegg_er_neg_con_classic_down.csv')
write.csv(as.data.frame(topKEGG(kerglnc,sort='Up',n = Inf)),'kegg_er_neg_con_glm_up.csv')
write.csv(as.data.frame(topKEGG(kerglnc,sort='Down',n = Inf)),'kegg_er_neg_con_glm_down.csv')

moder=model.matrix(~0+er3$samples$group)
moder
colnames(moder)
colnames(moder)=c('Control','neg','pos')
moder
ergl=estimateDisp(er3,design = moder)
ergl$common.dispersion
plotBCV(ergl)
erglm_pnc=glmQLFit(ergl,design = moder)
ergl_pn=glmQLFTest(erglm_pnc,contrast = c(0,-1,1))
ergl_pc=glmQLFTest(erglm_pnc,contrast = c(-1,0,1))
ergl_nc=glmQLFTest(erglm_pnc,contrast = c(-1,1,0))
summary(decideTests(ergl_pn))
summary(decideTests(ergl_pc))
summary(decideTests(ergl_nc))
topTags(ergl_pn)
topTags(ergl_pc)
topTags(ergl_pc)

dev.off()
png('de_er.png',width = 600,height = 972)
par(mfrow=c(3,2))
plotMD(erpn)
abline(h=c(-1,1),col='green',lwd=2)
text(0.5,0.5,'DE in ER',cex=2,font = 2)
plotMD(ergl_pn)
abline(h=c(-1,1),col='green',lwd=2)
plotMD(erpc)
abline(h=c(-1,1),col='green',lwd=2)
plotMD(ergl_pc)
abline(h=c(-1,1),col='green',lwd=2)
plotMD(ernc)
abline(h=c(-1,1),col='green',lwd=2)
plotMD(ergl_nc)
abline(h=c(-1,1),col='green',lwd=2)
dev.off()

#PR related
pr1$samples
pr1$samples$group=pr1$samples$PR
pr1$samples
keep=rowSums(cpm(pr1)>1) >=3 #Remove low count genes in case and control
pr2=pr1[keep, , keep.lib.sizes=F]
dim(pr2$counts)
dim(pr1$counts)
pr3=calcNormFactors(pr2)
pr3$samples
pr3=estimateDisp(pr3)
pr3$common.dispersion
plotBCV(pr3)
prpn=exactTest(pr3,pair = c('neg','pos'))
prpc=exactTest(pr3,pair = c('Control','pos'))
prnc=exactTest(pr3,pair = c('Control','neg'))
summary(decideTests(prpn))
summary(decideTests(prpc))
summary(decideTests(prnc))
topTags(prpn)
topTags(prpc)
topTags(prnc)

modpr=model.matrix(~0+pr3$samples$group)
modpr
colnames(modpr)
colnames(modpr)=c('Control','neg','pos')
modpr
prgl=estimateDisp(pr3,design = modpr)
prgl$common.dispersion
plotBCV(prgl)
prglm_pnc=glmQLFit(prgl,design = modpr)
prgl_pn=glmQLFTest(prglm_pnc,contrast = c(0,-1,1))
prgl_pc=glmQLFTest(prglm_pnc,contrast = c(-1,0,1))
prgl_nc=glmQLFTest(prglm_pnc,contrast = c(-1,1,0))
summary(decideTests(prgl_pn))
summary(decideTests(prgl_pc))
summary(decideTests(prgl_nc))
topTags(prgl_pn)
topTags(prgl_pc)
topTags(prgl_nc)

write.csv(as.data.frame(topTags(prpn)),'domain_pr_pos_neg_classic.csv')
write.csv(as.data.frame(topTags(prpc,n=2500,sort.by = 'PValue')),'domain_pr_pos_con_classic.csv')
write.csv(as.data.frame(topTags(prnc,n=200,sort.by = 'PValue')),'domain_pr_neg_con_classic.csv')
write.csv(as.data.frame(topTags(prgl_pn)),'domain_pr_pos_neg_glm.csv')
write.csv(as.data.frame(topTags(prgl_pc,n=2500,sort.by = 'PValue')),'domain_pr_pos_con_glm.csv')
write.csv(as.data.frame(topTags(prgl_nc,n=200,sort.by = 'PValue')),'domain_pr_neg_con_glm.csv')

goprpc=goana(prpc,geneid = prpc$genes$entrezgene)
goprnc=goana(prnc,geneid = prnc$genes$entrezgene)
goprglpc=goana(prgl_pc,geneid = prgl_pc$genes$entrezgene)
goprglnc=goana(prgl_nc,geneid = prgl_nc$genes$entrezgene)
write.csv(as.data.frame(topGO(goprpc, ont='BP',sort = 'Up',n=21392)),'go_classic_pr_pos_con_up.csv')
write.csv(as.data.frame(topGO(goprpc, ont='BP',sort = 'Down',n=21392)),'go_classic_pr_pos_con_down.csv')
write.csv(as.data.frame(topGO(goprglpc, ont='BP',sort = 'Up',n=21392)),'go_glm_pr_pos_con_up.csv')
write.csv(as.data.frame(topGO(goprglpc, ont='BP',sort = 'Down',n=21392)),'go_glm_pr_pos_con_down.csv')
write.csv(as.data.frame(topGO(goprglnc, ont='BP',sort = 'Down',n=21392)),'go_glm_pr_neg_con_down.csv')
write.csv(as.data.frame(topGO(goprglnc, ont='BP',sort = 'Up',n=21392)),'go_glm_pr_neg_con_up.csv')
write.csv(as.data.frame(topGO(goprnc, ont='BP',sort = 'Down',n=21392)),'go_classic_pr_neg_con_down.csv')
write.csv(as.data.frame(topGO(goprnc, ont='BP',sort = 'Up',n=21392)),'go_classic_pr_neg_con_up.csv')

kprpc=kegga(prpc,geneid = prpc$genes$entrezgene)
kprnc=kegga(prnc,geneid = prnc$genes$entrezgene)
kprglpc=kegga(prgl_pc,geneid = prgl_pc$genes$entrezgene)
kprglnc=kegga(prgl_nc,geneid = prgl_nc$genes$entrezgene)
write.csv(as.data.frame(topKEGG(kprpc,sort='Up',n = Inf)),'kegg_pr_pos_con_classic_up.csv')
write.csv(as.data.frame(topKEGG(kprpc,sort='Down',n = Inf)),'kegg_pr_pos_con_classic_down.csv')
write.csv(as.data.frame(topKEGG(kprglpc,sort='Up',n = Inf)),'kegg_pr_pos_con_glm_up.csv')
write.csv(as.data.frame(topKEGG(kprglpc,sort='Down',n = Inf)),'kegg_pr_pos_con_glm_down.csv')
write.csv(as.data.frame(topKEGG(kprnc,sort='Up',n = Inf)),'kegg_pr_neg_con_classic_up.csv')
write.csv(as.data.frame(topKEGG(kprnc,sort='Down',n = Inf)),'kegg_pr_neg_con_classic_down.csv')
write.csv(as.data.frame(topKEGG(kprglnc,sort='Up',n = Inf)),'kegg_pr_neg_con_glm_up.csv')
write.csv(as.data.frame(topKEGG(kprglnc,sort='Down',n = Inf)),'kegg_pr_neg_con_glm_down.csv')

png('de_pr.png',width = 600,height = 972)
par(mfrow=c(3,2))
plotMD(prpn)
abline(h=c(-1,1),col='green',lwd=2)
plotMD(prgl_pn)
abline(h=c(-1,1),col='green',lwd=2)
plotMD(prpc)
abline(h=c(-1,1),col='green',lwd=2)
plotMD(prgl_pc)
abline(h=c(-1,1),col='green',lwd=2)
plotMD(prnc)
abline(h=c(-1,1),col='green',lwd=2)
plotMD(prgl_nc)
abline(h=c(-1,1),col='green',lwd=2)
dev.off()

#HER2
hr1$samples
keep=rowSums(cpm(hr1)>1) >=3 #Remove low count genes in case and control
her2=hr1[keep, , keep.lib.sizes=F]
dim(hr1$counts)
dim(her2$counts)
her3=calcNormFactors(her2)
her3$samples
her3=estimateDisp(her3)
her3$common.dispersion
plotBCV(her3)
herpn=exactTest(her3,pair = c('neg','pos'))
herpc=exactTest(her3,pair = c('Control','pos'))
hernc=exactTest(her3,pair = c('Control','neg'))
summary(decideTests(herpn))
summary(decideTests(herpc))
summary(decideTests(hernc))
topTags(herpn)
topTags(herpc)
topTags(hernc)

modher=model.matrix(~0+her3$samples$group)
modher
colnames(modher)
colnames(modher)=c('Control','neg','pos')
modher
hergl=estimateDisp(her3,design = modher)
hergl$common.dispersion
plotBCV(hergl)
herglm_pnc=glmQLFit(hergl,design = modher)
hergl_pn=glmQLFTest(herglm_pnc,contrast = c(0,-1,1))
hergl_pc=glmQLFTest(herglm_pnc,contrast = c(-1,0,1))
hergl_nc=glmQLFTest(herglm_pnc,contrast = c(-1,1,0))
summary(decideTests(hergl_pn))
summary(decideTests(hergl_pc))
summary(decideTests(hergl_nc))
topTags(hergl_pn)
topTags(hergl_pc)
topTags(hergl_nc)

write.csv(as.data.frame(topTags(herpn)),'domain_her_pos_neg_classic.csv')
write.csv(as.data.frame(topTags(herpc,n=2500,sort.by = 'PValue')),'domain_her_pos_con_classic.csv')
write.csv(as.data.frame(topTags(hernc,n=200,sort.by = 'PValue')),'domain_her_neg_con_classic.csv')
write.csv(as.data.frame(topTags(hergl_pn)),'domain_her_pos_neg_glm.csv')
write.csv(as.data.frame(topTags(hergl_pc,n=2500,sort.by = 'PValue')),'domain_her_pos_con_glm.csv')
write.csv(as.data.frame(topTags(hergl_nc,n=200,sort.by = 'PValue')),'domain_her_neg_con_glm.csv')

gohernc=goana(hernc,geneid = hernc$genes$entrezgene)
goherglnc=goana(hergl_nc,geneid = hergl_nc$genes$entrezgene)
write.csv(as.data.frame(topGO(goherglnc, ont='BP',sort = 'Down',n=21392)),'go_glm_her_neg_con_down.csv')
write.csv(as.data.frame(topGO(goherglnc, ont='BP',sort = 'Up',n=21392)),'go_glm_her_neg_con_up.csv')
write.csv(as.data.frame(topGO(gohernc, ont='BP',sort = 'Down',n=21392)),'go_classic_her_neg_con_down.csv')
write.csv(as.data.frame(topGO(gohernc, ont='BP',sort = 'Up',n=21392)),'go_classic_her_neg_con_up.csv')

khernc=kegga(hernc,geneid = hernc$genes$entrezgene)
kherglnc=kegga(hergl_nc,geneid = hergl_nc$genes$entrezgene)
write.csv(as.data.frame(topKEGG(khernc,sort='Up',n = Inf)),'kegg_her_neg_con_classic_up.csv')
write.csv(as.data.frame(topKEGG(khernc,sort='Down',n = Inf)),'kegg_her_neg_con_classic_down.csv')
write.csv(as.data.frame(topKEGG(kherglnc,sort='Up',n = Inf)),'kegg_her_neg_con_glm_up.csv')
write.csv(as.data.frame(topKEGG(kherglnc,sort='Down',n = Inf)),'kegg_her_neg_con_glm_down.csv')

png('de_her.png',width = 600,height = 972)
par(mfrow=c(3,2))
plotMD(herpn)
abline(h=c(-1,1),col='green',lwd=2)
plotMD(hergl_pn)
abline(h=c(-1,1),col='green',lwd=2)
plotMD(herpc)
abline(h=c(-1,1),col='green',lwd=2)
plotMD(hergl_pc)
abline(h=c(-1,1),col='green',lwd=2)
plotMD(hernc)
abline(h=c(-1,1),col='green',lwd=2)
plotMD(hergl_nc)
abline(h=c(-1,1),col='green',lwd=2)
dev.off()

#AFR
sm3=data.frame(lapply(sm2, as.character), stringsAsFactors=FALSE) #Convert factor into character
afr1=DGEList(counts = cf1[,sm3$Sample_ID],genes = cf1[,(52:57)])
afr1$samples$AFR=sm3$afrq[match(rownames(afr1$samples),sm3$Sample_ID)]
afr1$samples
table(afr1$samples$AFR)
afr1$samples$group=afr1$samples$AFR
keep=rowSums(cpm(afr1)>1) >=4 #Remove low count genes in case and control
afr2=afr1[keep, , keep.lib.sizes=F]
dim(afr2$counts)
dim(afr1$counts)
afr3=calcNormFactors(afr2)
afr3$samples
afr3=estimateDisp(afr3)
afr3$common.dispersion
plotBCV(afr3)
afrq4q1=exactTest(afr3,pair = c('Q1','Q4'))
afr4=afr3
afr4$samples$group=fct_collapse(afr4$samples$group, Q1=c('Q1','Q2'),Q3=c('Q3','Q4'))
afr4$samples
afr4=calcNormFactors(afr4)
afr4=estimateDisp(afr4)
afr4$common.dispersion
afrq34q12=exactTest(afr4,pair = c('Q1','Q3'))
summary(decideTests(afrq4q1))
summary(decideTests(afrq34q12))
topTags(afrq4q1)
topTags(afrq34q12)

modafr=model.matrix(~0+afr3$samples$group)
modafr
colnames(modafr)
colnames(modafr)=c('Q1','Q2','Q3','Q4')
modafr
afrgl=estimateDisp(afr3,design = modafr)
afrgl$common.dispersion
plotBCV(afrgl)
afrglm=glmQLFit(afrgl,design = modafr)
afrglmq4q1=glmQLFTest(afrglm,contrast = c(-1,0,0,1))
afrglmq34q12=glmQLFTest(afrglm,contrast = c(-1,-1,1,1))
summary(decideTests(afrglmq4q1))
summary(decideTests(afrglmq34q12))
topTags(afrglmq4q1)
topTags(afrglmq34q12)

write.csv(as.data.frame(topTags(afrq4q1,n=500,sort.by = 'PValue')),'domain_afr_q4_q1_classic.csv')
write.csv(as.data.frame(topTags(afrq34q12,n=500,sort.by = 'PValue')),'domain_afr_q34_q12_classic.csv')
write.csv(as.data.frame(topTags(afrglmq4q1,n=500,sort.by = 'PValue')),'domain_afr_q4_q1_glm.csv')
write.csv(as.data.frame(topTags(afrglmq34q12,n=500,sort.by = 'PValue')),'domain_afr_q34_q12_glm.csv')

#EUR
eur1=DGEList(counts = cf1[,sm3$Sample_ID],genes = cf1[,(52:57)])
eur1$samples$EUR=sm3$eurq[match(rownames(eur1$samples),sm3$Sample_ID)]
eur1$samples
table(eur1$samples$EUR)
eur1$samples$group=eur1$samples$EUR
keep=rowSums(cpm(eur1)>1) >=4 #Remove low count genes in case and control
eur2=eur1[keep, , keep.lib.sizes=F]
dim(eur2$counts)
dim(eur1$counts)
eur3=calcNormFactors(eur2)
eur3$samples
eur3=estimateDisp(eur3)
eur3$common.dispersion
plotBCV(eur3)
eurq4q1=exactTest(eur3,pair = c('Q1','Q4'))
eur4=eur3
eur4$samples$group=fct_collapse(eur4$samples$group, Q1=c('Q1','Q2'),Q3=c('Q3','Q4'))
eur4$samples
eur4=calcNormFactors(eur4)
eur4=estimateDisp(eur4)
eur4$common.dispersion
eurq34q12=exactTest(eur4,pair = c('Q1','Q3'))
summary(decideTests(eurq4q1))
summary(decideTests(eurq34q12))
topTags(eurq4q1)
topTags(eurq34q12)

modeur=model.matrix(~0+eur3$samples$group)
modeur
colnames(modeur)
colnames(modeur)=c('Q1','Q2','Q3','Q4')
modeur
eurgl=estimateDisp(eur3,design = modeur)
eurgl$common.dispersion
plotBCV(eurgl)
eurglm=glmQLFit(eurgl,design = modeur)
eurglmq4q1=glmQLFTest(eurglm,contrast = c(-1,0,0,1))
eurglmq34q12=glmQLFTest(eurglm,contrast = c(-1,-1,1,1))
summary(decideTests(eurglmq4q1))
summary(decideTests(eurglmq34q12))
topTags(eurglmq4q1)
topTags(eurglmq34q12)

write.csv(as.data.frame(topTags(eurq4q1,n=500,sort.by = 'PValue')),'domain_eur_q4_q1_classic.csv')
write.csv(as.data.frame(topTags(eurq34q12,n=500,sort.by = 'PValue')),'domain_eur_q34_q12_classic.csv')
write.csv(as.data.frame(topTags(eurglmq4q1,n=500,sort.by = 'PValue')),'domain_eur_q4_q1_glm.csv')
write.csv(as.data.frame(topTags(eurglmq34q12,n=500,sort.by = 'PValue')),'domain_eur_q34_q12_glm.csv')


#NAT
nat1=DGEList(counts = cf1[,sm3$Sample_ID],genes = cf1[,(52:57)])
nat1$samples$NAT=sm3$natq[match(rownames(nat1$samples),sm3$Sample_ID)]
nat1$samples
table(nat1$samples$NAT)
nat1$samples$group=nat1$samples$NAT
keep=rowSums(cpm(nat1)>1) >=4 #Remove low count genes in case and control
nat2=nat1[keep, , keep.lib.sizes=F]
dim(nat2$counts)
dim(nat1$counts)
nat3=calcNormFactors(nat2)
nat3$samples
nat3=estimateDisp(nat3)
nat3$common.dispersion
plotBCV(nat3)
natq4q1=exactTest(nat3,pair = c('Q1','Q4'))
nat4=nat3
nat4$samples$group=fct_collapse(nat4$samples$group, Q1=c('Q1','Q2'),Q3=c('Q3','Q4'))
nat4$samples
nat4=calcNormFactors(nat4)
nat4=estimateDisp(nat4)
nat4$common.dispersion
natq34q12=exactTest(nat4,pair = c('Q1','Q3'))
summary(decideTests(natq4q1))
summary(decideTests(natq34q12))
topTags(natq4q1)
topTags(natq34q12)

modnat=model.matrix(~0+nat3$samples$group)
modnat
colnames(modnat)
colnames(modnat)=c('Q1','Q2','Q3','Q4')
modnat
natgl=estimateDisp(nat3,design = modnat)
natgl$common.dispersion
plotBCV(natgl)
natglm=glmQLFit(natgl,design = modnat)
natglmq4q1=glmQLFTest(natglm,contrast = c(-1,0,0,1))
natglmq34q12=glmQLFTest(natglm,contrast = c(-1,-1,1,1))
summary(decideTests(natglmq4q1))
summary(decideTests(natglmq34q12))
topTags(natglmq4q1)
topTags(natglmq34q12)

write.csv(as.data.frame(topTags(natq4q1,n=500,sort.by = 'PValue')),'domain_nat_q4_q1_classic.csv')
write.csv(as.data.frame(topTags(natq34q12,n=500,sort.by = 'PValue')),'domain_nat_q34_q12_classic.csv')
write.csv(as.data.frame(topTags(natglmq4q1,n=500,sort.by = 'PValue')),'domain_nat_q4_q1_glm.csv')
write.csv(as.data.frame(topTags(natglmq34q12,n=500,sort.by = 'PValue')),'domain_nat_q34_q12_glm.csv')

#DNA Repair System
dnar=read.csv('../dna_repair_hk_classic.csv',row.names = 1)
dna1=dnar[,1:2]
dna2=dna1
dna1$Case_vs_Control=etcc$table$logFC[match(rownames(dna1),rownames(etcc$table))]
head(dna1)
dna1$DRC_Q2_vs_Q1=clq2q1$table$logFC[match(rownames(dna1),rownames(clq2q1$table))]
head(dna1)
dna1=dna1[,-(1:2)]
head(dna1)
dna1$DRC_Q3_vs_Q1=clq3q1$table$logFC[match(rownames(dna1),rownames(clq3q1$table))]
head(dna1)
dna1$DRC_Q4_vs_Q1=clq4q1$table$logFC[match(rownames(dna1),rownames(clq4q1$table))]
head(dna1)
dna1$DRC_Q3_vs_Q2=clq3q2$table$logFC[match(rownames(dna1),rownames(clq3q2$table))]
head(dna1)
dna1$DRC_Q4_vs_Q2=clq4q2$table$logFC[match(rownames(dna1),rownames(clq4q2$table))]
head(dna1)
dna1$DRC_Q4_vs_Q3=clq4q3$table$logFC[match(rownames(dna1),rownames(clq4q3$table))]
head(dna1)
dna1$DRC_Q34_vs_Q12=clq34q12$table$logFC[match(rownames(dna1),rownames(clq34q12$table))]
head(dna1)
dna1$DRC_Cases_Q34_vs_Q12=clq34q12c$table$logFC[match(rownames(dna1),rownames(clq34q12c$table))]
head(dna1)
dna1$ER_pos_vs_neg=erpn$table$logFC[match(rownames(dna1),rownames(erpn$table))]
head(dna1)
dna1$ER_pos_vs_Control=erpc$table$logFC[match(rownames(dna1),rownames(erpc$table))]
head(dna1)
dna1$ER_neg_vs_Control=ernc$table$logFC[match(rownames(dna1),rownames(ernc$table))]
head(dna1)
dna1$PR_pos_vs_neg=prpn$table$logFC[match(rownames(dna1),rownames(prpn$table))]
head(dna1)
dna1$PR_pos_vs_Control=prpc$table$logFC[match(rownames(dna1),rownames(prpc$table))]
head(dna1)
dna1$PR_neg_vs_Conrol=prnc$table$logFC[match(rownames(dna1),rownames(prnc$table))]
head(dna1)
dna1$HER2_pos_vs_neg=herpn$table$logFC[match(rownames(dna1),rownames(herpn$table))]
head(dna1)
dna1$HER2_pos_vs_Control=herpc$table$logFC[match(rownames(dna1),rownames(herpc$table))]
head(dna1)
dna1$HER2_neg_vs_Control=hernc$table$logFC[match(rownames(dna1),rownames(hernc$table))]
head(dna1)
dna1$AFR_Q4_vs_Q1=afrq4q1$table$logFC[match(rownames(dna1),rownames(afrq4q1$table))]
head(dna1)
dna1$AFR_Q34_vs_Q12=afrq34q12$table$logFC[match(rownames(dna1),rownames(afrq34q12$table))]
head(dna1)
dna1$EUR_Q4_vs_Q1=eurq4q1$table$logFC[match(rownames(dna1),rownames(eurq4q1$table))]
head(dna1)
dna1$EUR_Q34_vs_Q12=eurq34q12$table$logFC[match(rownames(dna1),rownames(eurq34q12$table))]
head(dna1)
dna1$NAT_Q4_vs_Q1=natq4q1$table$logFC[match(rownames(dna1),rownames(natq4q1$table))]
head(dna1)
dna1$NAT_Q34_vs_Q12=natq34q12$table$logFC[match(rownames(dna1),rownames(natq34q12$table))]
head(dna1)
dna3=as.matrix(dna1)
dim(dna3)
head(dna3)
dna4=dna3[complete.cases(dna3),]
dim(dna4)
png('heatmap_classic_2.png',height = 3560,width = 2200)
pheatmap(dna4,color = rainbow(3),fontsize_row = 18,fontsize_col = 22)
dev.off()

dna2$Case_vs_Control=qlcc$table$logFC[match(rownames(dna2),rownames(qlcc$table))]
head(dna2)
dna2$DRC_Q2_vs_Q1=qlq2q1$table$logFC[match(rownames(dna2),rownames(qlq2q1$table))]
head(dna2)
dna2=dna2[,-(1:2)]
head(dna2)
dna2$DRC_Q3_vs_Q1=qlq3q1$table$logFC[match(rownames(dna2),rownames(qlq3q1$table))]
head(dna2)
dna2$DRC_Q4_vs_Q1=qlq4q1$table$logFC[match(rownames(dna2),rownames(qlq4q1$table))]
head(dna2)
dna2$DRC_Q3_vs_Q2=qlq3q2$table$logFC[match(rownames(dna2),rownames(qlq3q2$table))]
head(dna2)
dna2$DRC_Q4_vs_Q2=qlq4q2$table$logFC[match(rownames(dna2),rownames(qlq4q2$table))]
head(dna2)
dna2$DRC_Q4_vs_Q3=qlq4q3$table$logFC[match(rownames(dna2),rownames(qlq4q3$table))]
head(dna2)
dna2$DRC_Q34_vs_Q12=qlq34q12$table$logFC[match(rownames(dna2),rownames(qlq34q12$table))]
head(dna2)
dna2$DRC_Cases_Q34_vs_Q12=qlq34q12c$table$logFC[match(rownames(dna2),rownames(qlq34q12c$table))]
head(dna2)
dna2$ER_pos_vs_neg=ergl_pn$table$logFC[match(rownames(dna2),rownames(ergl_pn$table))]
head(dna2)
dna2$ER_pos_vs_Control=ergl_pc$table$logFC[match(rownames(dna2),rownames(ergl_pc$table))]
head(dna2)
dna2$ER_neg_vs_Control=ergl_nc$table$logFC[match(rownames(dna2),rownames(ergl_nc$table))]
head(dna2)
dna2$PR_pos_vs_neg=prgl_pn$table$logFC[match(rownames(dna2),rownames(prgl_pn$table))]
head(dna2)
dna2$PR_pos_vs_Control=prgl_pc$table$logFC[match(rownames(dna2),rownames(prgl_pc$table))]
head(dna2)
dna2$PR_neg_vs_Conrol=prgl_nc$table$logFC[match(rownames(dna2),rownames(prgl_nc$table))]
head(dna2)
dna2$HER2_pos_vs_neg=hergl_pn$table$logFC[match(rownames(dna2),rownames(hergl_pn$table))]
head(dna2)
dna2$HER2_pos_vs_Control=hergl_pc$table$logFC[match(rownames(dna2),rownames(hergl_pc$table))]
head(dna2)
dna2$HER2_neg_vs_Control=hergl_nc$table$logFC[match(rownames(dna2),rownames(hergl_nc$table))]
head(dna2)
dna2$AFR_Q4_vs_Q1=afrglmq4q1$table$logFC[match(rownames(dna2),rownames(afrglmq4q1$table))]
head(dna2)
dna2$AFR_Q34_vs_Q12=afrglmq34q12$table$logFC[match(rownames(dna2),rownames(afrglmq34q12$table))]
head(dna2)
dna2$EUR_Q4_vs_Q1=eurglmq4q1$table$logFC[match(rownames(dna2),rownames(eurglmq4q1$table))]
head(dna2)
dna2$EUR_Q34_vs_Q12=eurglmq34q12$table$logFC[match(rownames(dna2),rownames(eurglmq34q12$table))]
head(dna2)
dna2$NAT_Q4_vs_Q1=natglmq4q1$table$logFC[match(rownames(dna2),rownames(natglmq4q1$table))]
head(dna2)
dna2$NAT_Q34_vs_Q12=natglmq34q12$table$logFC[match(rownames(dna2),rownames(natglmq34q12$table))]
head(dna2)
dna5=dna2[complete.cases(dna2),]
dim(dna5)
png('heatmap_glm_1.png',height = 3560,width = 2200)
pheatmap(dna5,color = rainbow(3),fontsize_row = 18,fontsize_col = 22)
dev.off()
write.csv(as.data.frame(dna5),'dna_repair_glm.csv')
write.csv(as.data.frame(dna4),'dna_repair_classic.csv')
dna4$InterPro=natglmq34q12$genes$interpro_description[match(rownames(dna4),rownames(natglmq34q12$genes))]
dna5$InterPro=natglmq34q12$genes$interpro_description[match(rownames(dna5),rownames(natglmq34q12$genes))]
head(dna4)
write.csv(as.data.frame(dna5),'dna_repair_interpro.csv')
