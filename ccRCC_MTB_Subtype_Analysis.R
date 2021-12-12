
####Bach effect removing 

library(cluster)
library(oompaBase)
library(sva)
batch <- data.frame(batch = rep(c("CPTAC","Tokyo","TCGA"), times = c(ncol(CPTAC),ncol(Tokyo),ncol(TCGA))))
modcombat = model.matrix(~1, data = batch)
combined.expr.combat <- as.data.frame(ComBat(dat=as.matrix(combined.expr), batch=batch$batch, mod=modcombat))

#####Quantifying 100 KEGG metabolic pathways
library(IOBR)
sig_metabolism<-calculate_sig_score(pdata    = NULL,
                                    eset            = combined.expr.combat,
                                    signature       = signature_metabolism,
                                    method          = "ssgsea",
                                    mini_gene_count = 2)

#####Consensus Clustering for distinct ccRCC subtypes
library(ConsensusClusterPlus)
indata <- t(scale(sig_metabolism))
subtype <- ConsensusClusterPlus(d = as.matrix(indata),
                                maxK = 4, 
                                pItem = 0.8, 
                                pFeature = 1, 
                                reps = 1000, 
                                clusterAlg = "pam", 
                                innerLinkage = "ward.D", 
                                finalLinkage = "ward.D", 
                                distance = "euclidean", 
                                seed = 123456,
                                plot = "pdf", 
                                writeTable = TRUE,
                                title = "ConsensusCluster-Metablism") 

######Calculating cluster specific pathways
library(clusterProfiler)
library(GSVA)
library(pheatmap)
library(gplots)

colnames(combined.expr.combat) <- rownames(subt)
expr <- combined.expr.combat
subt <- ClusterSub

colnames(subt) <- 'TCGA_Subtype'
rownames(subt)  <- rownames(design)
table(subt$TCGA_Subtype)

display.progress = function ( index, totalN, breakN=20) {
  if ( index %% ceiling(totalN/breakN)  ==0  ) {
    cat(paste(round(index*100/totalN), "% ", sep=""))
  }
} 

n.sub <- length(table(subt$TCGA_Subtype)) 
n.sub.label <- unique(subt$TCGA_Subtype) 
treat_list <- ctrl_list <- degs.list <- list() 
for (i in 1:n.sub) {
  cat(paste0(n.sub.label[i], " vs. Others starts!\n"))
  treat_list[[i]] <- rownames(subt)[which(subt$TCGA_Subtype == n.sub.label[i])] 
  ctrl_list[[i]] <- rownames(subt)[-which(subt$TCGA_Subtype == n.sub.label[i])]
  
  meanA <- meanB <- p <- fc <- lgfc <- c() 
  for (k in 1:nrow(expr)) {
    display.progress(index = k,totalN = nrow(expr))
    a <- as.numeric(expr[k,treat_list[[i]]])
    b <- as.numeric(expr[k,ctrl_list[[i]]])
    p <- c(p,t.test(a,b,na.rm=T)$p.value) 
    meanA <- c(meanA,mean(a)) 
    meanB <- c(meanB,mean(b)) 
    fc <- c(fc,mean(a)/mean(b)) 
    lgfc <- c(lgfc,log2(mean(a)/mean(b))) 
  }
  fdr <- p.adjust(p,method = "fdr") 
  degs <- data.frame(mean_treat=meanA,
                     mean_ctrl=meanB,
                     FoldChange=fc,
                     log2FoldChange=lgfc,
                     pvalue=p,
                     padj=fdr,
                     row.names = rownames(expr),
                     stringsAsFactors = F)
  
  write.table(degs,paste0(n.sub.label[[i]],"_degs.txt"),sep = "\t",row.names = T,col.names = NA,quote = F)
  
  degs.list[[n.sub.label[i]]] <- as.data.frame(na.omit(degs))
  
  cat("\n")
}
###Up regulated
subtype_specific_gsea <- function(msigdb=NULL,n.top=10,mode=c("up","down"),degs.list=NULL,subtype.label=NULL,nPerm.gsea=1000,minGSSize.gsea=10,maxGSSize.gsea=500,pvalueCutoff.gsea=1){
  
  MSigDB <- read.gmt(msigdb)
  GSEA.list <- top.gs <- list() 
  
  if(!is.element(mode, c("up", "dn"))) { stop("mode must be up or dn!\n") }
  
  for (i in 1:n.sub) {
    degs <- degs.list[[n.sub.label[i]]]
    geneList <- degs$log2FoldChange; names(geneList) <- rownames(degs)
    geneList <- sort(geneList,decreasing = T) 
    
  
    cat(paste0("GSEA for ",subtype.label[i]," starts!\n"))
    GSEA.list[[subtype.label[i]]] <- GSEA(geneList = geneList,
                                          TERM2GENE=MSigDB,
                                          nPerm = nPerm.gsea,
                                          minGSSize = minGSSize.gsea,
                                          maxGSSize = maxGSSize.gsea,
                                          seed = T,
                                          verbose = F,
                                          pvalueCutoff = pvalueCutoff.gsea) 
    
    GSEA.dat <- as.data.frame(GSEA.list[[subtype.label[i]]])
    
    if(mode == "up") {
      GSEA.dat <- GSEA.dat[order(GSEA.dat$NES,decreasing = T),] 
    } else {
      GSEA.dat <- GSEA.dat[order(GSEA.dat$NES,decreasing = F),] 
    }
    
    
    write.table(GSEA.dat,paste0(subtype.label[[i]],"_degs_",mode,"_gsea.txt"),sep = "\t",row.names = T,col.names = NA,quote = F)
    
    
    top.gs[[subtype.label[i]]] <- rownames(GSEA.dat)[1:n.top] 
  }
  
  gs <- list()
  for (i in as.character(unlist(top.gs))) {
    gs[[i]] <- MSigDB[which(MSigDB[,1] %in% i),"gene"]
  }
  
  return(list(mode=mode,top.gs=top.gs,gs=gs))
}

msigdfFile = "c5.all.v6.2.symbols.gmt"
n.top = 10
mode = "up" 
gs.up <- subtype_specific_gsea(msigdb = msigdfFile,
                               n.top = n.top,
                               degs.list = degs.list,
                               subtype.label = n.sub.label,
                               mode = mode)

gsva_gs.up <- gsva(as.matrix(expr), gs.up$gs, method="gsva") 
dim(gsva_gs.up)

gsva_gs.up_mean <- data.frame(row.names = rownames(gsva_gs.up)) 
for (i in n.sub.label) {
  gsva_gs.up_mean <- cbind.data.frame(gsva_gs.up_mean,
                                      data.frame(rowMeans(gsva_gs.up[,rownames(subt)[which(subt$TCGA_Subtype == i)]])))
}
colnames(gsva_gs.up_mean) <- n.sub.label



### Calculating DEGs correlated to the corresponding MTB Subtypes
library(clusterProfiler)
library(gplots)
library(limma)
library(Glimma)
library(edgeR)

ClusterSub<- read.csv(file = 'ConsensusCluster-Metablism.k=3.consensusClass.csv',
                      row.names = 1)
ClusterSub[which(ClusterSub$class==1),] <- 'C1'
ClusterSub[which(ClusterSub$class==2),] <- 'C2'
ClusterSub[which(ClusterSub$class==3),] <- 'C3'
exprSet <- combined.expr.combat[,rownames(ClusterSub)]

group <- factor(ClusterSub$class)
design <- model.matrix(~0+factor(group))
colnames(design)=levels(factor(group))
rownames(design)= colnames(exprSet)

fit <- lmFit(exprSet,design)

contrast.matrix <- makeContrasts(C1vsC2 = C1-C2, 
                                 C1vsC3 = C1-C3, 
                                 C2vsC3 = C2-C3,levels=design)
fit2 <- contrasts.fit(fit,contrast.matrix)

fit2 <- eBayes(fit2)



#####Construction of MTB signatures A and B
library(Boruta)
library(ggplot2)
library(ggpubr)
library(pheatmap)
geneclust <- subtype[[3]]
samorder <- sort(geneclust$consensusClass)

outTab <- NULL
for (i in rownames(indata)) {
  tmp <- as.numeric(indata[i,names(samorder)])
  cor.res <- cor.test(tmp, as.numeric(samorder), method = "pearson")
  outTab <- rbind.data.frame(outTab,
                             data.frame(gene = i,
                                        r = cor.res$estimate,
                                        p = cor.res$p.value,
                                        stringsAsFactors = F),
                             stringsAsFactors = F)
}

# Correlation
outTab$direct <- ifelse(outTab$r > 0, "A","B")
outTab <- outTab[order(outTab$r, decreasing = T),]
table(outTab$direct) #
write.table(outTab,"ouput_MTB signature Gene.txt",sep = "\t", row.names = F, col.names = T, quote = F)



set.seed(20201024)
dat.boruta <- as.data.frame(t(indata[outTab$gene,rownames(annCol)]))
borutafit <- Boruta(x = as.matrix(dat.boruta), 
                    y = as.factor(annCol$`Gene cluster`), # multiclassification
                    doTrace = 2,
                    maxRuns = 100,
                    ntree = 500)
boruta_fea <- attStats(borutafit)
boruta_fea <- rownames(boruta_fea[which(boruta_fea$decision == "Confirmed"),])
boruta.all <- outTab[which(outTab$gene %in% boruta_fea),]
table(boruta.all$direct)

boruta.A <- boruta.all[1:119,]
save(boruta.A,file = 'boruta.A.rda')
boruta.B <- boruta.all[120:206,]
save(boruta.B,file = 'boruta.B.rda')
annRow <- data.frame("MTB signature gene" = rep(c("A","B"), c(nrow(boruta.A), nrow(boruta.B))),
                     row.names = c(boruta.A$gene,boruta.B$gene),
                     check.names = F,
                     stringsAsFactors = F)



####Survival analysis with optimal cutoff value
library(survminer)
library(survival)
res.cut <- surv_cutpoint(annCol, time = "Time", 
                         event = "Status", 
                         variables ="MTB_Score", 
                         minprop = 0.3) 
res.cat <- surv_categorize(res.cut)
table(res.cat$MTB_Score)

fitd <- survdiff(Surv(Time, Status) ~ MTB_Score ,
                 data      = res.cat,
                 na.action = na.exclude)

p.val <- 1 - pchisq(fitd$chisq, length(fitd$n) - 1)
fit <- survfit(Surv(Time, Status) ~ MTB_Score ,
               data      = res.cat,
               type      = "kaplan-meier",
               error     = "greenwood",
               conf.type = "plain",
               na.action = na.exclude)

names(fit$strata) <- gsub("MTB_Score=", "", names(fit$strata))
# kaplan-meier curve

blue <- "#2874C5"
yellow <- "#EABF00"
green <- "#008B8A"
red <- "#E21F26"


p <- ggsurvplot(fit               = fit,
                conf.int          = FALSE,
                risk.table        = TRUE,
                risk.table.col    = "strata",
                palette           = c("#0099B4CC","#AD002ACC" ),
                data              = res.cat,
                size              = 1,
                xlim              = c(0,100),
                break.time.by     = 20,
                legend.title      = "MTB_Score",
                pval              = FALSE,
                surv.median.line  = "hv",
                xlab              = "Time (month)",
                ylab              = "Survival probability",
                risk.table.y.text = FALSE)
p.lab <- paste0("Log rank test P", 
                ifelse(p.val < 0.001, " < 0.001",
                       paste0(" = ",round(p.val, 3))))

p$plot <- p$plot + annotate("text", 
                            x = 50, 
                            y = 0.55,
                            hjust = 0,
                            fontface = 4,
                            label = p.lab)








