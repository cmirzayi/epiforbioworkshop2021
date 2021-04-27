#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install("curatedMetagenomicData")


library(curatedMetagenomicData)
library(dagitty)
library(DESeq2)
library(phyloseq)
library(ggdag)

#Step 1
v_dag <- dagitty( "dag {
  crc -> microbiome
  crc <- age -> microbiome
  crc <- gender -> microbiome
  crc <- diet -> microbiome
  crc <- healthstatus -> microbiome
  crc -> hospitalization <- microbiome
}")
ggdag(v_dag, text_col="red")+theme_dag()+theme_dag_gray()
ggdag_adjustment_set(v_dag, exposure = "crc", outcome = "microbiome", text_col="purple")+theme_dag_gray()

#Step 2
vogtmann <- curatedMetagenomicData("VogtmannE_2016.metaphlan_bugs_list.stool", dryrun = FALSE)
vm.eset <-vogtmann[[1]]
table(vm.eset$study_condition, useNA="ifany")
vm.sub <- vm.eset[,!is.na(vm.eset$study_condition)]
table(vm.sub$study_condition, useNA="ifany")
colnames(phenoData(vm.sub))

# Step 3
vm.pseq <- ExpressionSet2phyloseq(vm.sub)
dx_only <- phyloseq_to_deseq2(vm.pseq, ~ study_condition)
dx_only <- DESeq(dx_only, test="Wald", fitType="parametric")
res <- results(dx_only, cooksCutoff = FALSE)
alpha <- 0.10
sigtab <- res[which(res$padj < alpha), ]
sigtab <- cbind(as(sigtab, "data.frame"), as(tax_table(vm.pseq)[rownames(sigtab), ], "matrix"))
sigtab

# Step 4
dx_age_gender <- phyloseq_to_deseq2(vm.pseq, ~ study_condition + age_category + gender)
dx_age_gender <- DESeq(dx_age_gender, test="Wald", fitType="parametric")
res = results(dx_age_gender, cooksCutoff = FALSE, contrast=c("study_condition", "CRC", "control"))
alpha = 0.10
sigtab = res[which(res$padj < alpha), ]
sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(vm.pseq)[rownames(sigtab), ], "matrix"))
sigtab
