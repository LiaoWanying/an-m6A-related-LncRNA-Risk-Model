# an-m6A-related-LncRNA-Risk-Model
1.getClinical.pl：extract clinical information from downloaded TCGA data
2.symbol.pl: translate ENSG id to gene name
3.biotype.pl: differentiate mRNA and lncRNA
4.m6aExp.R: extract the expression level of m6a gene
5.cor.R: identify m6a-related lncRNA
6.mergeTime.R: merge the expression of lncRNA and survival data
7.uniCox.R: identify prognostic m6a lncRNA
8.geneCor.R: correlation analysis of m6a gene and 29 prognostic lncRNA
9.heatmap.R: Heatmap of the expression levels of the 29 m6A-related prognostic lncRNAs in COAD samples and adjacent normal samples
10.riskPlot.R: Risk score distribution of m6A-LPS
11.survial.R: Kaplan–Meier survival analysis
12.ROC.R: receiver operating characteristic curve analysis of the 1-, 2-, and 3-year survival rates
13.cliGroupSur.R: Stratification analysis(investigate the differences in OS of all COAD patients stratified by the clinicopathological variables)
14.TCGAriskDiff_20220505.R: Identify differentially expressed genes(DEGs) between the high- and low-risk groups With sex and age applied as confounding factors
15.GO.R: Gene Ontology (GO) pathway enrichment analysis
16.KEGG.R: Kyoto Encyclopedia of Genes and Genomes (KEGG) pathway enrichment analysis
17.volcanoV2.Rmd: visualize the DEGs
18.indep.R: evaluate the independence of m6A-LPS in predicting overall survival (OS)
19.Nomogram.R: establish a nomogram combined two independent risk factors, m6A-LPS and pathologic stage
20.ROC曲线.R: ROC analysis to validate the nomogram in terms of its predictive accuracy
21.Cplot.R: Calibration plots of nomogram
23.TIDE.R: visualize the TIDE prediction scores
24.IPS.R: Immunophenoscore of the high- and low-risk groups
25.pRRophetic.R: IC50 of potential chemotherapeutic drugs
