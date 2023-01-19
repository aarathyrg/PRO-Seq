# PRO-Seq
Comparing early and delayed transcriptional response to interferons
The early and delayed transcriptional response to interferons are analyzed using nascent transript data.
Steps
1) Preprocessing and alignment was performed using proseq2.0 pipeline (https://github.com/Danko-Lab/proseq2.0)
2) Count matrices were generated using bigWig files according to https://github.com/Danko-Lab/tfTarget/blob/master/tfTarget/R/diffTXN.R
3) Using dREG package (https://dreg.dnasequence.org/) transcriptional regulatory elements were selected
4) Bedtools intersect was used for further estimating overlap between regions used for estimating putative enhancers as well as binding sites of transcription factors
5) Downstream analyses including differential analysis of Knockouts vs WT , clustering, pol II pausing analysis were mostly performed using R.
