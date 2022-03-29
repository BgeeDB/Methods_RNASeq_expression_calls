#!/bin/sh
echo "All figures, tsv files as well as supplementary data/files will be generated by default"

cd /Copy_all_repository/

## Figure 1: generate the Gaussians destribution after getting the sum file from Bgee pipeline
Rscript /Copy_all_repository/scripts/Figure_1/Gaussian_Density_distribution_mouse.R

## Figure 2: generate final figure after getting the calls file for each library from the Bgee pipeline
Rscript /Copy_all_repository/scripts/Figure_2/TPR_FPR_Mouse_Liver_samples.R

## Figure 3: generate the Gaussians destribution after getting the sum file from Bgee pipeline + plot proportion of protein coding using different methods to call expressed genes
Rscript /Copy_all_repository/scripts/Figure_3/Gaussian_Density_distribution_esox.R
Rscript /Copy_all_repository/scripts/Figure_3/Proportion_coding_genes_all_approaches.R

## Figure 4: Generate the calls of expression of protein coding genes across all species in Bgee 15 by using all the methods (Deconvolution, without deconvolution and TPM threshold).
Rscript /Copy_all_repository/scripts/Figure_4/ProportionProteinCoding_52Species_AllMethods.R

## Figure 5: Generate plot by using pValue (p ≤ 0.05) method across all libraries + plot performing the calls using TPM threshold + Correlation pValue vy TPM threshold in GTEx samples 
Rscript /Copy_all_repository/scripts/Figure_5/GTEx_plots.R

## Figure 6: Proportion of coding genes call present in RNA-Seq and single cell
Rscript /Copy_all_repository/scripts/Figure_6/ProportionCoding_bulk-RNA_scRNA_FL.R

## Figure 7: Density distribuiton of one random cell or cell population for droplet-based protocols (10X)
Rscript /Copy_all_repository/scripts/Figure_7/Droplet_based_OneCell_CellPop.R

## Figure 8: compare proportion of coding genes in scRNA-Seq (full-length or droplet-based)
Rscript /Copy_all_repository/scripts/Figure_8/Compare_scRNA-seq_protocols.R

## Figure 9: enrichment analysis using different approaches in droplet-based protocols
Rscript /Copy_all_repository/scripts/Figure_9/barplot_panther_pathways.R

## Figure 10: ROC curve human lung data
Rscript /Copy_all_repository/scripts/Figure_10/ROC_curve_Human_lung.R

## Figure 11: ROC curve Drosophila testis data
Rscript /Copy_all_repository/scripts/Figure_11/ROC_curve_Drosophila_testis.R

## Figure 12: ROC curve Mouse liver data
Rscript /Copy_all_repository/scripts/Figure_12/ROC_curve_Mouse_liver.R

## Figure 13: compare proportion of coding genes present by using pValue or qValue method
Rscript /Copy_all_repository/scripts/Figure_13/comparing_RefInt_methods_callExpression.R

## Figure 14: was generated using https://app.diagrams.net/?src=about

## Figure 15: generate intermediary files after mapping/pseudo-mapping (analysis and sum over all libraries) + final figure
# analysis per library: Rscript /Copy_all_repository/scripts/Figure_15/rna_seq_analysis_UPDATED.R
# sum over libraries of the species: Rscript /Copy_all_repository/scripts/Figure_15/rna_seq_sumSpecies_UPDATED.R
Rscript /Copy_all_repository/scripts/Figure_15/plotVeenDiagram_RefInt_using_all_mapping.R

########### export suppl figures

## Suppl figure 1: upset plot to show the ref intergenic regions selected per anatomical entity
# first sum libraries across anatomical entity for Mus musculus data: Rscript /Copy_all_repository/scripts/Suppl_figure_1/sum_libraries_by_anatomicalEntity.R
Rscript /Copy_all_repository/scripts/Suppl_figure_1/UpsetPlot_Mouse_RefInt_anatomicalEntity.R

## Suppl figure 2: correlation plot for mouse data using different methods to call expressed genes
Rscript /Copy_all_repository/scripts/Suppl_figure_2/Correlation_between_methods_mouseData.R

## Suppl figure 3: TPR and FPR for mouse data using BH correction per library
Rscript /Copy_all_repository/scripts/Suppl_figure_3/TPR_FPR_with_BH_correction_mouse.R

## Suppl figure 4: TPR and FPR with and without BH correction for human and drosophila data
Rscript /Copy_all_repository/scripts/Suppl_figure_4/TPR_FPR_Human_Drosophila_with_and_without_BH_correction.R

## Suppl figure 5: Proportion of coding genes call present with different methods and cutoffs
Rscript /Copy_all_repository/scripts/Suppl_figure_5/ProportionCoding_per_method_cutoff_AllSpecies.R

## Suppl figure 6: plot the density destribution of regions for all species
Rscript /Copy_all_repository/scripts/Suppl_figure_6/plot_densities_after_sum_all_libraries_per_species.R

## Suppl figure 7: Pearson correlation between pValue vs TPM threshold for GTEx data without blood
Rscript /Copy_all_repository/scripts/Suppl_figure_7/plot_GTEx_without_blood.R

## Suppl figure 8: Proportion of coding genes called present in droplet-based with CPM ≥ 1
Rscript /Copy_all_repository/scripts/Suppl_figure_8/callsExpression_CPM_threshold.R

## Suppl figure 9: Density and variance of Human lung samples
Rscript /Copy_all_repository/scripts/Suppl_figure_9/density_variance_human_lung_samples.R

## Suppl figure 10: Density and variance of Drosophila testis samples
Rscript /Copy_all_repository/scripts/Suppl_figure_10/density_variance_drosophila_testis_samples.R

## Suppl figure 11: Proportion of non coding versus proportion of reference intergenic region per species
Rscript /Copy_all_repository/scripts/Suppl_figure_11/ProportionNonCoding_vs_ProportionRefIntergenic.R

## Suppl figure 12: CG content 
Rscript /Copy_all_repository/scripts/Suppl_figure_12/CG_content_all_species.R




























