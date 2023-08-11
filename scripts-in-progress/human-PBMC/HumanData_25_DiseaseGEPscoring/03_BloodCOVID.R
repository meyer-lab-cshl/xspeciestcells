###
# Purpose: Plot GEP score by disease state
# Date: June 2023
# Author: Salomé Carcy
###


# **************
# 1. IMPORT ####
# **************

## 1.1. Libraries ####
library(ggplot2)
library(tidyverse)
library(dplyr)
# library(Seurat)
library(RColorBrewer)
library(ggrepel)
library(cowplot)
source("~/Projects/HumanThymusProject/scripts-final/colors_universal.R") # get color palettes

## 1.2. Data ####
df <- read.csv("~/Projects/HumanThymusProject/data/human-PBMC/HumanData_25_DiseaseGEPscoring/covid_metadata_withGEPs.csv", row.names=1)


# Remove data we're not interested in
df <- df %>%
  filter(Annotation_minor_subset != "nan")




# ********************
# 2. CLEANUP DATA ####
# ********************

# useful metadata columns:
# "Annotation_major_subset" - CD4, CD8, MAIT, iNKT, GD: https://ars.els-cdn.com/content/image/1-s2.0-S0092867422000708-mmc2.pdf
# "Annotation_minor_subset" - cluster annotation 
# "Annotation_cluster_name" - very detailed cluster annotation 
# "COMBAT_ID" - participant ID
# "scRNASeq_sample_ID", "COMBAT_participant_timepoint_ID" - same thing (sometimes there are multiple time point per participantID)
# "Source" - condition (with different levels of care)
#          => COVID_CRIT     - Oxford hospitalized patients (critical)
#          => COVID_HCW_MILD - Mild COVID from healthcare workers (in same hospital)
#          => COVID_LDN      - London ICU patients with COVID
#          => COVID_MILD     - Oxford hospitalized patients (mild)
#          => COVID_SEV      - Oxford hospitalized patients (severe)
#          => Flu            - London ICU patients with Influenza
#          => Sepsis         - Oxford patients hospitalized for sepsis (non-COVID)
#          => HV             - Healthy volunteers
# "DiseaseClassification" is simplified condition
# "Age", "Sex", "BMI", "Hospitalstay" (from -12 to 15)
# "SARSCoV2PCR" - 0 or 1
# "Tissue" - all blood

# look up some stuff
# table(df[, c("Source", "DiseaseClassification")], useNA="ifany") %>%
#   as_tibble() %>%
#   filter(n>0) %>%
#   group_by(COMBAT_ID) %>% filter(n_distinct(COMBAT_participant_timepoint_ID)>1) %>% print(n=29)
# 
# df %>%
#   select(Source, COMBAT_participant_timepoint_ID) %>%
#   distinct() %>%
#   group_by(Source) %>%
#   count()



# Get average proportions of each cell cluster per person & time point
df_plot <- df %>%
  # total nb of cells per sample
  group_by(COMBAT_participant_timepoint_ID) %>%
  mutate(totalcells_per_participantimept=n()) %>%
  # freq of each cluster per sample
  group_by(COMBAT_participant_timepoint_ID, Annotation_minor_subset) %>%
  mutate(cluster_freq = n()/totalcells_per_participantimept) %>%
  ungroup() %>%
  # keep only rows/columns of interest
  select(Annotation_major_subset, Annotation_minor_subset, COMBAT_participant_timepoint_ID, Source, cluster_freq) %>%
  distinct() %>%
  # mean+sd freq of each cluster in each disease
  group_by(Source, Annotation_minor_subset) %>%
  mutate(cluster_freq_mean=mean(cluster_freq),
         cluster_freq_sd=sd(cluster_freq)) %>%
  ungroup()
  # sanity check
  # select(Source, Annotation_minor_subset, cluster_freq_mean) %>%
  # distinct() %>%
  # summarise(sumfreq = sum(cluster_freq_mean), .by=Source) # pb withinfluenza




# *********************************************
# 3. PLOT ALL GEPs PER CLUSTER AND DISEASE ####
# *********************************************

# FLU
plot_grid(
  # Proportions
  df_plot %>%
    filter(!Annotation_minor_subset %in% c("DN", "DP", "CD8.TREG")) %>%
    filter(Source%in% c("HV", "Flu")) %>%
    ggplot(aes(fill=factor(Source, levels=c("HV", "Flu"))))+
      # geom_bar(aes(x=Annotation_minor_subset, y=cluster_freq_mean, fill=Source), stat="identity", position=position_dodge()) +
      geom_boxplot(aes(x=Annotation_minor_subset, y=cluster_freq), outlier.shape = NA)+
      geom_point(aes(x=Annotation_minor_subset, y=cluster_freq), position=position_jitterdodge(jitter.width = 0, dodge.width = 0.7), size=0.5)+
      # facet_wrap(~Source)+
      scale_fill_manual(values=c("#b2abd2", "#fdb863"), name="Condition")+
      theme_cowplot()+
      theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.5))+
      labs(x="", title="Hospitalized flu vs Healthy"),
  # GEP1
  df %>%
    filter(!Annotation_minor_subset %in% c("DN", "DP", "CD8.TREG")) %>%
    filter(Source%in% c("HV", "Flu")) %>%
    # plot
    ggplot(aes(x=Annotation_minor_subset, y=GEP1, fill=factor(Source, levels=c("HV", "Flu"))))+
      geom_hline(yintercept = 0, linetype="dashed", color="darkgrey")+
      geom_violin(position=position_dodge(width=0.7))+
      geom_point(position=position_jitterdodge(jitter.width = 0.1, dodge.width=0.7), size=0.1)+
      scale_fill_manual(values=c("#b2abd2", "#fdb863"), name="Condition")+
      theme_cowplot()+
      theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.5))+
      labs(x=""),
  # GEP4
  df %>%
    filter(!Annotation_minor_subset %in% c("DN", "DP", "CD8.TREG")) %>%
    filter(Source%in% c("HV", "Flu")) %>%
    # plot
    ggplot(aes(x=Annotation_minor_subset, y=GEP4, fill=factor(Source, levels=c("HV", "Flu"))))+
      geom_hline(yintercept = 0, linetype="dashed", color="darkgrey")+
      geom_violin(position=position_dodge(width=0.7))+
      geom_point(position=position_jitterdodge(jitter.width = 0.1, dodge.width=0.7), size=0.1)+
      scale_fill_manual(values=c("#b2abd2", "#fdb863"), name="Condition")+
      theme_cowplot()+
      theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.5))+
      labs(x=""),
  # GEP6
  df %>%
    filter(!Annotation_minor_subset %in% c("DN", "DP", "CD8.TREG")) %>%
    filter(Source%in% c("HV", "Flu")) %>%
    # plot
    ggplot(aes(x=Annotation_minor_subset, y=GEP6, fill=factor(Source, levels=c("HV", "Flu"))))+
      geom_hline(yintercept = 0, linetype="dashed", color="darkgrey")+
      geom_violin(position=position_dodge(width=0.7))+
      geom_point(position=position_jitterdodge(jitter.width = 0.1, dodge.width=0.7), size=0.1)+
      scale_fill_manual(values=c("#b2abd2",  "#fdb863"), name="Condition")+
      theme_cowplot()+
      theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.5))+
      labs(x=""),
  # PLOTGRID PARAMS
  ncol=1)
# ggsave("./scripts-in-progress/human-PBMC/HumanData_25_DiseaseGEPscoring/plots/covid/heathyVSflu.jpeg", width=10, height=30)





# SEPSIS
plot_grid(
  # Proportions
  df_plot %>%
    filter(!Annotation_minor_subset %in% c("DN", "DP", "CD8.TREG")) %>%
    filter(Source%in% c("HV", "Sepsis")) %>%
    ggplot()+
    # geom_bar(aes(x=Annotation_minor_subset, y=cluster_freq_mean, fill=Source), stat="identity", position=position_dodge()) +
    geom_boxplot(aes(x=Annotation_minor_subset, y=cluster_freq, fill=Source), outlier.shape = NA)+
    geom_point(aes(x=Annotation_minor_subset, y=cluster_freq, fill=Source), position=position_jitterdodge(jitter.width = 0, dodge.width = 0.7), size=0.5)+
    # facet_wrap(~Source)+
    scale_fill_manual(values=c("#b2abd2", "#fc9272"), name="Condition")+
    theme_cowplot()+
    theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.5))+
    labs(x="", title="Hospitalized sepsis vs Healthy"),
  # GEP1
  df %>%
    filter(!Annotation_minor_subset %in% c("DN", "DP", "CD8.TREG")) %>%
    filter(Source%in% c("HV", "Sepsis")) %>%
    # plot
    ggplot(aes(x=Annotation_minor_subset, y=GEP1, fill=Source))+
      geom_hline(yintercept = 0, linetype="dashed", color="darkgrey")+
      geom_violin(position=position_dodge(width=0.7))+
      geom_point(position=position_jitterdodge(jitter.width = 0.1, dodge.width=0.7), size=0.1)+
      scale_fill_manual(values=c("#b2abd2", "#fc9272"), name="Condition")+
      theme_cowplot()+
      theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.5))+
      labs(x=""),
  # GEP4
  df %>%
    filter(!Annotation_minor_subset %in% c("DN", "DP", "CD8.TREG")) %>%
    filter(Source%in% c("HV", "Sepsis")) %>%
    # plot
    ggplot(aes(x=Annotation_minor_subset, y=GEP4, fill=Source))+
      geom_hline(yintercept = 0, linetype="dashed", color="darkgrey")+
      geom_violin(position=position_dodge(width=0.7))+
      geom_point(position=position_jitterdodge(jitter.width = 0.1, dodge.width=0.7), size=0.1)+
      scale_fill_manual(values=c("#b2abd2", "#fc9272"), name="Condition")+
      theme_cowplot()+
      theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.5))+
      labs(x=""),
  # GEP6
  df %>%
    filter(!Annotation_minor_subset %in% c("DN", "DP", "CD8.TREG")) %>%
    filter(Source%in% c("HV", "Sepsis")) %>%
    # plot
    ggplot(aes(x=Annotation_minor_subset, y=GEP6, fill=Source))+
      geom_hline(yintercept = 0, linetype="dashed", color="darkgrey")+
      geom_violin(position=position_dodge(width=0.7))+
      geom_point(position=position_jitterdodge(jitter.width = 0.1, dodge.width=0.7), size=0.1)+
      scale_fill_manual(values=c("#b2abd2", "#fc9272"), name="Condition")+
      theme_cowplot()+
      theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.5))+
      labs(x=""),
  # GEP12
  df %>%
    filter(!Annotation_minor_subset %in% c("DN", "DP", "CD8.TREG")) %>%
    filter(Source%in% c("HV", "Sepsis")) %>%
    # plot
    ggplot(aes(x=Annotation_minor_subset, y=GEP12, fill=Source))+
    geom_hline(yintercept = 0, linetype="dashed", color="darkgrey")+
    geom_violin(position=position_dodge(width=0.7))+
    geom_point(position=position_jitterdodge(jitter.width = 0.1, dodge.width=0.7), size=0.1)+
    scale_fill_manual(values=c("#b2abd2", "#fc9272"), name="Condition")+
    theme_cowplot()+
    theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.5))+
    labs(x=""),
  # PLOTGRID PARAMS
  ncol=1)
# ggsave("./scripts-in-progress/human-PBMC/HumanData_25_DiseaseGEPscoring/plots/covid/heathyVSsepsis.jpeg", width=10, height=40)





# COVID
covid_list <- c("HV", "COVID_HCW_MILD", "COVID_MILD", "COVID_SEV", "COVID_CRIT", "COVID_LDN")
plot_grid(
  # Proportions
  df_plot %>%
    filter(!Annotation_minor_subset %in% c("DN", "DP", "CD8.TREG")) %>%
    filter(Source %in% covid_list) %>%
    ggplot(aes(fill=factor(Source, levels=covid_list)))+
      # geom_bar(aes(x=Annotation_minor_subset, y=cluster_freq_mean, fill=Source), stat="identity", position=position_dodge()) +
      geom_boxplot(aes(x=Annotation_minor_subset, y=cluster_freq), outlier.shape = NA)+
      geom_point(aes(x=Annotation_minor_subset, y=cluster_freq), position=position_jitterdodge(jitter.width = 0, dodge.width = 0.75), size=0.5)+
      # facet_wrap(~Source)+
      scale_fill_manual(values=c("#b2abd2", "#c7e9c0", "#41ab5d", "#238b45", "#006d2c", "#00441b"), name="Condition")+
      theme_cowplot()+
      theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.5))+
      labs(x="", title="COVID vs Healthy"),
  # GEP1
  df %>%
    filter(!Annotation_minor_subset %in% c("DN", "DP", "CD8.TREG")) %>%
    filter(Source%in% covid_list) %>%
    # plot
    ggplot(aes(x=Annotation_minor_subset, y=GEP1, fill=factor(Source, levels=covid_list)))+
    geom_hline(yintercept = 0, linetype="dashed", color="darkgrey")+
    geom_violin(position=position_dodge(width=0.7))+
    geom_point(position=position_jitterdodge(jitter.width = 0.1, dodge.width=0.7), size=0.1)+
    scale_fill_manual(values=c("#b2abd2", "#c7e9c0", "#41ab5d", "#238b45", "#006d2c", "#00441b"), name="Condition")+
    theme_cowplot()+
    theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.5))+
    labs(x=""),
  # GEP4
  df %>%
    filter(!Annotation_minor_subset %in% c("DN", "DP", "CD8.TREG")) %>%
    filter(Source%in% covid_list) %>%
    # plot
    ggplot(aes(x=Annotation_minor_subset, y=GEP4, fill=factor(Source, levels=covid_list)))+
    geom_hline(yintercept = 0, linetype="dashed", color="darkgrey")+
    geom_violin(position=position_dodge(width=0.7))+
    geom_point(position=position_jitterdodge(jitter.width = 0.1, dodge.width=0.7), size=0.1)+
    scale_fill_manual(values=c("#b2abd2", "#c7e9c0", "#41ab5d", "#238b45", "#006d2c", "#00441b"), name="Condition")+
    theme_cowplot()+
    theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.5))+
    labs(x=""),
  # GEP6
  df %>%
    filter(!Annotation_minor_subset %in% c("DN", "DP", "CD8.TREG")) %>%
    filter(Source%in%covid_list) %>%
    # plot
    ggplot(aes(x=Annotation_minor_subset, y=GEP6, fill=factor(Source, levels=covid_list)))+
    geom_hline(yintercept = 0, linetype="dashed", color="darkgrey")+
    geom_violin(position=position_dodge(width=0.7))+
    geom_point(position=position_jitterdodge(jitter.width = 0.1, dodge.width=0.7), size=0.1)+
    scale_fill_manual(values=c("#b2abd2", "#c7e9c0", "#41ab5d", "#238b45", "#006d2c", "#00441b"), name="Condition")+
    theme_cowplot()+
    theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.5))+
    labs(x=""),
  # PLOTGRID PARAMS
  ncol=1)
# ggsave("./scripts-in-progress/human-PBMC/HumanData_25_DiseaseGEPscoring/plots/covid/heathyVScovid.jpeg", width=15, height=30)


# end plot by gep and disease ####




# **************************
# 4. COMBINE CELL TYPES ####
# **************************

# Function for easier plotting
vlnplot_gep <- function(xvar, yvar, disease=c("HV", "Flu"), colors=c("#b2abd2", "#fdb863"), plotitle="Hospitalized flu vs Healthy"){
  # if plotting proportions
  if(yvar=="prop"){
    p <- df_plot %>%
      filter(!Annotation_minor_subset %in% c("DN", "DP", "CD8.TREG")) %>%
      filter(Source%in% disease) %>%
      ggplot(aes(x=.data[[xvar]], y=cluster_freq, fill=factor(Source, levels=disease)))+
        geom_boxplot(outlier.shape = NA)+
        geom_point(position=position_jitterdodge(jitter.width = 0, dodge.width = 0.7), size=0.5)+
        scale_fill_manual(values=colors, name="Condition")+
        theme_cowplot()+
        theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.5))+
        labs(x="", title=plotitle)
  }
  # if plotting gep score
  else if(grepl("GEP", yvar) == TRUE){
  p <- df %>%
    filter(!Annotation_minor_subset %in% c("DN", "DP", "CD8.TREG")) %>%
    filter(Source%in% disease) %>%
    # plot
    ggplot(aes(x=.data[[xvar]], y=.data[[yvar]], fill=factor(Source, levels=disease)))+
      geom_hline(yintercept = 0, linetype="dashed", color="darkgrey")+
      geom_violin(position=position_dodge(width=0.7))+
      geom_point(position=position_jitterdodge(jitter.width = 0.1, dodge.width=0.7), size=0.1)+
      scale_fill_manual(values=colors, name="Condition")+
      theme_cowplot()+
      theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.5))+
      labs(x="")
  }
  return(p)
}



# FLU
plot_grid(
  # Proportions
  vlnplot_gep(xvar="Annotation_major_subset", yvar="prop", plotitle = "Hospitalized flu vs Healthy"),
  # GEPs
  vlnplot_gep(xvar="Annotation_major_subset", yvar="GEP1"),
  vlnplot_gep(xvar="Annotation_major_subset", yvar="GEP4"),
  vlnplot_gep(xvar="Annotation_major_subset", yvar="GEP5"),
  vlnplot_gep(xvar="Annotation_major_subset", yvar="GEP6"),
  vlnplot_gep(xvar="Annotation_major_subset", yvar="GEP12"),
  # PLOTGRID PARAMS
  ncol=3)
# ggsave("./scripts-in-progress/human-PBMC/HumanData_25_DiseaseGEPscoring/plots/covid/heathyVSflu_annotationmajortypes.jpeg", width=15, height=10)

# FLU
plot_grid(
  # Proportions
  vlnplot_gep(xvar="Annotation_major_subset", yvar="prop", plotitle = "Hospitalized flu vs Healthy"),
  # GEPs
  vlnplot_gep(xvar="Annotation_major_subset", yvar="GEP1"),
  vlnplot_gep(xvar="Annotation_major_subset", yvar="GEP4"),
  vlnplot_gep(xvar="Annotation_major_subset", yvar="GEP5"),
  vlnplot_gep(xvar="Annotation_major_subset", yvar="GEP6"),
  vlnplot_gep(xvar="Annotation_major_subset", yvar="GEP12"),
  # PLOTGRID PARAMS
  ncol=3)
# ggsave("./scripts-in-progress/human-PBMC/HumanData_25_DiseaseGEPscoring/plots/covid/heathyVSflu_annotationmajortypes.jpeg", width=15, height=10)


# SEPSIS
plot_grid(
  # Proportions
  vlnplot_gep(xvar="Annotation_major_subset", yvar="prop", disease=c("HV", "Sepsis"), colors=c("#b2abd2", "#fc9272"), plotitle = "Hospitalized sepsis vs Healthy"),
  # GEPs
  vlnplot_gep(xvar="Annotation_major_subset", yvar="GEP1", disease=c("HV", "Sepsis"), colors=c("#b2abd2", "#fc9272")),
  vlnplot_gep(xvar="Annotation_major_subset", yvar="GEP4", disease=c("HV", "Sepsis"), colors=c("#b2abd2", "#fc9272")),
  vlnplot_gep(xvar="Annotation_major_subset", yvar="GEP5", disease=c("HV", "Sepsis"), colors=c("#b2abd2", "#fc9272")),
  vlnplot_gep(xvar="Annotation_major_subset", yvar="GEP6", disease=c("HV", "Sepsis"), colors=c("#b2abd2", "#fc9272")),
  vlnplot_gep(xvar="Annotation_major_subset", yvar="GEP12", disease=c("HV", "Sepsis"), colors=c("#b2abd2", "#fc9272")),
  # PLOTGRID PARAMS
  ncol=3)
# ggsave("./scripts-in-progress/human-PBMC/HumanData_25_DiseaseGEPscoring/plots/covid/heathyVSsepsis_annotationmajortypes.jpeg", width=15, height=10)


# COVID
plot_grid(
  # Proportions
  vlnplot_gep(xvar="Annotation_major_subset", yvar="prop", disease=covid_list, colors=c("#b2abd2", "#c7e9c0", "#41ab5d", "#238b45", "#006d2c", "#00441b"), plotitle = "Hospitalized sepsis vs Healthy"),
  # GEPs
  vlnplot_gep(xvar="Annotation_major_subset", yvar="GEP1",  disease=covid_list, colors=c("#b2abd2", "#c7e9c0", "#41ab5d", "#238b45", "#006d2c", "#00441b")),
  vlnplot_gep(xvar="Annotation_major_subset", yvar="GEP4",  disease=covid_list, colors=c("#b2abd2", "#c7e9c0", "#41ab5d", "#238b45", "#006d2c", "#00441b")),
  vlnplot_gep(xvar="Annotation_major_subset", yvar="GEP5",  disease=covid_list, colors=c("#b2abd2", "#c7e9c0", "#41ab5d", "#238b45", "#006d2c", "#00441b")),
  vlnplot_gep(xvar="Annotation_major_subset", yvar="GEP6",  disease=covid_list, colors=c("#b2abd2", "#c7e9c0", "#41ab5d", "#238b45", "#006d2c", "#00441b")),
  vlnplot_gep(xvar="Annotation_major_subset", yvar="GEP12", disease=covid_list, colors=c("#b2abd2", "#c7e9c0", "#41ab5d", "#238b45", "#006d2c", "#00441b")),
  # PLOTGRID PARAMS
  ncol=3)
# ggsave("./scripts-in-progress/human-PBMC/HumanData_25_DiseaseGEPscoring/plots/covid/heathyVScovid_annotationmajortypes.jpeg", width=20, height=10)




# *******************************************
# 5. CELL CLUSTERS OF INTEREST FIGURE 4B ####
# *******************************************


# Function for easier plotting
vlnplot_gep_fig4b <- function(xvar="Annotation_cluster_name", yvar, disease=c("HV", "Flu"), colors=c("#b2abd2", "#fdb863"), plotitle="Hospitalized flu vs Healthy"){
  clusters_fig4b <- c("CD4.TREG", "CD4.TEM.GZMK", "CD4.TEFF.GZMK", "CD4.TREG.CCR4hi", "CD4.TEFF.TCF7", "CD4.Th1.1",
                      "CD4.TEFF.prolif.1", "CD4.TEFF.prolif.2", "CD4.TEFF.prolif.3",
                      "CD8.TCM.CCL5.1", "CD8.NAIVE.1", "CD8.TCM.mitohi", "CD8.TEMRA.2", "CD8.TEFF.prolif.1", "CD8.TEFF.prolif.2",
                      "MAIT.CD8.CCL5.cytox_lo.1")
  # if plotting proportions
  if(yvar=="prop"){
    p <- df_plot %>%
      filter(Annotation_cluster_name %in% clusters_fig4b) %>%
      filter(Source%in% disease) %>%
      ggplot(aes(x=.data[[xvar]], y=cluster_freq, fill=factor(Source, levels=disease)))+
      geom_boxplot(outlier.shape = NA)+
      geom_point(position=position_jitterdodge(jitter.width = 0, dodge.width = 0.7), size=0.5)+
      scale_fill_manual(values=colors, name="Condition")+
      theme_cowplot()+
      theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.5))+
      labs(x="", title=plotitle)
  }
  # if plotting gep score
  else if(grepl("GEP", yvar) == TRUE){
    p <- df %>%
      filter(Annotation_cluster_name %in% clusters_fig4b) %>%
      filter(Source%in% disease) %>%
      # plot
      ggplot(aes(x=.data[[xvar]], y=.data[[yvar]], fill=factor(Source, levels=disease)))+
      geom_hline(yintercept = 0, linetype="dashed", color="darkgrey")+
      geom_violin(position=position_dodge(width=0.7))+
      geom_point(position=position_jitterdodge(jitter.width = 0.1, dodge.width=0.7), size=0.1)+
      scale_fill_manual(values=colors, name="Condition")+
      theme_cowplot()+
      theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.5))+
      labs(x="")
  }
  return(p)
}


# Get average proportions of each cell cluster per person & time point
df_plot_clustername <- df %>%
  # total nb of cells per sample
  group_by(COMBAT_participant_timepoint_ID) %>%
  mutate(totalcells_per_participantimept=n()) %>%
  # freq of each cluster per sample
  group_by(COMBAT_participant_timepoint_ID, Annotation_cluster_name) %>%
  mutate(cluster_freq = n()/totalcells_per_participantimept) %>%
  ungroup() %>%
  # keep only rows/columns of interest
  select(Annotation_major_subset, Annotation_cluster_name, COMBAT_participant_timepoint_ID, Source, cluster_freq) %>%
  distinct() %>%
  # mean+sd freq of each cluster in each disease
  group_by(Source, Annotation_cluster_name) %>%
  mutate(cluster_freq_mean=mean(cluster_freq),
         cluster_freq_sd=sd(cluster_freq)) %>%
  ungroup() %>%
  # sanity check
  select(Source, Annotation_cluster_name, cluster_freq_mean) %>%
  distinct() %>%
  summarise(sumfreq = sum(cluster_freq_mean), .by=Source) # pb withinfluenza



# FLU
plot_grid(
  # Proportions
  # vlnplot_gep_fig4b(yvar="prop", plotitle = "Hospitalized flu vs Healthy"),
  # GEPs
  vlnplot_gep_fig4b(yvar="GEP1"),
  vlnplot_gep_fig4b(yvar="GEP4"),
  vlnplot_gep_fig4b(yvar="GEP5"),
  vlnplot_gep_fig4b(yvar="GEP6"),
  vlnplot_gep_fig4b(yvar="GEP12"),
  # PLOTGRID PARAMS
  ncol=3)
# ggsave("./scripts-in-progress/human-PBMC/HumanData_25_DiseaseGEPscoring/plots/covid/heathyVSflu_clustersfig4b.jpeg", width=20, height=10)

# SEPSIS
plot_grid(
  # GEPs
  vlnplot_gep_fig4b(yvar="GEP1",  disease=c("HV", "Sepsis"), colors=c("#b2abd2", "#fc9272")),
  vlnplot_gep_fig4b(yvar="GEP4",  disease=c("HV", "Sepsis"), colors=c("#b2abd2", "#fc9272")),
  vlnplot_gep_fig4b(yvar="GEP5",  disease=c("HV", "Sepsis"), colors=c("#b2abd2", "#fc9272")),
  vlnplot_gep_fig4b(yvar="GEP6",  disease=c("HV", "Sepsis"), colors=c("#b2abd2", "#fc9272")),
  vlnplot_gep_fig4b(yvar="GEP12", disease=c("HV", "Sepsis"), colors=c("#b2abd2", "#fc9272")),
  # PLOTGRID PARAMS
  ncol=3)
# ggsave("./scripts-in-progress/human-PBMC/HumanData_25_DiseaseGEPscoring/plots/covid/heathyVSsepsis_clustersfig4b.jpeg", width=20, height=10)


# COVID
plot_grid(
  # GEPs
  vlnplot_gep_fig4b(yvar="GEP1",  disease=covid_list, colors=c("#b2abd2", "#c7e9c0", "#41ab5d", "#238b45", "#006d2c", "#00441b")),
  vlnplot_gep_fig4b(yvar="GEP4",  disease=covid_list, colors=c("#b2abd2", "#c7e9c0", "#41ab5d", "#238b45", "#006d2c", "#00441b")),
  vlnplot_gep_fig4b(yvar="GEP5",  disease=covid_list, colors=c("#b2abd2", "#c7e9c0", "#41ab5d", "#238b45", "#006d2c", "#00441b")),
  vlnplot_gep_fig4b(yvar="GEP6",  disease=covid_list, colors=c("#b2abd2", "#c7e9c0", "#41ab5d", "#238b45", "#006d2c", "#00441b")),
  vlnplot_gep_fig4b(yvar="GEP12", disease=covid_list, colors=c("#b2abd2", "#c7e9c0", "#41ab5d", "#238b45", "#006d2c", "#00441b")),
  # PLOTGRID PARAMS
  ncol=2)
ggsave("./scripts-in-progress/human-PBMC/HumanData_25_DiseaseGEPscoring/plots/covid/heathyVScovid_clustersfig4b.jpeg", width=20, height=20)


