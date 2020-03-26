# Zymomona

Nitrogen Fixation Experiment (March 2020) - lead with Julia Martien 

OBJECTIVE  
This project investigates the proteomic remodeling associated with N2 fixation and NH4+ assimilation in Zymomona mobilis. Previous work has elucidated the metabolic and proteomic response of Zymomona mobilis during steady state N2 fixation only. These experiments examine the proteomic profile during the dynamic shift of N2 and NH4+ provide an additional dimension for studying metabolic remodeling during the increased energy and redox demands imposed by expression of biofuel-producing pathways.  

RESULTS  
Code takes in MaxQuant LFQ output table that has had contaminants and the fasta header overspill previously removed. A series of functions and data wrangling set up the dataframe for imputation using DRB Imputation code in Visual Studios (C# code). Imputation was applied to separately to upshift and down shift experiments. 

Simple statistics of average, standard deviation, and Log2 fold change are calculated and printed into individual .csv sheets that were merged and shared with collaborator. 


