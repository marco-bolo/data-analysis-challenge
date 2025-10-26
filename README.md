
<img src="https://www.eu4oceanobs.eu/wp-content/uploads/2023/03/MARCO-BOLO_logo_col.png" alt="drawing" width="500" align="center"/>
</p>

# The Marco-Bolo eDNA Data Analysis Challenge

The [EU Horizon MARCO-BOLO project](https://marcobolo-project.eu/) launched the eDNA Data Analysis Challenge, 
inviting the global eDNA community to download the provided datasets, 
run their pipeline on it, and submit their results for comparison across all participating pipelines and entries.

Information about the datasets, and instructions to participate can be found in the [challenge participation page.](https://github.com/marco-bolo/wp2-wp5-workshop/tree/main/)

Here, we present the entries, code and results of the eDNA data analyis challenge challenge. 

## Goal

1. :floppy_disk:  We provide the eDNA datasets along with standardised reference libraries
    - The 18S plankton time series [described here](https://github.com/marco-bolo/wp2-wp5-workshop/tree/main/?tab=readme-ov-file#plankton-18s-time-series)
    -  The 12S/16S/COI aquarium samples [described here](https://github.com/marco-bolo/wp2-wp5-workshop/tree/main/?tab=readme-ov-file#fish-12s16scoi-aquarium-dataset)

2. :wrench:  Participants come with their preferred pipeline
3. :computer:  They run the pipeline on one or more of our plankton or fish eDNA datasets
4. :globe_with_meridians:  They send us their resulting OTU/ASV/ZOTU and Taxonomy tables

5. :bar_chart:  We compile and compare the resulting tables across participants' pipelines
6. :memo:  We report and share these results back with the participants


## Content

This repository is structured as follow:

- :file_folder: &nbsp;[**data/**](https://github.com/marco-bolo/data-analysis-challenge/tree/master/data):
 contains the compilation of the different pipeline outputs, harmonised into two large data files 
   - challengadata_18s.csv: 
   - challengedata_12s16sCOI.csv: 

   This folder also contains the metadata files for each of the large data compilation files. 
   - metadata_18s.csv:
   - metadata_12s16sCOI.csv

   For now, the pipeline outputs are anonymised but the main pipeline steps are categorised in the metadata file.
- :file_folder: &nbsp;[**scripts/**](https://github.com/marco-bolo/data-analysis-challenge/tree/master/scripts):
 contains R scripts to reproduce comparative analyses and figures.
- :file_folder: &nbsp;[**figures/**](https://github.com/marco-bolo/data-analysis-challenge/tree/master/figures):
 produced figures (.png) will be saved to this folder

 ## Comparative analyses

For the fish 12S, 16S and COI entries:
 - Detectability: Target species identification, sensitivity
 - Alpha diversity indices: (target) species richness, (target) Shannon diversity
 - Community composition: PCoA and Simper analyses.

 For the plankton 18S entries:
 - Fluctuation patterns: GAM model of Shannon diversity ~ month to reveal seasonal patterns
 - Community composition: NMDS on the Bray-Curtis distance of Hellinger-tranformed relative abundances

 
 ## Interactive exploration

*Here, we host and share the shiny web app to interactively explore the results.*

 ## Report

*Here, we share some insights of the results, and/or link to the online deliverables.*

 ## Acknowledgements

The MARCO-BOLO project is funded by the European Union under the Horizon Europe Programme, Grant Agreement No. 101082021 (MARCO-BOLO). Views and opinions expressed are however those of the author(s) only and do not necessarily reflect those of the European Union or the European Research Executive Agency (REA). Neither the European Union nor the granting authority can be held responsible for them.
UK participants in MARCO-BOLO are supported by the UKRI’s Horizon Europe Guarantee under the Grant No. 10068180 (MS); No. 10063994 (MBA); No. 10048178 (NOC).