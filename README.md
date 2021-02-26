CytoTalk 
================
Version 2.0 (February 22nd, 2021) 


## Overview

Signal transduction is the primary mechanism for cell-cell communication. scRNA-Seq technology holds great promise for studying cell-cell communication at much higher resolution. Signaling pathways are highly dynamic and cross-talk among them is prevalent. Due to these two features, simply examining expression levels of ligand and receptor genes cannot reliably capture the overall activities of signaling pathways and interactions among them. 

We have developed the CytoTalk algorithm for *de novo* construction of a signaling network (union of multiple signaling pathways emanating from the ligand-receptor pairs) between two cell types using single-cell transcriptomics data. The algorithm first constructs an integrated network consisting of intracellular and intercellular functional gene interactions. It then identifies the signaling network by solving a prize-collecting Steiner forest (PCSF) problem based on appropriately defined node prize (i.e. cell-specific gene activity) and edge cost (i.e. probability of functional interaction between two genes). The objective of the PCSF problem is to find an optimal subnetwork in the integrated network that includes **genes with high levels of cell-type-specific expression and close connection to highly active ligand-receptor pairs**. CytoTalk is currently implemented using a combination of MATLAB (version  R2018a), R (version  3.5.0) and Python (version  3.7.0). 


<div align=center><img src="https://github.com/huBioinfo/CytoTalk/blob/master/CytoTalk_schematic.png" width="60%" height="60%" /></div>
<br />

## Packages and Environment configurations      

 * The following four **R** packages should be installed (R version  3.5.0 is recommended). 

    - “entropy”: https://cran.r-project.org/web/packages/entropy/index.html

    - “infotheo”: https://cran.r-project.org/web/packages/infotheo/index.html

    - “doParallel”: https://cran.r-project.org/web/packages/doParallel/index.html

    - “parmigene”: https://cran.rstudio.com/web/packages/parmigene/index.html

    "doParallel" is used for parallel computation of mutual information for all gene pairs. Before running CytoTalk, set the number of logical cores available (default=14) in **Line 13** of the R scripts: **comp_MIcoexp_TypA_WinPara.R** and **comp_MIcoexp_TypB_WinPara.R**.

    "parmigene" is used for parallel computation of indirect edge-filtered gene networks. Before running CytoTalk, set system environment variable as following. n is the number of logical cores available.
  
    ```Bash
        export OMP_NUM_THREADS=n
    ```

 * The following three **Python** package should be installed (Python version  3.7.0 is recommended).

    - “pcsf_fast”: https://github.com/fraenkel-lab/pcst_fast

    - “numpy”: https://numpy.org/install/

    - “datetime”: https://pypi.org/project/DateTime/
 
    "pcsf_fast" is used for fast identification of a rooted prize-collecting Steiner tree in a network. Before running CytoTalk, set the system environment variable as following: 
  
    ```Bash
        export PYTHONPATH=$PYTHONPATH:/your installed pcsf_fast folder/
    ```
  
 * Set system environment variable to include the absolute path of the executable **MATLAB** program. An example in the macOS system is as following:

    ```Bash
        export PATH=/Applications/MATLAB_R2018a.app/bin/:$PATH
    ```

## Input files   

* A **comma-delimited “.csv”** file containing scRNA-Seq data for **each cell type** under study. Each file contains the **ln-transformed normalized scRNA-Seq data** for a cell type with rows as genes (GENE SYMBOL) and columns as cells. The files should be named as: **scRNAseq_Fibroblasts.csv**, **scRNAseq_Macrophages.csv**, **scRNAseq_EndothelialCells.csv**, **scRNAseq_CellTypeName.csv** …

* A **“TwoCellTypes.txt”** file indicating the two cell types between which the signaling network is predicted. Please make sure that the cell type names should be consistent with scRNA-Seq data files above.

* A **“LigandReceptor_Human.txt” or "LigandReceptor_Mouse.txt"** file listing all known ligand-receptor pairs. The first column (ligand) and the second column (receptor) are separated by a tab (\t). Currently, 1942 and 1855 ligand-receptor pairs are provided for human and mouse, respectively.

* A **“Species.txt”** file indicating the species from which the scRNA-Seq data are generated. Currently, “Human” and “Mouse” are supported.

* A **“Cutoff_GeneFilter.txt”** file indicating the cutoff for removing lowly-expressed genes in the processing of scRNA-Seq data. The default cutoff value is 0.1, which means that genes expressed in less than 10% of all cells of a given type are removed.

* A **“BetaUpperLimit.txt”** file indicating the upper limit of the test values of the algorithm parameter β, which is inversely proportional to the total number of genes in a given cell-type pair after removing lowly-expressed genes in the processing of scRNA-Seq data. Based on preliminary tests, the upper limit of β value is suggested to be 100 (default) if the total number of genes in a given cell-type pair is above 10,000. However, if the total number of genes is below 5000, it is necessary to increase the upper limit of β value to 500.

⚠ Please download **"CytoTalk_package_v2.0.zip"**. All example input files are in the **/Input/** folder and should be customized and copied into the **/CytoTalk/** folder before running. The /CytoTalk/ folder can only be used **ONCE** for a given cell-type pair. Please use a new /CytoTalk/ folder for analysis of other cell-type pairs.

## Run CytoTalk  

Copy the **input file-added “/CytoTalk/”** folder to your working directory and execute the following script:

```Bash
bash InferSignalingNetwork.sh
```

**[Alternative way]** The whole computation above may take 5.5 hours (2.3 GHz 8-Core Intel Core i9, 14 logical cores for parallel computation), of which 4 hours are used for computing pair-wise mutual information between genes in the construction of intracellular networks for the given two cell types. Considering that users may have alternative ways for constructing cell-type-specific intracellular networks, we divide the whole computation into two steps below.

```Bash
bash InferIntracellularNetwork_part1.sh  # around 4 hours
```

```Bash
bash InferIntercellularNetwork_part2.sh  # around 1.5 hours
```

The outputs of the script "part1.sh" are two comma-delimited files "IntracellularNetwork_TypeA.txt" and "IntracellularNetwork_TypeB.txt", containing the adjacency matrices of two intracellular networks for the given two cell types, respectively. These two files are the inputs of the script "part2.sh", which can generate the final predicted signaling network.

## CytoTalk output  

The output folder, “/CytoTalk/IllustratePCSF/”, contains a network topology file and six attribute files that are ready for import into Cytoscape for visualization and further analysis of the predicted signaling network between the given two cell types.

| Network topology | Edge attribute | Node attribute |
| :------ | :----- | :------ |
| PCSF_edgeSym.sif | PCSF_edgeCellType.txt<br>PCSF_edgeCost.txt| PCSF_geneCellType.txt，PCSF_geneExp.txt<br>PCSF_genePrize.txt，PCSF_geneRealName.txt |

   
     
## Cite CytoTalk  


* Hu Y, Peng T, Gao L, Tan K. CytoTalk: *De novo* construction of signal transduction networks using single-cell transcriptomics data. *Science Advances* (2021). Accepted.

## Reference  

* Shannon P, et al. Cytoscape: a software environment for integrated models of biomolecular interaction networks. *Genome Research* 13, 2498-2504 (2003).

## Contact  

Kai Tan, tank1@email.chop.edu<br />


