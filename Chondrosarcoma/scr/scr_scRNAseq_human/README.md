# Chondrosarcoma single-cell RNA-seq data integration

M1 Internship project to integrate differents datasets of single-cell RNA-seq (https://github.com/learogue/scRNAseq_project) adapted to my M2 intership project. These scripts process 10X data (3 files : barcodes.tsv, genes.tsv, matrix.mtx) from dataset GSE184118. We use it to analyze conventional chondrosarcoma and observe the expression of CMH 1 and CTAs.

## Installation and Requirements

These scripts use multiple packages. For installing, you can use these command lines in your terminal (on vs code terminal for exemple) or you can use conda (see below).
They also use packages pywin32 (**only if you have Windows**).
Use Python version >= 3.11

```
pip install scanpy==1.10.2
pip install igraph==0.10.13
pip install tkfilebrowser==2.3.2
pip install pywin32==306
```

>[!NOTE]
>It's possible that when you execute `main.py`, an error like `_tkinter.TclError: Item ... already exists`, you can refere here : https://github.com/hwstar/BOMtools/issues/2


## Contents

- `apply_filters.py`: function to apply filters on anndata object on the minimum of genes per cell, the minimum of cells per genes and the maximum of % mitochondrial genes, generate a violin plot and ask the user if he want to save the filtered object
- `create_anndata_object.py`: function to create anndata object to save data and generate a violin plot
- `create_umaps.py`: function to create UMAPs from a merged object
- `main.py`: main script which call function and interact with the user, use `apply_filters.py`, `create_anndata_object.py`, `merge.py`
- `merge.py`: function to merge multiple anndata objects

## Data Processing Script

This script interacts with the user to process and analyze single-cell RNA sequencing datasets (10X). Follow the process to create objects, filter data and merge objects.

## Usage
If you use VS code, you can Open Folder to be in the right folder. Execute `main.py` script and follow the processus to select the desired operation.

```
python main.py
```

## Main Menu Options

At the beginning, a window appear and show folders in Processing. Here, you can create a folder or select an existing folder to save processing files for a better tracability. (More details in the tutorial folder)

1. **Create objects**
3. **Apply filters**
4. **Merge objects**

### 10X format objects creation

- Choose one or more folders for process 10X
- Creates objects using the `create_anndata_object` function (similar to Smart-seq2).

### Filtering Data

- Choose one or more objects to add filters
- Requests filter parameters (minimum genes per cell, minimum cells per gene, and maximum percentage of mitochondrial genes).
- Executes the `apply_filters` function from the `apply_filters.py` script.
- Displays the number of cells and genes before and after filtering 
- Asks if the user wants to save the filtered object (e.g. `object_1_Deng_filtered_1.h5ad`) in `Objects/Objects_ori` and saves a violin plot with the filter parameters in the name. Moreover, a log file `filters_appliled.tab` is created with all the parameters choosen to a better tracability.

### Merge Objects

- Displays objects folder
- Choose one or more objects in `Objects/Objects_filtered` to merge (the best is to merge one object by one, not an already merge with another because the script take all genes and some genes can miss if you merged an already merge object with an other moreover it gives errors because of 0)
- Merges objects using the `merge` function from the `merge.py` script and saves the merged object in `Objects/Objects_merged` (e.g. `object_merged_1.h5ad`) and a log file `objects_merged.tab` with the object number and all the object names merged.
