## get package description
conda skeleton cran https://github.com/HCGB-IGTP/XICRA.stats

## build package
conda build r-xicra.stats

## upload to anaconda
anaconda upload $miniconda3_path/conda-bld/linux-64/r-xicra.stats-0.0.1-r35_0.tar.bz2
