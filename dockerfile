FROM r-base:4.1.3
# INSTALLATION
RUN apt-get update \
    && DEBIAN_FRONTEND=noninteractive apt-get install -y --no-install-recommends  build-essential libssl-dev libncurses5-dev libcurl4-openssl-dev liblzma-dev libbz2-dev libboost-all-dev sqlite3 libsqlite3-0 libsqlite3-dev libgsl0-dev zlib1g-dev libxml2-dev libgmp-dev libgmp10 libmpfr-dev cmake libnlopt-dev \
    && echo '# Install Bioconductor' \
    && R -e 'install.packages("BiocManager")' \
    && echo '## Use the argument ask=FALSE to update old packages without being prompted' \
    && R -e 'BiocManager::install(ask=FALSE)' \
    \
    \
    && echo '# Install R libraries used during the analysis' \
    && R -e 'BiocManager::install("BgeeDB")' \
    && R -e 'BiocManager::install("data.table")' \
    && R -e 'BiocManager::install("here")' \
    && R -e 'BiocManager::install("ggplot2")' \
    && R -e 'BiocManager::install("dplyr")' \
    && R -e 'BiocManager::install("mclust")' \
    && R -e 'BiocManager::install("RColorBrewer")' \
    && R -e 'BiocManager::install("reshape2")' \
    && R -e 'BiocManager::install("tidyr")' \
    && R -e 'BiocManager::install("gridExtra")' \
    && R -e 'BiocManager::install("VennDiagram")' \
    && echo '# Install R libraries from CRAN' \
    && echo '# First install dependencies: nloptr for ggpubr package' \
    && R -e 'install.packages("ggpubr")' \
    && R -e 'install.packages("UpSetR")' \
    && R -e 'install.packages("stringr")' \
    && R -e 'install.packages("skimr")' \
    && R -e 'install.packages("inflection")' \
    && R -e 'install.packages("stringi")'
# copy all directory (Methods_RNASeq_expression_calls) to docker destination
COPY data/ /Copy_all_repository/data/
COPY scripts/ /Copy_all_repository/scripts/
# Start with batch script by default
CMD ["/Copy_all_repository/scripts/call_all_Rscripts.sh"]
