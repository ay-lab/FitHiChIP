FROM r-base:4.1.0

RUN apt-get update
RUN apt-get -y upgrade
RUN apt-get -y install pkg-config libssl-dev curl libcurl4-openssl-dev

# Install required R packages
RUN R -e "install.packages(c('BiocManager', 'optparse', 'ggplot2', 'data.table', 'splines', 'fdrtool', 'parallel', 'tools', 'plyr', 'dplyr'), quitely=TRUE, repos='http://cran.rstudio.com/')"

RUN R -e 'BiocManager::install()'
RUN R -e "BiocManager::install('GenomicRanges');"
RUN R -e "BiocManager::install('edgeR');"

# Install python 3
RUN apt-get -y install python3
RUN apt-get -y install python-pip
RUN pip install networkx

# Install Bedtools
RUN wget https://github.com/arq5x/bedtools2/releases/download/v2.30.0/bedtools-2.30.0.tar.gz
RUN tar -zxvf bedtools-2.30.0.tar.gz
RUN cd bedtools2 && make && cp -r ./bin/* /usr/local/bin/


# Install Samtools and Htslib
RUN cd / && wget https://github.com/samtools/samtools/releases/download/1.16.1/samtools-1.16.1.tar.bz2 && wget https://github.com/samtools/htslib/releases/download/1.16/htslib-1.16.tar.bz2
RUN tar -xjf /samtools-1.16.1.tar.bz2
RUN tar -xjf /htslib-1.16.tar.bz2

RUN cd /htslib-1.16 && ./configure --prefix=/usr/local/ && make && make install
RUN cd /samtools-1.16.1 && ./configure --prefix=/usr/local/ && make && make install

# Install miniconda
RUN wget https://repo.anaconda.com/miniconda/Miniconda2-latest-Linux-x86_64.sh
RUN sh Miniconda2-latest-Linux-x86_64.sh -b -f

ENV PATH="/root/miniconda2/bin/:${PATH}"

# Add bioconda
RUN conda config --add channels defaults
RUN conda config --add channels bioconda
RUN conda config --add channels conda-forge

# Install bowtie2
RUN conda install bowtie2


# Install HiCPro dependencies
RUN pip install pysam bx-python numpy scipy
RUN conda install -y -c anaconda scipy 
RUN conda install -y -c anaconda numpy 
RUN conda install -y -c bcbio bx-python 
RUN conda install -y -c bioconda pysam 

RUN R -e "install.packages(c('RColorBrewer'), quitely=TRUE, repos='http://cran.rstudio.com/')"

# Install HiCPro
RUN cd / && wget https://github.com/nservant/HiC-Pro/archive/refs/tags/v3.1.0.tar.gz
RUN tar -zxvf v3.1.0.tar.gz

RUN cd HiC-Pro-3.1.0/ && make configure && make install
ENV PATH="/HiC-Pro-3.1.0/bin/:${PATH}"
ENV PATH="/HiC-Pro-3.1.0/bin/utils/:${PATH}"

# Install Macs2
RUN pip install MACS2

# Get FitHiChIP
RUN apt-get -y install git
RUN cd / && git clone https://github.com/ay-lab/FitHiChIP
RUN cd /FitHiChIP && sed -i 's/\/home\/sourya\/packages\/HiCPro\/HiC-Pro_2.9.0\//\/HiC-Pro-3.1.0\//g' configfile_*
RUN pip install networkx


# Cleanup
RUN rm -rf /*tar*
RUN rm -rf /bedtools2/
RUN rm -rf /htslib-1.16/
RUN rm -rf /samtools-1.16.1/

RUN chmod -R 777 /root
RUN conda install -y iced 

RUN conda clean --all -y

RUN apt -y install bc

ENTRYPOINT /bin/bash

