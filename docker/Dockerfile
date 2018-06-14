FROM ubuntu:trusty

# get linux fundamentals
RUN apt-get update \
    && apt-get -y install git python2.7 python-pip build-essential libboost-all-dev wget vim

# install required python modules
RUN pip install pgen cython pysam==0.9 Biopython==1.62

# clone the repo
RUN echo "alias python=\"python2.7\"" >> ~/.bashrc \
    && git clone https://github.com/teng-gao/BreaKmer.git ~/breaKmer \
    && mv ~/breaKmer/bin/blat ~/breaKmer/bin/faToTwoBit \
    ~/breaKmer/bin/gfClient ~/breaKmer/bin/gfServer /bin \
    && chmod +x bin/blat /bin/faToTwoBit \
    /bin/gfClient /bin/gfServer

# install cutadapt
RUN cd ~/breaKmer/bin \
    && tar -xzvf cutadapt-1.9.1.tar.gz

# install jellyfish
RUN cd ~/breaKmer/bin \
    && tar -xzvf jellyfish-1.1.11.tar.gz \
    && cd jellyfish-1.1.11 \
    && ./configure \
    && make \
    && make install

# put jellyfish in path
ENV LD_LIBRARY_PATH /root/breaKmer/bin/jellyfish-1.1.11/.libs

# prepare reference fasta
RUN mkdir /root/breaKmer/ref \
    && cd /root/breaKmer/ref \
    && wget http://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/chromFa.tar.gz \
    && tar -xzf chromFa.tar.gz \
    && cat $(ls | egrep "chr[0-9,X,Y,M][0-9]?.fa" | sort -V) > hg19_ref_genome.primary-assembly.fa \
    && sed -i 's/chr//g' hg19_ref_genome.primary-assembly.fa

# set up breakmer config file
COPY breakmer.cfg /root/breaKmer/breakmer.cfg

# get gene model
COPY ucsc_hg19_refgene.txt /root/breaKmer/ref/ucsc_hg19_refgene.txt

RUN echo "alias ll='ls -alFrt'" >> /root/.bashrc
