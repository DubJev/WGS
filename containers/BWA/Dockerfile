# specifying base image
FROM ubuntu:22.04

# set a directory for the app
WORKDIR /usr/src/app

# copy all the files to the container
COPY . .

# install dependencies
RUN apt-get update
RUN apt-get install -yqq autoconf automake make gcc perl 
RUN apt-get install -yqq zlib1g-dev libbz2-dev liblzma-dev libcurl4-gnutls-dev libssl-dev libncurses5-dev
RUN apt-get install -yqq git
RUN git clone https://github.com/samtools/samtools
RUN cd samtools
RUN apt install -yqq samtools
RUN cd ..
RUN git clone https://github.com/lh3/bwa.git
RUN cd bwa
RUN apt install -yqq bwa

# specify the port number that needs to be exposed
EXPOSE 5000

# command for running the application

