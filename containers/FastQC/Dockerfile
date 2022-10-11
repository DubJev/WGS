# specifying base image
FROM ubuntu:22.04

# set a directory for the app
WORKDIR /usr/src/app

# copy all the files to the container
COPY . .

# install dependencies
RUN apt-get update && \
    apt-get install -yqq samtools bwa fastqc

# specify the port number that needs to be exposed
EXPOSE 5000

# command for running the application

