# specifying base image
FROM ensemblorg/ensembl-vep

# set a directory for the app
WORKDIR /opt/vep/

# copy all the files to the container
COPY . .

# install dependencies
RUN chmod a+rwx /opt/vep/
RUN curl -O http://ftp.ensembl.org/pub/release-107/variation/vep/homo_sapiens_vep_107_GRCh37.tar.gz
RUN tar xzf homo_sapiens_vep_107_GRCh37.tar.gz
