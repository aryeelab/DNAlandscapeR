FROM rocker/shiny
RUN apt-get update && apt-get install -y git=1:2.8.1-1 libxml2-dev=2.9.3+dfsg1-1

# Install git lfs
RUN wget https://github.com/github/git-lfs/releases/download/v1.2.0/git-lfs-linux-amd64-1.2.0.tar.gz && \
    tar zxf git-lfs-linux-amd64-1.2.0.tar.gz && \
    cd git-lfs-1.2.0 && \
    ./install.sh

# Install DNAlandscapeR
RUN cd /srv/shiny-server && \
    git clone --depth=1 https://github.com/aryeelab/DNAlandscapeR.git

WORKDIR /srv/shiny-server/DNAlandscapeR

# Install packrat packages
RUN Rscript -e 'install.packages("packrat"); \
                packrat::restore()'

# Temporary permissions hack
RUN chown shiny:shiny -R packrat && \
    chown shiny:shiny .gitignore 

# Start and expose shiny server
EXPOSE 3838
CMD ["/usr/bin/shiny-server.sh"]

