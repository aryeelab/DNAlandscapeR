FROM rocker/shiny
RUN apt-get update && apt-get install -y libxml2-dev

# Install R package dependencies
RUN Rscript -e 'source("https://bioconductor.org/biocLite.R"); \
             biocLite(c( "shiny", \
                "shinythemes", \
                "ggplot2", \
                "GenomicRanges", \
                "diffloop", \
                "rtracklayer", \
                "shinyFiles", \
                "DT", \
                "bumphunter"))'

# Install DNAlandscapeR
RUN rm -fr /srv/shiny-server/*
COPY . /srv/shiny-server/

# Temporary permissions hack
RUN chmod -R 777 /srv/shiny-server

# Start and expose shiny server
EXPOSE 3838
CMD ["/usr/bin/shiny-server.sh"]

