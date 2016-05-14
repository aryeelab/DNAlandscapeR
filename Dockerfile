FROM rocker/shiny
RUN apt-get update && apt-get install -y libxml2-dev

# Install DNAlandscapeR
COPY . /srv/shiny-server/DNAlandscapeR

# Install R package dependencies
RUN cd /srv/shiny-server/DNAlandscapeR && \
    Rscript -e 'source("https://bioconductor.org/biocLite.R"); \
             biocLite(c( "shiny", \
                "shinythemes", \
                "ggplot2", \
                "GenomicRanges", \
                "diffloop", \
                "rtracklayer", \
                "shinyFiles", \
                "DT", \
                "bumphunter"))'


# Temporary permissions hack
RUN chmod -R 777 /srv/shiny-server/DNAlandscapeR

# Start and expose shiny server
EXPOSE 3838
CMD ["/usr/bin/shiny-server.sh"]

