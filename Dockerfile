FROM rocker/shiny
RUN apt-get update && apt-get install -y libxml2-dev

# Install DNAlandscapeR
COPY . /srv/shiny-server/DNAlandscapeR

# Install packrat packages
RUN cd /srv/shiny-server/DNAlandscapeR && \
    Rscript -e 'install.packages("packrat"); \
                packrat::restore()'

# Temporary permissions hack
RUN chmod -R 777 /srv/shiny-server/DNAlandscapeR

# Start and expose shiny server
EXPOSE 3838
CMD ["/usr/bin/shiny-server.sh"]

