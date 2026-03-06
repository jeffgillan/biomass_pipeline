FROM rocker/r-ver:4.4

# Install geospatial system libraries
RUN apt-get update && apt-get install -y \
    libgdal-dev libgeos-dev libproj-dev \
    libudunits2-dev libfontconfig1-dev \
    && rm -rf /var/lib/apt/lists/*

# Install R packages
RUN R -e "install.packages(c('lidR','terra','sf','itcSegment','BIOMASS','ForestTools','data.table','future'), repos='https://cloud.r-project.org')"

WORKDIR /pipeline
COPY R/ R/

ENTRYPOINT ["Rscript", "R/run_pipeline.R"]
