
FROM bioconductor/shiny

# system libraries of general use
## install debian packages
RUN apt-get update -qq && apt-get -y --no-install-recommends install \
    libxml2-dev \
    libcairo2-dev \
    libsqlite3-dev \
    libmariadbd-dev \
    libpq-dev \
    libssh2-1-dev \
    unixodbc-dev \
    libcurl4-openssl-dev \
    libssl-dev

## update system libraries
RUN apt-get update && \
    apt-get upgrade -y && \
    apt-get clean
    

RUN R -e "install.packages(pkgs=c('shiny','tidyverse','shinydashboard', 'DT', 'stringr', 'dplyr', 'shinycssloaders'), repos='https://cran.rstudio.com/')" 

RUN R -e 'BiocManager::install(ask = F)' && R -e 'BiocManager::install(c("Biostrings", \
    "BiocGenerics", "BSgenome.Hsapiens.UCSC.hg38", "Rsamtools", ask = F))'

RUN mkdir /home/shiny_app

# Copy the Shiny app code
COPY *.R CDS-human.csv /home/shiny_app/

EXPOSE 3838

# clean up
# RUN rm -rf /tmp/downloaded_packages/ /tmp/*.rds

# RUN dos2unix /usr/bin/shiny-server.sh && apt-get --purge remove -y dos2unix && rm -rf /var/lib/apt/lists/*
CMD ["R", "-e", "shiny::runApp('/home/shiny_app', host='0.0.0.0', port=3838)"]
