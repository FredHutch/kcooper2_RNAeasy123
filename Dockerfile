# build me as fredhutch/rnaeasy123
FROM fredhutch/r-shiny-server-base:3.6.0

RUN R -e "install.packages('BiocManager')"
RUN R -e "BiocManager::install(c('readr','limma','Glimma','edgeR','Mus.musculus','RColorBrewer','gplots','msigdbr','DT'))"

EXPOSE 3838

RUN rm -rf /srv/shiny-server/
ADD . /srv/shiny-server/01_hello

CMD /usr/bin/shiny-server.sh
