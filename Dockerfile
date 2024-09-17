FROM mtbls/safer-base
RUN Rscript -e "install.packages('apcluster')"
RUN Rscript -e "install.packages('DEoptimR')"
RUN Rscript -e "install.packages('sass')"
RUN Rscript -e "install.packages('openssl')"
RUN Rscript -e "install.packages('remotes')"
RUN Rscript -e "devtools::install_github('EBI-Metabolights/SAFERnmr')"
CMD ["R"]
