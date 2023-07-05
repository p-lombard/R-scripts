#https://github.com/FertigLab/SpaceMarkers
#https://www.biorxiv.org/content/10.1101/2022.06.02.490672v1.full
devtools::install_github("FertigLab/CoGAPS")
install.packages("remotes")
remotes::install_github("FertigLab/SpaceMarkers", dependencies = TRUE, build_vignettes = TRUE)