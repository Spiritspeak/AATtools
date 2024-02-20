#Generate help files
devtools::document()
#Check package
devtools::check(args="--as-cran")
#Build manual
devtools::build_manual(path=".")
#build pkgdown manual
devtools::build_site()
file.copy(list.files("W:/Sercan/package/AATtools/AATtools/docs/",full.names = TRUE),
          to="W:\\Sercan\\homepage\\spiritspeak.github.io\\content\\packages\\",recursive = TRUE)
#devtools::check_win_devel()
#build package
devtools::build()


