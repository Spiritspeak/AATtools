#Generate help files
devtools::document()
#Build manual
devtools::build_manual(path=".")
#build pkgdown manual
devtools::build_site()
file.copy(list.files("C:/Users/b1066151/Nextcloud/AAT/multilevel/package/AATtools/docs/",full.names = TRUE),
          to="Y:\\Sercan\\homepage\\spiritspeak.github.io\\content\\packages\\",recursive = TRUE)
#Check package
devtools::check(args="--as-cran")
#devtools::check_win_devel()
#build package
devtools::build()


