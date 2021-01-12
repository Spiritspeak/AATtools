#Generate help files
devtools::document()
#Build manual
devtools::build_manual(path=".")
#build pkgdown manual
devtools::build_site()
pkgdown::build_site(path="Y:\\Sercan\\homepage\\spiritspeak.github.io\\content\\packages")
#Check package
devtools::check(args="--as-cran")
#devtools::check_win_devel()
#build package
devtools::build()
