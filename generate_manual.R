#Generate help files
devtools::document()
#Build manual
devtools::build_manual(path=".")
#build pkgdown manual
pkgdown::build_site()
#Check package
devtools::check(args="--as-cran")
#devtools::check_win_devel()
#build package
devtools::build()
