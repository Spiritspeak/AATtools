#Generate help files
devtools::document()
#Build manual
devtools::build_manual(path=".")
#Check package
devtools::check(args="--as-cran")
#build package
devtools::build()
