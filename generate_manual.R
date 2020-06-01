#Generate help files
devtools::document()
#Build manual
devtools::build_manual(path=".")
#Check package
devtools::check()
