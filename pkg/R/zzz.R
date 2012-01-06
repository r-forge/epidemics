.First.lib <- function (lib, pkg){
#.initAdegenetClasses()
#.initAdegenetUtils()
    library.dynam("epidemics", pkg, lib)
    startup.txt <- "   ===========================\n    epidemics 1.0-0 is loaded\n   ===========================\n\n"

    packageStartupMessage(startup.txt)
}
