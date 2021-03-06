'.onLoad' <- function(libname, pkgname="lava") {
    addhook("heavytail.init.hook","init.hooks")
    addhook("glm.estimate.hook","estimate.hooks")
    addhook("ordinal.estimate.hook","estimate.hooks")
    addhook("cluster.post.hook","post.hooks")
    addhook("ordinal.sim.hook","sim.hooks")
    addhook("color.ordinal","color.hooks")
    addhook("ordinal.remove.hook","remove.hooks")
}

'.onAttach' <- function(libname, pkgname="lava") {
    #desc <- utils::packageDescription(pkgname)
    #packageStartupMessage(desc$Package, " version ",desc$Version)
    lava.options(cluster.index=versioncheck("mets",c(0,2,7)),
                 tobit=versioncheck("lava.tobit",c(0,5)))
}
