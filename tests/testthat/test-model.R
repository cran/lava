context("Model specification")

test_that("Basic model building blocks", {
    m <- lvm(y[m]~x)
    covariance(m) <- y~z
    testthat::expect_true(covariance(m)$rel["z","y"]==1)
    testthat::expect_true(regression(m)$rel["x","y"]==1)

    ## Children parent,nodes
    testthat::expect_true(children(m,~x)=="y")
    testthat::expect_true(parents(m,~y)=="x")
    testthat::expect_equivalent(parents(m),vars(m))
    testthat::expect_equivalent(children(m),vars(m))

    ## Remove association
    cancel(m) <- y~z+x
    testthat::expect_true(covariance(m)$rel["z","y"]==0)
    testthat::expect_true(regression(m)$rel["x","y"]==0)

    ## Remove variable
    kill(m) <- ~x
    testthat::expect_equivalent(vars(m),c("y","z"))
    testthat::expect_true(intercept(m)["y"]=="m")

    m <- lvm(c(y1,y2,y3)~x)
    d <- sim(m,50)
    e <- estimate(m,d)
    ## Equivalence
    ##equivalence(e,silent=TRUE)


    ## formula
    f <- formula(m,all=TRUE)
    testthat::expect_true(length(f)==length(vars(m)))
    testthat::expect_true(all(unlist(lapply(f,function(x) inherits(x,"formula")))))

    ## Parametrization
    m <- lvm(c(y1,y2,y3)~u)
    latent(m) <- ~u
    m2 <- fixsome(m,param=NULL)
    testthat::expect_true(all(is.na(regression(m2)$values)))
    m2 <- fixsome(m,param="relative")
    testthat::expect_true(regression(m2)$values["u","y1"]==1)
    testthat::expect_true(intercept(m2)[["y1"]]==0)
    m2 <- fixsome(m,param="hybrid")
    testthat::expect_true(regression(m2)$values["u","y1"]==1)
    testthat::expect_true(intercept(m2)[["u"]]==0)
    m2 <- fixsome(m,param="absolute")
    testthat::expect_true(all(is.na(regression(m2)$values)))
    testthat::expect_true(intercept(m2)[["u"]]==0)
    testthat::expect_true(covariance(m2)$values["u","u"]==1)

    ## Merge
    m1 <- lvm(c(y1,y2,y3)~1*u1[m1:v1])
    latent(m1) <- ~u1
    m2 <- lvm(c(y1,y2,y3)~2*u2[m2:v2])
    latent(m2) <- ~u2
    mm <- m1%++%m2

    testthat::expect_true(covariance(mm)$labels["u1","u1"]=="v1")
    testthat::expect_true(intercept(mm)[["u2"]]=="m2")

    ## LISREL
    mm <- fixsome(mm)
    L <- lisrel(mm,rep(1,length(coef(mm))))
    testthat::expect_equivalent(L$B,matrix(0,2,2))
    testthat::expect_equivalent(L$Theta,diag(3))
    testthat::expect_equivalent(L$Psi,diag(2))

})


test_that("Linear constraints", {
    m <- lvm(c(y[m:v]~b*x))
    constrain(m,b~a) <- base::identity
    d <- sim(m,100,seed=1)
    l <- lm(y~x, d)
    e <- estimate(m, d)
    err <- sum((coef(l)-coef(e)[c('y','a')])^2)
    testthat::expect_true(err<1e-12)
})


if (requireNamespace("Rgraphviz",quietly = TRUE))
test_that("Graph attributes", {
    m <- lvm(y~x)
    suppressMessages(g1 <- graph::updateGraph(plot(m,noplot=TRUE)))
    m1 <- graph2lvm(g1)
    testthat::expect_equivalent(m1$M, m$M)

    col <- "blue"
    v <- "y"
    g1 <- lava::addattr(g1, "fill", v, col)
    testthat::expect_true(col == graph::nodeRenderInfo(g1)$fill[[v]])
    nodecolor(m, v) <- "blue"

    g2 <- Graph(m, add=TRUE)
    testthat::expect_true(inherits(g2, "graph"))
    testthat::expect_true(col == graph::nodeRenderInfo(g2)$fill[[v]])
    testthat::expect_true(addattr(g2, "fill")[[v]] == "blue")
    graph::graphRenderInfo(g2)$rankdir <- "LR"
    Graph(m) <- g2
    testthat::expect_true(graph::graphRenderInfo(Graph(m))$rankdir=="LR")

    ## Labels
    labels(m) <- c(y = "Y")
    addattr(Graph(m, add=TRUE), "label")
    testthat::expect_true(addattr(finalize(m), "label")[["y"]]=="Y")
    labels(g2) <- c(y = "Y")
    testthat::expect_true(!is.null(graph::nodeRenderInfo(g2)$label["y"]))

    edgelabels(m, y~x) <- "a"
    testthat::expect_true(!is.null(edgelabels(finalize(m))))
})


test_that("Categorical variables", {
    m <- lvm()
    categorical(m,K=3,p=c(0.1,0.5)) <- ~x
    d1 <- simulate(m,10,seed=1)
    categorical(m,K=3) <- ~x
    d2 <- simulate(m,10,seed=1)
    testthat::expect_false(identical(d1,d2))

    regression(m,additive=FALSE,y~x) <- c(0,-5,5)
    d <- simulate(m,100,seed=1)
    l <- lm(y~factor(x),d)
    testthat::expect_true(sign(coef(l))[2]==-sign(coef(l))[3])

})
