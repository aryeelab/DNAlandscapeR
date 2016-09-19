### R code from vignette source 'siggenes.Rnw'

###################################################
### code chunk number 1: siggenes.Rnw:60-61
###################################################
library(siggenes)


###################################################
### code chunk number 2: siggenes.Rnw:64-65
###################################################
library(siggenes)


###################################################
### code chunk number 3: siggenes.Rnw:74-75
###################################################
data(golub)


###################################################
### code chunk number 4: siggenes.Rnw:149-150
###################################################
args(sam)


###################################################
### code chunk number 5: siggenes.Rnw:155-156
###################################################
args(d.stat)


###################################################
### code chunk number 6: siggenes.Rnw:190-191
###################################################
args.sam(summary)


###################################################
### code chunk number 7: siggenes.Rnw:264-266
###################################################
n <- 10
rep(1, n)


###################################################
### code chunk number 8: siggenes.Rnw:280-282
###################################################
n1 <- n2 <- 5
rep(c(0, 1), c(n1, n2))


###################################################
### code chunk number 9: siggenes.Rnw:308-310
###################################################
K <- 5
c((-1:-5), 1:5)


###################################################
### code chunk number 10: siggenes.Rnw:322-324
###################################################
K <- 5
rep(1:K, e = 2) * rep(c(-1 ,1), K)


###################################################
### code chunk number 11: siggenes.Rnw:335-337
###################################################
K <- 5
cbind(rep(c(-1, 1), 5), rep(1:5, e = 2))


###################################################
### code chunk number 12: siggenes.Rnw:374-376
###################################################
sam.out <- sam(golub, golub.cl, rand = 123, gene.names = golub.gnames[,3])
sam.out


###################################################
### code chunk number 13: siggenes.Rnw:400-401
###################################################
sam.out2 <- sam(golub, golub.cl, method = wilc.stat, rand = 123)


###################################################
### code chunk number 14: siggenes.Rnw:406-407
###################################################
summary(sam.out)


###################################################
### code chunk number 15: siggenes.Rnw:420-421
###################################################
print(sam.out, seq(1.5, 2.4, 0.1))


###################################################
### code chunk number 16: siggenes.Rnw:470-472
###################################################
sum.sam.out <- summary(sam.out, 3.3)
sum.sam.out


###################################################
### code chunk number 17: siggenes.Rnw:485-486
###################################################
print(sum.sam.out, varNames = "Proteins")


###################################################
### code chunk number 18: siggenes.Rnw:494-495
###################################################
sum.sam.out@row.sig.genes


###################################################
### code chunk number 19: siggenes.Rnw:500-501
###################################################
sum.sam.out@mat.fdr


###################################################
### code chunk number 20: siggenes.Rnw:506-507
###################################################
sum.sam.out@mat.sig


###################################################
### code chunk number 21: siggenes.Rnw:512-513
###################################################
list.siggenes(sam.out, 3.3)


###################################################
### code chunk number 22: siggenes.Rnw:520-521
###################################################
findDelta(sam.out, fdr = 0.05)


###################################################
### code chunk number 23: siggenes.Rnw:534-535
###################################################
findDelta(sam.out, genes = 200)


###################################################
### code chunk number 24: siggenes.Rnw:553-554
###################################################
find.out <- find.a0(golub, golub.cl, rand = 123)


###################################################
### code chunk number 25: siggenes.Rnw:560-561
###################################################
find.out


###################################################
### code chunk number 26: siggenes.Rnw:569-570
###################################################
print(find.out, 0.95)


###################################################
### code chunk number 27: siggenes.Rnw:604-605
###################################################
ebam(find.out)


###################################################
### code chunk number 28: siggenes.Rnw:612-613
###################################################
ebam(find.out, which.a0 = 2)


###################################################
### code chunk number 29: siggenes.Rnw:625-626
###################################################
ebam(golub, golub.cl, a0 = 0, fast = TRUE, rand = 123)


###################################################
### code chunk number 30: siggenes.Rnw:640-641
###################################################
ebam.out <- ebam(golub, golub.cl, a0 = 0, rand = 123)


###################################################
### code chunk number 31: siggenes.Rnw:651-652
###################################################
print(ebam.out, seq(0.91, 0.99, 0.01))


###################################################
### code chunk number 32: siggenes.Rnw:675-676
###################################################
summary(ebam.out, 0.99997)


###################################################
### code chunk number 33: siggenes.Rnw:693-694
###################################################
ebam(golub, golub.cl, a0 = 0, var.equal = TRUE, rand = 123)


###################################################
### code chunk number 34: siggenes.Rnw:701-702
###################################################
ebam(golub, golub.cl, quan.a0 = 0.5, rand = 123)


###################################################
### code chunk number 35: siggenes.Rnw:709-710
###################################################
ebam(golub, golub.cl, method = wilc.ebam, rand =123)


###################################################
### code chunk number 36: siggenes.Rnw:777-795
###################################################
t.stat <- function(data, cl){
    require(genefilter) ||
        stop("genefilter required.")
    cl <- as.factor(cl)
    row.out <- rowttests(data, cl)
    d <- row.out$statistic
    m <- length(na.exclude(d))
    d.bar <- qt(((1:m) - 0.5)/m, length(cl) - 2)
    p.value <- row.out$p.value
    vec.false <- m * p.value/2
    s <- row.out$dm/d
    msg <- paste("SAM Two-Class Analysis",
         "Assuming Normality\n\n")
    list(d = -d, d.bar = d.bar, p.value = p.value,
        vec.false = vec.false, s = s, s0 = 0,
        mat.samp = matrix(numeric(0)),
        msg = msg, fold = numeric(0))
}


###################################################
### code chunk number 37: siggenes.Rnw:808-809
###################################################
sam(golub, golub.cl, method = t.stat)


###################################################
### code chunk number 38: siggenes.Rnw:854-871
###################################################
t.find <- function(data, cl, B = 50){
    require(genefilter)
    z.fun <- function(data, cl){
        cl <- as.factor(cl)
        out <- rowttests(data, cl)
        r<- out$dm
        s<- r / out$statistic
        return(list(r = -r, s = s))
    }
    mat.samp <- matrix(0, B, length(cl))
    for(i in 1:B)
        mat.samp[i, ] <- sample(cl)
    z.out <- z.fun(data, cl)
    msg <- paste("EBAM Analysis with a Moderated t-Statistic\n\n")
    list(r = z.out$r, s = z.out$s,
        mat.samp = mat.samp, z.fun = z.fun, msg = msg)
}


###################################################
### code chunk number 39: siggenes.Rnw:881-883
###################################################
t.out <- find.a0(golub, golub.cl, method = t.find, B = 100, rand =123)
t.out


###################################################
### code chunk number 40: siggenes.Rnw:888-889
###################################################
find.a0(golub, golub.cl, var.equal = TRUE, rand =123)


###################################################
### code chunk number 41: siggenes.Rnw:894-895
###################################################
ebam(t.out)


###################################################
### code chunk number 42: siggenes.Rnw:928-941
###################################################
t.ebam<-function(data, cl){
    require(genefilter)
    cl <- as.factor(cl)
    out <- rowttests(data, cl)
    z <- -out$statistic
    z.dens <- denspr(z)$y
    m <- length(z)
    vec.pos <- m * out$p.value / 2
    z.null <- dt(z, length(cl) - 2)
    msg<-paste("EBAM Analysis with t-Statistic Assuming Normality.\n\n")
    list(z = z, ratio = z.null/z.dens, vec.pos = vec.pos,
        vec.neg = vec.pos, msg = msg)
}


###################################################
### code chunk number 43: siggenes.Rnw:951-952
###################################################
ebam(golub, golub.cl, method = t.ebam)


