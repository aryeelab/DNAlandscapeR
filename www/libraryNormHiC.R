
# Function that normalizes Hi-C data based on library complexity

.libraryNormHiC <- function(losm){
    options(scipen=999)
    
    # Set up long matrix to get the differences
    idx <- lapply(losm, function(m) summary(m)[,c(1,2)])
    all <- unique(Reduce("rbind", idx))
    counts <- sapply(losm, function(m){ as.matrix(m[cbind(all$i, all$j)])})
    long <- data.matrix(cbind(as.numeric(colnames(losm[[1]])[all$i]), as.numeric(colnames(losm[[1]])[all$j]), counts))
    diff <- abs(long[,1] - long[,2])
    longdiff <- cbind(long, diff)
    colnames(longdiff) <- c("idx1", "idx2", paste0("s", seq(1, length(losm), 1)), "diff")

    # Infer resolution and max size
    res <- min(diff[diff > 0])
    ma <- res * (dim(losm[[1]])[1]-1)
    
    # Aggregate and append means
    m <- aggregate(x = longdiff, by = list(longdiff[,3+length(losm)]), FUN = "mean")
    scaled <- m[,4:(3 + length(losm))]/rowMeans(m[,4:(3 + length(losm))])
    mlo <- data.frame(cbind(diff = m[,4 + length(losm)], scaled))
    ldm <- merge(longdiff, mlo, by.x = c("diff"), by.y = c("diff"))
    
    # Make long zeros matrix
    bins <- seq(0, ma, res)
    zeros.long <- cbind(t(combn(bins, 2)), 0)
    zeros.long <- rbind(zeros.long, cbind(bins, bins, 0))
    colnames(zeros.long) <- c("idx1", "idx2", "val")
    
    # Apply transform and make new list
    dat <- lapply(1:(length(losm)), function(i){
        a <- data.frame(cbind(idx1=ldm[,2], idx2=ldm[,3], val=ldm[,i+3]/ldm[,i+3+length(losm)]))
        Matrix(reshape2::acast(rbind(a, zeros.long), formula = idx1 ~ idx2, value.var = "val", fill = 0, fun.aggregate = sum))
    })
    return(dat)
}