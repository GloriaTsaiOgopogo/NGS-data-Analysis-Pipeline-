pierre.unpair.combine.2 =
function (h, controls, expts, num.permutations=100, winsize, conf, mat.sample=NULL, minrep)
{
  num.datasets <- length(h)
  # calculate the d statistics for all the datasets in h without permutation
  d.orig <- pierre.unpair.combine.elt.2(h,controls,expts,winsize,conf,minrep)
  d.stat <- d.orig$stat

  # now do the permutations (same permutation to each dataset).
  d.perm <- c()
  if(is.null(mat.sample)) {
    mat.sample <- sample.matrix(x=length(controls),y=length(expts),B=num.permutations,paired=FALSE)
  }
  for(i in 1:nrow(mat.sample)) {
    tmp.controls <- which(mat.sample[i,] == -1)
    tmp.expts <- which(mat.sample[i,] == 1)
    d <- pierre.unpair.combine.elt(h,controls=tmp.controls,expts=tmp.expts,minrep=minrep,winsize=winsize,conf=conf)
    d.perm <- cbind(d.perm,d$stat)
    cat("Done with ",i,", Doing permutation\n",sep="")
    print(mat.sample[i,])
  }
  q.val <- fdr.cal(d.stat,as.vector(d.perm))
  d.orig <- c(d.orig,list(q=q.val,dstat=d.stat))
  return(list(d.orig=d.orig,mat.sample=mat.sample,d.perm=d.perm,q=q.val))
}


pierre.unpair.combine.elt.2=
function (h, controls, expts, winsize, conf, minrep) 
{
    num.datasets <- length(h)
    all.n <- c()
    all.df <- c()
    all.var <- c()
    all.mean <- c()
    all.d <- c()
    all.sign <- c()
	ctrl.mean = c()
	expt.mean = c()
   for (i in 1:num.datasets) {
        X <- as.matrix(h[[i]][, c(controls, expts)])
        mode(X) <- "numeric"
        tmp.ctls <- c(1:length(controls))
        tmp.expts <- (length(controls) + 1):(length(expts) + 
            length(controls))
        d.orig <- pierre.unpair.element(X, tmp.ctls, tmp.expts, 
            minrep = minrep, conf = conf, winsize = winsize)
        all.d <- cbind(all.d, d.orig$Bay.t)
        all.df <- cbind(all.df, (d.orig$C.N + d.orig$E.N - 2))
        all.n <- cbind(all.n, (d.orig$C.N + d.orig$E.N))
        tmp.var <- sqrt(((d.orig$C.N - 1) * d.orig$C.var * d.orig$C.var + 
            (d.orig$E.N - 1) * d.orig$E.var * d.orig$E.var)/(d.orig$C.N + 
            d.orig$E.N - 2))
        all.var <- cbind(all.var, tmp.var)
        tmp.mean <- d.orig$E.mean - d.orig$C.mean
        ctrl.mean=	cbind(ctrl.mean, d.orig$C.mean) # to get all the means together in one table for part II of the analysis
        expt.mean= cbind(expt.mean, d.orig$E.mean)
        all.mean <- cbind(all.mean, d.orig$E.mean - d.orig$C.mean)
    }
    colnames(all.d) <- names(h)
    colnames(all.df) <- names(h)
    colnames(all.n) <- names(h)
    colnames(all.var) <- names(h)
    colnames(all.mean) <- names(h)
    all.n[is.na(all.d)] <- NA
    all.df[is.na(all.d)] <- NA
    all.mean[is.na(all.d)] <- NA
    all.var[is.na(all.d)] <- NA
 	row.mean.ctrl=apply(ctrl.mean,1,mean.na)
	row.mean.expt=apply(expt.mean,1,mean.na)
	fraction.mean=cbind(row.mean.ctrl, row.mean.expt)
	colnames(fraction.mean)=c("sum.mean.ctrl","sum.mean.treated")
    combined.mean <- apply(all.n * all.mean, 1, sum.na)/apply(all.n, 
        1, sum.na)
    var.a <- apply(all.df * all.var * all.var, 1, sum.na)
    var.b <- apply(all.n * (all.mean - combined.mean) * (all.mean - 
        combined.mean), 1, sum.na)
    stat.tmp <- combined.mean/sqrt(var.a + var.b)
    return(list(all.d = all.d, all.N = all.n, all.var = all.var, 
        all.mean = all.mean, var.a = var.a, var.b = var.b, combined.mean = combined.mean, 
        stat = stat.tmp, mean.ctrl=fraction.mean[,1], mean.expt=fraction.mean[,2]))
}

pierre.unpair.element.2=
function (h, controls, expts, winsize, conf, minrep) 
{
    C.N <- apply(h[, controls], 1, function(x) sum(x != 0 & !is.na(x))) # count the number on not-NA in the array
    E.N <- apply(h[, expts], 1, function(x) sum(x != 0 & !is.na(x)))	  # count the number on not-NA in the array
    C.mean <- apply(h[, controls], 1, function(x) if (sum(x[x !=        # calculate the mean of each triplet of experiment groups (contorl)
        0 & !is.na(x)] != 0)) 
        mean(x[x != 0 & !is.na(x)])
    else NA)
    E.mean <- apply(h[, expts], 1, function(x) if (sum(x[x !=           # calculate the mean of each triplet of experiment groups (treatment)
        0 & !is.na(x)] != 0)) 
        mean(x[x != 0 & !is.na(x)])
    else NA)
    C.sd <- apply(h[, controls], 1, function(x) if (sum(x[x !=          # calculate the SD of each triplet of experiment groups (contorl)
        0 & !is.na(x)] != 0) > 1) 
        sqrt(var(x[x != 0 & !is.na(x)]))
    else NA)
    E.sd <- apply(h[, expts], 1, function(x) if (sum(x[x != 0 &         # calculate the SD of each triplet of experiment groups (treatment) 
        !is.na(x)] != 0) > 1) 
        sqrt(var(x[x != 0 & !is.na(x)]))
    else NA)
    xxx <- C.sd[!is.na(C.sd)][order(C.mean[!is.na(C.sd)])]
    xxx <- runavg(xxx, winsize)
    xxx <- xxx[rank(C.mean[!is.na(C.sd)])]
    C.rasd <- rep(NA, nrow(h))
    C.rasd[!is.na(C.sd)] <- xxx
    xxx <- E.sd[!is.na(E.sd)][order(E.mean[!is.na(E.sd)])]
    xxx <- runavg(xxx, winsize)
    xxx <- xxx[rank(E.mean[!is.na(E.sd)])]
    E.rasd <- rep(NA, nrow(h))
    E.rasd[!is.na(E.sd)] <- xxx
    C.var <- sqrt((conf * C.rasd^2 + (C.N - 1) * C.sd^2)/(conf +        # calculate the variance of the experiment groups 
        C.N - 2))
    E.var <- sqrt((conf * E.rasd^2 + (E.N - 1) * E.sd^2)/(conf + 
        E.N - 2))
    Bay.result <- t(apply(cbind(C.N, E.N, C.mean, E.mean, C.var, 
        E.var), 1, function(x) tstat.general(x, 1, 2, 3, 4, 5, 
        6, 0.1, minrep)))
    colnames(Bay.result) <- c("t.log", "n.log", "vr.log")
    Bay.p <- log(2) + pt(abs(Bay.result[, 1]), df = Bay.result[, 
        2] + 2 * conf - 2, lower.tail = FALSE, log.p = TRUE)
    return(data.frame(C.mean = C.mean, E.mean = E.mean, C.sd = C.sd, 
        E.sd = E.sd, C.rasd = C.rasd, E.rasd = E.rasd, C.N = C.N, 
        E.N = E.N, C.var = C.var, E.var = E.var, Bay.vr = Bay.result[, 
            3], Bay.t = Bay.result[, 1], Bay.n = Bay.result[, 
            2], Bay.p = Bay.p))
}



