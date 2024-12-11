# # psuedo EBLUP estimation of domain means and model parameters
pseudo_EBLUP <- function(ys,ws,Xs,Xmean,areas,sigmae2est,sigmau2est,nd){

    D <- length(unique(areas))
    mean_PEBLUP<-rep(0,D)
    p2 <- 0; p4 <- 0

    # Calculations needed for pseudo-EBLUP 
    p <- ncol(Xs)
    lambda <- ifelse(all(Xs[,1] == 1), 1, 0)
    wsum <-rep(0,D)
    meanydsw<-rep(0,D)
    meanXsw<-matrix(0,nr=D,nc=p)
    meanXs<-matrix(0,nr=D,nc=p)
    meanXsCons<-0
    delta2d<-rep(0,D)
    gammadw<-rep(0,D)
    nums0 <- matrix(0,nr=p,nc=1)
    dens0 <- matrix(0,nr=p,nc=p)

    for (d in 1:D){    
        # Domain means of y 
        yds <-ys[areas==d]
        wds<-ws[areas==d]
        wsum[d]<-sum(wds)
        Xds<-Xs[areas==d,] |> matrix(ncol = p)

        meanydsw[d] <- sum(wds*yds)/wsum[d] 
        for (k in 1:p){
            meanXsw[d,k] <- sum(wds*Xds[,k])/wsum[d]
            meanXs[d,k] <- sum(Xds[,k])/nd[d]
        }

        delta2d[d]<-sum(wds^2)/wsum[d]^2
        gammadw[d]<-sigmau2est/(sigmau2est+sigmae2est*delta2d[d])
        dwds<-diag(wds, nrow = length(wds))
        yds_ast <- yds-gammadw[d]*meanydsw[d]
        Xds_Wds <- t(Xds)%*%dwds
        Xds_ast <- Xds-matrix(rep(gammadw[d]*meanXsw[d,],nd[d]),nr=nd[d],byrow=TRUE)

        nums0 <- nums0+(Xds_Wds%*%yds_ast)
        dens0 <- dens0+(Xds_Wds%*%Xds_ast)

        # Auxiliary computes for MSE estimator
        zd <- (dwds %*% Xds_ast)
        p2 <- p2 + (t(zd) %*% zd)
        p4 <- p4 + (colSums(zd) %*% t(colSums(zd)))

        meanXsCons <- meanXsCons + (meanXs[d,] %*% t(meanXs[d,])*nd[d]^2)
    }

    betaw <- solve(dens0)%*%nums0
    udw<-gammadw*(meanydsw-meanXsw%*%betaw)
    mean_PEBLUP<-Xmean%*%betaw+udw[,1]
 
    return(list(mean_PEBLUP = mean_PEBLUP, betaw = betaw))
}

# Parametric Bootstrap (PB) of MSE for PEBLUP
mse_PB <- function(sam, beta.est, sigmau2.est, sigmae2.est, df_pop_means,
                    PopnSegments, Xs, nd=samsize_dom$nd, formula.mod, weights.calib ){
  
  y.pred <- (Xs %*% beta.est)[,1]
  sam.boot <- sam |> 
    group_by(dom) |> 
    mutate(ed = rnorm(n(), 0, sqrt(sigmae2.est)),
           ud = rnorm(1, 0, sqrt(sigmau2.est))) |> 
    ungroup()
  sam.boot$y <- y.pred + sam.boot$ud + sam.boot$ed
  Xmean <- cbind(1:nrow(df_pop_means), df_pop_means[,-1])

  form.temp <- formula(paste0('y~', formula.mod))
  mod.est.boot <- lmer(form.temp, data = sam.boot)
  beta.est.boot <- fixef(mod.est.boot)
  sigmau2.est.boot <- as.numeric(VarCorr(mod.est.boot)); if(sigmau2.est.boot==0) sigmau2.est.boot <- 1e-16
  sigmae2.est.boot <- summary(mod.est.boot)$sigma^2
  
  ypred.boot <- pseudo_EBLUP(ys = sam.boot$y, ws = weights.calib,
      areas = sam$dom,
      Xs = Xs, 
      Xmean = df_pop_means, 
      sigmau2est = sigmau2.est.boot, 
      sigmae2est = sigmae2.est.boot, 
      nd=nd)

  ud.boot <- sam.boot |> select(dom, ud) |> unique() |> pull(ud)

  mse_PB <- tibble(dom = unique(sam$dom), 
    REAL.boot = (df_pop_means %*% beta.est)[,1] + ud.boot,
    est.PEBLUP.cal.boot = ypred.boot$mean_PEBLUP[,1]) |> 
    mutate(mse.PB = (est.PEBLUP.cal.boot - REAL.boot)^2) |> 
    select(dom, mse.PB)

  return(mse_PB)
}

# Parametric Bootstrap (PB) of MSE for UA
mse_FHT <- function(df_estimates, beta, sigmau2, Xmean, psid.est.var, y_var_dir){
  
    xnames <- colnames(Xmean)[-1]

    boot.u <- df_estimates |> 
      select(dom, all_of(c(xnames, psid.est.var, y_var_dir)), term2Psid) |> 
      unique() |> 
      rename(all_of(c(psid.est=psid.est.var, y_var =y_var_dir))) |> 
      group_by(dom) |> 
      mutate(ed = rnorm(1, 0, sqrt(psid.est)),
             ud = rnorm(1, 0, sqrt(sigmau2)))

    boot.u <- boot.u |> 
      ungroup() |> 
      mutate(Yd.b = (Xmean %*% beta + ud)[,1],
             Ydw.b = Yd.b + ed) |> 
      as.data.frame()
    
    Yd.b.real <- boot.u$Yd.b
    
    # UA predictor
    form.mod.temp <- paste0("Ydw.b~", paste0(xnames, collapse="+")) |> formula()
    library(nlme)
    vf1Fixed <- varFixed(~term2Psid)
    fit.b <- try(lme(fixed = form.mod.temp, 
                    random = ~1|dom, weights = vf1Fixed, data = boot.u))
    
    if(class(fit.b) == "try-error"){
      fit.b <- try(lme(fixed = form.mod.temp, 
                       random = ~1|dom, weights = vf1Fixed, data = boot.u, method = "ML"))
    }
    sigmae2.fit.b <- fit.b$sigma^2
    boot.u$psid.est.FH.b <- sigmae2.fit.b * boot.u$term2Psid
    
    est.FHT.b <- eblupFH(form.mod.temp, vardir=psid.est.FH.b, data = boot.u)
    
    if(class(est.FHT.b) == "list"){
      if(est.FHT.b$fit$convergence == FALSE){
        est.FHT.b <- try(eblupFH(form.mod.temp, vardir=psid.est.FH.b, data = boot.u,
                                 method = "ML"))
      }
    }
    
    if(class(est.FHT.b) == "list"){
      if(est.FHT.b$fit$convergence == FALSE){
        est.FHT.b <- try(eblupFH(form.mod.temp, vardir=psid.est.FH.b, data = boot.u,
                                 method = "FH"))
      }
    }
    
    if(class(est.FHT.b) == "try-error"){
      est.FHT.b <- try(eblupFH(form.mod.temp, vardir=psid.est.FH.b, data = boot.u,
                               method = "ML"))
      if(class(est.FHT.b) == "try-error"){
        est.FHT.b <- try(eblupFH(form.mod.temp, vardir=psid.est.FH.b, data = boot.u,
                                 method = "FH"))
      }
    }
    
    res <- ((est.FHT.b$eblup[,1] - Yd.b.real)^2)

    res3 <- tibble(dom = df_estimates$dom, mse_FHT = res)
    
    return(res3)  
}
