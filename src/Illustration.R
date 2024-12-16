###############################################################################
# # Illustration of U and UA predictor
# # Author: William Acero - wacero@ucm.es
###############################################################################
rm(list=ls())

# # Directories
inPath	<- "../input"

# # Libraries
library(pacman)
p_load(sae,tidyverse,furrr,progressr,ggh4x,srvyr,survey,nlme,fastDummies,
       data.table,ggpubr)
select <- dplyr::select

# Auxiliary functions
source("auxFunctions.R", encoding = "UTF-8")

# Number of replicates of PB MSE estimator of U and UA predictors
n.replicas.boot <- 500
seed <- 0305

# ICFES database
path <- file.path(inPath, "SB11_20182.RData")
load(path) # df

vars_pruebas <- c("PUNT_MATEMATICAS","PUNT_INGLES","PUNT_LECTURA_CRITICA","PUNT_C_NATURALES","PUNT_SOCIALES_CIUDADANAS")
df <- df |> select(ESTU_CONSECUTIVO, COLE_COD_DEPTO_UBICACION, COLE_DEPTO_UBICACION, all_of(vars_pruebas))

id_deptos <- df |> 
  group_by(COLE_COD_DEPTO_UBICACION, COLE_DEPTO_UBICACION) |> 
  count() |> 
  arrange(n) |> 
  ungroup() |> 
  mutate(dom = 1:n()) |> 
  select(COLE_COD_DEPTO_UBICACION, COLE_DEPTO_UBICACION, dom) |> 
  unique()

df <- df |> left_join(id_deptos |> select(-COLE_DEPTO_UBICACION), by = "COLE_COD_DEPTO_UBICACION")

df <- df |> select(-COLE_COD_DEPTO_UBICACION, -COLE_DEPTO_UBICACION) |> arrange(dom)

sizes <- df |> 
  group_by(dom) |> 
  count() |> 
  mutate(n_mue = round(n*0.015)) |> 
  select(dom, n_mue)

df_pop_totals <- df |> 
  group_by(dom) |> 
  summarise(Nd = n(), LC = sum(PUNT_LECTURA_CRITICA), CN = sum(PUNT_C_NATURALES) )

# Sample
set.seed(seed)
sam <- df %>% 
  group_by(dom) %>%
  group_split() %>%
  map2_dfr(sizes$n_mue, ~slice_sample(.x, n = .y)) %>%
  group_by(dom) %>%
  mutate(nd = n()) |> 
  left_join(df_pop_totals |> select(dom, Nd)) |> 
  mutate(dk = Nd/nd)


# Survey object
stsi.dsgn <- sam |> 
  as_survey_design(strata = dom, weights = dk, fpc = Nd)

sam.estimates <- stsi.dsgn |> 
  group_by(dom) |> 
  summarise(mate = survey_mean(PUNT_MATEMATICAS, vartype = c("var")))

sam.estimates <- sam.estimates |> 
  left_join(df_pop_totals |> 
    rowwise() |> 
    mutate(LC = LC/Nd, CN = CN/Nd) |> 
    select(-Nd), by = "dom") |> 
  data.frame()

################################################################################
# # Calibration
################################################################################

# Population totals
domains <- unique(sam$dom)
sam.estimates$mat_Cal <- NA
sam.estimates$mat_var_Cal <- NA
sam$dkCal <- NA

for(d in domains){
  stsi.dsgn <- sam |>
    filter(dom == d) |>
    as_survey_design(fpc = Nd, weights = dk)
  
  pop.totals <- df_pop_totals |> filter(dom  == d) |> 
    select(-dom, PUNT_LECTURA_CRITICA=LC, PUNT_C_NATURALES = CN) |> 
    as.vector() |> unlist()
  names(pop.totals)[1] <- '(Intercept)'
  
  greg.dsgn <- calibrate(design = stsi.dsgn,
                         formula = ~PUNT_LECTURA_CRITICA+PUNT_C_NATURALES,
                         calfun="linear",
                         population = pop.totals)

  sam <- sam |> mutate(dkCal = ifelse(dom == d, weights(greg.dsgn), dkCal))

  df_temp <- greg.dsgn |> summarise(mat = survey_mean(PUNT_MATEMATICAS, vartype = c("var")))
  sam.estimates <- sam.estimates |> mutate(mat_Cal = ifelse(dom == d, df_temp$mat, mat_Cal),
                                           mat_var_Cal = ifelse(dom == d, df_temp$mat, mat_var_Cal))
  cat("Finish ", d, "\n")
}

################################################################################
# # FHD predictor
################################################################################
est.FHD <- mseFH(mat_Cal ~ LC+CN, vardir = mat_var_Cal, data = sam.estimates)

################################################################################
# # U predictor and PB estimator
################################################################################
mod.est <- lmer(PUNT_MATEMATICAS ~ PUNT_LECTURA_CRITICA+PUNT_C_NATURALES + (1|dom), data = sam)
beta_est    <- fixef(mod.est)
sigmau2_est <- as.numeric(VarCorr(mod.est))
sigmae2_est <- summary(mod.est)$sigma^2

Xs <- sam |> ungroup() |> mutate(cons=1) |> select(cons,PUNT_LECTURA_CRITICA,PUNT_C_NATURALES) |> as.matrix()
Xd <- df_pop_totals |> mutate(cons=1) |> select(-dom) |> rowwise() |> mutate(LC = LC / Nd, CN = CN/Nd) |> ungroup() |> select(cons, LC, CN) |>  as.matrix()
nd <- sam |> select(dom, nd) |> unique() |> pull(nd)

# Diagonistic qqplot for residuals of U
ee <- residuals(mod.est, type='pearson')
sam$res <- ee

p.dm <- ggplot(sam, aes(x = res)) + 
  geom_histogram(aes(y = after_stat(density)),
                 colour = 1, fill = "white") +
  geom_density(alpha = .2, fill = "antiquewhite3") +
  theme_bw() + 
  xlab("Residuals")

df.residuals <- tibble(resid = ee)
p.qqRes <- df.residuals |> 
  ggplot(aes(sample = resid)) + 
    stat_qq() + 
    stat_qq_line() +
    theme_bw() +
    ylab("Sample quantiles") + 
    xlab("Theorical quantiles")

df.residuals$pred <- predict(mod.est)

p.fitRes <- df.residuals |> 
  ggplot(aes(x = pred, y = resid)) + 
  geom_point()+
  theme_bw() +
  ylab("Residuals") + 
  xlab("Fitted values")

u.pred <- ranef(mod.est)
df.ranef <- tibble(ranef = u.pred$dom[,1])
p.qqRanef <- df.ranef |> 
  ggplot(aes(sample = ranef)) + 
  stat_qq() + 
  stat_qq_line() +
  theme_bw() +
  ylab("Sample quantiles") + 
  xlab("Theorical quantiles")

# Diagnostic model plots
p3 <- ggarrange(ggarrange(p.dm, p.qqRes, ncol = 2, labels = c("A", "B")),
                ggarrange(p.fitRes, p.qqRanef, ncol = 2, labels = c("C", "D")),
                nrow = 2)

p3

# U predictor
est.U <- pseudo_EBLUP(ys=sam$PUNT_MATEMATICAS,
             ws=sam$dkCal,
             Xs=Xs,
             Xmean=Xd,
             areas=sam$dom,
             sigmae2est=sigmae2_est,
             sigmau2est=sigmau2_est,
             nd=nd)

# # MSE bootstrap estimator of U predictor
PopnSegments <- sam |> select(dom, Nd) |> unique() |> as.matrix()
  
results.PB <- NULL
set.seed(seed)
for(ii in 1:n.replicas.boot){
  results.PB[[ii]] <- mse_PB(sam = sam, 
         beta.est=beta_est,
         sigmau2.est=sigmau2_est,
         sigmae2.est=sigmae2_est,
         df_pop_means=Xd,
         Xs=Xs,
         PopnSegments=PopnSegments,
         nd= nd,
         formula.mod = 'PUNT_LECTURA_CRITICA+PUNT_C_NATURALES+(1|dom)',
         weights.calib = sam$dkCal)
  cat('Replicate ', ii, 'done \n')
}

mse.est.U.PB <- results.PB |> 
  bind_rows() |> 
  group_by(dom) |> 
  summarise_all(mean)

################################################################################
# # UA predictor and MSE bootstrap estimator 
################################################################################
df.UA <- sam.estimates

psid.est.UA <- sam |> 
  group_by(dom) |> 
  summarise(psid.est = sum(dkCal^2)/unique(Nd^2)) |> 
  pull(psid.est)

df.UA <- df.UA |> mutate(term2Psid = psid.est.UA)

vf1Fixed <- varFixed(~term2Psid)
fit1 <- try(lme(fixed = mat_Cal ~ LC + CN, 
  random = ~1|dom, weights = vf1Fixed, data = df.UA))

sigmae2.fit1 <- fit1$sigma^2

df.UA$psid.est.UA <- sigmae2.fit1 * df.UA$term2Psid

eblupUA <- eblupFH(mat_Cal ~ LC + CN, vardir=psid.est.UA, data = df.UA)

Xmean <- model.matrix(mat_Cal ~ LC + CN, df.UA)

results.PB <- NULL
set.seed(seed)
for(ii in 1:n.replicas.boot){
  results.PB[[ii]] <- mse_FHT(df_estimates = df.UA,
          beta = eblupUA$fit$estcoef[,1],
          sigmau2 = eblupUA$fit$refvar, 
          Xmean = Xmean,
          y_var_dir = "mat_var_Cal",
          psid.est.var = "psid.est.UA")
}

mse_PB_UA <- results.PB |> 
  bind_rows() |> 
  group_by(dom) |> 
  summarise_all(mean)

df.res <- tibble(dom = sam.estimates$dom,
                 Dir = sam.estimates$mat_Cal, 
                 Dir.var = sam.estimates$mat_var_Cal,
                 FHD = est.FHD$est$eblup[,1], 
                 mseFHD = est.FHD$mse,
                 U = est.U$mean_PEBLUP[,1],
                 mseU = mse.est.U.PB$mse.PB,
                 UA = eblupUA$eblup[,1], 
                 mseUA = mse_PB_UA$mse_FHT)

df.res <- df.res |> 
  left_join(id_deptos |> select(-COLE_COD_DEPTO_UBICACION), by = "dom") |> 
  relocate(COLE_DEPTO_UBICACION, .after = dom) |> 
  rename(Depto = COLE_DEPTO_UBICACION)

df.res
