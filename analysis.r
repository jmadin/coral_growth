# Growth 2 analysis

source("R/functions.R")

data <- read.csv("data/rspb20170053supp2_circularity.csv", as.is=TRUE)
data$x <- data$area_timestep / 10000
data$y <- data$area_nexttimestep / 10000
data$y_radius <- r_func(data$y)
data$x_radius <- r_func(data$x)
data$added_radius <- data$y_radius - data$x_radius
data <- data[data$added_radius < 0.2,]
data$x_log10 <- log10(data$x)
data$colony_id <- factor(data$colony_id)

gfs <- c("AH", "AC", "AI", "AR", "AD", "AS", "AL", "AN", "AM", "GR", "GP")
gfs_names <- c("Acropora hyacinthus", "Acropora cytherea", "Acropora intermedia", "Acropora robusta", "Acropora cf digitifera", "Acropora humilis", "Acropora spathulata", "Acropora nasuta", "Acropora millepora", "Goniastrea retiformis", "Goniastrea pectinata")

# GROWTH

png("output/figS1.png", width = 8.3, height = 6.5, units = 'in', res = 250)

  par(mfrow=c(3, 4), oma=c(2, 2, 0, 0), mar=c(3, 3, 2, 1))
  
  gfs_store <- data.frame()
  stat_growth0 <- data.frame()
  stat_growth1 <- data.frame()

  for (i in 1:length(gfs)) {
    if (i == 8) {
      plot(0,type='n',axes=FALSE,ann=FALSE)
    }
    dat1 <- data[data$species_code == gfs[i],]
    wholecolony <- sum(dat1$y == 0) / length(dat1$y)
    
    dat <- data[data$species_code == gfs[i] & data$y != 0,]
    x <- seq(log10(min(dat$x)), log10(max(dat$x)), 0.01)

    qr0 <- rq(added_radius ~ 1, 0.95, data=dat, ci=FALSE)
    qr1 <- rq(added_radius ~ log10(x), 0.95, data=dat, ci=FALSE)

    # qr2 <- lqmm(added_radius ~ log10(x), random = ~1, group=colony_id, data=dat, tau=0.999, nK = 11, type = "normal")
    # qr3 <- lqmm(added_radius ~ 1, random = ~1, group=colony_id, data=dat, tau=0.999, nK = 11, type = "normal")

    stat_growth0 <- rbind(stat_growth0, cbind(gfs_names[i], summary(qr0, se = "boot", bsmethod= "xy", R=10000)$coef))
    stat_growth1 <- rbind(stat_growth1, cbind(gfs_names[i], summary(qr1, se = "boot", bsmethod= "xy", R=10000)$coef))

    plot(added_radius ~ log10(x), dat, axes=FALSE, xlab="", ylab="Added radius by t+1 (m)", ylim=c(-0.3, 0.2), xlim=c(-3, 0.1), col="grey")
    mtext(make.italic(paste0(LETTERS[i], ". ", gfs_names[i])), adj=0)
    abline(h=0, lty=2)
    axis(1, at=c(-3, -2, -1, 0), labels=10^c(-3, -2, -1, 0))
    axis(2, las=2)
    lines(x, rep(qr0$coef, length(x)), lty=1, col="black")
    lines(x, coef(qr1)[1] + coef(qr1)[2] * x, lty=1, col="blue")
    lines(x, -r_func(10^x), lty=1)

    if (i == 1) {
      arrows(-1.9, -0.18, -1.55, -0.1, length=0.1)
      arrows(-2, 0.16, -2, 0.09, length=0.1)
      text(-2, -0.2, "Whole-colony death", cex=0.9)
      text(-2, 0.18, "Maximum radial growth", cex=0.9)
    }
    
    dat$potential_growth <- a_func(dat$x_radius + coef(qr0)[[1]])
    dat$pmort <- 1 - (dat$y / dat$potential_growth )
    dat <- dat[dat$pmort > 0.001,]
    dat$pmort_logit <- logit(dat$pmort)
    
    mod <- lm(pmort_logit ~ x_log10, dat)
    
    gfs_store <- rbind(gfs_store, data.frame(qr0=coef(qr0)[[1]], bint=coef(mod)[1], bslp=coef(mod)[2], bsig=sigma(mod), wholecolony))
  }
  mtext("Area at t (m^2)", 1, 1, outer=TRUE)

dev.off()

# Mortality

png("output/figS2.png", width = 8.3, height = 6.5, units = 'in', res = 250)
  
  par(mfrow=c(3, 4), oma=c(2, 2, 0, 0), mar=c(3, 3, 2, 1))
  
  stat_mort0 <- data.frame()
  stat_mort1 <- data.frame()
  
  for (i in 1:length(gfs)) {
    if (i == 8) {
      plot(0,type='n',axes=FALSE,ann=FALSE)
    }

    dat <- data[data$species_code == gfs[i] & data$y != 0,]
    x <- seq(log10(min(dat$x)), log10(max(dat$x)), 0.01)
    
    # Survival
    dat$potential_growth <- a_func(dat$x_radius + gfs_store$qr0[i])
    dat$pmort <- 1 - (dat$y / dat$potential_growth )
    dat <- dat[dat$pmort > 0.001,]
    dat$pmort_logit <- logit(dat$pmort)
    
    mod <- lm(pmort_logit ~ x_log10, dat)
    mod2 <- lm(pmort_logit ~ 1, dat)
    
    inv.logit(predict(mod, list(x_log10=log10(0.01))))
    stat_mort0 <- rbind(stat_mort0, cbind(species=gfs_names[i], m_pred=inv.logit(predict(mod, list(x_log10=log10(0.01))))))
    stat_mort1 <- rbind(stat_mort1, cbind(gfs_names[i], summary(mod)$coef))
    
    plot(pmort ~ x_log10, dat, xlab="", ylab="Proportion lost by t+1", axes=FALSE, ylim=c(0, 1), xlim=c(-3, 0), col="grey")
    axis(1, at=c(-3, -2, -1, 0), labels=10^c(-3, -2, -1, 0))
    axis(2, las=2)
    pp <- inv.logit(predict(mod, data.frame(x_log10=x), interval="confidence"))
    lines(x, pp[,1])
    x0 <- 10^x
    r0 <- r_func(x0)
    r1 <- r0 + gfs_store$qr0[i]
    x1 <- a_func(r1)

    lines(x, (x1 - x0) / x1, lty=2)
    if (i == 1) {
      arrows(-1.7, 0.95, -2, 0.82, length=0.1)
      text(-1.6, 1, "Colony stasis", cex=0.9)
    }
    mtext(make.italic(paste0(LETTERS[i], ". ", gfs_names[i])), adj=0)
  }
  mtext("Area at t (m^2)", 1, 1, outer=TRUE)
  
dev.off()

# Trade-off

trade <- data.frame(species=stat_growth0$V1, growth=(as.numeric(stat_growth0$Value)), mort=as.numeric(stat_mort0$m_pred), spp=LETTERS[1:11], spp2=c("Ahya", "Acyt", "Aint", "Arob", "Adig", "Ahum", "Aspa", "Anas", "Amil", "Gret", "Gpec"), whole_colony=gfs_store$wholecolony)

png("output/fig1.png", width=6, height=5.5, units = 'in', res = 250)

  par(mar=c(6, 6, 2, 1))

  plot(mort ~ growth, trade, xlab="", ylab="", axes=FALSE, xlim=c(0, 0.11), ylim=c(0.3, 0.8), pch=NA, cex=0.75)
  text(trade$growth, trade$mort, make.italic(trade$spp2), cex=0.7)
  axis(1, at=seq(0, 0.1, 0.02), las=2)
  axis(2, las=2)
  mtext("Growth", 1, 3, cex=1.2)
  mtext("(Added radius by t+1, m)", 1, 4, cex=1)
  mtext("Partial mortality", 2, 4, cex=1.2)
  mtext("(Proportion lost by t+1 at 0.01 m^2)", 2, 3, cex=1)
  
  rasterImage(readPNG("data/silh_massive.png"), 0.002, 0.325, 0.017, 0.37) 
  rasterImage(readPNG("data/silh_digitate.png"), 0.012, 0.38, 0.028, 0.44) 
  rasterImage(readPNG("data/silh_tabular.png"), 0.071, 0.73, 0.091, 0.78) 
  rasterImage(readPNG("data/silh_corymbose.png"), 0.025, 0.48, 0.038, 0.54) 
  rasterImage(readPNG("data/silh_branching.png"), 0.048, 0.62, 0.058, 0.72) 
 
dev.off()

# RELATIVE

png("output/fig2.png", width = 8.3, height = 6.5, units = 'in', res = 250)

  par(mfrow=c(3, 4), oma=c(2, 2, 0, 0), mar=c(3, 3, 2, 1))

  xi <- yi <- seq(-3, 0.1, 0.02)
  for (i in 1:length(gfs)) {

    # dat <- data2[data2$group.5 == gfs[i],]
    dat <- data[data$species_code == gfs[i],]

    store <- matrix(NA, length(xi), length(yi))
    for (j in 1:length(xi)) {
      store[j,] <- g_func(xi[j], yi, qr0=gfs_store$qr0[i], p1=gfs_store$bint[i], p2=gfs_store$bslp[i], p3=gfs_store$bsig[i])
    }
    image(store, x = xi, y = yi, col=rev(heat.colors(50)), axes=FALSE, xlab="Area at t (m^2)", ylab="", asp=1)
    points(log10(y) ~ log10(x), dat, col="grey")
    mtext(make.italic(paste0(LETTERS[i], ". ", gfs_names[i])), adj=0)
    abline(0, 1, lty=2)
    axis(1, at=c(-3, -2, -1, 0), labels=10^c(-3, -2, -1, 0))
    axis(2, at=c(-3, -2, -1, 0), labels=10^c(-3, -2, -1, 0), las=2)

    lines(xi, log10(a_func(r_func(10^xi) + gfs_store$qr0[i])))

    if (gfs[i] == "AL") {
      plot(0, 0, type = "n", axes = FALSE, xlab="", ylab="")
    }
  }
  mtext("Area at t (m^2)", 1, 0, outer=TRUE)
  mtext("Area at t+1 (m^2)", 2, 0, outer=TRUE)

dev.off()


# SAVE STATS

stat_growth0[] <- lapply(stat_growth0, as.character)
stat_growth0[,2:5] <- lapply(stat_growth0[,2:5], as.numeric)
stat_growth0[,2:5] <- lapply(stat_growth0[,2:5], round, digits=3)
stat_growth0$param <- c("Intercept")

write.csv(stat_growth0[,c(1, 6, 2, 3, 4, 5)], "output/stat_growth0.csv", row.names=FALSE, quote=FALSE)

stat_growth1[] <- lapply(stat_growth1, as.character)
stat_growth1[,2:5] <- lapply(stat_growth1[,2:5], as.numeric)
stat_growth1[,2:5] <- lapply(stat_growth1[,2:5], round, digits=3)
stat_growth1$param <- c("Intercept", "Slope")

write.csv(stat_growth1[,c(1, 6, 2, 3, 4, 5)], "output/stat_growth1.csv", row.names=FALSE, quote=FALSE)

stat_mort1[] <- lapply(stat_mort1, as.character)
stat_mort1[,2:5] <- lapply(stat_mort1[,2:5], as.numeric)
stat_mort1[,2:5] <- lapply(stat_mort1[,2:5], round, digits=3)
stat_mort1$param <- c("Intercept", "Slope")

write.csv(stat_mort1[,c(1, 6, 2, 3, 4, 5)], "output/stat_mort.csv", row.names=FALSE, quote=FALSE)

# Circularity

png("output/figS3.png", width = 9, height = 6.5, units = 'in', res = 250)
  par(mfrow=c(3, 4), oma=c(2, 2, 0, 0), mar=c(3, 5, 2, 1))

  for (i in 1:length(gfs)) {

    dat <- data[data$species_code == gfs[i],]

    qr0 <- rq(added_radius ~ 1, 0.95, data=dat, ci=FALSE)
    plot(dat$circularity, residuals(qr0), axes=FALSE, xlab="", ylab="", ylim=c(-0.3, 0.1), xlim=c(0.5, 1), xpd=TRUE, col="grey")
    axis(1)
    axis(2, las=2)

    ct <- cor.test(dat$circularity, residuals(qr0), method = "spearman")
    text(0.6, -0.2, paste0("rho=", round(ct$estimate, 3)))

    mtext(make.italic(paste0(LETTERS[i], ". ", gfs_names[i])), adj=0)
    abline(h=0, lty=2)

    if (gfs[i] == "AL") {
      plot(0, 0, type = "n", axes = FALSE, xlab="", ylab="")
    }
  }
  mtext("Circularity", 1, 0, outer=TRUE)
  mtext("Residuals, 95% quantile for added radius", 2, 0, outer=TRUE)

dev.off()
