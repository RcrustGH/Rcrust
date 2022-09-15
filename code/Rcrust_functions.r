correctMnzSat<-function (kd, c0, pm, cmins, min.props, melt.arg = list(), dont = character(0)) 
{
browser()
    c0 <- .sanitize(c0)
    m.pr <- .sanitize(min.props, normalize.TO = 1)
    lree.nm <- c("Th", "La", "Ce", "Pr", 
        "Nd", "Sm", "Gd")
    mw <- c(232.04, 138.91, 140.12, 140.91, 144.24, 150.3, 157.2, 
        30.974, 15.999)
    names(mw) <- c(lree.nm, "P", "O")
    if (any(toupper(melt.arg$mjrs) == melt.arg$mjrs)) {
        melt.arg$mjrs <- .TrueNames(melt.arg$mjrs)
    }
    W.R <- c(-520, -260, -116, 31, 177, 460, 750)
    names(W.R) <- lree.nm
    MW.REE <- mw[lree.nm]
    REE <- melt.arg$trc[lree.nm]
    M <- rep(0, length(lree.nm))
    Xmz <- 0.99
    milcats <- millications(melt.arg$mjrs)
    S <- REE[lree.nm]/MW.REE[lree.nm]
    ST <- sum(S)
    ee <- mzSaturationWithTh(cats = milcats, REE = REE, H2O = melt.arg$H2O, 
        Xmz = Xmz, Temp = melt.arg$TT - 273.15)
    AiT <- ee[1, "Sum.A.REE.sat"]
    if (ST <= AiT) {
        cm <- cmins
        cL <- melt.arg$trc
    }
    else {
        gamma <- exp(W.R/(melt.arg$TT))
        .compute <- function(MiT) {
            A <- gamma * AiT * S/(MiT + gamma * AiT)
            ee <- (sum(A) - AiT)^2
            return(ee)
        }
        m.guess <- 0
        out <- optim(m.guess, .compute, method = "Brent", 
            lower = 0, upper = 1)
        R2 <- out$value
        MiT <- out$par
        A <- gamma * AiT * S/(MiT + gamma * AiT)
        M <- S - A
        mf.mnz <- sum((MW.REE + mw["P"] + 4 * mw["O"]) * 
            M)/1e+06
        mnz.prop <- mf.mnz * pm/100
        lree.c <- A * MW.REE
        x <- M/sum(M)
        MW.mnz <- sum(x * MW.REE) + mw["P"] + 4 * mw["O"]
        cmnz <- x * MW.REE/MW.mnz * 1e+06
        kd["Mnz", lree.nm] <- cmnz/lree.c
        m.pr["Mnz"] <- mnz.prop
        m.pr <- m.pr/sum(m.pr)
        ee <- BatchPM(kd = kd, c0 = c0, pm = pm, min.props = m.pr)
        cL <- ee$cL
        cm <- ee$cmins
        cL[dont] <- melt.arg$trc[dont]
        all.but.mnz <- setdiff(rownames(cm), "Mnz")
        cm[all.but.mnz, dont] <- cmins[, dont]
        cL[lree.nm] <- lree.c
        cm["Mnz", lree.nm] <- cmnz
    }
    cS <- .sanitize(m.pr %*% cm)
    DD <- cS/cL
    dont <- c(dont, lree.nm)
    invisible(list(cL = .sanitize(cL), cmins = cm, kd = kd, DD = DD, 
        min.props = m.pr, dont = dont))
}