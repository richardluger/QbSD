rm(list = ls())  ## clear the workspace

load("../../Data/Data.Rdata")

library(MCS)

bnam <- c("SP500", "DJIA", "NASDAQ", "STOXX50", "FTSE", "DAX", "CAC40", "TSX")

nb <- length(bnam)
bdate <- matrix(
  c(
    2002, 250 * 10 / 12 + 18,
    2002, 250 * 10 / 12 + 18,
    2002, 250 * 10 / 12 + 18,
    2002, 250 * 10 / 12 + 18,
    2002, 250 * 10 / 12 + 18,
    2002, 250 * 10 / 12 + 18,
    2002, 250 * 10 / 12 + 18,
    2002, 250 * 10 / 12 + 18
  ),
  nb, 2, byrow = TRUE
)

datt <- list()
for (n in 1:nb) {
  datt[[n]] <- Data[[n]]
}

rw  <- c(250, 1250, 2500)
ALP <- c(0.01, 0.025, 0.05)

RLOC <- c("zero", "QAR")
CAV  <- c("SAV", "AS")
RFO  <- c("VaRf", "ESf")
LOC  <- c("zero", "AR")
ESD  <- c("basic", "AR")
gmo  <- c("GARCH-n", "GARCH-t", "GARCH-skewt", "GJR-n", "GJR-t", "GJR-skewt", "EGARCH-n")

####### load VaR and ES forecasts
GARCHFF <- list()
load("../Full_forecasting_output/GARCH.SP500.Rdata");   GARCHFF[[1]] <- GARCHF[[1]]
load("../Full_forecasting_output/GARCH.DJIA.Rdata");    GARCHFF[[2]] <- GARCHF[[2]]
load("../Full_forecasting_output/GARCH.NASDAQ.Rdata");  GARCHFF[[3]] <- GARCHF[[3]]
load("../Full_forecasting_output/GARCH.STOXX50.Rdata"); GARCHFF[[4]] <- GARCHF[[4]]
load("../Full_forecasting_output/GARCH.FTSE.Rdata");    GARCHFF[[5]] <- GARCHF[[5]]
load("../Full_forecasting_output/GARCH.DAX.Rdata");     GARCHFF[[6]] <- GARCHF[[6]]
load("../Full_forecasting_output/GARCH.CAC40.Rdata");   GARCHFF[[7]] <- GARCHF[[7]]
load("../Full_forecasting_output/GARCH.TSX.Rdata");     GARCHFF[[8]] <- GARCHF[[8]]

PattonFF <- list()
load("../Full_forecasting_output/Patton.SP500.Rdata");  PattonFF[[1]] <- PattonF[[1]]
load("../Full_forecasting_output/Patton.DJIA.Rdata");   PattonFF[[2]] <- PattonF[[2]]
load("../Full_forecasting_output/Patton.NASDAQ.Rdata"); PattonFF[[3]] <- PattonF[[3]]
load("../Full_forecasting_output/Patton.STOXX50.Rdata");PattonFF[[4]] <- PattonF[[4]]
load("../Full_forecasting_output/Patton.FTSE.Rdata");   PattonFF[[5]] <- PattonF[[5]]
load("../Full_forecasting_output/Patton.DAX.Rdata");    PattonFF[[6]] <- PattonF[[6]]
load("../Full_forecasting_output/Patton.CAC40.Rdata");  PattonFF[[7]] <- PattonF[[7]]
load("../Full_forecasting_output/Patton.TSX.Rdata");    PattonFF[[8]] <- PattonF[[8]]

RobustFF <- list()
load("../Full_forecasting_output/QbSD.SP500.Rdata");  RobustFF[[1]] <- RobustF[[1]]
load("../Full_forecasting_output/QbSD.DJIA.Rdata");   RobustFF[[2]] <- RobustF[[2]]
load("../Full_forecasting_output/QbSD.NASDAQ.Rdata"); RobustFF[[3]] <- RobustF[[3]]
load("../Full_forecasting_output/QbSD.STOXX50.Rdata");RobustFF[[4]] <- RobustF[[4]]
load("../Full_forecasting_output/QbSD.FTSE.Rdata");   RobustFF[[5]] <- RobustF[[5]]
load("../Full_forecasting_output/QbSD.DAX.Rdata");    RobustFF[[6]] <- RobustF[[6]]
load("../Full_forecasting_output/QbSD.CAC40.Rdata");  RobustFF[[7]] <- RobustF[[7]]
load("../Full_forecasting_output/QbSD.TSX.Rdata");    RobustFF[[8]] <- RobustF[[8]]

TaylorFF <- list()
load("../Full_forecasting_output/Taylor.SP500.Rdata");  TaylorFF[[1]] <- TaylorF[[1]]
load("../Full_forecasting_output/Taylor.DJIA.Rdata");   TaylorFF[[2]] <- TaylorF[[2]]
load("../Full_forecasting_output/Taylor.NASDAQ.Rdata"); TaylorFF[[3]] <- TaylorF[[3]]
load("../Full_forecasting_output/Taylor.STOXX50.Rdata");TaylorFF[[4]] <- TaylorF[[4]]
load("../Full_forecasting_output/Taylor.FTSE.Rdata");   TaylorFF[[5]] <- TaylorF[[5]]
load("../Full_forecasting_output/Taylor.DAX.Rdata");    TaylorFF[[6]] <- TaylorF[[6]]
load("../Full_forecasting_output/Taylor.CAC40.Rdata");  TaylorFF[[7]] <- TaylorF[[7]]
load("../Full_forecasting_output/Taylor.TSX.Rdata");    TaylorFF[[8]] <- TaylorF[[8]]

GARCHF  <- GARCHFF
PattonF <- PattonFF
RobustF <- RobustFF
TaylorF <- TaylorFF

GARCHnam <- c()
for (s in 1:length(LOC)) {
  for (g in 1:length(gmo)) {
    GARCHnam <- c(GARCHnam, paste(gmo[g], LOC[s]))
  }
}
GARCHmat <- matrix(1:(length(gmo) * length(LOC)), length(gmo), length(LOC))

Robustnam <- c()
for (s in 1:length(RLOC)) {
  for (v in 1:length(CAV)) {
    Robustnam <- c(Robustnam, paste("Robust", RLOC[s], CAV[v]))
  }
}
Robustmat <- array(1:4, c(length(CAV), length(RLOC)))

Taylornam <- c()
for (s in 1:length(LOC)) {
  for (v in 1:length(CAV)) {
    for (d in 1:length(ESD)) {
      Taylornam <- c(Taylornam, paste("Taylor", LOC[s], CAV[v], ESD[d]))
    }
  }
}
Taylormat <- array(0, c(length(ESD), length(CAV), length(LOC)))
Taylormat[,,1] <- matrix(1:4,  length(ESD), length(CAV))
Taylormat[,,2] <- matrix(5:8,  length(ESD), length(CAV))

QSS_GARCH  <- list()
QSS_Patton <- list()
QSS_Robust <- list()
QSS_Taylor <- list()

ALLS_GARCH  <- list()
ALLS_Patton <- list()
ALLS_Robust <- list()
ALLS_Taylor <- list()

for (n in 1:nb) {

  dat <- as.numeric(datt[[n]])
  TT  <- length(dat)

  QSS_GARCH[[n]]  <- list()
  QSS_Patton[[n]] <- list()
  QSS_Robust[[n]] <- list()
  QSS_Taylor[[n]] <- list()

  ALLS_GARCH[[n]]  <- list()
  ALLS_Patton[[n]] <- list()
  ALLS_Robust[[n]] <- list()
  ALLS_Taylor[[n]] <- list()

  for (w in 3) {  # <— run only w = 3  (= 2500)

    win <- rw[w]
    fp  <- TT - win
    yf  <- dat[(win + 1):TT]

    QSS_GARCH[[n]][[w]] <- array(
      0, c(fp, length(gmo) * length(LOC), length(ALP)),
      dimnames = list(c(), GARCHnam, ALP)
    )
    QSS_Patton[[n]][[w]] <- array(0, c(fp, length(ALP)), dimnames = list(c(), ALP))
    QSS_Robust[[n]][[w]] <- list()
    QSS_Taylor[[n]][[w]] <- array(
      0, c(fp, length(LOC) * length(ESD) * length(CAV), length(ALP)),
      dimnames = list(c(), Taylornam, ALP)
    )

    ALLS_GARCH[[n]][[w]] <- array(
      0, c(fp, length(gmo) * length(LOC), length(ALP)),
      dimnames = list(c(), GARCHnam, ALP)
    )
    ALLS_Patton[[n]][[w]] <- array(0, c(fp, length(ALP)), dimnames = list(c(), ALP))
    ALLS_Robust[[n]][[w]] <- list()
    ALLS_Taylor[[n]][[w]] <- array(
      0, c(fp, length(LOC) * length(ESD) * length(CAV), length(ALP)),
      dimnames = list(c(), Taylornam, ALP)
    )

    for (a in 1:length(ALP)) {

      print(c(bnam[n], rw[w], ALP[a]))
      garch  <- GARCHF[[n]][[w]][,,,,a]
      patton <- PattonF[[n]][[w]][,,a]
      robust <- RobustF[[n]][[w]][[a]]
      taylor <- TaylorF[[n]][[w]][,,,,,a]

      ## GARCH models
      for (s in 1:length(LOC)) {
        for (g in 1:length(gmo)) {
          VAR    <- garch[, g, 1, s]
          ES     <- garch[, g, 2, s]
          qscore <- (yf - VAR) * (ALP[a] - (yf <= VAR) * 1)
          QSS_GARCH[[n]][[w]][, GARCHmat[g, s], a]  <- qscore
          ALLS_GARCH[[n]][[w]][, GARCHmat[g, s], a] <- -log(abs((ALP[a] - 1) / ES)) - qscore / (ALP[a] * ES)
        }
      }

      ### Patton et al (2019) PZC
      VAR    <- patton[, 1]
      ES     <- patton[, 2]
      qscore <- (yf - VAR) * (ALP[a] - (yf <= VAR) * 1)
      QSS_Patton[[n]][[w]][, a]  <- qscore
      ALLS_Patton[[n]][[w]][, a] <- -log(abs((ALP[a] - 1) / ES)) - qscore / (ALP[a] * ES)

      ######## Robust
      pp <- seq(from = 0.05, to = 0.25, by = 0.05)
      QSS_Robust[[n]][[w]][[a]] <- array(
        0, c(fp, length(CAV) * length(RLOC)),
        dimnames = list(c(), Robustnam)
      )
      ALLS_Robust[[n]][[w]][[a]] <- array(
        0, c(fp, length(CAV) * length(RLOC)),
        dimnames = list(c(), Robustnam)
      )

      for (s in 1:length(RLOC)) {
        for (v in 1:length(CAV)) {
          VAR    <- robust[, 1, v, s]
          ES     <- robust[, 2, v, s]
          dr     <- Robustmat[v, s]
          qscore <- (yf - VAR) * (ALP[a] - (yf <= VAR) * 1)
          QSS_Robust[[n]][[w]][[a]][, dr]  <- qscore
          ALLS_Robust[[n]][[w]][[a]][, dr] <- -log(abs((ALP[a] - 1) / ES)) - qscore / (ALP[a] * ES)
        }
      }

      ############## Taylor
      for (s in 1:length(LOC)) {
        for (v in 1:length(CAV)) {
          for (d in 1:length(ESD)) {
            VAR    <- taylor[, 1, d, v, s]
            ES     <- taylor[, 2, d, v, s]
            dn     <- Taylormat[d, v, s]
            qscore <- (yf - VAR) * (ALP[a] - (yf <= VAR) * 1)
            QSS_Taylor[[n]][[w]][, dn, a]  <- qscore
            ALLS_Taylor[[n]][[w]][, dn, a] <- -log(abs((ALP[a] - 1) / ES)) - qscore / (ALP[a] * ES)
          }
        }
      }
    }
  }
}

########################################################################################
########################################################################################
########################################################################################
########################### MCS

MCSS_QS <- list()
MCSS_AL <- list()

for (n in 1:nb) {
  MCSS_QS[[n]] <- list()
  MCSS_AL[[n]] <- list()
  for (w in 3) {  # <— run only w = 3
    MCSS_QS[[n]][[w]] <- list()
    MCSS_AL[[n]][[w]] <- list()

    for (a in 1:length(ALP)) {
      print(c(bnam[n], rw[w], ALP[a]))

      QSS  <- cbind(
        QSS_GARCH[[n]][[w]][,,a],
        QSS_Patton[[n]][[w]][,a],
        QSS_Taylor[[n]][[w]][,,a],
        QSS_Robust[[n]][[w]][[a]]
      )
      dimnames(QSS)[[2]][dim(QSS_GARCH[[n]][[w]][,,a])[2] + 1] <- "Patton"

      ALLS <- cbind(
        ALLS_GARCH[[n]][[w]][,,a],
        ALLS_Patton[[n]][[w]][,a],
        ALLS_Taylor[[n]][[w]][,,a],
        ALLS_Robust[[n]][[w]][[a]]
      )
      dimnames(ALLS)[[2]][dim(ALLS_GARCH[[n]][[w]][,,a])[2] + 1] <- "Patton"

      k <- dim(QSS)[2]
      MCSS_QS[[n]][[w]][[a]] <- MCSprocedure(Loss = QSS,  alpha = 0.1, statistic = "TR", B = 1000)
      MCSS_AL[[n]][[w]][[a]] <- MCSprocedure(Loss = ALLS, alpha = 0.1, statistic = "TR", B = 1000)
    }
  }
}

############################################

RLOC <- c("zero", "QAR")
CAV  <- c("SAV", "AS")
RFO  <- c("VaRf", "ESf")
LOC  <- c("zero", "AR")
ESD  <- c("basic", "AR")
gmo  <- c("GARCH-n", "GARCH-t", "GARCH-skewt", "GJR-n", "GJR-t", "GJR-skewt", "EGARCH-n")

GARCHnam <- c()
for (s in 1:length(LOC)) {
  for (g in 1:length(gmo)) {
    GARCHnam <- c(GARCHnam, paste(gmo[g], LOC[s]))
  }
}

Robustnam <- c()
for (s in 1:length(RLOC)) {
  for (v in 1:length(CAV)) {
    Robustnam <- c(Robustnam, paste("Robust", RLOC[s], CAV[v]))
  }
}

Taylornam <- c()
for (s in 1:length(LOC)) {
  for (v in 1:length(CAV)) {
    for (d in 1:length(ESD)) {
      Taylornam <- c(Taylornam, paste("Taylor", LOC[s], CAV[v], ESD[d]))
    }
  }
}

Mnam <- c(GARCHnam, "Patton", Taylornam, Robustnam)
nm   <- length(Mnam)

Mnam1 <- c(
  "GARCH-Normal", "GARCH-t", "GARCH-Skew t", "GJR-Normal", "GJR-t", "GJR-Skew t",
  "EGARCH-Normal", "AR-GARCH-Normal", "AR-GARCH-t", "AR-GARCH-Skew t",
  "AR-GJR-Normal", "AR-GJR-t", "AR-GJR-Skew t", "AR-EGARCH-Normal",
  "GAS",
  "ALMult-SAV", "ALAR-SAV", "ALMulti-AS", "ALAR-AS", "AR-ALMult-SAV",
  "AR-ALAR-SAV", "AR-ALMulti-AS", "AR-ALAR-AS",
  "QbSD-gSAV", "QbSD-gAS", "QAR-QbSD-gSAV", "QAR-QbSD-gAS"
)

bnam <- c("SP500", "DJIA", "NASDAQ", "STOXX50", "FTSE", "DAX", "CAC40", "TSX")

nb  <- length(bnam)
rw  <- c(250, 1250, 2500)
ALP <- c(0.01, 0.025, 0.05)

########### Produce Tables 8,9,10

MCSR_AL <- array(
  0, c(length(Mnam1), nb, length(ALP), length(rw)),
  dimnames = list(Mnam1, bnam, ALP, rw)
)

for (n in 1:nb) {
  for (w in 3) {  # <— run only w = 3
    for (a in 1:length(ALP)) {
      mc   <- MCSS_AL[[n]][[w]][[a]]@show
      snam <- dimnames(mc)[[1]]
      for (i in 1:length(snam)) {
        loc <- which(snam[i] == Mnam)
        MCSR_AL[loc, n, a, w] <- mc[i, 4]
      }
    }
  }
}

MCS_AL_Rank <- list()
w <- 3
for (a in 1:length(ALP)) {
  mm <- MCSR_AL[,,a,w]
  MA <- matrix(
    0, nrow = dim(mm)[1], ncol = dim(mm)[2] + 3,
    dimnames = list(Mnam1, c(bnam, "#", "Avg.rank", "Final.rank"))
  )

  MA[, 1:nb] <- mm
  for (i in 1:nm) {
    MA[i, nb + 1] <- sum(mm[i, ] != 0)
    ind           <- (mm[i, ] == 0) * 1
    MA[i, nb + 2] <- round((sum(ind) * nm + sum((1 - ind) * mm[i, ])) / nb, 1)
  }

  od <- order(MA[, nb + 2])
  MCS_AL_Rank[[a]] <- MA[od, ]
  MCS_AL_Rank[[a]][, nb + 3] <- 1:nm
}

Table8  <- MCS_AL_Rank[[1]]
Table9  <- MCS_AL_Rank[[2]]
Table10 <- MCS_AL_Rank[[3]]
for (i in 1:nm) {
  for (j in 1:nb) {
    if (Table8[i, j]  == 0) Table8[i, j]  <- NA  ## "NA" represents empty cells
    if (Table9[i, j]  == 0) Table9[i, j]  <- NA
    if (Table10[i, j] == 0) Table10[i, j] <- NA
  }
}

print("Table 8");  Table8
print("Table 9");  Table9
print("Table 10"); Table10

########### Produce Tables 5,6,7

MCSR_QS <- array(
  0, c(length(Mnam1), nb, length(ALP), length(rw)),
  dimnames = list(Mnam1, bnam, ALP, rw)
)

for (n in 1:nb) {
  for (w in 3) {  # <— run only w = 3
    for (a in 1:length(ALP)) {
      mc   <- MCSS_QS[[n]][[w]][[a]]@show
      snam <- dimnames(mc)[[1]]
      for (i in 1:length(snam)) {
        loc <- which(snam[i] == Mnam)
        MCSR_QS[loc, n, a, w] <- mc[i, 4]
      }
    }
  }
}

MCS_QS_Rank <- list()
nnma <- ls()
w <- 3
for (a in 1:length(ALP)) {
  mm <- MCSR_QS[,,a,w]
  MA <- matrix(
    0, nrow = dim(mm)[1], ncol = dim(mm)[2] + 3,
    dimnames = list(Mnam1, c(bnam, "#", "Avg.rank", "Final.rank"))
  )

  MA[, 1:nb] <- mm
  for (i in 1:nm) {
    MA[i, nb + 1] <- sum(mm[i, ] != 0)
    ind           <- (mm[i, ] == 0) * 1
    MA[i, nb + 2] <- round((sum(ind) * nm + sum((1 - ind) * mm[i, ])) / nb, 1)
  }

  od <- order(MA[, nb + 2])
  MCS_QS_Rank[[a]] <- MA[od, ]
  MCS_QS_Rank[[a]][, nb + 3] <- 1:nm
}

Table5 <- MCS_QS_Rank[[1]]
Table6 <- MCS_QS_Rank[[2]]
Table7 <- MCS_QS_Rank[[3]]
for (i in 1:nm) {
  for (j in 1:nb) {
    if (Table5[i, j] == 0) Table5[i, j] <- NA  ## "NA" represents empty cells
    if (Table6[i, j] == 0) Table6[i, j] <- NA
    if (Table7[i, j] == 0) Table7[i, j] <- NA
  }
}

print("Table 5"); Table5
print("Table 6"); Table6
print("Table 7"); Table7
