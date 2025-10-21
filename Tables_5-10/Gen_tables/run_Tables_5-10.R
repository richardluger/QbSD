rm(list = ls())  ## clear the workspace

p <- "MCS"
if (!requireNamespace(p, quietly = TRUE)) {
  install.packages(p)
}
library(p, character.only = TRUE)

set.seed(123456)

load("../../Data/Data.Rdata")

bnam <- names(Data)
nb <- length(bnam) 

datt <- list()
for (n in 1:nb) {
  datt[[n]] <- Data[[n]]
}

rw  <- c(250, 1250, 2500)
ALP <- c(0.01, 0.025, 0.05)

RLOC <- c("zero", "QAR")
CAV  <- c("SAV", "AS")
LOC  <- c("zero", "AR")
ESD <- c("mult", "ar")
gmo  <- c("GARCH-n", "GARCH-t", "GARCH-skewt", "GJR-n", "GJR-t", "GJR-skewt", "EGARCH-n")

GARCHFF <- vector("list", nb)
GASFF   <- vector("list", nb)
QbSDFF  <- vector("list", nb)
ALFF    <- vector("list", nb)

####### load VaR and ES forecasts
## GARCH
load("../Full_forecasting_output/GARCH.SP500.Rdata");   GARCHFF[[1]] <- GARCH_oos[[1]]
load("../Full_forecasting_output/GARCH.DJIA.Rdata");    GARCHFF[[2]] <- GARCH_oos[[2]]
load("../Full_forecasting_output/GARCH.NASDAQ.Rdata");  GARCHFF[[3]] <- GARCH_oos[[3]]
load("../Full_forecasting_output/GARCH.STOXX50.Rdata"); GARCHFF[[4]] <- GARCH_oos[[4]]
load("../Full_forecasting_output/GARCH.FTSE.Rdata");    GARCHFF[[5]] <- GARCH_oos[[5]]
load("../Full_forecasting_output/GARCH.DAX.Rdata");     GARCHFF[[6]] <- GARCH_oos[[6]]
load("../Full_forecasting_output/GARCH.CAC40.Rdata");   GARCHFF[[7]] <- GARCH_oos[[7]]
load("../Full_forecasting_output/GARCH.TSX.Rdata");     GARCHFF[[8]] <- GARCH_oos[[8]]

## GAS
load("../Full_forecasting_output/GAS.SP500.Rdata");     GASFF[[1]] <- GAS_oos[[1]]
load("../Full_forecasting_output/GAS.DJIA.Rdata");      GASFF[[2]] <- GAS_oos[[2]]
load("../Full_forecasting_output/GAS.NASDAQ.Rdata");    GASFF[[3]] <- GAS_oos[[3]]
load("../Full_forecasting_output/GAS.STOXX50.Rdata");   GASFF[[4]] <- GAS_oos[[4]]
load("../Full_forecasting_output/GAS.FTSE.Rdata");      GASFF[[5]] <- GAS_oos[[5]]
load("../Full_forecasting_output/GAS.DAX.Rdata");       GASFF[[6]] <- GAS_oos[[6]]
load("../Full_forecasting_output/GAS.CAC40.Rdata");     GASFF[[7]] <- GAS_oos[[7]]
load("../Full_forecasting_output/GAS.TSX.Rdata");       GASFF[[8]] <- GAS_oos[[8]]

## QbSD 
load("../Full_forecasting_output/QbSD.SP500.Rdata");    QbSDFF[[1]] <- QbSD_oos[[1]]
load("../Full_forecasting_output/QbSD.DJIA.Rdata");     QbSDFF[[2]] <- QbSD_oos[[2]]
load("../Full_forecasting_output/QbSD.NASDAQ.Rdata");   QbSDFF[[3]] <- QbSD_oos[[3]]
load("../Full_forecasting_output/QbSD.STOXX50.Rdata");  QbSDFF[[4]] <- QbSD_oos[[4]]
load("../Full_forecasting_output/QbSD.FTSE.Rdata");     QbSDFF[[5]] <- QbSD_oos[[5]]
load("../Full_forecasting_output/QbSD.DAX.Rdata");      QbSDFF[[6]] <- QbSD_oos[[6]]
load("../Full_forecasting_output/QbSD.CAC40.Rdata");    QbSDFF[[7]] <- QbSD_oos[[7]]
load("../Full_forecasting_output/QbSD.TSX.Rdata");      QbSDFF[[8]] <- QbSD_oos[[8]]

## AL 
load("../Full_forecasting_output/AL.SP500.Rdata");      ALFF[[1]] <- AL_oos[[1]]
load("../Full_forecasting_output/AL.DJIA.Rdata");       ALFF[[2]] <- AL_oos[[2]]
load("../Full_forecasting_output/AL.NASDAQ.Rdata");     ALFF[[3]] <- AL_oos[[3]]
load("../Full_forecasting_output/AL.STOXX50.Rdata");    ALFF[[4]] <- AL_oos[[4]]
load("../Full_forecasting_output/AL.FTSE.Rdata");       ALFF[[5]] <- AL_oos[[5]]
load("../Full_forecasting_output/AL.DAX.Rdata");        ALFF[[6]] <- AL_oos[[6]]
load("../Full_forecasting_output/AL.CAC40.Rdata");      ALFF[[7]] <- AL_oos[[7]]
load("../Full_forecasting_output/AL.TSX.Rdata");        ALFF[[8]] <- AL_oos[[8]]


GARCHF  <- GARCHFF
GASF    <- GASFF
QbSDF <- QbSDFF
ALF     <- ALFF

GARCHnam <- NULL
for (s in 1:length(LOC)) {
  for (g in 1:length(gmo)) {
    GARCHnam <- c(GARCHnam, paste(gmo[g], LOC[s]))
  }
}
GARCHmat <- matrix(1:(length(gmo) * length(LOC)), length(gmo), length(LOC))

QbSDnam <- NULL
for (s in 1:length(RLOC)) {
  for (v in 1:length(CAV)) {
    QbSDnam <- c(QbSDnam, paste("QbSD", RLOC[s], CAV[v]))
  }
}
QbSDmat <- array(1:4, c(length(CAV), length(RLOC)))

ALnam <- NULL
for (s in 1:length(LOC)) {
  for (v in 1:length(CAV)) {
    for (d in 1:length(ESD)) {
      ALnam <- c(ALnam, paste("AL", LOC[s], CAV[v], ESD[d]))
    }
  }
}
ALmat <- array(0, c(length(ESD), length(CAV), length(LOC)))
ALmat[,,1] <- matrix(1:4,  length(ESD), length(CAV))
ALmat[,,2] <- matrix(5:8,  length(ESD), length(CAV))

QSS_GARCH  <- list()
QSS_GAS <- list()
QSS_QbSD <- list()
QSS_AL <- list()

ALLS_GARCH  <- list()
ALLS_GAS <- list()
ALLS_QbSD <- list()
ALLS_AL <- list()



# Build model name lists ONCE, before MCSS
GASnam  <- paste("GAS", LOC)  
Mnam    <- c(
  (function(){
    out <- character()
    for (s in LOC) for (g in gmo) out <- c(out, paste(g, s))
    out
  })(),
  GASnam,
  (function(){
    out <- character()
    for (s in LOC) for (v in CAV) for (d in ESD) out <- c(out, paste("AL", s, v, d))
    out
  })(),
  (function(){
    out <- character()
    for (s in RLOC) for (v in CAV) out <- c(out, paste("QbSD", s, v))
    out
  })()
)
nm   <- length(Mnam)
Mnam1 <- Mnam  


for (n in 1:nb) {

  dat <- as.numeric(datt[[n]])
  TT  <- length(dat)

  QSS_GARCH[[n]]  <- list()
  QSS_GAS[[n]] <- list()
  QSS_QbSD[[n]] <- list()
  QSS_AL[[n]] <- list()

  ALLS_GARCH[[n]]  <- list()
  ALLS_GAS[[n]] <- list()
  ALLS_QbSD[[n]] <- list()
  ALLS_AL[[n]] <- list()

  for (w in 2) {  # <— run only w = 2  (= 1250)

    win <- rw[w]
    fp  <- TT - win
    yf  <- dat[(win + 1):TT]
    
    QSS_GARCH[[n]][[w]] <- array(0, c(fp, length(gmo) * length(LOC), length(ALP)), dimnames = list(NULL, GARCHnam, ALP))
    QSS_GAS[[n]][[w]]  <- array(0, c(fp, length(LOC), length(ALP)), dimnames = list(NULL, LOC, ALP))
    QSS_QbSD[[n]][[w]] <- list()
    QSS_AL[[n]][[w]] <- array(0, c(fp, length(LOC) * length(ESD) * length(CAV), length(ALP)), dimnames = list(NULL, ALnam, ALP))

    ALLS_GARCH[[n]][[w]] <- array(0, c(fp, length(gmo) * length(LOC), length(ALP)), dimnames = list(NULL, GARCHnam, ALP))
    ALLS_GAS[[n]][[w]] <- array(0, c(fp, length(LOC), length(ALP)), dimnames = list(NULL, LOC, ALP))
    ALLS_QbSD[[n]][[w]] <- list()
    ALLS_AL[[n]][[w]] <- array(0, c(fp, length(LOC) * length(ESD) * length(CAV), length(ALP)), dimnames = list(NULL, ALnam, ALP))

    for (a in 1:length(ALP)) {

      print(c(bnam[n], rw[w], ALP[a]))
      garch  <- GARCHF[[n]][[w]][,,,,a]
      QbSD <- QbSDF[[n]][[w]][,,,, a] 
      AL <- ALF[[n]][[w]][,,,,,a]

      ## GARCH 
      for (s in 1:length(LOC)) {
        for (g in 1:length(gmo)) {
          VAR    <- garch[, g, 1, s]
          ES     <- garch[, g, 2, s]
          qscore <- (yf - VAR) * (ALP[a] - (yf <= VAR) * 1)
          QSS_GARCH[[n]][[w]][, GARCHmat[g, s], a]  <- qscore
          ES <- pmin(ES, -1e-10)          
          ALLS_GARCH[[n]][[w]][, GARCHmat[g, s], a] <- -log((ALP[a] - 1) / ES) - qscore / (ALP[a] * ES)
        }
      }

      ## GAS 

      GAS <- GASF[[n]][[w]][,,,a]   
      for (s in 1:length(LOC)) {
        VAR    <- GAS[, 1, s]
        ES     <- GAS[, 2, s]
        qscore <- (yf - VAR) * (ALP[a] - (yf <= VAR) * 1)
        QSS_GAS[[n]][[w]][, s, a]  <- qscore
        ES <- pmin(ES, -1e-10)
        ALLS_GAS[[n]][[w]][, s, a] <- -log((ALP[a] - 1) / ES) - qscore / (ALP[a] * ES)
      }
      
      ## QbSD
      QSS_QbSD[[n]][[w]][[a]] <- array(0, c(fp, length(CAV) * length(RLOC)), dimnames = list(NULL, QbSDnam) )
      ALLS_QbSD[[n]][[w]][[a]] <- array(0, c(fp, length(CAV) * length(RLOC)), dimnames = list(NULL, QbSDnam))

      for (s in 1:length(RLOC)) {
        for (v in 1:length(CAV)) {
          VAR    <- QbSD[, 1, v, s]
          ES     <- QbSD[, 2, v, s]
          dr     <- QbSDmat[v, s]
          qscore <- (yf - VAR) * (ALP[a] - (yf <= VAR) * 1)
          QSS_QbSD[[n]][[w]][[a]][, dr]  <- qscore
          ES <- pmin(ES, -1e-10)
          ALLS_QbSD[[n]][[w]][[a]][, dr] <- -log((ALP[a] - 1) / ES) - qscore / (ALP[a] * ES)
        }
      }

      ## AL
      for (s in 1:length(LOC)) {
        for (v in 1:length(CAV)) {
          for (d in 1:length(ESD)) {
            VAR    <- AL[, 1, d, v, s]
            ES     <- AL[, 2, d, v, s]
            dn     <- ALmat[d, v, s]
            qscore <- (yf - VAR) * (ALP[a] - (yf <= VAR) * 1)
            QSS_AL[[n]][[w]][, dn, a]  <- qscore
            ES <- pmin(ES, -1e-10)
            ALLS_AL[[n]][[w]][, dn, a] <- -log((ALP[a] - 1) / ES) - qscore / (ALP[a] * ES)
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
  for (w in 2) {  # <— run only w = 2
    MCSS_QS[[n]][[w]] <- list()
    MCSS_AL[[n]][[w]] <- list()

    for (a in 1:length(ALP)) {
      print(c(bnam[n], rw[w], ALP[a]))

      QSS <- cbind(
        QSS_GARCH[[n]][[w]][,,a],
        QSS_GAS[[n]][[w]][, , a],       
        QSS_AL[[n]][[w]][,,a],
        QSS_QbSD[[n]][[w]][[a]]
      )
      colnames(QSS) <- Mnam
      
      ALLS <- cbind(
        ALLS_GARCH[[n]][[w]][,,a],
        ALLS_GAS[[n]][[w]][, , a],      
        ALLS_AL[[n]][[w]][,,a],
        ALLS_QbSD[[n]][[w]][[a]]
      )
      colnames(ALLS) <- Mnam
      
      k <- dim(QSS)[2]
      MCSS_QS[[n]][[w]][[a]] <- MCSprocedure(Loss = QSS,  alpha = 0.1, statistic = "TR", B = 1000)
      MCSS_AL[[n]][[w]][[a]] <- MCSprocedure(Loss = ALLS, alpha = 0.1, statistic = "TR", B = 1000)
    }
  }
}

############################################


bnam <- c("SP500", "DJIA", "NASDAQ", "STOXX50", "FTSE", "DAX", "CAC40", "TSX")

nb  <- length(bnam)
rw  <- c(250, 1250, 2500)
ALP <- c(0.01, 0.025, 0.05)

########### Produce Tables 8,9,10

MCSR_AL <- array( 0, c(length(Mnam1), nb, length(ALP), length(rw)), dimnames = list(Mnam1, bnam, ALP, rw))

for (n in 1:nb) {
  for (w in 2) {  # <— run only w = 2
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
w <- 2
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

write.table(Table8,  file = "Table8.txt", sep = "\t", quote = FALSE, na = "", col.names = NA)
write.table(Table9,  file = "Table9.txt", sep = "\t", quote = FALSE, na = "", col.names = NA)
write.table(Table10, file = "Table10.txt", sep = "\t", quote = FALSE, na = "", col.names = NA)



########### Produce Tables 5,6,7

MCSR_QS <- array(
  0, c(length(Mnam1), nb, length(ALP), length(rw)),
  dimnames = list(Mnam1, bnam, ALP, rw)
)

for (n in 1:nb) {
  for (w in 2) {  # <— run only w = 2
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
w <- 2
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

write.table(Table5, file = "Table5.txt", sep = "\t", quote = FALSE, na = "", col.names = NA)
write.table(Table6, file = "Table6.txt", sep = "\t", quote = FALSE, na = "", col.names = NA)
write.table(Table7, file = "Table7.txt", sep = "\t", quote = FALSE, na = "", col.names = NA)

