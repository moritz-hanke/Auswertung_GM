#############################
### Laden der Datensätze/ des Datensatzes,
  # der ausgewertet werden soll

load()
# Daten müssen "erg" genannt werden



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##############################
### Benötigte Funktionen:

mu <- function( edges, vertices, threshold=0.5, epsilon=1e-100)
  ###### edges ist Vektor von Kantenh?ufigkeiten
  ######       falls notwendig kann aus einer Kantenmatrix dieser so erzeugt werden:
  ######       edges <- t(edgematrix)[lower.tri(t(edgematrix))]
  ###### vertices ist Anzahl von Knoten  
{ 
  e <- cbind( threshold*edges, (1-threshold)*(1-edges))
  edgeS <- apply(e,1,min)
  edgeS
  
  meu  <- sum(edgeS/ (threshold*(1-threshold)))/choose(vertices,2)
  mseu <- sum((edgeS/(threshold*(1-threshold)))^2)/choose(vertices,2)
  H    <- -1/log(2)*( edges*log(edges+epsilon) + (1-edges)*log(1-edges+epsilon) )
  mee  <- sum(H) / choose(vertices,2)
  out  <- list(meu=meu, mseu=mseu, mee=mee)
  out
}

######################################CLUSTER: extrahiere Cliquen (Modell)
cli <- function(x)
  #### Cluster: extrahiere Cliquen des Originalmodels
  #### cli( erg[[1]]$original )
  #### cli( erg[[1]]$bootmodel[[1]] )
{
  out <- lapply( x$glist, paste, collapse="")
  out <- sort(unlist( out ))
  out <- tolower( paste( out, collapse=",") )
  out
}


############################################################# Matrixparse
############ Erstelle aus einem Stringvektor eine Adjazenzmatrix
`.gm.matrixparse` <- function (result) 
  ### results = stringvektor mit Cliquen wie "ab,acd,bcd"
{
  m = strsplit(result, "")[[1]]
  elements = unique(m)
  if (any(elements == ",")) 
    elements = elements[-which(elements == ",")]
  elements = sort(elements)
  dep.table = matrix(0, nrow = length(elements), ncol = length(elements))
  dimnames(dep.table) = list(elements, elements)
  for (i in 1:length(m)) {
    if (length(which(elements == m[i])) == 0) {
    }
    else {
      j = 1
      while (i + j <= length(m) && length(which(elements == 
                                                  m[i + j])) > 0) {
        pos1 = which(elements == m[i])
        pos2 = which(elements == m[i + j])
        dep.table[pos1, pos2] = 1
        j = j + 1
      }
    }
  }
  dep.table
}

############################################################# Modelparse
############ Erstelle aus einer Adjazenzmatrix einen Stringvektor
`.gm.modelparse` <-  function (data) 
  ### data ist eine Matrix
{
  v <- dim(data)[1]
  for (i in 1:(v - 1)) {
    for (j in (i + 1):v) {
      data[j, i] <- data[i, j]
    }
  }
  maxvar <- v
  SetSize <- function(X) {
    i <- 0
    while ((X[(i + 1)] > 0) && (i < maxvar)) {
      i <- i + 1
    }
    out <- i
  }
  Neighbours <- function(u, G) {
    j <- 0
    Y <- rep(0, maxvar)
    for (i in 1:maxvar) {
      if (G[i, u] == 1) {
        j <- j + 1
        Y[j] <- i
      }
    }
    out <- Y
  }
  AndSet <- function(X, Y) {
    nx <- SetSize(X)
    ny <- SetSize(Y)
    k <- 0
    Z <- rep(0, maxvar)
    if ((nx > 0) && (ny > 0)) {
      for (i in 1:nx) {
        for (j in 1:ny) {
          if (X[i] == Y[j]) {
            k <- k + 1
            Z[k] <- X[i]
          }
        }
      }
    }
    out <- Z
  }
  BK <- function(R, P, X, G) {
    result <- NULL
    if ((SetSize(P) == 0) && (SetSize(X) == 0)) {
      result <- rbind(result, R)
    }
    else {
      k <- SetSize(P)
      if (k == 0) {
      }
      else {
        for (i in k:1) {
          u <- P[i]
          P <- c(P[-i], 0)
          nr <- SetSize(R)
          if (nr == 0) {
            Rnew <- c(u, rep(0, (maxvar - 1)))
          }
          else {
            Rnew <- c(R[1:nr], u, rep(0, (maxvar - nr - 
                                            1)))
          }
          N <- Neighbours(u, G)
          Pnew <- AndSet(P, N)
          Xnew <- AndSet(X, N)
          result = rbind(result, BK(Rnew, Pnew, Xnew, 
                                    G))
          nx <- SetSize(X)
          if (nx == 0) {
            X <- c(u, rep(0, (maxvar - 1)))
          }
          else {
            X <- c(X[1:nx], u, rep(0, (maxvar - nx - 
                                         1)))
          }
        }
      }
    }
    result
  }
  X <- BK(rep(0, maxvar), 1:maxvar, rep(0, maxvar), data)
  result <- NULL
  for (i in 1:dim(X)[1]) {
    result <- c(result, ",", letters[sort(X[i, ])])
  }
  result <- paste(result[-1], collapse = "")
  result <- paste(sort(strsplit(result, ",")[[1]]), collapse = ",")
  result
}






#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#############################################################
#############################################################
### Funktion zur Auswertung der Ergebnisse:
  # Für die simulierten Datensatz wird eine Liste erstellt, 
  # in der das erste Element ein Dataframe ist mit den 
  # Ergebnissen aus den Differenzen für alle x Simulationen.
  # Das zweite Element ist ebenfalls ein Data.Frame, der die
  # Ergebnisse für alle x Simulationen bzgl. der Funktion mu()
  # beinhaltet.
  # Das dritte Element ist wiederum eine Liste, das für jede
  # der x Simulationen eine Vector mit den Samples der n 
  # Bootstraps beinhaltet.


Results <- function(input.data, MODEL){
  ### input.data sind die Ergenisse aus den x Simulationen
    # und n Bootstraps
  ### MODEL ist das wahre Model auf dessen Grundlage simuliert
    # wurde (truemodel)
  
  
  ### Hier ist noch keine Fehlerabfrage, weil die Ergebnisse in Listen 
    # relativ komplex bzw. tief geschachtelt sind und die passende Abfrage
    # längere Zeit zum Programmieren bräuchte um alle möglichen Fehler
    # abzufangen. Für den aktuellen Fall (Stand 2014-05-09) funktioniert die
    # Funktion
  
  Auswertung <- lapply(seq_along(input.data), FUN=function(x){
    
    
    
    edge.rel.freq <- input.data[[x]]$edgefreq/input.data[[x]]$R
    meanmodel <- ifelse(edge.rel.freq >=0.5, 1, 0)
    meanmodel[lower.tri(meanmodel)] <- 0
    
    cli.original <- cli(input.data[[x]]$original)
    origmodel <- .gm.matrixparse(cli.original)
    
    
    truemodel <- .gm.matrixparse(MODEL)
    
    
    bm <- sapply(c(1:input.data[[x]]$R), FUN=function(r){
      cli(input.data[[x]]$bootmodel[[r]])
    }
    )
    
    
    
    
    mu.etc <-NULL
    edges <- edge.rel.freq[lower.tri(edge.rel.freq)]
    mu.etc <- unlist(mu(edges, nrow(edge.rel.freq)))
    
    
    
    
    
    ### Differenzen
    diffmeanbm <- difforigbm <- difftruebm <- NULL
    for(k in 1:length(bm))
    {
      diffmeanbm[k] <- sum(abs( meanmodel - .gm.matrixparse(bm[k])))  
      difforigbm[k] <- sum(abs( origmodel - .gm.matrixparse(bm[k])))  
      difftruebm[k] <- sum(abs( truemodel - .gm.matrixparse(bm[k])))  
    } 
    Dmeanbm <- sum(diffmeanbm)/length(bm)
    Dorigbm <- sum(difforigbm)/length(bm)
    Dtruebm <- sum(difftruebm)/length(bm)
    PI95mm  <- quantile(diffmeanbm, prob=c(.05,.95))
    PI95orig<- quantile(difforigbm, prob=c(.05,.95))
    PI95true<- quantile(difftruebm, prob=c(.05,.95))
    Dmean = sum(abs( meanmodel - truemodel))
    Dorig = sum(abs( origmodel - truemodel))
    
    OUT <- list(data.frame("Truemodel"= MODEL, "Originalmodel"= .gm.modelparse(origmodel), "Meanmodel"= .gm.modelparse(meanmodel), 
                           "Dtrue_boot"=Dtruebm, "PI95ut"=PI95true[1], "PI95ot"=PI95true[2], 
                           "Dmean_boot"=Dmeanbm, "PI95um"=PI95mm[1], "PI95om"=PI95mm[2], 
                           "Dorig_boot"=Dorigbm, "PI95uo"=PI95orig[1], "PI95oo"=PI95orig[2], 
                           "Dmean"=Dmean, "Dorig"=Dorig),
                "MU"=mu.etc,modelle=bm)
    OUT
    
    
    
  })
  
  
  ### Data.Frame der für alle x Simulationen bzgl Differenzen
  Differences <- lapply(seq_along(Auswertung), FUN=function(var){
    Auswertung[[var]][[1]]
  })
  Differences <- do.call("rbind", Differences)
  row.names(Differences) <- NULL   #sonst sind die rownames 5%
  
  ### Data.Frame der für alle x Simulationen bzgl. Ergebnisse 
  #der Funktion mu()
  MUs <- lapply(seq_along(Auswertung), FUN=function(var){
    Auswertung[[var]][[2]]
  })
  MUs <- do.call("rbind", MUs)
  
  ### Liste mit den n Bootstrap-Samples für die x Simulationen
  BS_Samples <- lapply(seq_along(Auswertung), FUN=function(var){
    Auswertung[[var]][[3]]
  })
  
  
  ergebnisse <- list(Differences, MUs, BS_Samples)
  names(ergebnisse) <- c("Differences", "MUs", "BS_Samples")
  return(ergebnisse)



}


Test <- Results(erg, "ace,acg,adg,bdh,be,cef,cfg,fjk,i")

