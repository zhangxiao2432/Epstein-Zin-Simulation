#--------------------------------------------------------------------
# Tauchen and Hussey(1991)
# author: Zhonghui Zhang
# This code is a translation of the Fortran77 code by Tauchen
#
#
# REFFERENCE:
#
#    TAUCHEN, GEORGE, 1987, "QUADRATURE-BASED METHODS FOR OBTAINING
#    APPROXIMATE SOLUTIONS TO NONLINEAR ASSET PRICING MODELS,"
#    DUKE UNIVERSITY.
#
# THIS PROGRAM CALCULATES A DISCRETE APPROXIMATION TO A CONTINUOUS 
# GAUSSIAN VECTOR AUTREGRESSION
#
# Y(T) = B + A(1)*Y(T-1) + A(2)*Y(T-2) + ... + A(NLAG)*Y(T-NLAG) + E(T)
#
# WHERE
#
#    Y(T) IS AN (NVAR X 1) VECTOR
#    B IS AN (NVAR X 1) VECTOR
#    A(I) ARE (NVAR X NVAR) MATRICES
#    E(T) IS AN (NVAR X 1) VECTOR OF I.I.D. MULTIVARIATE NORMAL 
#       RANDOM VARIBLES WITH MEAN VECTOR ZERO AND VARIANCE SIGMA 

# INPUTS:
#
#  NVAR - NUMBER OF VARIABLES IN THE VAR 
#  NLAG - NUMBER OF LAGS IN THE VAR 
#  NVAL - NUMBER OF DISCRETE POINTS FOR EACH VARIABLE 
#         NVAR ROWS 
#     B - VECTOR OF INTERCEPTS IN THE VAR 
#         Matrix(nrow = NVAR, 1)  
#  AUTO - HORIZONTALLY CONCATENATED MATRICES 
#         Matrix(nrow = NLAG, ncol = NVAR*NVAR)
# SIGMA - VARIANCE MATRIX OF E(T) READ COLUMNWISE
#         Matrix(nrow = NVAR, ncol = NVAR)
#
#    NS - NUMBER OF DISCRETE VALUES THAT Y(T) CAN TAKE ON;          
#         EQUAL TO THE PRODUCT OF THE ELEMENTS OF THE MATRIX 
#         NVAL.  
# 
#NSTATE - NUMBER OF STATES IN THE SYSTEM; EQUAL TO NS**NLAG.
#         THE STATE IS DEFINED BY THE DISCRETE VALUES OF
#         Y(T-1), Y(T-2), ...,Y(T-NLAG); SINCE EACH LAG CAN
#         TAKE ON NS DIFFERENT VALUES, THERE ARE NS**NLAG
#         STATES IN THE SYSTEM.
#
# EXAMPLE:  IF THERE ARE TWO VARIABLES (Y1 AND Y2) AND TWO
#            LAGS IN THE VAR AND TWO DISCRETE POINTS ARE USED FOR 
#            EACH VARIABLE, THE STATES ARE ARRANGED AS FOLLOWS:
#
#                                   DISCRETE VALUES
#                      
#            STATE NO.        TIME T-1            TIME T-2
#            ---------     ---------------     ---------------
#                           
#                1         Y1(1)     Y2(1)     Y1(1)     Y2(1)
#                2         Y1(2)     Y2(1)     Y1(1)     Y2(1)
#                3         Y1(1)     Y2(2)     Y1(1)     Y2(1)
#                4         Y1(2)     Y2(2)     Y1(1)     Y2(1)
#                5         Y1(1)     Y2(1)     Y1(2)     Y2(1)
#                6         Y1(2)     Y2(1)     Y1(2)     Y2(1)
#                7         Y1(1)     Y2(2)     Y1(2)     Y2(1)
#                8         Y1(2)     Y2(2)     Y1(2)     Y2(1)
#                9         Y1(1)     Y2(1)     Y1(1)     Y2(2)
#               10         Y1(2)     Y2(1)     Y1(1)     Y2(2)
#               11         Y1(1)     Y2(2)     Y1(1)     Y2(2)
#               12         Y1(2)     Y2(2)     Y1(1)     Y2(2)
#               13         Y1(1)     Y2(1)     Y1(2)     Y2(2)
#               14         Y1(2)     Y2(1)     Y1(2)     Y2(2)
#               15         Y1(1)     Y2(2)     Y1(2)     Y2(2)
#               16         Y1(2)     Y2(2)     Y1(2)     Y2(2)
#
# THREE OUTPUT MATRICES ARE RELEVANT
#
#  YMAT - (NS X NVAR) MATRIX OF DISCRETE VARIABLE VALUES
#
# PISTA - (NSTATE X 1) VECTOR THAT IS THE STATIONARY
#         PROBABILITY DISTRIBUTION OF DISCRETE STATES
#
# PIMAT - (NSTATE X NS) MATRIX OF TRANSITION PROBABILITIES.  
#         IF NLAG>1, THIS MATRIX IS NOT SQUARE BECAUSE FROM
#         ANY GIVEN STATE, ONLY NS OTHER STATES ARE REACHABLE
#         IN THE NEXT PERIOD.
#
# ALSO AVAILABLE IN THE SUBROUTINE MCHAIN IS THE (NSTATE X NS)
# MATRIX IREACH WHICH CONTAINS THE REACHABLE STATE NUMBERS.  THUS
# ROW I OF IREACH CONTAINS THE NUMBERS OF THE NS STATES THAT ARE 
# REACHABLE IN THE NEXT PERIOD WHEN THE STATE NUMBER IN THE CURRENT
# PERIOD IS I.
#====================================================================
#====================================================================

Tauchen <- function(NVAR,NLAG,NVAL,B,AUTO,SIGMA) {

library("MASS")
  
#CALCULATE THE NUMBER OF STATES
  NS <- 1
  for(J in 1:NVAR){
    NS <- NS*NVAL[J]
  }
  
#ERROR TRAPS
  if (NVAR > 2){
    stop("Choose NVAR = 1 or 2.") 
  }
  if (NLAG > 2) {
    stop("Choose NLAG = 1 or 2.") 
  }
  if (NS > 20) {
    stop("States can not be more than 20.") 
  }
  
#QUADRATURE DATA POINTS
#--------------------------------------------------------------------
  QUADAT <- c(1,    0.00000000000000,    1.00000000000000,
              2,   -1.00000000000000,    0.50000000000000,
              2,    1.00000000000000,    0.50000000000000,
              3,   -1.73205080756887,    0.16666666666667,
              3,   -0.00000000000000,    0.66666666666667,
              3,    1.73205080756888,    0.16666666666667,
              4,   -2.33441421833897,    0.04587585476807,
              4,   -0.74196378430273,    0.45412414523193,
              4,    0.74196378430272,    0.45412414523193,
              4,    2.33441421833898,    0.04587585476807,
              5,   -2.85697001387280,    0.01125741132772,
              5,   -1.35562617997427,    0.22207592200561,
              5,   -0.00000000000000,    0.53333333333333,
              5,    1.35562617997426,    0.22207592200561,
              5,    2.85697001387280,    0.01125741132772,
              6,   -3.32425743355212,    0.00255578440206,
              6,   -1.88917587775371,    0.08861574604191,
              6,   -0.61670659019259,    0.40882846955603,
              6,    0.61670659019259,    0.40882846955603,
              6,    1.88917587775371,    0.08861574604191,
              6,    3.32425743355212,    0.00255578440206,
              7,   -3.75043971772574,    0.00054826885597,
              7,   -2.36675941073454,    0.03075712396759,
              7,   -1.15440539473997,    0.24012317860501,
              7,   -0.00000000000000,    0.45714285714286,
              7,    1.15440539473997,    0.24012317860501,
              7,    2.36675941073454,    0.03075712396759,
              7,    3.75043971772574,    0.00054826885597,
              8,   -4.14454718612589,    0.00011261453838,
              8,   -2.80248586128754,    0.00963522012079,
              8,   -1.63651904243511,    0.11723990766176,
              8,   -0.53907981135138,    0.37301225767908,
              8,    0.53907981135137,    0.37301225767908,
              8,    1.63651904243510,    0.11723990766176,
              8,    2.80248586128754,    0.00963522012079,
              8,    4.14454718612589,    0.00011261453838,
              9,   -4.51274586339978,    0.00002234584401,
              9,   -3.20542900285647,    0.00278914132123,
              9,   -2.07684797867782,    0.04991640676522,
              9,   -1.02325566378913,    0.24409750289494,
              9,   -0.00000000000000,    0.40634920634921,
              9,    1.02325566378913,    0.24409750289494,
              9,    2.07684797867783,    0.04991640676522,
              9,    3.20542900285647,    0.00278914132123,
              9,    4.51274586339977,    0.00002234584401,
              10,   -4.85946282833230,    0.00000431065263,
              10,   -3.58182348355192,    0.00075807093431,
              10,   -2.48432584163895,    0.01911158050077,
              10,   -1.46598909439116,    0.13548370298027,
              10,   -0.48493570751550,    0.34464233493202,
              10,    0.48493570751550,    0.34464233493202,
              10,    1.46598909439115,    0.13548370298027,
              10,    2.48432584163895,    0.01911158050077,
              10,    3.58182348355192,    0.00075807093431,
              10,    4.85946282833231,    0.00000431065263,
              11,   -5.18800122437487,    0.00000081218498,
              11,   -3.93616660712997,    0.00019567193027,
              11,   -2.86512316064364,    0.00672028523554,
              11,   -1.87603502015484,    0.06613874607106,
              11,   -0.92886899738106,    0.24224029987397,
              11,   -0.00000000000000,    0.36940836940837,
              11,    0.92886899738106,    0.24224029987397,
              11,    1.87603502015484,    0.06613874607106,
              11,    2.86512316064364,    0.00672028523554,
              11,    3.93616660712998,    0.00019567193027,
              11,    5.18800122437486,    0.00000081218498,
              12,   -5.50090170446775,    0.00000014999272,
              12,   -4.27182584793227,    0.00004837184923,
              12,   -3.22370982877009,    0.00220338068753,
              12,   -2.25946445100079,    0.02911668791236,
              12,   -1.34037519715162,    0.14696704804533,
              12,   -0.44440300194414,    0.32166436151283,
              12,    0.44440300194414,    0.32166436151283,
              12,    1.34037519715161,    0.14696704804533,
              12,    2.25946445100079,    0.02911668791236,
              12,    3.22370982877009,    0.00220338068753,
              12,    4.27182584793227,    0.00004837184923,
              12,    5.50090170446774,    0.00000014999272,
              13,   -5.80016725238649,    0.00000002722628,
              13,   -4.59139844893651,    0.00001152659653,
              13,   -3.56344438028163,    0.00068123635044,
              13,   -2.62068997343221,    0.01177056050600,
              13,   -1.72541837958824,    0.07916895586045,
              13,   -0.85667949351945,    0.23787152296414,
              13,   -0.00000000000000,    0.34099234099234,
              13,    0.85667949351945,    0.23787152296414,
              13,    1.72541837958824,    0.07916895586045,
              13,    2.62068997343221,    0.01177056050600,
              13,    3.56344438028163,    0.00068123635044,
              13,    4.59139844893651,    0.00001152659653,
              13,    5.80016725238649,    0.00000002722628,
              14,   -6.08740954690129,    0.00000000486816,
              14,   -4.89693639734556,    0.00000266099134,
              14,   -3.88692457505976,    0.00020033955376,
              14,   -2.96303657983866,    0.00442891910695,
              14,   -2.08834474570194,    0.03865010882425,
              14,   -1.24268895548546,    0.15408333984251,
              14,   -0.41259045795460,    0.30263462681302,
              14,    0.41259045795460,    0.30263462681302,
              14,    1.24268895548546,    0.15408333984251,
              14,    2.08834474570194,    0.03865010882425,
              14,    2.96303657983866,    0.00442891910695,
              14,    3.88692457505976,    0.00020033955376,
              14,    4.89693639734555,    0.00000266099134,
              14,    6.08740954690128,    0.00000000486816,
              15,   -6.36394788882983,    0.00000000085896,
              15,   -5.19009359130477,    0.00000059754196,
              15,   -4.19620771126901,    0.00005642146405,
              15,   -3.28908242439876,    0.00156735750355,
              15,   -2.43243682700975,    0.01736577449214,
              15,   -1.60671006902873,    0.08941779539984,
              15,   -0.79912906832455,    0.23246229360973,
              15,   -0.00000000000000,    0.31825951825952,
              15,    0.79912906832454,    0.23246229360973,
              15,    1.60671006902873,    0.08941779539984,
              15,    2.43243682700975,    0.01736577449214,
              15,    3.28908242439876,    0.00156735750355,
              15,    4.19620771126901,    0.00005642146405,
              15,    5.19009359130477,    0.00000059754196,
              15,    6.36394788882983,    0.00000000085896,
              16,   -6.63087819839313,    0.00000000014978,
              16,   -5.47222570594934,    0.00000013094732,
              16,   -4.49295530252000,    0.00001530003216,
              16,   -3.60087362417154,    0.00052598492657,
              16,   -2.76024504763069,    0.00726693760118,
              16,   -1.95198034571633,    0.04728475235401,
              16,   -1.16382910055496,    0.15833837275095,
              16,   -0.38676060450056,    0.28656852123801,
              16,    0.38676060450056,    0.28656852123801,
              16,    1.16382910055496,    0.15833837275095,
              16,    1.95198034571633,    0.04728475235401,
              16,    2.76024504763069,    0.00726693760118,
              16,    3.60087362417154,    0.00052598492657,
              16,    4.49295530252001,    0.00001530003216,
              16,    5.47222570594934,    0.00000013094732,
              16,    6.63087819839312,    0.00000000014978,
              17,   -6.88912243989532,    0.00000000002584,
              17,   -5.74446007865941,    0.00000002808016,
              17,   -4.77853158962998,    0.00000401267945,
              17,   -3.90006571719800,    0.00016849143155,
              17,   -3.07379717532819,    0.00285894606228,
              17,   -2.28101944025299,    0.02308665702571,
              17,   -1.50988330779674,    0.09740637116272,
              17,   -0.75184260070390,    0.22670630846898,
              17,   -0.00000000000000,    0.29953837012661,
              17,    0.75184260070389,    0.22670630846898,
              17,    1.50988330779674,    0.09740637116272,
              17,    2.28101944025298,    0.02308665702571,
              17,    3.07379717532818,    0.00285894606228,
              17,    3.90006571719801,    0.00016849143155,
              17,    4.77853158962997,    0.00000401267945,
              17,    5.74446007865940,    0.00000002808016,
              17,    6.88912243989532,    0.00000000002584,
              18,   -7.13946484914647,    0.00000000000442,
              18,   -6.00774591135959,    0.00000000590549,
              18,   -5.05407268544274,    0.00000102155240,
              18,   -4.18802023162940,    0.00005179896144,
              18,   -3.37473653577809,    0.00106548479629,
              18,   -2.59583368891123,    0.01051651775194,
              18,   -1.83977992150864,    0.05489663248022,
              18,   -1.09839551809150,    0.16068530389351,
              18,   -0.36524575550770,    0.27278323465429,
              18,    0.36524575550770,    0.27278323465429,
              18,    1.09839551809150,    0.16068530389351,
              18,    1.83977992150864,    0.05489663248022,
              18,    2.59583368891124,    0.01051651775194,
              18,    3.37473653577809,    0.00106548479629,
              18,    4.18802023162940,    0.00005179896144,
              18,    5.05407268544274,    0.00000102155240,
              18,    6.00774591135959,    0.00000000590549,
              18,    7.13946484914647,    0.00000000000442,
              19,   -7.38257902403042,    0.00000000000075,
              19,   -6.26289115651325,    0.00000000122037,
              19,   -5.32053637733603,    0.00000025322200,
              19,   -4.46587262683102,    0.00001535114595,
              19,   -3.66441654745063,    0.00037850210941,
              19,   -2.89805127651575,    0.00450723542034,
              19,   -2.15550276131693,    0.02866669103012,
              19,   -1.42887667607837,    0.10360365727614,
              19,   -0.71208504404238,    0.22094171219914,
              19,   -0.00000000000000,    0.28377319275152,
              19,    0.71208504404238,    0.22094171219914,
              19,    1.42887667607837,    0.10360365727614,
              19,    2.15550276131693,    0.02866669103012,
              19,    2.89805127651575,    0.00450723542034,
              19,    3.66441654745063,    0.00037850210941,
              19,    4.46587262683103,    0.00001535114595,
              19,    5.32053637733603,    0.00000025322200,
              19,    6.26289115651324,    0.00000000122037,
              19,    7.38257902403042,    0.00000000000075,
              20,   -7.61904854167975,    0.00000000000013,
              20,   -6.51059015701364,    0.00000000024821,
              20,   -5.57873880589319,    0.00000006127490,
              20,   -4.73458133404605,    0.00000440212109,
              20,   -3.94396735065731,    0.00012882627996,
              20,   -3.18901481655338,    0.00183010313108,
              20,   -2.45866361117236,    0.01399783744710,
              20,   -1.74524732081412,    0.06150637206398,
              20,   -1.04294534880275,    0.16173933398400,
              20,   -0.34696415708136,    0.26079306344955,
              20,    0.34696415708135,    0.26079306344955,
              20,    1.04294534880275,    0.16173933398400,
              20,    1.74524732081412,    0.06150637206398,
              20,    2.45866361117236,    0.01399783744710,
              20,    3.18901481655338,    0.00183010313108,
              20,    3.94396735065731,    0.00012882627996,
              20,    4.73458133404605,    0.00000440212109,
              20,    5.57873880589319,    0.00000006127490,
              20,    6.51059015701366,    0.00000000024821,
              20,    7.61904854167974,    0.00000000000013)
  
  QUADAT <- t(matrix(QUADAT,3,210))
#--------------------------------------------------------------------

#COMPUTE THE INVERSE OF SIGMA
SIGINV <- solve(SIGMA)
  
#--------------------------------------------------------------------
#FUNCTION TO PUT ALL POSSIBLE COMBINATIONS OF QUADRATURE ABSCISSA 
#AND CORRESPONDING PRODUCTS OF QUADRATURE WEIGHTS IN RETURNED MATRIX                                                 
#--------------------------------------------------------------------
QUAD <- function(NVAL,NVAR,NS,QUADAT,NRQUAD){
  
  WVEC <- matrix(1,NS,1)
  ZMAT <- matrix(5,NS,NVAR)
  NPCUM <- 1
  for (I in 1:NVAR){
    NP <- NVAL[I] 
    for (J in 1:NP){
      Z <- QUADAT[((NP-1)*NP)/2+J,2]
      W <- QUADAT[((NP-1)*NP)/2+J,3]
      temp <- (NS/NP)/NPCUM
      for (K in 1: temp){
        for (L in 1:NPCUM){
          ZMAT[(K-1)*NPCUM*NP+(J-1)*NPCUM+L,I] <- Z
          WVEC[(K-1)*NPCUM*NP+(J-1)*NPCUM+L] <- WVEC[(K-1)*NPCUM*NP+(J-1)*NPCUM+L]*W
        }
      }
    }
    NPCUM <- NPCUM*NP 
  } 
  remove("temp")
  return(list(ZMAT=ZMAT, WVEC=WVEC))
}  
  
QUAD <- QUAD(NVAL,NVAR,NS,QUADAT,NRQUAD)  
ZMAT <- QUAD$ZMAT
WVEC <- QUAD$WVEC

#FIND THE CHOLESKY FACTORIZATION OF SIGMA 
ARMEAN <- function(AUTO,B,NVAR,NLAG){
  
  SUMAUT <- matrix(0,NVAR,NVAR)
  for (L in 1:NLAG){
    for (I in 1:NVAR){
      for (J in 1:NVAR){
        SUMAUT[I,J] <- SUMAUT[I,J]+AUTO[I,J+(L-1)*NVAR] 
      }
    }
  }
  for (M in 1:NVAR){
    for (N in 1:NVAR){
      SUMAUT[M,N] <- -SUMAUT[M,N]
      if (M == N) SUMAUT[M,N] <- 1+SUMAUT[M,N]
    }
  }
  
  temp <- ginv(SUMAUT)
  EXVAL <- temp %*% B

  return(EXVAL=EXVAL)
}

EXVAL <- ARMEAN(AUTO,B,NVAR,NLAG)

ZTOY <- function(SIGMA,ZMAT,NVAR,NS){
  
  S <- chol(SIGMA)
  YMAT <- ZMAT %*% S  
  
  for (I in 1:NS){
    for (J in 1:NVAR){
      YMAT[I,J] <- YMAT[I,J]+EXVAL[J]
    }
  }
  DETSIG=1
  for (K in 1:NVAR){
    DETSIG=DETSIG*S[K,K]
  }
  return(list(YMAT=YMAT,DETSIG=DETSIG))
}

YMAT <- ZTOY(SIGMA,ZMAT,NVAR,NS)$YMAT
DETSIG <- ZTOY(SIGMA,ZMAT,NVAR,NS)$DETSIG

#--------------------------------------------------------------------
#FUNCTION TO CODE EACH ROW OF IMAT AS                             
#X = N1 + (N2-1)*NBASE + ... + (NL-1)*NBASE(L-1)               
#WHERE (N1,N2,...,NL) IS A ROW OF IMAT AND L = COLS(IMAT).          
#--------------------------------------------------------------------
CODE <- function(IMAT,NR,NC,NBASE){
  
  ICODE <- matrix(0,NR,1)
  if (NBASE < 1){
    stop("CODE ERROR: BAD NBASE!")
  }
  for (I in 1:NR){
    for (J in 1:NC){
      if (IMAT[I,J] > NBASE | IMAT[I,J] < 0){
        stop("CODE ERROR: AN ELEMENT OF THE INPUT ARRAY EITHER EXCEEDS THE NBASE OR IS LESS THAN OR EQUAL TO ZERO")
      }
    }
  }
  for (M in 1:NR){
    ICODE[M] <- IMAT[M,1]
    if(NC > 1){
      for (N in 2:NC){
        ICODE[M] <- ICODE[M]+(IMAT[M,N]-1)*(NBASE^(N-1))
      }
    }
  }
  return( ICODE= ICODE)
}
#--------------------------------------------------------------------
#FUNCTION TO DECODE THE ROWS OF IX WHERE A TYPICAL ROW IS OF THE 
#FORM
#IX(K,.) =  N1 + (N2-1)*NBASE**1 + ... (NL-1)*NBASE**(L-1)        
#--------------------------------------------------------------------
DECODE <- function(IX,NR,NBASE,L){
  if (NBASE < 0 | L < 1){
    stop("DECODE ERROR BAD NBASE OR BAD LENGTH!")
  }
  for (I in 1:NR){
    if (IX[I] < 0){
      stop("DECODE WARNING: AT LEAST ONE ELEMENT OF INPUT VECTOR IS LESS THAN OR EQUAL TO ZERO!")
    }
  }
  ISTATE <-matrix(0,NR,L)
  for (M in 1:NR){
    for (N in 1:L){
      ISTATE[M,N] <- as.integer(as.integer((IX[M]-1)/(NBASE^(N-1)))%%NBASE)+1
    }
  }
  return(ISTATE=ISTATE)
}
#--------------------------------------------------------------------
#FUNCTION REACHR                                     
#--------------------------------------------------------------------
REACHR <- function(IVEC,NBASE,NLAG){
  IX <- matrix(0,NS,NLAG)
  for (I in 1:NBASE){
    IX[I,1] <- I
    if (NLAG>1){
      for (J in 2:NLAG){
      IX[I,J] <- IVEC[J-1]
      }
    }
  }
  IY <- CODE(IX,NBASE,NLAG,NBASE)
  return(IY=IY)
}
#--------------------------------------------------------------------
#--------------------------------------------------------------------
#FUNCTION FOR ROW-WISE MULTIVARIATE NORMAL DISTRIBUTION                  
#NOTES: THE ELEMENTS OF RETURNED VECTOR ARE THE MULTIVARIATE     
#NORMAL DENSITY N(MU,SIGMA) EVALUATED AT EACH ROW OF XMAT 
#--------------------------------------------------------------------
PHIMAT <- function(XMAT,MU,SIGINV,DETSIG,N,NVAR){
  XCTR <- matrix(1,N,NVAR)
  for (I in 1:N){
    for (J in 1:NVAR){
      XCTR[I,J] <- XMAT[I,J]-MU[J]
    }
  }
  XSIG <- XCTR %*% SIGINV
  FZMAT <- matrix(0,N,1)
  for (L in 1:N){
    for (M in 1:NVAR){
      FZMAT[L] <- FZMAT[L]+XSIG[L,M]*XCTR[L,M]
    }
    FZMAT[L] <- exp(-0.5*FZMAT[L])/((2*pi*DETSIG)^.5)
  }
  return(FZMAT=FZMAT)
}
#--------------------------------------------------------------------
#FUNCTION FOR TRANSITION PROBABILITES AS A FUNCTION OF X                
#--------------------------------------------------------------------
PIFN <- function(XMAT,NRXMAT,YMAT,NS,NVAR,WVEC,EXVAL,SIGINV,DETSIG,AUTO,B,NLAG){
  I1 <- 1
  I2 <- I1+NS*NVAR
  
  FZMAT <- PHIMAT(YMAT,EXVAL,SIGINV,DETSIG,NS,NVAR)  
  
  TEMP <- matrix(0,NS,1)
  PIVAL <- matrix(0,NRXMAT,NS)
  CONEXP <- matrix(0,NVAR,1)
  
  for (I in 1:NS){
    TEMP[I] <- WVEC[I]/FZMAT[I]
  }
  
  for (J in 1:NRXMAT){
    
    CONEXP <- t(B) + XMAT%*%t(AUTO)
    FZMAT <- PHIMAT(YMAT,CONEXP,SIGINV,DETSIG,NS,NVAR)
    
    for (N in 1:NS){
      PIVAL[J,N] <- FZMAT[N]*TEMP[N]
    }
  }
  return(PIVAL=PIVAL)
}  
#--------------------------------------------------------------------
#FUNCTION FOR SELECTED ROWS OF PROBABILITY TRANSITION MATRIX      
#--------------------------------------------------------------------
PISUB <- function(JVEC,NRJ,YMAT,NS,NVAR,WVEC,EXVAL,SIGINV,DETSIG,AUTO,B,NLAG){
  I1 <- 1
  I2 <- I1+NS
  I3 <- I2+NS
  I4 <- I3+NVAR*NLAG
  I5 <- I4+NVAR
  
  IST <- DECODE(JVEC,NRJ,NS,NLAG)
  PIOUT <- matrix(0,NRJ,NS)
  X <- matrix(0,1,NVAR*NLAG)
  
  for (I in 1:NRJ){
    for (L in 1:NLAG){
      for (K in 1:NVAR){
        X[(L-1)*NVAR+K] <- YMAT[IST[I,L],K]
      }
    }
    PIOUTT <- PIFN(X,1,YMAT,NS,NVAR,WVEC,EXVAL,SIGINV,DETSIG,AUTO,B,NLAG)  
    VSUM <- 0
    for (J in 1:NS){
      VSUM <- VSUM+PIOUTT[1,J]
    }
    for (M in 1:NS){
      PIOUT[I,M] <- PIOUTT[1,M]/VSUM
    }
  }
  return(PIOUT=PIOUT)
}
#--------------------------------------------------------------------
NRXMAT <- 1
NSTATE <- NS^NLAG
ISEQA <- matrix(0,NSTATE,1)
IREACH <- matrix(0,NSTATE,NS)


for (J in 1:NSTATE){
  IS <- DECODE(J,1,NS,NLAG)
  IR <- REACHR(IS,NS,NLAG)
  
  for (K in 1:NS){
    IREACH[J,K] <- IR[K]
  }
  ISEQA[J] <- J
}

PIOUT <- PISUB(ISEQA,NSTATE,YMAT,NS,NVAR,WVEC,EXVAL,SIGINV,DETSIG,AUTO,B,NLAG)
#--------------------------------------------------------------------
#COMPUTE THE CHAIN'S STATIONARY DISTRIBUTION                        
#--------------------------------------------------------------------
MARKST <- function(PIMAT,NSTATE,NS,NLAG){
  
  #INITIALIZE
  ESPTA <- 0.0001
  MAXIT <- 50
  CRIT <- ESPTA + 1
  NBASE <- NS
  PISTA <- matrix(1/NSTATE,NSTATE,1)
  ISTOLD <- matrix(0,NS,NLAG)
  JOLD <- matrix(0,NS,1)
  TEMP <- matrix(0,NSTATE,1)
  
  #MAIN LOOP TO COMPUTE THE CHAIN'S STATIONARY DISTRIBUTION   
  for (NUMIT in 1:MAXIT){
    
    if (CRIT > ESPTA){
      if (NLAG == 1){
        for (J in 1:NSTATE){
          TEMP[J] <- 0
          for (K in 1:NSTATE){
            TEMP[J] <- TEMP[J] + PISTA[K]*PIMAT[K,J]
          }
        }
      } else {
        for (J in 1:NSTATE){
          IST <- DECODE(J,1,NBASE,NLAG)
          for (K in 1:NS){
            for (L in 1:NLAG-1){
              ISTOLD[K,L] <- IST[L+1]
            }
            ISTOLD[K,NLAG] <- K
          }
          JOLD <- CODE(ISTOLD,NS,NLAG,NBASE)
          SUM <- 0
          for (K in 1:NS){
            SUM <- SUM + PIOUT[JOLD[K],IST[1]]*PISTA[JOLD[K]]
          }
          TEMP[J] <- SUM
        }
      }
      
      CIRT <- 0
      
      for (J in 1:NSTATE){
        DELTA <- abs(TEMP[J]-PISTA[J])
        if (DELTA > CRIT) {CIRT <- DELTA}
        PISTA[J] <- TEMP[J]
      }
    }
    
    if (NUMIT == MAXIT && CIRT > ESPTA){
      stop("COMPLETE!")
    }
  }
  return(PISTA=PISTA)  
}

PISTA <- MARKST(PIOUT,NSTATE,NS,NLAG)
#--------------------------------------------------------------------

return(list(YMAT=YMAT,PIOUT=PIOUT,PISTA=PISTA))

}
