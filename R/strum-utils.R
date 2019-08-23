#------------------------------------------------------------------------------
# Model analysis data Y from StrumData given StrumModel
# - Extract Y
#------------------------------------------------------------------------------
.getAnalysisY = function(myStrumModel, myStrumData, probands)
{
  yNames = myStrumModel@varList$name[myStrumModel@varList$inY==TRUE]

  # Error check if names not there
  if( !length(yNames) )
    stop("No y names exist in StrumModel")

  yAll = tryCatch(dataVals(myStrumData)[,yNames, drop=FALSE],
                  error=function(e)
                  {
                    stop(paste("Some of variables expected from the model ",
                               "are not in the data.  Please check the ",
                               "variable names in your model and data!",
                               sep=""))
                  })

  yAll = split(yAll, dataVals(myStrumData)$family)
  yAll = lapply(yAll, data.matrix)

  yPro = list()
  if( length(probands) > 0 )
    yPro = mapply(.filterMissing, yAll, probands, TRUE, SIMPLIFY=FALSE)

  return(list(yAll = yAll, yPro = yPro))
}

#------------------------------------------------------------------------------
# Model analysis data X from StrumData given StrumModel
# - Extract X from covariates
#------------------------------------------------------------------------------
.getAnalysisX = function(myStrumModel, myStrumData, probands)
{
  xNames = myStrumModel@varList$name[myStrumModel@varList$covariate==TRUE]

  xAll = list()

  # If names not there, then intercept only
  if( !length(xNames) )
  {
    xAll = rep(1, nrow(dataVals(myStrumData)))
    xAll = split(xAll, dataVals(myStrumData)$family)
    xAll = lapply(xAll, data.matrix)
  } else
  {
    xAll = tryCatch(dataVals(myStrumData)[,xNames, drop=FALSE],
                    error=function(e)
                    {
                      stop(paste("Some of covariates expected from the model ",
                                 "are not in the data.  Please check the ",
                                 "covariate names in your model and data!",
                                 sep=""))
                    })

    na_default = options("na.action")$na.action
    options(na.action='na.pass')
    xAll = lapply(split(xAll, dataVals(myStrumData)$family),
                  function(xi)
                  {
                    xi = data.matrix(xi)
                    model.matrix(~xi)
                  })
    options(na.action=na_default)
  }

  xPro = list()
  if( length(probands) > 0 )
    xPro = mapply(.filterMissing, xAll, probands, TRUE, SIMPLIFY=FALSE)

  return(list(xAll = xAll, xPro = xPro))
}

#------------------------------------------------------------------------------
# Model analysis data VC from StrumData given StrumModel
# - Use allRandomEffects to extract VC values
#------------------------------------------------------------------------------
.getAnalysisVC = function(myStrumModel, myStrumData, allRE)
{
  vcPEC = list()

  if( any(allRE == "p") )
    vcPEC = lapply(phi(myStrumData), function(p) return(list(p=p)))

  if( any(allRE == "e") )
  {
    vce = lapply(split(dataVals(myStrumData), dataVals(myStrumData)$family),
                 function(p) return(list(e=diag(nrow(p)))))

    if( length(vcPEC) )
      vcPEC = mapply(function(vci, vcei) return(c(vci, vcei)),
                     vcPEC, vce,
                     SIMPLIFY = FALSE)
    else
      vcPEC = vce
  }

  if( any(allRE == "c") )
  {
    vcc = lapply(phi(myStrumData),
                 function(p) return(list(c=matrix(rep(1,nrow(p)*nrow(p)),nrow=nrow(p)))))

    if( length(vcPEC) )
      vcPEC = mapply(function(vci, vcci) return(c(vci, vcci)),
                     vcPEC, vcc,
                     SIMPLIFY = FALSE)

    else
      vcPEC = vcc
  }

  return(vcPEC)
}

.getValidAnalysisData = function(y, x, vc)
{
  missingX = lapply(x$xAll, .findMissing, isX=TRUE)
  missingY = lapply(y$yAll, function(yk) return(rowSums(is.na(yk))!= ncol(yk)))

  missingXY = mapply("&", missingX, missingY, SIMPLIFY=FALSE)

  xa  = mapply(.filterMissing,   x$xAll,   missingXY, TRUE,      SIMPLIFY=FALSE)
  ya  = mapply(.filterMissing,   y$yAll,   missingXY, TRUE,      SIMPLIFY=FALSE)
  vca = mapply(.filterMissingVC, vc$vcAll, missingXY, missingXY, SIMPLIFY=FALSE)

  xan = Filter(function(x) nrow(x) > 0, xa)
  ped_names = names(xan)
  yan  = ya[ped_names]
  vcan = vca[ped_names]

  ypn = list()
  xpn = list()
  vcpn = list()

  if( length(y$yPro) > 0 )
  {
    missingXp = lapply(x$xPro, .findMissing, isX=TRUE)
    missingYp = lapply(y$yPro, function(yk) return(rowSums(is.na(yk))!= ncol(yk)))

    missingXYp = mapply("&", missingXp, missingYp, SIMPLIFY=FALSE)

    xp  = mapply(.filterMissing,   x$xPro,   missingXYp, TRUE,       SIMPLIFY=FALSE)
    yp  = mapply(.filterMissing,   y$yPro,   missingXYp, TRUE,       SIMPLIFY=FALSE)
    vcp = mapply(.filterMissingVC, vc$vcPro, missingXYp, missingXYp, SIMPLIFY=FALSE)

    ypn  = yp[ped_names]
    xpn  = xp[ped_names]
    vcpn = vcp[ped_names]
  }


  return(list(y  = list(yAll=yan,yPro=ypn),
              x  = list(xAll=xan,xPro=xpn),
              vc = list(vcAll=vcan,vcPro=vcpn)))
}

#------------------------------------------------------------------------------
# Find the missing data
#------------------------------------------------------------------------------
.findMissing = function(yk, isX)
{
  if( isX )
    return(complete.cases(yk))

  return(!is.na(yk))
}

#------------------------------------------------------------------------------
# Filter the missing entries in x and y.
#------------------------------------------------------------------------------
.filterMissing = function(yk, mk, isX)
{
  if( !isX )
    return(yk[mk, drop = FALSE])

  return(yk[mk,, drop = FALSE])
}

#------------------------------------------------------------------------------
# Filter the missing entries in vc matrix.
#------------------------------------------------------------------------------
.filterMissingVC = function(vck, mk1, mk2)
{
  return(lapply(vck, function(z) return(z[mk1,mk2, drop = FALSE])))
}

#------------------------------------------------------------------------------
# Calculate the model parameters in stage 1.
# - Estimate delta parameters, C and V in Equation 5 and 6.
#------------------------------------------------------------------------------
.estimateDeltaParameter = function(y, x, vc, step1OptimControl)
{
  numTrait = ncol(y$yAll[[1]])
  numCov   = ncol(x$xAll[[1]])
  numVC    = length(vc$vcAll[[1]])
  numPed   = length(y$yAll)

  betaMat  = matrix(0, nrow = numCov, ncol = numTrait, byrow = FALSE)
  sigmaMat = list()
  for( j in 1:numVC )
    sigmaMat[[j]] = matrix(0, numTrait, numTrait)

  ypt1 = list()
  ypt2 = list()
  xpt1 = list()
  xpt2 = list()
  vcpt1 = list()
  vcpt2 = list()
  vcp12 = list()

  usedN = rep(0, numTrait)
  for( t1 in 1:numTrait )
  {
    #
    # 1. Compute univariate log-likelihood
    #--------------------------------------
    # all
    #-----
    t1dat = .getYFilteredData(y$yAll, x$xAll, vc$vcAll, t1)
    yt1  = t1dat$y
    xt1  = t1dat$x
    vct1 = t1dat$vc
    usedN[t1] = usedN[t1] + length(unlist(yt1))

    if( length(y$yPro) )
    {
      # probands
      #----------
      pt1dat = .getYFilteredData(y$yPro, x$xPro, vc$vcPro, t1)
      ypt1  = pt1dat$y
      xpt1  = pt1dat$x
      vcpt1 = pt1dat$vc
    }

    temp1 = .computeUnivariateLL(yt1, xt1, vct1, ypt1, xpt1, vcpt1, step1OptimControl)

    betaMat[,t1] = temp1$beta
    for( v in 1:numVC )
      sigmaMat[[v]][t1,t1] = temp1$sigma[v]

    if( t1 == 1 )
      next

    #
    # 2. Compute bivariate log-likelihood
    #-------------------------------------
    for( t2 in (t1-1):1 )
    {
      # all
      #-----
      t2dat = .getYFilteredData(y$yAll, x$xAll, vc$vcAll, t2)
      yt2  = t2dat$y
      xt2  = t2dat$x
      vct2 = t2dat$vc
      vc12 = .getVCFilteredData(vc$vcAll, t1dat$missing, t2dat$missing)

      if( length(y$yPro) )
      {
        # probands
        #----------
        pt2dat = .getYFilteredData(y$yPro, x$xPro, vc$vcPro, t2)
        ypt2  = pt2dat$y
        xpt2  = pt2dat$x
        vcpt2 = pt2dat$vc
        vcp12 = .getVCFilteredData(vc$vcPro, pt1dat$missing, pt2dat$missing)
      }

      bet12 = matrix(betaMat[,c(t1,t2)], ncol=2)
      sig12 = lapply(sigmaMat, function(s) return(diag(s)[c(t1,t2)]))

      temp2 = .computeBivariateLL(yt1,  xt1,  vct1,  yt2,  xt2,  vct2,  vc12,
                                  ypt1, xpt1, vcpt1, ypt2, xpt2, vcpt2, vcp12,
                                  bet12, sig12, step1OptimControl)

      for( v in 1:numVC )
      {
        sigmaMat[[v]][t1,t2] = temp2$sigma[v]
        sigmaMat[[v]][t2,t1] = temp2$sigma[v]
      }
    }
  }

  return(list(C=betaMat, V=sigmaMat, usedN=usedN))
}

#------------------------------------------------------------------------------
# Get filtered list of data Y, X VC
#------------------------------------------------------------------------------
.getYFilteredData = function(y, x, vc, t)
{
  tmp_y = lapply(y, function(yk) return(yk[,t]))

  missingY = lapply(tmp_y, .findMissing, isX=FALSE)

  yall  = mapply(.filterMissing,   tmp_y,  missingY, FALSE,    SIMPLIFY=FALSE)
  xall  = mapply(.filterMissing,   x,      missingY, TRUE,     SIMPLIFY=FALSE)
  vcall = mapply(.filterMissingVC, vc,     missingY, missingY, SIMPLIFY=FALSE)

  return(list(y=yall,x=xall,vc=vcall, missing=missingY))
}

#------------------------------------------------------------------------------
# A function which utilizes the C likelihood function to estimate
#  the log-likelihood of "saturated model" - univariate .
#------------------------------------------------------------------------------
.computeUnivariateLL = function(yt1, xt1, vct1, ypt1, xpt1, vcpt1, step1OptimControl)
{
  numCov = ncol(xt1[[1]])
  numVC  = length(vct1[[1]])

  yb = do.call("c",yt1)
  xb = do.call(rbind,xt1)

  fit = lm(yb ~ xb-1)
  betaStart = coef(fit)
  sigmaTot = summary(fit)$sigma

  startVal = c(rep(sigmaTot,numVC)/numVC, betaStart)

  logLikelihoodP = 0.0;

  Q = function(theta)
  {
    sigma = .makeSigmaList(theta[1:numVC], list())

    pd = sapply(sigma, .testPosDefMatrix)
    if( any(pd == FALSE) )
      return(-1*10^6)

    beta = matrix(theta[(numVC+1):(numVC + numCov)], ncol = 1)

    logLikelihood = .Call("computeLL",
                          yt1, xt1, vct1, list(), list(), list(), list(), beta, sigma)

    if( length(ypt1) )
      logLikelihoodP = .Call("computeLL",
                             ypt1, xpt1, vcpt1, list(), list(), list(), list(), beta, sigma)

    return(logLikelihood - logLikelihoodP)
  }

  optimOut = optim(startVal, Q, control=step1OptimControl)

  if( optimOut$convergence != 0 )
    warning("optim did not converge in stage 1.")

  ret = list(sigma = optimOut$par[1:numVC],
             beta  = optimOut$par[(numVC+1):(numVC + numCov)])

  return(ret)
}

#------------------------------------------------------------------------------
# Form a list of variacne-covariance matrices.
#------------------------------------------------------------------------------
.makeSigmaList = function(vari, cor)
{
  numVC = length(vari)
  sigmaL = list()

  if( length(cor) == 0 )
    for( v in 1:numVC )
      sigmaL[[v]] = matrix(c(vari[v]), nrow = 1, ncol = 1, byrow = FALSE)
  else
    for( v in 1:numVC )
    {
      covari = sqrt(vari[[v]][1]*vari[[v]][2]) * cor[v]

      sigmaL[[v]] = matrix(c(vari[[v]][1],covari,covari,vari[[v]][2]),
                           nrow = 2, ncol = 2, byrow = FALSE)
    }

  return(sigmaL)
}

#------------------------------------------------------------------------------
# Check whether a matrix is positive definite or not.
#------------------------------------------------------------------------------
.testPosDefMatrix = function(mat)
{
  if( nrow(mat) == 1 )
  {
    if( mat[1,1] > 0 )
      return(TRUE)
    else
      return(FALSE)
  }

  if( (mat[1,1]>0) & (mat[2,2]>0) & ((mat[1,1]*mat[2,2]) > mat[2,1]*mat[1,2]) )
    return(TRUE)
  else
    return(FALSE)
}

#------------------------------------------------------------------------------
# Get filtered list of data COV
#------------------------------------------------------------------------------
.getVCFilteredData = function(vc, missingY1, missingY2)
{
  return(mapply(.filterMissingVC, vc, missingY1, missingY2, SIMPLIFY=FALSE))
}

#------------------------------------------------------------------------------
# A function which utilizes the C likelihood function to estimate
#  the log-likelihood of "saturated model" - bivariate .
#------------------------------------------------------------------------------
.computeBivariateLL = function(yt1,  xt1,  vct1,  yt2,  xt2,  vct2,  vc12,
                               ypt1, xpt1, vcpt1, ypt2, xpt2, vcpt2, vcp12,
                               beta, sig, step1OptimControl)
{
  numVC = length(vct1[[1]])

  startVal1 = rep(0, numVC)
  startVal2 = .rdirichlet(1, rep(1, numVC))

  logLikelihoodP = 0.0;

  Q = function(theta)
  {
    sigma = .makeSigmaList(sig, theta)

    pd = sapply(sigma, .testPosDefMatrix)
    if( any(pd == FALSE) )
      return(-1*10^6)

    logLikelihood = .Call("computeLL",
                          yt1, xt1, vct1, yt2, xt2, vct2, vc12, beta, sigma)

    if( length(ypt1) )
      logLikelihoodP = .Call("computeLL",
                             ypt1, xpt1, vcpt1, ypt2, xpt2, vcpt2, vcp12, beta, sigma)

    return(logLikelihood - logLikelihoodP)
  }

  optimOut1 = optim(startVal1, Q, control=step1OptimControl)

  if( optimOut1$convergence != 0 )
    warning("optim did not converge in stage 1.")

  optimOut2 = optim(startVal2, Q, control=step1OptimControl)

  if( optimOut2$convergence != 0 )
    warning("optim did not converge in stage 1.")

  logLLs   = c(optimOut1$value, optimOut2$value)
  maxLLPos = which.max(logLLs)

  optimPars = rbind(optimOut1$par, optimOut2$par)

  sigVal0 = sqrt(sig[[1]][1] * sig[[1]][2]) * optimPars[maxLLPos,][1]
  if( numVC > 1 )
    for( v in 2:numVC )
      sigVal0 = c(sigVal0, sqrt(sig[[v]][1] * sig[[v]][2]) * optimPars[maxLLPos,][v])

  ret = list(sigma = sigVal0)

  return(ret)
}

#------------------------------------------------------------------------------
# A function for sampling from a Dirichlet distribution
# - improvement over Null {MASS}
#------------------------------------------------------------------------------
.rdirichlet = function(n, alpha)
{
  k = length(alpha)
  r = matrix(0, nrow=n, ncol=k)
  for( i in 1:k )
  {
    r[,i] = rgamma(n, alpha[i], 1)
  }
  r = matrix(mapply(function(r, s) {return (r/s)}, r, rowSums(r)), ncol=k)

  return (r)
}

#------------------------------------------------------------------------------
# Estimate the 1st and 2nd derivatives of C and V
#------------------------------------------------------------------------------
.estimateDeltaDerivative = function(y, x, vc, betaMat, sigmaMat)
{
  numTrt = ncol(y$yAll[[1]])
  numCov = ncol(x$xAll[[1]])
  numVC  = length(vc$vcAll[[1]])
  numPed = length(y$yAll)

  fDeriv = matrix(ncol = numPed, nrow = 0)
  sDeriv = list()
  s1Deriv = list()
  s2Deriv = list()

  for( t1 in 1:numTrt )
  {
    # 1. univariate derivatives
    #---------------------------

    bet1 = matrix(betaMat[,t1], ncol = 1)
    sig1 = lapply(sigmaMat, function(s) return(matrix(s[t1,t1])))

    # all
    #-----
    t1dat = .getYFilteredData(y$yAll, x$xAll, vc$vcAll, t1)
    yt1  = t1dat$y
    xt1  = t1dat$x
    vct1 = t1dat$vc

    dUni = .Call("computeDeriv",
                 yt1, xt1, vct1, list(), list(), list(), list(), bet1, sig1)

    if( length(y$yPro) )
    {
      # probands
      #----------
      pt1dat = .getYFilteredData(y$yPro, x$xPro, vc$vcPro, t1)
      ypt1  = pt1dat$y
      xpt1  = pt1dat$x
      vcpt1 = pt1dat$vc

      dUniP = .Call("computeDeriv",
                    ypt1, xpt1, vcpt1, list(), list(), list(), list(), bet1, sig1)

      fDeriv = rbind(fDeriv, dUni$dPar - dUniP$dPar)
      sDeriv[[length(sDeriv)+1]] = dUni$d2Par - dUniP$d2Par
    } else
    {
      fDeriv = rbind(fDeriv, dUni$dPar)
      sDeriv[[length(sDeriv)+1]] = dUni$d2Par
    }

    if( t1 == 1 )
      next

    # 2. bivariate derivatives
    #--------------------------
    for( t2 in (t1-1):1 )
    {
      bet12 = matrix(betaMat[,c(t1,t2)], ncol=2)
      sig12 = lapply(sigmaMat, function(s) return(s[c(t1,t2),c(t1,t2)]))

      # all
      #-----
      t2dat = .getYFilteredData(y$yAll, x$xAll, vc$vcAll, t2)
      yt2  = t2dat$y
      xt2  = t2dat$x
      vct2 = t2dat$vc
      vc12 = .getVCFilteredData(vc$vcAll, t1dat$missing, t2dat$missing)

      dBi = .Call("computeDeriv",
                  yt1, xt1, vct1, yt2, xt2, vct2, vc12, bet12, sig12)

      if( length(y$yPro) )
      {
        # probands
        #----------
        pt2dat = .getYFilteredData(y$yPro, x$xPro, vc$vcPro, t2)
        ypt2  = pt2dat$y
        xpt2  = pt2dat$x
        vcpt2 = pt2dat$vc
        vcp12 = .getVCFilteredData(vc$vcPro, pt1dat$missing, pt2dat$missing)

        dBiP = .Call("computeDeriv",
                     ypt1, xpt1, vcpt1, ypt2, xpt2, vcpt2, vcp12, bet12, sig12)

        fDeriv = rbind(fDeriv, dBi$dPar - dBiP$dPar)
        sDeriv[[length(sDeriv)+1]]   = dBi$d2Par - dBiP$d2Par
        s1Deriv[[length(s1Deriv)+1]] = dBi$dP1dS - dBiP$dP1dS
        s2Deriv[[length(s2Deriv)+1]] = dBi$dP2dS - dBiP$dP2dS
      } else
      {
        fDeriv = rbind(fDeriv, dBi$dPar)
        sDeriv[[length(sDeriv)+1]]   = dBi$d2Par
        s1Deriv[[length(s1Deriv)+1]] = dBi$dP1dS
        s2Deriv[[length(s2Deriv)+1]] = dBi$dP2dS
      }
    }
  }

  return(list(numTrt=numTrt, numCov=numCov, numVC=numVC,
              fDeriv=fDeriv, sDeriv=sDeriv, s1Deriv=s1Deriv, s2Deriv=s2Deriv))
}

#------------------------------------------------------------------------------
# Estimate Covariance(delta) = (J'FJ'), Equation 8.
#------------------------------------------------------------------------------
.estimateDeltaCovariance = function(deltaDeriv)
{
  #
  # 1. Estimate J and F
  #---------------------
  K = ncol(deltaDeriv$fDeriv) # equal to numPed

  # Use .bdiag to stack them on a diagonal
  # - a sparse matrix is constructed.
  # Needs library(Matrix)
  J = .bdiag(deltaDeriv$sDeriv)

  if( deltaDeriv$numTrt > 1 )
  {
    index = 1
    rowOffset = deltaDeriv$numCov + deltaDeriv$numVC

    for( t1 in 2:deltaDeriv$numTrt )
    {
      r1e = rowOffset * t1 + .getCovBlockCount(t1-2) * deltaDeriv$numVC
      r1s = r1e - rowOffset + 1

      for( t2 in (t1-1):1 )
      {
        r2e = rowOffset * t2 + .getCovBlockCount(t2-2) * deltaDeriv$numVC
        r2s = r2e - rowOffset + 1

        cs = r1e + 1 + (t1-t2-1) * deltaDeriv$numVC
        ce = cs + deltaDeriv$numVC - 1

        # upper matrix
        #J[r1s:r1e, cs:ce] = deltaDeriv$s1Deriv[[index]]
        #J[r2s:r2e, cs:ce] = deltaDeriv$s2Deriv[[index]]

        # lower matrix
        J[cs:ce, r1s:r1e] = t(deltaDeriv$s1Deriv[[index]])
        J[cs:ce, r2s:r2e] = t(deltaDeriv$s2Deriv[[index]])

        index = index + 1
      }
    }
  }

  J = J/K

  F = (deltaDeriv$fDeriv %*% t(deltaDeriv$fDeriv))/K

  #
  # 2. Calculate covariance matrix
  #--------------------------------
  #Ji = solve(J)
  Ji = tryCatch(solve(J),
                error=function(e)
                {
                  stop(paste("Covariance matrix of the model parameters ",
                             "is not invertable in step 1!  The model may ",
                             "not be identifiable, or the sample size ",
                             "is too small.",
                             sep=""))
                })

  JiFJit = Ji %*% F %*% t(Ji)

  w = 1/diag(JiFJit)

  return(list(K=K, deltaCov=JiFJit, w=w))
}

#------------------------------------------------------------------------------
# Estimate asymptotic covariance(theta) = (A'BA')/k of the model.
#------------------------------------------------------------------------------
.estimateThetaCovariance = function(model, s1, theta)
{
  dFdTheta = .funDeriv(fun=function(theta) .thetaToDelta(model,theta), theta)
  #dfOrthog = Null(t(dFdTheta)) #Library MASS
  dfOrthog = .strumNull(t(dFdTheta), tol=1e-07)

  # Calculate covariance matrix of theta.
  #---------------------------------------
  W = diag(s1$w)
  WJFJW = W %*% s1$deltaCov %*% W

  A = dFdTheta %*% W     %*% t(dFdTheta)
  B = dFdTheta %*% WJFJW %*% t(dFdTheta)

  #Ai = solve(A)
  Ai = tryCatch(solve(A), error=function(e) return(ginv(A)))

  thetaCov = (Ai %*% B %*% Ai) / s1$K

  # Calculate covariance matrix of diff=(delta-theta).
  #----------------------------------------------------
  diff     = s1$delta - .thetaToDelta(model, theta)
  dFAidFW  = t(dFdTheta) %*% Ai %*% dFdTheta %*% W
  IdFAidFW = diag(nrow(W)) - dFAidFW

  diffCov = (IdFAidFW %*% s1$deltaCov %*% t(IdFAidFW))

  U = W - (W %*% dFAidFW)

  return(list(thetaCov = thetaCov,
              diff     = diff,
              diffCov  = diffCov,
              dfOrthog = dfOrthog,
              U        = U))
}

#------------------------------------------------------------------------------
# Find the number of previous covariance blocks.
#------------------------------------------------------------------------------
.getCovBlockCount = function(N)
{
  if( N <= 0 )
    return(0)
  else
    return(sum(1:N))
}

#------------------------------------------------------------------------------
# Fitting the Model, stage 1
# - Form a limited information estimate for the saturated model.
#------------------------------------------------------------------------------
.fitModelStep1 = function(y, x, vc, step1OptimControl)
{
  CV = .estimateDeltaParameter(y, x, vc, step1OptimControl)
  De = .estimateDeltaDerivative(y, x, vc, CV$C, CV$V)
  JF = .estimateDeltaCovariance(De)

  delta = .Call("getDelta", t(CV$C), CV$V)

  .printInfoLine("Fitting model step 1", "Done", 50, 2)

  return(list(CV=CV, K=JF$K, delta=delta, deltaCov=JF$deltaCov, w=JF$w, Wi=De$W))
}

#------------------------------------------------------------------------------
# Fitting the Model, stage 2
# - Estimate parameters and asymptotic covariance matrix of the model.
#------------------------------------------------------------------------------
.fitModelStep2 = function(s1, model, startValueControl, step2OptimControl)
{
  startTheta = .generateStartValue(model, s1$delta, s1$w, startValueControl)
  optimTheta = .estimateThetaParameter(model, s1$delta, s1$w, startTheta, step2OptimControl)
  parCov     = .estimateThetaCovariance(model, s1, optimTheta)

  .printInfoLine("Fitting model step 2", "Done", 50, 2)

  return(list(theta    = optimTheta,
              thetaCov = parCov$thetaCov,
              diff     = parCov$diff,
              diffCov  = parCov$diffCov,
              dfOrthog = parCov$dfOrthog,
              U        = parCov$U))
}

#------------------------------------------------------------------------------
# Print analysis progress information,
#  a property name and value with padding dots between them.
#------------------------------------------------------------------------------
.printInfoLine = function(pInfo, value, width, leadSize=8)
{
  value = as.character(value)
  if( length(value) == 0 )
    value = "NA"

  padSize   = width - nchar(value) - nchar(pInfo)
  dotPad    = paste(rep(".", padSize), collapse="")
  leadSpace = paste(rep(" ", leadSize), collapse="")

  cat(leadSpace, pInfo, dotPad, value, "\n")
}

#------------------------------------------------------------------------------
# Generate the starting values, step 2
#------------------------------------------------------------------------------
.generateStartValue = function(model, delta, w, startValueControl)
{
  Q = function(theta)
  {
    diff = delta - .thetaToDelta(model, theta)
    return(sum(diff*diff*w))
  }

  positive = grep("<", paramNames(model))

  initialPopulation = length(paramNames(model))*120
  if( !is.null(startValueControl[["initPopulation"]]) )
    initialPopulation = startValueControl[["initPopulation"]]

  startVals = lapply(1:initialPopulation,
                     function(dummy)
                     {
                       startVal = rnorm(length(paramNames(model)))
                       startVal[positive] = abs(startVal[positive])
                       return(list(startVal=startVal, Q=Q(startVal)))
                     })

  nChildren         = length(paramNames(model))*2
  nGenerations      = 20
  selection1        = 15
  selection2        = 15

  if( !is.null(startValueControl[["nChildren"]]) )
    nChildren = startValueControl[["nChildren"]]

  if( !is.null(startValueControl[["nGenerations"]]) )
    nGenerations = startValueControl[["nGenerations"]]

  if( !is.null(startValueControl[["selection1"]]) )
    selection1 = startValueControl[["selection1"]]

  if( !is.null(startValueControl[["selection2"]]) )
    selection2 = startValueControl[["selection2"]]

  QS = sapply(startVals, function(x) return(x$Q))
  startVals = startVals[order(QS)[1:selection1]]

  Delta = 0 != .funDeriv(fun = function(theta) .thetaToDelta(model,theta),
                         startVals[[1]]$startVal )

  for( gen in 1:nGenerations )
  {
    newGen = list()
    keepTrack = rep(0,length(startVals))

    for( i in 1:length(startVals) )
    {
      deltaTemp = sign(.thetaToDelta(model, startVals[[i]]$startVal)) != sign(delta)
      keepTrack[i] = length(deltaTemp)-sum(deltaTemp)
      deltaTemp = (Delta%*%deltaTemp)
      gammai = which(deltaTemp>=1)
      gammai = gammai[!(gammai %in% positive)]
      stVali = startVals[[i]]$startVal
      stValiGammai = stVali[gammai]

      for( child in 1:nChildren )
      {
        newgenTmp = stVali
        newgenTmp[gammai] = stValiGammai * (rbinom(length(gammai), 1, 0.5)-.5)*2
        newgenTmp = newgenTmp+rnorm(length(stVali), sd=.01)
        deltaTemp = sign(.thetaToDelta(model, newgenTmp)) != sign(delta)
        newGen[[child]] = list(startVal=newgenTmp,
                               Q=Q(newgenTmp)+ sum(deltaTemp)*1000)
      }

      QS = sapply(newGen, function(x) return(x$Q))
      startVals[i] = newGen[order(QS)[1]]
    }

    #cat(gen, keepTrack, "\n")
    tbl = as.matrix(table(keepTrack))
    tmp = as.integer(rownames(tbl))
    maxrow = which(tmp==max(tmp))

    if( tbl[maxrow]/sum(tbl) > .85 )
      break
  }

  QS = sapply(startVals, function(x) return(x$Q))
  startVals = startVals[order(QS)[1:selection2]]

  return(startVals)
}

#------------------------------------------------------------------------------
# Calculate delta values (C and V as a vector) given theta values.
#  - equation 5 and 6
#  - use library(MASS)
#------------------------------------------------------------------------------
.thetaToDelta = function(model, theta)
{
  B  = model@B(theta)
  L  = model@L(theta)
  Gm = model@Gm(theta)
  Gs = model@Gs(theta)
  Z  = lapply(model@Z, function(v) v(theta))
  E  = lapply(model@E, function(v) v(theta))

  return(.Call("computeDelta", B, L, Gs, Gm, Z, E))
}

#------------------------------------------------------------------------------
# A function which calculates the vector of partial derivatives of a function
#  numerically
#------------------------------------------------------------------------------
.funDeriv = function(fun, theta, step = 1E-6)
{
  f0 = fun(theta)
  len_fun = length(f0)
  d = matrix(0, length(theta), len_fun)

  for( i in 1:length(theta) )
  {
    h     = rep(0, length(theta))
    h[i]  = step
    d[i,] = (fun(theta + h) - f0)/step
  }

  return(d)
}

#------------------------------------------------------------------------------
# Calculate the model parameters in stage 2.
# - Estimate reduced model parameters (theta), Equation 10.
#------------------------------------------------------------------------------
.estimateThetaParameter = function(model, delta, w, startTheta, step2OptimControl)
{
  positive = grep("<", paramNames(model))

  lengthNV = length(paramNames(model)) - length(positive)

  Qstar = function(thetaNV)
  {
    mMat = .getModelMats(model, thetaNV)

    B  = mMat$B
    L  = mMat$L
    Gs = mMat$Gs
    Gm = mMat$Gm

    Zs = mMat$Zs
    Es = mMat$Es

    q = .Call("computeQstar", B, L, Gs, Gm, Zs, Es, as.matrix(delta), as.matrix(w))

    return(q)
  }

  optimOut = list()
  mymin = 10^100

  for( i in 1:length(startTheta) )
  {
    optimOut[[i]] = tryCatch(
      optim(startTheta[[i]]$startVal[1:lengthNV],
            Qstar,
            control = step2OptimControl,
            method = "BFGS"),
      error = function(d)
      {
        sval = startTheta[[i]]$startVal[1:lengthNV] + rnorm(lengthNV, sd=.1)
        return(tryCatch(
          optim(sval,
                Qstar,
                control = step2OptimControl,
                method = "BFGS"),
          error = function(d) list(value=10^100,par=NA)))
      })

    if( abs((optimOut[[i]]$value-mymin)/mymin) < 10^-4 )
      break

    if( optimOut[[i]]$value < mymin )
      mymin = optimOut[[i]]$value
  }

  tmp = sapply(optimOut, function(ooi) ooi$value)
  minOut = c(which(tmp==min(tmp)))
  optimOut = optimOut[[minOut[1]]]

  if( optimOut$convergence != 0 )
    warning("optim did not converge in stage 2.")

  mMat = .getModelMats(model, optimOut$par)

  B  = mMat$B
  L  = mMat$L
  Gs = mMat$Gs
  Gm = mMat$Gm

  Zs = mMat$Zs
  Es = mMat$Es

  thetaV = .Call("computeThetaV", B, L, Gs, Gm, Zs, Es, as.matrix(delta), as.matrix(w))

  optimPar = c(optimOut$par, thetaV)
  names(optimPar) = paramNames(model)

  return(optimPar)
}

#------------------------------------------------------------------------------
# Estimate asymptotic covariance(theta) = (A'BA')/k of the model.
#------------------------------------------------------------------------------
.estimateThetaCovariance = function(model, s1, theta)
{
  dFdTheta = .funDeriv(fun=function(theta) .thetaToDelta(model,theta), theta)
  #dfOrthog = Null(t(dFdTheta)) #Library MASS
  dfOrthog = .strumNull(t(dFdTheta), tol=1e-07)

  # Calculate covariance matrix of theta.
  #---------------------------------------
  W = diag(s1$w)
  WJFJW = W %*% s1$deltaCov %*% W

  A = dFdTheta %*% W     %*% t(dFdTheta)
  B = dFdTheta %*% WJFJW %*% t(dFdTheta)

  #Ai = solve(A)
  Ai = tryCatch(solve(A), error=function(e) return(ginv(A)))

  thetaCov = (Ai %*% B %*% Ai) / s1$K

  # Calculate covariance matrix of diff=(delta-theta).
  #----------------------------------------------------
  diff     = s1$delta - .thetaToDelta(model, theta)
  dFAidFW  = t(dFdTheta) %*% Ai %*% dFdTheta %*% W
  IdFAidFW = diag(nrow(W)) - dFAidFW

  diffCov = (IdFAidFW %*% s1$deltaCov %*% t(IdFAidFW))

  U = W - (W %*% dFAidFW)

  return(list(thetaCov = thetaCov,
              diff     = diff,
              diffCov  = diffCov,
              dfOrthog = dfOrthog,
              U        = U))
}
