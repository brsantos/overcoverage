# Auxiliary function to fit model
#' @importFrom tseries read.matrix
bict2 <- function (formula, data = NULL, n.sample, prior = "SBH", cens = NULL,
          start.formula = NULL, start.beta = NULL, start.sig = NULL,
          start.y0 = NULL, save = 0, name = NULL, null.move.prob = 0.5,
          a = 0.001, b = 0.001, progress = FALSE)
{
  if (n.sample <= 0) {
    stop("n.sample must be positive")
  }
  if (prior != "UIP" & prior != "SBH") {
    stop("prior not found")
  }
  if (save < 0) {
    stop("save must be non-negative")
  }
  if (null.move.prob < 0 | null.move.prob > 1) {
    stop("null.move.prob is a probability and should be between 0 and 1")
  }
  if (a < 0 & a != (-1)) {
    stop("a and b must be non-negative")
  }
  if (b < 0) {
    stop("a and b must be non-negative")
  }
  ptm <- (proc.time())[3]
  # if (!is.null(data)) {
  #   if (attributes(data)$class == "table") {
  #     data <- data.frame(data)
  #   }
  # }
  if (save > 0) {
    if (is.null(name)) {
      name_RJACC <- "RJACC.txt"
      name_MHACC <- "MHACC.txt"
      name_BETA <- "BETA.txt"
      name_MODEL <- "MODEL.txt"
      name_SIG <- "SIG.txt"
      name_Y0 <- "Y0.txt"
    }
    else {
      name_RJACC <- paste(name, "RJACC.txt", sep = "")
      name_MHACC <- paste(name, "MHACC.txt", sep = "")
      name_BETA <- paste(name, "BETA.txt", sep = "")
      name_MODEL <- paste(name, "MODEL.txt", sep = "")
      name_SIG <- paste(name, "SIG.txt", sep = "")
      name_Y0 <- paste(name, "Y0.txt", sep = "")
    }
    if (file.exists(name_BETA)) {
      stop(paste("A file named ", name_BETA, " already exists in the working directory",
                 sep = ""))
    }
    if (file.exists(name_MODEL)) {
      stop(paste("A file named ", name_MODEL, " already exists in the working directory",
                 sep = ""))
    }
    if (file.exists(name_SIG)) {
      stop(paste("A file named ", name_SIG, " already exists in the working directory",
                 sep = ""))
    }
    if (file.exists(name_Y0)) {
      stop(paste("A file named ", name_Y0, " already exists in the working directory",
                 sep = ""))
    }
    if (file.exists(name_RJACC)) {
      stop(paste("A file named ", name_RJACC, " already exists in the working directory",
                 sep = ""))
    }
    if (file.exists(name_MHACC)) {
      stop(paste("A file named ", name_MHACC, " already exists in the working directory",
                 sep = ""))
    }
  }
  priortypes <- c("UIP", "SBH")
  priornum <- c(1, 2)[prior == priortypes]
  options(contrasts = c("contr.sum", "contr.poly"), warn = -1)
  if (!is.null(data)) {
    maximal.mod <- stats::glm(formula = formula, data = data, method = "model.frame",
                       na.action = na.pass, family = poisson, control = list(maxit = 1),
                       x = TRUE, y = TRUE)
  }
  else {
    maximal.mod <- stats::glm(formula = formula, method = "model.frame",
                       na.action = na.pass, family = poisson, control = list(maxit = 1),
                       x = TRUE, y = TRUE)
  }
  options(contrasts = c("contr.treatment", "contr.poly"),
          warn = 0)
  missing1 <- (1:length(maximal.mod[, 1]))[is.na(maximal.mod[,
                                                             1])]
  data <- maximal.mod
  data[missing1, 1] <- 0
  if (!is.null(cens)) {
    missing2 <- cens
  }
  else {
    missing2 <- c()
  }
  missing <- c(missing1, missing2)
  missing_details <- data[missing1, -1]
  censored_details <- data[missing2, -1]
  options(contrasts = c("contr.sum", "contr.poly"), warn = -1)
  if (!is.null(data)) {
    maximal.mod <- stats::glm(formula = formula, data = data, family = poisson,
                       control = list(maxit = 1), x = TRUE, y = TRUE)
  }
  else {
    maximal.mod <- stats::glm(formula = formula, family = poisson,
                       control = list(maxit = 1), x = TRUE, y = TRUE)
  }
  options(contrasts = c("contr.treatment", "contr.poly"), warn = 0)
  big.X <- maximal.mod$x
  y <- maximal.mod$y
  n <- dim(big.X)[1]
  IP <- t(big.X) %*% big.X/n
  IP[, 1] <- 0
  IP[1, 0] <- 0
  if (is.null(start.beta)){
    bmod <- conting::beta_mode(X = big.X[-missing, ],
                               y = y[-missing],
                               prior = prior, IP = IP, a = a, b = b)
    eta.hat <- as.vector(big.X %*% matrix(bmod, ncol = 1))
  }
  else {
    eta.hat <- as.vector(big.X %*% matrix(start.beta, ncol = 1))
  }

  if (is.null(start.formula)) {
    start.index <- rep(1, dim(big.X)[2])
  }
  else {
    start.index <- conting::formula2index(big.X = big.X,
                                          formula = start.formula,
                                          data = data)
  }
  if (is.null(start.beta)) {
    start.beta <- bmod[start.index == 1]
  }
  if (is.null(start.sig)) {
    start.sig <- 1
  }
  if (is.null(start.y0)) {
    start.y0 <- round(exp(eta.hat[missing]))
  }
  start.mod <- conting::index2model(start.index)
  runit <- bict.fit2(priornum = priornum, missing1 = missing1,
                    missing2 = missing2,
                    maximal.mod = maximal.mod, IP = IP,
                    eta.hat = eta.hat,
                    ini.index = start.index,
                    ini.beta = start.beta,
                    ini.sig = start.sig,
                    ini.y0 = start.y0, iters = n.sample,
                    save = save, name = name,
                    null.move.prob = null.move.prob,
                    a = a, b = b, progress = progress)
  BETA <- runit$BETA
  MODEL <- runit$MODEL
  SIG <- runit$SIG
  Y0 <- runit$Y0
  rj_acc <- runit$rj_acc
  mh_acc <- runit$mh_acc
  if (save > 0) {
    rj_acc <- tseries::read.matrix(file = name_RJACC, header = FALSE)
    mh_acc <- tseries::read.matrix(file = name_MHACC, header = FALSE)
    BETA <- tseries::read.matrix(file = name_BETA, header = FALSE)
    SIG <- tseries::read.matrix(file = name_SIG, header = FALSE)
    Y0 <- matrix(tseries::read.matrix(file = name_Y0, header = FALSE),
                 nrow = length(SIG))
    MODEL <- as.character(utils::read.table(file = name_MODEL,
                                            header = FALSE)[, 1])
  }
  time <- (proc.time())[3] - ptm
  est <- list(BETA = BETA, MODEL = MODEL, SIG = SIG, Y0 = Y0,
              missing1 = missing1, missing2 = missing2, missing_details = missing_details,
              censored_details = censored_details, rj_acc = rj_acc,
              mh_acc = mh_acc, priornum = priornum, maximal.mod = maximal.mod,
              IP = IP, eta.hat = eta.hat, save = save, name = name,
              null.move.prob = null.move.prob, time = time, a = a,
              b = b)
  class(est) <- "bict"
  est
}

# Auxiliary function to fit model
bict.fit2 <- function (priornum, missing1, missing2, maximal.mod, IP, eta.hat,
                      ini.index, ini.beta, ini.sig, ini.y0, iters, save, name,
                      null.move.prob, a, b, progress)
{
  missing <- c(missing1, missing2)
  if (is.null(name)) {
    name_RJACC <- "RJACC.txt"
    name_MHACC <- "MHACC.txt"
    name_BETA <- "BETA.txt"
    name_MODEL <- "MODEL.txt"
    name_SIG <- "SIG.txt"
    name_Y0 <- "Y0.txt"
  }
  else {
    name_RJACC <- paste(name, "RJACC.txt", sep = "")
    name_MHACC <- paste(name, "MHACC.txt", sep = "")
    name_BETA <- paste(name, "BETA.txt", sep = "")
    name_MODEL <- paste(name, "MODEL.txt", sep = "")
    name_SIG <- paste(name, "SIG.txt", sep = "")
    name_Y0 <- paste(name, "Y0.txt", sep = "")
  }
  big.X <- maximal.mod$x
  y <- maximal.mod$y
  data <- maximal.mod$data
  curr.index <- ini.index
  curr.X <- big.X[, curr.index == 1]
  curr.ivar <- IP[curr.index == 1, curr.index == 1]
  MODEL <- c()
  BETA <- c()
  curr.beta <- ini.beta
  SIG <- c()
  curr.sig <- ini.sig
  Y0 <- c()
  curr.y0 <- ini.y0
  rj_acc <- c()
  mh_acc <- c()
  counter <- 0
  if (progress) {
    pb <- utils::txtProgressBar(min = 0, max = iters, style = 3)
  }
  while (counter < iters) {
    curr.y <- y
    curr.y[missing] <- curr.y0
    prop_a_model <- conting::prop_mod(curr.index = curr.index, data = data,
                                      maximal.mod = maximal.mod,
                                      null.move.prob = null.move.prob)
    prop.index <- prop_a_model$new.index
    if (prop_a_model$type == "null") {
      new.beta <- conting::iwls_mh(curr.y = curr.y, curr.X = curr.X,
                                   curr.beta = curr.beta,
                                   iprior.var = curr.ivar/curr.sig)
      if (all(new.beta == curr.beta)) {
        mh_acc <- c(mh_acc, 0)
      }
      else {
        mh_acc <- c(mh_acc, 1)
      }
      new.index <- curr.index
    }
    if (prop_a_model$type != "null") {
      rho_bot <- (1 - prop_a_model$null.move.prob)/prop_a_model$total.choices
      prop_a_rev <- conting::prop_mod(curr.index = prop.index, data = data,
                                      maximal.mod = maximal.mod,
                                      null.move.prob = 0)
      rho_top <- (1 - prop_a_model$null.move.prob)/prop_a_rev$total.choices
      rj <- conting::RJ_update(prop.index = prop.index, curr.index = curr.index,
                               curr.beta = curr.beta, eta.hat = eta.hat,
                               curr.y = curr.y,
                               big.X = big.X,
                               proposal.probs = c(rho_top, rho_bot),
                               i.prop.prior.var =
                                 IP[prop.index == 1,
                                    prop.index == 1] / curr.sig,
                               i.curr.prior.var = curr.ivar/curr.sig)
      new.beta <- rj$new.beta
      new.index <- rj$new.index
      if (all(curr.index == new.index)) {
        rj_acc <- c(rj_acc, 0)
      }
      else {
        rj_acc <- c(rj_acc, 1)
      }
    }
    if (priornum == 2) {
      iR <- IP[new.index == 1, new.index == 1]
      curr.sig <- 1/stats::rgamma(n = 1, shape = 0.5 *
                                    (length(new.beta) - 1 + a), rate = 0.5 *
                                    (b + as.vector(matrix(new.beta[-1],
                                                          nrow = 1) %*%
                                                     iR[-1, -1] %*%
                                                     matrix(new.beta[-1],
                                                            ncol = 1))))
    }
    mutar <- exp(as.vector(big.X[, new.index == 1] %*% matrix(new.beta,
                                                              ncol = 1)))
    curr.y0 <- rep(0, length(missing))
    if (length(missing1) > 0){
      curr.y0[1:length(missing1)] <- stats::rpois(n = length(missing1),
                                                  lambda = mutar[missing1])
    }
    if (length(missing2) > 0) {
      for (i in 1:length(missing2)) {
        ppp <- stats::dpois(x = seq(1, y[missing2[i]], by = 1),
                            lambda = mutar[missing2[i]],
                            log = TRUE)
        ppp <- ppp - max(ppp)
        curr.y0[length(missing1) + i] <- sample(x =
                                                seq(1, y[missing2[i]], by = 1),
                                                size = 1, prob = exp(ppp))
      }
    }
    curry <- rep(0, dim(big.X)[2])
    curry[new.index == 1] <- new.beta
    curr.index <- new.index
    curr.beta <- new.beta
    curr.X <- big.X[, curr.index == 1]
    curr.ivar <- IP[curr.index == 1, curr.index == 1]
    BETA <- rbind(BETA, curry)
    SIG <- c(SIG, curr.sig)
    MODEL <- c(MODEL, conting::index2model(new.index))
    Y0 <- rbind(Y0, curr.y0)
    counter <- counter + 1
    if (progress) {
      utils::setTxtProgressBar(pb, counter)
    }
    if (save > 0) {
      if (counter%%save == 0) {
        utils::write.table(file = name_BETA, x = BETA, row.names = FALSE,
                    col.names = FALSE, append = TRUE)
        utils::write.table(file = name_MODEL, x = MODEL, row.names = FALSE,
                    col.names = FALSE, append = TRUE)
        utils::write.table(file = name_SIG, x = SIG, row.names = FALSE,
                    col.names = FALSE, append = TRUE)
        utils::write.table(file = name_Y0, x = Y0, row.names = FALSE,
                    col.names = FALSE, append = TRUE)
        utils::write.table(file = name_RJACC, x = rj_acc, row.names = FALSE,
                    col.names = FALSE, append = TRUE)
        utils::write.table(file = name_MHACC, x = mh_acc, row.names = FALSE,
                    col.names = FALSE, append = TRUE)
        rj_acc <- c()
        mh_acc <- c()
        BETA <- c()
        MODEL <- c()
        SIG <- c()
        Y0 <- c()
      }
    }
  }
  if (progress) {
    close(pb)
  }
  list(BETA = BETA, SIG = SIG, MODEL = MODEL, Y0 = Y0, rj_acc = rj_acc,
       mh_acc = mh_acc)
}
