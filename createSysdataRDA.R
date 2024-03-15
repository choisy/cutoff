dHash <- list(normal       = stats::dnorm,
              `log-normal` = stats::dlnorm,
              gamma        = stats::dgamma,
              weibull      = stats::dweibull)


pHash <- list(normal       = stats::pnorm,
              `log-normal` = stats::plnorm,
              gamma        = stats::pgamma,
              weibull      = stats::pweibull)


qHash <- list(normal       = stats::qnorm,
              `log-normal` = stats::qlnorm,
              gamma        = stats::qgamma,
              weibull      = stats::qweibull)

save(dHash, pHash, qHash, file = here::here("R/sysdata.rda"))
