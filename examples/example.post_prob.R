
H0 <- structure(list(logml = -20.8084543022433, niter = 4, method = "normal"),
                .Names = c("logml", "niter", "method"), class = "bridge")
H1 <- structure(list(logml = -17.9623077558729, niter = 4, method = "normal"),
                .Names = c("logml", "niter", "method"), class = "bridge")
H2 <- structure(list(logml = -19, niter = 4, method = "normal"),
                .Names = c("logml", "niter", "method"), class = "bridge")


post_prob(H0, H1, H2)
post_prob(H1, H0)

## all produce the same (only names differ):
post_prob(H0, H1, H2)
post_prob(H0$logml, H1$logml, H2$logml)
post_prob(c(H0$logml, H1$logml, H2$logml))
post_prob(H0$logml, c(H1$logml, H2$logml))


