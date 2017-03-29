
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
post_prob(H0$logml, c(H1$logml, H2$logml), model_names = c("H0", "H1", "H2"))


### with bridge list elements:
H0L <- structure(list(logml = c(-20.8088381186739, -20.8072772698116,
-20.808454454621, -20.8083419072281, -20.8087870541247, -20.8084887398113,
-20.8086023582344, -20.8079083169745, -20.8083048489095, -20.8090050811436
), niter = c(4, 4, 4, 4, 4, 4, 4, 4, 4, 4), method = "normal",
    repetitions = 10), .Names = c("logml", "niter", "method",
"repetitions"), class = "bridge_list")

H1L <- structure(list(logml = c(-17.961665507006, -17.9611290723151,
-17.9607509604499, -17.9608629535992, -17.9602093576442, -17.9600223300432,
-17.9610157118017, -17.9615557696561, -17.9608437034849, -17.9606743200309
), niter = c(4, 4, 4, 4, 4, 4, 4, 4, 3, 4), method = "normal",
    repetitions = 10), .Names = c("logml", "niter", "method",
"repetitions"), class = "bridge_list")

post_prob(H1L, H0L)
post_prob(H1L, H0L, H0) # last element recycled with warning.

