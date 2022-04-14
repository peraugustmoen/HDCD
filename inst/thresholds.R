## setup 1: 
n = 1000
p = 2000

ss = 1: floor(sqrt(p*log(n^4*log(p))))

rho = 18*c(pmax(ss*log(exp(1)*p*log(n^4)/ss^2), log(n^4)), sqrt(p*log(n^4)))
plot(rho)
as = c(log(exp(1)*p*log(n^4*log(p))/ss^2), 0)
plot(as)
t = c(pmax(ss/exp(1), log(n^4*log(p))), sqrt(p*log(n^4*p)))

plot(rho, type="l", ylim=c(min(c(min(t), min(rho))), max(c(max(t), max(rho)))))
lines(t, type="l", col=2)
lines(as*t, type="l", col=3)
sum(t>rho)
rho[1]
tt[1]


n = 1000
p = 2000

ss = 1: floor(sqrt(p*log(n^4)))

rho = 2*c(pmax(ss*log(exp(1)*p*log(n^4)/ss^2), log(n^4)), sqrt(p*log(n^4)))
plot(rho)
as = c(log(exp(1)*p*log(n^4)/ss^2), 0)
plot(as)
t = c(pmax(ss/exp(1), log(n^4*log(p))), sqrt(p*log(n^4*p)))

plot(rho, type="l", ylim=c(min(c(min(t), min(rho))), max(c(max(t), max(rho)))))
lines(t, type="l", col=2)
sum(t>rho)
rho[1]
tt[1]




n = 1000
p = 2000

ss = 1: floor(sqrt(p*log(n^4)))

rho = 2*pmax(ss*log(exp(1)*p*log(n^4)/ss^2), log(n^4))

t = pmax(ss^2*sqrt(log(n^4*log(p))/exp(2)/p/log(n^4)^2), log(n^4*log(p)))

plot(rho, type="l", ylim=c(min(c(min(t), min(rho))), max(c(max(t), max(rho)))))
lines(t, type="l", col=2)
sum(t>rho)
rho[1]
tt[1]



# 
# 
# # with a * 2:
# tt = pmax(ss*sqrt(log(n^4*log(p))/exp(1)/log(n^4)), log(n^4*log(p)))
# lines(tt, type="l")