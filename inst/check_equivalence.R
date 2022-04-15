# expression from pilliat and liu are equivalent :-)

p = 100000000
n = 10000


ss = 1: floor(sqrt(p*log(n^4)))

par(mfrow=c(1,1))
regular = ss*log(exp(1)*p*log(n^4)/ss^2)
new = ss*log(exp(1)-1 + sqrt(p*log(n^4))/ss)
regular[1]
new[1]
plot(ss, regular,type="l")
lines(ss, new,col=2)

lines(ss, exp(1)*ss*log(exp(1)-1 + sqrt(p*log(n^4))/ss),col=3)

sum(regular<new)
sum(regular >exp(1)*new)

