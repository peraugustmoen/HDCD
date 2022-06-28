erfc = function(z){
  return(2*pnorm(z*sqrt(2), lower.tail=FALSE))
}

erfclowb = function(z){
  return(2/sqrt(pi)*exp(-z^2)/(z+sqrt(z^2+2)))
}
erfcupb = function(z){
  return(1/sqrt(pi)*exp(-z^2)/(z))
}

mu = function(c, ss,aa){
  nu_aa  = 1 + aa*exp(dnorm(aa, log=TRUE)-pnorm(aa, lower.tail = FALSE, log=TRUE))
  EE1 = 1/2*(c -nu_aa +1)*erfc((aa-sqrt(c))/sqrt(2)) + exp(-(aa-sqrt(c))^2/2)*(aa+sqrt(c))/sqrt(2*pi)
  EE2 =  1/2  *(c-nu_aa + 1)*erfc((aa+sqrt(c))/sqrt(2)) + exp(-(aa+sqrt(c))^2/2)*(aa-sqrt(c))/sqrt(2*pi)
  EE = EE1+EE2
  return(list(ss*EE, EE1, EE2))
}

mudiff = function(c, ss,aa){
  nu_aa  = 1 + aa*exp(dnorm(aa, log=TRUE)-pnorm(aa, lower.tail = FALSE, log=TRUE))
  EE1 = sqrt(c) *( erfc((aa-sqrt(c))/sqrt(2)) + erfc((aa+sqrt(c))/sqrt(2)))
  EE2 = 1/sqrt(2*pi)*(aa^2 +2 - nu_aa)*(exp(-1/2*(aa-sqrt(c))^2) - exp(-1/2*(aa+sqrt(c))^2))
  EE = EE1+EE2
  return(EE)
}
a = 4

thetz = seq(a*0.000001,a/4, by=a*0.000001)

muthetz= thetz[]
mudiffs = thetz[]
mudiff_bound = thetz[]
bb = thetz[]
mudivthet = thetz[]
mudivthet2 = thetz[]
mudivthet3 = thetz[]
temp1 = thetz[]
temp2 = thetz[]
temp3 = thetz[]
temp4 = thetz[]
temp5 = thetz[]
temp6 = thetz[]
temp7 = thetz[]
temp8 = thetz[]
muthetdiff = thetz[]
muthetdifflower = thetz[]
muthetdiffupper = thetz[]
for (i in 1:length(thetz)) {
  rez = mu(thetz[i]^2, 1, a)
  muthetz[i] = rez[[1]]
  mudiffs[i] = mudiff(thetz[i]^2, 1, a)
  
  nu_a  = 1 + a*exp(dnorm(a, log=TRUE)-pnorm(a, lower.tail = FALSE, log=TRUE))
  thet = thetz[i]
  temp1[i] = 2*(nu_a - 1)*(pnorm(a-thet, lower.tail=FALSE) + pnorm(a+thet,lower.tail=FALSE))
  temp2[i] = dnorm(a-thet)*(thet*(a^2-nu_a) - 2*a)
  #temp2[i] = dnorm(a-thet)*(thet*(-2) - 2*a)
  temp3[i] = dnorm(a+thet)*(thet*(nu_a - a^2) - 2*a)
  muthetdiff[i] = temp1[i]+temp2[i]+temp3[i]
  #muthetdifflower[i] = (nu_a - 1)*(erfclowb((a-thet)/sqrt(2)) + erfclowb((a+thet)/sqrt(2))) + temp2[i] + temp3[i]
  muthetdifflower[i] = 2*(nu_a - 1)*(erfc((a)/sqrt(2))) + temp2[i] + temp3[i]
  muthetdiffupper[i] = (nu_a - 1)*(erfcupb((a-thet)/sqrt(2)) + erfcupb((a+thet)/sqrt(2))) + temp2[i] + temp3[i]
  
  
  temp4[i] = 2*(nu_a - 1) + exp(dnorm(a-thet,log=TRUE) - pnorm(a-thet,log=TRUE, lower.tail=FALSE))*
  (thet*(a^2-nu_a)-2*a )
  temp5[i] = 2*(nu_a - 1) + exp(dnorm(a+thet,log=TRUE) - pnorm(a+thet,log=TRUE, lower.tail=FALSE))*
    (thet*(nu_a-a^2 )-2*a )
  
  temp6[i] = (nu_a-1)*(erfc((a-thet)/sqrt(2)) - thet*dnorm(a-thet) + 
                         erfc((a+thet)/sqrt(2)) + thet*dnorm(a+thet)) + 
    thet*(a^2-1)*(dnorm(a-thet) -dnorm(a+thet)) - 
    2*a*(dnorm(a-thet)+ dnorm(a+thet))
  
  temp7[i] = (dnorm(a-thet) - dnorm(a+thet))*((nu_a-1)*(2-thet)+ a^2- 1-2*a) +
    thet*(a^2-1)*(dnorm(a-thet) + dnorm(a+thet))
  
  temp8[i] = (dnorm(a-thet) - dnorm(a+thet))*(a^2 - 1- thet*(nu_a-1)) +
    thet*(a^2-1)*(dnorm(a-thet) + dnorm(a+thet))
  
  #temp3[i] = dnorm(a+thet)*(thet*(2) - 2*a)
  #temp3[i] = dnorm(a+thet)*( - 2*a)
  # mudiff_bound[i] = (nu_a-1)*( erfc((a-thet)/sqrt(2)) - dnorm(a-thet) + 
  #                               erfc((a+thet)/sqrt(2)) - dnorm(a+thet))
  # mudiff_bound[i] = mudiff_bound[i] + dnorm(a-thet) *(thet*(a^2+3) - 2*a) + 
  #                   dnorm(a+thet)*(thet*(1-a^2)-2*a)
  #           
  # bb[i] = dnorm(a, log=TRUE) - pnorm(a, lower.tail=FALSE, log=TRUE)+
  #   log(pnorm(a - thetz[i], lower.tail = FALSE) + pnorm(a +thetz[i], lower.tail = FALSE)) 
  # bb[i] = exp(bb[i]) - pnorm(a-thetz[i])
  # mudivthet[i] = (1-a/thet^2*dnorm(a)/pnorm(a, lower.tail=FALSE))*(pnorm(a-thet,lower.tail=FALSE) + 
  #                                                                  pnorm(a+thet, lower.tail=FALSE)) +
  #   dnorm(a-thet)*(a+thet)/thet^2 + dnorm(a+thet)*(a-thet)/thet^2
  # mudivthet2[i] = (-a/thet^2*dnorm(a)/pnorm(a, lower.tail=FALSE))*(pnorm(a-thet,lower.tail=FALSE) + 
  #                                                                  pnorm(a+thet, lower.tail=FALSE)) +
  #   dnorm(a-thet)*(a)/thet^2 + dnorm(a+thet)*(a)/thet^2
  # 
  # mudivthet3[i] = 1/thet * (dnorm(a-thet) - dnorm(a+thet))
}

plot(temp1, type="l", ylim=c(-0.1,0.1))
lines(temp2)
lines(temp3)
temp1+temp2+temp3
bound = exp(-9/2 *a^2/16)*thetz^2

plot(muthetdiff,type="l")
lines(muthetdifflower,type="l",col=2)
lines(muthetdiffupper,type="l",col=3)

plot(muthetdiff[1:5000], type="l", ylim=c(-1e-12, 1e-12))
lines(muthetdifflower[1:1000], type="l", col=2)

plot(thetz, log(muthetz), ylim = c(max(-1000,min(c(log(muthetz), log(bound)))), max(c(log(muthetz), log(bound)))),type="l")
lines(thetz,log(bound),col=2)
sum(muthetz> bound)

cbind(log(muthetz), log(bound))

plot(thetz, log(muthetz)-log(bound) ,type="l")
gg = log(muthetz)-log(bound)
plot(thetz[1:100], gg[1:100],type="l")

plot(thetz, mudiffs,type="l")
plot(thetz, thetz*mudiffs- (2*muthetz),type="l")
plot(thetz, log(thetz*mudiffs)- log(2*muthetz),type="l")
gg = log(thetz*mudiffs/(2*muthetz))
plot(thetz[1:100], gg[1:100],type="l")

plot(thetz, log(mudiff_bound))
