as = c(4,6,8,10,15,20)
#as = c(2)
nu_as = 1 + as*exp(dnorm(as, log=TRUE)-pnorm(as, lower.tail = FALSE, log=TRUE))


ys = seq(0,30, length.out=1000)

rez = matrix(NA, nrow = length(ys), ncol = length(as))

for (aind in 1:length(as)) {
  a = as[aind]
  nu_a = nu_as[aind]
  
  for (yind in 1:length(ys)) {
    y = ys[yind]
    if(y<a){
      rez[yind, aind] = 1-y^2
    }else{
      rez[yind, aind] = 1 - nu_a
    }
  }
}

matplot(rez,type="l")


b= 2
a = 2.2
#as = c(2)
nu_a = 1 + a*exp(dnorm(a, log=TRUE)-pnorm(a, lower.tail = FALSE, log=TRUE))
nu_b = 1 + b*exp(dnorm(b, log=TRUE)-pnorm(b, lower.tail = FALSE, log=TRUE))


ys = seq(0,30, length.out=1000)

rez = matrix(NA, nrow = length(ys), 1)


for (yind in 1:length(ys)) {
  y = ys[yind]
  if(y<b){
    rez[yind] = 0
  }else if(y<a){
    rez[yind] = nu_b - y^2
  }
  else{
    rez[yind] = nu_b - nu_a
  }
}


matplot(rez,type="l")




b= 5
a = 7
#as = c(2)
nu_a = 1 + a*exp(dnorm(a, log=TRUE)-pnorm(a, lower.tail = FALSE, log=TRUE))
nu_b = 1 + b*exp(dnorm(b, log=TRUE)-pnorm(b, lower.tail = FALSE, log=TRUE))


#ys = seq(0,30, length.out=1000)
N = 100000
#shifts = c(0,1,2,3,4,10,100)
shifts = seq(b-10, a+10, length.out = 30)
#shifts = c(shifts, 3, 5,10,20,50)
rez = matrix(NA, nrow = N, ncol = length(shifts))
rez2 = matrix(NA, nrow = N, ncol = length(shifts))

for (yind in 1:N) {
  xx = rnorm(1)
  for (i in 1:length(shifts)) {
  shift = shifts[i]
  
    #y = ys[yind]
    y = xx + shift
    if(abs(y)<b){
      rez[yind,i] = 0
    }else if(y<a){
      rez[yind,i] = nu_b - y^2
    }
    else{
      rez[yind,i] = nu_b - nu_a
    }
    
    rez2[yind, i] = 0
    if(abs(y)> a){
      rez2[yind, i] =  (y^2 - nu_a)
    }
    if(abs(y)> b){
      rez2[yind, i] = rez2[yind, i] - (y^2 - nu_b)
    }
    
  }
  
}


plot(ecdf(rez2[,1]))
lines(ecdf(rez2[,2]),col=2)
lines(ecdf(rez2[,3]),col=3)


#mgf: 

erf = function(z){
  return(2*pnorm(sqrt(2)*z)-1)
}
erfc = function(z){
  return(1-erf(z))
}


erfclog = function(z){
  return(log(2) + pnorm(sqrt(2)*z,lower.tail=FALSE, log=TRUE))
}
erfcloglowb = function(z){
  return(log(2) -1/2*log(pi)- z^2 - log(z + sqrt(z^2+2)))
}
erfcloglowbmod = function(z,k){
  return(log(2) - z^2 - log(k + sqrt(k^2+2)))
}


mgf_thet = function(a, b,nu_a, nu_b, t, l){
  val = 1 - 1/2*(erfc((b-t)/sqrt(2)) + erfc((b+t)/sqrt(2)))
  val = val + 1/2*exp(l*(nu_b - nu_a) +erfclog((a+t)/sqrt(2)))
  val = val + 1/2*exp(l*(nu_b - nu_a) +erfclog((a-t)/sqrt(2)))
  
  val = val - 1/sqrt(2*l+1) * 1/2 *exp(l*nu_b -l*t^2/(2*l+1) +erfclog((a*(2*l+1)+t)/sqrt(2*(2*l+1))) )
  val = val - 1/sqrt(2*l+1) * 1/2 *exp(l*nu_b -l*t^2/(2*l+1) +erfclog((a*(2*l+1)-t)/sqrt(2*(2*l+1))) )
  val = val + 1/sqrt(2*l+1) * 1/2 *exp(l*nu_b -l*t^2/(2*l+1) +erfclog((b*(2*l+1)-t)/sqrt(2*(2*l+1))) )
  val = val + 1/sqrt(2*l+1) * 1/2 *exp(l*nu_b -l*t^2/(2*l+1) +erfclog((b*(2*l+1)+t)/sqrt(2*(2*l+1))) )
  
  
  # val = 1 - 1/2*(erfc((b-t)/sqrt(2)) + erfc((b+t)/sqrt(2)))
  #   val = val + 1/2*exp(l*(nu_b - nu_a))*1/2*(erfc((a+t)/sqrt(2)) + erfc((a-t)/sqrt(2)))
  # val = val + sqrt(2*l+1) * exp(l*nu_b -l*t^2/(2*l+1)) * 1/2*(
  #   erf((a*(2*l+1)+t)/sqrt(2*(2*l+1))) + erf((a*(2*l+1)-t)/sqrt(2*(2*l+1))) -
  #     erf((b*(2*l+1)+t)/sqrt(2*(2*l+1))) - erf((b*(2*l+1)-t)/sqrt(2*(2*l+1)))
  #   
  # )
  
  return(val)
  
}

mgf_thet22 = function(a, b,nu_a, nu_b, t, l){
  val = 0
  val = val + 1/sqrt(2*l+1) * 1/2 *exp(l*nu_b -l*t^2/(2*l+1) +erfclog((b*(2*l+1)-t)/sqrt(2*(2*l+1))) )
  #val = val + 1/sqrt(2*l+1) * 1/2 *exp(l*nu_b -l*t^2/(2*l+1) +erfclog((b*(2*l+1)+t)/sqrt(2*(2*l+1))) )
  
  
  # val = 1 - 1/2*(erfc((b-t)/sqrt(2)) + erfc((b+t)/sqrt(2)))
  #   val = val + 1/2*exp(l*(nu_b - nu_a))*1/2*(erfc((a+t)/sqrt(2)) + erfc((a-t)/sqrt(2)))
  # val = val + sqrt(2*l+1) * exp(l*nu_b -l*t^2/(2*l+1)) * 1/2*(
  #   erf((a*(2*l+1)+t)/sqrt(2*(2*l+1))) + erf((a*(2*l+1)-t)/sqrt(2*(2*l+1))) -
  #     erf((b*(2*l+1)+t)/sqrt(2*(2*l+1))) - erf((b*(2*l+1)-t)/sqrt(2*(2*l+1)))
  #   
  # )
  
  return(val)
  
}

ls = seq(0,0.5, length.out = 100)


mgfz = matrix(NA, nrow = length(ls), ncol = length(shifts))
mgfz_theory = matrix(NA, nrow = length(ls), ncol = length(shifts))
mgfz_theory22 = matrix(NA, nrow = length(ls), ncol = length(shifts))
maxes = rep(NA, length(ls))
for (ll in 1:length(ls)) {
  maxx = 0
  maxind = 1
  l = ls[ll]
  for (j in 1:length(shifts)) {
    mgfz[ll, j] = mean(exp(l*(rez[,j])))
    mgfz_theory[ll,j] = mgf_thet(a,b, nu_a, nu_b, shifts[j], l)
    mgfz_theory22[ll,j] = mgf_thet22(a,b, nu_a, nu_b, shifts[j], l)
    if(mgfz_theory[ll,j]>maxx){
      maxx = mgfz_theory[ll,j]
      maxind = j
    }
  }
  maxes[ll] = maxind
}

# plot(mgfz[,1])
# lines(mgfz[,2], col=2)
# lines(mgfz[,3], col=2)
# lines(mgfz[,4], col=2)

matplot(ls, mgfz,type="l")
matplot(ls, mgfz_theory,type="l")

matplot(ls, mgfz_theory22,type="l")
maxes

for (ll in 1:length(ls)) {

  l = ls[ll]
  for (j in 2:length(shifts)) {
    if(mgfz_theory[ll, j] > mgfz_theory[ll, j-1]){
      print(sprintf("ERROR on pos j=%d\n", j))
      print(mgfz_theory[ll, j] - mgfz_theory[ll, j-1])
    }
  }
  maxes[ll] = maxind
}



plot(mgfz_theory[length(ls),])
plot(mgfz_theory[1,])




b= 500
a =510
#as = c(2)
nu_a = 1 + a*exp(dnorm(a, log=TRUE)-pnorm(a, lower.tail = FALSE, log=TRUE))
nu_b = 1 + b*exp(dnorm(b, log=TRUE)-pnorm(b, lower.tail = FALSE, log=TRUE))


ts = seq(max(b-5,0), a+5, length.out= 1000)
ts = seq(0, a+5, length.out= 100000)


term1 = rep(NA, length(ts))
term2 = rep(NA, length(ts))
term3 = rep(NA, length(ts))
term4 = rep(NA, length(ts))
term5 = rep(NA, length(ts))
term6 = rep(NA, length(ts))
term7 = rep(NA, length(ts))
term8 = rep(NA, length(ts))
term9 = rep(NA, length(ts))
term10 = rep(NA, length(ts))
term11 = rep(NA, length(ts))
term12 = rep(NA, length(ts))
term13 = rep(NA, length(ts))
term14= rep(NA, length(ts))
term15= rep(NA, length(ts))
term16= rep(NA, length(ts))
term17= rep(NA, length(ts))
term18= rep(NA, length(ts))
term19= rep(NA, length(ts))
term20= rep(NA, length(ts))
l = 0.2

for (i in 1:length(ts)) {
  t = ts[i]
  term1[i] = 1/sqrt(2*pi)*(exp(-(b+t)^2/2) - exp(-(b-t)^2/2))
  
  term2[i] = 1/sqrt(2*pi) * exp(l*(nu_b -nu_a) -(a-t)^2/2) 
  term2[i] = term2[i] -1/sqrt(2*pi) * exp(l*(nu_b -nu_a) -(a+t)^2/2) 
  
  term3[i] = -(-2*l*t) * 1/2* 1/(2*l+1)^(3/2) *exp(l*(nu_b - t^2/(2*l+1)) + erfclog((a*(2*l+1)+t)/sqrt(2*(2*l+1))))
  term3[i] = term3[i] -(-2*l*t) * 1/2* 1/(2*l+1)^(3/2) *exp(l*(nu_b - t^2/(2*l+1)) + erfclog((a*(2*l+1)-t)/sqrt(2*(2*l+1))))
  term3[i] = term3[i] +(-2*l*t) * 1/2* 1/(2*l+1)^(3/2) *exp(l*(nu_b - t^2/(2*l+1)) + erfclog((b*(2*l+1)-t)/sqrt(2*(2*l+1))))
  term3[i] = term3[i] +(-2*l*t) * 1/2* 1/(2*l+1)^(3/2) *exp(l*(nu_b - t^2/(2*l+1)) + erfclog((b*(2*l+1)+t)/sqrt(2*(2*l+1))))
  
  term19[i] = -(-2*l*t) * 1/2* 1/(2*l+1)^(3/2) *exp(l*(nu_b - t^2/(2*l+1)) + erfclog((a*(2*l+1)+t)/sqrt(2*(2*l+1))))
  term19[i] = term19[i] -(-2*l*t) * 1/2* 1/(2*l+1)^(3/2) *exp(l*(nu_b - t^2/(2*l+1)) + erfclog((a*(2*l+1)-t)/sqrt(2*(2*l+1))))
  term18[i] = (-2*l*t) * 1/2* 1/(2*l+1)^(3/2) *exp(l*(nu_b - t^2/(2*l+1)) + erfclog((b*(2*l+1)-t)/sqrt(2*(2*l+1))))
  term18[i] = term18[i] +(-2*l*t) * 1/2* 1/(2*l+1)^(3/2) *exp(l*(nu_b - t^2/(2*l+1)) + erfclog((b*(2*l+1)+t)/sqrt(2*(2*l+1))))

  
  term5[i] =1/sqrt(2*pi)  * 1/(2*l+1) * exp(l*(nu_b-t^2/(2*l+1)) -(a*(2*l+1) + t)^2/(2*l+1)/2)  
  term5[i] =term5[i] - 1/sqrt(2*pi)  * 1/(2*l+1) * exp(l*(nu_b-t^2/(2*l+1)) -(a*(2*l+1)- t)^2/(2*l+1)/2)  
  term6[i] = -1/sqrt(2*pi)  * 1/(2*l+1) * exp(l*(nu_b-t^2/(2*l+1)) -(b*(2*l+1)+t)^2/(2*l+1)/2)  
  term6[i] =term6[i] +1/sqrt(2*pi)  * 1/(2*l+1) * exp(l*(nu_b-t^2/(2*l+1)) -(b*(2*l+1)-t)^2/(2*l+1)/2)  
  
  term4[i] = term5[i] + term6[i]
  
  term7[i] =1/sqrt(2*pi)  * 1/(2*l+1) * exp(l*(nu_b-a^2) -(a + t)^2/2)  
  term7[i] =term7[i] - 1/sqrt(2*pi)  * 1/(2*l+1) * exp(l*(nu_b-a^2) -(a- t)^2/2)
  
  term8[i] =-1/sqrt(2*pi)  * 1/(2*l+1) * exp(l*(nu_b-b^2) -(b + t)^2/2)  
  term8[i] =term8[i] +1/sqrt(2*pi)  * 1/(2*l+1) * exp(l*(nu_b-b^2) -(b- t)^2/2)  

  
  term9[i] = 1/sqrt(2*pi)*(exp(-1/2*(b-t)^2) - exp(-1/2*(b+t)^2))*(exp(l*(nu_b - b^2))/(2*l+1) -1)
  
  term10[i] = 1/sqrt(2*pi)*(exp(-1/2*(a-t)^2) - exp(-1/2*(a+t)^2))*( exp(l*(nu_b - nu_a))- exp(l*(nu_b - a^2))/(2*l+1) )
  
  term11[i] = term3[i]
  
  term12[i] = +(-2*l*t) * 1/2* 1/(2*l+1)^(3/2) *exp(l*(nu_b - t^2/(2*l+1)) + erfclog((b*(2*l+1)-t)/sqrt(2*(2*l+1))))
  term12[i] = term12[i] -(-2*l*t) * 1/2* 1/(2*l+1)^(3/2) *exp(l*(nu_b - t^2/(2*l+1)) + erfclog((b*(2*l+1)+t)/sqrt(2*(2*l+1))))
  
  
  term13[i] = -1/sqrt(2*pi)*exp(log(t) + log(l) - 3/2*log(2*l+1) +
                               l*(nu_b-b^2) -1/2*(b-t)^2) 
  term13[i] = term13[i] - 1/sqrt(2*pi)*exp(log(t) + log(l) - 3/2*log(2*l+1) +
                                           l*(nu_b-b^2) -1/2*(b+t)^2) 
  term13[i] = term13[i] + 1/sqrt(2*pi)*exp(log(t) + log(l) - 3/2*log(2*l+1) +
                                           l*(nu_b-a^2) -1/2*(a-t)^2)
                                         
  term13[i] = term13[i] + 1/sqrt(2*pi)*exp(log(t) + log(l) - 3/2*log(2*l+1) +
                                           l*(nu_b-a^2) -1/2*(a+t)^2)
  
  term14[i] =1/sqrt(2*pi)  * 1/(2*l+1) * exp(l*(nu_b) -(a*(2*l+1) + t)^2/(2*l+1)/2)  
  term14[i] =term14[i] - 1/sqrt(2*pi)  * 1/(2*l+1) * exp(l*(nu_b) -(a*(2*l+1)- t)^2/(2*l+1)/2)  
  term15[i] = -1/sqrt(2*pi)  * 1/(2*l+1) * exp(l*(nu_b) -(b*(2*l+1)+t)^2/(2*l+1)/2)  
  term15[i] =term15[i] +1/sqrt(2*pi)  * 1/(2*l+1) * exp(l*(nu_b) -(b*(2*l+1)-t)^2/(2*l+1)/2)  

  
  
  term16[i] = -(-2*l*t) * 1/2* 1/(2*l+1)^(3/2) *exp(l*(nu_b - t^2/(2*l+1)) + dnorm((a*(2*l+1)+t)/sqrt((2*l+1)),log=TRUE))
  term16[i] = term16[i] -(-2*l*t) * 1/2* 1/(2*l+1)^(3/2) *exp(l*(nu_b - t^2/(2*l+1)) +dnorm((a*(2*l+1)-t)/sqrt((2*l+1)),log=TRUE))
  term16[i] = term16[i] +(-2*l*t) * 1/2* 1/(2*l+1)^(3/2) *exp(l*(nu_b - t^2/(2*l+1)) + erfcloglowb((b*(2*l+1)-t)/sqrt(2*(2*l+1))))
  term16[i] = term16[i] +(-2*l*t) * 1/2* 1/(2*l+1)^(3/2) *exp(l*(nu_b - t^2/(2*l+1)) + erfcloglowb((b*(2*l+1)+t)/sqrt(2*(2*l+1))))
  
  oldt = t[]
  factor = 1
  if((a*(2*l+1)+t)<0){
    t = -t
    factor = -1
  }
  term17[i] = -factor*(-2*l*t) * 1/2* 1/(2*l+1)^(3/2) *exp(l*(nu_b - t^2/(2*l+1)) + dnorm((a*(2*l+1)+t)/sqrt((2*l+1)),log=TRUE))
  term17[i] = term17[i] -factor*(-2*l*t) * 1/2* 1/(2*l+1)^(3/2) *exp(l*(nu_b - t^2/(2*l+1)) +dnorm((a*(2*l+1)-t)/sqrt((2*l+1)),log=TRUE))
  t = oldt
  factor = 1
  if((b*(2*l+1)+t)<0){
    t = -t
    factor = -1
  }
  term17[i] = term17[i] +factor*(-2*l*t) * 1/2* 1/(2*l+1)^(3/2) *exp(l*(nu_b - t^2/(2*l+1)) + erfcloglowbmod((b*(2*l+1)-t)/sqrt(2*(2*l+1)),(b*(2*l+1)+t)/sqrt(2*(2*l+1))))
  #term17[i] = term17[i] +factor*(-2*l*t) * 1/2* 1/(2*l+1)^(3/2) *exp(l*(nu_b - t^2/(2*l+1)) + erfcloglowb((b*(2*l+1)+t)/sqrt(2*(2*l+1))))
    #term13[i] = exp(l*(nu_b - b^2))/(2*l+1) -1 - exp(log(t) - log(t*b) +l*(nu_b -b^2) -1/2*log(2*l+1))
  # term1[i] = 1/sqrt(2*pi)*(exp(-(b+t)^2/2) - exp(-(b-t)^2/2))
  # term2[i] = exp(l*(nu_b -nu_a )) * 1/sqrt(2*pi) * (exp(-(a-t)^2/2) - exp(-(a+t)^2/2))
  # term3[i] = 1/(2*l+1)^(3/2) *exp(l*(nu_b - t^2/(2*l+1))) * (-2*l*t) * 
  #   1/2*(erf((a*(2*l+1)+t)/sqrt(2*(2*l+1))) + erf((a*(2*l+1)-t)/sqrt(2*(2*l+1))) -
  #          erf((b*(2*l+1)+t)/sqrt(2*(2*l+1))) - erf((b*(2*l+1)-t)/sqrt(2*(2*l+1))))
  # term4[i] = 1/(2*l+1) * exp(l*(nu_b-t^2/(2*l+1))) * 1/sqrt(2*pi) * 
  #   ( exp(-(a*(2*l+1) + t)^2/(2*l+1)/2)  -exp(-(a*(2*l+1)- t)^2/(2*l+1)/2)    +
  #       -exp(-(b*(2*l+1) + t)^2/(2*l+1)/2)  +exp(-(b*(2*l+1)- t)^2/(2*l+1)/2)   )
  # term5[i] = 1/(2*l+1) * exp(l*(nu_b-t^2/(2*l+1))) * 1/sqrt(2*pi) * 
  #   ( exp(-(a*(2*l+1) + t)^2/(2*l+1)/2)  -exp(-(a*(2*l+1)- t)^2/(2*l+1)/2))
  # 
  # term6[i] = 1/(2*l+1) * exp(l*(nu_b-t^2/(2*l+1))) * 1/sqrt(2*pi) * 
  #   ( -exp(-(b*(2*l+1) + t)^2/(2*l+1)/2)  +exp(-(b*(2*l+1)- t)^2/(2*l+1)/2)   )
  
  
  term20[i] = dnorm(b+t) - dnorm(b-t) + exp(0.5*(nu_b - nu_a))*(dnorm(a-t) - dnorm(a+t))
  
}

# plot(term1+term2+term3+term4)
# lines(term1+term2+term17+term4)
sum(term1+term2+term17+term4>0)

sum(term1+term2+term3+term4>0)



matplot(ts,cbind(term1, term2, term3,term4, term1+term2+term3+term4),type="l")
legend(x="bottomright",legend=c("term1", "term2", "term3", "term4", "sum"), col=1:5, lty=1:5)


matplot(ts,cbind(term9, term10, term11, term9+term10+term11),type="l")
legend(x="bottomright",legend=c("term9", "term10", "term11","sum"), col=1:4, lty=1:5)


plot(ts, term1+term2+term3+term4)

sum(term1+term2+term3+term4>0)
sum(term1+term4>0)
sum(term2+term3>0)




a = c(4)
#as = c(2)
nu_a = 1 + a*exp(dnorm(a, log=TRUE)-pnorm(a, lower.tail = FALSE, log=TRUE))


ys = seq(-10,10, length.out=1000)
shifts = c(0,2,-2)
rez = matrix(NA, nrow = length(ys), ncol = length(shifts))

#shifts = seq(0, 4, by = 0.5)

for (aind in 1:length(shifts)) {
  shift = shifts[aind]
  
  for (yind in 1:length(ys)) {
    y = ys[yind]
    if(abs(y+shift)<a){
      rez[yind, aind] = 1-(y+shift)^2
    }else{
      rez[yind, aind] = 1 - nu_a
    }
  }
}

matplot(ys,rez,type="l")


# idea: 
# for all values of a (or t)
# ((y+thet)^2 - nu_a ) ind( |y+thet| \geq a) - ((y+thet)^2 -1)  
# \leq (y^2 - nu_a)ind(|y| \geq a) - (y^2-1)
#.. and the above sums can be controlled well. 
# the \leq at least holds in a probabalistic sense ?





a = c(4)
#as = c(2)
nu_a = 1 + a*exp(dnorm(a, log=TRUE)-pnorm(a, lower.tail = FALSE, log=TRUE))


ys = seq(-10,10, length.out=1000)
shifts = c(0,2,-2)
rez = matrix(NA, nrow = length(ys), ncol = length(shifts))

#shifts = seq(0, 4, by = 0.5)

for (aind in 1:length(shifts)) {
  shift = shifts[aind]
  
  for (yind in 1:length(ys)) {
    y = ys[yind]
    if(abs(y+shift)<a){
      rez[yind, aind] = 0
    }else{
      rez[yind, aind] = (y+shift)^2-nu_a
    }
    rez[yind, aind] = rez[yind, aind]  - shift^2
  }
}

matplot(ys,rez,type="l")

a = c(2)
N = 1000000
rez = matrix(NA, nrow = N, ncol = length(shifts))

#shifts = seq(0, 4, by = 0.5)

for (aind in 1:length(shifts)) {
  shift = shifts[aind]
  
  for (yind in 1:N) {
    y = rnorm(1) + shift
    if(abs(y+shift)<a){
      rez[yind, aind] = 0
    }else{
      rez[yind, aind] = (y+shift)^2-nu_a
    }
    rez[yind, aind] = rez[yind, aind]  - shift^2
  }
}

plot(ecdf(rez[,1]))
lines(ecdf(rez[,2]))



mfg_clean = function(l,a){
  val = 0
  nu_a = 1 + a*exp(dnorm(a, log=TRUE)-pnorm(a, lower.tail = FALSE, log=TRUE))
  val = val + exp(-l*nu_a + erfclog(a/sqrt(2)))
  val = val + 1/sqrt(2*l+1) * erf((sqrt(2*l+1)*a)/sqrt(2))
  val = val*exp(l)
  return(val)
}
ls = seq(0, 0.5, length.out=100)
rezz = rep(NA, length(ls))

for (i in 1:length(ls)) {
  l = ls[i]
  rezz[i] = mfg_clean(l, a)
}
plot(ls, rezz)
