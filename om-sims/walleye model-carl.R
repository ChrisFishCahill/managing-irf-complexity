vbk=.23
s=.86
age=c(1:20)
n=array(0.,20)  #numbers at age
n[1]=2.
for (a in 2:20){n[a]=n[a-1]*s}  #initial numbers at age
n[20]=n[20]/(1-s)  #set n in age 20 as plus group
w=array(0,20) #weight at age
vul=array(0,20)
w=(1-exp(-vbk*age))^3
mwt=array(0,20)
vul=1/(1+exp(-.5*(age-5)))
mwt=w/(1+exp(-.5*(age-7)))  #maturity x weight for SSB calculation
spro=sum(n*mwt)
cr=6
reca=cr/spro
ro=1
recb=(cr-1)/(ro*spro)
recmult=array(1,200)
by=c(10,20,30,40,50,60,70,80,90,100,110,120,130,140,150,160,170,180,190)
for (t in by){recmult[t]=20}  #big recruitment every 10 yrs
ut=array(.1,200)
yield=array(0,200)
abar=array(0,200)
for (t in 1:200){
  yield[t]=ut[t]*sum(vul*n*w)
  ssb=sum(n*mwt)
  abar[t]=sum(age*n)/sum(n)
  n=s*n*(1-vul*ut[t])
  n[20]=n[20]+n[19]   #update plus group n
  for (a in 19:2){n[a]=n[a-1]}  #move fish up one age
  n[1]=reca*ssb/(1+recb*ssb)*recmult[t]   #put in new recruits
 # if (t/10-int(t/10)=0){n[1]=n[1]*20} #big cohort every 10 years
}
totalyield=sum(yield)
totalutility=sum(yield^.6)
plot(yield,type="l")
plot(abar,type="l")