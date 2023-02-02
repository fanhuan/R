argument_generator <- function(...) {
    argument <-list(...)
    arguments <- expand.grid(argument)
    return(arguments)
}

#USAGE: argument_generator(ID=c('1604','1605'),r=c('0.0',0.05,0.1),id=c(1,2,3,4))
ID=seq(16052200,16052299,1)
r='0.0'
k=15
n=1
G=18
t=6
L=seq(1000,2000,200)
args_May22<-argument_generator(ID=ID,r=r,k=k,n=n,G=G,t=t,L=L)
write.table(args_May22,file='Arguments_SBA_single_160522.txt',row.names = F,quote = F)

###############160526###################
ID=seq(16052600,16052699,1)
r=c('0.0','0.01')
k=15
n=1
G=30
t=5
L=seq(1000,5000,1000)
args_May26<-argument_generator(ID=ID,r=r,k=k,n=n,G=G,t=t,L=L)
write.table(args_May26,file='Arguments_SBA_pairend_160526.txt',row.names = F,quote = F)

###############160527###################
ID=seq(16052300,16052399,1)
r=c('0.0','0.01')
k=15
n=1
G=30
t=5
L=seq(1000,5000,1000)
args_May27<-argument_generator(ID=ID,r=r,k=k,n=n,G=G,t=t,L=L)
write.table(args_May27,file='Arguments_SBA_pairwise_160527.txt',row.names = F,quote = F)

###############160530###################
ID=seq(16052600,16052699,1)
r=c('0.0','0.01')
k=15
n=1
G=30
t=5
L=c(500,1500,2500)
args_May30<-argument_generator(ID=ID,r=r,k=k,n=n,G=G,t=t,L=L)
write.table(args_May30,file='Arguments_SBA_pairend_160530.txt',row.names = F,quote = F)

###############160613###################
#RAD with d=0.1
#simrrls
ID=seq(16061300,16061399,1)
L=10000
dm=1
ds=0
t=c('simrrls_default.0.1.tre')
args_June13<-argument_generator(ID=ID,L=L,dm=dm,ds=ds,Tree=t)
write.table(args_June13,file='Arguments_simrrls_RAD_160613.txt',row.names = F,quote = F)
#SBA_single
ID=seq(16061300,16061399,1)
r=c('0.0','0.01','0.05')
k=15
n=1
G=18
t=6
L=seq(2000,10000,2000)
args_June13<-argument_generator(ID=ID,r=r,k=k,n=n,G=G,t=t,L=L)
write.table(args_June13,file='Arguments_SBA_single_160613.txt',row.names = F,quote = F)


#RAD with d=0.01
#simrrls
ID=seq(160613100,160613199,1)
L=10000
dm=1
ds=0
t=c('simrrls_default.0.01.tre')
args_June13<-argument_generator(ID=ID,L=L,dm=dm,ds=ds,Tree=t)
write.table(args_June13,file='Arguments_simrrls_RAD_160613_1.txt',row.names = F,quote = F)
#SBA_single
ID=seq(160613100,160613199,1)
r=c('0.0','0.01','0.05')
k=15
n=1
G=18
t=6
L=seq(2000,10000,2000)
args_June13<-argument_generator(ID=ID,r=r,k=k,n=n,G=G,t=t,L=L)
write.table(args_June13,file='Arguments_SBA_single_160613_1.txt',row.names = F,quote = F)


#ddRAD with d=0.1
ID=seq(160613200,160613299,1)
L=5000
dm=1
ds=0
t=c('simrrls_default.0.1.tre')
args_June13<-argument_generator(ID=ID,L=L,dm=dm,ds=ds,Tree=t)
write.table(args_June13,file='Arguments_simrrls_ddRAD_160613.txt',row.names = F,quote = F)
#SBA_pairend
ID=seq(160613200,160613299,1)
r=c('0.0','0.01','0.05')
k=15
n=1
G=30
t=6
L=seq(1000,5000,1000)
args_June13<-argument_generator(ID=ID,r=r,k=k,n=n,G=G,t=t,L=L)
write.table(args_June13,file='Arguments_SBA_pairend_160613.txt',row.names = F,quote = F)
#ddRAD with d=0.01
#simrrls
ID=seq(160613300,160613399,1)
L=5000
dm=1
ds=0
t=c('simrrls_default.0.01.tre')
args_June13<-argument_generator(ID=ID,L=L,dm=dm,ds=ds,Tree=t)
write.table(args_June13,file='Arguments_simrrls_ddRAD_160613_1.txt',row.names = F,quote = F)
#SBA_pairend
ID=seq(160613200,160613299,1)
r=c('0.0','0.01','0.05')
k=15
n=1
G=5
t=5
L=seq(1000,5000,1000)
args_June13<-argument_generator(ID=ID,r=r,k=k,n=n,G=G,t=t,L=L)
write.table(args_June13,file='Arguments_SBA_pairend_160613.txt',row.names = F,quote = F)
#pairwise for RAD d=0.1
ID=seq(16061300,16061399,1)
r=c('0.0','0.01','0.05')
k=15
n=1
G=30
t=5
L=seq(2000,10000,2000)
args_June13<-argument_generator(ID=ID,r=r,k=k,n=n,G=G,t=t,L=L)
write.table(args_June13,file='Arguments_SBA_pairwise_160613.txt',row.names = F,quote = F)
#pairwise for RAD d=0.01
ID=seq(160613100,160613199,1)
r=c('0.0','0.01','0.05')
k=15
n=1
G=30
t=5
L=seq(2000,10000,2000)
args_June13<-argument_generator(ID=ID,r=r,k=k,n=n,G=G,t=t,L=L)
write.table(args_June13,file='Arguments_SBA_pairwise_160613_1.txt',row.names = F,quote = F)
#pairwise for ddRAD d=0.1
ID=seq(160613200,160613299,1)
r=c('0.0','0.01','0.05')
k=15
n=1
G=5
t=5
L=seq(1000,5000,1000)
args_June13<-argument_generator(ID=ID,r=r,k=k,n=n,G=G,t=t,L=L)
write.table(args_June13,file='Arguments_SBA_pairwise_160613_2.txt',row.names = F,quote = F)
#pairwise for ddRAD d=0.01
ID=seq(160613300,160613399,1)
r=c('0.0','0.01','0.05')
k=15
n=1
G=5
t=5
L=seq(1000,5000,1000)
args_June13<-argument_generator(ID=ID,r=r,k=k,n=n,G=G,t=t,L=L)
write.table(args_June13,file='Arguments_SBA_pairwise_160613_3.txt',row.names = F,quote = F)

#coverage for RAD d=0.01
ID=seq(160711000,160711099,1)
L=20000
dm=2
ds=1
t='simrrls_default.0.01.tre'
e=0.001
args_July11<-argument_generator(ID=ID,L=L,dm=dm,ds=ds,tree=t,e=e)
write.table(args_July11,file='Arguments_simrrls_RAD_160711.txt',row.names = F,quote = F)

##SBA for RAD d=0.01
ID=seq(160711000,160711099,1)
r=c(0.01,0.05)
k=15
n=2
G=18
t=4
L=c(10000,15000,20000)
args_July11<-argument_generator(ID=ID,r=r,k=k,n=n,G=G,t=t,L=L)
write.table(args_July11,file='Arguments_SBA_single_160711.txt',row.names = F,quote = F)
