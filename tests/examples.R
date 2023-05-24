library(TFunHDDC)
set.seed(1027)
#simulataed univariate data
data = genModelFD(ncurves=300, nsplines=35, alpha=c(0.9,0.9,0.9),
                  eta=c(10, 7, 17))
plot(data$fd, col = data$groupd)
clm = data$groupd
model1=c("AkjBkQkDk", "AkjBQkDk", "AkBkQkDk", "ABkQkDk", "AkBQkDk", "ABQkDk")
t1<-tfunHDDC(data$fd,K=3,threshold=0.2,init="kmeans",nb.rep=1,
             dfconstr="no", dfupdate="numeric", model=model1[1], itermax = 10)
table(clm, t1$class)
###############example when some classifications are known
if (FALSE) { # ommited due to long run times
  known1=rep(NA,1,300)
  known1[1]=clm[1]
  known1[103]=clm[103]
  known1[250]=clm[250]
  t2<-tfunHDDC(data$fd,K=3,threshold=0.2,init="kmeans",nb.rep=10,dfconstr="no", 
               dfupdate="numeric", model=model1[1],known=known1)
  table(clm, t2$class)
  ################### example when some classifications are known and given in training
  known1=rep(NA,1,300)
  known1[1:100]=rep(3,1,50)
  t3<-tfunHDDC(data$fd,K=3,threshold=0.2,init="kmeans",nb.rep=10,dfconstr="no", 
               dfupdate="numeric", model=model1[1],known=known1)
  table(clm, t3$class)
}
####################classification example with predictions
training=c(1:50,101:150, 201:250)
test=c(51:100,151:200, 251:300)
known1=clm[training]
t4<-tfunHDDC(data$fd[training],K=3,threshold=0.2,init="kmeans",nb.rep=1,
             dfconstr="no", dfupdate="numeric", model=model1[1],known=known1, 
             itermax = 10)
table(clm[training], t4$class)
p1<-predict.tfunHDDC(t4,data$fd[test] )
table(clm[test], p1$class)
###########################NOX data
data1=fitNOxBenchmark(15)
plotNOx(data1)

if (FALSE) { # ommited due to long run times
  t1<-tfunHDDC(data1$fd,K=2,threshold=0.6,init="kmeans",nb.rep=20,dfconstr="no", 
               model=model1)
  #t2<-tfunHDDC(data1$fd,K=2,threshold=0.4,init="kmeans",nb.rep=20, model=c("AkjBkQkDk", "AkjBQkDk", "AkBkQkDk", "ABkQkDk", "AkBQkDk", "ABQkDk"))
  #t3<-tfunHDDC(data1$fd,K=2,threshold=0.2,init="kmeans",nb.rep=20, model=c("AkjBkQkDk", "AkjBQkDk", "AkBkQkDk", "ABkQkDk", "AkBQkDk", "ABQkDk"))
  #t3<-tfunHDDC(data1$fd,K=2,threshold=0.05,init="kmeans",nb.rep=20, model=c("AkjBkQkDk", "AkjBQkDk", "AkBkQkDk", "ABkQkDk", "AkBQkDk", "ABQkDk"))
  
  table(data1$groupd, t1$class)
  #table(data1$groupd, t2$class)
  #table(data1$groupd, t2$class)
  #table(data1$groupd, t3$class)
  #table(data1$groupd, t4$class)
  ###example for prediction
  training=c(1:50)
  test=c(51:115)
  known1=data1$groupd[training]
  t1<-tfunHDDC(data1$fd[training],K=2,threshold=0.6,init="kmeans",nb.rep=10,
               dfconstr="no", model=model1,known=known1) 
  table(data1$groupd[training], t1$class)
  p1<-predict.tfunHDDC(t1,data1$fd[test] )
  table(data1$groupd[test], p1$class)
}
############################multivariate simulated data
set.seed(2341)
conTrig <- genTriangles()
# plotTriangles(conTrig)
cls = conTrig$groupd # groups 5 and 6 (contaminated) go into 1 and 3 respectively
res_s = tfunHDDC(conTrig$fd, K=4, dfconstr="no", dfupdate="numeric", 
                 model="ABKQKDK", init="kmeans", threshold=0.2, nb.rep=1, 
                 itermax=10)
table(cls, res_s$class)


