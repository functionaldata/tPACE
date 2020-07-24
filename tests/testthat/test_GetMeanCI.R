setwd("/Users/xinerzhou/Downloads/tPACE-meanCI/R")
library(tidyverse)
library(stringr)
library(readr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(gridExtra)
library(fdapace)
library(mgcv)
library(cowplot)
library(wesanderson)
 
files<-list.files()
for(i in 1:length(files)){source(files[i])}
 
#####################################3
setwd("/Users/xinerzhou/Dropbox/COVID19/Week5")
source("process.R")
#View(covid_df)
# log transformation for outcome under consideration
df<-covid_df%>%mutate(case=log(1000000*total_cases/pop,base=10),
                      death=log(1000000*(1+total_deaths)/pop,base=10))%>%
  filter(t<=60) 

country.list<-df%>%group_by(country)%>%dplyr::summarise(tmax=max(t))%>%ungroup()
df<-df%>%filter(country %in% country.list$country)




L3 <- MakeFPCAInputs(IDs = df[["country"]], tVec=as.double(df[["t"]]), df[["case"]])
FPCAfit <- FPCA(L3$Ly, L3$Lt, optns = list(usergrid=TRUE,methodMuCovEst = "smooth", kernel = "epan"))
#methodMuCovEst = "smooth", kernel = "epan"
plot(FPCAfit)
confi<-as.data.frame(GetMeanCI(L3$Ly, L3$Lt)$CI)
mean.dat<-data.frame(mean=FPCAfit$mu,
                     t=FPCAfit$workGrid)


confi%>%ggplot(aes(x=CIgrid))+
  geom_ribbon(data=confi,aes(x=CIgrid,ymin = lower, ymax = upper),fill="yellow")+
  geom_line(data=mean.dat,aes(x=t,y=mean ),col="red")+
  geom_hline(yintercept = 0)+
  labs(title="COVID 19 Cumulative Case")+
  xlab("t")+
  ylab("")+
  geom_line(data=df,aes(x=t,y=case,group=country),alpha=0.2)





###########################################
judge<-readr::read_csv('/Users/xinerzhou/Dropbox/FDA/Judge/2018_MQscores/justices.csv',guess_max = 1000)
View(judge)
president<-readr::read_csv('/Users/xinerzhou/Dropbox/FDA/Judge/year_president.csv',guess_max = 1000)%>%
  mutate(senate_r_ratio=senate_r/senate_tt,house_r_ratio=house_r/house_tt)%>%select(year,president,party,senate_r_ratio,house_r_ratio)

# How many justices ? 47
length(unique(judge[["justice"]]))

# rename judges
unique(judge$justiceName)
names<-data.frame(justiceName=c("AFortas","AJGoldberg","AMKennedy","AScalia","BMKavanaugh","BNCardozo",   
                                "BRWhite","CEHughes2","CEWhittaker","CThomas","DHSouter","EKagan",     
                                "EWarren","FFrankfurter","FMurphy","FMVinson","GSutherland","HABlackmun",  
                                "HFStone","HHBurton","HLBlack","JCMcReynolds","JFByrnes","JGRoberts",   
                                "JHarlan2","JPStevens","LDBrandeis","LFPowell","NMGorsuch","OJRoberts",   
                                "PButler","PStewart","RBGinsburg", "RHJackson","SAAlito","SDOConnor",  
                                "SFReed", "SGBreyer","SMinton","SSotomayor","TCClark","TMarshall",  
                                "WBRutledge","WEBurger","WHRehnquist","WJBrennan","WODouglas"),
                  fullName=c("Abe Fortas","Arthur Goldberg","Anthony Kennedy","Antonin Scalia","Brett Kavanaugh","Benjamin N. Cardozo",
                             "Byron White","Charles E. Hughes","Charles E. Whittaker","Clarence Thomas","David Souter","Elena Kagan",
                             "Earl Warren","Felix Frankfurter","Frank Murphy","Fred M. Vinson","George Sutherland","Harry Blackmun",
                             "Harlan F. Stone","Harold H. Burton","Hugo Black","James C. McReynolds","James F. Byrnes","John Roberts",
                             "John M. Harlan II","John P. Stevens","Louis Brandeis","Lewis F. Powell","Neil Gorsuch","Owen J. Roberts",
                             "Pierce Butler","Potter Stewart","Ruth B. Ginsburg","Robert H. Jackson","Samuel Alito","Sandra D. O’Connor",
                             "Stanley F. Reed","Stephen Breyer","Sherman Minton","Sonia Sotomayor","Tom C. Clark","Thurgood Marshall",
                             "Wiley B. Rutledge","Warren E. Burger","William Rehnquist","William J. Brennan","William O. Douglas"),
                  president=c("Lyndon B. Johnson","John F. Kennedy","Ronald W. Reagan","Ronald W. Reagan","Donald Trump","Herbert Hoover",
                              "John F. Kennedy","William H. Taft","Dwight D. Eisenhower","George H. W. Bush","George H. W. Bush","Barack Obama",
                              "Dwight D. Eisenhower","F. D. Roosevelt","F. D. Roosevelt","Harry S. Truman","Warren G. Harding","Richard M. Nixon",
                              "Calvin Coolidge","Harry S. Truman","F. D. Roosevelt","Woodrow Wilson","F. D. Roosevelt","George W. Bush",
                              "Dwight D. Eisenhower","Gerald Ford","Woodrow Wilson","Richard M. Nixon","Donald Trump","Herbert Hoover",
                              "Warren G. Harding","Dwight D. Eisenhower","Bill Clinton","F. D. Roosevelt","George W. Bush","Ronald W. Reagan",
                              "F. D. Roosevelt","Bill Clinton","Harry S. Truman","Barack Obama","Harry S. Truman","Lyndon B. Johnson",
                              "F. D. Roosevelt","Richard M. Nixon","Richard M. Nixon","Dwight D. Eisenhower","F. D. Roosevelt"),
                  president.party=c("D","D","R","R","R","R",
                                    "D","R","R","R","R","D",
                                    "R","D","D","D","R","R",
                                    "R","D","D","D","D","R",
                                    "R","R","D","R","R","R",
                                    "R","R","D","D","R","R",
                                    "D","D","D","D","D","D",
                                    "D","R","R","R","D"),
                  birth=c(1910,1908,1936,1936,1965,1870,1917,1862,1901,1948,1939,1960,1891,1882,1890,1890,
                          1862,1908,1872,1888,1886,1862,1882,1955,1833,1920,1856,1907,1967,1875,1866,1915,1933,1892,
                          1950,1930,1884,1938,1890,1954,1899,1908,1894,1907,1924,1906,1898),
                  tenure=c(1965,1962,1988,1986,2018,1932,1962,1910,1957,1991,1990,2020,1953,1939,1940,1946,1922,1970,1925,1945,1937,1914,1941,2005,
                           1877,1975,1916,1972,2017,1930,1923,1958,1993,1941,2006,1981,1938,1994,1949,2009,1949,1967,1943,1969,1972,1956,1939))

judge<-judge%>%left_join(names,by="justiceName")%>%mutate(age=tenure-birth)

### First, keep all samples longer than 5 years to estimate mean function and covariance function

# number of terms per justice?
term<-judge %>% 
  arrange(justice,term)%>%
  group_by(justice)%>%
  count()
first.term<-judge %>% 
  arrange(justice,term)%>%
  group_by(justice)%>%
  dplyr::summarise(first_term=min(term))
last.term<-judge %>% 
  arrange(justice,term)%>%
  group_by(justice)%>%
  dplyr::summarise(last_term=max(term))

judge<-judge%>%
  left_join(first.term,by="justice")%>%
  left_join(last.term,by="justice")%>%
  mutate(delta_term=as.double(term-first_term),duration=last_term-first_term+1)

judge.dense<-judge%>%filter(delta_term<=20 )%>%filter(delta_term<=20 )%>%filter(duration>=20)
 
denseobj <- MakeFPCAInputs(IDs = judge.dense[["fullName"]], 
                           tVec=as.double(judge.dense[["delta_term"]]), judge.dense[["post_mn"]])

densefit <- FPCA(denseobj$Ly, denseobj$Lt,
                 optns = list(kernel = "epan",methodMuCovEst='cross-sectional'))
CreateCovPlot(densefit, 'Fitted')
plot(densefit)
confi<-as.data.frame(GetMeanCI(denseobj$Ly, denseobj$Lt)$CI)
mean.dat<-data.frame(mean=densefit$mu,
                     t=densefit$workGrid)
 
 
confi%>%ggplot(aes(x=CIgrid))+
  geom_ribbon(data=confi,aes(x=CIgrid,ymin = lower, ymax = upper),fill="yellow")+
  geom_line(data=mean.dat,aes(x=t,y=mean ),col="red")+
  geom_hline(yintercept = 0)+
  labs(title="Ideology")+
  xlab("t")+
  ylab("")+
  geom_line(data=judge.dense,aes(x=delta_term,y=post_mn,group=fullName),alpha=0.2)




########################################
library(MASS)

gaussprocess <- function(from = 0, to = 1, K = function(s, t) {min(s, t)},
                         start = 0, m = 100) {
  # Simulates a Gaussian process with a given kernel
  #
  # args:
  #   from: numeric for the starting location of the sequence
  #   to: numeric for the ending location of the sequence
  #   K: a function that corresponds to the kernel (covariance function) of
  #      the process; must give numeric outputs, and if this won't produce a
  #      positive semi-definite matrix, it could fail; default is a Wiener
  #      process
  #   start: numeric for the starting position of the process
  #   m: positive integer for the number of points in the process to simulate
  #
  # return:
  #   A data.frame with variables "t" for the time index and "xt" for the value
  #   of the process
  
  t <- seq(from = from, to = to, length.out = m)
  Sigma <- sapply(t, function(s1) {
    sapply(t, function(s2) {
      K(s1, s2)
    })
  })
  
  path <- mvrnorm(mu = rep(0, times = m), Sigma = Sigma)
  path <- path - path[1] + start  # Must always start at "start"
  
  return(data.frame("t" = t, "xt" = path))
}
temp1<-c()
temp2<-c()
temp3<-c()
 
for(i in 1:100){
  #i<-1
  temp<-gaussprocess()
  temp1<-append(temp1,temp$t)
  temp2<-append(temp2,temp$xt)
  temp3<-append(temp3,rep(i,100))
}
dat<-data.frame(t=temp1,y=temp2,id=temp3)

dat%>%ggplot(aes(x=t,y=y,group=id))+
  geom_line()

L3 <- MakeFPCAInputs(IDs = dat[["id"]], tVec=as.double(dat[["t"]]), dat[["y"]])
FPCAfit <- FPCA(L3$Ly, L3$Lt, optns = list(usergrid=TRUE,methodMuCovEst = "smooth", kernel = "epan"))
#methodMuCovEst = "smooth", kernel = "epan"
plot(FPCAfit)
confi<-as.data.frame(GetMeanCI(L3$Ly, L3$Lt)$CI)
mean.dat<-data.frame(mean=FPCAfit$mu,
                     t=FPCAfit$workGrid)

confi%>%ggplot(aes(x=CIgrid))+
  geom_ribbon(data=confi,aes(x=CIgrid,ymin = lower, ymax = upper),fill="yellow")+
  geom_line(data=mean.dat,aes(x=t,y=mean ),col="red")+
  geom_hline(yintercept = 0)+
  labs(title="Gaussian Process (Brownian Motion)")+
  xlab("t")+
  ylab("")+
  geom_line(data=dat,aes(x=t,y=y,group=id),alpha=0.2)

 



 