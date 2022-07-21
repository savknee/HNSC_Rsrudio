##### RUN THIS CODE FOR MAKING KM PLOTS ############################
###### use packages "survival" and "survminer" before you run this code #######
ab<-read.csv("path/filename.csv")
fit101<-survfit(Surv(X_TIME_TO_EVENT,OS)~class,data=ab)
ggsurv101<-ggsurvplot(fit101,data=zz,pval=TRUE,title="title you want to give" )
ggsurv101



########## FOR BOXPLOTS ##########
###### LOAD PACKAGES "ggplot2" , "ggpubr" #######
##### IF YOU WANT NORMAL BOXPLOTS NO NEED FOR ANY PACKAGES THEN ##########

#forG,S,T,N,M stages
bp<-read.csv("path/filename.csv")
p<-ggboxplot(bp,x = "name in excelsheet, name must match the x-axis and name in excelsheet", y = "name in excelsheet, name must match the y-axis and name in excelsheet",color = "name of x-axis",add = "jitter",legend="right", palette = "Dark2", title = "title you want")
p
p + stat_compare_means()

#####if you want a normal boxplot ###########

#(eg: can be used for previous labels as well ie above ggplot code, y-axis: mrna, x-axis:neoplasm_histologic_grade)
## G1-G4: they are sample sizes
boxplot(mrna~neoplasm_histologic_grade,data = bp,font.lab=2,col=c("red","blue","green","orange"),xlab="Histological Grade ", ylab = "gene expression", main="title you want to give")
legend("topright",cex = 0.55, legend=c("G1=18","G2=104","G3=32","G4=4"), title="Samples",fill = c("red","blue","green","orange","cyan","silver"))