setwd("C:/Users/Dr. Anita A. Joshi/Desktop/HiSeqV2/HNSC_clinicalMatrix")
a2<-read.csv(file.choose())
#ggplot(data = a2, aes(x = clinical_stage))geom_bar()
ggplot(data = a2 , aes(x=a2$clinical_stage))geom_bar()
cls<-a2$clinical_stage
cls
ggplot(data = a2 , aes(x=cls))geom_bar()
a2 %>% ggplot(aes(cls))+geom_bar()
a2 %>% mutate(highlight_flag = ifelse(cls=='Stage IVA',T,F))%>%ggplot(aes(x=cls))+geom_bar(aes(fill=highlight_flag))
a2 %>% mutate(highlight_flag = ifelse(cls=='Stage IVA',T,F))%>%ggplot(aes(x=cls))+geom_bar(aes(fill=highlight_flag))+labs(x='stages',y='count', title =c("Clinical stages in TCGA HNSC"))+coord_flip()
ncs<-a2$person_neoplasm_cancer_status
view(ncs)
plot(ncs)
f<-read.csv(ncs,na.strings =c("","NA"), sep = "\t")
#f<-read.table(ncs,na.strings =c("","NA"), sep = "\t")
#f<-read.csv(file.choose(),na.strings =c("","NA"), sep = "\t")
f1<-read.csv("C:/Users/Dr. Anita A. Joshi/Desktop/HiSeqV2/HNSC_clinicalMatrix/ncs.csv")
f1
##we want read.csv to treat empty rows as NA
f3<-read.csv("C:/Users/Dr. Anita A. Joshi/Desktop/HiSeqV2/HNSC_clinicalMatrix/ncs.csv",na.strings = c("","NA"))
f3
#and then we filter them out using na.omit()function
f3<- f3 %>% na.omit()
plot(f3)
mrna<-read.csv("C:/Users/Dr. Anita A. Joshi/Desktop/HiSeqV2/HiSeqV2")
head(mrna)
glimpse(mrna)
#rem<-function(mrna){
#x<-as.matrix(mrna)
#x<-t(apply(mrna,1,as.numeric))
#r<-as.numeric(apply(x,1,function(i)sum(i==0)))
 #remove<-which(r > dim(mrna)[2]*0.5)
 #return(remove)
#}
rem<-function(mrna)
x<-as.matrix(mrna)
x<-t(apply(mrna,1,as.numeric))
r<-as.numeric(apply(x,1,function(i)sum(i==0)))
remove<-which(r > dim(mrna)[2]*0.5)
return(remove)
os<-file.choose()
require("survival")
os1<-file.choose()
require("survival")
fit<-survfit(Surv(time)~1, data = os1)
if(!require(devtools))install.packages("devtools")
devtools::install_github("kasambara/survminer")
require("survival")
fit<-survfit(Surv(time)~1,data = file.choose())
p1<-file.choose()
View(p1)
p1<-read.csv("C:/Users/Dr. Anita A. Joshi/Desktop/HiSeqV2/HNSC_clinicalMatrix/gender_pri_site")
fit<-survfit(Surv(time)~1,data = p1)
vr<-file.choose()
View(vr)
vr
vr<-read.csv("C:/Users/Dr. Anita A. Joshi/Desktop/HiSeqV2/HNSC_clinicalMatrix/parameters/parameter1")
vr1<-read.csv("D:/IIT/IIT_r/IIT/parameter1")
dim(vr1)
is.na(vr1$X_TIME_TO_EVENT)
count(is.na(vr1$anatomic_neoplasm_subdivision))
vr1<-read.csv("D://IIT//IIT_r//IIT//parameter1.csv")
bn<-lm(vr1$X_TIME_TO_EVENT~vr1$age_at_initial_pathologic_diagnosis)
summary(bn)
vr1<-read.csv("D://IIT//IIT_r//IIT//parameter1.csv")
 mm<-lm(vr1$X_TIME_TO_EVENT~vr1$age_at_initial_pathologic_diagnosis)
 summary(mm)
 bn<-lm(vr1$X_TIME_TO_EVENT~vr1$age_at_initial_pathologic_diagnosis+vr1$anatomic_neoplasm_subdivision)
  summary(bn)
  vb<-lm(vr1$X_TIME_TO_EVENT~vr1$anatomic_neoplasm_subdivision+vr1$age_at_initial_pathologic_diagnosis+vr1$alcohol_history_documented)
  summary(vb)
  vb<-lm(vr1$X_TIME_TO_EVENT~vr1$age_at_initial_pathologic_diagnosis)
  summary(vb)
  library("dplyr", lib.loc="~/R/win-library/3.5")
  cv1<-filter(vr1,(vr1$alcohol_history_documented=="NO"))
  cv2<-filter(vr1,(vr1$alcohol_history_documented=="YES"))
 t.test(cv1$X_TIME_TO_EVENT,cv2$X_TIME_TO_EVENT)
 vb<-lm(vr1$X_TIME_TO_EVENT~vr1$age_at_initial_pathologic_diagnosis)
  summary(vb)
  vb<-lm(vr1$X_TIME_TO_EVENT~vr1$anatomic_neoplasm_subdivision)
  summary(vb)
  c3<-filter(vr1,(vr1$anatomic_neoplasm_subdivision=="anatomic_neoplasm_subdivisionTonsil"))
  c4<-filter(vr1,(vr1$anatomic_neoplasm_subdivision=="anatomic_neoplasm_subdivisionOral Tongue"))
 t.test(c3$X_TIME_TO_EVENT,c4$X_TIME_TO_EVENT)
 c3<-filter(vr1,(vr1$anatomic_neoplasm_subdivision=="anatomic_neoplasm_subdivisionFloor of mouth"))
 t.test(c3$X_TIME_TO_EVENT,c4$X_TIME_TO_EVENT)
 c4<-filter(vr1,(vr1$anatomic_neoplasm_subdivision=="anatomic_neoplasm_subdivisionOralTonsil"))
 t.test(c3$X_TIME_TO_EVENT,c4$X_TIME_TO_EVENT)
 c5<-filter(vr1,(vr1$anatomic_neoplasm_subdivision=="Tonsil"))
 t.test(c5$X_TIME_TO_EVENT,c4$X_TIME_TO_EVENT)
 c6<-filter(vr1,(vr1$anatomic_neoplasm_subdivision=="Buccal Mucosa"))
 t.test(c6$X_TIME_TO_EVENT,c5$X_TIME_TO_EVENT)
 file1<-file.choose()
 View(file1)
  file1<-read.csv("C:/Users/Dr. Anita A. Joshi/Desktop/HiSeqV2/HNSC_clinicalMatrix/hiseq.csv")
  View(file1)
  c10<-filter(file1,(file1$pathologic_M=="M0"))
  dim(c10)
  c11<-filter(file1,(file1$pathologic_M=="M1"))
  dim(c11)
  t.test(c10$X_TIME_TO_EVENT,c11$X_TIME_TO_EVENT)
  c11<-filter(file1,(file1$other_dx=="NO"))
  c21<-filter(file1,(file1$other_dx=="YES"))
   t.test(c11,c21)
   t.test(c11$X_TIME_TO_EVENT,c21$X_TIME_TO_EVENT)
   vb<-lm(file1$X_TIME_TO_EVENT~file1$age_at_initial_pathologic_diagnosis+afile1$pathologic_M)
   file2<-read.csv("C:/Users/Dr. Anita A. Joshi/Desktop/HiSeqV2/HNSC_clinicalMatrix/parameters/parameter3.csv")
   vb<-lm(file2$X_TIME_TO_EVENT~file2$age_at_initial_pathologic_diagnosis+file2$pathologic_M)
   summary(vb)
   file2$pathologic_M<-as.factor(file2$pathologic_M)
   vb<-lm(file2$X_TIME_TO_EVENT~file2$age_at_initial_pathologic_diagnosis+file2$pathologic_M)
   summary(vb)
   plot(file2$age_at_initial_pathologic_diagnosis,file2$X_TIME_TO_EVENT,col="firebrick",pch=9)
   plot(file2$age_at_initial_pathologic_diagnosis,file2$X_TIME_TO_EVENT,col="firebrick",pch=9)
    regs<-lm(file2$X_TIME_TO_EVENT~file2$age_at_initial_pathologic_diagnosis)
    abline(regs,col="steelblue",lwd=3)
    regs<-lm(file2$X_TIME_TO_EVENT~file2$age_at_initial_pathologic_diagnosis,xlim=" age at initial pathological diagnosis", ylim="time to event")
    regs<-lm(file2$X_TIME_TO_EVENT~file2$age_at_initial_pathologic_diagnosis,xlab=" age at initial pathological diagnosis", ylab="time to event")
    plot(file2$age_at_initial_pathologic_diagnosis,file2$X_TIME_TO_EVENT,col="firebrick",pch=9,xlab=" age at initial pathological diagnosis", ylab="time to event")
    regs<-lm(file2$X_TIME_TO_EVENT~file2$age_at_initial_pathologic_diagnosis)
    abline(regs,col="steelblue",lwd=3)
   file3<-read.csv("C:/Users/Dr. Anita A. Joshi/Desktop/HiSeqV2/HNSC_clinicalMatrix/parameters/parameter1.csv") 
head(file3)  
 View(file3)
#col1=X_TIME_TO_EVENT, col2=anatomic_neoplasm_subdivision,col3=age_at_initial_pathologic_diagnosis,col4=alcohol_history_documented
 library(survival)
#data1<-read.csv("file3.csv",header = TRUE,row.names = 1)
x1<-survfit(Surv(file3[,2],file3[,3])~file3[,1])
install.packages("dplyr")
s<-read.csv("E:/head_and_neck_analysis/transfered/HiSeqV2/HNSC_clinicalMatrix/parameters/tcgaparameter1.csv")
View(s)
require("survival")
fit<-survfit(Surv(s$X_TIME_TO_EVENT+s$OS)~s$mrna,data = s)
fit
ggsurvplot(fit,data = s)
s1<-read.csv("E:/head_and_neck_analysis/transfered/HiSeqV2/HNSC_clinicalMatrix/parameters/tcgaparagender.csv")
View(s1)
fit1<-survfit(Surv(s1$X_TIME_TO_EVENT+dat$gen.n)~s1$mrna,data = s1)
gen<-s1$gender
View(gen)
gen1<-factor(c("MALE","FEMALE"))
#View(gen1)
#levels(gen1)
#nlevels(gen1)
dat<-data.frame(gen=sample(c("MALE","FEMALE"),139,replace = TRUE))
View(dat)                    
dat$gen.n<-as.numeric(as.character(factor(dat$gen,levels = c("MALE","FEMALE"),labels = c("2","1"))))
View(dat$gen.n)
ggsurvplot(fit1,data = s1~dat$gen.n)
dat<-data.frame(gen=sample(c("MALE","FEMALE"),139,replace = TRUE))
 View(dat)                    
dat$gen.n<-as.numeric(as.character(factor(dat$gen,levels = c("MALE","FEMALE"),labels = c("2","1"))))
View(dat$gen.n)
View(dat$gen.n)
fit1<-survfit(Surv(s1$X_TIME_TO_EVENT+s1$OS)~s1$gender,data = s1)
View(fit1)
View(s1)
ggsurvplot(fit1,data = s1)
ggsurv<-ggsurvplot(fit1,data = s1,risk.table = TRUE,pval = TRUE)
ggsurv
fit2<-survfit(Surv(s1$mrna+s1$OS)~s1$gender,data = s1)
ggsurv1<-ggsurvplot(fit2,data = s1,risk.table = TRUE,pval = TRUE,legend.lab=c("FEMALE","MALE"))
ggsurv1
s2<-read.csv("E:/head_and_neck_analysis/transfered/HiSeqV2/sox2/tcga_larynx.csv")
View(s2)
fit3<-survfit(Surv(s2$mrna+s2$OS)~s2$gender_num,data=s2)
ggsurv2<-ggsurvplot(fit3,data = s2,pval = TRUE,pval.method = TRUE,legend.lab=c("FEMALE","MALE"),title="Larynx")
ggsurv2
fit4<-survfit(Surv(s2$mrna+s2$OS)~s2$class,data=s2)
ggsurv3<-ggsurvplot(fit4,data = s2,pval = TRUE,pval.method = TRUE,title="Larynx_sox2")
ggsurv3
s3<-read.csv("E:/head_and_neck_analysis/transfered/HiSeqV2/sox2/tcga_oralt.csv")
View(s3)
fit5<-survfit(Surv(s3$X_TIME_TO_EVENT+s3$OS)~s3$class,data=s3)
ggsurv4<-ggsurvplot(fit5,data = s3,pval = TRUE,pval.method = TRUE,title="oraltongue_sox2")
ggsurv4
s4<-read.csv("E:/head_and_neck_analysis/transfered/HiSeqV2/sox2/tcga_bm.csv")
View(s4)
fit6<-survfit(Surv(s4$OS+s4$X_TIME_TO_EVENT)~s4$class,data=s4)
ggsurv5<-ggsurvplot(fit6,data = s4,pval = TRUE,pval.method = TRUE,title="buccal_mucosa_sox2")
ggsurv5
ssn<-read.csv("E:/head_and_neck_analysis/ihc_2_med_low.csv")
View(ssn)
fit7<-survfit(Surv(ssn$OS+ssn$Death)~ssn$Expression,data=ssn)
ggsuv6<-ggsurvplot(fit7,data = ssn,pval = TRUE,pval.method = TRUE,title="sushants_ihc")
ggsuv6
ssn1<-read.csv("E:/head_and_neck_analysis/ifa.csv")
View(ssn1)
fit8<-survfit(Surv(ssn1$OS+ssn1$Death)~ssn1$Expression,data=ssn)
ggsuv7<-ggsurvplot(fit8,data = ssn1,pval = TRUE,pval.method = TRUE,title="sushants_ifa")
ggsuv7
s5<-read.csv("E:/head_and_neck_analysis/transfered/HiSeqV2/sox2/tcga_bm.csv")
View(s5)
fit9<-survfit(Surv(s5$OS+s5$X_TIME_TO_EVENT)~s5$class,data=s5)
ggsurv8<-ggsurvplot(fit9,data = s5,pval = TRUE,pval.method = TRUE,title="buccal_mucosa")
ggsurv8
#fit9<-survfit(Surv(s5$X_TIME_TO_EVENT+s5$class)~s5$OS,data=s5)
s6<-read.csv("E:/head_and_neck_analysis/transfered/HiSeqV2/sox2/tcga_boft.csv")
view(s6)
View(s6)
fit10<-survfit(Surv(s6$OS+s6$X_TIME_TO_EVENT)~s6$class,data=s6)
ggsurv9<-ggsurvplot(fit10,data = s6,pval = TRUE,pval.method = TRUE,title="base of tongue")
ggsurv9
#fit10<-survfit(Surv(s6$OS+s6$class)~s6$X_TIME_TO_EVENT,data=s6)
s7<-read.csv("E:/head_and_neck_analysis/transfered/HiSeqV2/sox2/tcga_ar.csv")
View(s7)
fit11<-survfit(Surv(s7$OS+s7$X_TIME_TO_EVENT)~s7$class,data=s7)
ggsurv10<-ggsurvplot(fit11,data = s7,pval = TRUE,pval.method = TRUE,title="alveolar_ridge")
ggsurv10
s8<-read.csv("E:/head_and_neck_analysis/transfered/HiSeqV2/sox2/tcga_fom.csv")
View(s8)
fit12<-survfit(Surv(s8$OS+s8$X_TIME_TO_EVENT)~s8$class,data=s8)
ggsurv11<-ggsurvplot(fit12,data = s8,pval = TRUE,pval.method = TRUE,title="floor of mouth")
ggsurv11
s9<-read.csv("E:/head_and_neck_analysis/transfered/HiSeqV2/sox2/tcga_hypoph.csv")
View(s9)
fit14<-survfit(Surv(s9$OS+s9$X_TIME_TO_EVENT)~s9$class,data=s9)
ggsurv12<-ggsurvplot(fit14,data = s9,pval = TRUE,pval.method = TRUE,title="hypopharynx")
ggsurv12
s10<-read.csv("E:/head_and_neck_analysis/transfered/HiSeqV2/sox2/tcga_larynx.csv")
View(s10)  
fit15<-survfit(Surv(s10$OS+s10$X_TIME_TO_EVENT)~s10$class,data=s10)
ggsurv13<-ggsurvplot(fit15,data = s10,pval = TRUE,pval.method = TRUE,title="larynx")
ggsurv13
s11<-read.csv("E:/head_and_neck_analysis/transfered/HiSeqV2/stat3/tcga_ar.csv")
fit12<-survfit(Surv(s11$OS+s11$X_TIME_TO_EVENT)~s11$class,data=s11)
ggsurv12<-ggsurvplot(fit12,data = s11,pval = TRUE,pval.method = TRUE,title="alveolar_ridge stat3")
ggsurv12
s11<-read.csv("E:/head_and_neck_analysis/transfered/HiSeqV2/stat3/tcga_bm.csv")
fit121<-survfit(Surv(s11$OS+s11$X_TIME_TO_EVENT)~s11$class,data=s11)
ggsurv12<-ggsurvplot(fit121,data = s11,pval = TRUE,pval.method = TRUE,title="buccal_mucosa stat3")
ggsurv12
s19<-read.csv("E:/head_and_neck_analysis/transfered/HiSeqV2/stat3/tcga_fofm.csv")
fit122<-survfit(Surv(s19$OS+s19$X_TIME_TO_EVENT)~s19$class,data=s19)
ggsurv12<-ggsurvplot(fit122,data = s19,pval = TRUE,pval.method = TRUE,title="florr_of_mouth stat3")
ggsurv12
s20<-read.csv("E:/head_and_neck_analysis/transfered/HiSeqV2/stat3/tcga_hp.csv")
fit122<-survfit(Surv(s20$OS+s20$X_TIME_TO_EVENT)~s20$class,data=s20)
ggsurv12<-ggsurvplot(fit122,data = s20,pval = TRUE,pval.method = TRUE,title="hardpalate stat3")
ggsurv12
s21<-read.csv("E:/head_and_neck_analysis/transfered/HiSeqV2/stat3/tcga_hypoph.csv")
fit122<-survfit(Surv(s21$OS+s21$X_TIME_TO_EVENT)~s21$class,data=s21)
ggsurv12<-ggsurvplot(fit122,data = s21,pval = TRUE,pval.method = TRUE,title="hypopharynx stat3")
ggsurv12
s22<-read.csv("E:/head_and_neck_analysis/transfered/HiSeqV2/stat3/tcga_larynx.csv")
fit122<-survfit(Surv(s22$OS+s22$X_TIME_TO_EVENT)~s22$class,data=s22)
ggsurv12<-ggsurvplot(fit122,data = s22,pval = TRUE,pval.method = TRUE,title="larynx stat3")
ggsurv12
s23<-read.csv("E:/head_and_neck_analysis/transfered/HiSeqV2/stat3/tcga_lip.csv")
fit122<-survfit(Surv(s23$OS+s23$X_TIME_TO_EVENT)~s23$class,data=s23)
ggsurv12<-ggsurvplot(fit122,data = s23,pval = TRUE,pval.method = TRUE,title="lip stat3")
ggsurv12
s24<-read.csv("E:/head_and_neck_analysis/transfered/HiSeqV2/stat3/tcga_oc.csv")
fit122<-survfit(Surv(s24$OS+s24$X_TIME_TO_EVENT)~s24$class,data=s24)
ggsurv12<-ggsurvplot(fit122,data = s24,pval = TRUE,pval.method = TRUE,title="oral cavity stat3")
ggsurv12
s24<-read.csv("E:/head_and_neck_analysis/transfered/HiSeqV2/stat3/tcga_.csv")
fit122<-survfit(Surv(s24$OS+s24$X_TIME_TO_EVENT)~s24$class,data=s24)
ggsurv12<-ggsurvplot(fit122,data = s24,pval = TRUE,pval.method = TRUE,title="oropharynx  stat3")
ggsurv12

s24<-read.csv("E:/head_and_neck_analysis/transfered/HiSeqV2/stat3/tcga_tonsil.csv")
fit122<-survfit(Surv(s24$OS+s24$X_TIME_TO_EVENT)~s24$class,data=s24)
ggsurv12<-ggsurvplot(fit122,data = s24,pval = TRUE,pval.method = TRUE,title="tonsil stat3")
ggsurv12


d1<-read.csv("E:/head_and_neck_analysis/transfered/HiSeqV2/sox2/tcga_oc.csv")
fit31<-survfit(Surv(d1$OS+d1$X_TIME_TO_EVENT)~d1$class,data = d1)
ggsurv15<-ggsurvplot(fit31,data = d1,pval = TRUE,pval.method = TRUE,title="oral cavity sox2")
ggsurv15


d1<-read.csv("E:/head_and_neck_analysis/transfered/HiSeqV2/sox2/tcga_oropharynx.csv")
fit31<-survfit(Surv(d1$OS+d1$X_TIME_TO_EVENT)~d1$class,data = d1)
ggsurv15<-ggsurvplot(fit31,data = d1,pval = TRUE,pval.method = TRUE,title="oropharynx sox2")
ggsurv15


d1<-read.csv("E:/head_and_neck_analysis/transfered/HiSeqV2/sox2/tcga_boft.csv")
fit31<-survfit(Surv(d1$OS+d1$X_TIME_TO_EVENT)~d1$class,data = d1)
ggsurv15<-ggsurvplot(fit31,data = d1,pval = TRUE,pval.method = TRUE,title="base of tongue sox2")
ggsurv15

pj<-read.csv("E:/head_and_neck_analysis/transfered/HiSeqV2/bmi1/tcga_tonsil.csv")

View(pj)
pj
fit123<-survfit(Surv(pj$OS+pj$X_TIME_TO_EVENT)~pj$class,data = pj)
ggsurv64<-ggsurvplot(fit123,data = pj,pval = TRUE,pval.method = TRUE,title="tonsil BMI1")
ggsurv


#############################################################################################

pj1<-read.csv("E:/head_and_neck_analysis/transfered/HiSeqV2/gins2/tcga_bm.csv")
View(pj1)
fit22<-survfit(Surv(pj1$OS+pj1$X_TIME_TO_EVENT)~pj1$class,data = pj1)
ggsurv22<-ggsurvplot(fit22,data = pj1,pval = T,pval.method = T,title="buccal mucosa gins2")
ggsurv22



pj2<-read.csv("E:/head_and_neck_analysis/transfered/HiSeqV2/bmi1/tcga_lip.csv")
fit23<-survfit(Surv(pj2$OS+pj2$X_TIME_TO_EVENT)~pj2$class,data=pj2)
ggsurv23<-ggsurvplot(fit23,data=pj2,pval=T, pval.method = T,title="Lip bmi1")
ggsurv23

zz4<-read.csv("F:/transfered/Genes/cd44/tcga_boft.csv")
fit23<-survfit(Surv(X_TIME_TO_EVENT,OS)~class,data=zz4)
ggsurv23<-ggsurvplot(fit23,data=zz2,pval=TRUE,title="cd44 base of tongue" )
ggsurv23


#######################################################################################                                                                                                   
#######################################################################################

#forG
bp<-read.csv("E:\\head_and_neck_analysis\\transfered\\HiSeqV2\\new_groups\\newgroups_boxplots\\SMPD4\\grp3_G.csv")
p<-ggboxplot(bp,x = "neoplasm_histologic_grade", y = "mrna",color = "neoplasm_histologic_grade",add = "jitter",legend="right", palette = "Dark2", title = "smpd4 grp 3 tumor grade")
p
p + stat_compare_means()
#boxplot(mrna~neoplasm_histologic_grade,data = bp,font.lab=2,col=c("red","blue","green","orange"),xlab="Histological Grade ", ylab = "gene expression", main="stat3 group2 Tumour Grade")
#legend("topright",cex = 0.55, legend=c("G1=18","G2=104","G3=32","G4=4"), title="Samples",fill = c("red","blue","green","orange","cyan","silver"))


#forS
bp1<-read.csv("E:/head_and_neck_analysis/transfered/HiSeqV2/new_groups/newgroups_boxplots/SMPD4/grp1_S.csv")
p<-ggboxplot(bp1,x = "pathologic_stage", y = "mrna",color = "pathologic_stage",add = "jitter",legend="right", palette = "Dark2", title = "smpd4 grp 1 stage")
pgg
p + stat_compare_means()
boxplot(mrna~pathologic_stage,data = bp1,font.lab=2,col=c("pink","purple","yellow","darkgreen"),xlab="Pathological stage ", ylab = "gene expression", main="stat3 group3 stage")
legend("topleft",cex = 0.55, legend=c("I=6","II=24","III=18","IV=103"), title="Samples",fill = c("pink","purple","yellow","darkgreen"))

#forT 
bp2<-read.csv("E:/head_and_neck_analysis/transfered/HiSeqV2/new_groups/newgroups_boxplots/SMPD4/grp3_T.csv")
p<-ggboxplot(bp2,x = "pathologic_T", y = "mrna",color = "pathologic_T",add = "jitter",legend="right", palette = "Dark2", title = "smpd4 grp 3 tumor")
p
p + stat_compare_means()
#boxplot(mrna~pathologic_T,data = bp2,font.lab=2,col=c("chocolate1","grey38","firebrick3","cyan","navyblue","pink"),xlab="Pathological T ", ylab = "gene expression", main="stat3 group1 T")
#legend("topright",  legend=c("T1=8","T2=59","T3=27" ,"T4=11", "T4a=84","T4b=3"), title="Samples",cex = 0.55,fill = c("chocolate1","grey38","firebrick3","cyan","navyblue","pink") )


#forN
bp3<-read.csv("E:/head_and_neck_analysis/transfered/HiSeqV2/new_groups/newgroups_boxplots/SMPD4/grp1_N.csv")
p<-ggboxplot(bp3,x = "pathologic_N", y = "mrna",color = "pathologic_N",add = "jitter",legend="right", palette = "Dark2", title = "smpd4 grp 1 neoplasm")
p
p + stat_compare_means()
#boxplot(mrna~pathologic_N,data = bp3,font.lab=2,col=terrain.colors(7),xlab="pathologic node ", ylab = "gene expression", main="stat3 group1 node")
#legend(  "top" ,  legend=c("N0=75","N1=30","N2=6","N2a=1","N2b=31","N2c=20","N3=3"), title="Samples",cex = 0.45,fill = c(terrain.colors(7)))


View(sfrp)
sfrp<-read.csv("F:/transfered/revised_boxplots/stat3/grp3T.csv")
p<-ggboxplot(sfrp,x = "Tumor_stages", y = "mrna",color = "Tumor_stages",add = "jitter",legend="right", palette = "Dark2", title = "Grp1 T stat3")
p
p + stat_compare_means()

boxplot(mrna~pathologic_stage,data = sfrp1,font.lab=2,col=c("pink","purple","yellow","darkgreen"),xlab="pathological stage ", ylab = "gene expression", main="sfrp1 stage")
legend("topleft",cex = 0.55, legend=c("I=6","II=24","III=18","IV=103"), title="Samples",fill = c("pink","purple","yellow","darkgreen"))








