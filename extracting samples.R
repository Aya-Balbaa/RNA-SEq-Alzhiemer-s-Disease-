cts3=txi_lengthScaled$counts
cts3=as.data.frame(cts3)

##extractin WT samples
##WT aging 9m , 18m, 24m
wt2m=cbind(cts3$GEN181107_S_721877_lib70400_1337_5, cts3$GEN190822_E_896359_lib93587_1689_3, cts3$GEN190822_E_896360_lib93603_1689_4, cts3$GEN181107_S_170824_3_lib70422_1337_7, cts3$GEN181107_S_170874_2_lib70411_1337_6, cts3$GEN190822_E_193953_1_lib93601_1689_4, cts3$GEN190822_E_901251_lib93579_1689_2, cts3$GEN190822_E_905258_lib93566_1689_1, cts3$GEN190822_E_905263_lib93572_1689_2, cts3$GEN181107_S_170525_1_lib70439_1337_8, cts3$GEN181107_S_178525_2_lib70405_1337_5, cts3$GEN181107_S_195899_2_lib70377_1337_3)
batchwt2m=factor(c(1,2,2,1,1,2,2,2,2,1,1,1))
dim(wt2m)
colnames(wt2m)=c("wt2.1", "wt2.2", "wt2.3", "wt2.4", "wt2.5", "wt2.6", "wt2.7", "wt2.8", "wt2.9", "wt2.10", "wt2.11", "wt2.12")


wt4m=cbind(cts3$GEN190822_E_846061_lib93573_1689_2, cts3$GEN190822_E_846062_lib93597_1689_4, cts3$GEN190822_E_869576_lib93612_1689_5, cts3$GEN181107_S_164700_5_lib70409_1337_5, cts3$GEN181107_S_164700_6_lib70421_1337_6, cts3$GEN181107_S_187153_1_lib70438_1337_8, cts3$GEN190822_E_887611_lib93640_1689_7, cts3$GEN190822_E_887612_lib93618_1689_6, cts3$GEN190822_E_897369_lib93605_1689_4, cts3$GEN181107_S_182828_1_lib70397_1337_4,cts3$GEN181107_S_W545_1_lib70391_1337_4, cts3$GEN181107_S_W545_2_lib70429_1337_7)
batchwt4m=factor(c(2,2,2,1,1,1,2,2,2,1,1,1))
colnames(wt4m)=c("wt4.1", "wt4.2", "wt4.3", "wt4.4", "wt4.5", "wt4.6", "wt4.7", "wt4.8", "wt4.9", "wt4.10", "wt4.11", "wt4.12")

wt9m=cbind(cts3$GEN181107_S_532886_lib70419_1337_6, cts3$GEN181107_S_532887_lib70408_1337_5, cts3$GEN181107_S_532888_lib70364_1337_2, cts3$GEN190822_E_661572_lib93561_1689_1, cts3$GEN190822_E_691785_lib93578_1689_2, cts3$GEN190822_E_691785_lib93578_1689_2, cts3$GEN190822_E_806516_lib93584_1689_3, cts3$GEN190822_E_806518_lib93583_1689_3, cts3$GEN181107_S_125173_4_lib70413_1337_6, cts3$GEN181107_S_125173_5_lib70427_1337_7, cts3$GEN181107_S_W523_1_lib70381_1337_3, cts3$GEN181107_S_W526_4_lib70436_1337_8 )
batchw9m=factor(c(2,2,2,2,2,2,2,2,1,1,1,1))
wt92=round(wt9m)
dim(wt14)
colnames(wt9m)=c("wt9.1", "wt9.2", "wt9.3", "wt9.4", "wt9.5", "wt9.6", "wt9.7", "wt9.8", "wt9.9", "wt9.10", "wt9.11", "wt9.12")
genderw9m=factor(c(rep("f", 6), rep("m", 6)))

wt14=cbind(cts3$GEN190822_E_625485_lib93642_1689_8, cts3$GEN190822_E_625486_lib93617_1689_5, cts3$GEN190822_E_663090_lib93611_1689_5,cts3$GEN190822_E_706622_lib93575_1689_2, cts3$GEN181107_S_W545_2_lib70429_1337_7, cts3$GEN190822_E_127620_4_lib93577_1689_2, cts3$GEN190822_E_533997_lib93586_1689_3, cts3$GEN190822_E_647846_lib93563_1689_1, cts3$GEN190822_E_675512_lib93580_1689_2, cts3$GEN181107_S_W525_1_lib70417_1337_6, cts3$GEN181107_S_W525_2_lib70399_1337_5, cts3$GEN181107_S_W526_2_lib70393_1337_4)
batchwt14=factor(c(2,2,2,2,2,3,2,2,2,1,1,1))
colnames(wt14)=c("wt14.1", "wt14.2", "wt14.3", "wt14.4", "wt14.5", "wt14.6", "wt14.7", "wt14.8", "wt14.9", "wt14.10", "wt14.11", "wt14.12")
genderw14m=factor(c(rep("f", 6), rep("m", 6)))
dim(wt14)
wt18m=cbind(cts3$GEN190822_E_585935_lib93648_1689_8, cts3$GEN190822_E_585939_lib93599_1689_4, cts3$GEN190822_E_597014_lib93602_1689_4, cts3$GEN190822_E_624681_lib93581_1689_2, cts3$GEN190822_E_633832_lib93568_1689_1, cts3$`NG-24492_695434_lib386555_6738_1`, cts3$GEN190822_E_526841_lib93625_1689_6, cts3$GEN190822_E_526843_lib93588_1689_3, cts3$GEN190822_E_HTau142_1_lib93598_1689_4, cts3$GEN181107_S_Htau142_LR_lib70431_1337_7, cts3$GEN181107_S_Htau148_R_lib70444_1337_8, cts3$GEN181107_S_Lau145_L_lib70384_1337_3 )
batchw18m=factor(c(1,2,2,2,1,1,1,1,1,1,2,2))
dim(wt18m)
colnames(wt18m)=c("wt18.1", "wt18.2", "wt18.3", "wt18.4", "wt18.5", "wt18.6", "wt18.7", "wt18.8", "wt18.9", "wt18.110", "wt18.11", "wt18.12")
genderw18m=factor(c(rep("f", 6), rep("m", 6)))

wt24m=cbind(cts3$`NG-24492_585932_lib386556_6738_1`, cts3$`NG-27329_733441_lib478299_7389_2`,cts3$`NG-27329_808957_lib478298_7389_2`, cts3$`NG-27329_838802_lib478297_7389_2`,cts3$GEN190822_E_416963_lib93645_1689_8, cts3$`NG-24492_416964_lib386552_6738_1`, cts3$GEN190822_E_421934_lib93604_1689_4, cts3$GEN190822_E_421936_lib93567_1689_1, cts3$`NG-24492_453986_lib386551_6738_1`, cts3$`NG-24492_473954_lib386553_6738_1`)
batchw24m=factor(c(3,4,4,4,2,3,2,2,3,3))
colnames(wt24m)=c("wt24.1", "wt24.2", "wt24.3", "wt24.4", "wt24.5", "wt24.6", "wt24.7", "wt24.8", "wt24.9", "wt24.10")
genderw24m=factor(c(rep("f", 4), rep("m", 6)))
dim(wt24m)

wtage=cbind(wt9m, wt14, wt18m, wt24m)

wtage3=cbind(wt9m, wt18m, wt24m)
wtage3=round(wtage3)
allwt=cbind(wt2m, wt4m, wt9m, wt14, wt18m, wt24m)
dim(allwt)
batchallwt=factor(c(1,2,2,1,1,2,2,2,2,1,1,1,2,2,2,1,1,1,2,2,2,1,1,1,2,2,2,2,2,2,2,2,1,1,1,1,2,2,2,2,2,3,2,2,2,1,1,1,1,2,2,2,1,1,1,1,1,1,2,2,3,4,4,4,2,3,2,2,3,3))


##NLF
##4months
nlf4m=cbind(cts3$GEN181107_S_719581_lib70435_1337_8, cts3$GEN190822_E_869576_lib93612_1689_5, cts3$GEN190822_E_880491_lib93634_1689_7, cts3$GEN190822_E_907227_lib93559_1689_1, cts3$`NG-24492_972964_lib386561_6738_1`, cts3$GEN181107_S_176825_1_lib70443_1337_8, cts3$GEN181107_S_176825_2_lib70374_1337_3, cts3$GEN181107_S_713818_lib70390_1337_4, cts3$GEN181107_S_730373_lib70380_1337_3, cts3$GEN190822_E_869949_lib93616_1689_5, cts3$GEN190822_E_880490_lib93650_1689_8, cts3$GEN190822_E_882184_lib93569_1689_1)
colnames(nlf4m)=c("nlf41", "nlf42", "nlf43", "nlf44", "nlf45", "nlf46", "nlf47", "nlf48", "nlf49", "nlf410", "nlf411", "nlf412")
dim(nlf4m)
batchnlf4=factor(c(1,2,2,2,3,1,1,1,1,2,2,2))

##9 month
nlf9=cbind(cts3$GEN181107_S_635250_lib70366_1337_2, cts3$GEN181107_S_639498_lib70402_1337_5, cts3$`NG-24492_812244_lib386623_6738_6`, cts3$GEN190822_E_822441_lib93628_1689_6, cts3$GEN190822_E_822446_lib93626_1689_6, cts3$GEN181107_S_156396_F_lib70375_1337_3, cts3$GEN181107_S_660410_lib70376_1337_3, cts3$GEN181107_S_660411_lib70426_1337_7, cts3$GEN190822_E_812440_lib93610_1689_5, cts3$GEN190822_E_825116_lib93600_1689_4, cts3$GEN190822_E_825216_lib93562_1689_1, cts3$GEN181107_S_AB24_M_lib70398_1337_5)
dim(nlf9)
nlf9ps=nlf9
batchn9=factor(c(1,1,2,2,2,1,1,1,2,2,2,1))
gendernf9m=factor(c(rep("f", 6), rep("m", 6)))
psn9= as.numeric(c("2.25", "2.5", "1.5", "1", "1.5", "2.25", "1.5",	"1.5",	"2.25",	"1.67",	"1.5",	"1"))
row.names(nlf9)=row.names(cts3)
colnames(nlf9)=c("nlf91", "nlf92", "nlf93", "nlf94", "nlf95", "nlf96", "nlf97", "nlf98", "nlf99", "nlf910", "nlf911", "nlf912")

##14-months
nlf14=cbind(cts3$GEN190822_E_587283_lib93638_1689_7, cts3$GEN190822_E_587284_lib93596_1689_4, cts3$`NG-24492_781945_lib386562_6738_2`, cts3$`NG-24492_781946_lib386563_6738_2`, cts3$`NG-24492_864376_lib386564_6738_2`, cts3$GEN190822_E_A70_1_lib93636_1689_7, cts3$GEN190822_E_545852_lib93590_1689_3, cts3$GEN190822_E_635278_lib93609_1689_5, cts3$GEN190822_E_639497_lib93614_1689_5, cts3$`NG-24492_781944_lib386557_6738_1`, cts3$GEN181107_S_AB15M_lib70416_1337_6, cts3$GEN181107_S_AB20M_lib70394_1337_4)
nlf14ps=nlf14[,-4]
batchn14=factor(c(2,2,3,3,3,2,2,2,2,3,1,1))
batchn14ps=factor(c(2,2,3,3,2,2,2,2,3,1,1))
gendernf14m=factor(c(rep("f", 6), rep("m", 6)))
dim(nlf14)
##for plaque score remove sample ##4 as plaque score is NA
psn14=as.numeric(c("3.25",	"3.5",	"2","2.17","3.5","3",	"3",	"3.25",	"2.67",	"2.5","1"))
row.names(nlf14)=row.names(cts3)
colnames(nlf14)=c("nlf141", "nlf142", "nlf143", "nlf144", "nlf145", "nlf146", "nlf147", "nlf148", "nlf149", "nlf1410", "nlf1411", "nlf1412")

##18 months
nlf18=cbind(cts3$GEN190822_E_481242_lib93613_1689_5, cts3$GEN190822_E_603744_lib93639_1689_7, cts3$GEN190822_E_662539_lib93570_1689_2, cts3$`NG-24492_691776_lib386565_6738_2`, cts3$`NG-27329_952124_lib478286_7389_2`, cts3$GEN190822_E_551590_lib93594_1689_4, cts3$`GEN181107_S_A49-L_lib70362_1337_2`, cts3$GEN181107_S_A50_1_lib70385_1337_3, cts3$GEN181107_S_A50_2_lib70440_1337_8, cts3$GEN181107_S_A50_3_lib70395_1337_4, cts3$GEN190822_E_A54_3_lib93560_1689_1 )
batchnlf18=factor(c(2,2,2,3,4,2,1,1,1,1,2))
##plaque score samples
nlf18ps=nlf18[,-11]
nlf18ps=nlf18ps[,-4]
nlf18ps=nlf18ps[,-1]
batchnlf18ps=factor(c(2,2,4,2,1,1,1,1))
gendernf18m=factor(c(rep("f", 5), rep("m", 6)))
dim(nlf18)
##for plaque score remove samples 1 and 4 and 11 (the last one)
psn18=as.numeric(c("NA","3.17","3.5","NA", "3", "3.75","3.17","3.67","3.83","3.67", "NA"))
row.names(nlf18)=row.names(cts3)
colnames(nlf18)=c("nlf181", "nlf182", "nlf183", "nlf184", "nlf185", "nlf186", "nlf187", "nlf188", "nlf189", "nlf1810", "nlf1811")

##24 months
nlf24=cbind(cts3$`NG-27329_814334_lib478284_7389_2`, cts3$`NG-27329_822445_lib478285_7389_2`, cts3$`NG-27329_836314_lib478287_7389_2`, cts3$`NG-27329_836315_lib478288_7389_2`, cts3$`NG-24492_408248_lib386558_6738_1`, cts3$GEN190822_E_417161_lib93582_1689_3, cts3$`NG-24492_423942_lib386559_6738_1`, cts3$GEN190822_E_423945_lib93653_1689_8, cts3$GEN190822_E_426651_lib93593_1689_3)
nlf24ps=nlf24
batchn24=factor(c(4,4,4,4,3,2,3,2,2))
gendernf24m=factor(c(rep("f", 4), rep("m", 5)))
psn24=as.numeric(c( "3.5",	"3.5",	"4",	"3.75","3.67",	"4.33",	"3",	"4",	"3.67"))
row.names(nlf24)=row.names(cts3)
colnames(nlf24)=c("nlf241", "nlf242", "nlf243", "nlf244", "nlf245", "nlf246", "nlf247", "nlf248", "nlf249")
dim(nlf24)
allnlf=cbind(nlf4m,nlf9,nlf14,nlf18,nlf24)
dim(allnlf)
batchallnlf=factor(c(1,2,2,2,3,1,1,1,1,2,2,2,1,1,2,2,2,1,1,1,2,2,2,1,2,2,3,3,3,2,2,2,2,3,1,1,2,2,2,3,4,2,1,1,1,1,2,4,4,4,4,3,2,3,2,2))

nlfage=cbind(nlf9,nlf14, nlf18, nlf24)
dim(nlfage)
nlfageps=cbind(nlf9ps, nlf14ps, nlf18ps, nlf24ps)
nlfageps2=round(nlfageps)
nlfage3=cbind(nlf9, nlf18, nlf24)
dim(nlfage3)
nlfage3=round(nlfage3)
nlfagepsf=as.numeric(c("2.25", "2.5", "1.5", "1", "1.5", "2.25", "1.5",	"1.5",	"2.25",	"1.67",	"1.5",	"1", "3.25",	"3.5",	"2","2.17", "3.5", "3",	"3",	"3.25",	"2.67",	"2.5","1", "3.17",	"3.5", "3", "3.75",	"3.17",	"3.67",	"3.83",	"3.67", "3.5",	"3.5",	"4",	"3.75","3.67",	"4.33",	"3",	"4",	"3.67"))
batchnlfage=factor(c(batchn9, batchn14, batchnlf18, batchn24))
batchnlfage=factor(c(1,1,2,2,2,1,1,1,2,2,2,1,2,2,3,3,3,2,2,2,2,3,1,1,2,2,2,3,4,2,1,1,1,1,2,4,4,4,4,3,2,3,2,2))
batchnlfage3=factor(c(batchn9, batchnlf18, batchn24))
batchnlfageps=factor(c(1,1,2,2,2,1,1,1,2,2,2,1,2,2,3,3,2,2,2,2,3,1,1,2,2,4,2,1,1,1,1,4,4,4,4,3,2,3,2,2 ))
gendernf=factor(c(gendernf9m, gendernf14m, gendernf18m, gendernf24m))
gendernf3= factor(c(gendernf9m, gendernf18m, gendernf24m))

nlfage3=cbind(nlf9, nlf18, nlf24)
dim(nlfage3)
nlfage3=round(nlfage3)


##nlgf

nlgf2=cbind(cts3$GEN190822_E_892669_lib93651_1689_8, cts3$GEN190822_E_892670_lib93565_1689_1, cts3$GEN190822_E_931176_lib93637_1689_7, cts3$GEN181107_S_191332_1_lib70445_1337_8, cts3$GEN181107_S_197129_3_lib70368_1337_2, cts3$GEN190822_E_892665_lib93643_1689_8, cts3$GEN190822_E_892666_lib93607_1689_5, cts3$GEN190822_E_897331_lib93606_1689_5, cts3$GEN190822_E_931175_lib93571_1689_2, cts3$GEN181107_S_191329_1_lib70410_1337_6, cts3$GEN181107_S_191330_1_lib70430_1337_7)
batchnlgf2=factor(c(2,2,2,1,1,2,2,2,2,1,1))
dim(nlgf2)
colnames(nlgf2)=c("nlgf2.1", "nlgf2.2", "nlgf2.3", "nlgf2.4", "nlgf2.5", "nlgf2.6", "nlgf2.7", "nlgf2.8", "nlgf2.9", "nlgf2.10", "nlgf2.11")
psnlgf2=as.numerical(c("2.83",	"2.33",	"2.83",	"2.17",	"4", "2",	"3.17",	"1.67",	"2.83",	"2.5",	"2.17"))
dim(nlgf2)

##4month
dim(cts3)
nlgf4=cbind(cts3$GEN190822_E_846025_lib93622_1689_6, cts3$GEN190822_E_846026_lib93608_1689_5, cts3$GEN190822_E_846027_lib93558_1689_1, cts3$GEN181107_S_183647_3_lib70418_1337_6 , cts3$GEN181107_S_186319_2_lib70396_1337_4 , cts3$GEN181107_S_186319_3_lib70442_1337_8, cts3$GEN181107_S_745232_lib70404_1337_5, cts3$`NG-24492_805516_lib386566_6738_2`, cts3$GEN190822_E_843819_lib93633_1689_7, cts3$GEN190822_E_896314_lib93629_1689_6, cts3$GEN181107_S_N57_2_lib70428_1337_7, cts3$GEN181107_S_N57_3_lib70369_1337_2) 
batchnlgf4=factor(c(2,2,2,1,1,1,1,3,2,2,1,1))
dim(nlgf4)
colnames(nlgf4)=c("nlgf4.1", "nlgf4.2", "nlgf4.3", "nlgf4.4", "nlgf4.5", "nlgf4.6", "nlgf4.7", "nlgf4.8", "nlgf4.9", "nlgf4.10", "nlgf4.11", "nlgf4.12")
nlgf4ps=as.numeric(c("4.5",	"4.75",	"4.5",	"4.5", "5",	"4.75", "4",	"4",	"4.5",	"4.25",	"4.5",	"4.5"))

nlgf9=cbind(cts3$GEN190822_E_746848_lib93630_1689_7, cts3$GEN190822_E_803956_lib93564_1689_1, cts3$GEN190822_E_803957_lib93621_1689_6, cts3$GEN190822_E_803958_lib93576_1689_2, cts3$GEN190822_E_805530_lib93623_1689_6, cts3$GEN190822_E_825465_lib93627_1689_6, cts3$GEN190822_E_723920_lib93635_1689_7, cts3$GEN190822_E_723921_lib93615_1689_5, cts3$GEN190822_E_825290_lib93574_1689_2, cts3$GEN181107_S_N35_1_lib70401_1337_5, cts3$`GEN181107_S_N35-2_lib70372_1337_2`, cts3$GEN181107_S_N37_1_lib70403_1337_5)
batchng9=factor(c(2,2,2,2,2,2,2,2,2,1,1,1))
dim(nlgf9)
colnames(nlgf9)=c("nlgf91", "nlgf92", "nlgf93", "nlgf94", "nlgf95", "nlgf96", "nlgf97", "nlgf98", "nlgf99", "nlgf910", "nlgf911", "nlgf912")
genderng9m=factor(c(rep("f", 6), rep("m", 6)))
nlgf9ps= as.numeric(c("5.75",	"5",	"5.25",	"5", "5.5",	"5", "5.5", "5.75",	"5",	"5",	"4.75",	"5"))

nlgf18=cbind( cts3$GEN190822_E_492037_lib93631_1689_7, cts3$GEN190822_E_571489_lib93649_1689_8, cts3$GEN190822_E_589853_lib93620_1689_6, cts3$GEN190822_E_590567_lib93632_1689_7, cts3$`NG-24492_805531_lib386567_6738_2`, cts3$GEN190822_E_496030_lib93647_1689_8, cts3$GEN190822_E_508442_lib93641_1689_7, cts3$GEN190822_E_509979_lib93589_1689_3, cts3$GEN190822_E_524732_lib93652_1689_8, cts3$GEN190822_E_540627_lib93595_1689_4, cts3$GEN190822_E_551246_lib93646_1689_8) 
batchnlgf18=factor(c(2,2,2,2,3,2,2,2,2,2,2))
dim(nlgf18)
colnames(nlgf18)=c("nlgf181", "nlgf182", "nlgf183", "nlgf184", "nlgf185", "nlgf186", "nlgf187", "nlgf188", "nlgf189", "nlgf1810", "nlgf1811")
genderng18m=factor(c(rep("f", 5), rep("m", 6)))
nlgf18ps=as.numeric(c("4.83",	"5.67",	"5.5",	"5.83",	"5.67", "5.5",	"5.83",	"5.33",	"5.67",	"5.67",	"5.17"))

nlgf24=cbind(cts3$`NG-24492_509981_lib386573_6738_2`, cts3$`NG-24492_511780_lib386572_6738_2`, cts3$`NG-24492_511781_lib386574_6738_3`, cts3$`NG-24492_527103_lib386569_6738_2`, cts3$`NG-24492_527104_lib386570_6738_2`, cts3$`NG-24492_527105_lib386571_6738_2`)         
batchnlgf24=factor(c(3,3,3,3,3,3))
dim(nlgf24)
colnames(nlgf24)=c("nlgf241", "nlgf242", "nlgf243", "nlgf244", "nlgf245", "nlgf246")
genderng24m=factor(c( rep("m", 5)))
nlgf24ps=as.numeric(c("5.67",	"5.67",	"5.83",	"5.75",	"5.75",	"5.25"))
## 2,2,2,2,2,2,2,2,2,1,1,1,2,2,2,2,3,2,2,2,2,2,2, 3,3,3,3,3,3
##5.75",	"5",	"5.25",	"5", "5.5",	"5", "5.5", "5.75",	"5",	"5",	"4.75",	"5", "4.83",	"5.67",	"5.5",	"5.83",	"5.67", "5.5",	"5.83",	"5.33",	"5.67",	"5.67",	"5.17", "5.67",	"5.67",	"5.83",	"5.75",	"5.75",	"5.25"
nlgfage=cbind(nlgf9, nlgf18, nlgf24)
nlgfage3=cbind(nlgf9, nlgf18)
dim(nlgfage)
row.names(nlgfage)=row.names(cts3)
nlgfageps=as.numeric(c("5.75",	"5",	"5.25",	"5", "5.5",	"5", "5.5", "5.75",	"5",	"5",	"4.75",	"5", "4.83",	"5.67",	"5.5",	"5.83",	"5.67", "5.5",	"5.83",	"5.33",	"5.67",	"5.67",	"5.17", "5.67",	"5.67",	"5.83",	"5.75",	"5.75",	"5.25"))
nlgfageps2=scale(nlgfageps, center = TRUE)
row.names(nlgfage)=row.names(cts3)
batchnlgfage=factor(c(batchng9, batchnlgf18, batchnlgf24 ))
batchnlgfage=factor(c(2,2,2,2,2,2,2,2,2,1,1,1,2,2,2,2,3,2,2,2,2,2,2,3,3,3,3,3,3))
gendernlgf= factor(c(genderng9m, genderng18m))
dim(nlgfage)

allnlgf=cbind(nlgf2, nlgf4, nlgf9, nlgf18, nlgf24)
dim(allnlgf)
#htau

htau9= cbind(cts3$`NG-24492_806526_lib386622_6738_6`, cts3$`NG-24492_812244_lib386623_6738_6`, cts3$`NG-24492_831427_lib386624_6738_6`, cts3$`NG-24492_835328_lib386634_6738_7`, cts3$`NG-24492_845420_lib386625_6738_6`, cts3$GEN190822_E_743931_lib93585_1689_3, cts3$GEN190822_E_831424_lib93619_1689_6, cts3$`NG-24492_835330_lib386620_6738_6`, cts3$`NG-24492_856365_lib386619_6738_6`, cts3$`NG-24492_882506_lib386621_6738_6`)
dim(htau9)
batchhtau9= factor(c(3,3,3,3,3,2,2,3,3,3))
genderhtau9=factor(c(rep("f", 5), rep("m", 5)))
rownames(htau9)=rownames(cts3)
colnames(htau9)=c("htau9.1", "htau9.2", "htau9.3", "htau9.4", "htau9.5", "htau9.6", "htau9.7", "htau9.8", "htau9.9", "htau9.10")

htau18=cbind(cts3$`NG-24492_645217_lib386636_6738_7`, cts3$`NG-24492_728833_lib386637_6738_8`, cts3$`NG-24492_744283_lib386633_6738_7`, cts3$`NG-24492_681466_lib386629_6738_7`, cts3$`NG-24492_713817_lib386627_6738_7`, cts3$`NG-24492_719582_lib386628_6738_7`, cts3$`NG-24492_726801_lib386626_6738_7`, cts3$`NG-27329_871902_lib478279_7389_2`)
batchhtau18=factor(c(3,3,3,3,3,3,3,4))
dim(htau18)
genderhtau18= factor(c(rep("f", 3), rep("m", 5)))
rownames(htau18)=rownames(cts3)
colnames(htau18)=c("htau18.1", "htau18.2", "htau18.3", "htau18.4", "htau18.5", "htau18.6", "htau18.7", "htau18.8")


htau24=cbind(cts3$`NG-24492_645218_lib386638_6738_8`, cts3$`NG-27329_719741_lib478283_7389_2`, cts3$`NG-27329_731934_lib478282_7389_2`, cts3$`NG-24492_599928_lib386630_6738_7`, cts3$`NG-24492_599929_lib386631_6738_7`, cts3$`NG-24492_646726_lib386632_6738_7`, cts3$`NG-27329_835319_lib478281_7389_2`, cts3$`NG-27329_812245_lib478280_7389_2`)
batchhtau24=factor(c(4,4,4,4,4,4,4,4))
dim(htau24)
genderhtau24= factor(c(rep("f", 3), rep("m", 5)))
rownames(htau24)=rownames(cts3)
colnames(htau24)=c("htau24.1", "htau24.2", "htau24.3", "htau24.4", "htau24.5", "htau24.6", "htau24.7", "htau24.8")

htauage=cbind(htau9, htau18, htau24)
batchallhtau=factor(c(3,3,3,3,3,2,2,3,3,3,3,3,3,3,3,3,3,4,4,4,4,4,4,4,4,4))
dim(htauage)

##nlfhtau

nfht4= cbind(cts3$GEN181107_S_719739_lib70424_1337_7, cts3$GEN181107_S_719766_lib70425_1337_7, cts3$GEN181107_S_730297_lib70407_1337_5, cts3$`NG-24492_860287_lib386590_6738_4`, cts3$`NG-24492_865160_lib386591_6738_4`, cts3$`NG-24492_922226_lib386589_6738_4`, cts3$GEN181107_S_713816_lib70441_1337_8, cts3$GEN181107_S_719581_lib70435_1337_8, cts3$GEN181107_S_719756_lib70420_1337_6, cts3$`NG-24492_860293_lib386576_6738_3`, cts3$`NG-24492_865155_lib386577_6738_3`, cts3$`NG-24492_865158_lib386575_6738_3`)
batchhtaunfht4=factor(c(1,1,1,3,3,3,1,1,1,3,3,3))
dim(nfht4)
gendernfht4= factor(c(rep("f", 6), rep("m", 6)))
rownames(nfht4)=rownames(cts3)
colnames(nfht4)=c("nfht4.1", "nfht4.2", "nfht4.3", "nfht4.4", "nfht4.5", "nfht4.6", "nfht4.7", "nfht4.8","nfht4.9", "nfht4.10", "nfht4.11", "nfht4.12")


nfht9=cbind(cts3$`NG-24492_938735_lib386592_6738_4`, cts3$`NG-24492_922222_lib386578_6738_3`, cts3$`NG-24492_922223_lib386579_6738_3`, cts3$`NG-24492_937261_lib386582_6738_3`, cts3$`NG-24492_938733_lib386580_6738_3`, cts3$`NG-24492_938982_lib386581_6738_3`, cts3$`NG-24492_951257_lib386583_6738_3`)
batchhtaunfht9=factor(c(3,3,3,3,3,3,3))
dim(nfht9)
gendernfht9= factor(c(rep("f", 1), rep("m", 6)))
rownames(nfht9)=rownames(cts3)
colnames(nfht9)=c("nfht9.1", "nfht9.2", "nfht9.3", "nfht9.4", "nfht9.5", "nfht9.6", "nfht9.7")
nfht9ps=as.numeric(c("1.5", "1.166666667",	"1.5",	"1.5",	"1.5", "1.333333333", 	"1.5"))

##"1.5", "1.166666667",	"1.5",	"1.5",	"1.5", "1.333333333", 	"1.5", "4",	"3",	"3.333333333", "3.833333333",	"2.333333333",	"3.5",	"3.333333333",	"2.5",	"4", "NA",	"3.25", "3.833333333", "3.666666667"

nfht18=cbind(cts3$`NG-24492_810235_lib386594_6738_4`, cts3$`NG-27329_938734_lib478294_7389_2`, cts3$`NG-27329_951260_lib478295_7389_2`, cts3$`NG-24492_732585_lib386584_6738_3`, cts3$`NG-24492_744454_lib386585_6738_4`, cts3$`NG-24492_744457_lib386586_6738_4`, cts3$`NG-27329_937263_lib478292_7389_2`, cts3$`NG-27329_937264_lib478293_7389_2`, cts3$`NG-27329_951258_lib478296_7389_2`)
batchnfht18=factor(c(3,4,4,3,3,3,4,4,4))
dim(nfht18)
gendernfht18= factor(c(rep("f", 3), rep("m", 6)))
rownames(nfht18)=rownames(cts3)
colnames(nfht18)=c("nfht18.1", "nfht18.2", "nfht18.3", "nfht18.4", "nfht18.5", "nfht18.6", "nfht18.7", "nfht18.8", "nfht18.9")
nfht18ps=as.numeric(c("4",	"3",	"3.333333333", "3.833333333",	"2.333333333",	"3.5",	"3.333333333",	"2.5",	"4"))

nfht24=cbind(cts3$`NG-27329_732586_lib478289_7389_2`, cts3$`NG-27329_841870_lib478300_7389_2`, cts3$`NG-27329_744453_lib478290_7389_2`, cts3$`NG-27329_744455_lib478291_7389_2`)
batchhtaunfht24=factor(c(4,4,4,4))
dim(nfht24)
gendernfht24= factor(c(rep("f", 2), rep("m", 2)))
rownames(nfht24)=rownames(cts3)
colnames(nfht24)=c("nfht24.1", "nfht24.2", "nfht24.3", "nfht24.4")
nfht24ps=as.numeric(c("NA",	"3.25", "3.833333333", "3.666666667"))
nfht=cbind(nfht9, nfht18, nfht24)
allnfht=cbind(nfht4,nfht9,nfht18,nfht24)
batchallnfht=factor(c(1,1,1,3,3,3,1,1,1,3,3,3,3,3,3,3,3,3,3,3,4,4,3,3,3,4,4,4,4,4,4,4))
b=as.matrix(batchallnfht)
dim(allnfht)
nfhtage=cbind(nfht9, nfht18, nfht24)
dim(nfhtage)
nfhtage2=round(nfhtage)
dim(nfhtage)
batchnhtage=factor(c(3,3,3,3,3,3,3,3,4,4,3,3,3,4,4,4,4,4,4,4))


##NLGFHtau
nght4=cbind(cts3$`NG-24492_856399_lib386605_6738_5`, cts3$`NG-24492_856453_lib386606_6738_5`, cts3$`NG-24492_961589_lib386602_6738_5`, cts3$`NG-24492_961593_lib386601_6738_5`, cts3$`NG-24492_979476_lib386603_6738_5`, cts3$`NG-24492_979488_lib386604_6738_5`, cts3$`NG-24492_856451_lib386600_6738_5`, cts3$`NG-24492_856570_lib386599_6738_5`, cts3$`NG-24492_939009_lib386598_6738_4`, cts3$`NG-24492_942017_lib386597_6738_4`, cts3$`NG-24492_979474_lib386595_6738_4`, cts3$`NG-24492_979482_lib386596_6738_4`)
batchnght4=factor(c(3,3,3,3,3,3,3,3,3,3,3,3))
dim(nght4)
gendernght4= factor(c(rep("f", 6), rep("m", 6)))
rownames(nght4)=rownames(cts3)
colnames(nght4)=c("nght4.1", "nght4.2", "nght4.3", "nght4.4", "nght4.5", "nght4.6", "nght4.7", "nght4.8", "nght4.9", "nght4.10", "nght4.11", "nght4.12")


nght9=cbind(cts3$`NG-24492_942025_lib386612_6738_6`, cts3$`NG-24492_945527_lib386614_6738_6`, cts3$`NG-24492_937267_lib386607_6738_5`, cts3$`NG-24492_937269_lib386608_6738_5`, cts3$`NG-24492_942020_lib386609_6738_5`, cts3$`NG-24492_945524_lib386610_6738_6`, cts3$`NG-24492_945535_lib386611_6738_6`)
batchnght9=factor(c(3,3,3,3,3,3,3))
dim(nght9)
gendernght9= factor(c(rep("f", 2), rep("m", 5)))
rownames(nght9)=rownames(cts3)
colnames(nght9)=c("nght9.1", "nght9.2", "nght9.3", "nght9.4", "nght9.5", "nght9.6", "nght9.7")
nght9ps=as.numeric(c("5.166667",	"4.833333333", "5",	"5.166666667",	"5.5",	"5",	"4.66"))

nght18=cbind(cts3$`NG-27329_937271_lib478302_7389_2`, cts3$`NG-27329_939014_lib478303_7389_2`, cts3$`NG-27329_942024_lib478305_7389_2`, cts3$`NG-27329_897330_lib478306_7389_2`, cts3$`NG-27329_918830_lib478301_7389_2`, cts3$`NG-27329_937270_lib478304_7389_2`)
batchnght18=factor(c(4,4,4,4,4,4))
dim(nght18)
gendernght18= factor(c(rep("f", 3), rep("m", 3)))
rownames(nght18)=rownames(cts3)
colnames(nght18)=c("nght18.1", "nght18.2", "nght18.3", "nght18.4", "nght18.5", "nght18.6")
nght18ps=as.numeric(c("5.5",	"6",	"5.833333333", "6",	"5.833333333",	"6"))
## "5.166667",	"4.833333333", "5",	"5.166666667",	"5.5",	"5",	"4.66", "5.5",	"6",	"5.833333333", "6",	"5.833333333",	"6"
##all nght
nghtage= cbind(nght4, nght9, nght18)
nghtage= cbind(nght9, nght18)
dim(nghtage)
batchnght=factor(c(3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,4,4,4,4,4,4))

