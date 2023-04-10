*草稿测试
使用plink运行

```sh
#!/bin/bash
#CSUB -J PCA14
#CSUB -o %J.out 
#CSUB -e %J.error
#CSUB -q fat4
#CSUB -n 24
#CSUB -R "span[hosts=1]"
source /share/apps/anaconda3/bin/activate wgs

#输入多月份表型数据文件
phe_file=S6.tab
#输入基因型文件前綴位置
snp_prefix=/share/home/hnu_wangy/project/wgs/297sample/popfilter/bcftools/snp4/310/snp
#主成分个数
pca_number=14

#创建目录并生成每个月份的表型文件
sed -i 's/NA/\-9/g' $phe_file
mkdir -p ${phe_file%.*}
NF=$(head -1 $phe_file |awk '{print NF}')
for i in $(seq 2 $NF); do cut -f $i $phe_file > tmp; id=$(head -1 tmp); sed -i '1d' tmp; paste $snp_prefix.nosex tmp > ${phe_file%.*}/$id.phe; rm tmp; done

#计算生成分构建协变量
mkdir -p pca$pca_number;cd pca$pca_number
plink --bfile $snp_prefix --pca $pca_number --out pca
paste pca.nosex <(cut -d' ' -f3- pca.eigenvec) > cov.txt

#linear关联
for i in `ls ../${phe_file%.*}`
do
	plink --bfile $snp_prefix --pheno ../${phe_file%.*}/$i --linear --covar cov.txt --out ${i%.*} --allow-no-sex --hide-covar 
	awk '!/'NA'/' ${i%.*}.assoc.linear > ${i%.*}.linear.result
done




#CSUB -J plot4
#CSUB -o %J.out 
#CSUB -e %J.error
#CSUB -q fat4
#CSUB -n 1
#CSUB -R "span[hosts=1]"
source /share/apps/anaconda3/bin/activate rstudio

for i in `ls *linear.result`;do Rscript linear_cmplot.R ${i%%.*};done
mkdir figure
mv *jpg figure
```





```R
#linear_plot.R
library(qqman)
args <- commandArgs(trailingOnly=TRUE)
i <- as.character(args[1])

result <- read.table(paste(i,"linear.result",sep = "."), head=TRUE)

cat("START",i,"\n")
jpeg(paste(i,"manhattan.jpeg",sep = "_"))
manhattan(result,chr="CHR",bp="BP",p="P",snp="SNP", main = i)
dev.off()
    
jpeg(paste(i,"QQ-Plot.jpeg",sep = "_"))
qq(result$P, main = i)
dev.off()
cat("END",i,"\n")



#linear_cmplot.R
library(CMplot)
library(data.table)
args <- commandArgs(trailingOnly=TRUE)
i <- as.character(args[1])

dat<-fread(paste(i,"linear.result",sep = "."),head=T)
dat2<-data.frame(dat$SNP,dat$CHR,dat$BP,dat$P)
CMplot(dat2,plot.type = "q",memo = i,main = "")
CMplot(dat2,plot.type = "m",threshold = 1e-5,threshold.col="red",
       threshold.lty = 2,threshold.lwd = 1,memo = i)

```



```sh
#!/bin/bash
#CSUB -J plink_S8D10
#CSUB -o %J.out 
#CSUB -e %J.error
#CSUB -q fat4
#CSUB -n 24
#CSUB -R "span[hosts=1]"
source /share/apps/anaconda3/bin/activate wgs

#输入多月份表型数据文件
phe_file=S8D10.tab
#输入基因型文件前綴位置
snp_prefix=/share/home/hnu_wangy/project/wgs/297sample/popfilter/bcftools/snp4/310/snp
#主成分个数
pca_number=4

#创建目录并生成每个月份的表型文件
sed -i 's/NA/\-9/g' $phe_file
mkdir -p ${phe_file%.*}
NF=$(head -1 $phe_file |awk '{print NF}')
for i in $(seq 2 $NF); do cut -f $i $phe_file > tmp; id=$(head -1 tmp); sed -i '1d' tmp; paste $snp_prefix.nosex tmp > ${phe_file%.*}/$id.phe; rm tmp; done

#计算生成分构建协变量
plink --bfile $snp_prefix --pca $pca_number --out pca
paste pca.nosex <(cut -d' ' -f3- pca.eigenvec) > cov.txt

#linear关联
for i in `ls ${phe_file%.*}`
do
	plink --bfile $snp_prefix --pheno ${phe_file%.*}/$i --linear --covar cov.txt --out ${i%.*} --allow-no-sex --hide-covar 
	awk '!/'NA'/' ${i%.*}.assoc.linear > ${i%.*}.linear.result
done

#画图
source /share/apps/anaconda3/bin/activate rstudio
for i in `ls ${phe_file%.*}`;do Rscript linear_plot.R ${i%.*};done
```









GEMMA

```sh
#!/bin/bash
#CSUB -J test
#CSUB -o %J.out 
#CSUB -e %J.error
#CSUB -q fat8
#CSUB -n 48
#CSUB -R "span[hosts=1]"
source /share/apps/anaconda3/bin/activate wgs

#输入多月份表型数据文件
phe_file=test.tab
#输入基因型文件前綴位置
snp_prefix=/share/home/hnu_wangy/project/wgs/297sample/popfilter/bcftools/snp4/310/snp
#主成分个数
pca_number=4

#创建目录并生成每个月份的表型文件
mkdir -p ${phe_file%.*};cd ${phe_file%.*}
NF=$(head -1 ../$phe_file |awk '{print NF}')
for i in $(seq 2 $NF); do cut -f $i ../$phe_file > tmp; id=$(head -1 tmp); sed -i '1d' tmp; mv tmp $id.phe; done
cd ..

#计算主成分构建协变量
mkdir -p pca$pca_number;cd pca$pca_number
plink --bfile $snp_prefix --pca $pca_number --out pca
cut -f 3- -d ' ' pca.eigenvec |awk '{print 1,$0}' > cov.txt

#生成K矩阵
gemma -bfile $snp_prefix -gk 2 -p ../${phe_file%.*}/*10*.phe -o kin

#linear关联
for i in `ls ../${phe_file%.*}`
do
	echo "###START ${i%.*}###"
	gemma -bfile $snp_prefix -k output/kin.sXX.txt -lmm 1 -p ../${phe_file%.*}/${i%.*}.phe -c cov.txt -o ${i%.*}
	cat <(head -1 output/${i%.*}.assoc.txt) <(sed '1d' output/${i%.*}.assoc.txt |sed 's/chr//') > output/${i%.*}.assoc.result
	echo "###END ${i%.*}###"
done

#画图
source /share/apps/anaconda3/bin/activate rstudio
for i in `ls ../${phe_file%.*}`;do Rscript ../gemma_plot.R output/${i%.*};done



#只用一个PCA
#GROWD5.assoc.txt
chr     rs      ps      n_miss  allele1 allele0 af      beta    se      logl_H1 l_remle p_wald
chr1       SNP2    8352    11      T       C       0.410   9.945367e-01    5.709739e-01    -9.866381e+02   1.000000e-05    8.261041e-02
chr1       SNP3    8401    3       C       T       0.233   1.787456e+00    7.378280e-01    -9.852422e+02   1.000000e-05    1.602936e-02
chr1       SNP8    28079   2       T       C       0.113   -1.739455e-01   8.683910e-01    -9.881271e+02   1.000000e-05    8.413817e-01


```





```R
#gemma_plot.R
library(qqman)
args <- commandArgs(trailingOnly=TRUE)
i <- as.character(args[1])

result <- read.table(paste(i,"assoc.result",sep = "."), head=TRUE)

cat("START",i,"\n")
jpeg(paste(i,"manhattan.jpeg",sep = "_"))
manhattan(result,chr="chr",bp="ps",p="p_wald",snp="rs", main = i)
dev.off()
    
jpeg(paste(i,"QQ-Plot.jpeg",sep = "_"))
qq(result$p_wald, main = i)
dev.off()
cat("END",i,"\n")



p_value=result$p_wald
z = qnorm(p_value/ 2)
lambda = round(median(z^2, na.rm = TRUE) / 0.454, 3)
```







```sh
#!/bin/bash
#CSUB -J gemma_S8D10
#CSUB -o %J.out 
#CSUB -e %J.error
#CSUB -q cpu_cae
#CSUB -n 12
#CSUB -R "span[hosts=1]"
source /share/apps/anaconda3/bin/activate wgs

#输入多月份表型数据文件
phe_file=S8D10.tab
#输入基因型文件前綴位置
snp_prefix=/share/home/hnu_wangy/project/wgs/297sample/popfilter/bcftools/snp4/310/snp
#主成分个数
pca_number=4

#创建目录并生成每个月份的表型文件
mkdir -p ${phe_file%.*};cd ${phe_file%.*}
NF=$(head -1 ../$phe_file |awk '{print NF}')
for i in $(seq 2 $NF); do cut -f $i ../$phe_file > tmp; id=$(head -1 tmp); sed -i '1d' tmp; mv tmp $id.phe; done
cd ..

#计算主成分构建协变量
#plink --bfile $snp_prefix --pca $pca_number --out pca
#cut -f 3- -d ' ' pca.eigenvec |awk '{print 1,$0}' > cov.txt

#生成K矩阵
#gemma -bfile $snp_prefix -gk 2 -p ${phe_file%.*}/*9*.phe -o kin

#lmm关联
for i in `ls ${phe_file%.*}`
do
	echo "###START ${i%.*}###"
	gemma -bfile $snp_prefix -k output/kin.sXX.txt -lmm 1 -p ${phe_file%.*}/${i%.*}.phe -c cov.txt -o ${i%.*}
	cat <(head -1 output/${i%.*}.assoc.txt) <(sed '1d' output/${i%.*}.assoc.txt |sed 's/chr//') > output/${i%.*}.assoc.result
	echo "###END ${i%.*}###"
done

#画图
source /share/apps/anaconda3/bin/activate rstudio
mkdir fig
for i in `ls ${phe_file%.*}`;do Rscript gemma_plot.R output/${i%.*};done


#314
#CSUB -J plot2
#CSUB -o %J.out 
#CSUB -e %J.error
#CSUB -q fat4
#CSUB -n 1
#CSUB -R "span[hosts=1]"
source /share/apps/anaconda3/bin/activate rstudio
for i in `ls *assoc.result`;do Rscript gemma_cmplot.R ${i%%.*};done
mkdir figure
mv *jpg figure
```





```R
#gemma_plot.R
library(qqman)
args <- commandArgs(trailingOnly=TRUE)
i <- as.character(args[1])

result <- read.table(paste(i,"assoc.result",sep = "."), head=TRUE)

cat("START",i,"\n")
jpeg(paste(i,"manhattan.jpeg",sep = "_"))
manhattan(result,chr="chr",bp="ps",p="p_wald",snp="rs", main = i)
dev.off()
    
jpeg(paste(i,"QQ-Plot.jpeg",sep = "_"))
qq(result$p_wald, main = i)
dev.off()
cat("END",i,"\n")



#CMplot
library(CMplot)
CMplot(dat,plot.type = "q",memo = "test",main = "")
CMplot(dat,plot.type = "m",threshold = 1e-5,threshold.col="red",
       threshold.lty = 2,threshold.lwd = 1,memo = "test")


#gemma_cmplot.R
library(CMplot)
library(data.table)
args <- commandArgs(trailingOnly=TRUE)
i <- as.character(args[1])

dat<-fread(paste(i,"assoc.result",sep = "."),head=T)
dat2<-data.frame(dat$rs,dat$chr,dat$ps,dat$p_wald)
CMplot(dat2,plot.type = "q",memo = i,main = "")
CMplot(dat2,plot.type = "m",threshold = 1e-5,threshold.col="red",
       threshold.lty = 2,threshold.lwd = 1,memo = i)


#test
H10.assoc.result
dat<-fread("H10.assoc.result",head=T)

dat2<-data.frame(dat$rs,dat$chr,dat$ps,dat$p_wald)
CMplot(dat2,plot.type = "q",main = "")
CMplot(dat2,plot.type = "m",threshold = 1e-5,threshold.col="red",
       threshold.lty = 2,threshold.lwd = 1)



```



GEMMA算法中使用的技巧依赖于每个SNP的完整或估算基因型数据。也就是说，GEMMA要求用户在关联测试之前输入所有缺失的基因型。这一插补步骤可以说比简单地删除缺失基因型的个体更好，因为它可以提高检测关联的能力。因此，为了拟合LMM或BSLMM，建议首先插补缺失基因型。否则，任何缺失超过某一阈值（默认值为5%）的SNP将不会被分析，而缺失的SNP基因型如果没有超过该阈值，将被简单地替换为该SNP的估计平均基因型。然而，对于预测，所有SNP都将被使用，无论它们是否缺失。测试集中缺失的基因型将被训练集中的平均基因型替换







```sh
#!/bin/bash
#CSUB -J test
#CSUB -o %J.out 
#CSUB -e %J.error
#CSUB -q cpu
#CSUB -n 1
#CSUB -R "span[hosts=1]"
source /share/apps/anaconda3/bin/activate wgs
plink --bfile ../snp --indep-pairwise 50 10 0.95 --out ld95


#LD过滤（可选，SNP数量太大还是过滤不然运行太慢了）
#注意：map文件中SNP需要有名字，vcf用plink转换这列为.，LD找出来的是SNP名字
SNP=../snp
#50 10 0.2  
#滑动窗口是50个SNP  
#每次滑动的大小是10个SNP  
#r2 0.1表示R方小于0.1；50 10 0.2 ；50 5 0.4
#在50个SNP滑动窗口（每次前进10个SNP）内去除R2值大于0.1的每个SNP和任何其他SNP。
plink --bfile $SNP --indep-pairwise 50 10 0.95 --out ld95
#将剩余的（未定向的）SNP复制到修剪数据。
plink --bfile $SNP --extract ld.prune.in --make-bed --out LDfiltered

#50 10 0.1(admixture mannual example)，但是过滤太多可以用50 10 0.2过滤少些




半衰距离 5603bp
```





## SNP显著位点筛选

GEMMA结果

```bash

#awk '$12<1e-5 {print $0}' GROWD5.assoc.result > GROWD5.marked.result

id=D10
mkdir ${id}_marked
for i in `ls *$id.assoc.result`;do cat <(head -1 $i) <(awk '$12<1e-5 {print $0}' $i) > ${i%%.*}.marked.result;done
mv *$id.marked.result ${id}_marked

#取并集
cd ${id}_marked;cat [0-9]*$id.marked.result| sort | uniq > union.txt



id=H
mkdir ${id}_marked
for i in `ls $id*.assoc.result`;do cat <(head -1 $i) <(awk '$12<1e-5 {print $0}' $i) > ${i%%.*}.marked.result;done
mv $id*.marked.result ${id}_marked



grep -Fx -f 

cat file1.txt file2.txt file3.txt ... fileN.txt | sort | uniq > union.txt
grep -Fx -f file1.txt file2.txt | \
  grep -Fx -f - file3.txt | \
  grep -Fx -f - file4.txt | \
  ... | \
  grep -Fx -f - fileN.txt > intersection.txt
  

#phe=D5
grep -Fx -f H01 H02 | grep -Fx -f - H03 | grep -Fx -f - H04 | \
grep -Fx -f - H05 | grep -Fx -f - H06 | grep -Fx -f - H07 | \
grep -Fx -f - H08 | grep -Fx -f - H09 | grep -Fx -f - H10 | \
grep -Fx -f - H11 > intersection.txt



suffix=assoc.result
prefix=D5

grep -Fx -f 01$prefix.$suffix 02$prefix.$suffix | \
grep -Fx -f - 03$prefix.$suffix | \
grep -Fx -f - 04$prefix.$suffix | \
grep -Fx -f - 05$prefix.$suffix | \
grep -Fx -f - 06$prefix.$suffix | \
grep -Fx -f - 07$prefix.$suffix | \
grep -Fx -f - 08$prefix.$suffix | \
grep -Fx -f - 09$prefix.$suffix | \
grep -Fx -f - 10$prefix.$suffix | \
grep -Fx -f - 11$prefix.$suffix > intersection.txt






```



plink

```sh
id=D10
mkdir ${id}_marked
for i in `ls *$id.linear.result`;do cat <(head -1 $i) <(awk '$12<1e-5 {print $0}' $i) > ${i%%.*}.marked.result;done
mv *$id.marked.result ${id}_marked

#取并集
cd ${id}_marked; cat [0-9]*$id.marked.result| sort | uniq > union.txt


cat 01* 02*|sort |uniq -d




file=union.txt
#sed -i '$d' $file
sed 's/^/chr0&/g' $file |sed 's/chr010/chr10/' |awk '{OFS="\t"}{print $1,$3,$3,$2}' > ${file%.*}.snp_info.bed
bedtools intersect -a ${file%.*}.snp_info.bed -b gene_info.bed  -wa -wb > snp_gene.txt
cut -f 8 snp_gene.txt |sort -u > snp_gene.id
LD=5603
sed 's/^/chr0&/g' $file |sed 's/chr010/chr10/' |awk -va=$LD '{OFS="\t"}{print $1,$3-a,$3+a,$2}' > ${file%.*}.snp_region.bed
bedtools intersect -a ${file%.*}.snp_region.bed -b gene_info.bed -wa -wb > snp_region_gene.txt
cut -f 8 snp_region_gene.txt |sort -u > snp_region_gene.id


file=grow.txt
sed 's/^/chr0&/g' $file |sed 's/chr010/chr10/' |awk '{OFS="\t"}{print $1,$3,$3,$2}' > ${file%.*}.snp_info.bed
bedtools intersect -a ${file%.*}.snp_info.bed -b gene_info.bed  -wa -wb > snp_gene.txt
cut -f 8 snp_gene.txt |sort -u > snp_gene.id
LD=5603
sed 's/^/chr0&/g' $file |sed 's/chr010/chr10/' |awk -va=$LD '{OFS="\t"}{print $1,$3-a,$3+a,$2}' > ${file%.*}.snp_region.bed
bedtools intersect -a ${file%.*}.snp_region.bed -b gene_info.bed -wa -wb > snp_region_gene.txt
cut -f 8 snp_region_gene.txt |sort -u > snp_region_gene.id
```







Post-GWAS: eQTL、mQTL共定位分析(SMR)

https://zhuanlan.zhihu.com/p/337969763







论文可参考



群体结构分析structure原理

https://www.omicsclass.com/article/1366

https://mp.weixin.qq.com/s/x8I9zBaYin6yRQF4oQLW4w

https://mp.weixin.qq.com/s/wHbEirrBFEnIUwGBOwMs5g





GWAS eQTL 共定位

https://mp.weixin.qq.com/s/oghsIJLe9rYEIHuW3FL5lg



RNA-seq进行eQTL分析的一系列问题

https://www.omicshare.com/forum/forum.php?mod=viewthread&tid=704&highlight=eQTL

组学大讲堂

https://www.omicsclass.com/search/all?word=GO%E5%AF%8C%E9%9B%86



```R


source("../GAPIT.library.R.txt")
source("../gapit_functions.txt")

myY <- read.table("../height.txt", head = T)
myG <- read.table("mdp_genotype_test.hmp.txt", head = FALSE)
myKI <- read.table("GAPIT.Genotype.Kin_Zhang.csv", head = FALSE)

#Step 2: Run GAPIT
myGAPIT <- GAPIT(
Y=myY[,1:2],
G=myG,
KI=myKI,
PCA.total=4,
file.G="../snpf4_chr",
file.Ext.G="hmp.txt",
file.from=1,
file.to=10
)

GAPIT.Genotype.Kin_Zhang.csv



gene=scaffold_142.207

cat <(head -1 tpm.tab) <(grep $gene tpm.tab ) > gene.tpm

#转置
awk '
{
  for (i=1; i<=NF; i++) {
    a[NR,i] = $i
  }
}
NF>p { p = NF }
END {
  for(j=1; j<=p; j++) {
    str=a[1,j]
    for(i=2; i<=NR; i++){
      str=str" "a[i,j];
    }
    print str
  }
}' gene.tpm |sed '1d' |sort |cut -f 2 -d ' '> tmp

paste test.nosex tmp > test.phe


plink --file test --pheno test.phe --assoc -out $gene --allow-no-sex 

plink --file test --pheno test.phe --linear -out out2  --allow-no-sex --covar cov.txt --hide-covar

```



我有一个输入文件共有十三列，现在需要使用R对前十二列中的每一列绘制频数分布图并添加正态曲线，将十二个图放在一张图中输出



https://mp.weixin.qq.com/s/CEjrIyunRL6ndABjBGUSYQ



```R
# 读取文件A
df_a <- read.csv("A.csv", header = TRUE)

# 使用 diff() 函数计算每一列相邻元素的差值
df_b <- data.frame(diff(df_a))

# 将文件B写入CSV文件
write.csv(df_b, "B.csv", row.names = FALSE)

```



```R
# 读取数据文件
setwd("D:/毕设数据/result/1. phe/毕设")
data <- read.table("GROW.txt", header=TRUE,row.names=1)

xname="Height(cm)"
xname="Diameter(mm)"
# 绘制频数分布图并添加正态曲线
par(mfrow=c(3, 4)) # 将图像分为3行4列
for(i in 1:12){
  x<-data[,i]
  x <- x[!is.na(x)]
  h<-hist(x,main=colnames(data)[i], xlab=xname,cex.lab=1.5,cex.sub=1.5,cex.main=2)
  xfit<-seq(min(x),max(x),length=100) ##生成从X的最小值到最大值的等间距的40个数
  yfit<-dnorm(xfit,mean=mean(x),sd=sd(x)) ##使用dnorm()函数生成服从正态分布的概率密度函数值
  yfit <- yfit*diff(h$mids[1:2])*length(x) ##在这里diff()函数是计算两数之差，也即直方图的组距；这一行是计算出模拟的Y值，为后续绘图做准备。
  lines(xfit, yfit, col="red", lwd=2)
}

# 保存图像并输出
dev.copy(jpeg, "output.jpg")
dev.off()


par(mfrow=c(1, 3)) 
x<-data[,3]
x <- x[!is.na(x)]
h<-hist(x,main=colnames(data)[i], xlab=xname,cex.lab=1.5,cex.sub=1.5,cex.main=2)
xfit<-seq(min(x),max(x),length=100) ##生成从X的最小值到最大值的等间距的40个数
yfit<-dnorm(xfit,mean=mean(x),sd=sd(x)) ##使用dnorm()函数生成服从正态分布的概率密度函数值
yfit <- yfit*diff(h$mids[1:2])*length(x) ##在这里diff()函数是计算两数之差，也即直方图的组距；这一行是计算出模拟的Y值，为后续绘图做准备。
lines(xfit, yfit, col="red", lwd=2)

现有一个文件A，共有十二列，要获得文件B，B的第一列等于A的第二列减去第一列，B的第二列等于A的第三列减去第二列，依此类推，共得到十一列，请用R实现，注意A与B中可能含有NA值
library(dplyr)

# 读取A.csv文件
A <- read.csv("A.csv")

# 对A进行处理，生成B
B <- A %>%
  mutate_at(vars(2:12), ~ ifelse(is.na(.), NA, .)) %>%  # 处理NA值
  mutate(across(2:12, ~ . - lag(.))) %>%
  select(-1)  # 去除第一列

# 将结果写入文件B.csv
write.csv(B, "B.csv", row.names = FALSE)

```





SnpEff安装使用及报错解决

https://zhuanlan.zhihu.com/p/476561285







```


https://mp.weixin.qq.com/s/RLsTaAumNi9VfX8CXjzxSA



id=pop1
less $id.stat.gz|sort -k 2,2nr |head -1



perl ~/bin/Plot_MultiPop.pl -inList multi.list --output re3 -bin1 10 -bin2 500



```



```R
BiocManager::install("CMplot")
library(CMplot)
setwd("yourpath")
mydata<-read.table("snp_density.txt",header=TRUE,sep="\t")
CMplot(mydata,plot.type="d",bin.size=1e6,col=c("darkgreen","yellow", "red"),file="jpg",memo="snp_density",dpi=300)




dat2<-data.frame(dat$rs,dat$chr,dat$ps,dat$p_wald)
CMplot(dat2,plot.type = "q",main = "")
CMplot(dat2,plot.type = "m",threshold = 1e-5,threshold.col="red",
       threshold.lty = 2,threshold.lwd = 1)

setwd("/share/home/hnu_wangy/project/wgs/297sample/popfilter/bcftools/snp4/310/anno")

dat<-read.table("../snp.map")

CMplot(dat,plot.type = "d",main = "")



CMplot(mm,plot.type="d",bin.size=1e6,col=c("darkgreen", "yellow", "red"),
       file="tiff",memo="",dpi=300,file.output=TRUE, verbose=TRUE)

CMplot(mm,plot.type="d",bin.size=1e6,col=c("darkgreen", "yellow", "red"),
       file="tiff",memo="",dpi=300,file.output=FALSE, verbose=TRUE)





library(CMplot)
dat = read.table("../snp.map")
CMplot(dat, plot.type="d")


library(CMplot)
data("pig60K")
CMplot(pig60K, plot.type="d")

```





请问一下大家有什么服务器传文件速度快的软件，大批量数据上传太费时间了

邓飞：fileFTP，速度很快，比xshell子代的xftp还要快。





augustus和genemark预测的文件怎么转换为标准的gff3文件，为什么使用gffread会把gene行过滤掉。

iTools 里面有gff gtf互转的功能









在conda中，`remove`和`uninstall`都用于从环境中卸载包，但它们之间存在一些细微的差异：

- `conda remove`：这个命令会从环境中删除指定的包以及其依赖项，但不会删除任何未使用的依赖项。这意味着如果一个已安装的包是另一个已安装的包的依赖项，那么这个依赖项不会被删除，因为它仍然是另一个包的依赖项。如果您想完全删除一个包及其所有未使用的依赖项，您可以使用`--uninstall`选项。
- `conda uninstall`：这个命令也会从环境中卸载指定的包及其依赖项，但它还会删除未使用的依赖项。这是因为`conda uninstall`命令会检查环境中是否有其他包依赖于要卸载的包，如果没有，它会删除所有未使用的依赖项。因此，`conda uninstall`命令可能会比`conda remove`命令删除更多的包，但也会导致环境中存在较少的未使用的包。

综上所述，如果您想删除一个包及其所有未使用的依赖项，您应该使用`conda uninstall`命令；如果您只想删除指定的包及其直接依赖项，而不删除任何未使用的依赖项，您可以使用`conda remove`命令。





```bash
#TASSEL

#性状文件中缺失数据使用-999表示，可以多个性状在一个文件中


#!/bin/bash
#CSUB -J test_glm
#CSUB -o %J.out 
#CSUB -e %J.error
#CSUB -q cpu
#CSUB -n 12
#CSUB -R "span[hosts=1]"
source /share/apps/anaconda3/bin/activate wgs2
VCF=/share/home/hnu_wangy/project/wgs/297sample/popfilter/bcftools/snps.f4.vcf
run_pipeline.pl -Xmx6g -fork1 -vcf $VCF -fork2 -importGuess S6.tab -combine3 -input1 -input2 -intersect -FixedEffectLMPlugin -endPlugin -export test_glm



ADMITURE的P文件作为群体结构文件，注意去掉最后最一列，保证三列之和小于1
表头加上协变量英文 防止识别错误



-excludeLastTrait This removes last column of Phenotype data. For
example… Can be used to remove last column of
population structure for use with MLM or GLM.

run_pipeline.pl -Xms10g -Xmx100g -fork1 -h snpf4.hmp.txt -fork2 -importGuess S6.tab -fork3 -importGuess pca1.txt -excludeLastTrait -combine4 -input1 -input2 -input3 -intersect -FixedEffectLMPlugin -endPlugin -export test3_glm

<Covariate>
<Trait>
sed -i '1i\<Covariate>\n<Trait>' test



###########################
#测试不同方式加群体结构Q矩阵
./tassel-5-standalone/run_pipeline.pl -fork1 -h test.hmp.txt -fork2 -importGuess test.traits -combine3 -input1 -input2 -intersect -FixedEffectLMPlugin -endPlugin -export te
#直接用ADMIXTURE生成的Q矩阵（无表头信息）进行GLM，TASSEL软件并未读入该矩阵，结果tt和t完全一致
./tassel-5-standalone/run_pipeline.pl -fork1 -h test.hmp.txt -fork2 -importGuess test.traits -fork3 -importGuess LDfiltered.14.Q -combine4 -input1 -input2 -input3 -intersect -FixedEffectLMPlugin -endPlugin -export 2te
#给LDfiltered.14.Q增加表头信息（列名和行名）后生成test.Q矩阵正常读入GLM分析
./tassel-5-standalone/run_pipeline.pl -fork1 -h test.hmp.txt -fork2 -importGuess test.traits -fork3 -importGuess test.Q -combine4 -input1 -input2 -input3 -intersect -FixedEffectLMPlugin -endPlugin -export 3te
#加-excludeLastTrait
./tassel-5-standalone/run_pipeline.pl -fork1 -h test.hmp.txt -fork2 -importGuess test.traits -fork3 -importGuess test.Q -excludeLastTrait -combine4 -input1 -input2 -input3 -intersect -FixedEffectLMPlugin -endPlugin -export 4te


./tassel-5-standalone/run_pipeline.pl -fork1 -h test.hmp.txt -fork2 -importGuess test.traits -fork3 -importGuess test.Q -excludeLastTrait -combine4 -input1 -input2 -input3 -intersect -FixedEffectLMPlugin -endPlugin -export w




#加不加-excludeLastTrait结果一致





./tassel-5-standalone/run_pipeline.pl -Xms10g -Xmx100g -fork1 -h snpf4.hmp.txt -fork2 -importGuess S6.tab -fork3 -importGuess LDfiltered.14.Q -combine4 -input1 -input2 -input3 -intersect -FixedEffectLMPlugin -endPlugin -export test3_glm



#!/bin/bash
#CSUB -J kinship
#CSUB -o %J.out 
#CSUB -e %J.error
#CSUB -q cpu
#CSUB -n 12
#CSUB -R "span[hosts=1]"
source /share/apps/anaconda3/bin/activate wgs2
VCF=/share/home/hnu_wangy/project/wgs/297sample/popfilter/bcftools/snps.f4.vcf
run_pipeline.pl -Xmx100g -vcf $VCF -KinshipPlugin -method Centered_IBS -endPlugin -export kinship.txt -exportType SqrMatrix



#!/bin/bash
#CSUB -J pca
#CSUB -o %J.out 
#CSUB -e %J.error
#CSUB -q fat4
#CSUB -n 6
#CSUB -R "span[hosts=1]"
source /share/apps/anaconda3/bin/activate wgs2
VCF=/share/home/hnu_wangy/project/wgs/297sample/popfilter/bcftools/snps.f4.vcf
run_pipeline.pl -Xms10g -Xmx1000g -vcf $VCF -PrincipalComponentsPlugin -covariance true -endPlugin -export pca



#conda安装不同build的TASSEL运行PCA分析，小数据量测试样本可以正常运行，但大样本不报错自动中止；使用git clone安装TASSEL 5.2.88后运行后有报错信息
#For Tassel 5.0 (Creates installation directory tassel-5-standalone) (Requires Java 8)
#git clone https://bitbucket.org/tasseladmin/tassel-5-standalone.git
#cd tassel-5-standalone/
#git pull

[pool-1-thread-1] ERROR net.maizegenetics.plugindef.AbstractPlugin - Matrix size exceeds the size of an integer
https://groups.google.com/g/tassel/c/GAauDZen94k
umber of Taxa: 297
Number of Sites: 7914042
Sites x Taxa: 2350470474
#自动停止退出原因应该是分析的变异数据过大，可尝试使用LD过滤再TASSEL或使用PLINK

#TASSEL提取部分样本的数据
#提取部分位点的数据

-excludeLastTrait This removes last column of Phenotype data. For
example… Can be used to remove last column of
population structure for use with MLM or GLM



./tassel-5-standalone/run_pipeline.pl -ListPlugins -full TRUE
./tassel-5-standalone/run_pipeline.pl -ListPlugins -usage TRUE
./tassel-5-standalone/run_pipeline.pl -ListPlugins -usage TRUE|less

LINUX中以下任务系统提交命令有什么含义:
#!/bin/bash
#CSUB -J pca
#CSUB -o %J.out 
#CSUB -e %J.error
#CSUB -q fat4
#CSUB -n 4
#CSUB -R "span[hosts=1]"


我应该如何修改将标准输出重定向到作业名称.out文件中。
```





**conda build** 

在使用conda search命令查找包时，会显示多个版本的包及其不同的构建（build）。

Build是指同一版本的软件包的不同构建方式或编译方式。每个build都有一个唯一的标识符，用于区分不同的构建。

不同的build可能包含不同的功能、依赖项或优化选项，这些都取决于软件包的编译方式和使用的编译器。在选择软件包时，可以根据自己的需要选择适合的build版本。

安装特定build版本号的tassel`conda install -c bioconda tassel=5.2.40=1`

当conda search命令显示的软件包版本的build号为1或2时，这通常表示这是该软件包的默认构建版本。

默认的build版本是指在构建软件包时，如果没有指定特定的构建选项，则使用的构建选项。这通常是软件包维护者认为最稳定和最通用的选项。

在大多数情况下，您可以使用默认的构建版本，而不必指定特定的构建号。但是，如果您需要使用特定的构建选项或遇到与默认构建不兼容的问题，则可以指定特定的构建版本。









**你可以把位点做一下ld质控，用质控后的位点作为N，进行P值矫正**

#LD 0.1 369589 ld.prune.in 7544453 ld.prune.out
#LD 0.2 647358 ld.prune.in 7266684 ld.prune.out
#50 10 0.1(admixture mannual example)，但是过滤太多可以用50 10 0.2过滤少些

#1/n=2.71e-6

#1/n=1.54e-6



```bash
awk '$12<1e-5 {print $0}' GROWD5.assoc.result > GROWD5.marked.result

awk '$12<1.54e-6 {print $0}' GROWD5.assoc.result > GROWD5.marked.result


awk '$12<2.71e-6 {print $0}' GROWD5.assoc.result > GROWD5.marked.result
```

















