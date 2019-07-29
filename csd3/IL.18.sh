# 29-7-2019 JHZ

# IL.10R1 data

cut -f2 work/test/IL.18R1-chr2:102810080_A_G.jma.cojo | awk 'NR>1' | sort -k1,1 > l
awk /chr2:/ work/INTERVAL.rsid | sort -k1,1 | join - l > ll
cut -d' ' -f2 ll > lll
export s=/rds/project/jmmh2/rds-jmmh2-projects/olink_proteomics/scallop/INF/INTERVAL/o5000-inf1-outlier_in-r2.sample
qctool  -filetype bgen -g INTERVAL/per_chr/interval.imputed.olink.chr_2.bgen \
        -ofiletype dosage -og IL18.gen.gz -s ${s} \
        -incl-rsids lll
zcat IL18.gen.gz | \
awk '
{
  for (i=1;i<=NF;i++) a[NR,i]=$i
  if(big <= NF) big=NF
}
END {
   for(i=1;i<=big;i++)
   {
     for(j=1;j<=NR;j++) printf OFS a[j,i]; printf "\n"
   }
}' | \
awk 'NR != 1 && NR !=3 && NR != 4 && NR != 5 && NR != 6' > geno
awk 'NR != 2' $s > phencovar

R --no-save -q <<END
  g <- read.table("geno",as.is=TRUE,header=TRUE)
  pc <-read.table("phenocovar",as.is=TRUE,header=TRUE)
  gpc <- merge(g,pc,by.x="SNPID",by.y="ID_1")
  s <- read.table("lll",as.is=TRUE)
  nog <- c("IL.18R1___Q13478~", "age", "sexPulse", "season", "plate", "bleed_to_process_time", paste0("PC",1:20), t(s))
  m1 <- with(gpc, lm(as.formula(paste(nog,g,sep="+"))))
  summary(m1)
END

## SNPTEST v2.5.2

function test()
{
awk '' phenotype
snptest_v2.5.2 -data IL18.gen.gz ${s} \
        -condition_on rs1558649 rs78545931 rs77152652 rs11683213 rs141398063 rs76565432 rs2160203 \
                      rs56151044 rs4851005 rs12987260 rs113030214 rs116635243 rs57942946 rs115725744 rs112893345
        -condtion_on chr2:102810080_A_G chr2:102882352_A_G chr2:102892093_C_T chr2:102904244_A_C chr2:102927130_A_G \
                     chr2:102927649_A_G chr2:102960824_A_G chr2:103011329_A_G chr2:103011552_C_T chr2:103055634_G_T \
                     chr2:103064186_G_T chr2:103074113_A_G chr2:103076057_A_G chr2:103088622_A_G chr2:103236159_C_T
}
