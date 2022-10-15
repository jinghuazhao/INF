#!/usr/bin/bash

# count genomic regions for all sentinels
(
  awk -v OFS='\t' '(NR==1){print $1,$2,$3,$5,$6,$8,$9}' work/INF1.merge
  awk -vd=1e6 -v OFS='\t' '
    (NR>1){
    if($3-$2<=2) {$2=$2-d;$3=$3+d}
    if ($2<0) $2=0
    print $1,$2,$3,$5,$6,$8,$9
  }' work/INF1.merge | \
  sort -k6,6n -k2,2n
) | \
bedtools merge | \
wc -l

awk -vFS="," -vOFS="\t" 'NR>1{print $3,$22,$21}' INF1.jma-rsid.cis.vs.trans | \
sort -k1n -k2n | \
awk '{print "chr" $0}' | \
bedtools merge | \
wc -l

# ST2
R --no-save -q <<END
  cvt<- read.table("work/INF1.merge.cis.vs.trans",as.is=TRUE,header=TRUE)
  m <- read.delim("work/INF1.merge",as.is=TRUE)
  m <- merge(m[c("prot","MarkerName","Start","End")],cvt[c("uniprot","prot","SNP","cis.trans")],by.x=c("prot","MarkerName"),by.y=c("prot","SNP"))
  load("ds/latest/fp/INF1.rda")
  e <- with(tbl,prot=="CCL25" & MarkerName=="chr19:49206145_C_G")
  tbl <- subset(tbl,!e)
  m <- merge(tbl,m,by=c("prot","MarkerName"))
  m <- merge(rsid,m,by="MarkerName")
  s <- setdiff(names(m),c("FreqSE","MinFreq","MaxFreq"))
  write.table(m[s],file="work/INF1.METAL",quote=FALSE,row.names=FALSE,sep="\t")
END

# CNVplot
R --no-save -q <<END
  png("work/regions.png",width=10,height=8,units="in",pointsize=8,res=300)
  xy <- function(x) if (x<23) x else if (x==23) "X" else if (x==24) "Y";
  metal <- read.delim("work/INF1.METAL",as.is=TRUE)[c("Chromosome","Start","End", "prot", "cis.trans")]
  names(metal) <- c("chr","start","end","prot","cis.trans")
  d <- within(metal,{chr<-replace(chr,chr=="X",23); chr<-replace(chr,chr=="Y",24)})
  t <- table(with(metal,chr))
  n <- length(t)
  pos <- vector("numeric")
  for (x in 1:n) pos[x] <- with(subset(d,chr==paste(x)),max(end))
  CM <- cumsum(pos)
  par(xaxt = "n", yaxt = "n")
  xycoords <- xy.coords(c(0,CM), seq(1,5,by=4/n))
  with(xycoords,plot(x, y, type = "n", ann = FALSE, axes = FALSE))
  par(xaxt = "s", yaxt = "s", xpd = TRUE)
  for (x in 1:n) with(subset(d,chr==paste(x)), {
      l <- ifelse(x==1,0,CM[x-1])
      pos[x] <<- ifelse(x == 1, CM[x]/2, (CM[x-1] + CM[x])/2)
      h <- 1+1:5/t[x]
      segments(l+start,h,l+end,h,lwd="3",col=ifelse(cis.trans=="cis","red","blue"))
  })
  axis(1,labels=names(t),at=pos)
  title(main="Flanking regions of sentinels (red=cis, blue=trans)",xlab="Chromosome",ylab="",line=2.5)
  dev.off()
END

Rscript -e '
  cvt <- read.table("work/INF1.merge.out",as.is=TRUE,header=TRUE,nrows=70)
  H <- with(cvt,table(total))
  M <- names(H)
  png(file = "work/signals_by_protein.png",width=7,height=5,units="in",res=300)
  barplot(H,names.arg=M,xlab="No. of pQTL regions",ylab="No. of proteins",
          ylim=c(0,25),col="darkgrey",border="black",cex=0.8,cex.axis=1.5,cex.names=1.5)
  dev.off()
'

# twas_fusion coloc results
function get_fusion()
{
  export tissue=Whole_Blood
  export tag=${tissue}-coloc.dat
  cd ${INF}/work/coloc
  (
    cat *${tag} | head -1 | awk -vOFS="\t" '{print "prot", $0}'
    ls *${tag} | awk -vtag=-${tag} '{sub(tag,"")};1' | \
    parallel --env tag -C' ' 'grep -v PP4 {}-${tag} | awk -vf={} -vOFS="\t" "NR>1&&\$NF>0.5{print f,\$0}"'
  ) > INF1.merge.fusion
}

# fastenloc coloc results
function get_fastenloc()
{
  export tissue=Whole_Blood
  cd ${INF}/work/fastenloc
  (
    cat *${tissue}.sig.out | head -1 | awk -vOFS="\t" '{print "prot", $0}'
    ls *${tissue}.sig.out | awk -vtag=-Whole_Blood.sig.out '{sub(tag,"")};1' | \
    parallel -C' ' 'grep -v RCP {}-${tissue}.sig.out | awk -vf={} -vOFS="\t" "NR>1&&\$NF>0.5{print f,\$0}"'
  ) > INF1.merge.fastenloc-sig
  (
    cat *${tissue}.snp.out | head -1 | awk -vOFS="\t" '{print "prot", $0}'
    ls *${tissue}.snp.out | awk -vtag=-${tissue}.snp.out '{sub(tag,"")};1' | \
    parallel -C' ' 'grep -v SCP {}-${tissue}.snp.out | awk -vf={} -vOFS="\t" "NR>1&&\$NF>0.5{print f,\$0}"'
  ) > INF1.merge.fastenloc-snp
}

# eQTL overlap
awk '!/BDNF/ && NR > 1 {
  if ($3=="\"Q8NF90\"") $7="\"FGF5\""; else if ($3=="\"Q8WWJ7\"") $7="\"CD6\"";print}
' FS='\t' OFS='\t' doc/olink.inf.panel.annot.tsv | \
cut -f2,3,7 | \
cut -f3 | \
sed 's/"//g' > work/inf1.gene

sed '1d' work/INF1.merge.cis.vs.trans | \
cut -d, -f10 | \
uniq > work/INF1.gene

function ms()
{
  gunzip -c ~/rds/results/public/gwas/multiple_sclerosis/discovery_metav3.0.meta.gz | \
  awk -vchr=${chr} -vstart=${start} -vend=${end} 'NR>1 && $1==chr && $2>=start && $2<=end {$3="chr" $1 ":" $2;$6="";print $0,"MS"}' | \
  sort -k3,3 | \
  join -12 -23 ${INF}/work/snp_pos - | \
  awk 'a[$1]++==0' > ${INF}/work/${prot}-QTL-${rsid}.dat
}

function rsid()
{
  sort -k1,1 | \
  join ${INF}/work/INTERVAL.rsid - | \
  awk '{$1="";print}' | \
  awk '{$1=$1};1'

  ln -sf ${INF}/INTERVAL/cardio/INTERVAL.bed ${INF}/MS/INTERVAL.bed
  sort -k2,2 ${INF}/INTERVAL/cardio/INTERVAL.bim | \
  join -12 -21 - ${INF}/work/INTERVAL.rsid | \
  awk '{$1=$2" "$7; $7=""};1' | \
  sort -k1,1n -k4,4n | \
  tr ' ' '\t' > ${INF}/MS/INTERVAL.bim
  ln -sf ${INF}/INTERVAL/cardio/INTERVAL.fam ${INF}/MS/INTERVAL.fam
}

function ieu()
{
  md5sum *.vcf.gz* > vcf.md5sum
  md5sum --check vcf.md5sum
  ls  *.tsv.gz | parallel -C' ' 'echo {}; gunzip -c {} | wc -l' | paste - - | sort -k1,1 > INF1.size
}

# idmap
export gwas_id=${INF}/INTERVAL/per_chr/interval.imputed.olink.chr_9.fam
export inf_id=${INF}/INTERVAL/o5000-inf1-outlier_in-r2.sample
export scallop_inf_id=${HOME}/COVID-19/SCALLOP-Seq

wc -l ${gwas_id}
wc -l <(sed '1,2d' ${inf_id})
join <(cut -f1 ${gwas_id} | sort -k1,1n) <(sed '1,2d' ${inf_id} | cut -d' ' -f1 | sort -k1,1n) | wc -l

module load ceuadmin/stata

stata <<END
local scallop_inf_id : env scallop_inf_id
di "\`scallop_inf_id\'"
insheet using \`scallop_inf_id\'/omicsMap.csv, case clear comma
sort identifier
merge 1:1 identifier using \`scallop_inf_id\'/work/INTERVALdata_28FEB2020
format Olink_inf_QC_24m %15.0g
format Olink_inf_gwasQC_24m %15.0g
format Affymetrix_QC_bl %15.0g
format Affymetrix_gwasQC_bl %15.0g
format Affymetrix_QC_24m %15.0g
format Affymetrix_gwasQC_24m %15.0g

keep identifier Olink_inf_QC_24m Olink_inf_gwasQC_24m Affymetrix_gwasQC_bl ethnicPulse agePulse sexPulse
d,f
tab ethnicPulse
di _N
outsheet if Olink_inf_QC_24m!=. | Olink_inf_gwasQC_24m!=. using idmap.tsv,noquote replace
gen str20 ethnic=ethnicPulse
replace ethnic="EUR" if inlist(ethnicPulse,"Eng/W/Scot/NI/Brit","White Irish")==1
replace ethnic="EAS" if inlist(ethnicPulse,"Asian- Bangladeshi","Asian- Indian","Asian- Pakistani","Chinese")==1
replace ethnic="MID" if ethnicPulse=="Arab"
gen ethnic_NA=ethnic
replace ethnic_NA="NA" if inlist(ethnic,"EUR","EAS","MID")==0
gen FID=0
rename Affymetrix_gwasQC_bl IID
keep if Olink_inf_QC_24m!=. | Olink_inf_gwasQC_24m!=.
outsheet FID IID ethnic using ethnic.txt if IID!=., noquote replace
outsheet FID IID ethnic_NA using ethnic_NA.txt if IID!=., noquote replace
END

cut -f1 ${gwas_id} | grep -f - idmap.tsv | wc -l

# Scaled phenotypes
R --no-save <<END
  suppressMessages(library(dplyr))
  INF <- Sys.getenv("INF")
  ethnic <- read.delim(file.path(INF,"work","ethnic.txt"))
  eur <- filter(ethnic,ethnic %in% c("EUR"))
  f <- file.path(INF,"INTERVAL","o5000-inf1-outlier_in-r2.sample")
  ph <- read.delim(f,sep=" ",nrows=1)
  p <- read.table(f,col.names=names(ph),skip=2)
  eur <- with(p,ID_1 %in% with(eur,IID))
  p[!eur,4:ncol(p)] <- NA
  scatterplot3d::scatterplot3d(p[c("PC1","PC2","PC3")])
  rgl::plot3d(p[c("PC1","PC2","PC3")])
  PCs <- paste0("PC", 1:20)
  covars <- c(c("sexPulse", "age", "season", "plate", "bleed_to_process_time"), PCs)
  proteins <- p[setdiff(names(ph),c("ID_1","ID_2","missing",covars))]
  regress <- function(x)
  {
    fmla <- as.formula(paste(names(proteins)[x]," ~ ", paste(covars, collapse= "+")))
    l <- lm(fmla,data=p)
    r <- proteins[,x]-predict(l,na.action=na.pass)
    scale(r)
  }
  z <- sapply(1:92,regress)
  colnames(z) <- names(proteins)
  rownames(z) <- with(p,ID_1)
  z <- as.data.frame(z)
  write.table(p[eur,c("ID_1")],file.path(INF,"finemapping","eur.id"),row.names=FALSE,col.names=FALSE,quote=FALSE)
  gap::snptest_sample(p,file.path(INF,"finemapping","INTERVAL.sample"),P=names(proteins))
END

function scaled_assoc()
{
  head -1 finemapping/INTERVAL.sample | cut -d' ' -f1-3 --complement | tr ' ' '\n' > ${INF}/finemapping/INTERVAL.prot
  export M=1e6
  awk 'NR>1{print $5,$6,$8,$9}' ${INF}/work/INF1.merge | \
  parallel -C' ' --env M '
  export prot={1}
  export MarkerName={2}
  export chr={3}
  export pos={4}
  export start=$(awk -vpos=${pos} -vflanking=${M} "BEGIN{start=pos-flanking;if(start<0) start=0;print start}")
  export end=$(awk -vpos=${pos} -vflanking=${M} "BEGIN{print pos+flanking}")
  function bgen()
  {
    qctool -g ${INF}/INTERVAL/per_chr/interval.imputed.olink.chr_${chr}.bgen \
           -s ${INF}/finemapping/INTERVAL.sample \
           -incl-range ${start}-${end} \
           -og ${INF}/finemapping/${prot}-${MarkerName}.bgen -os ${INF}/finemapping/${prot}-${MarkerName}.sample
    bgenix -g ${INF}/finemapping/${prot}-${MarkerName}.bgen -index -clobber
  }
  bgen
  export phenocol=$(grep ${prot}_ ${INF}/finemapping/INTERVAL.prot)
  function assoc_test()
  {
    snptest \
            -data ${INF}/finemapping/${prot}-${MarkerName}.bgen ${INF}/finemapping/INTERVAL.sample \
            -log ${INF}/finemapping/${prot}-${MarkerName}-snptest.log \
            -filetype bgen \
            -frequentist 1 -hwe -missing_code NA,-999 -use_raw_phenotypes \
            -method score \
            -pheno ${phenocol} -printids \
            -o ${INF}/finemapping/${prot}-${MarkerName}.out
  }
  assoc_test
 '
}

function tbi_INF()
# to generate tbl.gz/tbl.gz.tbi from ${INF}/METAL
{
  export rt=${INF}/METAL
  ls ${rt}/*tbl.gz | \
  xargs -I{} basename {} .tbl.gz | \
  parallel -j3 --env INF --env rt -C' ' '
    (
      gunzip -c ${rt}/{}.tbl.gz | head -1
      gunzip -c ${rt}/{}.tbl.gz | \
      sed "1d" | \
      awk -f ${INF}/tryggve/metal.awk | \
      sort -k1,1n -k2,2n
    ) | \
    bgzip -f > {}.tbl.gz; \
    tabix -f -S1 -s1 -b2 -e2 {}.tbl.gz
  '
}

function tbi_INTERVAL()
{
  cd ${INF}/sumstats
  ls INTERVAL/ | sed 's/INTERVAL.//;s/.gz//' | parallel -j15 -C' ' 'gunzip -c INTERVAL/INTERVAL.{}.gz | bgzip -f > INTERVAL.{}.gz'
  mv INTERVAL.* INTERVAL
  ls INTERVAL/ | sed 's/INTERVAL.//;s/.gz//' | parallel -j15 -C' ' 'tabix -f -S1 -s2 -b3 -e3 INTERVAL/INTERVAL.{}.gz'
}

function fp()
{
R --no-save -q <<END
  METAL_forestplot <- function(tbl,all,rsid,package="meta",split=FALSE,...)
  {
  prot <- MarkerName <- NA
  requireNamespace("dplyr")
  dplyr_rsid <- function(df,rsid)
  {
  # d <- dplyr::nest_join(df,rsid)
    d <- dplyr::left_join(df,rsid)
  # dy <- d["y"]
    m <- within(d, {
  #   rsid <- ifelse(length(lapply(dy,"[[",1)) == 1, unlist(d[["y"]]), unlist(lapply(lapply(dy,"[[",1),"[",1)))
  #   rsid <- unlist(d[["y"]])
      isna <- is.na(rsid)
      rsid[isna] <- MarkerName[isna]
    })
  }
  t <- dplyr_rsid(tbl,rsid)
  a <- dplyr_rsid(all,rsid)
  for(i in 1:nrow(tbl))
  {
     p <- tbl[i,"prot"]
     target <- dplyr::filter(gap.datasets::inf1,prot==p) %>% select(target.short)
     m <- tbl[i,"MarkerName"]
     A1 <- toupper(tbl[i,"Allele1"])
     A2 <- toupper(tbl[i,"Allele2"])
     print(paste0(i,"-",p,":",m))
     with(subset(all,prot==p & MarkerName==m), {
       print(subset(all,prot==p & MarkerName==m))
       e <- toupper(EFFECT_ALLELE)
       r <- toupper(REFERENCE_ALLELE)
       a1 <- e
       a2 <- r
       c <- rep(1,length(e))
       j <- sapply(a1,'!=',A1)
       a1[j] <- r[j]
       a2[j] <- e[j]
       c[j] <- -1
       print(cbind(A1,A2,EFFECT_ALLELE,REFERENCE_ALLELE,a1,a2,format(BETA,digits=3),format(BETA*c,digits=3)))
       BETA <- BETA * c
       title <- sprintf("%s [%s (%s) (%s/%s) N=%.0f]",target,m,t[i,"rsid"],A1,A2,tbl[i,"N"])
       if (split) pdf(paste0(p,"-",m,".pdf"),...)
       requireNamespace("meta")
       mg <- meta::metagen(BETA,SE,sprintf("%s (%.0f)",study,N),title=title)
       meta::forest(mg,colgap.forest.left = "1cm")
       requireNamespace("grid")
       grid::grid.text(title,0.5,0.9)
       with(mg,cat("prot =", p, "MarkerName =", m, "Q =", Q, "df =", df.Q, "p =", pval.Q, "I2 =", I2, "lower.I2 =", lower.I2, "upper.I2 =", upper.I2, "\n"))
       if (split) dev.off()
     })
  }
  }
  INF <- Sys.getenv("INF")
  load(file.path(INF,"ds","latest","fp","INF1.rda"))
  drop <- names(tbl)[c(6:9,12:17)]
  library(dplyr)
  tbl_ord <- tbl %>%
             select(-any_of(drop)) %>%
             filter(!(prot=="CCL25" & MarkerName=="chr19:49206145_C_G")) %>%
             arrange(prot,Chromosome,Position)
# 180 forest plots
  METAL_forestplot(tbl_ord,all,rsid,split=TRUE,width=8.75,height=5)
END
}

function pdf()
{
  cd ${INF}/ds/latest
  if [ ! -d work ]; then mkdir work; fi
# forest/locuszoom top-down format
  qpdf --empty --pages $(ls lz/*.pdf | grep -v CCL25-chr19:49206145_C_G.lz.pdf) -- lz2.pdf
  qpdf -show-npages lz2.pdf
  qpdf --pages . 1-360:odd -- lz2.pdf lz.pdf
  rm lz2.pdf
  rm -f work/*pdf
  cd work; fp; cd -
  qpdf --empty --pages $(ls work/*pdf) -- fp.pdf
# left-right with very small file size
# Split files, note the naming scheme
  pdfseparate fp.pdf temp-%04d-fp.pdf
  pdfseparate lz.pdf temp-%04d-lz.pdf
# Combine the final pdf
  pdfjam temp-*-*.pdf --nup 2x1 --landscape --papersize '{5in,16in}' --outfile fp+lz.pdf
# Clean up
# Images for the GitHub page using output from qqman.sb (no border for Q-Q plot)
  convert -density 300 -resize 110% work/fp-lz-OPG-chr17:26694861_A_G.png OPG.png
  convert ${INF}/plots/work/OPG-qqman.png OPG.png -append -density 300 ~/INF/doc/OPG.png
  rm OPG.png
  cd ~/EWAS-fusion/IL.12B.tmp
  pdftopng -r 300 ewas-plot.pdf ewas-plot
  export rt=ewas-plot-00000
  convert \( ${rt}1.png ${rt}2.png +append \) \( ${rt}3.png ${rt}4.png +append \) -append ewas-plot.png
}

function pdf_test()
{
# OCR/resolution is poor (though layout is nice) with the following:
  qpdf fp-lz.pdf --pages . 1 -- 1.pdf
  qpdf fp-lz.pdf --pages . 2 -- 2.pdf
  convert 1.pdf 2.pdf -density 300 +append 12.pdf
  rm 1.pdf 2.pdf
# The .tiff format is possible with the tiff64 tag but too large
  qpdf --empty -collate --pages fp.pdf lz.pdf -- fp+lz.pdf
# forest/locuszoom side-by-side format, OCR via PDF-viewer and compressed by Adobe
  rm -f temp-*-*.pdf
  source ~/COVID-19/py37/bin/activate
# pip install img2pdf
  awk -vFS="\t" 'NR>1 {print $5,$6}' ${INF}/work/INF1.merge | \
  parallel -C' ' '
    export rt={1}-{2}
    export suffix=$(printf "%06d\n" 1)
    if [ ! -f work/fp-${rt}.png ]; then pdftopng -r 300 -f 1 -l 1 work/${rt}.pdf work/fp-${rt}; fi
    if [ ! -f work/lz-${rt}.png ]; then pdftopng -r 300 -f 1 -l 1 lz/${rt}.lz.pdf work/lz-${rt}; fi
    convert work/fp-${rt}-${suffix}.png work/lz-${rt}-${suffix}.png -density 300 +append work/fp-lz-${rt}.png
    convert work/fp-lz-${rt}.png -quality 0 work/fp-lz-${rt}.jp2
    img2pdf -o work/fp-lz-${rt}.pdf work/fp-lz-${rt}.jp2
    rm work/fp-${rt}-${suffix}.png work/lz-${rt}-${suffix}.png work/fp-lz-${rt}.jp2
  '
  qpdf --empty --pages $(ls work/fp-lz-*.pdf) -- fp+lz.pdf
# convert fp+lz.pdf -density 300 tiff64:fp+lz.tiff
# qml/
# 91 Q-Q/Manhattan (left+right collation dropping cis-locuszoom) and tif via PDF-viewer
  rm -f work/*pdf
  ls qml/*qq*png | xargs -l basename -s .qq.png | grep -v BDNF | \
  parallel -C' ' '
    convert qml/{}.manhattan.png qml/{}.qq.png -density 300 +append work/{}.png
    convert work/{}.png -quality 0 work/{}.jp2
    img2pdf -o work/{}.pdf work/{}.jp2
    rm work/{}.jp2
  '
  qpdf --empty --pages $(ls work/*.pdf) -- qq+manhattan.pdf
# Not working very well
# see https://legacy.imagemagick.org/Usage/layers/
  convert $(paste -d ' ' <(ls *qq*) <(ls *manhattan*) | xargs -l -I {} echo '\(' {} +append '\)' '\')
          -append qq_manhattan.pdf
# convert qq_manhattan.pdf -density 300 tiff64:qq+manhattan.tiff
# locuszoom plots for 91 cis-regions are possible with pdfunite but got complaints from qpdf
# pdfunite *.pdf ~/lz.pdf
}

function ppi()
{
Rscript -e '
  options(width=200)
  INF <- Sys.getenv("INF")
  library(dplyr)
  inf1 <- pQTLtools::inf1
  head(inf1)
  write.table(inf1,file=file.path(INF,"ppi","inf1.txt"),quote=FALSE,row.names=FALSE,sep="\t")
  library(openxlsx)
  sheet1 <- read.xlsx(file.path(INF,"ppi","humphreys21.xlsx")) %>%
            select(-c("ourID","pdbmono1","pdbmono2","KEGGpath1","KEGGpath2","Phob1","Phob2","UniFun1","UniFun2","Uniloc1","Uniloc2","UniPubsInt"))
  write.table(sheet1,file=file.path(INF,"ppi","ppi.txt"),quote=FALSE,row.names=FALSE,sep="\t")
  head(sheet1)
  filter(sheet1,Uniprot1 %in% inf1$uniprot | Uniprot2 %in% inf1$uniprot)
  '
  sed '1d' ${INF}/ppi/inf1.txt | cut -f5 | grep -f - -w ${INF}/ppi/ppi.txt
}

function googlesheet()
# https://www.hdfstutorial.com/blog/connect-r-google-sheet-use/
{
  Rscript -e '
    url <- "https://docs.google.com/spreadsheets/d/1nLqPlG5RwGcxAzolYI5bsctIVGl25f_W/edit?usp=sharing&ouid=102539576739161755759&rtpof=true&sd=true"
    library(googlesheets)
    gs_auth(new_user = TRUE)
    gs_ls()
    for_gs <- gs_title("mtcars")
    for_gs_sheet <- gs_read(for_gs)
    gs_new(title = "IRIS Dataset", ws_title = "first_sheet", input = iris)
  '
}

function f2()
# Figure 2 for the Google document
{
  convert signals_by_protein.png -resize 110% 2a.png
  convert hotspot-rs12075.png -resize 60% 2b.png
  convert IL.12B-mhtplot.trunc.png -resize 80% 2c.png
  convert TRAIL-mhtplot.trunc.png -resize 80% 2d.png
  convert +append 2a.png 2b.png f2-1.png
  convert +append 2c.png 2d.png f2-2.png
  convert -append f2-1.png f2-2.png f2.png
  rm 2a.png 2b.png 2c.png 2d.png f2-1.png f2-2.png
}

function pdf_final()
# Rework
{
  cd ${INF}/METAL/qqmanhattanlz
# qpdf --empty --pages $(ls *_rs*.pdf) -- lz2.pdf
  qpdf --empty --pages $(ls *_rs*.pdf | \
                         xargs -l basename -s .pdf | \
                         join - <(awk 'NR>1{print $3"_"$2,$1}' ${INF}/work/INF1.METAL | sort -k1,1) | \
                         sed 's/_/ /' | \
                         sort -k1,1 -k3,3 | \
                         awk '{print $1"_"$2".pdf"}') -- lz2.pdf
  qpdf -show-npages lz2.pdf
  qpdf --pages . 1-360:odd -- lz2.pdf lz.pdf
  rm lz2.pdf
  pdfseparate ${INF}/ds/latest/fp.pdf temp-%04d-fp.pdf
  pdfseparate lz.pdf temp-%04d-lz.pdf
# Combine the final pdf
  pdfjam temp-*-*.pdf --nup 2x1 --landscape --papersize '{5in,16in}' --outfile fp+lz.pdf
  rm temp*pdf
  ls *_qq.png | xargs -l basename -s _qq.png | \
  parallel -C' ' 'convert -resize 150% {}_qq.png {}_qq.pdf;convert {}_manhattan.png {}_manhattan.pdf'
  qpdf --empty --pages $(ls *_qq.pdf) -- qq.pdf
  qpdf --empty --pages $(ls *_manhattan.pdf) -- manhattan.pdf
  pdfseparate qq.pdf temp-%04d-qq.pdf
  pdfseparate manhattan.pdf temp-%04d-manhattan.pdf
  pdfjam temp-*-*.pdf --nup 2x1 --landscape --papersize '{5in,16in}' --outfile qq-manhattan.pdf
}

function cis_info()
{
  awk '$21=="cis" {print $3}' ${INF}/work/INF1.METAL | sort | uniq | grep -w -f - ${INF}/work/INF1.merge.genes | \
  awk -vM=1e6 '{$4=$4-M;$5=$5+M};1' | \
  sort -k1,1 | join <(cut -f3,7 doc/olink.inf.panel.annot.tsv | sed 's/"//g' | sort -k1,1) - > TNFB/cis.dat
}

function gsmr_png()
{
# https://bioinformatics.psb.ugent.be/webtools/Venn/
# https://shiny.cnsgenomics.com/LMOR/
  convert gsmr-5e-8-3/gsmr-efo.png gsmr-5e-8-10/gsmr-efo.png -append -density 300 -resize 20% gsmr.png
  convert +append gsmr-5e-8-3/gsmr-efo.png gsmr-5e-8-10/gsmr-efo.png -resize 10% gsmr-lr.png
  convert -resize 25% protein_venn_diagram.png latest/protein_venn_diagram.png
}

function chrpos_rsid()
{
  Rscript -e '
    suppressMessages(library(dplyr))
    suppressMessages(library(ieugwasr))
    INF <- Sys.getenv("INF")
    cis <- read.table(file.path(INF,"TNFB","cis.dat"),col.names=c("uniprot","gene","prot","chr","start","end")) %>%
           mutate(region=paste0(chr,":",start,"-",end))
    info <- sapply(cis$region, function(x) variants_chrpos(x) %>% select(name,chr,pos))
    sapply(1:59,function(x)
           write.table(info[,x],file=file.path(INF,"mr","gsmr","region",paste0(cis[x,"prot"],".rsid.txt")),quote=FALSE,row.names=FALSE,sep="\t"))
  '
}
