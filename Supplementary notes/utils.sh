#!/usr/bin/bash

# count genomic regions for all sentinels
cut -d, -f3,4 ${INF}/work/INF1.jma-rsid.cis.vs.trans | \
tr ',' '\t' | \
awk -vflanking=1e6 -vOFS='\t' 'NR>1{print $1,$2-flanking,$2+flanking}' | \
awk -vOFS='\t' '{if(NR>1&&$2<0) $2=0};1' | \
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
     gene <- dplyr::filter(gap.datasets::inf1,prot==p) %>% select(gene)
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
       title <- sprintf("%s (%s) [%s (%s) (%s/%s) N=%.0f]",target,gene,m,t[i,"rsid"],A1,A2,tbl[i,"N"])
       if (split) pdf(paste0(p,"-",m,".pdf"),...)
       requireNamespace("meta")
       mg <- meta::metagen(BETA,SE,sprintf("%s (%.0f)",study,N),title=title)
       meta::forest(mg,colgap.forest.left = "0.5cm",digits.TE=3,digits.se=2)
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
  METAL_forestplot(tbl_ord,all,rsid,split=TRUE,width=9,height=5)
END
}

function pdf()
{
  cd ${INF}/ds/latest
  if [ ! -d work ]; then mkdir work; fi
# forest/locuszoom top-down format
  qpdf --empty --pages $(ls *_*.pdf | grep -v CCL25-chr19:49206145_C_G.lz.pdf) -- lz2.pdf
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
  awk -vFS="\t" 'NR>1 {print $3,$1,$2}' ${INF}/work/INF1.METAL | sort -k1,1 -k2,2 | \
  parallel -C' ' '
    export rt={1}-{2}
    export r_t={1}_{3}
    export suffix=$(printf "%06d\n" 1)
    if [ ! -f fp-${rt}.png ]; then pdftopng -r 300 -f 1 -l 1 fp/${rt}.pdf fp-${rt}; fi
    if [ ! -f lz-${rt}.png ]; then
       pdftopng -r 300 -f 1 -l 1 qqmanhattanlz/${r_t}.pdf lz-${r_t}; convert -density 300 -size 5x9 lz-${r_t}-${suffix}.png lz-${rt}.png;
    fi
    convert +append fp-${rt}-${suffix}.png lz-${rt}.png -resize x500 -density 300 fp-lz-${rt}.png
    convert fp-lz-${rt}.png -quality 0 fp-lz-${rt}.jp2
    img2pdf -o fp-lz-${rt}.pdf fp-lz-${rt}.jp2
    rm fp-${rt}-${suffix}.png lz-${r_t}-${suffix}.png lz-${rt}.png fp-lz-${rt}.jp2
  '
  qpdf --empty --pages $(ls fp-lz-*.pdf) -- SF-fp-lz.pdf
  rm fp-lz-*.*
# convert fp+lz.pdf -density 300 tiff64:fp+lz.tiff
# qml/
# 91 Q-Q/Manhattan (left+right collation dropping cis-locuszoom) and tif via PDF-viewer
  ls qqmanhattanlz/*qq*png | xargs -l basename -s _qq.png | grep -v BDNF | \
  parallel -C' ' '
    convert +append qqmanhattanlz/{}_manhattan.png qqmanhattanlz/{}_qq.png -resize x500 -density 300 {}.png
    convert {}.png -quality 0 {}.jp2
    img2pdf -o {}.pdf {}.jp2
    rm {}.jp2
  '
  qpdf --empty --pages $(ls *.pdf) -- SF-manhattan-qq.pdf
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
# f2 hotspot-rs12075.png
{
  cd ${INF}/work
  export figure2a=${INF}/hotspots/$1
  convert ${figure2a} -resize 80% 2a.png
  convert signals_by_protein.png -resize 110% 2b.png
  convert IL.12B-mhtplot.trunc.png -resize 80% 2c.png
  convert TRAIL-mhtplot.trunc.png -resize 80% 2d.png
  convert +append 2a.png 2b.png f2-1.png
  ln -sf 2a.png f2-1.png
  convert +append 2c.png 2d.png f2-2.png
  convert -append f2-1.png f2-2.png f2.png
  rm 2a.png 2b.png 2c.png 2d.png f2-1.png f2-2.png
  cd -
}

f2 hotspot-rs3184504.png

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
  convert -resize 25% protein_MR_GSMR.png latest/protein_MR_GSMR.png
  convert -resize 20% gsmr-5e-8-3/gsmr-efo-reduce.png gsmr-efo.png
  convert -resize 20% mr-efo.png latest/mr-efo.png
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

function llod()
{
  Rscript -e '
    library(dplyr)
    annot <- readRDS(file.path("~","pQTLtools","tests","annot.RDS")) %>%
             left_join(pQTLdata::inf1[c("prot","gene")],by=c("ID"="prot")) %>%
             mutate(order=1:n()) %>%
             arrange(desc(order))
    xtick <- seq(1, nrow(annot))
    INF <- Sys.getenv("INF")
    attach(annot)
    png(file.path(INF,"work","llod.png"),width=10,height=8,units="in",pointsize=8,res=300)
    par(mar=c(5,4,1,1))
    plot(100-annot$pc.belowLOD.new, col=annot$pQTL, las=1,
         ylab = "100 - % samples with very low abundance per protein", xaxt= "n",
         xlab = "Ordered proteins", cex=0.8, pch=19)
    axis(side=1, at = seq(from=0, to=nrow(annot), by=20))
    legend("topright", legend= c("no pQTL", "pQTL"), col= c("red","blue"), pch = 19)
    dev.off()

    png(file.path(INF,"work","llod2.png"),width=10,height=8,units="in",pointsize=8,res=300)
    par(mar=c(5,5,1,1))
    plot(100-pc.belowLOD.new,pch=19,cex=2,col=pQTL,axes=FALSE,ann=FALSE)
    axis(1,cex.axis=1.5,lwd.tick=0.5,line=0)
    axis(2,cex.axis=1.5,lwd.tick=0.5,line=0)
    legend(x=1.5,y=-0.5,c("cis","trans"),box.lwd=0,cex=2,col=c("red","blue"),pch=19)
    mtext("% samples above LLOD",side=2,line=2.5,cex=1.5)
    mtext("Ordered proteins",side=1,line=2.5,cex=1.5)
    dev.off()

    png(file.path(INF,"work","SF-LLOD.png"),width=15,height=8,units="in",pointsize=8,res=300)
    par(mar=c(10,5,1,1))
    plot(100-pc.belowLOD.new,cex=2,pch=19,col=pQTL,xaxt="n",xlab="",ylab="",cex.axis=1.2)
    text(66,16,"IL-17C",offset=0,pos=2,cex=1.5,font=2,srt=0)
    arrows(67,16,71,16,lwd=2)
    axis(1, at=xtick, labels=gene, lwd.tick=0.5, lwd=0, las=2, hadj=1, cex.axis=1.2)
    mtext("% samples above LLOD",side=2,line=2.5,cex=1.5)
    mtext("Ordered protein",side=1,line=8.5,cex=1.5,font=1)
    legend(x=1,y=25,c("without pQTL","with pQTL"),box.lwd=0,cex=2,col=c("red","blue"),pch=19)
    dev.off()

    detach(annot)
  '
}

function pve()
{
R --no-save -q <<END
  INF <- Sys.getenv("INF")
  png(file.path(INF,"h2","SF-Effect-size-PVE.png"), res=300, units="cm", width=30, height=30)
  jma <- read.delim(file.path(INF,"work","INF1.jma-rsid"))
  jma.cistrans <- read.csv(file.path(INF,"work","INF1.jma-rsid.cis.vs.trans"))[c("prot","cis.trans")]
  INF1_jma <- merge(jma,jma.cistrans,by="prot")
  with(INF1_jma,
  {
    MAF <- freq
    repl <- MAF > 1-MAF
    MAF[repl] <- 1-MAF[repl]
    Effect <- bJ
    v <- 2*MAF*(1-MAF)*Effect^2
    col <- c("blue","red")[1+(cis.trans=="trans")]
#   plot(MAF,abs(Effect),cex.axis=1.3,cex.lab=1.3,pch=19,main="a",xlab="MAF", ylab="Effect size",col=col)
    library(scatterplot3d)
    scatterplot3d(MAF,abs(Effect),v,color=col, pch=16, type="h",
                  xlab="MAF", ylab="Effect", zlab=expression(italic(2*MAF(1-MAF)*Effect^2)), cex.axis=1.3, cex.lab=1.3)
    legend("right", legend=levels(as.factor(cis.trans)), box.lwd=0, col=c("red", "blue"), pch=16)
  })
  dev.off()
END
}

function IL17C()
# gunzip -c ${INF}/METAL/IL.17C-1.tbl.gz | cut -f1-3,10-11 | gzip -f > ${INF}/work/INF-IL.17C.gz
{
Rscript -e '
  suppressMessages(library(dplyr))
  suppressMessages(library(gap))
  INF <- Sys.getenv("INF")
  genes <- data.frame(MarkerName=c("chr16:88684495_G_T"),gene=c("IL17C"),color=c("red"))
  IL.17C <- read.delim(file.path(INF,"work","INF-IL.17C.gz"),as.is=TRUE) %>%
            mutate(Z=Effect/StdErr,P=as.numeric(pvalue(Z)),gene=NA,color=NA) %>%
            filter(!is.na(Z)) %>%
            select(Chromosome,Position,MarkerName,Z,P,color) %>%
            arrange(Chromosome,Position)
  loc <- with(IL.17C,MarkerName=="chr16:88684495_G_T")
  IL.17C[loc,"gene"] <- "IL17C"
  mhtdata <- select(IL.17C,Chromosome,Position,P,color,gene)
  cis <- with(mhtdata,Chromosome==16 & Position>=88704999-1e6 & Position<=88706881+1e6)
  mhtdata[cis,"color"] <- "IL17C"
  subset(mhtdata,!is.na(gene))
  png(file.path(INF,"work","INF-IL.17C-mhtplot.png"), res=300, units="in", width=12, height=8)
    par(cex=0.8, mar=c(6,6,3,1))
    ops <- mht.control(colors=rep(c("blue4","skyblue"),11),srt=0,yline=3)
    hops <- hmht.control(data=filter(mhtdata,!is.na(gene)),colors="red")
    mhtplot2(mhtdata,ops,hops,xlab="Chromosome",ylab="-log10(P)",srt=0, cex.axis=2)
    abline(h=-log10(5e-10),col="red")
    axis(2,at=0:15)
  dev.off()
'
}

function gene_annotation()
{
Rscrtip -e '
  library(dplyr)
  genes_trans <- read.delim("~/INF/work/genes_trans.tsv") %>%
                 setNames(c("rsid","causal_gene"))
  inf1_metal <- read.delim("~/INF/work/INF1.METAL",as.is=TRUE) %>%
                left_join(pQTLdata::inf1[c("prot","gene")]) %>%
                left_join(genes_trans) %>%
                mutate(causal_gene=if_else(cis.trans=="cis",gene,causal_gene))
  write.table(inf1_metal,file="inf1_metal.tsv",row.names=FALSE,quote=FALSE,sep="\t")
'
}

function gsmr3()
{
  export efo=ebi-a-GCST004133
  export prot=CXCL5
  gunzip -c ${INF}/mr/gsmr/out/${prot}-${efo}.eff_plot.gz | awk '/effect_begin/,/effect_end/' | grep -v effect | cut -d' ' -f1 | \
  grep -f - ${INF}/work/INTERVAL.rsid > ${INF}/${prot}-${efo}.rsid
  Rscript -e '
    options(width=200)
    INF <- Sys.getenv("INF")
    efo <- Sys.getenv("efo")
    p <- Sys.getenv("prot")
    trait <- "Ulcerative colitis"
    r <- "rs450373"
    gene <- "CXCL5"
    rt <- paste0(INF,"/",p,"-",efo)
    snpid_rsid <- read.table(paste0(rt,".rsid"),col.names=c("snpid","SNP"))
    suppressMessages(require(dplyr))
    suppressMessages(require(MendelianRandomization))
    source(file.path(INF,"rsid","gsmr_plot.r"))
    eff_dat <- paste0(INF,"/mr/gsmr/out/",p,"-",efo,".eff_plot.gz")
    gsmr_data <- read_gsmr_data(eff_dat)
    gsmr_summary(gsmr_data)
    attach(gsmr_data)
    gsmr_snp_effect <- as.data.frame(snp_effect[,-(2:3)]) %>%
                       setNames(c("snpid","effect_allele","other_allele","eaf","bx","bxse","by","byse")) %>%
                       mutate(bx=as.numeric(bx),bxse=as.numeric(bxse),by=as.numeric(by),byse=as.numeric(byse))
    detach(gsmr_data)
    gsmr_snp_effect <- snpid_rsid %>%
                       left_join(gsmr_snp_effect)
    attach(gsmr_snp_effect)
    png(paste0(INF,"/SF-",p,"-",trait,".png"),res=300,height=5,width=4,units="in")
    mr_dat <- mr_input(bx = bx, bxse = bxse, by = by, byse = byse, exposure = gene, outcome = trait, snps = SNP,
              effect_allele = effect_allele, other_allele = other_allele, eaf = eaf)
    mr_plot(mr_dat, interactive = FALSE, labels = FALSE)
    dev.off()
    detach(gsmr_snp_effect)
  # j <- left_join(gsmr_snp_effect,snpid_rsid) %>% pull(SNP)
  # rsids <- c(r,j)
  # ldr <- ieugwasr::ld_matrix(rsids,pop="EUR",with_alleles=FALSE)
  # write.table(ldr,file=paste0(p,"-",r,".r"),quote=FALSE,sep="\t")
  # pdf(paste0(INF,"/",p,"-",trait,".eff_plot.pdf"))
  # par(mar=c(6,6,5,1),mgp=c(4,1,0),xpd=TRUE)
  # plot_gsmr_effect(gsmr_data, p, efo, colors()[75])
  # dev.off()
  '
}

function GCST()
{
#!/usr/bin/bash

#SBATCH --job-name=gsmr
#SBATCH --account CARDIO-SL0-CPU
#SBATCH --partition cardio
#SBATCH --qos=cardio
#SBATCH --mem=28800
#SBATCH --time=5-00:00:00
#SBATCH --output=/rds/user/jhz22/hpc-work/work/GCST_%A_%a.out
#SBATCH --error=/rds/user/jhz22/hpc-work/work/GCST_%A_%a.err
#SBATCH --export ALL

export TMPDIR=/rds/user/jhz22/hpc-work/work

Rscript -e '
  library(dplyr)
  library(gwasrapidd)
  INF <- Sys.getenv("INF")
  GCST <- gwasrapidd::get_associations(variant_id = c(unique(ps_na_disease)$rsid))
  GCST_id <- pull(GCST@associations,association_id)
  assoc_study <- association_to_study(GCST_id)
  study_id <- pull(assoc_study, study_id)
  GCST_traits <- get_traits(study_id)
  GCST_studies <- get_studies(study_id)
  save(GCST,assoc_study,GCST_studies,GCST_traits,file=file.path(INF,"work","GCST.rda"))
'
}

function CXCL5_prep()
{
  export eQTLGen=~/rds/public_databases/eQTLGen
  export cis_full=${eQTLGen}/cis-eQTLs_full_20180905.txt.gz
  export eQTLGen_tabix=${eQTLGen}/tabix
  export GTEx=~/rds/public_databases/GTEx/csv
  export chr=4
  export start=74861359
  export end=74864496
  export M=1e6
# SCALLOP/INF
  for prot in CXCL5
  do
    gunzip -c ${INF}/METAL/${prot}-1.tbl.gz | \
    awk -vOFS="\t" -vchr=${chr} -vstart=${start} -vend=${end} -vM=${M} '
    {
      if ($1 == chr && $2 >= start - M && $2 <= end + M)
      {
        split($3,a,"_")
        print a[1],$1,$2,$10/$11,$3,toupper($4),toupper($5)
      }
    }' | \
    sort -k1,1 | \
    join -12 -21 <(grep chr${chr} ${INF}/work/snp_pos) - | \
    awk -vOFS="\t" '{print $6, $2, $3, $4, $7, $8, $5}' | \
    gzip -f > ${INF}/work/${prot}.tsv.gz
  done
  export region=4:73861359-75864496
  export CD=${INF}/OpenGWAS/ebi-a-GCST004132.vcf.gz
  cat <(tabix -H ${CD} ${region} | tail -n -1) \
      <(tabix ${CD} ${region}) | \
  cut -f1-5,10 | \
  awk -vOFS="\t" '{
         if (NR==1) print "snpid","rsid","chr","pos","a1","a2","z";
         else {
           if ($4<$5) snpid="chr"$1":"$2"_"$4"_"$5; else snpid="chr"$1":"$2"_"$5"_"$4;
           split($6,a,":"); z=a[1]/a[2];
           print snpid, $3, $1, $2, $5, $4, z
         }
       }' | \
  gzip -f > ${INF}/work/CD.tsv.gz
  export UC=${INF}/OpenGWAS/ebi-a-GCST004133.vcf.gz
  cat <(tabix -H ${UC} ${region} | tail -n -1) \
      <(tabix ${UC} ${region}) | \
  cut -f1-5,10 | \
  awk  -vOFS="\t" '{
         if (NR==1) print "snpid","rsid","chr","pos","a1","a2","z";
         else {
           if ($4<$5) snpid="chr"$1":"$2"_"$4"_"$5; else snpid="chr"$1":"$2"_"$5"_"$4;
           split($6,a,":"); z=a[1]/a[2];
           print snpid, $3, $1, $2, $5, $4, z
         }
       }' | \
  gzip -f > ${INF}/work/UC.tsv.gz
# done via eQTLGen.sh
# cat <(gunzip -c ${cis_full} | head -1) <(gunzip -c ${cis_full} | sed '1d' | sort -k3,3n -k4,4n) | \
# bgzip -f > ${eQTLGen_tabix}/cis_full.txt.gz
# tabix -S1 -s3 -b4 -e4 -f ${eQTLGen_tabix}/cis_full.txt.gz
  cat <(gunzip -c ${eQTLGen_tabix}/cis_full.txt.gz | head -1) \
      <(tabix ${eQTLGen_tabix}/cis_full.txt.gz  ${region}) |
  awk 'NR==1||/CXCL5/' | \
  cut -f2-7 | \
  awk -vOFS="\t" '{
         if (NR==1) print "snpid","rsid","chr","pos","a1","a2","z";
         else {
           if ($5<$6) snpid="chr"$2":"$3"_"$5"_"$6; else snpid="chr"$2":"$3"_"$6"_"$5;
           print snpid, $1, $2, $3, $5, $6, $4
         }
       }' | \
  gzip -f > ${INF}/work/eQTLGen.tsv.gz
# GTEx v8 data which requires liftOver
# 1 variant
# 2 r2
# 3 pvalue
# 4 molecular_trait_object_id
# 5 molecular_trait_id
# 6 maf
# 7 gene_id
# 8 median_tpm
# 9 beta
# 10 se
# 11 an
# 12 ac
# 13 chromosome
# 14 position
# 15 ref
# 16 alt
# 17 type
# 18 rsid
  Rscript -e '
    library(dplyr)
    chr <- Sys.getenv("chr") %>% as.numeric()
    start <- Sys.getenv("start") %>% as.numeric()
    end <- Sys.getenv("end") %>% as.numeric()
    flanking <- Sys.getenv("M") %>% as.numeric()
    region <- Sys.getenv("region")
    HPC_WORK <- Sys.getenv("HPC_WORK")
    f <- file.path(HPC_WORK,"bin","hg19ToHg38.over.chain")
    chain <- rtracklayer::import.chain(f)
    require(GenomicRanges)
    gr <- GenomicRanges::GRanges(seqnames=chr,IRanges::IRanges(start,end)+flanking)
    seqlevelsStyle(gr) <- "UCSC"
    gr38 <- rtracklayer::liftOver(gr, chain)
    chr <- gsub("chr","",colnames(table(seqnames(gr38))))
    start <- min(unlist(start(gr38)))
    end <- max(unlist(end(gr38)))
    list(chr=chr,start=start,end=end,region=paste0(chr,":",start,"-",end))
  '
  export region38=4:72995642-74939286
  export CT=${GTEx}/Colon_Transverse.tsv.gz
  cat <(tabix -H ${CT} ${region38} | tail -n -1) \
      <(tabix ${CT} ${region38}) | \
  awk -vOFS="\t" '/ENSG00000163735/{
         if (NR==1) print "snpid","rsid","chr","pos","a1","a2","z";
         else {
           if ($15<$16) snpid="chr"$13":"$14"_"$15"_"$16; else snpid="chr"$13":"$14"_"$16"_"$15;
           z=$9/$10;
           print snpid, $18, $13, $14, $16, $15, z
         }
       }' | \
  gzip -f > ${INF}/work/CT.tsv.gz
  export WB=${GTEx}/Whole_Blood.tsv.gz
  cat <(tabix -H ${WB} ${region38} | tail -n -1) \
      <(tabix ${WB} ${region38}) | \
  awk -vOFS="\t" '/ENSG00000163735/{
         if (NR==1) print "snpid","rsid","chr","pos","a1","a2","z";
         else {
           if ($15<$16) snpid="chr"$13":"$14"_"$15"_"$16; else snpid="chr"$13":"$14"_"$16"_"$15;
           z=$9/$10;
           print snpid, $18, $13, $14, $16, $15, z
         }
       }' | \
  gzip -f > ${INF}/work/WB.tsv.gz
  Rscript -e '
    library(dplyr)
    library(gap)
    library(grid)
    INF <- Sys.getenv("INF")
    gsmr <- read.delim(file.path(INF,"mr/gsmr/","gsmr-efo-reduce.txt")) %>%
            left_join(pQTLdata::inf1[c("target.short","gene")],by=c("protein"="target.short")) %>%
            select(Disease,bxy,se,p,p_qtl,gene) %>%
            filter(grepl("CXCL5",gene)&grepl("Crohn\'s disease|Ulcerative colitis",Disease))
    pdf(file.path(INF,"SF-CXCL5-MR.pdf"),height=3,width=7)
    mr_forestplot(gsmr,colgap.forest.left="0.05cm", fontsize=14,
                  leftcols=c("studlab"), leftlabs=c("Disease"),
                  plotwidth="3inch", sm="OR", sortvar=dat[["bxy"]],
                  rightcols=c("effect","ci","pval"), rightlabs=c("OR","95%CI","P"),
                  digits=2, digits.pval=2, scientific.pval=TRUE,
                  common=FALSE, random=FALSE, print.I2=FALSE, print.pval.Q=FALSE, print.tau2=FALSE,
                  addrow=TRUE, backtransf=TRUE, spacing=1.6, at=c(0.7,0.8,0.9,1,1.1),xlim=c(0.7,1.1))
    grid::grid.text("GSMR results", 0.5, 0.9)
    dev.off()
    vars <- c("snpid","rsid","chr","pos","a1","a2","z")
    cxcl5 <- read.delim(file.path(INF,"work","CXCL5.tsv.gz")) %>% setNames(vars)
    uc <- read.delim(file.path(INF,"work","UC.tsv.gz")) %>% setNames(vars)
    cd <- read.delim(file.path(INF,"work","CD.tsv.gz")) %>% setNames(vars)
    eqtl <- read.delim(file.path(INF,"work","eQTLGen.tsv.gz")) %>% setNames(vars)
    cxcl5_z <- select(cxcl5,snpid,rsid,z,chr,pos) %>% mutate(marker=rsid,cxcl5=z)
    eqtl_z <- mutate(eqtl,eqtl=z) %>% select(snpid,eqtl)
    uc_z <- mutate(uc,uc=z) %>% select(snpid,uc)
    cd_z <- mutate(cd,cd=z) %>% select(snpid,cd)
    traits <- cxcl5_z %>% left_join(eqtl_z) %>% left_join(uc_z) %>% left_join(cd_z) %>%
              select(chr,pos,marker,cxcl5,eqtl,uc,cd)
    library(data.table)
    dt <- data.table(traits)
    dup <- dt[duplicated(marker), cbind(.SD[1], number = .N), by = marker] %>% pull(marker)
    traits <- filter(traits,!marker %in% dup)
    plink_bin <- "/rds/user/jhz22/hpc-work/bin/plink"
    chr <- Sys.getenv("chr")
    bfile <- file.path(INF,"INTERVAL","per_chr",paste0("interval.imputed.olink.chr_",chr))
    r <- ieugwasr::ld_matrix(select(traits,marker),with_alleles=TRUE,pop="EUR",bfile=bfile,plink_bin=plink_bin)
    rnames <- gsub("_[A-Z]*","",colnames(r))
    traits <- subset(traits,marker %in% rnames)
    rsids <- intersect(rnames,with(traits,marker))
    ld <- r
    colnames(ld) <- rownames(ld) <- rnames
    ld <- ld[rsids,rsids]
    d <- subset(traits,marker %in% rsids)
    z <- d[c("cxcl5","eqtl","uc","cd")] %>% setNames(c("CXCL5","Gene expression","Ulcerative colitis","Crohn's disease"))
    rownames(z) <- with(d,marker)
    library(gassocplot)
    pdf(file.path(INF,"CXCL5-gassoc.pdf"),height=20,width=8)
    sap <- stack_assoc_plot(d[c("marker","chr","pos")], z, ld, traits=names(z), ylab="-log10(P)", top.marker="rs450373",legend=TRUE)
    grid::grid.draw(sap)
    dev.off()
    system("qpdf CXCL5-gassoc.pdf --pages . 2 -- --replace-input")
    ct <- read.delim(file.path(INF,"work","CT.tsv.gz")) %>% setNames(vars) %>%
          mutate(SNP=rsid,CHR=chr,POS=pos) %>% select(SNP,CHR,POS,a1,a2,z)
    wb <- read.delim(file.path(INF,"work","WB.tsv.gz")) %>% setNames(vars) %>%
          mutate(SNP=rsid,CHR=chr,POS=pos) %>% select(SNP,CHR,POS,a1,a2,z)
    library(catalogueR)
    wb.lifted <- liftover(gwas_data=wb, build.conversion="hg38.to.hg19") %>% data.frame()
    ct.lifted <- liftover(gwas_data=ct, build.conversion="hg38.to.hg19") %>% data.frame()
    wb_z <- mutate(wb.lifted,wb=z,snpid=gap::chr_pos_a1_a2(CHR,start,a1,a2)) %>% select(snpid,wb)
    ct_z <- mutate(ct.lifted,ct=z,snpid=gap::chr_pos_a1_a2(CHR,start,a1,a2)) %>% select(snpid,ct)
    traits <- cxcl5_z %>% left_join(wb_z) %>% left_join(ct_z) %>% left_join(uc_z) %>% left_join(cd_z) %>%
              select(chr,pos,marker,cxcl5,ct,wb,uc,cd)
    library(data.table)
    dt <- data.table(traits)
    dup <- dt[duplicated(marker), cbind(.SD[1], number = .N), by = marker] %>% pull(marker)
    traits <- filter(traits,!marker %in% dup)
    plink_bin <- "/rds/user/jhz22/hpc-work/bin/plink"
    chr <- Sys.getenv("chr")
    bfile <- file.path(INF,"INTERVAL","per_chr",paste0("interval.imputed.olink.chr_",chr))
    r <- ieugwasr::ld_matrix(select(traits,marker),with_alleles=TRUE,pop="EUR",bfile=bfile,plink_bin=plink_bin)
    rnames <- gsub("_[A-Z]*","",colnames(r))
    traits <- subset(traits,marker %in% rnames)
    rsids <- intersect(rnames,with(traits,marker))
    ld <- r
    colnames(ld) <- rownames(ld) <- rnames
    ld <- ld[rsids,rsids]
    d <- subset(traits,marker %in% rsids)
    z <- d[c("cxcl5","wb","ct","uc","cd")] %>% setNames(c("CXCL5","Whole blood","Colon transverse","Ulcerative colitis","Crohn's disease"))
    rownames(z) <- with(d,marker)
    library(gassocplot)
    pdf(file.path(INF,"CXCL5-gassoc2.pdf"),height=20,width=8)
    sap <- stack_assoc_plot(d[c("marker","chr","pos")], z, ld, traits=names(z), ylab="-log10(P)", top.marker="rs450373",legend=TRUE)
    grid::grid.draw(sap)
    dev.off()
  '
}

function CXCL5()
{
  export eQTLGen=~/rds/public_databases/eQTLGen
  export cis_full=${eQTLGen}/cis-eQTLs_full_20180905.txt.gz
  export eQTLGen_tabix=${eQTLGen}/tabix
  export GTEx=~/rds/public_databases/GTEx/csv
  export chr=4
  export start=74861359
  export end=74864496
  export M=1e6
  Rscript -e '
    library(dplyr)
    library(gap)
    library(grid)
    INF <- Sys.getenv("INF")
    vars <- c("snpid","rsid","chr","pos","a1","a2","z")
    cxcl5 <- read.delim(file.path(INF,"work","CXCL5.tsv.gz")) %>% setNames(vars)
    uc <- read.delim(file.path(INF,"work","UC.tsv.gz")) %>% setNames(vars)
    cd <- read.delim(file.path(INF,"work","CD.tsv.gz")) %>% setNames(vars)
    eqtl <- read.delim(file.path(INF,"work","eQTLGen.tsv.gz")) %>% setNames(vars)
    cxcl5_z <- select(cxcl5,snpid,rsid,z,chr,pos) %>% mutate(marker=rsid,cxcl5=z)
    eqtl_z <- mutate(eqtl,eqtl=z) %>% select(snpid,eqtl)
    ct <- read.delim(file.path(INF,"work","CT.tsv.gz")) %>% setNames(vars) %>%
          mutate(SNP=rsid,CHR=chr,POS=pos) %>% select(SNP,CHR,POS,a1,a2,z)
    wb <- read.delim(file.path(INF,"work","WB.tsv.gz")) %>% setNames(vars) %>%
          mutate(SNP=rsid,CHR=chr,POS=pos) %>% select(SNP,CHR,POS,a1,a2,z)
    uc_z <- mutate(uc,uc=z) %>% select(snpid,uc)
    cd_z <- mutate(cd,cd=z) %>% select(snpid,cd)
    library(catalogueR)
    wb.lifted <- liftover(gwas_data=wb, build.conversion="hg38.to.hg19") %>% data.frame()
    ct.lifted <- liftover(gwas_data=ct, build.conversion="hg38.to.hg19") %>% data.frame()
    wb_z <- mutate(wb.lifted,wb=z,snpid=gap::chr_pos_a1_a2(CHR,start,a1,a2)) %>% select(snpid,wb)
    ct_z <- mutate(ct.lifted,ct=z,snpid=gap::chr_pos_a1_a2(CHR,start,a1,a2)) %>% select(snpid,ct)
    start <- Sys.getenv("start") %>% as.numeric()
    end <- Sys.getenv("end") %>% as.numeric()
    traits <- cxcl5_z %>% left_join(eqtl_z) %>% left_join(ct_z) %>% left_join(uc_z) %>% left_join(cd_z) %>%
              select(chr,pos,marker,cxcl5,eqtl,ct,uc,cd) %>%
              filter(pos>=start-300000 & pos<=end+300000)
    library(data.table)
    dt <- data.table(traits)
    dup <- dt[duplicated(marker), cbind(.SD[1], number = .N), by = marker] %>% pull(marker)
    traits <- filter(traits,!marker %in% dup)
    plink_bin <- "/rds/user/jhz22/hpc-work/bin/plink"
    chr <- Sys.getenv("chr")
    bfile <- file.path(INF,"INTERVAL","per_chr",paste0("interval.imputed.olink.chr_",chr))
    r <- ieugwasr::ld_matrix(select(traits,marker),with_alleles=TRUE,pop="EUR",bfile=bfile,plink_bin=plink_bin)
    rnames <- gsub("_[A-Z]*","",colnames(r))
    traits <- subset(traits,marker %in% rnames)
    rsids <- intersect(rnames,with(traits,marker))
    ld <- r
    colnames(ld) <- rownames(ld) <- rnames
    ld <- ld[rsids,rsids]
    d <- subset(traits,marker %in% rsids)
    z <- d[c("cxcl5","eqtl","ct","uc","cd")] %>% setNames(c("","","","",""))
    rownames(z) <- with(d,marker)
    library(gassocplot)
    pdf(file.path(INF,"SF-CXCL5-WB-CT-UC-CD.pdf"),height=20,width=6)
    sap <- stack_assoc_plot(d[c("marker","chr","pos")], z, ld, traits=names(z), ylab="-log10(P)", top.marker="rs450373",legend=FALSE)
    grid::grid.draw(sap)
    dev.off()
  '
  # setNames(c("CXCL5","Whole blood","Colon transverse","Ulcerative colitis","Crohn's disease"))
}

function SCALLOP_MRC()
{
  module load ceuadmin/stata
  cd ~/INF
  stata <<\ \ END
  insheet using INTERVAL/o5000-inf1-outlier_in-r2.sample, delim(" ")
  format id* %15.0g
  gen AGE=real(age)
  sum AGE
  END
  cd ~/COVID-19/HGI
  stata <<\ \ END
  insheet using "20210317/INTERVALdata_17MAR2021.csv", case clear
  sort identifier
  keep identifier sexPulse agePulse ht_bl wt_bl
  rename sexPulse sex
  save a
  insheet using "20210317/INTERVAL_OmicsMap_20210317.csv", case clear
  keep identifier Affymetrix_QC_bl Affymetrix_gwasQC_bl Olink_inf_gwasQC_24m
  format Affymetrix_QC_bl %15.0g
  format Affymetrix_gwasQC_bl %15.0g
  format Olink_inf_gwasQC_24m %15.0g
  merge 1:1 identifier using a, gen(omics_data)
  l Affymetrix_gwasQC_bl sex age ht_bl wt_bl if wt_bl!=. & wt_bl>200, ab(15) linesize(150)
  replace wt_bl=. if wt_bl==777
  gen bmi=wt_bl/ht_bl/ht_bl
  sum bmi if Olink_inf_gwasQC_24m!=.
  l Affymetrix_gwasQC_bl sex age ht_bl wt_bl bmi if bmi!=. & bmi>50, ab(15) linesize(150)
  erase a.dta
  END
  cd -
}

function CXCL5_wb_ct()
{
  cat <(head -1 ${INF}/work/INF1.METAL) <(grep -w rs450373 ${INF}/work/INF1.METAL) | cut -f2,6,7,9,10
# beta se ref alt rsid
  cat <(gunzip -c ~/rds/public_databases/GTEx/csv/Stomach.tsv.gz | head -1) \
      <(zgrep -w rs450373 ~/rds/public_databases/GTEx/csv/Colon_Transverse.tsv.gz | grep ENSG00000163735) \
      <(zgrep -w rs450373 ~/rds/public_databases/GTEx/csv/Whole_Blood.tsv.gz | grep ENSG00000163735) | cut -f9,10,15,16,18
  R --no-save <<\ \ END
  CXCL5 <- '
       "Plasma protein" -0.5115   0.0183
     "Whole blood mRNA" -0.546318 0.0542992
           "Colon mRNA" -0.407191 0.0858433
  '
  CXCL5 <- as.data.frame(scan(file=textConnection(CXCL5),what=list("",0,0))) %>%
           setNames(c("outcome","Effect","StdErr")) %>%
           mutate(outcome=gsub("\\b(^[a-z])","\\U\\1",outcome,perl=TRUE))
  library(gap)
  pdf(file.path(INF,"SF-CXCL5-WB-CT.pdf"),height=3,width=9)
  mr_forestplot(CXCL5,colgap.forest.left="0.05cm", fontsize=14,
                leftcols=c("studlab"), leftlabs=c("CXCL5 measurement"),
                plotwidth="3inch", sm="OR",
                rightcols=c("effect","ci","pval"), rightlabs=c("OR","95%CI","P"),
                digits=2, digits.pval=2, scientific.pval=TRUE,
                common=FALSE, random=FALSE, print.I2=FALSE, print.pval.Q=FALSE, print.tau2=FALSE,
                addrow=TRUE, backtransf=TRUE, at=5:11/10, spacing=1.5, xlim=c(0.5,1.1))
  dev.off()
  END
}
# rs450373 (A/G)
#           Source    Effect    StdErr
#            CXCL5 -0.5115   0.0183
#      Whole blood -0.546318 0.0542992
# Colon transverse -0.407191 0.0858433

function CXCL4_rs450373()
# CXCL5___P42830
{
  cd ${INF}/CXCL5
  plink --bfile ${INF}/INTERVAL/per_chr/interval.imputed.olink.chr_4 --snp rs450373 --recode --out rs450374
  plink --bfile ${INF}/INTERVAL/per_chr/interval.imputed.olink.chr_4 --snp rs450373 --recode A --out rs450374
  sed 's/ID_1/FID/;s/ID_2/IID/;2d' ${INF}/INTERVAL/o5000-inf1-outlier_in-r2.sample > rs450374.pheno
  plink --bfile ${INF}/INTERVAL/per_chr/interval.imputed.olink.chr_4 \
        --pheno rs450374.pheno --pheno-name CXCL5___P42830 --snp rs450373 --recode A include-alt --out rs450374
  head rs450374.raw rs450374.ped
  Rscript -e '
    p <- read.delim("rs450374.pheno",sep=" ")
    head(p[c("FID","IID","CXCL5___P42830")])
    library(dplyr)
    r <- read.delim("rs450374.raw",na.strings="-9", sep=" ") %>%
         rename(rs450374=rs450373_G..A.,CXCL5=PHENOTYPE) %>%
         filter(rs450374!="NA") %>%
         mutate(Genotype=case_when(rs450374==0 ~ "AA", rs450374==1 ~ "AG", rs450374==2 ~ "GG", TRUE ~ "NA"))
    s <- group_by(r,rs450374) %>%
         summarise(
             Mean=mean(CXCL5,na.rm=TRUE)%>%round(digits=2),
             SD=sd(CXCL5,na.rm=TRUE)%>%round(digits=2),
             N=n()) %>%
         data.frame
    row.names(s) <- c("AA","AG","GG")
    print(s)
    require(ggplot2)
    require(ggpubr)
    v <- ggplot(r, aes(x=Genotype, y=CXCL5, fill=Genotype)) +
                geom_violin() +
                geom_boxplot(width=0.1) +
                xlab("rs450374 genotype") +
                ylab("CXCL5 level") +
                theme_bw(base_rect_size=0) + theme(legend.position = "none")
    m <- ggtexttable(s %>%
                     mutate(Genotype=case_when(rs450374==0 ~ "AA", rs450374==1 ~ "AG", rs450374==2 ~ "GG", TRUE ~ "NA"),
                            N=format(N,big.mark=",")) %>%
                     select(Genotype, Mean, SD, N),
                     rows = NULL, theme = ttheme("lBlueWhite")) + theme_bw(base_rect_size=0)
    p <- ggarrange(v,m,ncol=1,nrow=2)
    ggsave(p,file=file.path(INF,"CXCL5","SF-CXCL5-rs450374.png"),dpi=300,height=4,width=5,units="in")
  '
  cd -
}

function tbi()
# manageable sort
{
  ls *-1.tbl.gz | sed 's/.gz//' | grep -v BDNF | parallel -j10 -C' ' '
  gunzip {}.gz;
# cat <(head -1 {}) <(sed "1d" {} | sort -k1,1n -k2,2n) | bgzip -f > {}.gz
  (
    head -1 {}
    for chr in {1..22}
    do
      awk -vchr=${chr} "\$1==chr" {} | sort -k2,2n
    done
  ) | bgzip -f > {}.gz
  tabix -f -S1 -s1 -b2 -e2 {}.gz
  '
}
