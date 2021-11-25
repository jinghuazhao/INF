# 9-4-2019 JHZ

source tryggve/analysis.ini

echo "--> 1000Genomes reference data"

# individual genotypes
export EUR=/services/tools/locuszoom/1.4/data/1000G/genotypes/2014-10-14/EUR

seq 22 | \
parallel 'plink --bfile $EUR/chr{} --recode vcf bgz --out chr{}'
seq 22 | \
parallel 'tabix -f -p vcf chr{}.vcf.gz'
seq 22 | \
awk -vp=$PWD '{print p "/chr" $1 ".vcf.gz"}' > EUR.list
vcf-concat --files EUR.list | \
bcftools norm -d both - | \
bcftools convert - -O z -o EUR.vcf.gz
plink --vcf EUR.vcf.gz --make-bed --out EUR
awk -vFS="\t" '
{
  CHR=$1
  POS=$4
  a1=$5
  a2=$6
  if (a1>a2) snpid="chr" CHR ":" POS "_" a2 "_" a1;
  else snpid="chr" CHR ":" POS "_" a1 "_" a2
  print snpid, $2

# 1 CHR
# 2 SNPID
# 3 0
# 4 POS
# 5 A1
# 6 A2

}' EUR.bim > EUR.snpid

plink --bfile EUR --update-name EUR.snpid 1 2 --make-bed --out EUR1KG
plink --bfile EUR1KG --chr 6 --from-mb 25 --to-mb 35 --make-bed --out MHC
cut -f2 MHC.bim > MHC.snpid

plink --bfile EUR1KG \
      --exclude range tryggve/high-LD-regions-hg19.txt \
      --make-bed --out EUR

echo "--> SNP/gene databases and indices"

export hg19=/services/tools/locuszoom/1.4/data/database/locuszoom_hg19.db
sqlite3 $hg19 <<END
-- locuszoom 1.4 tables and indices
.tables
.indices
.separator "\t"
.header on
.output gencode.tsv
select * from gencode;
.output recomb_rate.tsv
select * from recomb_rate;
.output refFlat.tsv
select * from refFlat;
.output refsnp_trans.tsv
select * from refsnp_trans;
.output snp_pos.tsv
select * from snp_pos;
.output snp_set.tsv
select * from snp_set;
END

# to establish correspondence between snpid and rsid

gzip -f snp_pos.tsv > snp_pos.tsv.gz

zcat snp_pos.tsv.gz | \
awk 'NR>1' | \
awk '{chrpos="chr" $2 ":" $3; print $1,chrpos}' | \
sort -k2,2 > snp_pos
sort -k2,2 EUR.snpid | \
join -j2 -a2 - snp_pos | \
sort -k1,1 > EUR.rsid
