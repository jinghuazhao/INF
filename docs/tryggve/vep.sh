# 5-11-2018 JHZ

function vep()
# a1/a2 assertion is necessary
{
  (
  echo -e "##fileformat=VCFv4.0"
  echo -e "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO"
  awk -vOFS="\t" '(NR>1){
    snpid=$1
    rsid=$2
    split(snpid,a,":")
    chr=a[1]
    split(a[2],b,"_")
    pos=b[1]
    a1=b[2]
    a2=b[3]
    QUAL=".";FILTER=".";INFO="."
    print chr,pos,pos,a1,a2,QUAL,FILTER,INFO
  }' INTERVAL.jma.dat | \
  uniq > INTERVAL.vcf
  )

  module load perl/5.20.0 vep/88
  ln -sf /usr/local/Cluster-Apps/vep/release-88/cache $/HOME/.vep
# export vepfile=/usr/local/Cluster-Apps/vep/release-88/examples/homo_sapiens_GRCh37.vcf
  export vepfile=INTERVAL.vcf
  vep -i $vepfile --assembly GRCh37 -o $(basename $vepfile .vcf).out --force_overwrite -offline
}
