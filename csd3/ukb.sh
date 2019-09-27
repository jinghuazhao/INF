# 27-9-2019 JHZ

export INF=/rds/project/jmmh2/rds-jmmh2-projects/olink_proteomics/scallop/INF
export srcdir=${INF}/ukb

function individual()
{
  sed '1d' work/INF1.merge | \
  sed 's/chr//g' | \
  awk -v srcdir=${srcdir} -v flanking=1e6 -v INF=${INF} '
  {
    if ($2 >= flanking) start=$2-flanking;
    else start = 0;
    end = $3 + flanking
    range = $1 ":" start "-" end;
    if($1<=9) range0=0 range;
    else range0=range
    cmd=sprintf("qctool -g %s/ukb_imp_chr%d_v3.bgen -incl-range %s -ofiletype bgen -og ukb/%s-chr%s.bgen", srcdir, $1, range0, $5, $6)
    print cmd
  }'
}

function combined ()
{
  sed '1d' work/INF1.merge | \
  sortBed -i | \
  mergeBed -i - -d 1000000 | \
  sed 's/chr//g' | \
  awk -v srcdir=${srcdir} -v INF=${INF} '
  {
    range = $1 ":" $2 "-" $3;
    if($1<=9) range0=0 range;
    else range0=range
    cmd=sprintf("qctool -g %s/ukb_imp_chr%d_v3.bgen -incl-range %s -ofiletype bgen -og ukb/chr%s.bgen", srcdir, $1, range0, range)
    print cmd
  }'
}
