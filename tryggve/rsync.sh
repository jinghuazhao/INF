# 10-5-2020 JHZ

rsync -avrzP . jinhua@10.44.124.6:/home/jinhua/INF \
      --exclude=.git --exclude=ds --exclude=sumstats --exclude=METAL --exclude=work
