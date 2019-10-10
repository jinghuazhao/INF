# 10-10-2019 JHZ

rsync -avrzP . jinhua@10.44.124.6:/home/jinhua/INF \
      --exclude=.git --exclude=sumstats --exclude=METAL --exclude=work
