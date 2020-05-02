## Joint/Condtional analysis and finemapping

The INTERVAL data is used as reference panel. The logic of this specific directory is a simple solution of the delimma that the reference data, possibly
like others, uses reference sequence ID (rsid) whenever possible. However, during meta-analysis the practice of using rsid is undesirable so SNPID, i.e.,
chr:pos_A1_A2, (A1<=A2) is necessary.

After rather a long and laborious process involving many software, it turned out a simple wayout is to obtain sentinels using SNPID but return to rsid at
this stage and forward. The implementation here reflects this. The file INTERVAL.rsid contains SNPID-rsid mapping and could be generated from programs
such as qctool/bgenix/plink.

### A summary

File | Function
-----|------------------------------
NLRP2.sh | the exclusion list
ma.sh | INF1 sumstats
INTERVAL-ma.sh | INTERVAL sumstats
prune.sh | pruning
slct.sh  | GCTA --cojo-slct
finemap.sh | finemapping
st.sh | batch command file

### Steps

st.sh executes the following elements,

0. NLRP2.sh
1. prune.sh
2. ma.sh
3. slct.sh
4. finemap.sh

*Date last changed:* **2/5/2020**
