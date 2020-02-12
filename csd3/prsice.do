/* 12-2-2020 JHZ*/

foreach v in "chd" "cv" {
  insheet using work/crp-`v'.best, case clear delim(" ")
  sort FID
  merge 1:1 FID using work/ukb
  logit `v' sex ages PRS
  stset ages, failure(`v')
  stcox sex PRS if `v'!=.
}
