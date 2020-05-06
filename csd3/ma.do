// 6-5-2020 JHZ

local p : env p
insheet x1-x25 using `p'.ab, clear
save `p'.dta, replace
use `p'.dta if x13 < -323
gen index=.
if _N>0 {
  sort x1 x25 x13
  by x1 x25: replace index=_n
}
append using `p'.dta
drop if index==. & x13 < -323
replace index=0 if index==.
sort x25 x4
outsheet x1-x25 index using `p'.txt, replace noquote

// x1 = chr1
// x4 = SNPID
// x13= log10(p)

