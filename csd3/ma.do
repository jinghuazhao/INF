// 6-11-2019 JHZ

local p : env p
insheet x1-x25 using `p'.ab, clear
save `p'.dta, replace
use `p'.dta if x13 < -323
sort x1 x25 x13
by x1 x25: gen index=_n
save `p'0.dta, replace
use `p'.dta if x13 >= -323
append using `p'0.dta
replace index=0 if index==.
sort x25 x4
outsheet x1-x25 index using `p'.txt, replace noquote
