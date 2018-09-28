# Source http://users.nber.org/~rdehejia/data/nswdata2.html


nsw <- foreign::read.dta("data-raw/nsw.dta")
devtools::use_data(nsw)
