# Source http://users.nber.org/~rdehejia/data/nswdata2.html


lalonde <- haven::read_dta("data-raw/lalonde.dta")
devtools::use_data(lalonde)
