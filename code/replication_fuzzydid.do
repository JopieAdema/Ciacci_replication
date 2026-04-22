*Standalone runner for Table C1 (Fuzzy DiD). Split out of replication.do
*because `fuzzydid' with breps=100 takes several minutes and is rarely re-run.
*Exports: console/log output only (no polished .tex export yet).
*Requires: `ssc install fuzzydid'.

*Block 0: Setup ----
{
    clear all
    cd "C:\Users\c4041171\Dropbox\Ciacci replication\replication_220426\"
    set seed 1234
    * Number of bootstrap replications for `fuzzydid' standard errors.
    * Production value ~100; lower for faster iteration.
    global B_REPS = 99
}


*Table C1: Fuzzy DiD (Wald-DiD with and without time correction) ----
{
    cap which fuzzydid
    if _rc {
        di as error "Table C1 skipped: `fuzzydid' not installed. Run: ssc install fuzzydid"
        exit
    }

    use ".\data\Sweden_tables.dta", clear
    cap drop sum_sex_pur
    bys region (ym): gen sum_sex_pur = sum(sex_purchase_abs_)
    gen d_sum_sex_pur      = (sum_sex_pur>0)
    gen d_sex_purchase_abs_= (sex_purchase_abs_>0)
    bys region (ym): gen d_sum_sex_pur1 = sum(d_sex_purchase_abs_)
    xtset regionc ym
    gen fd_sum_sex_pur = d.sum_sex_pur
    gen sample1 = 0 if missing(fd_sum_sex_pur)
    replace  sample1 = 1 if sample1==.
    gen D = d_sum_sex_pur1
    gen T = ym
    gen G = regionc
    bys G T: egen mean_D = mean(D)
    by G: gen lag_mean_D = mean_D[_n-1] if G==G[_n-1] & T-1==T[_n-1]
    gen GT = sign(mean_D - lag_mean_D) if sample1==1
    gen GTplus1 = GT[_n+1] if G==G[_n+1] & T+1==T[_n+1]

    di _n as res "TableC1 col 1 & 4 (monthly, Wald-DiD + TC)"
    fuzzydid lrape GT GTplus1 year D, did tc breps($B_REPS) cluster(regionc) ///
        continuous(police) qualitative(regionc year month)
    matrix m_b_1  = e(b_LATE)
    matrix m_se_1 = e(se_LATE)
    scalar c1_b_DID  = m_b_1[1,1]
    scalar c1_se_DID = m_se_1[1,1]
    scalar c1_b_TC   = m_b_1[2,1]
    scalar c1_se_TC  = m_se_1[2,1]
    scalar c1_N      = e(N)

    di _n as res "TableC1 col 2 & 5 (monthly, time-as-T)"
    fuzzydid lrape GT GTplus1 T    D, did tc breps($B_REPS) cluster(regionc) ///
        continuous(police) qualitative(regionc year month)
    matrix m_b_2  = e(b_LATE)
    matrix m_se_2 = e(se_LATE)
    scalar c2_b_DID  = m_b_2[1,1]
    scalar c2_se_DID = m_se_2[1,1]
    scalar c2_b_TC   = m_b_2[2,1]
    scalar c2_se_TC  = m_se_2[2,1]
    scalar c2_N      = e(N)

    *annual collapse for cols 3 & 6
    use ".\data\Sweden_tables.dta", clear
    collapse (sum) rape sex_purchase_abs_ (mean) police, by(regionc year)
    gen lrape = log(rape+1)
    cap drop sum_sex_pur
    bys region (year): gen sum_sex_pur = sum(sex_purchase_abs_)
    gen d_sum_sex_pur      = (sum_sex_pur>0)
    gen d_sex_purchase_abs_= (sex_purchase_abs_>0)
    bys region (year): gen d_sum_sex_pur1 = sum(d_sex_purchase_abs_)
    xtset regionc year
    gen fd_sum_sex_pur = d.sum_sex_pur
    gen sample1 = 0 if missing(fd_sum_sex_pur)
    replace sample1 = 1 if sample1==.
    gen D = d_sum_sex_pur1
    gen T = year
    gen G = regionc
    bys G T: egen mean_D = mean(D)
    by G: gen lag_mean_D = mean_D[_n-1] if G==G[_n-1] & T-1==T[_n-1]
    gen GT      = sign(mean_D - lag_mean_D) if sample1==1
    gen GTplus1 = GT[_n+1] if G==G[_n+1] & T+1==T[_n+1]
    di _n as res "TableC1 col 3 & 6 (annual)"
    fuzzydid lrape GT GTplus1 year D, did tc breps($B_REPS) cluster(regionc) ///
        continuous(police) qualitative(regionc year)
    matrix m_b_3  = e(b_LATE)
    matrix m_se_3 = e(se_LATE)
    scalar c3_b_DID  = m_b_3[1,1]
    scalar c3_se_DID = m_se_3[1,1]
    scalar c3_b_TC   = m_b_3[2,1]
    scalar c3_se_TC  = m_se_3[2,1]
    scalar c3_N      = e(N)
}


*--- write tC1.tex -------------------------------------------------------
* Layout harmonised with main tables (t2.tex etc.):
*   columns (1)-(3) = Wald-DiD for (monthly fine-time, monthly calendar-time, annual)
*   columns (4)-(6) = Wald-DiD + time correction, same three specs.
{
    capture file close tC1
    file open tC1 using "./output/tC1.tex", write replace

    foreach est in DID TC {
        foreach c in 1 2 3 {
            local b   = c`c'_b_`est'
            local se  = c`c'_se_`est'
            local p   = 2*(1 - normal(abs(`b'/`se')))
            local stars = cond(`p'<.01,"***",cond(`p'<.05,"**",cond(`p'<.10,"*","")))
            local b_`c'_`est'_s  = string(`b',  "%9.3f") + "`stars'"
            local se_`c'_`est'_s = "(" + string(`se',"%9.3f") + ")"
            local p_`c'_`est'_s  = "\textit{" + string(`p',"%9.3f") + "}"
        }
    }
    local N1_s = string(c1_N,"%9.0f")
    local N2_s = string(c2_N,"%9.0f")
    local N3_s = string(c3_N,"%9.0f")

    file write tC1 "\begin{tabular}{lcccccc}" _n
    file write tC1 " & (1) & (2) & (3) & (4) & (5) & (6) \\" _n
    file write tC1 " & \multicolumn{3}{c}{Wald-DiD} & \multicolumn{3}{c}{Wald-DiD + time correction} \\" _n
    file write tC1 " & Monthly & Monthly & Annual & Monthly & Monthly & Annual \\" _n
    file write tC1 " & (fine time) & (calendar) & & (fine time) & (calendar) & \\" _n
    file write tC1 "Reform = 1 & `b_1_DID_s' & `b_2_DID_s' & `b_3_DID_s' & `b_1_TC_s' & `b_2_TC_s' & `b_3_TC_s' \\" _n
    file write tC1 " & `se_1_DID_s' & `se_2_DID_s' & `se_3_DID_s' & `se_1_TC_s' & `se_2_TC_s' & `se_3_TC_s' \\" _n
    file write tC1 "p-value & `p_1_DID_s' & `p_2_DID_s' & `p_3_DID_s' & `p_1_TC_s' & `p_2_TC_s' & `p_3_TC_s' \\" _n
    file write tC1 "Observations & `N1_s' & `N2_s' & `N3_s' & `N1_s' & `N2_s' & `N3_s' \\" _n
    file write tC1 "\end{tabular}" _n
    file close tC1
}
