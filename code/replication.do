*Block 0: Setup ----
{
    clear all
    cd "C:\Users\c4041171\Dropbox\Ciacci replication\replication_220426\"
    set seed 1234
    * Number of cluster-bootstrap replications for Table 2 cols 3 & 6.
    * Production value ~99-999; lower for faster iteration.
    global B_REPS = 99

    * ensure required SSC packages are present
    foreach pkg in rdrobust reghdfe outreg2 tuples fuzzydid {
        cap which `pkg'
        if _rc cap ssc install `pkg', replace
    }
}
*Block 1: Data prep (Sweden_tables.dta) ----
* Centralised cleaning used by Table 2 and the bootstrap.
* Table 3 reloads raw data because its sample (year<=2005) and season coding differ.
{
    * fail informatively if any required raw data file is missing
    foreach f in "Sweden_tables.dta" "Sweden_tables3.dta" ///
                 "yearly national data rape reports.dta" ///
                 "annual data nordic countries.dta" {
        cap confirm file "./data/`f'"
        if _rc {
            di as error "Raw data file missing: ./data/`f'"
            di as error "Please place all .dta files in ./data/ before running."
            exit 601
        }
    }

    use ".\data\Sweden_tables.dta", clear

    * running variable + treatment
    gen running         = ym - 468
    gen Treat           = (running >= 0)
    gen Treat_running   = Treat * running
    gen running_2       = running^2
    gen Treat_running_2 = Treat * running_2

    * cluster variable
    cap drop regionc_m
    egen regionc_m = group(regionc month)

    * fixed-effect dummies
    qui tab month,   gen(mdum)
    qui tab year,    gen(ydum)
    qui tab regionc, gen(cdum)

    * first differences
    xtset regionc ym
    gen dlrape  = d.lrape
    gen dlrape2 = d.lrape2
    gen dpolice = d.police

    * inverse-hyperbolic sine of rapes (for Tables B1/B2 col 2)
    gen ihs_rape  = asinh(rape)
    gen dihs_rape = d.ihs_rape

    * month-residualised log rape (full sample + restricted to year<2001)
    qui reg lrape i.month
    predict lrape_mres,   resid
    qui reg lrape i.month if year<2001
    predict lrape_mres_r, resid

    tempfile prepped
    save `prepped'
}

*==============================================================================
* MAIN PAPER
*==============================================================================

*Figure 1: yearly national rape reports ----
{
    use ".\data\yearly national data rape reports.dta", clear
    twoway (scatter reports year), ///
        ylabel(0(1000)10000, angle(horizontal)) ///
        xline(1997.5 1998.5 2004.5 2017.5) ///
        xlabel(1991(1)2022, angle(forty_five)) ///
        scheme(s1mono)
    graph export ".\output\F1.pdf", replace
}
*Figure 2: from Ciacci (2025) -- external, nothing to run ----
*Figure 3: Nordic countries comparison ----
{
    use ".\data\annual data nordic countries.dta", clear
    foreach c in dk swe fi no {
        gen temp = `c' if year==1998
        egen `c'_99 = max(temp)
        drop temp
    }
    gen dk_rel = dk/dk_99
    gen se_rel = swe/swe_99
    gen fi_rel = fi/fi_99
    gen no_rel = no/no_99

    twoway ///
        (connected dk_rel year, lpattern(shortdash_dot) msymbol(triangle) lcolor(black) mcolor(black) msize(vsmall)) ///
        (connected fi_rel year, lpattern(dash)          msymbol(diamond)  lcolor(gs6)   mcolor(gs6)   msize(vsmall)) ///
        (connected no_rel year, lpattern(dot)           msymbol(square)   lcolor(gs6)   mcolor(gs6)   msize(vsmall)) ///
        (connected se_rel year, lpattern(solid)         msymbol(circle)   lcolor(black) mcolor(black)) ///
        , xline(1998.5, lcolor(red)   lpattern(solid)) ///
          xline(1997.5, lcolor(black) lpattern(dash))  ///
          xline(2004.5, lcolor(black) lpattern(dash))  ///
          ylabel(0(0.5)3, angle(horizontal) grid) ///
          xlabel(1990(2)2010, angle(forty_five)) ///
          ytitle("Reported rapes per 100,000 population (1998 = 1)") xtitle("") ///
          legend(order(1 "Denmark" 2 "Finland" 3 "Norway" 4 "Sweden") position(11) ring(0) cols(1)) ///
          text(2.7 1998.7 "Criminalization" "of buying sex" "in Sweden", placement(e) size(small) color(red)) ///
          text(0.3 2004.7 "Definition of rape" "expanded" "in Sweden",  placement(e) size(small) color(black)) ///
          text(0.3 1993.5 "Definition of rape" "expanded" "in Sweden",  placement(e) size(small) color(black)) ///
          scheme(s1mono) graphregion(lcolor(none)) plotregion(lstyle(none))
    graph export ".\output\F3.pdf", replace
}
*Table 1: text-only, nothing to run ----
*Table 2: RDiT estimates of reform on rape ----
{
    *--- Cols 1-3 (restricted sample, year<2001) ---
    use `prepped', clear
    preserve
    keep if year<2001

    qui rdbwselect lrape ym, c(468)
    scalar t2_bw_c1       = e(b_mserd)
    global lrape_b_l_r    = 468-e(b_mserd)
    global lrape_b_u_r    = 468+e(b_mserd)
    global lrape_bwidth_r "if ym>${lrape_b_l_r} & ym<${lrape_b_u_r}"

    *col 1
    qui reg lrape Treat Treat_running Treat_running_2 running running_2 police ///
        i.regionc i.year i.month ${lrape_bwidth_r}, cl(regionc_m)
    scalar t2_b_c1  = _b[Treat]
    scalar t2_se_c1 = _se[Treat]
    scalar t2_p_c1  = 2*ttail(e(df_r), abs(_b[Treat]/_se[Treat]))
    qui count if Treat==0 & ym>${lrape_b_l_r} & ym<${lrape_b_u_r} & e(sample)
    scalar t2_nl_c1 = r(N)
    qui count if Treat==1 & ym>${lrape_b_l_r} & ym<${lrape_b_u_r} & e(sample)
    scalar t2_nr_c1 = r(N)

    *col 3
    qui rdrobust lrape_mres_r running, covs(police cdum2-cdum21) ///
        vce(cluster regionc_m) p(2) kernel(uniform)
    scalar t2_b_c3  = e(tau_cl)
    scalar t2_se_c3 = e(se_tau_cl)
    scalar t2_p_c3  = e(pv_rb)     // analytical p (will be overwritten by bootstrap below)
    scalar t2_bw_c3 = e(h_l)
    scalar t2_nl_c3 = e(N_h_l)
    scalar t2_nr_c3 = e(N_h_r)

    restore

    *--- Cols 4-6 (full sample) ---
    qui rdbwselect lrape ym, c(468)
    scalar t2_bw_c4    = e(b_mserd)
    global lrape_b_l   = 468-e(b_mserd)
    global lrape_b_u   = 468+e(b_mserd)
    global lrape_bwidth "if ym>${lrape_b_l} & ym<${lrape_b_u}"

    *col 4
    qui reg lrape Treat Treat_running Treat_running_2 running running_2 police ///
        i.regionc i.year i.month ${lrape_bwidth}, cl(regionc_m)
    scalar t2_b_c4  = _b[Treat]
    scalar t2_se_c4 = _se[Treat]
    scalar t2_p_c4  = 2*ttail(e(df_r), abs(_b[Treat]/_se[Treat]))
    qui count if Treat==0 & ym>${lrape_b_l} & ym<${lrape_b_u} & e(sample)
    scalar t2_nl_c4 = r(N)
    qui count if Treat==1 & ym>${lrape_b_l} & ym<${lrape_b_u} & e(sample)
    scalar t2_nr_c4 = r(N)

    *col 6
    qui rdrobust lrape_mres running, covs(police cdum2-cdum21) ///
        vce(cluster regionc_m) p(2) kernel(uniform)
    scalar t2_b_c6  = e(tau_cl)
    scalar t2_se_c6 = e(se_tau_cl)
    scalar t2_p_c6  = e(pv_rb)     // analytical p (will be overwritten by bootstrap below)
    scalar t2_bw_c6 = e(h_l)
    scalar t2_nl_c6 = e(N_h_l)
    scalar t2_nr_c6 = e(N_h_r)
}
*Block: cluster bootstrap for Table 2, cols 3 & 6 ----
* Resample clusters (regionc_m) with replacement, re-residualise lrape on month FE
* within each draw, refit rdrobust, and compute bootstrap SE and p-value.
{
    local B = $B_REPS

    *--- Col 3 bootstrap (restricted sample) ---
    use `prepped', clear
    keep if year<2001
    tempfile data_restr
    save `data_restr'

    qui levelsof regionc_m, local(clusters)
    local n_clusters : word count `clusters'

    matrix bs_c3 = J(`B', 1, .)
    forval b = 1/`B' {
        clear
        set obs `n_clusters'
        gen double regionc_m = .
        gen int    _rep      = .
        forval i = 1/`n_clusters' {
            local draw : word `= ceil(runiform() * `n_clusters')' of `clusters'
            qui replace regionc_m = `draw' in `i'
            qui replace _rep      = `i'    in `i'
        }
        tempfile keys
        qui save `keys'

        use `data_restr', clear
        qui joinby regionc_m using `keys'
        cap drop lrape_mres_b
        qui reg lrape i.month
        qui predict lrape_mres_b, resid
        qui rdrobust lrape_mres_b running, covs(police cdum2-cdum21) ///
            vce(cluster regionc_m) p(2) kernel(uniform)
        matrix bs_c3[`b', 1] = e(tau_cl)
    }
    clear
    matrix colnames bs_c3 = tau
    svmat bs_c3, names(col)
    qui summ tau
    scalar t2_se_c3_boot = r(sd)
    scalar t2_p_c3_boot  = 2*(1 - normal(abs(t2_b_c3/t2_se_c3_boot)))

    *--- Col 6 bootstrap (full sample) ---
    use `prepped', clear
    tempfile data_full
    save `data_full'

    qui levelsof regionc_m, local(clusters)
    local n_clusters : word count `clusters'

    matrix bs_c6 = J(`B', 1, .)
    forval b = 1/`B' {
        clear
        set obs `n_clusters'
        gen double regionc_m = .
        gen int    _rep      = .
        forval i = 1/`n_clusters' {
            local draw : word `= ceil(runiform() * `n_clusters')' of `clusters'
            qui replace regionc_m = `draw' in `i'
            qui replace _rep      = `i'    in `i'
        }
        tempfile keys
        qui save `keys'

        use `data_full', clear
        qui joinby regionc_m using `keys'
        cap drop lrape_mres_b
        qui reg lrape i.month
        qui predict lrape_mres_b, resid
        qui rdrobust lrape_mres_b running, covs(police cdum2-cdum21) ///
            vce(cluster regionc_m) p(2) kernel(uniform)
        matrix bs_c6[`b', 1] = e(tau_cl)
    }
    clear
    matrix colnames bs_c6 = tau
    svmat bs_c6, names(col)
    qui summ tau
    scalar t2_se_c6_boot = r(sd)
    scalar t2_p_c6_boot  = 2*(1 - normal(abs(t2_b_c6/t2_se_c6_boot)))

    di as txt "Bootstrap SE col 3: " t2_se_c3_boot "  bootstrap p: " t2_p_c3_boot
    di as txt "Bootstrap SE col 6: " t2_se_c6_boot "  bootstrap p: " t2_p_c6_boot
}
*Block: write t2.tex ----
{
    capture file close t2
    file open t2 using ".\output\t2.tex", write replace

    foreach c in c1 c3 c4 c6 {
        local b   = t2_b_`c'
        local se  = t2_se_`c'
        local p   = t2_p_`c'
        local stars = cond(`p'<.01,"***",cond(`p'<.05,"**",cond(`p'<.10,"*","")))
        local b_`c'_s  = string(`b',  "%9.3f") + "`stars'"
        local se_`c'_s = "(" + string(`se',"%9.3f") + ")"
        local p_`c'_s  = "\textit{" + string(`p',"%9.3f") + "}"
        local bw_`c'_s = string(t2_bw_`c',"%9.1f")
        local nl_`c'_s = string(t2_nl_`c',"%9.0f")
        local nr_`c'_s = string(t2_nr_`c',"%9.0f")
    }

    file write t2 "\begin{tabular}{lcccccc}" _n
    file write t2 " & (1) & (2) & (3) & (4) & (5) & (6) \\" _n
    file write t2 "Reform = 1 & `b_c1_s' & N/A & `b_c3_s' & `b_c4_s' & N/A & `b_c6_s' \\" _n
    file write t2 " & `se_c1_s' & N/A & `se_c3_s' & `se_c4_s' & N/A & `se_c6_s' \\" _n
    file write t2 "p-value & `p_c1_s' & \textit{N/A} & `p_c3_s' & `p_c4_s' & \textit{N/A} & `p_c6_s' \\" _n
    file write t2 "Bandwidth & `bw_c1_s' & N/A & `bw_c3_s' & `bw_c4_s' & N/A & `bw_c6_s' \\" _n
    file write t2 "Obs., left & `nl_c1_s' & N/A & `nl_c3_s' & `nl_c4_s' & N/A & `nl_c6_s' \\" _n
    file write t2 "Obs., right & `nr_c1_s' & N/A & `nr_c3_s' & `nr_c4_s' & N/A & `nr_c6_s' \\" _n
    local p_c3_boot_s = "\textit{" + string(t2_p_c3_boot,"%9.3f") + "}"
    local p_c6_boot_s = "\textit{" + string(t2_p_c6_boot,"%9.3f") + "}"
    file write t2 "Wild cluster bootstrap p-value &  & N/A & `p_c3_boot_s' &  & N/A & `p_c6_boot_s' \\" _n
    file write t2 "Year FE & x & x &  & x & x &  \\" _n
    file write t2 "Region FE & x & x & x & x & x & x \\" _n
    file write t2 "Month FE & x & x &  & x & x &  \\" _n
    file write t2 "Month resid. &  &  & x &  &  & x \\" _n
    file write t2 "\end{tabular}" _n
    file close t2
}
*Table 3: RD estimates through 2005 ----
{
    use ".\data\Sweden_tables.dta", clear
    keep if year<=2005

    gen running         = ym - 468
    gen Treat           = (running >= 0)
    gen Treat_running   = Treat * running
    gen running_2       = running^2
    gen Treat_running_2 = Treat * running_2
    gen pop             = pop_f_ + pop_m_

    * season dummies (3-month blocks) used in cols 1, 4
    gen season = .
    replace season = 1 if month < 4
    replace season = 2 if month > 3 & month < 7
    replace season = 3 if month > 6 & month < 10
    replace season = 4 if month > 9
    qui tab season, gen(season_)

    * centred/calendar seasons used in cols 2, 5
    gen season_c = .
    replace season_c = 1 if month<3 | month==12
    replace season_c = 2 if month>2 & month<6
    replace season_c = 3 if month>5 & month<9
    replace season_c = 4 if month>8 & month<12
    qui tab season_c, gen(season_c_)

    * month residualisation used in cols 3 and 6
    qui reg lrape i.month
    predict lrape_mres, resid

    *col 1 -- Ciacci spec (3-month seasons, vce(nn 3))
    qui rdrobust lrape ym, c(468) p(2) kernel(uni) vce(nn 3) ///
        covs(police season_1 season_2 season_3)
    scalar t3_b_c1  = e(tau_cl)
    scalar t3_se_c1 = e(se_tau_cl)
    scalar t3_p_c1  = e(pv_rb)
    scalar t3_bw_c1 = e(h_l)
    scalar t3_nl_c1 = e(N_h_l)
    scalar t3_nr_c1 = e(N_h_r)

    *col 2 -- centred/calendar seasons
    qui rdrobust lrape running, covs(police season_c_1 season_c_2 season_c_3) ///
        p(2) kernel(uniform)
    scalar t3_b_c2  = e(tau_cl)
    scalar t3_se_c2 = e(se_tau_cl)
    scalar t3_p_c2  = e(pv_rb)
    scalar t3_bw_c2 = e(h_l)
    scalar t3_nl_c2 = e(N_h_l)
    scalar t3_nr_c2 = e(N_h_r)

    *col 3 -- month residualisation
    qui rdrobust lrape_mres running, covs(police) p(2) kernel(uniform)
    scalar t3_b_c3  = e(tau_cl)
    scalar t3_se_c3 = e(se_tau_cl)
    scalar t3_p_c3  = e(pv_rb)
    scalar t3_bw_c3 = e(h_l)
    scalar t3_nl_c3 = e(N_h_l)
    scalar t3_nr_c3 = e(N_h_r)

    *col 4 -- Ciacci spec with population weights
    qui rdrobust lrape ym, c(468) p(2) weights(pop) kernel(uni) vce(nn 3) ///
        covs(police season_1 season_2 season_3)
    scalar t3_b_c4  = e(tau_cl)
    scalar t3_se_c4 = e(se_tau_cl)
    scalar t3_p_c4  = e(pv_rb)
    scalar t3_bw_c4 = e(h_l)
    scalar t3_nl_c4 = e(N_h_l)
    scalar t3_nr_c4 = e(N_h_r)

    *col 5 -- centred/calendar seasons with population weights
    qui rdrobust lrape running, covs(police season_c_1 season_c_2 season_c_3) ///
        weights(pop) p(2) kernel(uniform)
    scalar t3_b_c5  = e(tau_cl)
    scalar t3_se_c5 = e(se_tau_cl)
    scalar t3_p_c5  = e(pv_rb)
    scalar t3_bw_c5 = e(h_l)
    scalar t3_nl_c5 = e(N_h_l)
    scalar t3_nr_c5 = e(N_h_r)

    *col 6 -- month residualisation with population weights
    qui rdrobust lrape_mres running, covs(police) weights(pop) ///
        p(2) kernel(uniform)
    scalar t3_b_c6  = e(tau_cl)
    scalar t3_se_c6 = e(se_tau_cl)
    scalar t3_p_c6  = e(pv_rb)
    scalar t3_bw_c6 = e(h_l)
    scalar t3_nl_c6 = e(N_h_l)
    scalar t3_nr_c6 = e(N_h_r)

    *--- Write t3.tex ---
    capture file close t3
    file open t3 using ".\output\t3.tex", write replace

    foreach c in c1 c2 c3 c4 c5 c6 {
        local b   = t3_b_`c'
        local se  = t3_se_`c'
        local p   = t3_p_`c'
        local stars = cond(`p'<.01,"***",cond(`p'<.05,"**",cond(`p'<.10,"*","")))
        local b_`c'_s  = string(`b',"%9.3f") + "`stars'"
        local se_`c'_s = "(" + string(`se',"%9.3f") + ")"
        local p_`c'_s  = "\textit{" + string(`p',"%9.3f") + "}"
        local bw_`c'_s = string(t3_bw_`c',"%9.1f")
        local nl_`c'_s = string(t3_nl_`c',"%9.0f")
        local nr_`c'_s = string(t3_nr_`c',"%9.0f")
    }

    file write t3 "\begin{tabular}{lcccccc}" _n
    file write t3 " & (1) & (2) & (3) & (4) & (5) & (6) \\" _n
    file write t3 "Reform = 1 & `b_c1_s' & `b_c2_s' & `b_c3_s' & `b_c4_s' & `b_c5_s' & `b_c6_s' \\" _n
    file write t3 " & `se_c1_s' & `se_c2_s' & `se_c3_s' & `se_c4_s' & `se_c5_s' & `se_c6_s' \\" _n
    file write t3 "p-value & `p_c1_s' & `p_c2_s' & `p_c3_s' & `p_c4_s' & `p_c5_s' & `p_c6_s' \\" _n
    file write t3 " & & & & & & \\" _n
    file write t3 "Bandwidth & `bw_c1_s' & `bw_c2_s' & `bw_c3_s' & `bw_c4_s' & `bw_c5_s' & `bw_c6_s' \\" _n
    file write t3 "Obs., left & `nl_c1_s' & `nl_c2_s' & `nl_c3_s' & `nl_c4_s' & `nl_c5_s' & `nl_c6_s' \\" _n
    file write t3 "Obs., right & `nr_c1_s' & `nr_c2_s' & `nr_c3_s' & `nr_c4_s' & `nr_c5_s' & `nr_c6_s' \\" _n
    file write t3 "\end{tabular}" _n
    file close t3
}
*Table 4: sex-purchase prosecutions by year ----
{
    use ".\data\Sweden_tables.dta", clear
    collapse (sum) sex_purchase_abs_, by(year)

    capture file close t4
    file open t4 using ".\output\t4.tex", write replace

    local ncols = _N
    local colspec "l"
    forvalues i = 1/`ncols' {
        local colspec "`colspec'c"
    }

    local header "Year"
    local row    "Sex purchase prosecutions"
    forvalues i = 1/`ncols' {
        local y = year[`i']
        local v = sex_purchase_abs_[`i']
        local header "`header' & `y'"
        local row    "`row' & `v'"
    }

    file write t4 "\begin{tabular}{`colspec'}" _n
    file write t4 "`header' \\" _n
    file write t4 "`row' \\" _n
    file write t4 "\end{tabular}" _n
    file close t4
}


*==============================================================================
* APPENDIX EXHIBITS
*==============================================================================


*Figure B1: RDiT variation, whole sample (wide window, -24 to +28 months) ----
*  Maps to manuscript Appendix Figure B1 ("RDiT variation, whole sample"),
*  transplanted from Ciacci r&r do-file lines 257-284.
{
    use `prepped', clear
    preserve
    keep if ym>=468-24 & ym<=468+29
    qui reg Treat Treat_running Treat_running_2 running running_2 i.regionc ///
        ydum2 ydum3 ydum4 mdum2 mdum3 mdum4 mdum5 mdum6 mdum7 mdum8 mdum9 mdum10 mdum11, cl(regionc_m)
    cap drop residx
    predict residx, residuals
    qui reg lrape Treat_running Treat_running_2 running running_2 i.regionc ///
        ydum2 ydum3 ydum4 mdum2 mdum3 mdum4 mdum5 mdum6 mdum7 mdum8 mdum9 mdum10 mdum11, cl(regionc_m)
    cap drop residy
    predict residy, residuals

    collapse (mean) residx residy running rawx=Treat rawy=lrape, by(ym)
    la var residx "Treatment, net of FEs and controls"
    la var residy "Log(rapes+1), net of FEs and controls"
    la var rawx "Treatment"
    la var rawy "Log(rapes+1)"

    line residy rawy running, xline(-.5) legend(pos(6) row(1)) lwidth(medthick) ///
        title(A: Outcome) xtitle("Months to treatment") xlabel(-24(6)28,labsize(medsmall)) ///
        ylabel(,labsize(medsmall)) legend(size(small)) ///
        xline(-13 -1 11 23, lcolor(green%24) lp(solid) lw(2.25)) ///
        xline(-24 -12 0 12 24, lcolor(red%24) lp(solid) lw(2.25)) ///
        scheme(s1mono)
    graph export ".\output\FB1A.pdf", replace

    line residx rawx running, xline(-.5) legend(pos(6) row(1)) lwidth(medthick) ///
        title(B: Treatment) xtitle("Months to treatment") xlabel(-24(6)28,labsize(medsmall)) ///
        ylabel(,labsize(medsmall)) legend(size(small)) ///
        xline(-13 -1 11 23, lcolor(green%24) lp(solid) lw(2.25)) ///
        xline(-24 -12 0 12 24, lcolor(red%24) lp(solid) lw(2.25)) ///
        scheme(s1mono)
    graph export ".\output\FB1B.pdf", replace
    restore
}


*Figure B2: RDiT variation, restricted sample (narrow window +/- 9 months) ----
*  Maps to manuscript Appendix Figure B2 ("RDiT variation, restricted sample");
*  the narrow window itself limits the panel to ym 459-477 (Apr 1997 - Oct 1999).
*  FLAG: if B2 is meant to be year<2001 rather than narrow-window,
*  add `keep if year<2001` at the top of this block.
{
    use `prepped', clear
    preserve
    keep if ym>=468-9 & ym<=468+9
    qui reg Treat Treat_running Treat_running_2 running running_2 i.regionc ///
        mdum2 mdum3 mdum4 mdum5 mdum6 mdum7 mdum8 mdum9 mdum10, cl(regionc_m)
    cap drop residx
    predict residx, residuals
    qui reg lrape Treat_running Treat_running_2 running running_2 i.regionc ///
        mdum2 mdum3 mdum4 mdum5 mdum6 mdum7 mdum8 mdum9 mdum10, cl(regionc_m)
    cap drop residy
    predict residy, residuals

    collapse (mean) residx residy rawx=Treat rawy=lrape running, by(ym)
    la var residx "Treatment, net of FEs and controls"
    la var residy "Log(rapes+1), net of FEs and controls"
    la var rawx "Treatment"
    la var rawy "Log(rapes+1)"

    line residy rawy running, xline(-.5) xlabel(-9(3)9) title(A: Outcome) ///
        legend(pos(6) row(1)) lwidth(medthick) xtitle("Months to treatment") ///
        ylabel(,labsize(medsmall)) legend(size(small)) ///
        xline(-1, lcolor(green%24) lp(solid) lw(6.5)) ///
        xline(-2 0, lcolor(red%24) lp(solid) lw(6.5)) scheme(s1mono)
    graph export ".\output\FB2A.pdf", replace

    line residx rawx running, xline(-.5) xlabel(-9(3)9) title(B: Treatment) ///
        legend(pos(6) row(1)) lwidth(medthick) xtitle("Months to treatment") ///
        ylabel(,labsize(medsmall)) legend(size(small)) ///
        xline(-1, lcolor(green%24) lp(solid) lw(6.5)) ///
        xline(-2 0, lcolor(red%24) lp(solid) lw(6.5)) scheme(s1mono)
    graph export ".\output\FB2B.pdf", replace
    restore
}


*Figure B4: RDiT variation for Ciacci (2025), Table 2 col 1 ----
*  Maps to manuscript Appendix Figure B4. Uses Ciacci's 3-month season dummies
*  and narrow +/-7-month window, matching the Ciacci (2025) Table 2 col 1 spec.
{
    use `prepped', clear
    preserve
    keep if ym>=468-7 & ym<=468+7
    gen season = .
    replace season = 1 if month < 4
    replace season = 2 if month > 3 & month < 7
    replace season = 3 if month > 6 & month < 10
    replace season = 4 if month > 9
    qui tab season, gen(season_)

    qui reg Treat Treat_running Treat_running_2 running running_2 police ///
        season_1 season_2 season_3, cl(regionc_m)
    cap drop residx
    predict residx, residuals
    qui reg lrape Treat_running Treat_running_2 running running_2 police ///
        season_1 season_2 season_3, cl(regionc_m)
    cap drop residy
    predict residy, residuals

    collapse (mean) residx residy rawx=Treat rawy=lrape running, by(ym)
    la var residx "Treatment, net of FEs and controls"
    la var residy "Log(rapes+1), net of FEs and controls"
    la var rawx "Treatment"
    la var rawy "Log(rapes+1)"

    line residy rawy running, xline(-.5) xlabel(-6(2)6) title(A: Outcome) ///
        legend(pos(6) row(1)) lwidth(medthick) xtitle("Months to treatment") ///
        ylabel(,labsize(medsmall)) legend(size(small)) ///
        xline(-4 -1 2 5, lcolor(green%24) lp(solid) lw(8.5)) ///
        xline(-6 -3 0 3 6, lcolor(red%24) lp(solid) lw(8.5)) scheme(s1mono)
    graph export ".\output\FB4A.pdf", replace

    line residx rawx running, xline(-.5) xlabel(-6(3)6) title(B: Treatment) ///
        legend(pos(6) row(1)) lwidth(medthick) xtitle("Months to treatment") ///
        ylabel(,labsize(medsmall)) legend(size(small)) ///
        xline(-4 -1 2 5, lcolor(green%24) lp(solid) lw(8.5)) ///
        xline(-6 -3 0 3 6, lcolor(red%24) lp(solid) lw(8.5)) scheme(s1mono)
    graph export ".\output\FB4B.pdf", replace
    restore
}


*Figure B5: four RDiT rdplot panels (uniform vs triangular, raw vs month-resid) ----
* STATUS: READY -- exports ./output/FB5.pdf (combined 4-panel).
*  Source: "new rdd code 251029.do" lines 76-87. Uses the Ciacci 2025 sample
*  (year<=2005) with both `season_` (3-month blocks) and `season_c_`
*  (centered/calendar seasons) dummies.
{
    use ".\data\Sweden_tables.dta", clear
    keep if year<=2005
    gen running = ym - 468
    *3-month block seasons
    gen season = .
    replace season = 1 if month < 4
    replace season = 2 if month > 3 & month < 7
    replace season = 3 if month > 6 & month < 10
    replace season = 4 if month > 9
    qui tab season, gen(season_)
    *centered/calendar seasons
    gen season_c = .
    replace season_c = 1 if month<3 | month==12
    replace season_c = 2 if month>2 & month<6
    replace season_c = 3 if month>5 & month<9
    replace season_c = 4 if month>8 & month<12
    qui tab season_c, gen(season_c_)

    qui reg lrape i.month
    predict lrape_mres, resid

    tempfile fb5_1 fb5_2 fb5_3 fb5_4

    * Panel 1: lrape, uniform kernel
    rdplot lrape running if running>=-9.7 & running<=9.7,  h(9.7 9.7) ///
        covs(police season_c_1 season_c_2 season_c_3) p(2) kernel(uniform) ///
        graph_options(legend(off) title("Outcome: log(rape+1), uniform"))
    graph save "`fb5_1'", replace

    * Panel 2: lrape_mres, uniform kernel
    rdplot lrape_mres running if running>=-8.7 & running<=8.7, h(8.7 8.7) ///
        covs(police season_c_1 season_c_2 season_c_3) p(2) kernel(uniform) ///
        graph_options(legend(off) title("Outcome: month-resid log(rape+1), uniform"))
    graph save "`fb5_2'", replace

    * Panel 3: lrape, triangular kernel (rdplot default)
    rdplot lrape running if running>=-14.4 & running<=14.4, h(14.4 14.4) ///
        covs(police season_c_1 season_c_2 season_c_3) p(2) ///
        graph_options(legend(off) title("Outcome: log(rape+1), triangular"))
    graph save "`fb5_3'", replace

    * Panel 4: lrape_mres, triangular kernel
    rdplot lrape_mres running if running>=-15.7 & running<=15.7, h(15.7 15.7) ///
        covs(police season_c_1 season_c_2 season_c_3) p(2) ///
        graph_options(legend(off) title("Outcome: month-resid log(rape+1), triangular"))
    graph save "`fb5_4'", replace

    graph combine "`fb5_1'" "`fb5_3'" "`fb5_2'" "`fb5_4'", ///
        scheme(s1mono) col(2) ysize(13) xsize(20) iscale(.6)
    cap noi graph export ".\output\FB5.pdf", replace
}


*Figure B6: RDiT treatment effects at placebo thresholds ----
*  Source: "new rdd code 251029.do" lines 91-148. Sweeps the RD cutoff
*  from -12 to +11 months around the true threshold and stores
*  treatment-effect estimates and 90% CIs.
{
    use ".\data\Sweden_tables.dta", clear
    keep if year<=2005
    gen running = ym - 468
    qui reg lrape i.month
    predict lrape_mres, resid

    tempfile fb6_uni fb6_tri

    *--- uniform kernel ---
    preserve
    gen running_alt = 0
    matrix B = J(50,4,.)
    forvalues n = 1/24 {
        scalar bandwidth = `n' - 13
        qui replace running_alt = running - 13 + `n'
        qui rdrobust lrape_mres running_alt, c(0) covs(police) kernel(uniform) p(2)
        matrix B[`n',1] = e(tau_cl)
        matrix B[`n',2] = e(tau_cl) + 1.645*e(se_tau_cl)
        matrix B[`n',3] = e(tau_cl) - 1.645*e(se_tau_cl)
        matrix B[`n',4] = bandwidth
    }
    svmat B
    rename B1 estimate
    rename B2 upper
    rename B3 lower
    rename B4 bandwidth_pl

    graph twoway (scatter estimate bandwidth_pl, mc(gs1) msymbol(square)) ///
        (rbar upper lower bandwidth_pl, fcolor(gs0) fintensity(100) lwidth(none) barwidth(.05)), ///
        legend(off) xline(0, lp(dash)) yline(0) ///
        ytitle("Estimated treatment effect") ///
        xtitle("Months between placebo threshold and actual threshold") ///
        scheme(s1mono) ylabel(, angle(horizontal) grid) xlabel(-12(3)12) ///
        plotregion(lcolor(none))
    graph save "`fb6_uni'", replace
    restore

    *--- triangular kernel (default) ---
    preserve
    gen running_alt = 0
    matrix B = J(50,4,.)
    forvalues n = 1/24 {
        scalar bandwidth = `n' - 13
        qui replace running_alt = running - 13 + `n'
        qui rdrobust lrape_mres running_alt, c(0) covs(police) p(2)
        matrix B[`n',1] = e(tau_cl)
        matrix B[`n',2] = e(tau_cl) + 1.645*e(se_tau_cl)
        matrix B[`n',3] = e(tau_cl) - 1.645*e(se_tau_cl)
        matrix B[`n',4] = bandwidth
    }
    svmat B
    rename B1 estimate
    rename B2 upper
    rename B3 lower
    rename B4 bandwidth_pl

    graph twoway (scatter estimate bandwidth_pl, mc(gs1) msymbol(square)) ///
        (rbar upper lower bandwidth_pl, fcolor(gs0) fintensity(100) lwidth(none) barwidth(.05)), ///
        legend(off) xline(0, lp(dash)) yline(0) ///
        ytitle("Estimated treatment effect") ///
        xtitle("Months between placebo threshold and actual threshold") ///
        scheme(s1mono) ylabel(, angle(horizontal) grid) xlabel(-12(3)12) ///
        plotregion(lcolor(none))
    graph save "`fb6_tri'", replace
    restore

    graph combine "`fb6_uni'" "`fb6_tri'", ///
        scheme(s1mono) col(2) ysize(7) xsize(20) iscale(1.2) ycommon
    cap noi graph export ".\output\FB6.pdf", replace
}


*Tables B1 & B2: RD estimates on lrape and lrape2 with optimal bandwidth ----
*   Matches the structure of Appendix RR2 Tables B1 and B2:
*     - 2 columns: (1) Log(x+1),  (2) IHS(x)
*     - 2 panels:  A = full sample, B = year<2001
*     - Two rows per panel: first-order and second-order polynomial in the
*       running variable, each reporting the Treatment coefficient.
*   Table B1 uses the non-differenced outcomes with controls
*     `police i.regionc i.year i.month'  (clustering on regionc_m).
*   Table B2 uses the first-differenced outcomes with controls
*     `dpolice i.year'                   (clustering on regionc_m).
{
    *--- capture scalars for every {table x sample x outcome x polynomial} cell
    foreach table in B1 B2 {
        foreach sample in A B {
            use `prepped', clear
            if "`sample'"=="B" keep if year<2001

            * choose outcomes + controls by table
            if "`table'"=="B1" {
                local outcomes  lrape ihs_rape
                local controls "police i.regionc i.year i.month"
            }
            else {
                local outcomes  dlrape dihs_rape
                local controls "dpolice i.year"
            }

            local cols 1 2
            foreach col of local cols {
                local out : word `col' of `outcomes'
                qui rdbwselect `out' ym, c(468)
                local h_l = 468 - e(h_mserd)
                local h_u = 468 + e(h_mserd)

                *first-order polynomial
                qui reg `out' Treat Treat_running running `controls' ///
                    if ym>`h_l' & ym<`h_u', cl(regionc_m)
                scalar b_`table'_`sample'_`col'_1  = _b[Treat]
                scalar se_`table'_`sample'_`col'_1 = _se[Treat]
                scalar p_`table'_`sample'_`col'_1  = 2*ttail(e(df_r), abs(_b[Treat]/_se[Treat]))
                scalar N_`table'_`sample'_`col'_1  = e(N)

                *second-order polynomial
                qui reg `out' Treat Treat_running Treat_running_2 running running_2 `controls' ///
                    if ym>`h_l' & ym<`h_u', cl(regionc_m)
                scalar b_`table'_`sample'_`col'_2  = _b[Treat]
                scalar se_`table'_`sample'_`col'_2 = _se[Treat]
                scalar p_`table'_`sample'_`col'_2  = 2*ttail(e(df_r), abs(_b[Treat]/_se[Treat]))
                scalar N_`table'_`sample'_`col'_2  = e(N)
            }
        }
    }

    *--- write tB1.tex and tB2.tex -----------------------------------------
    foreach table in B1 B2 {
        local outpath "./output/t`table'.tex"
        capture file close t_out
        file open t_out using "`outpath'", write replace

        file write t_out "\begin{tabular}{lcc}" _n
        file write t_out " & (1) & (2) \\" _n
        file write t_out " & Log(x+1) & IHS(x) \\" _n

        foreach sample in A B {
            local panel_label = cond("`sample'"=="A","Panel A: full sample","Panel B: restricted sample (year<2001)")
            file write t_out "\multicolumn{3}{l}{\textit{`panel_label'}} \\" _n

            forvalues poly = 1/2 {
                local poly_label = cond(`poly'==1,"First order polynomial: Treatment","Second order polynomial: Treatment")
                local b1  = b_`table'_`sample'_1_`poly'
                local b2  = b_`table'_`sample'_2_`poly'
                local se1 = se_`table'_`sample'_1_`poly'
                local se2 = se_`table'_`sample'_2_`poly'
                local p1  = p_`table'_`sample'_1_`poly'
                local p2  = p_`table'_`sample'_2_`poly'
                local st1 = cond(`p1'<.01,"***",cond(`p1'<.05,"**",cond(`p1'<.10,"*","")))
                local st2 = cond(`p2'<.01,"***",cond(`p2'<.05,"**",cond(`p2'<.10,"*","")))
                local b1_s  = string(`b1', "%9.3f") + "`st1'"
                local b2_s  = string(`b2', "%9.3f") + "`st2'"
                local se1_s = "(" + string(`se1',"%9.3f") + ")"
                local se2_s = "(" + string(`se2',"%9.3f") + ")"
                local p1_s  = "\textit{" + string(`p1',"%9.3f") + "}"
                local p2_s  = "\textit{" + string(`p2',"%9.3f") + "}"

                file write t_out "`poly_label' & `b1_s' & `b2_s' \\" _n
                file write t_out " & `se1_s' & `se2_s' \\" _n
                file write t_out "p-value & `p1_s' & `p2_s' \\" _n
            }
            local N1 = N_`table'_`sample'_1_1
            local N2 = N_`table'_`sample'_2_1
            file write t_out "Observations & " (string(`N1',"%9.0f")) " & " (string(`N2',"%9.0f")) " \\" _n
            if "`sample'"=="A" file write t_out " & & \\" _n
        }

        file write t_out "\end{tabular}" _n
        file close t_out
    }
}


*Table C1: Fuzzy DiD ----
*  Runs in ./code/replication_fuzzydid.do (split out because fuzzydid
*  with breps=99 takes several minutes and is rarely re-run).


*Table C2: Event study with out-of-window dummy ----
*  Matches Appendix RR2 Table C2: two columns.
*    (1) main     baseline event study
*    (2) nw_dummy adds a dummy for observations outside the event window
{
    cap which reghdfe
    local _ok1 = (_rc==0)
    cap which outreg2
    local _ok2 = (_rc==0)
    if !`_ok1' | !`_ok2' {
        di as error "Table C2 skipped: need `reghdfe' and `outreg2'."
    }
    else {
        use ".\data\Sweden_tables.dta", clear
        cap drop sum_sex_pur dummy_cum_sex_pr rt_sp rt_sex_pr

        bys regionc (ym): gen sum_sex_pur = sum(sex_purchase_abs_)
        gen dummy_cum_sex_pr = (sex_purchase_abs_>0)
        gen time_sp = .
        replace time_sp = ym if dummy_cum_sex_pr==1
        replace time_sp = time_sp[_n+2] if dummy_cum_sex_pr[_n+2]==1 & time_sp==. & regionc==regionc[_n+2]
        replace time_sp = time_sp[_n+1] if dummy_cum_sex_pr[_n+1]==1 & time_sp==. & regionc==regionc[_n+1]
        replace time_sp = time_sp[_n-1] if dummy_cum_sex_pr[_n-1]==1 & time_sp==. & regionc==regionc[_n-1]
        replace time_sp = time_sp[_n-2] if dummy_cum_sex_pr[_n-2]==1 & time_sp==. & regionc==regionc[_n-2]
        bys regionc: gen running_time_sp = ym - time_sp
        gen rt_sp = running_time_sp
        recode rt_sp (.=-1)
        char rt_sp[omit] -1
        xi i.rt_sp, pref(R_sp)

        cap la var R_sprt_sp_1 "2 months pre-fine"
        cap la var R_sprt_sp_3 "Fine month"
        cap la var R_sprt_sp_4 "1 month post-fine"
        cap la var R_sprt_sp_5 "2 months post-fine"

        * Out-of-window flag (col 2)
        gen non_window = running==.

        cap erase ".\output\tC2.tex"
        * col 1 - main
        reghdfe lrape R_sp* police, a(regionc year month) cl(regionc year month) version(3)
        outreg2 using ".\output\tC2.tex", tex(frag) dec(3) replace ctitle(main) ///
            keep(R_sprt_sp_1 R_sprt_sp_3 R_sprt_sp_4 R_sprt_sp_5) label stats(coef se pval) nocons
        * col 2 - out-of-window dummy added
        reghdfe lrape R_sp* police non_window, a(regionc year month) cl(regionc year month) version(3)
        outreg2 using ".\output\tC2.tex", tex(frag) dec(3) append ctitle(nw_dummy) ///
            keep(R_sprt_sp_1 R_sprt_sp_3 R_sprt_sp_4 R_sprt_sp_5) label stats(coef se pval) nocons
    }
}


*Table C3: Event study, alternative clustering ----
*  Requires `reghdfe', `outreg2' (both on SSC).
{
    cap which reghdfe
    local _ok1 = (_rc==0)
    cap which outreg2
    local _ok2 = (_rc==0)
    if !`_ok1' | !`_ok2' {
        di as error "Table C2 skipped: need `reghdfe' and `outreg2'. Run: ssc install reghdfe; ssc install outreg2"
    }
    else {
        use ".\data\Sweden_tables.dta", clear
        cap drop sum_sex_pur dummy_cum_sex_pr running_time_sex_pr rt_sp rt_sex_pr
        bys regionc (ym): gen sum_sex_pur = sum(sex_purchase_abs_)
        gen dummy_cum_sex_pr = (sex_purchase_abs_>0)

        gen time_sp = .
        replace time_sp = ym if dummy_cum_sex_pr==1
        replace time_sp = time_sp[_n+2] if dummy_cum_sex_pr[_n+2]==1 & time_sp==. & regionc==regionc[_n+2]
        replace time_sp = time_sp[_n+1] if dummy_cum_sex_pr[_n+1]==1 & time_sp==. & regionc==regionc[_n+1]
        replace time_sp = time_sp[_n-1] if dummy_cum_sex_pr[_n-1]==1 & time_sp==. & regionc==regionc[_n-1]
        replace time_sp = time_sp[_n-2] if dummy_cum_sex_pr[_n-2]==1 & time_sp==. & regionc==regionc[_n-2]
        bys regionc: gen running_time_sp = ym - time_sp
        gen rt_sp = running_time_sp
        recode rt_sp (.=-1)
        char rt_sp[omit] -1
        xi i.rt_sp, pref(R_sp)
        cap la var R_sprt_sp_1 "2 months pre-fine"
        cap la var R_sprt_sp_3 "Fine month"
        cap la var R_sprt_sp_4 "1 month post-fine"
        cap la var R_sprt_sp_5 "2 months post-fine"

        egen r_y   = group(regionc year)
        egen r_m   = group(regionc month)
        egen r_y_m = group(regionc year month)

        cap erase ".\output\tC3.tex"
        reghdfe lrape R_sp* police, a(regionc year month) cl(regionc year month) version(3)
        outreg2 using ".\output\tC3.tex", tex(frag) dec(3) replace ctitle(main)              keep(R_sp*) label stats(coef se pval) nocons
        reghdfe lrape R_sp* police, a(regionc year month) cl(region)
        outreg2 using ".\output\tC3.tex", tex(frag) dec(3) append  ctitle(region)            keep(R_sp*) label stats(coef se pval) nocons
        reghdfe lrape R_sp* police, a(regionc year month) cl(r_y)
        outreg2 using ".\output\tC3.tex", tex(frag) dec(3) append  ctitle(region*year)       keep(R_sp*) label stats(coef se pval) nocons
        reghdfe lrape R_sp* police, a(regionc year month) cl(r_m)
        outreg2 using ".\output\tC3.tex", tex(frag) dec(3) append  ctitle(region*month)      keep(R_sp*) label stats(coef se pval) nocons
        reghdfe lrape R_sp* police, a(regionc year month) cl(r_y_m)
        outreg2 using ".\output\tC3.tex", tex(frag) dec(3) append  ctitle(region*year*month) keep(R_sp*) label stats(coef se pval) nocons noaster
        reghdfe lrape R_sp* police, a(regionc year month)
        outreg2 using ".\output\tC3.tex", tex(frag) dec(3) append  ctitle(no clustering)     keep(R_sp*) label stats(coef se pval) nocons

        *Keep the data for Table C3 below
        tempfile event_data
        save `event_data'
        global event_data "`event_data'"
    }
}

*Table C5: IV specifications ----
*  Manuscript Table C4 is the airports list (hardcoded, not produced here).
{
    cap which outreg2
    if _rc {
        di as error "Table C5 skipped: `outreg2' not installed."
    }
    else {
        use ".\data\Sweden_tables3.dta", clear
        xtset regionc ym
        foreach var in new_land_inter3 new_land_eintra4 {
            gen `var'_ues = `var'
            gen temp       = `var' if region==12
            bysort month year: egen `var'_stock = max(temp)
            replace `var'_ues = `var'_stock if region==14
            drop temp `var'_stock
            gen `var'_ns = `var'
            replace `var'_ns = 0 if region==12
        }

        cap la var sex_purchase_abs_ "Sex purchase prosecutions"
        cap erase ".\output\tC5.tex"

        *helper: first-stage F via a separate OLS + test (more robust than
        *`estat firststage', which can return empty r(mineig) here).
        *col 1
        xi: regress sex_purchase_abs_ new_land_inter3 new_land_eintra4 police ///
            i.regionc*year i.year i.month, cluster(regionc)
        test new_land_inter3 new_land_eintra4
        scalar F1 = r(F)
        xi: ivregress 2sls lrape (sex_purchase_abs_ = new_land_inter3 new_land_eintra4) ///
            police i.regionc*year i.year i.month, cluster(regionc)
        outreg2 using ".\output\tC5.tex", tex(frag) dec(3) replace ctitle(main) ///
            keep(sex_purchase_abs_) label addstat("First-stage F", F1) nocons

        *col 2
        xi: regress sex_purchase_abs_ new_land_inter3_ues new_land_eintra4_ues police ///
            i.regionc i.regionc*year i.year i.month, cluster(regionc)
        test new_land_inter3_ues new_land_eintra4_ues
        scalar F2 = r(F)
        xi: ivregress 2sls lrape (sex_purchase_abs_ = new_land_inter3_ues new_land_eintra4_ues) ///
            police i.regionc i.regionc*year i.year i.month, cluster(regionc)
        *freeze col 2's estimation sample so col 3 uses the same rows
        gen byte c5_samp2 = e(sample)
        outreg2 using ".\output\tC5.tex", tex(frag) dec(3) append  ctitle(upp_inc) ///
            keep(sex_purchase_abs_) label addstat("First-stage F", F2) nocons

        *col 3 -- Ciacci's "Stockholm excluded" spec: instrument values for
        *region 12 set to 0 via the *_ns variables. Sample restricted to
        *col 2's e(sample) so N matches across columns.
        xi: regress sex_purchase_abs_ new_land_inter3_ns new_land_eintra4_ns police ///
            i.regionc i.regionc*year i.year i.month if c5_samp2==1, cluster(regionc)
        test new_land_inter3_ns new_land_eintra4_ns
        scalar F3 = r(F)
        xi: ivregress 2sls lrape (sex_purchase_abs_ = new_land_inter3_ns new_land_eintra4_ns) ///
            police i.regionc i.regionc*year i.year i.month if c5_samp2==1, cluster(regionc)
        outreg2 using ".\output\tC5.tex", tex(frag) dec(3) append  ctitle(sto_excl) ///
            keep(sex_purchase_abs_) label addstat("First-stage F", F3) nocons
    }
}

*Cleanup: outreg2 always writes a .txt companion next to its .tex output.
*There is no option to suppress it, so erase them after the fact.
{
    foreach f in tC2 tC3 tC5 {
        cap erase "./output/`f'.txt"
    }
}
