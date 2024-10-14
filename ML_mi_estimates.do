# MI estimate trick - ML probability output
*************************
** Internal validation **
*************************
** Set directory **
cd "filepath\development"

forval i = 1/50 {
    import delimited "filepath\XGBClassifier_cluster_mi`i'_Proba_Test.csv", clear
	gen patid= _n	
	gen m=`i'	
	order id, before(test)
	save "probabilities\XGBClassifier_cluster_mi`i'_Proba_Test", replace
}


use "probabilities\XGBClassifier_cluster_mi1_Proba_Test.dta", clear
forval i = 2/50 {
	append using  "probabilities\XGBClassifier_cluster_mi`i'_Proba_Test"
	save ndipass_period1_miset, replace
}


* First imputation - say it is m=0 (complete case) 
keep if m==1 
replace m=0 

* Generate mi style missing indicator = predictions are missing 
* later use Stata's mi machinery to handle the 50 predictions for each id
gen i_proba = 1 
replace proba = . if i_proba==1 

append using ndipass_period1_miset

save ndipass_period1_miset_ready, replace
* Format as multiply imputed dataset - m=1 to m=50 
use ndipass_period1_miset_ready, clear
mi import flong, m(m) id(patid) imputed(proba)


gen cv_xb=log(proba / (1 - proba))		// generate linear predictor

************************************************
* Internal validation metrics *
************************************************
* Calibration slope overall * 
* Regress linear predictor * 
mi estimate, dots: logistic test cv_xb
* Calibration-in-the-large * 
mi estimate, dots: logistic test, offset(cv_xb)

* Calibration plot *
pmcalplot proba test if _mi_m==1, ci nospike nostat  ///
xtitle("Predicted event probability", size(medsmall) margin(medsmall)) ytitle("Observed probability", size(medsmall) margin(medsmall)) graphregion(color(white)) ///
title ("{bf:XGBoost (development data)}", color(black) position(11) yoffset(1) size(medsmall)) legend(off) ///
name(NDIdev_pmcalplot_xgb, replace) 
cd "filepath\pmcalplot"
graph save "NDIdev_pmcalplot_xgb", replace
graph export "NDIdev_pmcalplot_xgb", as(pdf) replace

cd "filepath\development"

* Check calibration plot for consistency across imputations *
/*
forval i = 1/50 {  
    quietly pmcalplot int_proba ndipass if _mi_m == `i', ci nospike ///
    xtitle("Predicted event probability", size(small)) ///
    ytitle("Observed probability", size(small)) graphregion(color(white)) ///
    title("Calibration plot - Logistic regression model (Imputation `i')", size(small))
    
    // Add a brief pause between plots
    sleep 1
}
*/

* Discrimination *
* Define eroctab to use all imputations *
mi xeq 1: roctab test cv_xb
return list
cap program drop eroctab
program eroctab, eclass
        version 12.0

        /* Step 1: perform ROC analysis */
        args refvar classvar
        roctab `refvar' `classvar'

        /* Step 2: save estimate and its variance in temporary matrices*/
        tempname b V
        mat `b' = r(area)
        mat `V' = r(se)^2
        local N = r(N)

        /* Step 3: make column names and row names consistent*/
        mat colnames `b' = AUC
        mat colnames `V' = AUC
        mat rownames `V' = AUC

        /*Step 4: post results to e()*/
        ereturn post `b' `V', obs(`N')
        ereturn local cmd "eroctab"
        ereturn local title "ROC area"
end
mi estimate, cmdok dots: eroctab test cv_xb

* Predicted values *
* Define program to use all imputations *
cap program drop calc_ppv_npv
program calc_ppv_npv, eclass
    version 12.0
    args truevar predvar

    // Explicitly limit to non-missing predicted values
    preserve
    qui: keep if !missing(`predvar')

    // Check if there are any remaining observations
    qui: count
    if r(N) == 0 exit

    // Continue with PPV and NPV calculations
    qui: count if `truevar' == 1 & `predvar' == 1
    local TP = r(N)
    qui: count if `truevar' == 0 & `predvar' == 1
    local FP = r(N)
    qui: count if `truevar' == 1 & `predvar' == 0
    local FN = r(N)
    qui: count if `truevar' == 0 & `predvar' == 0
    local TN = r(N)

    // Calculate PPV and NPV
    local PPV = `TP' / (`TP' + `FP')
    local NPV = `TN' / (`TN' + `FN')

    // Calculate variances for PPV and NPV
    local var_PP = (`PPV' * (1 - `PPV')) / (`TP' + `FP')
    local var_NP = (`NPV' * (1 - `NPV')) / (`TN' + `FN')

    // Handle edge cases for variance calculation
    if (`TP' + `FP') == 0 local var_PP = 0
    if (`TN' + `FN') == 0 local var_NP = 0

    matrix b = (`PPV', `NPV')
    matrix V = (`var_PP', 0 \ 0, `var_NP')

    matrix colnames b = PPV NPV
    matrix colnames V = PPV NPV
    matrix rownames V = PPV NPV

    ereturn post b V, obs(`=_N')
    ereturn local cmd "calc_ppv_npv"
    ereturn local title "PPV and NPV"
    restore
end
gen byte pred = (proba >=0.5) if proba !=.

mi estimate, cmdok dots: calc_ppv_npv test pred
*************************
** Temporal validation **
*************************
** Set directory **
cd "filepath\temporal"

# Repeat steps for temporal validation predictions

*************************
** External validation **
*************************
** Set directory **
cd "filepath\external"

# Repeat steps for external validation predictions