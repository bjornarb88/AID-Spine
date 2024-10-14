## Logistic regression model ##
cd "filepath\extval\data\"
*******************************************************************************
*					         1) EXTERNAL VALIDATION RUN						  *
*******************************************************************************
use "mice\ndipass_mi50.dta", clear

** Generate FP terms ** 
mi passive: gen double age_x = age/10
mi passive: gen age_fp1 = (age_x^0.5)-2.267084573 
mi passive: gen age_fp2 = (age_x^3)-135.7707852

mi passive: gen double eqvas_x = (eqvas+1)/100
mi passive: gen eqvas_fp1 = (eqvas_x^-2)-3.568306712  

** Register as imputed variables **
mi register imputed age_fp1 age_fp2 eqvas_fp1

global covariates = "i.sex age_fp1 age_fp2 bmi i.smoker i.nativespeaker i.education i.workreq i.dispension i.paindurneck i.paindurarms nrsheadache nrsneck nrsarm i.painarms i.emshand ndi eq5d eqvas_fp1 i.anxdep prevsurg i.comorbid asa i.analgesicsfreq i.posteriorapp i.private"

save "model\ndi\ndipass_imputed_validation.dta", replace 
* Convert to long format - to save all predictions * 
mi convert flong

******************
* Run validation *
******************
* Predict outcomes for SweSpine using the developed model
mi predict double val_xb using period1_ndi_model, xb storecomplete

sum val_xb, d

* Save updated dataset with predictions: use for the model evaluation*
save "model\ndi\ndipass_imputed_validation.dta", replace 

*******************************************************************************
*	          				2) PERFORMANCE ASSESSMENT                		  *
*******************************************************************************
use  "model\ndi\ndipass_imputed_validation.dta", clear 

**************************************************
**		External validation performance			**
**************************************************
* Calibration slope *  
mi estimate, dots: logistic ndipass val_xb
* Calibration-in-the-large * 
mi estimate, dots: logistic ndipass, offset(val_xb)
* Calibration plot *
* Calculated 'predicted' probs - use linear predictors
gen val_proba = invlogit(val_xb)

pmcalplot val_proba ndipass if _mi_m==2, ci nospike nostat  ///
xtitle("Predicted event probability", size(medsmall) margin(medsmall)) ytitle("Observed probability", size(medsmall) margin(medsmall)) graphregion(color(white)) ///
title ("{bf:Logistic regression}", color(black) position(11) yoffset(1) size(medsmall)) legend(off) ///
name(NDIextval_pmcalplot, replace) 
graph save "NDIextval_pmcalplot" "model\ndi\estimates\graphs\pmcalplot\NDIextval_pmcalplot.gph", replace
cd "filepath\data\model\ndi\estimates\graphs\pmcalplot"
graph export NDIextval_pmcalplot, as(pdf) replace
cd "filepath\data\"

* Check calibration plot for consistency across imputations *
/*
forval i = 1/50 {  
    quietly pmcalplot val_proba ndipass if _mi_m == `i', ci nospike ///
    xtitle("Predicted event probability", size(small)) ///
    ytitle("Observed probability", size(small)) graphregion(color(white)) ///
    title("Calibration plot - Logistic regression model (Imputation `i')", size(small))
    
    // Add a brief pause between plots
    sleep 1
}
*/

* Discrimination *
* Define eroctab to use all imputations *
mi xeq 0: roctab ndipass recal_val_xb
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
mi estimate, cmdok dots: eroctab ndipass recal_val_xb
**
* Predicted values *
gen byte predicted_val = (recal_val_proba >=0.5) if recal_val_proba !=.
tab ndipass predicted_val

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

mi estimate, cmdok dots: calc_ppv_npv ndipass predicted_val

*******************************************************************************
*******************************************************************************