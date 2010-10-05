#!/usr/bin/env python

# Definitions of all custom functions to be used by PYRAMDS

import numpy as np

def calc_det_limit(LM, RM, spec_array):
    """
    Calculates detection limits by first averaging counts about left/right
    markers and then using that to determine average counts per channel in
    background, mu_b. From there, critical level (el_c) and detection limit
    (el_d) are determined to a 95% confidence level.
    
    """
    
    avg_pm = 3
    N_chn = (LM - RM) + 1
    LM_avg = np.average(spec_array[(LM - avg_pm):(LM + avg_pm + 1)])
    RM_avg = np.average(spec_array[(RM - avg_pm):(RM + avg_pm + 1)])
    
    mu_b = N_chn * (LM_avg + RM_avg) * 0.5
    mu_b_sqroot = np.sqrt(mu_b)
    
    el_c = 2.325 * mu_b_sqroot
    el_d = 2.706 + 4.653 * mu_b_sqroot
    
    return el_c, el_d

def calc_det_limit_sel(chn, marker, spec_array):
    
    cent_en = en_coeff[chn][0] + en_coeff[chn][1] * marker
    width = fwhm_coeff[chn][0] + fwhm_coeff[chn][1] * marker + fwhm_coeff[chn][2] * marker * marker
    
    left_marker = int(round(cent_en - 1.75*width))
    right_marker = int(round(cent_en + 1.75*width))

    return calc_det_limit(left_marker, right_marker, spec_array)
    
def bin_count(t_start, t_stop, group):
    
    cnt_array = np.zeros( (energy_max + 1), dtype=np.int32)
    
    for chan in group:
        t_pass = [i for i in chan if (t_start <= chan[i] <= t_stop)]
        counts = len(t_pass)
        cnt_array[chan] = counts

    return cnt_array