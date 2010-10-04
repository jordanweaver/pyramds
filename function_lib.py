#!/usr/bin/env python

# Definitions of all custom functions to be used by PYRAMDS

def calc_det_limit(chn, zaid, spec_array):
    """
    Calculates detection limits by first averaging counts about left/right
    markers and then using that to determine average counts per channel in
    background, mu_b. From there, critical level (el_c) and detection limit
    (el_d) are determined to a 95% confidence level.
    
    """
    
    avg_pm = 3
    LM = sig_lookup[chn][zaid][0]
    RM = sig_lookup[chn][zaid][1]
    N_chn = (LM - RM) + 1
    LM_avg = np.average(spec_array[(LM - avg_pm):(LM + avg_pm + 1)])
    RM_avg = np.average(spec_array[(RM - avg_pm):(RM + avg_pm + 1)])
    
    mu_b = N_chn * (LM_avg + RM_avg) * 0.5
    mu_b_sqroot = np.sqrt(mu_b)
    
    el_c = 2.325 * mu_b_sqroot
    el_d = 2.706 + 4.653 * mu_b_sqroot
    
    return el_c, el_d

def bin_count(t_start, t_stop, type, det):
    
    cnt_array = np.zeros( (energy_max + 1), dtype=np.int32)
    
    for bin in range(energy_max + 1):
        array_obj = '/root/spectra/%s/%s/t_arrays/chan%d' % (type, det, bin)
        t_pass = [evt for evt in array_obj if (t_start <= evt <= t_stop)]
        counts = len(t_pass)
        cnt_array[bin] = counts

    return cnt_array 