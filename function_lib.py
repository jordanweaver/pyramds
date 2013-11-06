#!/usr/bin/env python

# Definitions of all custom functions to be used by PYRAMDS

import numpy as np
import datetime
import time
import tables as tb

from pyramds_cfg import *


def get_bin_info(file_path):

    with open(file_path + '.ifm','rU') as finfo:
        # Read in the entire .ifm file as a series of lines. From info_str_list, the
        # necessary sections can be acquired for various data parameters. This works
        # granted the .ifm file never changes its format.
        info_str_list = finfo.readlines()

        date_str = info_str_list[1][23:-2]
        run_start_time = datetime.datetime(*time.strptime(date_str, \
                            '%I:%M:%S %p %a, %b %d, %Y')[0:6])

        # Total time = real time
        total_time = info_str_list[6].split()[3]

        # Live time for each detector channel
        live_time =[0]*4
        for channel in range(4):
            live_time[channel] = info_str_list[9 + channel].split()[2]

        times = {
            'start' :   run_start_time,
            'total' :   total_time,
            'live'  :   live_time
        }

        bufheadlen = int(info_str_list[33].split()[1])
        eventheadlen = int(info_str_list[34].split()[1])
        chanheadlen = 2 #int(info_str_list[35].split()[1])

    return times, bufheadlen, eventheadlen, chanheadlen

def calc_det_limit(LM, RM, spec_array):
    """
    Calculates detection limits by first averaging counts about left/right
    markers and then using that to determine average counts per channel in
    background, mu_b. From there, critical level (el_c) and detection limit
    (el_d) are determined to a 95% confidence level.

    """

    avg_pm = 3
    N_chn = abs(RM - LM) + 1
    LM_avg = np.average(spec_array[(LM - avg_pm):(LM + avg_pm + 1)])
    RM_avg = np.average(spec_array[(RM - avg_pm):(RM + avg_pm + 1)])

    mu_b = N_chn * (LM_avg + RM_avg) * 0.5
    mu_b_sqroot = np.sqrt(mu_b)

    el_c = 2.325 * mu_b_sqroot
    el_d = 2.706 + 4.653 * mu_b_sqroot

    return el_c, el_d

def calc_det_limit_sel(chn, marker, en_coeff, fwhm_coeff, spec_array):

    cent_en = en_coeff[0] + en_coeff[1] * marker
    width = fwhm_coeff[0] + fwhm_coeff[1] * marker + fwhm_coeff[2] * marker * marker

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

def marker2energy(marker, en_coeff):
    energy = en_coeff[0] + en_coeff[1] * marker
    return energy

def write_spec(file_path, data_title):

    f = tb.openFile(file_path, 'r')
    spectra = f.root.spectra

    stl = f.root.times.start.read()
    startt = datetime.datetime(stl[0],stl[1],stl[2],stl[3],stl[4],stl[5])
    times = {
        'live'      :   list(f.root.times.live.read()),
        'start'     :   startt,
        'total'     :   f.root.times.total.read()[0]
    }

    for spec_type in spectra:
        for group in [x for x in spec_type if (x.name[-4:] == 'spec')]:
            title_list = group.title.split()
            spec_type = title_list[0]
            det_no = title_list[-1]
            group_title = spec_type + '-Det' + det_no
            file_string = 'data/%s-%s_PYRAMDS.Spe' % (data_title, group_title)

            spec_markers = group[-1]

            with open(file_string, 'w') as specout:
                field_width = 8

                # see: "ORTEC-Sofware-File-Structure-Manual.pdf" for more info
                spec_id = 'PYRAMDS %s %s' % (data_title, group_title)
                spec_rem = 'DET# %s\r\nDETDESC# PYRAMDS\r\nAP# Maestro Version 6.04' % det_no
                spe_date = times['start'].strftime("%m/%d/%Y %H:%M:%S")
                spe_meas = str(int(float(times['live'][int(det_no)]))) + ' ' + str(int(float(times['total'])))

                specout.write('$SPEC_ID:\r\n' + spec_id + '\r\n$SPEC_REM:\r\n' + spec_rem +
                              '\r\n$DATE_MEA:\r\n' + spe_date + '\r\n$MEAS_TIM:\r\n' + spe_meas +
                              '\r\n$DATA:\r\n0 ' + '8191' + '\r\n')

                # Begin writing out the values for each channel in this user-defined
                # set of markers.
                for en in range(8192):
                    str_number = str(spec_markers[en])
                    specout.write('%s\r\n' % (str_number.rjust(field_width)))

                roi_mark_list = []

                specout.write('$ROI:\r\n' + str(len(roi_mark_list)) +
                              '\r\n$PRESETS:\r\nNone\r\n0\r\n0\r\n$ENER_FIT:\r\n' + enerfit[det_no] +
                              '\r\n$MCA_CAL:\r\n' + mca_cal[det_no] +
                              '\r\n$SHAPE_CAL:\r\n' + shape_cal[det_no] + '\r\n')

    f.close()
