# Define all functions

def build_gamma(lib_name, fwhm_coeff, en_coeff):
    with open(lib_name, 'r') as ginput:
        gamma_lib = {}
        
        line_data = ginput.readline().split()
        while line_data:
            ZAID = line_data[0]
            name_str = line_data[1]
            energy = line_data[2]
            
            gamma_lib[ZAID] = [name_str, float(energy)]
            
            line_data = ginput.readline().split()
        
        sig_lookup = [{}]*4
        for c,chn in enumerate(en_coeff.keys()):
            for z,zaid in enumerate(gamma_lib.keys()):
                top_val = gamma_lib[zaid][1] - en_coeff[chn][0]
                bottom_val = en_coeff[chn][1]
                cent_chn = top_val/bottom_val
                width = fwhm_coeff[chn][0] + fwhm_coeff[chn][1] * cent_chn + fwhm_coeff[chn][2] * cent_chn * cent_chn
                
                left_marker = int(round(cent_chn - 1.75*width))
                right_marker = int(round(cent_chn + 1.75*width))
                
                sig_lookup[int(chn)][zaid] = [left_marker, right_marker]
            
    return gamma_lib, sig_lookup

def write_spec(spe_choices, file_series, spec_start, live_t, total_t, detector_choices, energy_max, spec_markers, enerfit, mca_cal, shape_cal, roi_list):
    
    # Write out the files for the desired detectors and coincidence patterns chosen
    # by the user. Using 'spe_choices' as spec_choices, the 'for' loops cycle through
    # all the user defined spectra to be generated. The format of the .Spe file is
    # described in "ORTEC-Sofware-File-Structure-Manual.pdf".
    
    for i,dc in enumerate(spe_choices.keys()):
        for j,ht in enumerate(spe_choices[dc]):
            
            # File names for the individual spectra contain the name of the .bin
            # file that was parsed followed by an identifier showing it was handled
            # by PYRAMDS and then which spectra it is in the choices made by the
            # user. The SPEC_ID found in the .Spe file lists the hit patterns used
            # in its generation.
            file_string = '%s_PYRAMDS_Det%s_Spec%s.Spe' % (file_series, dc, str(j+1))
            with open(file_string, 'w') as specout:
                field_width = 8
                
                # see: "ORTEC-Sofware-File-Structure-Manual.pdf" for more info
                spec_id = 'PYRAMDS %s [%s]' % (file_series, ' '.join(ht))
                spec_rem = 'DET# %s\r\nDETDESC# PYRAMDS\r\nAP# Maestro Version 6.04' % dc
                spe_date = spec_start.strftime("%m/%d/%Y %H:%M:%S")
                spe_meas = str(int(float(live_t[detector_choices[i]]))) + \
                            ' ' + str(int(float(total_t)))
                
                specout.write('$SPEC_ID:\r\n' + spec_id + '\r\n$SPEC_REM:\r\n' + spec_rem + \
                              '\r\n$DATE_MEA:\r\n' + spe_date + '\r\n$MEAS_TIM:\r\n' + spe_meas + \
                              '\r\n$DATA:\r\n0 ' + str(energy_max - 1) + '\r\n')
                
                # Begin writing out the values for each channel in this user-defined
                # set of markers.
                for en in range(energy_max):
                    str_number = str(spec_markers[int(dc)][j][en])
                    specout.write('%s\r\n' % (str_number.rjust(field_width)))
                
                roi_mark_list = ''
                for x in roi_list[dc]:
                    roi_mark_list += '%d %d\r\n' % (x[0],x[1])
                
                specout.write('$ROI:\r\n' + str(len(roi_list[dc])) + '\r\n' + roi_mark_list + \
                              '\r\n$PRESETS:\r\nNone\r\n0\r\n0\r\n$ENER_FIT:\r\n' + enerfit[dc] + \
                              '\r\n$MCA_CAL:\r\n' + mca_cal[dc] + \
                              '\r\n$SHAPE_CAL:\r\n' + shape_cal[dc] + '\r\n')