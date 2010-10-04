def detection_limits_to_html(detection_limits):
    s = ""

    if detection_limits['lookup_selected'] != {}:
        pass

    if detection_limits['clicked_selected'] != {}:
        s += "<h3>Peaks Selected from Plot:</h3>\n"

        for key in sorted(detection_limits['clicked_selected'].keys()):
            s += "<b>Peak Index:</b> {0}<br>\n".format(key)
            s += "<b>Detection Limit:</b> {0}<br><br>\n".format( detection_limits['clicked_selected'][key])
    return s
