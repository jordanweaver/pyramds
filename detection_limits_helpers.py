def detection_limits_to_html(detection_limits):
    """Convert detection limits to HTML for viewer."""
    s = ""

    if detection_limits['lookup_selected'] != []:
        s += "<h3>Lookup Selected Peaks:</h3>\n"

        for row in sorted(detection_limits['lookup_selected']):
            s += "<b>Isotope:</b> {0}<br>\n".format(row[0])
            s += "<b>Peak Number:</b> {0}<br>\n".format(row[1])
            s += "<b>LD:</b> {0}<br><br>\n".format( row[2] )
            s += "<b>LC:</b> {0}<br><br>\n".format( row[3] )

        s += "<br>\n"

    if detection_limits['clicked_selected'] != {}:
        s += "<h3>Peaks Selected from Plot:</h3>\n"

        for key in sorted(detection_limits['clicked_selected'].keys()):
            s += "<b>Peak Index:</b> {0}<br>\n".format(key)
            s += "<b>LD:</b> {0}<br><br>\n".format( detection_limits['clicked_selected'][key][0] )
            s += "<b>LC:</b> {0}<br><br>\n".format( detection_limits['clicked_selected'][key][1] )

    return s


def detection_limits_to_tsv(detection_limits):
    """Convert detection limits to tab-separated value text file for saving."""
    s = ""

    if detection_limits['lookup_selected'] != {}:
        pass

    if detection_limits['clicked_selected'] != {}:
        s += ("Peaks Selected from Plot:\n"
             "Index\tLimit\n")

        for key in sorted(detection_limits['clicked_selected'].keys()):
            s += "{0}\t{1}\n".format(key, detection_limits['clicked_selected'][key])

    return s
