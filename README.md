# Polygon-Dilation
Truncates and Dilates portions of soil map unit polygons
--------------------------------------------------------


This is a workflow to truncate deep and narrow corners and dilate narrow segments 
of soil polygons. This workflow is intened to assist in the QA/QC process of SSURGO

Five genveral steps:
    1.  Smooth and Generalize features
    2.  Identify deep corners and short narrow ends
    3.  Truncate features which only share boundary with two MU's
    4.  Identify narrow segements of polygons
    5.  Dilate narrow segments and update

This is a working prototype. I am hoping folks and try it out and provide feedback. Needs to be formally parameterized and set-up for usage as an ArcToolbox script. Perhaps for inclusion in the SSURGO_QA tool.
