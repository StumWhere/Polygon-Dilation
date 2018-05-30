# -*- coding: utf-8 -*-
"""
This is a workflow to truncate deep and narrow corners and dilate narrow segments 
of soil polygons. This workflow is intened to assist in the QA/QC process of SSURGO

Five genveral steps:
    1) Smooth and Generalize features
    2) Identify deep corners and short narrow ends
    3) Truncate features which only share boundary with two MU's
    4) Identify narrow segements of polygons
    5) Dilate narrow segments and update

@author: alexander.stum

Last revised 5/30/18

PrintMsg and errorMsg were written by Steve Peaslee. At the momment these functions
are not utilized. When set-up to be used as an ArcToolbox Script, they will become
useful


"""

import arcpy,traceback,sys

#======= User Defined Variables  ==========
#The geodatabase with MUPOLYGON and SAPOLYGON
gdb = r'G:\GIS_specialist\Data\Geospatial\soils\FY18\Lubbock\Borden\workspace\Lubbock18_aks051618.gdb'
arcpy.env.workspace = gdb
#The main soil polygon layer
mupolygon = "MUPOLYGON"
#The survey area boundary, NOTE that no features touching boundary will be modified
sapolygon = "SAPOLYGON"
#selection query, specify if the entire survey area is not being evaluated
#otherwise set to None
query = "Update_Name IS NOT NULL"
#query = None

#======= Set Parameters  ==========
#Parameter for line smooting of polylines. Subjective, I like 25
smooth_tol = "20 Meters"
#Parameter for the generalization of polylines. Subjective, I like 2.
gen_tol = "2 Meters"
#Half the miniumum width standard, currently national standard is 38 m
min_width = 19
#All deep/narrow corners with Shape_Length less than this parameter will not 
#be treated. Leave it to the human
min_len = "300"
#Look ahead parameter, used by centerline function
f = 5


#======= Functions  ==========
def PrintMsg(msg, severity=0):
    # Adds tool message to the geoprocessor and prints to console window

    #Split the message on \n first, so that if it's multiple lines, a GPMessage will be added for each line
    try:
        for string in msg.split('\n'):
            #Add a geoprocessing message (in case this is run as a tool)
            if severity == 0:
                arcpy.AddMessage(string)

            elif severity == 1:
                arcpy.AddWarning(string)

            elif severity == 2:
                arcpy.AddMessage("    ")
                arcpy.AddError(string)

    except:
        pass


## ===================================================================================
def errorMsg():
    try:
        tb = sys.exc_info()[2]
        tbinfo = traceback.format_tb(tb)[0]
        theMsg = tbinfo + "\n" + str(sys.exc_type)+ ": " + str(sys.exc_value)
        PrintMsg(theMsg, 2)

    except:
        PrintMsg("Unhandled error in errorMsg method", 2)
        pass
    
    
## ===================================================================================   
def centerline(narrows, centerline, extend_to, f=5):
    """Calculates centerline between two more or less parrallel polylines (multipart)
    which represent the outline of a polygon.
    
    It finds the centerline by starting at the endpoint of a feature. It first
    looks at the first five vertices of the opposing part. If the nearest vertex
    is the last of these five, it will look at the next five. It does this for
    each vertex of the first part. At the end it verifies that the enddpoints of
    the second part were involved, if not it similarly finds the nearest vertices
    in part one to the uninvolved endbpoints of part 2.
    As the two part polyline outline is produced by buffering, the algorithm extends
    the centerline to the nearest negative buffer feature to extend the centerline
    beyond the inward bowing minicus.

    - **Parameters**::
    
        Name:      Type:           Description:
        narrows     String          Name of the polyline outline feature
        centerline  String          Name of the centerline feature created by function
        extend_to   String          Name of the polygon feature to which centerline will be stretched to
        f           Integer         Look ahead parameter, how many verticies to consider at a time
    
    - **Returns**::
        
        True        if successful
        
    .. note:: Requires ArcPy
`
    """  
    
    s_r = arcpy.Describe(narrows).spatialReference                  
    arcpy.CreateFeatureclass_management(arcpy.env.workspace, centerline, 'POLYLINE',\
                                        spatial_reference= s_r)
    
    def midpoint(pt1,pt2):
        x = pt1.X-(pt1.X-pt2.X)/2.0
        y = pt1.Y-(pt1.Y-pt2.Y)/2.0
        return arcpy.Point(x,y)
    
    sCur = arcpy.da.SearchCursor(narrows,'SHAPE@')
    iCur = arcpy.da.InsertCursor(centerline,'SHAPE@')
    for poly, in sCur:
        if poly and poly.partCount==2:
            line  = arcpy.Array()
            pi = 0   #current position in part 2
            p1 = poly.getPart(0)
            p2 = [p for p in poly.getPart(1)] #put in list to be able to slice
            p1c = p1.count
            dlast0 = ((p1[p1c-1].X-p2[0].X)**2+(p1[p1c-1].Y-p2[0].Y)**2)**.5
            d00= ((p1[0].X-p2[0].X)**2+(p1[0].Y-p2[0].Y)**2)**.5
            if dlast0 < d00:
                index=  range(p1c-1,-1,-1)
                line.add(midpoint(p2[0],p1[p1.count-1]))
            else:
                index=  range(p1c)
                line.add(midpoint(p2[0],p1[0]))           
    
            for i in index:
                pf = f-1 #relative forward position
                while pf == f-1:        #if the relative forward position is the same as the last relative position, continue looking
                    d = [((p1[i].X-p.X)**2+(p1[i].Y-p.Y)**2)**.5 \
                         for p in p2[pi:(pi+f) or None]]
                    pf = d.index(min(d)) #nearest forward point
                    pi += pf
    
                line.add(midpoint(p1[i],p2[pi]))
            if pi < len(p2)-1: #if last point of part 2 not included, extend
                line.add(midpoint(p1[i],p2[-1]))
            if (line[0].X==line[1].X) and (line[0].Y==line[1].Y):
                line.remove(0) #If first two points are identical, remove one
            #split line
    
            if line.count == 2:
                mid = midpoint(line[0],line[1])
                partA = arcpy.Array([line[0],mid])
                partB = arcpy.Array([mid,line[1]])
            elif line.count <2: print line.count
            else:
                partA = arcpy.Array([line[i] for i in range(line.count/2+line.count%2)])
                partB = arcpy.Array([line[i] for i in range(line.count/2,line.count)])
                
            iCur.insertRow([arcpy.Polyline(partA)])
            iCur.insertRow([arcpy.Polyline(partB)])
    
    del sCur, iCur

    #Tie centerlines to the nearest non-narrow segment of polygon
    arcpy.Near_analysis(centerline, extend_to, "", "LOCATION", "NO_ANGLE")
    
    uCur = arcpy.da.UpdateCursor(centerline,['SHAPE@','NEAR_X','NEAR_Y'])
    try:
        seg = uCur.next()
        while seg:
            partA = seg[0].getPart(0)
            if ((partA[0].X-seg[1])**2+(partA[0].Y-seg[2])**2)**.5 <= 18:
                partA.insert(0,arcpy.Point(seg[1],seg[2]))
            uCur.deleteRow()
            seg = uCur.next()
            partA.extend(seg[0].getPart(0))
            if ((partA[partA.count-1].X-seg[1])**2+ \
                (partA[partA.count-1].Y-seg[2])**2)**.5 <= 18:
                partA.add(arcpy.Point(seg[1],seg[2]))
            uCur.updateRow([arcpy.Polyline(partA),seg[1],seg[2]])
            seg = uCur.next()
    except StopIteration:
        return True

###******************* Part I *****************
#Step 0
arcpy.FeatureToLine_management(mupolygon, "MU_lines", "", "NO_ATTRIBUTES")
arcpy.Dissolve_management("MU_lines", "MU_lines_dis", "",\
                          "", "SINGLE_PART", "DISSOLVE_LINES")

arcpy.FeatureToPoint_management(mupolygon, "MU_point", "INSIDE")

arcpy.MakeFeatureLayer_management('MU_lines_dis', 'MU_lines_select')
arcpy.Delete_management('MU_lines_dis')
arcpy.Delete_management("MU_lines")
if query:
    arcpy.MakeFeatureLayer_management(mupolygon, 'MU_select', query)
    arcpy.SelectLayerByLocation_management("MU_lines_select","SHARE_A_LINE_SEGMENT_WITH",\
    'MU_select',"#","NEW_SELECTION")
    arcpy.SelectLayerByLocation_management("MU_lines_select","SHARE_A_LINE_SEGMENT_WITH",sapolygon,\
                                       '#',"REMOVE_FROM_SELECTION")
else:
    arcpy.SelectLayerByLocation_management("MU_lines_select","SHARE_A_LINE_SEGMENT_WITH",sapolygon,\
                                       '#',"NEW_SELECTION","INVERT")   

arcpy.SmoothLine_cartography("MU_lines_select", "MU_lines_gen", "PAEK",\
                             "25 Meters","FIXED_CLOSED_ENDPOINT", "NO_CHECK")
arcpy.Generalize_edit("MU_lines_gen", "2 Meters")
arcpy.SelectLayerByAttribute_management("MU_lines_select","SWITCH_SELECTION")
arcpy.Merge_management("MU_lines_select;MU_lines_gen", "MU_lines_gen_merge")

arcpy.FeatureToPolygon_management("MU_lines_gen_merge", "MU_gen",\
                                      "", "ATTRIBUTES", "MU_point")
arcpy.Delete_management("MU_lines_gen_merge")
arcpy.Delete_management("MU_point")
arcpy.Delete_management("MU_lines_gen")
if query:
    arcpy.MakeFeatureLayer_management('MU_gen', 'MU_gen_select', query)
    #Step 1
    #Expect to see several warnings of features dissappearing
    arcpy.Buffer_analysis("MU_gen_select", "MU_negbuff19", "-"+str(min_width)+" Meters")
        #Step 2
    arcpy.Buffer_analysis("MU_negbuff19", "MU_rebuff19", str(min_width+.1)+" Meters")
    
    #Step 3
    arcpy.Erase_analysis("MU_gen_select", "MU_rebuff19", \
                         "too_narrow", cluster_tolerance="0.1 Meters")
else:
    #Step 1
    #Expect to see several warnings of features dissappearing
    arcpy.Buffer_analysis("MU_gen", "MU_negbuff19", "-"+str(min_width)+" Meters")
    
    #Step 2
    arcpy.Buffer_analysis("MU_negbuff19", "MU_rebuff19",\
                          str(min_width+.1)+" Meters")
    
    #Step 3
    arcpy.Erase_analysis("MU_gen", "MU_rebuff19", "too_narrow", \
                         cluster_tolerance="0.1 Meters")

arcpy.Delete_management("MU_negbuff19")
arcpy.Delete_management("MU_rebuff19")

#Step 4
arcpy.MultipartToSinglepart_management("too_narrow", "too_narrow_sing")

#Step 5
arcpy.MakeFeatureLayer_management('too_narrow_sing','too_narrow_select',\
                                  'Shape_Length >=30 AND Shape_Length <=300')
    
#Step 6
arcpy.Buffer_analysis("too_narrow_select", "corner_buff_4", "0.4 Meters")

#Step 7 Intersect is necessary to produce a polygon county on both sides of narrows to equal 3
#Tolerance set to clean out polygons with 0 area
arcpy.Intersect_analysis("MU_gen #;corner_buff_4 #", \
                         "corner_inter", "ALL", "0.05 Meters")

#Step 8
arcpy.MultipartToSinglepart_management("corner_inter", "corner_inter_sing")

#Step 9
arcpy.SpatialJoin_analysis("corner_buff_4", "corner_inter_sing",\
            "corner_sum","JOIN_ONE_TO_ONE", "KEEP_ALL", "", "INTERSECT")

#Step 10
arcpy.MakeFeatureLayer_management('corner_sum','corner_sum_select',"Join_Count =2")

#Step 11
arcpy.SelectLayerByLocation_management('corner_sum_select',"CROSSED_BY_THE_OUTLINE_OF",\
                                       sapolygon,'#',"NEW_SELECTION",'INVERT')

#Step 12
arcpy.SelectLayerByLocation_management('too_narrow_select',"HAVE_THEIR_CENTER_IN"\
                                       ,'corner_sum_select','#',"NEW_SELECTION")

#Step 13
arcpy.Update_analysis("MU_gen", "too_narrow_select", "MU_cornered", "BORDERS", "0.1 Meters")

#Step 14
arcpy.MakeFeatureLayer_management("MU_cornered","MU_cornered_select")
arcpy.SelectLayerByLocation_management("MU_cornered_select","COMPLETELY_WITHIN",\
                                       'corner_sum_select')
arcpy.CopyFeatures_management("too_narrow_select",'truncated')
#Step 15 Eliminate
arcpy.Eliminate_management("MU_cornered_select", "MU_decornered", "LENGTH")

arcpy.Delete_management("too_narrow")
arcpy.Delete_management("too_narrow_sing")
arcpy.Delete_management("corner_buff_4")
arcpy.Delete_management("corner_inter")
arcpy.Delete_management("corner_inter_sing")
arcpy.Delete_management("corner_sum")
arcpy.Delete_management("MU_cornered")


###******************* Part II *****************


#*********************************************************
if query:
    arcpy.MakeFeatureLayer_management("MU_decornered", 'MU_decornered_select', query)
    #Step 1
    arcpy.Buffer_analysis("MU_decornered_select", "MU_2negbuff19","-19 Meters")
    
    #Step 2
    arcpy.Buffer_analysis("MU_2negbuff19", "MU_2rebuff19", "19.4 Meters")
    
    #Step 3
    arcpy.Erase_analysis('MU_decornered_select', "MU_2rebuff19", \
                         "too_narrow2", "0.1 Meters")
else:
    #Step 1
    arcpy.Buffer_analysis("MU_decornered", "MU_2negbuff19","-19 Meters")
    
    #Step 2
    arcpy.Buffer_analysis("MU_2negbuff19", "MU_2rebuff19", "19.4 Meters")
    
    #Step 3
    arcpy.Erase_analysis("MU_decornered", "MU_2rebuff19",\
                         "too_narrow2", "0.1 Meters")
arcpy.Delete_management("MU_2rebuff19")

#Step 4
arcpy.MultipartToSinglepart_management("too_narrow2", "too_narrow_sing2")

#Step 5
arcpy.FeatureToLine_management("too_narrow_sing2", "narrow_outlines")

#Step 6
arcpy.Intersect_analysis("narrow_outlines #;MU_decornered #", \
                         "narrow_outlines_inter", "ALL", "0.5 Meters")

#Step 7
arcpy.MakeFeatureLayer_management("narrow_outlines_inter","narrow_outlines_select",\
                                  "MUSYM <> MUSYM_1")
#Step 8 
arcpy.Dissolve_management("narrow_outlines_select", "casings", "FID_too_narrow_sing2")

#Step 9
arcpy.Densify_edit("casings", "DISTANCE", "30 Meters")

#Step 10 
centerline('casings','centerlines','MU_2negbuff19')

#Step 11
arcpy.Identity_analysis("centerlines", "MU_decornered", "centerlines_attributed", "ALL")

#Step 12
arcpy.Buffer_analysis("centerlines_attributed","inserts",\
                      "19 Meters", "FULL", "ROUND", "LIST", "FID_MU_decornered;MUSYM")

#Step 13
arcpy.MakeFeatureLayer_management('inserts','inserts_select')
arcpy.SelectLayerByLocation_management('inserts_select',"CROSSED_BY_THE_OUTLINE_OF",\
                                       sapolygon,'#',"NEW_SELECTION",'INVERT')
arcpy.CopyFeatures_management("inserts_select",'dilations')

#Step 14
arcpy.Update_analysis("MU_decornered", "dilations", "MU_inserted", \
                      "BORDERS", "0.1 Meters")

#Step 15
arcpy.Dissolve_management("MU_inserted", "MU_dilated", "MUSYM", "", \
                          "SINGLE_PART", "DISSOLVE_LINES")

arcpy.Delete_management("MU_2negbuff19")
arcpy.Delete_management("too_narrow2")
arcpy.Delete_management("too_narrow_sing2")
arcpy.Delete_management("narrow_outlines")
arcpy.Delete_management("narrow_outlines_inter")
arcpy.Delete_management("centerlines_attributed")
arcpy.Delete_management("inserts")
arcpy.Delete_management("MU_inserted")


### Densify 31.75 degrees to replace arcs with appropriate density
arcpy.Densify_edit("MU_dilated", "ANGLE", max_angle="31.75")

arcpy.FeatureToPoint_management("MU_decornered","MUpoints","INSIDE")

arcpy.SpatialJoin_analysis("MU_dilated", "MUpoints", "MU_final","JOIN_ONE_TO_ONE")
if query:
    arcpy.MakeFeatureLayer_management("MU_final", 'MU_final_select', query)
    arcpy.Buffer_analysis("MU_final_select", "MU_3negbuff19","-18.85 Meters")
    
    #Step
    arcpy.Buffer_analysis("MU_3negbuff19", "MU_3rebuff19", "19 Meters")
    
    #Step
    arcpy.Erase_analysis("MU_final_select", "MU_3rebuff19",\
                         "too_narrow3", "0.1 Meters")
else:
    arcpy.Buffer_analysis("MU_final", "MU_3negbuff19","-18.85 Meters")
    
    #Step
    arcpy.Buffer_analysis("MU_3negbuff19", "MU_3rebuff19", "19.25 Meters")
    
    #Step
    arcpy.Erase_analysis("MU_final", "MU_3rebuff19",\
                         "too_narrow3", "0.1 Meters")
    #Step
    arcpy.MultipartToSinglepart_management("too_narrow3", "Remaining_Narrows")

arcpy.Delete_management("MU_3negbuff19")
arcpy.Delete_management("MU_3rebuff19")
arcpy.Delete_management("too_narrow3")
arcpy.Delete_management("too_narrow3")

#https://www.ian-ko.com/ET_GeoWizards/WhitePapers/gw_smooth_generalize.htm
#http://desktop.arcgis.com/en/arcmap/10.3/analyze/modelbuilder/the-in-memory-workspace.htm#ESRI_SECTION1_1BDEBF62493E489FA7A1CD7E4E951D5A