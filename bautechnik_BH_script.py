# -*- coding: utf-8 -*-

# Copyright (C) 2024 Andreas Wiegel (andreas.wiegel@tum.de)
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU Affero General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU Affero General Public License for more details.

## Attribution
# Parts of this code are derived from the original source at 
#    http://academy.ifcopenshell.org/creating-a-simple-wall-with-property-set-and-quantity-information/
#    https://gitlab.com/CHGEOL/geol_bim-geo2ifc

import numpy as np
import pandas as pd
import ifcopenshell
from bautechnik_ifc_functions import (create_template, 
                                       create_guid, 
                                       create_ifcsite_with_survey_point,
                                       make_IfcBoreholes,
                                       make_plane,
                                       open_ifc, 
                                       template2tempfile)



geo_list = [["Ton","s_type",[0,1.5],[204, 0, 204]],
            ["Sand","s_type",[1.5,2.5],[255, 128, 0]],
            ["Kies","s_type",[2.5,3.5],[255, 255, 0]]]

geo_vertices_list = []
geo_faces_list = []

    
# -----------------------------------------------------------------------------

unc_list = [["Good","var",[0.0,0.15],[0, 255, 0]],
            ["Fair","var",[0.15,0.3],[255,255,0]],
            ["Poor","var",[0.3,0.45],[255,165,0]],
            ["Very Poor","var",[0.45,1],[255, 0, 0]]]

unc_vertices_list = []
unc_faces_list = []

# -----------------------------------------------------------------------------

SURVEY = pd.DataFrame([])
FROM_TO = pd.DataFrame([])

diameter_in_mm = 0.5
verticaloption = True


nbor = len(SURVEY) # number of Boreholes
IDS = []
x,y,z = [],[],[]
azm= 0
dip= 90
IDF =[]
lowerdep = []
soiltype =[]

# -----------------------------------------------------------------------------

gw_list =  [["p5m_hw40.shp",1940,'HOEHE',  [204, 255, 255]]]

gw_vertices_list = []
gw_faces_list = []


###############################################################################
 # # # # # # # # # # # # # # # CREATE IFC FILE # # # # # # # # # # # # # # # #
###############################################################################


file_name = "my_IFC_file.ifc"
creator = "Andreas Wiegel"
organization = "TUM_ZENTRUM_GEOTECHNIK"
project_description = "New Town Hall (Munich)"
project_name = "TUM_IFC"
project_origin = [0,0,0]
epsg = 32632
rotation = 0 

schema_version = "IFC4"
georeferencing = "IfcSite" 
# alternatively "IfcGeometricRepresentationContext" or "IfcMapConversion"
IfcSiteName = "New Town Hall"


###############################################################################

print("--setup file")
# create IFC-header
template = create_template(file_name, creator, project_name, organization,
    schema_version, project_description,  
    project_description, "", 
    project_origin, rotation, georeferencing, minimal = False)

# write template to temp-file
temp_handle, temp_filename = template2tempfile(template)
ifcfile, project, owner_history, context = open_ifc(temp_filename)

# create IfcSite: Site base point (IfcSite.ObjectPlacement) at the project 
# origin (0, 0, 0)
site, site_placement = create_ifcsite_with_survey_point(
    ifcfile, project, owner_history, context, project_origin, epsg,
    rotation, georeferencing
    )
## Importing .xls File and naming the column names   

geometric_representation_contexts = ifcfile.by_type(
    "IfcGeometricRepresentationContext")
subcontext = ifcfile.createIfcGeometricRepresentationSubContext(
    "Body",            # Context Identifier
    "Model",           # Target View (e.g., 'Body', 'Model', 'Plan')
    None, None, None, None,
    geometric_representation_contexts[0],        # Parent Context
    None, "MODEL_VIEW", None
    )

###############################################################################
# GEO MODEL
###############################################################################

print("--creating layers")

building = ifcfile.createIfcBuilding(
    create_guid(), owner_history,
    Name="Baugrundschichten",  # Set the name of the building
    Description="Baugrundschichten abgeleitet aus probabilistischen Baugrundmodell",  
    # Set the description of the building
    ObjectType="Building",  # Set the object type ("Building" in this case)
    ObjectPlacement=site_placement,  # Optional: IfcLocalPlacement to define 
                                     # the placement of the building
)
ifcfile.createIfcRelAggregates(
    create_guid(), owner_history, "Baugrundschichten", None, site, [building])
geo_layer = ifcfile.createIfcBuildingStorey(
    create_guid(), owner_history, "Baugrundschichten", None, None, 
    site_placement, None, None, 'COMPLEX', None
    )
ifcfile.createIfcRelAggregates(create_guid(), owner_history, 
                               "Baugrundschichten", None, building, 
                               [geo_layer])
    
proxy_list =[]
for vertices, faces, geo_sub_list in zip(
        geo_vertices_list, geo_faces_list, geo_list):
    vertices = [tuple(map(float, row)) for row in vertices - project_origin]
    faces = [tuple(map(int, row)) for row in faces + 1] #!!
    
    col_RGB = geo_sub_list[3]
    properties = {"Homogenbereich" : geo_sub_list[0]}
    layer_name = geo_sub_list[0]
    
    proxy = make_plane(
        ifcfile, layer_name, vertices, faces, col_RGB, properties,
        site_placement, context, owner_history, schema_version,
        geo_layer, site)
    proxy_list.append(proxy)
container_space = ifcfile.createIfcRelAggregates(
    create_guid(), owner_history, None, None, geo_layer, proxy_list)

    
unc_model_option = True

if unc_model_option:
    unc_layer = ifcfile.createIfcBuildingStorey(
        create_guid(), owner_history, "Prognosesicherheit", None, None, 
        site_placement, None, None, 'COMPLEX', None
        )
    ifcfile.createIfcRelAggregates(create_guid(), owner_history, 
                                   "Prognosesicherheit", None, building, 
                                   [unc_layer])
    proxy_list =[]
    for vertices, faces, unc_sub_list in zip(unc_vertices_list, unc_faces_list, 
                                             unc_list):
        vertices = [tuple(map(float, row)) for row in vertices - project_origin]
        faces = [tuple(map(int, row)) for row in faces + 1] #!!
        
        col_RGB = unc_sub_list[3]
        properties = {"Homogenbereich" : unc_sub_list[0]}
        layer_name = unc_sub_list[0]
        
        proxy = make_plane(
            ifcfile, layer_name, vertices, faces, col_RGB, properties,
            site_placement, context, owner_history, schema_version,
            unc_layer, site)
        proxy_list.append(proxy)
    container_space = ifcfile.createIfcRelAggregates(
        create_guid(), owner_history, None, None, unc_layer, proxy_list)

print("--finished creating layers")
###############################################################################   
# GROUNDWATER MODEL
###############################################################################  
print("--creating groundwatertables")
building = ifcfile.createIfcBuilding(
    create_guid(), owner_history,
    "Grundwasser",  # Set the name of the building
    None,  # Set the description of the building
    "Building",  # Set the object type ("Building" in this case)
    site_placement,  # Optional: IfcLocalPlacement to define the placement of the building
    None, None, 'COMPLEX', None, None, None)
ifcfile.createIfcRelAggregates(
    create_guid(), owner_history, "Grundwasser", 
    "Grundwasserstichtagsmessungen verschiedener Jahre", site, [building])
gw_layer = ifcfile.createIfcBuildingStorey(
    create_guid(), owner_history, "Grundwasserhöhen", None, None, 
    site_placement, None, None, 'COMPLEX', None
    )
ifcfile.createIfcRelAggregates(
    create_guid(), owner_history, "Grundwasserhöhen", None, building, 
    [gw_layer])

proxy_list =[]
for vertices, faces, gw_sub_list in zip(
        gw_vertices_list, gw_faces_list, gw_list):
    vertices = [tuple(map(float, row)) for row in vertices - project_origin]
    faces = [tuple(map(int, row)) for row in faces + 1] #!!
    col_RGB = gw_sub_list[3]
    properties = {"Stichtagsmessung [Jahr]" : gw_sub_list[1]}
    layer_name = gw_sub_list[1]
    proxy = make_plane(
        ifcfile, str(layer_name), vertices, faces, col_RGB, properties,
        site_placement, context, owner_history, schema_version, gw_layer, site)
    proxy_list.append(proxy)
    break
container_space = ifcfile.createIfcRelAggregates(
    create_guid(), owner_history, None, None, gw_layer, proxy_list)

print("--finished creating groundwatertables")
###############################################################################
# BOREHOLE MODEL
###############################################################################
print("--creating boreholes")
diameter_in_mm = 0.5

direction = ifcfile.createIfcDirection((0., 0., 1.))
cartesian_point_template = ifcfile.createIfcCartesianPoint((0., 0.))
axis_placement_template = ifcfile.createIfcAxis2Placement2D(
    cartesian_point_template, None)
circle_profile_template= ifcfile.createIfcCircleProfileDef(
    'AREA', None, axis_placement_template, 0.5 * diameter_in_mm)
colour_rgb_template = ifcfile.createIfcColourRgb(None, 0.7, 0.7, 0.7)
surface_style_shading_template = ifcfile.createIfcSurfaceStyleShading(
    colour_rgb_template)
surface_style_template = ifcfile.createIfcSurfaceStyle(
    None, 'BOTH', [surface_style_shading_template])
presentation_style_assignmen_template = ifcfile.createIfcPresentationStyleAssignment(
    (surface_style_template,))

building = ifcfile.createIfcBuilding(
    create_guid(), owner_history,
    "Baugrundaufschlüsse",  # Set the name of the building
    "Bohrdaten erhalten von LfU",  # Set the description of the building
    "Building",  # Set the object type ("Building" in this case)
    site_placement,  # Optional: IfcLocalPlacement 
    None, None, 'COMPLEX', None, None, None)

make_IfcBoreholes(
    x,y,z,nbor,IDS,FROM_TO,IDF,soiltype,lowerdep,verticaloption,
    ifcfile,site_placement,owner_history, circle_profile_template,
    axis_placement_template, direction,subcontext,SURVEY,
    site, building  )
print("--finished creating boreholes")
print("--writing to disc")

ifcfile.write(str(file_name))
