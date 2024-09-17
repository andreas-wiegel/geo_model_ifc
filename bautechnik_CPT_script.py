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
import ifcopenshell
from bautechnik_ifc_functions import (create_template, 
                           create_guid, 
                           create_ifcsite_with_survey_point,
                           make_plane,
                           open_ifc, 
                           template2tempfile)
import colorsys # The percentage will map to a hue value in the HSV color space

nz_tot_list = [299]

Hi_Vis_col = [  [170, 85, 0],    # Humos:
                [170, 0, 255],   # Clay:
                [170, 157, 13],   # Silt:
                [255, 193, 8],   # SandM:
                [255, 116, 3],   # Sand:
                [255, 238, 1]]   # Gravel:

geo_list = [["Organic soils – clay","sbt_zone_etype",
             [1.5,2.5],Hi_Vis_col[0],"ORGANIC"],    # TSD 5-2-6
            ["Clays – silty clay to clay","sbt_zone_etype",
             [2.5,3.5],Hi_Vis_col[1],"xxxCLAY"],    # TSD 14-5-4 -> 14-4-4
            ["Silt mixtures – clayey silt to silty clay","sbt_zone_etype",
             [3.5,4.5],Hi_Vis_col[2],"xxxSILT"],    # TSD 14-5-4
            ["Sand mixtures – silty sand to sandy silt","sbt_zone_etype",
             [4.5,5.5],Hi_Vis_col[3],"SIxSAND"], # TSD NA
            ["Sands – clean sand to silty sand","sbt_zone_etype",
             [5.5,6.5],Hi_Vis_col[4],"xxxSAND"], # TSD 6-6-2
            ["Gravelly sand to dense sand","sbt_zone_etype",
             [6.5,7.5],Hi_Vis_col[5],"GRxSAND"]] # TSD 2-6-1

geo_vertices_list = []
geo_faces_list = []

unc_list = [["Fair","sbt_zone_var",[0,0.3],[255,255,0],"xxxFAIR"],
            ["Poor","sbt_zone_var",[0.3,0.5],[255,165,0],"xxxPOOR"],
            ["Very Poor","sbt_zone_var",[0.5,1],[255, 0, 0],"VRYPOOR"]]

unc_vertices_list = []
unc_faces_list = []
cpt_list = []
thresh = 20
azm = 180
###############################################################################
 # # # # # # # # # # # # # # # CREATE IFC FILE # # # # # # # # # # # # # # # #
###############################################################################


file_name = "my_IFC_file.ifc"
creator = "Andreas Wiegel"
organization = "TUM_ZENTRUM_GEOTECHNIK"
project_description = "CPT"
project_name = "TUM_IFC"
project_origin = [0,0,0]
epsg = 32632
rotation = 0 
schema_version = "IFC4"
georeferencing = "IfcSite" 
# alternatively "IfcGeometricRepresentationContext" or "IfcMapConversion"
IfcSiteName = "CPT_TEST"

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

# create IfcSite: 
# Site base point (IfcSite.ObjectPlacement) at the project origin (0, 0, 0)
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

GIM = ifcfile.createIfcBuilding(
    create_guid(), owner_history,
    Name="Geotechnical Information Model",  # Set the name of the building
    Description="Geotechnical Information",  # Set the description of the 
                                             # building
    ObjectType="Building",  # Set the object type ("Building" in this case)
    ObjectPlacement=site_placement,  # Optional: IfcLocalPlacement to define 
                                     # the placement of the building
)
ifcfile.createIfcRelAggregates(create_guid(), owner_history, "Factual Data",
                               None, site, [GIM])

# -----------------------------------------------------------------------------

GDM = ifcfile.createIfcBuildingStorey(
    create_guid(), owner_history, "GeoDocu-Model", None, None, site_placement, 
    None, None, 'COMPLEX', None
    )
ifcfile.createIfcRelAggregates(create_guid(), owner_history, "GDM", None, GIM, 
                               [GDM])

# -----------------------------------------------------------------------------

for cpt,k in zip(cpt_list,range(len(cpt_list))):
    CPT = ifcfile.createIfcSpace(
        create_guid(), owner_history, "CPT", None, None, site_placement, None, 
        None, 'COMPLEX', None
        )
    ifcfile.createIfcRelAggregates(create_guid(), owner_history, "CPT", None, 
                                   GDM, [CPT])
    
    proxy_list =[]
    vertical_resolution = cpt[0,2]-cpt[1,2]
    for j, interval in zip(range(len(cpt)),cpt):
        interval = interval.astype(float)
        x0,y0,z,qc = interval[0], interval[1], interval[2], interval[3]
        var = qc*0.5
        upperdepth = z + vertical_resolution
        lowerdepth = z 
        
        fx = np.sin(azm/180*np.pi)
        fy = np.cos(azm/180*np.pi)
        uppery0 = tuple([x0,y0,upperdepth])
        lowery0 = tuple([x0,y0,lowerdepth])
        uppery1 = tuple([x0 + fx*var,y0 + fy*var,upperdepth])
        lowery1 = tuple([x0 + fx*var,y0 + fy*var,lowerdepth])
        vertices = [uppery0,lowery0,uppery1,lowery1]
        vertices = np.round(vertices,3)
        vertices = [
            tuple(map(float, row)) for row in vertices - project_origin]
        faces = [(1,2,3),(2,3,4)]
        
        percentage = qc / thresh
        if ~(percentage > 0): continue
        if percentage > 1: percentage = 1
        
        hue = (1-percentage) * 0.7  
        # 0.83 to cover the range of colors in a rainbow 
        # (excluding pinks at the end of the spectrum)
        
        r, g, b = colorsys.hsv_to_rgb(hue, 1.0, 1.0) # Convert HSV to RGB
        col_RGB = [int(r*255),int(g*255),int(b*255)]
        properties = {"q_c [MPa]" : interval[3],
                      "f_s [MPa]" : interval[4],
                      "u_2 [Mpa]" : interval[5]
                      }

        layer_name = f"{np.round(z,2):06.2f}"
        
        proxy = make_plane(ifcfile, layer_name, vertices, faces, col_RGB, 
                           properties, site_placement, context, owner_history, 
                           schema_version, CPT, site)
        proxy_list.append(proxy)

    property_values = [ifcfile.createIfcPropertySingleValue(
                    "DAUB ObjectID", 
                    None, 
                    ifcfile.create_entity(
                        "IfcText", f"BAUxx-TUMZG-PROJ1-xxx-GIM-GDM-xxx-xxx-CPT-CPT-{int(k+1):03d}-KM{int(x0):05d}-{int(y0):04d}-0001"), 
                    None)
                    ] 
    property_set = ifcfile.createIfcPropertySet(
        create_guid(), owner_history, "pset_bhl_common", None, property_values)
    ifcfile.createIfcRelDefinesByProperties(
        create_guid(), owner_history, None, None, [CPT], property_set) 

    RelContainedInSpatialStructure = ifcfile.createIfcRelContainedInSpatialStructure(
        create_guid(), owner_history, None,
        None, proxy_list, CPT
        )
    
print("--finished creating CPTs")

###############################################################################

IPM = ifcfile.createIfcBuildingStorey(
    create_guid(), owner_history, "Interpretative Model", None, None, 
    site_placement, None, None, 'COMPLEX', None
    )
ifcfile.createIfcRelAggregates(create_guid(), owner_history, 
                               "Interpretative Model", None, GIM, [IPM])

# -----------------------------------------------------------------------------

GEO = ifcfile.createIfcSpace(
    create_guid(), owner_history, "Geo", None, None, site_placement, None,
    None, 'COMPLEX', None
    )
ifcfile.createIfcRelAggregates(create_guid(), owner_history, "Geo", None, IPM, 
                               [GEO])
    
proxy_list =[]
for vertices, faces, sub_list in zip(geo_vertices_list, geo_faces_list, 
                                     geo_list):
    vertices = [tuple(map(float, row)) for row in vertices - project_origin]
    faces = [tuple(map(int, row)) for row in faces + 1] #!!
    
    col_RGB = sub_list[3]
    properties = {"Soil Behavior Type" : sub_list[0],
                  "DAUB ObjectID": f"BAUxx-TUMZG-PROJ1-xxx-GIM-IPM-xxx-xxx-GEO-GEO-SBT-{sub_list[4]}-xxxx-xxxx"}
    layer_name = sub_list[0]
    
    proxy = make_plane(ifcfile, layer_name, vertices, faces, col_RGB, 
                       properties, site_placement, context, owner_history, 
                       schema_version, IPM, site)
    proxy_list.append(proxy)
RelContainedInSpatialStructure = ifcfile.createIfcRelContainedInSpatialStructure(
    create_guid(), owner_history, None, None, proxy_list, GEO)

print("--finished creating geo layers")

# -----------------------------------------------------------------------------

UNC = ifcfile.createIfcSpace(
    create_guid(), owner_history, "Uncertainty", None, None, site_placement,
    None, None, 'COMPLEX', None
    )
ifcfile.createIfcRelAggregates(create_guid(), owner_history, "Uncertainty", 
                               None, IPM, [UNC])

proxy_list =[]
for vertices, faces, unc_sub_list in zip(unc_vertices_list, unc_faces_list, 
                                         unc_list):
    vertices = [tuple(map(float, row)) for row in vertices - project_origin]
    faces = [tuple(map(int, row)) for row in faces + 1] #!!
    
    col_RGB = unc_sub_list[3]
    properties = {
        "Prediction Variability" : unc_sub_list[0],
        "DAUB ObjectID": f"BAUxx-TUMZG-PROJ1-xxx-GIM-IPM-xxx-xxx-GEO-UNC-VAR-{unc_sub_list[4]}-xxxx-xxxx"}
    layer_name = unc_sub_list[0]
    
    proxy = make_plane(ifcfile, layer_name, vertices, faces, col_RGB, 
                       properties, site_placement, context, owner_history, 
                       schema_version, UNC, site)
    proxy_list.append(proxy)
RelContainedInSpatialStructure = ifcfile.createIfcRelContainedInSpatialStructure(
    create_guid(), owner_history, None, None, proxy_list, UNC)

print("--finished creating unc layers")

print("--writing to disc")

ifcfile.write(str(file_name))
