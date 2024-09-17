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

from math import radians, sin, cos
import ifcopenshell
import time
import uuid
import logging
import tempfile
import numpy as np
from tqdm import tqdm
from pyproj import CRS, Transformer

module_logger = logging.getLogger('geol_bim.ifc_utils.read_write')

# Standard parameters if not specified
O = (0., 0., 0.)  # Location
X = (1., 0., 0.)  # Direction of X-Axis
Y = (0., 1., 0.)  # Direction of Y-Axis, for True North
Z = (0., 0., 1.)  # Direction of Z-Axis

# Creates GUID
create_guid = lambda: ifcopenshell.guid.compress(uuid.uuid1().hex)

def template2tempfile(template):
    """
    Writes the template to a tempfile.
    :param template: str, IFC header and first lines of the DATA section
    :return:
        temp_handle:    handle to the open file
        temp_filename:  filename
    """
    temp_handle, temp_filename = tempfile.mkstemp(suffix='.ifc')
    with open(temp_filename, 'w') as file:
        file.write(template)
    module_logger.debug('Writing tempfile...')

    return temp_handle, temp_filename

def create_ifcaxis2placement3D(ifcfile, point=O, dirZ=Z, dirX=X):
    """
    Creates an IfcAxis2Placement3D from Location, Axis and RefDirection 
    specified as Python tuples
    :param ifcfile: Instance of ifcopenshell
    :param point: tuple, (x,y,z), Location
    :param dirZ: tuple, Direction of Z-axis
    :param dirX: tuple, Direction of X-Axis
    :return: IfcAxis2Placement3D
    """
    # Convert integer coordinates to float
    point = tuple(float(i) for i in point)
    dirZ = tuple(float(i) for i in dirZ)
    dirX = tuple(float(i) for i in dirX)

    point = ifcfile.createIfcCartesianPoint(point)
    # creates less IfcCartesianPoints
    if dirZ == (0., 0., 1.) and dirX == (1., 0., 0.):  
        axis2placement = ifcfile.createIfcAxis2Placement3D(point)
    else:
        dirZ = ifcfile.createIfcDirection(dirZ)
        dirX = ifcfile.createIfcDirection(dirX)
        axis2placement = ifcfile.createIfcAxis2Placement3D(point, dirZ, dirX)

    return axis2placement

def create_ifclocalplacement(ifcfile, point=O, dirZ=Z, dirX=X, 
                             relative_to=None):
    """
    Creates an IfcLocalPlacement from Location, Axis and RefDirection, and 
    relative placement
    :param ifcfile: Instance of ifcopenshell
    :param point: tuple, (x,y,z), Location
    :param dirZ: tuple, Direction of Z-axis
    :param dirX: tuple, Direction of X-Axis
    :param relative_to: IfcLocalPlacement
    :return: IfcLocalPlacement
    """
    axis2placement = create_ifcaxis2placement3D(ifcfile, point, dirZ, dirX)
    ifclocalplacement = ifcfile.createIfcLocalPlacement(relative_to, 
                                                        axis2placement)

    return ifclocalplacement

def open_ifc(filename, minimal=False):
    """
    Open file with ifcopenshell, get references
    :param filename: str, for example temp_filename
    :param minimal: Bool
                    if True, open an IFC without OwnerHistory, 
                    GeometricRepresentationContext, Units
    :return: entities of ofcopenshell
        ifcfile:        representing the file
        project:        representing IfcProject (minimal=False)
        owner_history:  representing IfcOwnerHistory (minimal=False)
        context:        representing IfcContext (minimal=False)
    """
    # Open file with ifcopenshell, get references
    ifcfile = ifcopenshell.open(filename)
    module_logger.debug('Loading tempfile, obtain references...')
    if minimal:
        project = ifcfile.by_type("IfcProject")[0]
        return ifcfile, project
    else:
        # Obtain references to instances defined in template
        owner_history = ifcfile.by_type("IfcOwnerHistory")[0]
        project = ifcfile.by_type("IfcProject")[0]
        context = ifcfile.by_type("IfcGeometricRepresentationContext")[0]

        return ifcfile, project, owner_history, context

def rotation2xdir(rot, mode='deg'):
    """
    Creates a 3D coordinate-tuple of the X-direction
    :param rot: Rotation (counterclockwise)
    :param mode: 'deg' or 'rad'
    :return: Coordinate-tuple of x-direction, (1,0,0) for rot=0
    """
    if mode == 'deg':
        rotation = radians(rot)
    elif mode == 'rad':
        rotation = rot
    else:
        print("Rotation: in degrees ('deg') or radians ('rad')")
        raise ValueError

    # Cosmetics: if no rotation -> 1,0,0
    if rot == 0:
        xdirection = (1., 0., 0.)
    else:
        xdirection = (round(cos(rotation), 5), round(sin(rotation), 5), 0.)

    return xdirection


def rotation2ydir(rot, mode='deg'):
    """
    Creates a 3D Coordinate-Tuple of the Y-Direction
    Only used to calculate True North (for Axis2Placement, X-/Z-Directions are used)
    :param rot: Rotation (counterclockwise)
    :param mode: 'deg' or 'rad'

    :return: Coordinate-tuple of Y-direction (1,0,0) for rot=0
    """
    if mode == 'deg':
        rotation = radians(rot)
    elif mode == 'rad':
        rotation = rot
    else:
        print("Rotation: in degrees ('deg') or radians ('rad')")
        raise ValueError

    # Cosmetics: if no rotation -> 0,1,0
    if rot == 0:
        tn = (0., 1., 0.)
    else:
        tn = (round(sin(rotation), 5), round(cos(rotation), 5), 0.)

    return tn



def create_template(
        filename, author, project_name, organization, schema_version, 
        file_descr='', project_descr=None, authorization='', origin=O, rotation=0.,
        georeferencing='IfcGeometricRepresentationContext', minimal = False
        ):
    """
    Creates the template for an ifc_utils-file
    :param filename: str
    :param author: str
    :param project_name: str
    :param organization: str
    :param schema_version: 'IFC2X3' or 'IFC4'
    :param minimal: Bool
                        if True, write a minimal Header without OwnerHistory, 
                        GeometricRepresentationContext, Units
    :param file_descr: str
                        Optional File Description for STEP-Header
    :param project_descr: str
                        Optional Project Description for IfcProject
    :param authorization: str
                        Optional Auth-orization for STEP-Header
    :param origin: Coordinate-Tuple of the project reference point / survey 
                   point in the global reference system.
    :param rotation: Rotation FROM true North TO positive y-Axis
                     Degrees, counterclockwise (right-handed Cartesian 
                                                coordinate system)
                     Rotation center == origin
                     Rotation of the underlying project coordinate system,
                     relative to the true north (geographic northing direction).

    :param georeferencing:  '': No Georef
                            'IfcSite': Georef with IfcSite 
                                        (Not within Template)
                            'IfcGeometricRepresentationContext': Georef with 
                                        IfcGeometricRepresentationContext
                            'IfcMapConversion': Georef with 
                                        IfcGeometricRepresentationContext AND 
                                        IfcMapConversion



    :return: template
                        String with STEP-header and first section of an 
                        IFC-File
    """

    app, app_version = 'IfcOpenShell', ifcopenshell.version
    project_globalid = create_guid()
    # project_globalid = 12345
    t = time.time()
    t_string = time.strftime("%Y-%m-%dT%H:%M:%S", time.localtime(t))
    x_direction = rotation2xdir(rotation)
    true_north = rotation2ydir(rotation)

    module_logger.debug('Check schema versions...')


    # STEP Header https://en.wikipedia.org/wiki/ISO_10303-21
    header = f"""ISO-10303-21;
HEADER;
FILE_DESCRIPTION(
/* description */ ('{file_descr}'),
/* implementation_level */ '2;1');
FILE_NAME(
/* name */ '{filename}',
/* time_stamp */ '{t_string}',
/* author */ ('{author}'),
/* organization */ ('{organization}'),
/* preprocessor_version */ '{app}',
/* originating_system */ '{app}',
/* authorization */ '{authorization}');
FILE_SCHEMA(('{schema_version}'));
ENDSEC;
"""
#FILE_SCHEMA(('{ifcopenshell.schema_identifier}'));  --> this line shoud be in 100


    # ifc_utils header
    if minimal:  # Only IfcProject TODO -> IfcProjectLibrary? For Classification / Pset / Types
        data = f"""
DATA;
#1 = IFCPROJECT('{project_globalid}',$,'{project_name}','{project_descr}',$,$,$,$,$);
ENDSEC;
END-ISO-10303-21;
"""

    else:  # Full Template with Representation Context, Units, etc.
        data = f"""
DATA;
#1=IFCPERSON($,$,'{author}',$,$,$,$,$);
#2=IFCORGANIZATION($,'{organization}',$,$,$);
#3=IFCPERSONANDORGANIZATION(#1,#2,$);
#4=IFCAPPLICATION(#2,'{app_version}','{app}','');
#5=IFCOWNERHISTORY(#3,#4,$,.ADDED.,{int(t)},#3,#4,{int(t)});
#6=IFCDIMENSIONALEXPONENTS(0,0,0,0,0,0,0);
#7=IFCSIUNIT(*,.LENGTHUNIT.,$,.METRE.);
#8=IFCSIUNIT(*,.AREAUNIT.,$,.SQUARE_METRE.);
#9=IFCSIUNIT(*,.VOLUMEUNIT.,$,.CUBIC_METRE.);
#10=IFCSIUNIT(*,.PLANEANGLEUNIT.,$,.RADIAN.);
#11=IFCMEASUREWITHUNIT(IFCPLANEANGLEMEASURE(0.017453292519943295),#10);
#12=IFCCONVERSIONBASEDUNIT(#6,.PLANEANGLEUNIT.,'DEGREE',#11);
#13=IFCUNITASSIGNMENT((#7,#8,#9,#12));
"""
        if georeferencing == 'IfcGeometricRepresentationContext':
            # Georeferencing (Translation, Rotation with the placement of the GeometricRepresentationContext
            # If no georef needed, x_direction=X, origin=O, true_north=Z[0:2]
            module_logger.debug('Georeferencing with IfcGeometricRepresentationContext')
            georef = f"""
#14=IFCDIRECTION({x_direction});
#15=IFCDIRECTION((0.,0.,1.));
#16=IFCCARTESIANPOINT({origin});
#17=IFCAXIS2PLACEMENT3D(#16,#15,#14);
#18=IFCDIRECTION({true_north[0:2]});
#19=IFCGEOMETRICREPRESENTATIONCONTEXT($,'Model',3,1.E-05,#17,#18);
#20=IFCPROJECT('{project_globalid}',#5,'{project_name}','{project_descr}',$,$,$,(#19),#13);
ENDSEC;
END-ISO-10303-21;
"""
        elif georeferencing == 'IfcMapConversion':
            module_logger.debug('Georeferencing with IfcMapConversion')
            georef = f"""
#14=IFCDIRECTION({X});
#15=IFCDIRECTION((0.,0.,1.));
#16=IFCCARTESIANPOINT({O});
#17=IFCAXIS2PLACEMENT3D(#16,#15,#14);
#18=IFCDIRECTION({true_north[0:2]});
#19=IFCGEOMETRICREPRESENTATIONCONTEXT($,'Model',3,1.E-05,#17,#18);
#20=IFCPROJECT('{project_globalid}',#5,'{project_name}','{project_descr}',$,$,$,(#19),#13);
#21=IFCPROJECTEDCRS('EPSG:2056','CH1903+ / LV95 -- Swiss CH1903+ / LV95','CH1903+','LN02','CH1903+ / LV95',$,#7);
#22=IFCMAPCONVERSION(#19,#21,{origin[0]},{origin[1]},{origin[2]}, {x_direction[0]},{x_direction[1]},$);
ENDSEC;
END-ISO-10303-21;
"""
        else:
            module_logger.debug('No georeferencing in IfcGeometricRepresentationContext or IfcMapConversion')
            georef = f"""
#14=IFCDIRECTION({X});
#15=IFCDIRECTION((0.,0.,1.));
#16=IFCCARTESIANPOINT({O});
#17=IFCAXIS2PLACEMENT3D(#16,#15,#14);
#18=IFCDIRECTION({Y[0:2]});
#19=IFCGEOMETRICREPRESENTATIONCONTEXT($,'Model',3,1.E-05,#17,#18);
#20=IFCPROJECT('{project_globalid}',#5,'{project_name}','{project_descr}',$,$,$,(#19),#13);
ENDSEC;
END-ISO-10303-21;
"""

        data = data + georef

    template = header + data

    return template

def to_wgs84(easting, northing, epsg = 32632):
    # Create UTM CRS and WGS84 CRS objects
    utm_crs = CRS.from_epsg(epsg)
    wgs84_crs = CRS.from_epsg(4326)

    # Create the transformer object for conversion
    transformer = Transformer.from_crs(utm_crs, wgs84_crs)

    # Convert UTM coordinates to WGS84 latitude and longitude
    lon, lat = transformer.transform(easting, northing)

    # Convert decimal degrees to degrees, minutes, seconds, and millionth 
    # seconds
    def decimal_degrees_to_dmsm(decimal_degrees):
        degrees = int(decimal_degrees)
        decimal_minutes = (decimal_degrees - degrees) * 60
        minutes = int(decimal_minutes)
        decimal_seconds = (decimal_minutes - minutes) * 60
        seconds = int(decimal_seconds)
        millionth_seconds = int((decimal_seconds - seconds) * 1000000)
        return degrees, minutes, seconds, millionth_seconds

    # Convert latitude and longitude to desired format
    dmsm_latitude = decimal_degrees_to_dmsm(lat)
    dmsm_longitude = decimal_degrees_to_dmsm(lon)

    return dmsm_latitude, dmsm_longitude


def create_ifcsite_with_survey_point(ifcfile, project, owner_history, context, 
                                     projectorigin,epsg = 32632, rotation=0, 
                                     georef='', descr='Grundstueck'):
    """
    IfcSite,
    ShapeRepresentation as the Product Base Point (0/0/0)
    :param ifcfile: Instance of ifcopenshell
    :param owner_history: IfcOwnerHistory
    :param project: IfcProject
    :param context: IfcGeometricRepresentationContext
    :param origin: tuple, coordinate of project origin for georeferencing
    :param rotation: project rotation to world coordinate system for georeferencing
    :param georef: georeferencing method
                if georef == 'IfcSite
                -> Georeferencing with the placement of IfcSite
    :return: IfcSite, IfcLocalPlacement
    """
    lon, lat = to_wgs84(projectorigin[0], projectorigin[1])
    elevation = projectorigin[2]

    # Survey Point at the Project Origin
    survey_point = ifcfile.createIfcCartesianPoint((0., 0., 0.))
    # GeometricCurveSet
    geom_curve_set = ifcfile.createIfcGeometricCurveSet([survey_point])
    site_representation = ifcfile.createIfcShapeRepresentation(
        context, 'ProductBasePoint', 'GeometricCurveSet', [geom_curve_set])
    site_product_shape = ifcfile.createIfcProductDefinitionShape(
        'Georeferenzierung', None, [site_representation])
    
    # Site placement, all other entities should be relative to the site placement
    if georef == 'IfcSite':  # Do georeferencing with IfcSite
        module_logger.debug('Georeferencing with Placement of IfcSite')
        site_placement = create_ifclocalplacement(ifcfile, projectorigin, 
                                                  dirZ=(0., 0., 1.), 
                                                  dirX=rotation2xdir(rotation))
    else:  # No Georef or Done with IfcMapConversion OR IfcGeometricRepresentationContext
        site_placement = create_ifclocalplacement(ifcfile)
    
    # create IfcSite
    site = ifcfile.createIfcSite(create_guid(), owner_history, descr, 
                                 'ProjectBasePoint = SurveyPoint', None, 
                                 site_placement, site_product_shape, None, 
                                 'ELEMENT', lat, lon, elevation, None,
                                 None)
    container_project = ifcfile.createIfcRelAggregates(
        create_guid(), owner_history, "Project Container", None, project, 
        [site])
    return site, site_placement


def colorFromXlsx(val, rKey, gKey, bKey):
    def readOne(key):
        try:
            n = (abs(int(val[key])) % 256) / 256
        except:
            n = 0
        return n
    return [readOne(rKey),readOne(gKey),readOne(bKey)]

def read_STL(stl_file):
    # Conversion to nested list of triangles:
    # [triangle [vertex [x,y,z]],[vertex [x,y,z]],[vertex [x,y,z]]]
    triangles = []
    with open(stl_file) as f:
        while True:
            line = f.readline().strip()
            if not line:
                break
            if line == "outer loop":
                triangle = [f.readline().strip() for i in range(3)]
                vertices = []
                for i in range(3):
                    vertex = triangle[i].rstrip().split(" ")[1:]
                    vertex = [float(i) for i in vertex]
                    vertices.append(vertex)
                triangles.extend([vertices])

    # Create CartesianPointList Tuple
    vertices = []
    for i in range(len(triangles)):
        for j in range(3):
            vertex = tuple(triangles[i][j])
            vertices.append(vertex)

    # Create IfcTriangulatedFaceSet Tuple
    faces = []
    for i in range(len(triangles)):
        face = (i * 3 + 1, i * 3 + 2, i * 3 + 3)
        faces.append(face)

    return vertices, faces

def make_plane(ifcfile, layer_name, vertices, faces, col_RGB, properties,
               site_placement, context, owner_history, schema_version,
               building, site):
    """
    Create a 3D Plane in IFC format with associated properties and geometry.

    :param ifcfile: Instance of ifcopenshell, representing the IFC file.
    :param layer_name: Name of the layer to assign to the created IFC object.
    :param vertices: List of tuples representing the coordinates of the plane's
        vertices.
    :param faces: List of tuples representing the face indices for creating the
        plane geometry.
    :param col_RGB: Color representation in RGB format, either as a list of 
        three values or a predefined IFC color object.
    :param properties: Dictionary of key-value pairs representing custom 
        properties to be added to the IFC object.
    :param site_placement: IfcLocalPlacement, defining the placement of the 
        object relative to the IfcSite.
    :param context: IfcGeometricRepresentationContext, representing the context
        for geometric representations.
    :param owner_history: IfcOwnerHistory, tracking the ownership and 
        modification history of the IFC object.
    :param schema_version: String defining the IFC schema version, e.g., 
        'IFC2X3' or 'IFC4'.
    :param building: IfcBuilding, specifying the building structure to which 
        the plane object belongs.
    :param site: IfcSite, specifying the site structure for spatial 
        organization.
    :return: IfcBuildingElementProxy, representing the proxy object created for
        the plane with geometric and property data.
    """


    proxy_placement = create_ifclocalplacement(
        ifcfile, relative_to=site_placement
        )
    # writeGeometryToProxy(
    #     val_df.iloc[0], ifcfile, context, owner_history, placement, site, 
    #     vertices, faces, file_name, val_df.iloc[0], "val", col_RGB
    # )
    if schema_version == "IFC2X3":
        cords = [ifcfile.createIfcCartesianPoint(p) for p in vertices]
        outer_bound_list = [ifcfile.createIfcFaceOuterBound(
            ifcfile.createIfcPolyLoop(
                [cords[i-1] for i in face]), True) for face in faces]
        
        face_list = [ifcfile.createIfcFace(
            [outer_bound]) for outer_bound in outer_bound_list]
        
        connected_face_set = ifcfile.createIfcConnectedFaceSet(
            [face for face in face_list])
        face_based_surface_model = ifcfile.createIfcFaceBasedSurfaceModel(
            [connected_face_set])

        representation_item = ifcfile.createIfcShapeRepresentation(
            context,"Body", "SweptSolid", [face_based_surface_model])

        product_shape = ifcfile.createIfcProductDefinitionShape(
            None, None, [representation_item]
            )
        proxy = ifcfile.createIfcBuildingElementProxy(
            create_guid(), owner_history, str(layer_name), None, None, 
            proxy_placement, product_shape, None
            )
        
   
        if type(col_RGB) == list:
            color = ifcfile.createIfcColourRgb(
                None,col_RGB[0]/256,col_RGB[1]/256,col_RGB[2]/256)
            styleShading = ifcfile.createIfcSurfaceStyleShading(color)
            surfaceStyle = ifcfile.createIfcSurfaceStyle(
                None,"BOTH",[styleShading])
            presentation_style_assignment = ifcfile.createIfcPresentationStyleAssignment(
                [surfaceStyle])
            styledItem = ifcfile.createIfcStyledItem(
                connected_face_set,[presentation_style_assignment],None)
        else:
            styled_item_col = ifcfile.createIfcStyledItem(
                connected_face_set, [col_RGB], None)


    else:
        cords = ifcfile.createIfcCartesianPointList3d(vertices)
    
        faces = ifcfile.createIfcTriangulatedFaceSet(
            cords, None, None, faces, None)
    
        # Colour geometries
        if type(col_RGB) == list:
            color = ifcfile.createIfcColourRgb(
                None,col_RGB[0]/256,col_RGB[1]/256,col_RGB[2]/256)
            styleShading = ifcfile.createIfcSurfaceStyleShading(color, None)
            surfaceStyle = ifcfile.createIfcSurfaceStyle(
                None,"BOTH",(styleShading,))
            styledItem = ifcfile.createIfcStyledItem(
                faces,(surfaceStyle,),None)
        else:
            styled_item_col = ifcfile.createIfcStyledItem(
                faces, [col_RGB], None)
        bodyRepresentation = ifcfile.createIfcShapeRepresentation(
            context, "Body", "SweptSolid", [faces]
            )
        product_shape = ifcfile.createIfcProductDefinitionShape(
            None, None, [bodyRepresentation]
            )
        proxy = ifcfile.createIfcBuildingElementProxy(
            create_guid(), owner_history, str(layer_name), None, None, 
            proxy_placement, product_shape, None
            )    

    if properties != {}:   
        propertySingleValues = []
        for property_label, property_value in properties.items():
            valueStr = '{}'.format(property_value)
            propertySingleValueWriter = ifcfile.createIfcPropertySingleValue(
                property_label, property_label,
                ifcfile.create_entity('IfcText', valueStr),
                None
                )
            propertySingleValues.append(propertySingleValueWriter)
        propertySet = ifcfile.createIfcPropertySet(
            create_guid(), owner_history, '{}'.format("Eigenschaften"),
            None, propertySingleValues
            )
        RelDefinesByProperties = ifcfile.createIfcRelDefinesByProperties(
            create_guid(), owner_history, None, None, [proxy], propertySet
            )
    return proxy
    

def make_IfcBoreholes(x,y,z,nbor,IDS,FROM_TO,IDF,soiltype,lowerdep,
                      verticaloption, ifcfile,site_placement,owner_history, 
                      circle_profile_template, axis_placement_template, 
                      direction,subcontext,SURVEY, site, building):
    """
    Create a 3D representation of boreholes with associated geospatial and 
    property data in IFC format.

    :param x: List of X coordinates for the borehole locations.
    :param y: List of Y coordinates for the borehole locations.
    :param z: List of Z coordinates for the borehole locations.
    :param nbor: Number of boreholes (overwritten as the length of the x, y, z 
                                      lists).
    :param IDS: List of unique borehole identifiers.
    :param FROM_TO: DataFrame or array containing borehole layer data.
    :param IDF: Array indicating borehole IDs corresponding to layers in 
        FROM_TO.
    :param soiltype: Array representing the type of soil encountered at each 
        layer of the boreholes.
    :param lowerdep: Array of lower depth values for each borehole layer.
    :param verticaloption: Boolean flag to indicate if the boreholes are 
        vertical (True) or angled (False).
    :param ifcfile: Instance of ifcopenshell, representing the IFC file.
    :param site_placement: IfcLocalPlacement, defining the placement of 
        boreholes relative to the site.
    :param owner_history: IfcOwnerHistory, tracking the ownership and 
        modification history of the IFC object.
    :param circle_profile_template: IfcCircleProfileDef, defining the 
        cross-sectional profile of the borehole.
    :param axis_placement_template: IfcAxis2Placement3D, defining the 
        orientation of the borehole.
    :param direction: IfcDirection, specifying the extrusion direction for the 
        borehole geometry.
    :param subcontext: IfcGeometricRepresentationSubContext, representing the 
        geometric context for borehole representation.
    :param SURVEY: DataFrame containing general survey properties for the 
        boreholes.
    :param site: IfcSite, specifying the site structure for spatial 
        organization.
    :param building: IfcBuilding, specifying the building structure associated 
        with the boreholes.
    :return: None, modifies the IFC file by adding borehole elements and their 
        associated properties.
    """

    building_story_list = []
   
    for i_bhl in tqdm(np.arange(0,nbor,1)):
        H_Id = IDS[i_bhl] 
        if sum(IDF==H_Id) == 0: 
            print(f"No layer descriptions available for: {H_Id}")
            continue

        borFT_pset = FROM_TO[IDF==H_Id]
        borFT_soiltype = np.array(soiltype[IDF==H_Id])
        borFT_lowerdep = np.array(lowerdep[IDF==H_Id])
        
        DIN_arr = []
        HBA = []
        for DIN in borFT_soiltype:
            #print(ba)
            DIN_temp = DIN.split(sep=",")
            HBA_temp = 10
            if 'G' in DIN_temp[0]:
                HBA_temp = 3
                r,g,b = 255, 255, 0
            if 'A' in DIN_temp[0]:
                HBA_temp = 10
                r,g,b = 255, 255, 255
            if 'S' in DIN_temp[0]: 
                HBA_temp = 2
                r,g,b = 255, 128, 0
            if 'U' in DIN_temp[0]: 
                HBA_temp = 1
                r,g,b = 102, 102, 0
            if 'T' in DIN_temp[0]: 
                HBA_temp = 1
                r,g,b = 204, 0, 204        
            if 'Z' in DIN_temp[0]: 
                HBA_temp = 4
                r,g,b = 0, 0, 255
            HBA.append(HBA_temp)
            DIN_arr.append(DIN_temp[0])

        xi = x[i_bhl]
        yi = y[i_bhl]
        zi = z[i_bhl]
    
        drill_location = ifcfile.createIfcCartesianPoint(
            (float(xi), float(yi), float(zi)))
        axis_placement_3d_2 = ifcfile.createIfcAxis2Placement3D(
            drill_location, None, None)
        drill_placement = ifcfile.createIfcLocalPlacement(
            site_placement, axis_placement_3d_2)
        building_storey = ifcfile.createIfcBuildingStorey(
            create_guid(), owner_history, H_Id, None, None, drill_placement, 
            None, None, 'COMPLEX', None
            )
        building_story_list.append(building_storey)
        
        property_values = [] 
        for single_property in SURVEY.columns:
            property_values.append(
                ifcfile.createIfcPropertySingleValue(
                    "{}".format(single_property), 
                    "{}".format(single_property), 
                    ifcfile.create_entity("IfcText", "{}".format(
                        SURVEY[single_property].iloc[i_bhl])), 
                    None
                    )
                )

        property_set = ifcfile.createIfcPropertySet(
            create_guid(), owner_history, "pset_bhl_common", None, 
            property_values)
        ifcfile.createIfcRelDefinesByProperties(
            create_guid(), owner_history, None, None, [building_storey], 
            property_set) 
    
        segment_len = np.insert(
            np.diff(borFT_lowerdep), 0, borFT_lowerdep[0] - 0)

        IfcSpaces = []
        for i_bh_segment in range(len(borFT_lowerdep)):
            HBA[i_bh_segment]
            # Get the length of the line segment for the extrusion
            extrusion = segment_len[i_bh_segment]  
            
            r,g,b = 255, 255, 255
            if 'G' in DIN_arr[i_bh_segment]: r,g,b = 255, 255, 0
            if 'A' in DIN_arr[i_bh_segment]: r,g,b = 255, 255, 255
            if 'S' in DIN_arr[i_bh_segment]: r,g,b = 255, 128, 0
            if 'U' in DIN_arr[i_bh_segment]: r,g,b = 102, 102, 0
            if 'T' in DIN_arr[i_bh_segment]: r,g,b = 204, 0, 204        
            if 'Z' in DIN_arr[i_bh_segment]: r,g,b = 0, 0, 255
            cartesian_point_col = ifcfile.createIfcCartesianPoint(
                (0., 0., float(-borFT_lowerdep[i_bh_segment])))
            axis_placement_col = ifcfile.createIfcAxis2Placement3D(
                cartesian_point_col, None, None)
            local_placement_col = ifcfile.createIfcLocalPlacement(
                drill_placement, axis_placement_col)
            
            
            extruded_area_solid = ifcfile.createIfcExtrudedAreaSolid(
                circle_profile_template, axis_placement_template, direction, 
                float(extrusion))
            #################################################################################
            shape_representation = ifcfile.createIfcShapeRepresentation(
                subcontext, 'Body', 'SweptSolid', [extruded_area_solid])
            product_definition_shape = ifcfile.createIfcProductDefinitionShape(
                None, None, [shape_representation])
            
            colour_rgb_col = ifcfile.createIfcColourRgb(
                None, r/256, g/256, b/256)
            surface_style_shading_col = ifcfile.createIfcSurfaceStyleShading(
                colour_rgb_col)
            surface_style_col = ifcfile.createIfcSurfaceStyle(
                None,'BOTH', [surface_style_shading_col])
            presentation_style_assignment_col = ifcfile.createIfcPresentationStyleAssignment(
                (surface_style_col,))
            styled_item_col = ifcfile.createIfcStyledItem(
                extruded_area_solid, (presentation_style_assignment_col,), 
                None)
            building_element_proxy = ifcfile.createIfcBuildingElementProxy(
                create_guid(), owner_history, str(i_bh_segment), None, None, 
                local_placement_col, 
                product_definition_shape, None, None)
            space = building_element_proxy
    
            property_values = [] 
            for single_property in borFT_pset.columns:
                property_values.append(
                    ifcfile.createIfcPropertySingleValue(
                        "{}".format(single_property), 
                        "{}".format(single_property), 
                        ifcfile.create_entity("IfcText", "{}".format(
                            borFT_pset[single_property].iloc[i_bh_segment])), 
                        None
                        )
                    )
    
            property_set = ifcfile.createIfcPropertySet(
                create_guid(), owner_history, "pset_bhl_common", None, 
                property_values)
            ifcfile.createIfcRelDefinesByProperties(
                create_guid(), owner_history, None, None, [space], 
                property_set) 
        
            
            IfcSpaces.append(space)
    
        container_space = ifcfile.createIfcRelAggregates(
            create_guid(), owner_history, None, None, building_storey, 
            IfcSpaces)
    
    ifcfile.createIfcRelAggregates(create_guid(), owner_history, "Bohrungen", 
                                   None, building, building_story_list)
    ifcfile.createIfcRelAggregates(create_guid(), owner_history, "Bohrungen", 
                                   None, site, [building])