""" Functions used for the pre-processing """

import os
from geopy.distance import vincenty
from boltons.iterutils import pairwise
import geopandas as gpd
import numpy as np
from time import sleep
import pandas as pd

import matplotlib as mpl
import matplotlib.pyplot as plt


def spatial_overlays(df1, df2, how='intersection'):
    """
    Compute overlay intersection of two GeoPandasDataFrames df1 and df2. This
    function is much faster compared to the original geopandas overlay method.
    
    Adapted from https://github.com/geopandas/geopandas/issues/400 
    
    Args:
        *df1* : The dataframe to be either intersected or 'erased', using the difference function.
            
        *df2* : The dataframe to intersect or to erase, using the difference function.
        
    Return:
        *df1*: an either intersected or (partly) erased geopandas dataframe.
    """

    df1 = df1.copy()
    df2 = df2.copy()
    df1['geometry'] = df1.geometry.buffer(0)
    df2['geometry'] = df2.geometry.buffer(0)

# =============================================================================
#     """ the part of the function which does the intersection analysis """
# =============================================================================
    if how=='intersection':
        # Spatial Index to create intersections
        spatial_index = df2.sindex
        df1['bbox'] = df1.geometry.apply(lambda x: x.bounds)
        df1['histreg']=df1.bbox.apply(lambda x:list(spatial_index.intersection(x)))
        pairs = df1['histreg'].to_dict()
        nei = []
        for i,j in pairs.items():
            for k in j:
                nei.append([i,k])
        
        pairs = gpd.GeoDataFrame(nei, columns=['idx1','idx2'], crs=df1.crs)
        pairs = pairs.merge(df1, left_on='idx1', right_index=True)
        pairs = pairs.merge(df2, left_on='idx2', right_index=True, suffixes=['_1','_2'])
        pairs['Intersection'] = pairs.apply(lambda x: (x['geometry_1'].intersection(x['geometry_2'])).buffer(0), axis=1)
        pairs = gpd.GeoDataFrame(pairs, columns=pairs.columns, crs=df1.crs)
        cols = pairs.columns.tolist()
        cols.remove('geometry_1')
        cols.remove('geometry_2')
        cols.remove('histreg')
        cols.remove('bbox')
        cols.remove('Intersection')
        dfinter = pairs[cols+['Intersection']].copy()
        dfinter.rename(columns={'Intersection':'geometry'}, inplace=True)
        dfinter = gpd.GeoDataFrame(dfinter, columns=dfinter.columns, crs=pairs.crs)
        dfinter = dfinter.loc[dfinter.geometry.is_empty==False]
        return dfinter
# =============================================================================
#     """ the part of the function which does the difference/erase analysis """
# =============================================================================
    elif how=='difference':
        spatial_index = df2.sindex
        df1['bbox'] = df1.geometry.apply(lambda x: x.bounds)
        df1['histreg']=df1.bbox.apply(lambda x:list(spatial_index.intersection(x)))
        df1['new_g'] = df1.apply(lambda x: reduce(lambda x, y: x.difference(y).buffer(0), [x.geometry]+list(df2.iloc[x.histreg].geometry)) , axis=1)
        df1.geometry = df1.new_g
        df1 = df1.loc[df1.geometry.is_empty==False].copy()
        df1.drop(['bbox', 'histreg', 'new_g'], axis=1, inplace=True)
        return df1


def get_country(country,continent_osm,base_path,overwrite=False,RAI=False):  

    """
    Extraction of the road data for the specified country from OSM. We use the continental OSM file, downloaded from http://download.geofabrik.de.
    
    Args:
        *country* : The country for which we calculate the RAI.
        
        *continent_osm* : The continent the country 'belongs' to. This is required for the osm extraction.
                        
        *base_path* : The base path to location of all files and scripts.
        
        *overwrite* : This is set on True by default, set on False if you are sure the input files are correct (it will save a lot of computation time).
        *RAI* : This is set on False by default. set on True if this country extracting is used for the RAI analysis. It will skip the road length calculation (saving computation time).
        
    Returns:
        *load_country* : A geodataframe containing all the roads and their length.
    """    

# =============================================================================
#     """ First set all paths for output dirs"""
# =============================================================================
    ctry_out = os.path.join(base_path,'country_data')
    osm_path_out = os.path.join(base_path,'osm_country')
    poly_dir = os.path.join(base_path,'poly_files')

    
# =============================================================================
#     """ Set the paths for the files we are going to use and write"""
# =============================================================================
    country_poly = os.path.join(poly_dir,country+'.poly') 
    country_shp =  os.path.join(ctry_out,country+'.shp') 
    country_pbf = os.path.join(osm_path_out,country+'.osm.pbf') 
    
# =============================================================================
#     # extract osm file for the country and write to shapefile
# =============================================================================

    if os.path.exists(country_pbf) is not True:
        clip_osm(continent_osm,country_poly,country_pbf)
        
    # get shapefile as output
#    if   (os.path.getsize(country_shp) == 143):
    if (os.path.exists(country_shp) is not True) or overwrite == True:
        extract_osm(country_shp,country_pbf)

# =============================================================================
#     Load the shapefile, remove road tags which only occur less than 15 times 
#     and estimate the length of the remaining roads    
# =============================================================================
    try:
        load_country = gpd.read_file(country_shp)
    except:
        sleep(30)
        try:
            load_country = gpd.read_file(country_shp)
        except:
            sleep(60)
            load_country = gpd.read_file(country_shp)

    uniq = load_country['highway'].value_counts()
    uniq = list(uniq[uniq > 20].index)
    uniq.extend(['primary','secondary','trunk','motorway'])
    uniq = list(set(uniq))
    load_country = load_country[load_country['highway'].isin(uniq)]
    
    if RAI is False:
        load_country['distance'] = load_country['geometry'].apply(line_length)
        load_country = load_country[load_country['distance'] < 500]        
        

# =============================================================================
#   Add a new column to the dataframe with the aggregated road classification
# =============================================================================
    load_country = map_roads(load_country)
    
# =============================================================================
#     And return the geodataframe
# =============================================================================
    return load_country



def clip_osm(continent_osm,country_poly,country_pbf):
    """ Clip the country osm file from the larger continent (or planet) file and save to a new osm.pbf file. 
    This is much faster compared to clipping the osm.pbf file while extracting through ogr2ogr.
    
    This function uses the osmconvert tool, which can be found at http://wiki.openstreetmap.org/wiki/Osmconvert. 
    
    Either add the directory where this executable is located to your environmental variables or just put it in the **scripts** directory.
    
    Args:
        *continent_osm* : path string to the osm.pbf file of the continent associated with the country.
        
        *country_poly* : path string to the .poly file, made through the 'create_poly_files' function.
        
        *country_pbf* : path string indicating the final output dir and output name of the new .osm.pbf file.
        
    Returns:
        A clipped .osm.pbf file.
    """ 

    os.system('osmconvert64 %s -B=%s --complete-ways -o=%s' %(continent_osm,country_poly,country_pbf))


def extract_osm(country_shp,country_pbf):
    """Extract a shapefile with all the road information from the openstreetmap file.
    
    Args:
        *country_shp* : The path string indicating the final output directory and output name of the new shapefile.
        
        *country_pbf* : The path string indicating the directory and name of the .osm.pbf file.
        
    Returns:
        A shapefile with all the roads of the clipped country. The shapefile will be in *WGS84* (epsg:4326). This is the same coordinate system as Openstreetmap.
    """
    
    os.system('ogr2ogr -overwrite -skipfailures -f "ESRI Shapefile" -progress \
              -sql "SELECT highway FROM lines WHERE highway IS \
              NOT NULL" -lco ENCODING=UTF-8 '+country_shp+" " + country_pbf)    
        


def line_length(line, ellipsoid='WGS-84'):
    """Length of a line in meters, given in geographic coordinates.
    Adapted from https://gis.stackexchange.com/questions/4022/looking-for-a-pythonic-way-to-calculate-the-length-of-a-wkt-linestring#answer-115285
    Args:
        *line* : A shapely LineString object with WGS-84 coordinates.
        
        *ellipsoid* : The string name of an ellipsoid that `geopy` understands (see http://geopy.readthedocs.io/en/latest/#module-geopy.distance).
    Returns:
        The length of the line in meters.
    """
    if line.geometryType() == 'MultiLineString':
        return sum(line_length(segment) for segment in line)

    try:
        return sum(
            vincenty(a, b, ellipsoid=ellipsoid).kilometers
            for a, b in pairwise(line.coords)
        )
    except:
        return sum(
            vincenty(a, b, ellipsoid=ellipsoid).kilometers
            for a, b in pairwise(list([t[::-1] for t in list(line.coords)]))
        )



def map_roads(load_country):

    """ 
    To create a new column with an aggregated list of road types. 
    
    Args:
        *load_country* : A geodataframe containing all the roads of a country.
        
    Returns:
        *load_country* : The same geodataframe but with an additional 'roads' column containing the aggregated road types.
    """

    dict_map = {
"disused" : "other",
"dummy" : "other",
"planned" : "other",
"platform" : "other",
"unsurfaced" : "track",
"traffic_island" : "other",
"razed" : "other",
"abandoned" : "other",
"services" : "track",
"proposed" : "other",
"corridor" : "track",
"bus_guideway" : "other",
"bus_stop" : "other",
"rest_area" : "other",
"yes" : "other",
"trail" : "other",
"escape" : "track",
"raceway" : "other",
"emergency_access_point" : "track",
"emergency_bay" : "track",
"construction" : "track",
"bridleway" : "track",
"cycleway" : "other",
"footway" : "other",
"living_street" : "tertiary",
"path" : "track",
"pedestrian" : "other",
"primary" : "primary",
"primary_link" : "primary",
"residential" : "tertiary",
"road" : "secondary",
"secondary" : "secondary",
"secondary_link" : "secondary",
"service" : "tertiary",
"steps" : "other",
"tertiary" : "tertiary",
"tertiary_link" : "tertiary",
"track" : "track",
"unclassified" : "tertiary",
"trunk" : "primary",
"motorway" : "primary",
"trunk_link" : "primary",
"motorway_link" : "primary",
"via_ferrata": "other",
"elevator": "other",
"crossing": "other",
"seasonal": "other",
"traffic_signals":"other",
"piste":"other",
"dismantled": "other",
"winter_road":"other",
"access":"other",
"ohm:military:Trench":"other",
"no":"other",
"byway":"other",
"unmarked_route":"other",
"track_grade1":"track",
"track_grade2":"track",
"track_grade3":"track",
"track_grade4":"track",
"track_grade5":"track",
"unknown":"other"
}
    
    load_country['roads'] = load_country['fclass'].map(lambda x: (dict_map[x])) 
    
    return load_country


def explode(gdf):
    """ 
    Explodes a geodataframe 
    
    Will explode muti-part geometries into single geometries. Original index is
    stored in column level_0 and zero-based count of geometries per multi-
    geometry is stored in level_1
    
    Args:
        gdf (gpd.GeoDataFrame) : input geodataframe with multi-geometries
        
    Returns:
        gdf (gpd.GeoDataFrame) : exploded geodataframe with a new index 
                                 and two new columns: level_0 and level_1
        
    """
    gs = gdf.explode()
    gdf2 = gs.reset_index().rename(columns={0: 'geometry'})
    gdf_out = gdf2.merge(gdf.drop('geometry', axis=1), left_on='level_0', right_index=True)
    gdf_out = gdf_out.set_index(['level_0', 'level_1']).set_geometry('geometry')
    gdf_out.crs = gdf.crs
    return gdf_out

def extract_osm_rail(country_shp,country_pbf):
    """Extract a shapefile with all the road information from the openstreetmap file.
    
    Args:
        *country_shp* : The path string indicating the final output directory and output name of the new shapefile.
        
        *country_pbf* : The path string indicating the directory and name of the .osm.pbf file.
        
    Returns:
        A shapefile with all the roads of the clipped country. The shapefile will be in *WGS84* (epsg:4326). This is the same coordinate system as Openstreetmap.
    """
    
    os.system('ogr2ogr -overwrite -skipfailures -f "ESRI Shapefile" -progress \
              -sql "SELECT railway FROM lines WHERE railway IS \
              NOT NULL" -lco ENCODING=UTF-8 '+country_shp+" " + country_pbf)    


#%% Method to Create buffer
def Create_Buffer(gdf,epsg1,size):
    """ Create buffer of size "size" from GeoDataFrame in the projection of epsg1
    Args:
        *gdf* : GeoDataFrame
        *epsg1* : projection number
        *size* : size of buffer
    Returns:
        A GeoDataFrame with new buffered geometries
    """

    gdf1 = gdf.copy()
    gdf1 = gdf1.to_crs(epsg=epsg1)
    gdf1['geometry'] = gdf1.buffer(size)
    gdf1 = gdf1.to_crs(epsg=4326)
    return gdf1


def geom_within_country(x,geo_country):
    """ Given a geometry, return the geometry of its intersection with geo_country
    Args:
        *x* : Geopandas row
        *geo_country* : Shapely geometry
    Returns:
        A shapely geometry of the intersection
    """

    if not x['within_country']:
        return x['geometry'].intersection(geo_country)
    else:
        return x['geometry']

    
def delete_roads_urb(x,geo_urban):
    """ Given a geometry, returns the geometry of its difference with geo_country
    Args:
        *x* : Geopandas row
        *geo_urban* : Shapely geometry
    Returns:
        A shapely geometry of the difference
    """
    if x['inter_urb']:
        return x['geometry'].difference(geo_urban)
    else:
        return x['geometry']