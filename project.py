import os

# PostGIS sets GDAL_DATA and PROJ_LIB on the system environments to its own folders, and due to that fiona stopped working properly.
# The next two lines will remove those variables from the environment dictionary. 
# This is not documented anywhere, it's my own solution
os.environ.pop('GDAL_DATA', None)
os.environ.pop('PROJ_LIB', None)

import shutil
import fiona
from fiona.crs import from_epsg
from pyproj import Transformer, CRS
import urllib.request
import geojson
from geojson import FeatureCollection
from shapely.geometry import shape, mapping
from shapely.ops import unary_union
import itertools
from operator import itemgetter

import psycopg2
import psycopg2.errors

from osgeo import gdal
from osgeo import ogr
from osgeo import osr

from rasterstats import zonal_stats

# General workflow

# Create output and temporary folders
# Reproject CAOP into WGS 84, which is the coordinate system of ipma's data
# Dissolve CAOP into counties
# Obtain IPMA's data
# Save to the database
# Remove stations data outside the boundary of Portugal Continental
# Interpolate to get data on counties with no data
# Create a raster of the interpolation and clip it to the boundary of Portugal Continental
# Create a choropleth map

def create_folders(paths):
    # Always recreate the folders, even if it existed previously
    for path in paths:
        if os.path.exists(path):
            shutil.rmtree(path)
        os.mkdir(path)
        print("\t Folder created: {0}".format(path))
    
def delete_folders(paths):
    for path in paths:
        if os.path.exists(path):
            shutil.rmtree(path)

def reproject(infile, outfile, epsg):
    # Only reproject if the reprojected file does not exist.
    # This is an exception because it's still an input data, not an output and will probably never change.
    # As such, speeds up the processing time.
    if os.path.exists(outfile):
        print ("\t Reprojection not needed. Skipping step.")
        return

    # The reprojection works by creating a new output file with the new coordinate system.
    # Then every coordinate needs to be transformed from the input coordinate system to the output coordinate system.
    try:
        print("\t Reprojecting {0} to EPSG {1}...".format(infile, epsg))
        with fiona.open(infile) as input:
            output_schema = input.schema.copy()
            with fiona.open(outfile, 'w', input.driver, output_schema, from_epsg(epsg)) as output:
                transformer = Transformer.from_crs(input.crs_wkt, output.crs_wkt)
                for feature in input:
                    new_coords = []
                    for coordinate in feature['geometry']['coordinates']:
                        x2, y2 = transformer.transform(*zip(*coordinate)) # Swap the coordinates. EPSG 4326 uses lon / lat, not lat / lon
                        new_coords.append(zip(y2, x2))
                    feature['geometry']['coordinates'] = new_coords
                    output.write(feature)
    except:
        print("\t Reprojection failed, cleanup any created files.")

# Implementation of the dissolve inspired by the dissolve in ArcGIS
# If no fields are specified, it will dissolve every geometry. Useful to have a country boundary.
# Else it will dissolve by the parametres specified
# This solution was inspired by an answer on the GIS StackExchange: https://gis.stackexchange.com/a/274516
# It will group every item in 
def dissolve(infile, outfile, fields=[], print_tab_count = 1):
    print('{0} Start dissolve geoprocessing on {1}...'.format('\t' * print_tab_count, infile))
    with fiona.open(infile) as input:
        with fiona.open(outfile, 'w', **input.meta) as output:
            if (len(fields) != 0):
                grouper = itemgetter(*fields)
                key = lambda k: grouper(k['properties'])
                for k, group in itertools.groupby(sorted(input, key=key), key):
                    properties, geom = zip(*[(feature['properties'], shape(feature['geometry'])) for feature in group])
                    output.write({'geometry': mapping(unary_union(geom)), 'properties': properties[0]})
            else:
                properties, geom = zip(*[(feature['properties'], shape(feature['geometry'])) for feature in input])
                output.write({'geometry': mapping(unary_union(geom)), 'properties': properties[0]})
    print('{0} Finished dissolve geoprocessing. File: {1}'.format('\t' * print_tab_count, outfile))

def save_observation_data(data):
    # https://www.psycopg.org/docs/usage.html
    # Connect to the postgres DB and open a cursor to perform operations
    # Password hardcoded, but it was kept simple as this script is not meant for production release.
    with psycopg2.connect('dbname=geoprog user=postgres password=b') as conn:
        with conn.cursor() as cur:
            # Execute queries
            for feature in data.features:
                geom = geojson.dumps(feature.geometry)
                properties = feature.properties
                try:
                    cur.execute("""INSERT INTO main.stationdata (geom, intensidadeVentoKM, temperatura, idEstacao, pressao, humidade, localEstacao, precAcumulada, idDireccVento, radiacao, time, intensidadeVento, descDirVento)
                                VALUES (ST_GeomFromGeoJSON(%s), %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s);"""
                                , (geom, properties['intensidadeVentoKM'], properties['temperatura'], properties['idEstacao'], properties['pressao'], properties['humidade'], properties['localEstacao'],properties['precAcumulada'],properties['idDireccVento'],properties['radiacao'],properties['time'], properties['intensidadeVento'], properties['descDirVento']))
                except psycopg2.errors.UniqueViolation:
                    # Ignore this error, happens many times in tests. 
                    # Should never happen with IPMA's own data because the script should run every X hours on a schedule.
                    # In any case, if it happens, it's covered here and the data won't be inserted again.
                    continue 
                finally:
                    # Make the changes to the database permanent, close the transaction
                    conn.commit()            

def import_ipma_data():
    # This function fully utilizes GDAL instead of fiona and shapely
    # as it is very efficient to clip station data and will be used to generate
    # the interpolation and then also to clip the interpolation raster.
    # GDAL itself was already installed by fiona, so it's not an extra library.

    print("\t Importing data from the last 3 hours from IPMA's API ...")
    
    # First obtain the data
    url = 'https://api.ipma.pt/open-data/observation/meteorology/stations/obs-surface.geojson'
    r = urllib.request.urlopen(url)
    observation_data = geojson.loads(r.read())

    # Save data to the database
    print ("\t Saving data to the database ...")
    save_observation_data(observation_data)

    print("\t Import finished")

def interpolate_temperature(shapefile, output_path, tmp_path):

    print ("\t Start interpolation process:")
    # First, let's create a boundary of Portugal Continental
    print("\t\t Creating boundary of Portugal Continental...")
    boundary_file = os.path.join(output_path, "boundary_portugal_continental.shp")
    dissolve(shapefile, boundary_file, print_tab_count = 2)
    
    driver = ogr.GetDriverByName("ESRI Shapefile")
    dataSource = driver.Open(boundary_file, 0)
    shpLyr = dataSource.GetLayer()

    for feature in shpLyr:
        boundary = feature.GetGeometryRef()
        bounds = boundary.GetEnvelope()
    
    print ("\t\t Connecting to the database to retrieve the observed data...")
    # Putting the password directly on the code is not a best practice, but it's simplified to speed up development
    connString = 'PG: host={0} dbname={1} user={2} password={3}'.format('localhost', 'geoprog', 'postgres', 'b')
    conn = ogr.Open(connString)
    statement = """ SELECT geom, temperatura, localestacao
                    FROM main.stationdata
                    WHERE time = 
                    (
                        SELECT max(time)
                        FROM main.stationdata
                    )
                    AND temperatura <> -99.0 """
    lyr = conn.ExecuteSQL(statement)

    print ("\t\t Data retrieved, clip data to remove data outside the boundary of Portugal Continental...")
    lyr.SetSpatialFilter(boundary)

    print ("\t\t Generating the shapefile that has the stations data. Needed for the interpolation...")
    stations_shapefile = os.path.join(output_path, 'stationpoints.shp')

    # Generate the point shapefile 
    data_source = driver.CreateDataSource(stations_shapefile)

    # Create the Spatial Reference, WGS84
    srs = osr.SpatialReference()
    srs.ImportFromEPSG(4326)

    # create the layer
    layer = data_source.CreateLayer('points', srs, ogr.wkbPoint)

    # Add the fields we're interested in
    layer.CreateField(ogr.FieldDefn('Temp', ogr.OFTReal))
    field_local = ogr.FieldDefn('Local', ogr.OFTString)
    field_local.SetWidth(80)
    layer.CreateField(field_local)

    for lyr_feature in lyr:
        # create the feature
        feature = ogr.Feature(layer.GetLayerDefn())
        # Set the attributes using the values from the postgis layer
        feature.SetField('Temp', lyr_feature.GetField('temperatura'))
        feature.SetField('Local', lyr_feature.GetField('localestacao'))

        feature.SetGeometry(lyr_feature.GetGeometryRef())

        layer.CreateFeature(feature)

        feature = None

    data_source = None
    conn = None

    print ("\t\t Shapefile generated, interpolating temperature with Inverse Distance Weighted (IDW)...")

    # The bounds used by gdal grid function have a different ordering from the envelope given by gdal itself. 
    # outputBounds have the following order [minx, miny, maxx, maxy], gdal's envelope is [minx, maxx, miny, maxy]
    bounds_corrected = [bounds[0], bounds[2], bounds[1], bounds[3]]
    
    # Create an interpolation with IDW and save an uncliped file to the tmp folder.
    idw_file = os.path.join(tmp_path, "invdist.tif") 
    gdal.Grid(idw_file, stations_shapefile, zfield="Temp", algorithm = "invdist", outputBounds=bounds_corrected, width=800, height = 800)

    print ("\t\t Clipping raster to the boundary of Portugal Continental...")
    idw_clipped_file = os.path.join(output_path, "invdist_clip.tif")
    ds = gdal.Open(idw_file)
    gdal.Warp(idw_clipped_file, ds, cutlineDSName = boundary_file, cropToCutline= True, dstNodata = "-999")
    ds = None

    print ("\t Finished interpolation process.")
    return idw_clipped_file

def generate_zonal_stats(shapefile, idw, output_path):
    print ("Calculate zonal stats to get a mean value of the temperature in every county.")
    # Calculate zonal statistics to get the mean temperature in every Country
    stats = zonal_stats(shapefile, idw, stats="mean", geojson_out = True, nodata = -999)

    # Generate shapefile to allow for further analysis. The geojson generated by the zonal_stats is too big and slow to process on GIS tools
    schema = {
        'geometry': 'Polygon',
        'properties': {
            'meantemp': 'float'
        }
    }
    zonal_stats_path = os.path.join(output_path, "zonal_stats.shp")
    with fiona.open(zonal_stats_path, 'w', 'ESRI Shapefile', schema, from_epsg(4326)) as output:
        for stat in stats:
            feature = {
                'geometry': stat['geometry'],
                'properties': {
                    'meantemp': stat['properties']['mean']
                }
            }
            output.write(feature)
    
def start_processing():
    BASE_PATH = os.getcwd()
    DATA_PATH = os.path.join(BASE_PATH, "dados")
    OUTPUT_PATH = os.path.join(BASE_PATH, "output")
    TMP_PATH = os.path.join(BASE_PATH, "tmp")
    CAOP_FILE = os.path.join(DATA_PATH, 'Cont_AAD_CAOP2021', 'Cont_AAD_CAOP2021.shp')
    REPROJECTED_CAOP_FILE = os.path.join(DATA_PATH, 'Cont_AAD_CAOP2021', 'Cont_AAD_CAOP2021_Reprojected.shp')
    COUNTY_CAOP_FILE = os.path.join(OUTPUT_PATH, 'Cont_AAD_CAOP2021_Concelhos.shp')

    # First create or recreate the output folder to store the relevant processed data
    # Then create a temporary folder to store files that are needed for intermediate processing
    create_folders([OUTPUT_PATH, TMP_PATH])

    # This script prefers WGS 84, so first the CAOP needs to be reprojected from EPSG 3763 to 4326
    reproject(CAOP_FILE, REPROJECTED_CAOP_FILE, 4326)

    # Dissolve CAOP, remove every parish (Freguesia)
    dissolve(REPROJECTED_CAOP_FILE, COUNTY_CAOP_FILE, ['Concelho'])

    # Retrieve the data stations data from IPMA
    import_ipma_data()

    # Interpolate temperature
    temp_idw_raster = interpolate_temperature(COUNTY_CAOP_FILE, OUTPUT_PATH, TMP_PATH)

    generate_zonal_stats(COUNTY_CAOP_FILE, temp_idw_raster, OUTPUT_PATH)

    # Delete the tmp folder
    delete_folders([TMP_PATH])
            
# https://www.freecodecamp.org/news/if-name-main-python-example/
if __name__ == '__main__':
    start_processing()
    