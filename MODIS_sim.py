# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""

import sys
sys.path.append('C:/Users/hzm0045/AppData/Local/Programs/Python/Python36/Lib/site-packages')
import pymodis

#from pymodis import convertmodis_gdal
import glob, os
import re
import gdal, ogr, osr, sys, numpy
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
import argparse

parser = argparse.ArgumentParser()

###Parameters with default values
parser.add_argument('-t', '--tiles', nargs='+',help="Enter a list the MODIS tiles; see https://modis-land.gsfc.nasa.gov/MODLAND_grid.html", default=['h13v11,h13v10'])
parser.add_argument('-p', '--products', nargs='+', help="Enter a list of the MODIS Land Product; see https://lpdaac.usgs.gov/dataset_discovery/modis/modis_products_table",default=['MOD13Q1.006'])
parser.add_argument('-c', '--epsg', help="Enter the new coordinate system",default='29192')
parser.add_argument('-hd', '--hdfd', help="Enter the directory where the hdf files will be downloaded",default='C:/MODIS/RawImages/')
parser.add_argument('-td', '--tifd', help="Enter directory wher the projected tif files will be located",default='C:/MODIS/NDVI_EVI_LAI/')

###Mandatory parameters
parser.add_argument('-y', '--years', nargs='+', help="Enter a list of years")
parser.add_argument('-s', '--start', help="Enter the start day and month using format mm-dd")
parser.add_argument('-e', '--end', help="Enter the end day and month using format mm-dd")
parser.add_argument('-sf', '--shp', help="Enter the shapefile",default='C:/MODIS/Shape_files/canasat_2011_few.shp')

args = parser.parse_args()

print(args.tiles)
def make_list_txt_files(ext, Directory):
    '''
    This is to create a file with the name of the files in an specific directory and with some extension
    '''
    name_file=Directory + 'list_' + ext + '_files.txt'
    os.chdir(Directory)
    text_file = open(name_file, "w")
    for file in glob.glob('*.' + ext):
        file=Directory + file + ' \n'
        text_file.write(file)
    text_file.close()
    return name_file

def get_the_name_identifiers(infile,regex,caso):
    '''
    Get the unique string identifiers from a list 
    '''
    list_files=open(infile,'r')
    my_strings=[]
    for my_string in list_files:        
        cn_strings=re.findall(regex,my_string) #cn of components of the name
        if (caso==1):
            my_strings.append(cn_strings[0])
        if (caso==2):
            my_strings.append(''.join(cn_strings[0]))
    u_strings=list(set(my_strings)) #u unique
    return u_strings, my_strings


def mosaic_hdf_files(listfile):
    regex=r'(MOD.*?).h\d{2}.'
    mrtpath='C:/MODIS_Viejo/MRT/'
    u_strings, my_strings =get_the_name_identifiers(listfile,regex,1)
    len_tab = len(u_strings)
    t_file=args.hdfd + 'temporal_file.txt' # de temporal
    for j in range(len_tab): ## as many outputs (figures) as unique product identifiers
        i_tab=[i for i, e in enumerate(my_strings) if e == u_strings[j]]               
        with open(listfile) as f:
            lines = f.readlines()
            lines_sub = [lines[i] for i in i_tab]
            with open(t_file, "w") as f1:
                f1.writelines(lines_sub)
        print(u_strings[j])
        mos=pymodis.convertmodis.createMosaic(listfile=t_file,outprefix = u_strings[j],mrtpath=mrtpath)
        mos.run()
        with open(t_file, "r") as f2:
            for item in f2:
                path_file=args.hdfd+item[:-1]
                os.remove(path_file)
  
def down_modis():
    '''
    Here I download the images
    '''
    print(args.tiles)
    print(args.years)
    ###Define dir paths
    
    for i in range(int(args.years[0]),(int(args.years[-1])+1)):
        for tiles in args.tiles:
            for product in args.products:
                today=str(i)+'-'+str(args.start)
                enddate = str(i)+'-'+str(args.end)
                user='Fabiormarin'
                passwd='Esalq2018'
                dest=args.hdfd
                listfile=dest + 'listfile' + product + '.txt'
                print('product is: ',product)
                print('from:', today)
                print('to: ', enddate)
                print(type(args.tiles))
                print('dest={}, tiles={}, today={}, enddate={}, user = {}, password = {}, product = {}'.format(dest, args.tiles ,today,enddate,user,passwd,product))
     
                modisDown = pymodis.downmodis.downModis(destinationFolder=dest, tiles=tiles, today=today, enddate=enddate, user = user, password = passwd, product = product)
                modisDown.connect()
                modisDown.downloadsAllDay()
                if (len(tiles)>8):
                    mosaic_hdf_files(listfile)

    file_list=make_list_txt_files('hdf',dest)
    return file_list

def convert_modis(list_txt_file):
    '''
    Here I change change the coordinate system of the images 
    '''
    sub_MOD13Q1='(1 1 0 0 0 0 0 0 0 0 0 0)'
    sub_MOD15A2H='(0 1 0 0 0 0)' 
    list_files=open(list_txt_file,'r')
    regex_prod=r'(MOD\d{2}.+).A'
    regex_pref=r'(MOD\d{2}.+).hdf'
    for file in list_files:
        inf=file[:-1]
        print('inf is: ', inf)
        product=re.findall(regex_prod,file)
        print('product is: ', product)
        prefix=re.findall(regex_pref,file)
        prefix=args.tifd+ prefix[0]
        print('prefix is: ',prefix)
        if product[0]=='MOD13Q1':
            out=pymodis.convertmodis_gdal.convertModisGDAL(hdfname=inf, prefix=prefix, subset=sub_MOD13Q1, res=250, outformat='GTiff', epsg=args.epsg)
            out.run(quiet=False)
        elif product[0]=='MOD15A2H':
            out=pymodis.convertmodis_gdal.convertModisGDAL(hdfname=inf, prefix=prefix, subset=sub_MOD15A2H, res=500, outformat='GTiff', epsg=args.epsg)
            out.run(quiet=False)                


def spatially_averaged_vegindex():
    '''
    Here I get and spatially averaged output
    '''
    regex=r'(MOD\d.+).A(\d{4}).+ (.+).tif'
    u_tablas, tablas = get_the_name_identifiers(make_list_txt_files('tif',args.tifd),regex,2)
    len_tab = len(u_tablas)
    myshapefile = args.shp
    cuenta=0
    for j in range(len_tab): ## as many outputs (figures) as unique product identifiers
        i_tab=[i for i, e in enumerate(tablas) if e == u_tablas[j]]
        name_txt_file=make_list_txt_files('tif',args.tifd) 
        list_files=open(name_txt_file,'r')
        lines=list_files.readlines()
        lines_sub=[lines[i] for i in i_tab]
        julians=[]
        ave_output=[]
        cuenta+=1
        regex=r'\.A\d{4}(\d{3})' # search the julian day
        regex1=r'\d{4}(.+)' #serach the vegetation index
        veg_ind=re.findall(regex1,u_tablas[j])
        
        for filein in lines_sub:
            print(filein)
            inf=filein[:-1]
            julian=re.findall(regex,filein)
            julians.append(int(julian[0]))
            raster = gdal.Open(inf)
            shp = ogr.Open(myshapefile)
            lyr = shp.GetLayer()
            featList = range(lyr.GetFeatureCount())
            statDict = {}        
            for FID in featList:
                feat = lyr.GetFeature(FID)
                # Get raster georeference info
                transform = raster.GetGeoTransform()
                xOrigin = transform[0]
                yOrigin = transform[3]
                pixelWidth = transform[1]
                pixelHeight = transform[5]
                
                # Reproject vector geometry to same projection as raster
                sourceSR = lyr.GetSpatialRef()
                targetSR = osr.SpatialReference()
                targetSR.ImportFromWkt(raster.GetProjectionRef())
                coordTrans = osr.CoordinateTransformation(sourceSR,targetSR)
            # feat = lyr.GetNextFeature()
                feat = lyr.GetFeature(FID)
                geom = feat.GetGeometryRef()
                geom.Transform(coordTrans)
                
                # Get extent of feat
                geom = feat.GetGeometryRef()
                if (geom.GetGeometryName() == 'MULTIPOLYGON'):
                    count = 0
                    pointsX = []; pointsY = []
                    for polygon in geom:
                        geomInner = geom.GetGeometryRef(count)
                        ring = geomInner.GetGeometryRef(0)
                        numpoints = ring.GetPointCount()
                        for p in range(numpoints):
                                lon, lat, z = ring.GetPoint(p)
                                pointsX.append(lon)
                                pointsY.append(lat)
                        count += 1
                elif (geom.GetGeometryName() == 'POLYGON'):
                    ring = geom.GetGeometryRef(0)
                    numpoints = ring.GetPointCount()
                    pointsX = []; pointsY = []
                    for p in range(numpoints):
                            lon, lat, z = ring.GetPoint(p)
                            pointsX.append(lon)
                            pointsY.append(lat)
                
                else:
                    sys.exit("ERROR: Geometry needs to be either Polygon or Multipolygon")
                
                xmin = min(pointsX)
                xmax = max(pointsX)
                ymin = min(pointsY)
                ymax = max(pointsY)
                
                # Specify offset and rows and columns to read
                xoff = int((xmin - xOrigin)/pixelWidth)
                yoff = int((yOrigin - ymax)/pixelWidth)
                xcount = int((xmax - xmin)/pixelWidth)+1
                ycount = int((ymax - ymin)/pixelWidth)+1
                
                # Create memory target raster
                target_ds = gdal.GetDriverByName('MEM').Create('', xcount, ycount, 1, gdal.GDT_Byte)
                target_ds.SetGeoTransform((
                    xmin, pixelWidth, 0,
                    ymax, 0, pixelHeight,
                ))
                
                # Create for target raster the same projection as for the value raster
                raster_srs = osr.SpatialReference()
                raster_srs.ImportFromWkt(raster.GetProjectionRef())
                target_ds.SetProjection(raster_srs.ExportToWkt())
                
                # Rasterize zone polygon to raster
                gdal.RasterizeLayer(target_ds, [1], lyr, burn_values=[1])
                
                # Read raster as arrays
                banddataraster = raster.GetRasterBand(1)
                dataraster = banddataraster.ReadAsArray(xoff, yoff, xcount, ycount).astype(numpy.float)
                
                bandmask = target_ds.GetRasterBand(1)
                datamask = bandmask.ReadAsArray(0, 0, xcount, ycount).astype(numpy.float)
                
                # Mask zone of raster
                zoneraster = numpy.ma.masked_array(dataraster,  numpy.logical_not(datamask))
                
                # Calculate statistics of zonal raster
                meanValue=zoneraster[zoneraster.nonzero()].mean()
                if meanValue>0:
                    statDict[FID] = meanValue
            ave_output.append(sum(statDict.values()) )#/float(len(statDict)))
            print(ave_output)
            if (len_tab % 2 >0):
                plt.subplot(2,(len_tab-1)/2+1,cuenta)
                plt.scatter(julians,ave_output)
            # plt.title(u_tablas[j])
                plt.xlabel('Julian day')
                plt.ylabel(veg_ind)
            if (len_tab % 2 == 0):
                plt.subplot(2,(len_tab)/2,cuenta)
                plt.scatter(julians,ave_output)
            # plt.title(u_tablas[j])
                plt.xlabel('Julian day')
                plt.ylabel(veg_ind)
                
    plt.savefig('C:/MODIS/Shape_files/Veg_index_dyn.pdf')


#file_with_list=down_modis()
#print('file_with_list is: ',file_with_list)
file_with_list='C:/MODIS/RawImages/list_hdf_files.txt'
convert_modis(file_with_list)
spatially_averaged_vegindex()

#python C:/Pythoncourse_summer18/MODIS.py -t h13v11 -p MOD13A1.006 -y 2017 -s 03-01 -e 09-12 -c 29192 -hd C:/MODIS/RawImages/ -td C:/MODIS/NDVI_EVI_LAI/ -sf C:/MODIS/Areas_canna/canasat_2011_few.shp