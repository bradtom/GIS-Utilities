import sys
import arcpy, sys, os, math, glob
from arcpy import env
import dbfpy
from dbfpy import dbf
import struct, datetime, decimal, itertools
# Check out the ArcGIS 3D Analyst extension license
from arcgisscripting import ExecuteError
arcpy.CheckOutExtension("3D")
arcpy.CheckOutExtension("Spatial")
print "Check extension OK"
from shutil import copyfile


class CalcGISVolume:
    """
    This is the current version as of 3/12/2019. This is the version that was used.
    Mostly written by Ines, with some modifications by Brad
    Calculates volume and the 2D area using Cut & Fill.
    Inputs are a shapefile containing polygons representing DSM2 Channels, and a number of DEMs
    Polygons are selected one at a time, then for the selcted polygon:
        a raster is created using Polygon To Raster
        Volume is calculated using CutFill
    Results are written to ASCII and dbf files.
    """

    INPUT_POLYGON_DIR = 'D:/dsm2GisReference/data_source_dem_v4'
    DEM_INPUT_DIR = 'E:/dsm2GisReference/bay_delta_dem_v4.gdb'
    GEOPROCESSING_OUTPUT_ROOT_DIR = 'E:/dsm2GisReference/dsm2ChanVolCalc'
    RESULTS_OUTPUT_DIR = GEOPROCESSING_OUTPUT_ROOT_DIR+'/GeoProOutput'
    POLYGON_RASTERS_GDB = GEOPROCESSING_OUTPUT_ROOT_DIR+'/PolygonRasters.gdb'
    CUT_AND_FILL_GDB = GEOPROCESSING_OUTPUT_ROOT_DIR+'/VolumeCalcCutFill.gdb'

    def __init__(self):
        # stores rasters here, where there is more space
        os.environ["TMP"] = "E:/ArcGISTempFolder/"

        polygonShapefilePath = CalcGISVolume.INPUT_POLYGON_DIR + "/modelGridChannelPolygons.shp"

        # polygonShapefilePath = CalcGISVolumeModified.GISbasedata_dir + "/testChan150.shp"
        field_name = "Shape"
        # dbfPath = CalcGISVolumeModified.RESULTS_OUTPUT_DIR + "/" + field_name + "DataTable.dbf"
        # dbfFile = open(dbfPath, 'wb')
        # asciiFilename = CalcGISVolumeModified.RESULTS_OUTPUT_DIR + "/ChannelConveyanceCharacteristics.csv"
        # asciiFile = open(asciiFilename, 'w')

        # creating polygon rasters takes a long time, and should not need to be done more than once. Only set this
        # to true if they need to be created again.
        repeatPolygonRasterCreation = False

        # print "Computing volumes for ", field_name
        # # Make feature layer from Polygon (either channel or reservoir).
        # # The layer has the same data as the in-features, but is required for the next step.
        # out_layer = field_name + "_lyr"
        # arcpy.MakeFeatureLayer_management(polygonShapefilePath, out_layer)
        # print "Feature Layer OK"

        if repeatPolygonRasterCreation:
            if arcpy.Exists(CalcGISVolume.POLYGON_RASTERS_GDB):
                arcpy.Delete_management(CalcGISVolume.POLYGON_RASTERS_GDB)
            # Execute CreateFileGDB
            if not arcpy.Exists("PolygonRasters.gdb"):
                arcpy.CreateFileGDB_management(CalcGISVolume.GEOPROCESSING_OUTPUT_ROOT_DIR, "PolygonRasters.gdb")

        DEMFiles = [
            # 'dem_ccfb_south_delta_san_joaquin_rvr_2m_20180612',
            # 'dem_columbia_cut_2m_20120911',
            # 'dem_false_rvr_piper_sl_fisherman_cut_2m_20171109',
            # 'dem_montezuma_sl_2m_20180604',
            # 'dem_north_delta_2m_20171220',
            # 'dem_Sanjoaquin_bradfords_2m_20151204_clip',
            # 'dem_turner_cut_2m_20120907',
            # 'dem_yolo_merge3_mod_20161202',
            'dem_delta_10m_20180615'
        ]
        for DEMFile in DEMFiles:
            dbfPath = CalcGISVolume.RESULTS_OUTPUT_DIR + "/" + DEMFile + "_DataTable.dbf"
            dbfFile = open(dbfPath, 'wb')
            asciiFilename = CalcGISVolume.RESULTS_OUTPUT_DIR + "/" + DEMFile + "_CutFillResults.csv"
            asciiFile = open(asciiFilename, 'w')
            # Delete the cut/fill geodatabase
            # Check if Volume calc file geodatabases and file folders exist.  If they do, delete old and create new one.
            if arcpy.Exists(CalcGISVolume.CUT_AND_FILL_GDB):
                arcpy.Delete_management(CalcGISVolume.CUT_AND_FILL_GDB)
            # Execute CreateFileGDB
            # out_folder_path = CalcGISVolumeModified.root_dir
            if not arcpy.Exists("VolumeCalcCutFill.gdb"):
                arcpy.CreateFileGDB_management(CalcGISVolume.GEOPROCESSING_OUTPUT_ROOT_DIR, "VolumeCalcCutFill.gdb")

            # DEMFile = 'dem_ccfb_south_delta_san_joaquin_rvr_2m_201806121.tif'
            DEMPath = CalcGISVolume.DEM_INPUT_DIR + "/" + DEMFile
            print "about to call computeGISVolume: DEMPath="+DEMPath
            print "about to call computeGISVolume: polygonShapefilePath="+polygonShapefilePath
            print "about to call computeGISVolume: field_name="+field_name
            print "about to call computeGISVolume: asciiFilename="+asciiFilename
            self.computeGISVolume(repeatPolygonRasterCreation, DEMPath, field_name, polygonShapefilePath, asciiFile, dbfFile)
            # now that the rasters are created, turn off raster creation
            repeatPolygonRasterCreation = False
            asciiFile.close()
            dbfFile.close()
        print "done"

    def computeGISVolume(self, repeatPolygonRasterCreation, DEMfile, field_name, polygonShapefilePath, asciiFile, dbfFile):
        """
        Adapted from Ines's code
        :param field_name:
        :param channelPolygonLayerPath:
        :param demPath:

        :return:
        """

        # Start index for channels
        chanIndex = -1
        polygonNames = {}
        polygonAreas = {}
        cutFillAreas = {}

        max_vol = {}
        channel = []

        print "Computing volumes for ", field_name
        # Make feature layer from Polygon (either channel or reservoir).
        # The layer has the same data as the in-features, but is required for the next step.
        out_layer = field_name + "_lyr"
        # this line will throw an exception because the layer already exists. I tried enclosing it in a try/except,
        # but the resulted in the layer only containing one polygon, so still need to run this script one dem at a time
        arcpy.MakeFeatureLayer_management(polygonShapefilePath, out_layer)
        print "Feature Layer made"

        # Set cursor to loop through each channel top area in the layer just created
        # (polygon to raster must be done one-by-one, otherwise it creates a single raster)
        polygonCursor = arcpy.da.SearchCursor(out_layer, ["FID", "id", "Area"])

        for row in polygonCursor:
            print "==========================================================================="
            chanIndex += 1
            print "chanIndex = " + str(chanIndex)
            print row[0]
            pIndex = row[0]
            pName = row[1]
            pArea = row[3]
            print "polygon index, name, area="+str(pIndex)+","+str(pName)+","+str(pArea)
            print "polygonIndex = " + str(pIndex)

            # Select current channel (or reservoir) top area
            InLayer = out_layer
            Select = "NEW_SELECTION"
            query = 'FId=%s' % pIndex
            print 'FId=%s' % pIndex
            # query = 'id=%s' % polygonIndex
            # print 'id=%s' % polygonIndex
            # Execute SelectLayerByAttribute
            arcpy.SelectLayerByAttribute_management(InLayer,Select, query)
            print "Selection ok"

            # Convert selected Polygon in Polygons to Raster (Polygon to Raster with cell size 10).
            if repeatPolygonRasterCreation:
                FCname = field_name + str(pIndex)
                # ValueField = "Z"
                ValueField = "Elevation"
                OutRaster = CalcGISVolume.POLYGON_RASTERS_GDB + '/' + FCname
                print "OutRaster", OutRaster
                cell_assignment="CELL_CENTER"
                priority_field="NONE"
                cellsize="2"
                # Execute PolygonToRaster
                # print "InLayer=", InLayer, "  ValueField=", ValueField, "  OutRaster=", OutRaster,"cell_assignment=", cell_assignment, "priority_field & cellsize=", priority_field, cellsize
                arcpy.PolygonToRaster_conversion(InLayer, ValueField, OutRaster, cell_assignment, priority_field, cellsize)
                # print "polygon to raster ok"

            # arcpy.SurfaceVolume_3d(in_surface=OutRaster, out_text_file=root_dir+"/surfaceVolume"+str(chanIndex)+".txt",
            #                        reference_plane="ABOVE", base_z="0.0", z_factor="1", pyramid_level_resolution="0")
            # Now that we have rasters for each channel, we want to find the volume between the ChannelPolygons raster and the DEM (topo/bathymetry)
            # use the arcpy.CutFill tool and then select the highest positive volume
            in_before_surface = CalcGISVolume.POLYGON_RASTERS_GDB + '/' + field_name + str(pIndex)
            out_raster= CalcGISVolume.CUT_AND_FILL_GDB + "/CutFill" + str(pIndex)
            print in_before_surface, out_raster

            #  Use CutFill tool
            print "Cut and fill arguments: in_before_surface = " + in_before_surface
            # print "Cut and fill arguments: in_after_surface = " + DEMfile
            # print "Cut and fill arguments: output_raster = " + out_raster
            try:
                arcpy.CutFill_3d(in_before_surface, DEMfile, out_raster, z_factor="1")
                print "CutFill ok"

                field_names = ["VOLUME", "AREA"]
                in_raster = out_raster

                max = 0
                temp_area = 0
                lines = arcpy.da.SearchCursor(in_raster, field_names)
                for line in lines:
                    if line[0] > max:
                        max = line[0]
                        temp_area = line[1]
                polygonNames[pIndex] = pName
                polygonAreas[pIndex] = pArea
                max_vol[pIndex] = max
                cutFillAreas[pIndex] = temp_area
                print "CutFill results (vol, area) for "+pName+": "+str(max_vol[pIndex])+","+\
                      str(cutFillAreas[pIndex])
                channel.append(pIndex)
            except ExecuteError:
                polygonNames[pIndex] = pName
                polygonAreas[pIndex] = pArea
                max_vol[pIndex] = -999999
                cutFillAreas[pIndex] = -999999
                channel.append(pIndex)
                print "skipping shape "+in_before_surface + " because CutFill operation failed"
            finally:
                pass

        # Write data to a DBF file
        if field_name == 'Reservoir':
            res_names = {525:"SanRafael", 526:"CMadera", 574:"RichPt", 571:"ElCerrito", 572:"Albany", 573:"Emeryville",
                            576:"Oakland", 580:"SouthBay", 510:"Pinole", 520:"SanPablo", 504:"NapaSlough",
                            516:"PetalumaR", 522:"Ignacio", 482:"NapaR", 582: "Richardson"   }
            fieldnames = [field_name, 'Volume', 'Area', 'ResName']
            fieldspecs = [('N', 9, 0), ('N', 19, 10), ('N', 19, 10), ('C', 15, 0)]
        else:
            # include name of DEM
            fieldnames = [field_name, 'Polygon', 'Volume', 'Area', 'PolygonArea', 'DEM']
            fieldspecs = [('N', 9, 0), ('C', 15, 0), ('N', 19, 10), ('N', 19, 10), ('N', 19, 10), ('C', 50, 0)]
            # fieldnames = [field_name, 'Polygon', 'Volume', 'Area', 'DEM']
            # fieldspecs = [('N', 9, 0), ('C', 15, 0), ('N', 19, 10), ('N', 19, 10), ('C', 50, 0)]

        records = [[] for i in range(len(channel))]

        # Delete old file
        # if arcpy.Exists(dbfPath):
        #     arcpy.Delete_management(dbfPath)

        # IOError: [Errno 2] No such file or directory: 'D:/dsm2GisReference/data_source_dem_v4/GeoProOutput/ShapeDataTable.dbf'
        # dbfFile = dbf.Dbf(dbfPath, new=False, readOnly=False)
        for i in range(len(channel)):
            chan_no = channel[i]
            records[i].append(chan_no)
            records[i].append(polygonNames[chan_no])
            records[i].append(max_vol[chan_no])
            records[i].append(cutFillAreas[chan_no])
            records[i].append(polygonAreas[chan_no])
            records[i].append(DEMfile)
            if field_name == 'Reservoir':
                records[i].append(res_names[chan_no])
    #        print 'records = ', records
            if max_vol[chan_no] > 0 and cutFillAreas[chan_no] >0:
                asciiFile.write(str(chan_no)+","+str(polygonNames[chan_no])+","+str(max_vol[chan_no])+","+
                                str(cutFillAreas[chan_no])+","+str(polygonAreas[chan_no])+","+DEMfile+"\n")
        # dbfFile.addField(fieldspecs, records)
        self.dbfwriter(dbfFile, fieldnames, fieldspecs, records)

    def dbfreader(self, f):
        """
        Copied from Ines

        Returns an iterator over records in a Xbase DBF file.

        The first row returned contains the field names.
        The second row contains field specs: (type, size, decimal places).
        Subsequent rows contain the data records.
        If a record is marked as deleted, it is skipped.

        File should be opened for binary reads.

        """
        # See DBF format spec at:
        #     http://www.pgts.com.au/download/public/xbase.htm#DBF_STRUCT
        print "f= ", f
        numrec, lenheader = struct.unpack('<xxxxLH22x', f.read(32))
        numfields = (lenheader - 33) // 32
        print 'numrec, lenheader, numfields', numrec, lenheader, numfields

        fields = []
        for fieldno in xrange(numfields):
            name, typ, size, deci = struct.unpack('<11sc4xBB14x', f.read(32))
            name = name.replace('\0', '')       # eliminate NULs from string
            fields.append((name, typ, size, deci))
        yield [field[0] for field in fields]
        yield [tuple(field[1:]) for field in fields]

        terminator = f.read(1)
        assert terminator == '\r'

        fields.insert(0, ('DeletionFlag', 'C', 1, 0))
        fmt = ''.join(['%ds' % fieldinfo[2] for fieldinfo in fields])
        fmtsiz = struct.calcsize(fmt)
        for i in xrange(numrec):
            record = struct.unpack(fmt, f.read(fmtsiz))
            if record[0] != ' ':
                continue                        # deleted record
            result = []
            for (name, typ, size, deci), value in itertools.izip(fields, record):
                if name == 'DeletionFlag':
                    continue
                if typ == "N":
                    value = value.replace('\0', '').lstrip()
                    if value == '':
                        value = 0
                    elif deci:
                        value = decimal.Decimal(value)
                    else:
                        value = int(value)
                elif typ == 'D':
                    y, m, d = int(value[:4]), int(value[4:6]), int(value[6:8])
                    value = datetime.date(y, m, d)
                elif typ == 'L':
                    value = (value in 'YyTt' and 'T') or (value in 'NnFf' and 'F') or '?'
                elif typ == 'F':
                    value = float(value)
                result.append(value)
            yield result

    def dbfwriter(self, f, fieldnames, fieldspecs, records):
        """
        Copied from Ines
        Return a string suitable for writing directly to a binary dbf file.

        File f should be open for writing in a binary mode.

        Fieldnames should be no longer than ten characters and not include \x00.
        Fieldspecs are in the form (type, size, deci) where
            type is one of:
                C for ascii character data
                M for ascii character memo data (real memo fields not supported)
                D for datetime objects
                N for ints or decimal objects
                L for logical values 'T', 'F', or '?'
            size is the field width
            deci is the number of decimal places in the provided decimal object
        Records can be an iterable over the records (sequences of field values).

        """
        # header info
        ver = 3
        now = datetime.datetime.now()
        yr, mon, day = now.year-1900, now.month, now.day
        numrec = len(records)
        numfields = len(fieldspecs)
        lenheader = numfields * 32 + 33
        lenrecord = sum(field[1] for field in fieldspecs) + 1
        hdr = struct.pack('<BBBBLHH20x', ver, yr, mon, day, numrec, lenheader, lenrecord)
        f.write(hdr)

        # field specs
        for name, (typ, size, deci) in itertools.izip(fieldnames, fieldspecs):
            name = name.ljust(11, '\x00')
            fld = struct.pack('<11sc4xBB14x', name, typ, size, deci)
            f.write(fld)

        # terminator
        f.write('\r')

        # records
        for record in records:
            f.write(' ')                        # deletion flag
            for (typ, size, deci), value in itertools.izip(fieldspecs, record):
                if typ == "N":
                    value = str(value).rjust(size, ' ')
                elif typ == 'D':
                    value = value.strftime('%Y%m%d')
                elif typ == 'L':
                    value = str(value)[0].upper()
                else:
                    value = str(value)[:size].ljust(size, ' ')
                assert len(value) == size
                f.write(value)

        # End of file
        f.write('\x1A')

CalcGISVolume()