#
#  Copyright (c) 2001-2014, Scott D. Peckham
#
#
# -----------------------------------------------------------------------
#  NOTES:  This file defines a "base class" for permafrost
#          components as well as functions used by most or
#          all permafrost methods.  The methods of this class
#          should be over-ridden as necessary for different
#          methods of modeling snowmelt.

#          update_snow_vars() in precip.py sets values here.
# -----------------------------------------------------------------------
#
#  class snow_component    (inherits from BMI_base.py)
#
#      set_constants()
#      -----------------------
#      initialize()
#      update()
#      finalize()
#      set_computed_input_vars()
#      --------------------------
#      check_input_types()
#      initialize_computed_vars()
#      ----------------------------
#      update_meltrate()
#      enforce_max_meltrate()
#      update_SM_integral()
#      update_swe()
#      update_depth()
#      -----------------------
#      open_input_files()
#      read_input_files()
#      close_input_files()
#      -----------------------
#      update_outfile_names()
#      open_output_files()
#      write_output_files()
#      close_output_files()
#      save_grids()
#      save_pixel_values()
#
# -----------------------------------------------------------------------
"""
*The MIT License (MIT)*
Copyright (c) 2016 permamodel
Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:
The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.
THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
*
"""
from __future__ import print_function

import os

import numpy as np

from permamodel.utils import BMI_base, model_input

from .. import data_directory

# import gdal
# from gdalconst import *  # Import GDAL constants, eg GA_ReadOnly
# import osr
# from pyproj import Proj, transform


# from permamodel.utils import model_output


# -----------------------------------------------------------------------
class PermafrostComponent(BMI_base.BMI_component):
    # ------------------------------------------------------------
    # Notes: rho_H2O, Cp_snow, rho_air and Cp_air are currently
    #        hardwired (not adjustable with GUI).
    # ------------------------------------------------------------
    def set_constants(self):
        # -----------------------------------
        # Constants not changeable by user
        # ---------------------------------------------------------
        # Cp_snow = mass-specific isobaric heat capacity of snow
        #           value from:  NCAR CSM Flux Coupler web page
        # ---------------------------------------------------------
        # Lf = latent heat of fusion for water [J kg -1]
        # ---------------------------------------------------------
        self.Cp_snow = np.float64(2090.0)  # [J kg-1 K-1]
        self.Lf = np.float64(334000)  # [J kg-1]

        # --------------------------------------
        # Not a constant; read from CFG file.
        # --------------------------------------
        ## self.rho_snow = np.float64(300)
        ## self.rho_H2O  = np.float64(1000)  # (See initialize() method.)

    #   set_constants()
    # -------------------------------------------------------------------

    def open_input_files(self):
        # ------------------------------------------------------
        # Each component that inherits from snow_base.py must
        # implement its own versions of these.
        # ------------------------------------------------------
        print("ERROR: open_input_files() for permafrost component")
        print("       has not been implemented for this component.")

    #   open_input_files()
    # -------------------------------------------------------------------

    def extract_grid_value_from_GSD(
        self, input_file, lonname, lon_grid_scale, latname, lat_grid_scale, varname
    ):
        """
        The function is to extract the grid value from NetCDF file,
        according to input of latitude and longitude;

        INPUTs:
            input_lat: Latitude;
            input_lon: Longitude;
            input_file: grid data file (NetCDF file)
            lonname: name of longitude in "input_file"
            latname: name of latitude in "input_file"
            lon_grid_scale: grid size of longitude
            lat_grid_scale: grid size of latitude
            varname: name of variable should be extracted.

        OUTPUTs:
            p_data: grid value
        """

        from netCDF4 import Dataset

        # import numpy as np
        # Read the nc file

        fh = Dataset(input_file, mode="r")

        # Get the lat and lon
        #   Set the grid size for lat. and lon. (here is 0.5 degree)

        lon_grid = fh.variables[lonname][:]
        lat_grid = fh.variables[latname][:]

        # Get boundary of each grid
        #    Including top and bottom of latitude
        #              top and bottom of longitude

        lon_grid_top = lon_grid + lon_grid_scale / 2.0
        lat_grid_top = lat_grid + lat_grid_scale / 2.0

        lon_grid_bot = lon_grid - lon_grid_scale / 2.0
        lat_grid_bot = lat_grid - lat_grid_scale / 2.0

        # Get the index of input location acccording to lat and lon inputed
        idx_lon = np.where((self.lon <= lon_grid_top) & (self.lon > lon_grid_bot))
        idx_lat = np.where((self.lat <= lat_grid_top) & (self.lat > lat_grid_bot))

        idx_lon = np.array(idx_lon)
        idx_lat = np.array(idx_lat)

        p_data = fh.variables[varname][idx_lat[0, 0], idx_lon[0, 0]]

        # fh.close()

        return p_data

    #   extract_grid_value_from_GSD()
    # -------------------------------------------------------------------

    def get_param_nc4_filename(self, filename_root):
        filename = os.path.join(data_directory, "%s.nc4" % filename_root)
        if not os.path.isfile(filename):
            print(("Error: File %s does not exist" % filename))
            exit(-1)
        return filename

    def initialize_from_config_file(self, cfg_file=None):
        # ---------------------------------------------------------
        # Notes:  Need to make sure than h_swe matches h_snow ?
        #         User may have entered incompatible values.
        # ---------------------------------------------------------

        # SILENT and mode were original optional arguments
        #   they should be removed, but are simply defined here for brevity

        SILENT = True

        if not SILENT:
            print("Permafrost component: Initializing...")

        self.status = "initializing"  # (OpenMI 2.0 convention)

        # Set the cfg file if it exists, otherwise, a default
        if os.path.isfile(cfg_file):
            if not SILENT:
                print(("passed cfg_file: %s" % cfg_file))
            self.cfg_file = cfg_file
        else:
            cfg_file = (
                "./permamodel/examples/Frostnumber_example_singlesite_singleyear.cfg"
            )
            print("No valid configuration file specified, trying: ")
            print(("   %s" % cfg_file))
            self.cfg_file = cfg_file

            if os.path.isfile(cfg_file):
                print(("Default config file exists: %s" % cfg_file))
            else:
                print("Default config file does not exist: ")
                print(("  %s" % cfg_file))
                raise ValueError(
                    "Default frostnumber config file {fname} does not exist".format(
                        fname=cfg_file
                    )
                )

        # print mode, cfg_file

        # -----------------------------------------------
        # Load component parameters from a config file
        # -----------------------------------------------
        self.set_constants()  # in this file, unless overridden
        self.initialize_config_vars()  # in BMI_base.py, unless overridden
        # self.initialize_config_vars() reads the configuration file
        # using read_config_file() in BMI_base.py

        # ---------------------------------------------
        # Extract soil texture from Grid Soil Database (Netcdf files)
        # according to locations
        # ---------------------------------------------
        # ScottNote: this isn't needed for all models,
        #   so I have commented it out here
        # self.initialize_soil_texture_from_GSD()

        self.initialize_time_vars()  # These time values refer to clock time

        if self.comp_status == "Disabled":
            #########################################
            #  DOUBLE CHECK THIS; SEE NOTES ABOVE
            #########################################
            ####### and (ep.method != 2):  ??????
            if not (SILENT):
                print("Permafrost component: Disabled.")
            self.lat = self.initialize_scalar(0, dtype="float64")
            self.lon = self.initialize_scalar(0, dtype="float64")
            self.T_air = self.initialize_scalar(0, dtype="float64")
            self.h_snow = self.initialize_scalar(0, dtype="float64")
            self.vwc_H2O = self.initialize_scalar(0, dtype="float64")
            self.Hvgf = self.initialize_scalar(0, dtype="float64")
            self.Hvgt = self.initialize_scalar(0, dtype="float64")
            self.Dvf = self.initialize_scalar(0, dtype="float64")
            self.Dvt = self.initialize_scalar(0, dtype="float64")
            self.start_year = self.initialize_scalar(0, dtype="float64")
            self.DONE = True
            self.status = "initialized"
            return

        # ---------------------------------------------
        # Open input files needed to initialize vars
        # ---------------------------------------------
        self.open_input_files()
        self.read_input_files()

        # ---------------------------
        # Initialize computed vars
        # ---------------------------
        self.check_input_types()  # (maybe not used yet)

        self.status = "initialized"

    #   initialize()
    # -------------------------------------------------------------------
    def initialize_soil_texture_from_GSD(self):
        # ScottNote: this should be moved to the component that needs it

        Clay_file = self.get_param_nc4_filename("T_CLAY")
        Sand_file = self.get_param_nc4_filename("T_SAND")
        Silt_file = self.get_param_nc4_filename("T_SILT")
        Peat_file = self.get_param_nc4_filename("T_OC")

        # print("Clay_file: %s" % Clay_file)
        # print("Sand_file: %s" % Sand_file)
        # print("Silt_file: %s" % Silt_file)
        # print("Peat_file: %s" % Peat_file)

        lonname = "lon"
        lon_grid_scale = 0.05
        latname = "lat"
        lat_grid_scale = 0.05

        self.p_clay = self.extract_grid_value_from_GSD(
            Clay_file, lonname, lon_grid_scale, latname, lat_grid_scale, "T_CLAY"
        )
        self.p_sand = self.extract_grid_value_from_GSD(
            Sand_file, lonname, lon_grid_scale, latname, lat_grid_scale, "T_SAND"
        )
        self.p_peat = self.extract_grid_value_from_GSD(
            Peat_file, lonname, lon_grid_scale, latname, lat_grid_scale, "T_OC"
        )
        self.p_silt = self.extract_grid_value_from_GSD(
            Silt_file, lonname, lon_grid_scale, latname, lat_grid_scale, "T_SILT"
        )
        self.p_peat = 0

    #   initialize_soil_texture_from_GSD()
    # -------------------------------------------------------------------

    def read_input_files(self):
        print("ERROR: read_input_files() for permafrost component")
        print("       has not been implemented.")
        print(" Ignoring this message for now. Can not find the rti class")

    #   read_input_files()
    # -------------------------------------------------------------------
    def check_input_types(self):
        # ----------------------------------------------------
        # Notes: Usually this will be overridden by a given
        #        method of computing snow meltrate.
        # ----------------------------------------------------
        # Notes: rho_H2O, Cp_snow, rho_air and Cp_air are
        #        currently always scalars.
        # ----------------------------------------------------
        are_scalars = np.array(
            [
                self.is_scalar("lat"),
                self.is_scalar("lon"),
                self.is_scalar("T_air"),
                self.is_scalar("A_air"),
                self.is_scalar("h_snow"),
                self.is_scalar("rho_snow"),
                # ----------------------------------
                self.is_scalar("vwc_H2O"),
                self.is_scalar("Hvgf"),
                self.is_scalar("Hvgt"),
                self.is_scalar("Dvf"),
                self.is_scalar("Dvt"),
            ]
        )

        self.ALL_SCALARS = np.all(are_scalars)

    #   check_input_types()
    # -------------------------------------------------------------------

    def update_ALT(self):
        # ---------------------------------------------------------
        # Note: We don't need to update any variables if
        #       the soil_thermal_conductivity method is None.  But we need
        #       to make sure that thermal_conductivity = 0.0.
        #       This "method" will be over-ridden by a
        #       particular snowmelt method.
        # --------------------------------------------------
        print("ERROR: update_ALT() method for Permafrost component")
        print("       has not been implemented.")

    #   update_ALT()
    # -------------------------------------------------------------------
    def update_ground_temperatures(self):
        # ---------------------------------------------------------
        # Note: We don't need to update any variables if
        #       the soil_thermal_conductivity method is None.  But we need
        #       to make sure that thermal_conductivity = 0.0.
        #       This "method" will be over-ridden by a
        #       particular snowmelt method.
        # --------------------------------------------------
        print("ERROR: update_ground_temperatures method for Permafrost component")
        print("       has not been implemented.")

    #   update_ground_temperatures()
    # -------------------------------------------------------------------

    ## def update(self, dt=-1.0, time_seconds=None):
    def update(self, dt=-1.0):
        # ----------------------------------------------------------
        # Note: The read_input_files() method is first called by
        #       the initialize() method.  Then, the update()
        #       method is called one or more times, and it calls
        #       other update_*() methods to compute additional
        #       variables using input data that was last read.
        #       Based on this pattern, read_input_files() should
        #       be called at end of update() method as done here.
        #       If the input files don't contain any additional
        #       data, the last data read persists by default.
        # ----------------------------------------------------------

        # -------------------------------------------------
        # Note: self.SM already set to 0 by initialize()
        # -------------------------------------------------
        if self.comp_status == "Disabled":
            return
        self.status = "updating"  # (OpenMI)

        # -------------------------
        # Update computed values
        # -------------------------
        self.update_ground_temperatures()
        self.update_ALT()

        # -----------------------------------------
        # Read next perm vars from input files ? NOTE: does not work see the read_input_files()
        # -------------------------------------------
        # Note that read_input_files() is called
        # by initialize() and these values must be
        # used for "update" calls before reading
        # new ones.
        # -------------------------------------------
        if self.time_index > 0:
            self.read_input_files()

        # ----------------------------------------------
        # Write user-specified data to output files ?
        # ----------------------------------------------
        # Components use own self.time_sec by default.
        # -----------------------------------------------
        self.write_output_files()

        # -----------------------------
        # Update internal clock
        # after write_output_files()
        # -----------------------------
        self.update_time(dt)
        self.status = "updated"  # (OpenMI)

    #   update()
    # -------------------------------------------------------------------

    def write_output_files(self, time_seconds=None):
        # -----------------------------------------
        # Allows time to be passed from a caller
        # -----------------------------------------
        if time_seconds is None:
            time_seconds = self.time_sec
        eodel_time = int(time_seconds)

        # ----------------------------------------
        # Save computed values at sampled times
        # ----------------------------------------
        if model_time % int(self.save_grid_dt) == 0:
            self.save_grids()
        if model_time % int(self.save_pixels_dt) == 0:
            self.save_pixel_values()

        # ----------------------------------------
        # Save computed values at sampled times
        # ----------------------------------------

    ##        if ((self.time_index % self.grid_save_step) == 0):
    ##             self.save_grids()
    ##        if ((self.time_index % self.pixel_save_step) == 0):
    ##             self.save_pixel_values()

    #   write_output_files()
    # -------------------------------------------------------------------
    def save_grids(self):
        # Saves the grid values based on the prescribed ones in cfg file

        # if (self.SAVE_MR_GRIDS):
        #    model_output.add_grid( self, self.T_air, 'T_air', self.time_min )

        if self.SAVE_HS_GRIDS:
            model_output.add_grid(self, self.h_snow, "h_snow", self.time_min)

        # if (self.SAVE_SW_GRIDS):
        #    model_output.add_grid( self, self.Tps, 'Tps', self.time_min )

        # if (self.SAVE_CC_GRIDS):
        #    model_output.add_grid( self, self.Zal, 'Zal', self.time_min )

    #   save_grids()
    # -------------------------------------------------------------------
    def save_pixel_values(self):
        # IDs  = self.outlet_IDs  : not sure what is this
        time = self.time_min  ###

        # if (self.SAVE_MR_PIXELS):
        #    model_output.add_values_at_IDs( self, time, self.SM, 'mr', IDs )

        if self.SAVE_HS_PIXELS:
            model_output.add_values_at_IDs(self, time, self.h_snow, "h_snow", IDs)

        # if (self.SAVE_SW_PIXELS):
        #    model_output.add_values_at_IDs( self, time, self.h_swe, 'sw', IDs )

        # if (self.SAVE_CC_PIXELS):
        #    model_output.add_values_at_IDs( self, time, self.Ecc, 'cc', IDs )

    #   save_pixel_values()
    # -------------------------------------------------------------------
    def finalize(self):
        self.status = "finalizing"  # (OpenMI)
        self.close_input_files()  ##  TopoFlow input "data streams"
        self.close_output_files()
        self.status = "finalized"  # (OpenMI)

        self.print_final_report(comp_name="Permafrost component")

        # ---------------------------
        # Release all of the ports
        # ----------------------------------------
        # Make this call in "finalize()" method
        # of the component's CCA Imple file
        # ----------------------------------------
        # self.release_cca_ports( port_names, d_services )

    #   finalize()
    # -------------------------------------------------------------------
    def close_output_files(self):
        # Note: these are component specific, and so are commented out here
        # if (self.SAVE_MR_GRIDS): model_output.close_gs_file( self, 'mr')
        # if (self.SAVE_HS_GRIDS): model_output.close_gs_file( self, 'hs')
        # if (self.SAVE_SW_GRIDS): model_output.close_gs_file( self, 'sw')
        # if (self.SAVE_CC_GRIDS): model_output.close_gs_file( self, 'cc')
        # -----------------------------------------------------------------
        # if (self.SAVE_MR_PIXELS): model_output.close_ts_file( self, 'mr')
        # if (self.SAVE_HS_PIXELS): model_output.close_ts_file( self, 'hs')
        # if (self.SAVE_SW_PIXELS): model_output.close_ts_file( self, 'sw')
        # if (self.SAVE_CC_PIXELS): model_output.close_ts_file( self, 'cc')

        # Since there are no commands, we need a 'pass' statement here
        pass

    # -------------------------------------------------------------------
    def close_input_files(self):
        print("ERROR: close_input_files() for Snow component")
        print("       has not been implemented.")

    #   close_input_files()
    # -------------------------------------------------------------------

    # ----------------------------------------
    # CRU geotiff file interpretation routines
    # ----------------------------------------

    """

    WARNING: These don't work "automagically" because GDAL is busted in Python

    def get_cru_indexes_from_lon_lat(self, lon, lat, month, year):
        # Inspired by:
        #   http://gis.stackexchange.com/questions/122335/using-gdals-getprojection-information-to-make-a-coordinate-conversion-in-pyproj
        #   http://geoinformaticstutorial.blogspot.com/2012/09/
        #          reading-raster-data-with-python-and-gdal.html

        temp_filename = self.get_temperature_tiff_filename(year, month)

        ds = gdal.Open(temp_filename, GA_ReadOnly)
        tiff_proj_wkt = ds.GetProjection()
        proj_converter = osr.SpatialReference()
        proj_converter.ImportFromWkt(tiff_proj_wkt)
        tiff_Proj4_string = proj_converter.ExportToProj4()
        p1 = Proj(tiff_Proj4_string)

        # (xm, ym) is the point on the projected grid in meters
        (xm, ym) = p1(lon, lat)

        xdim = ds.RasterXSize
        ydim = ds.RasterYSize
        geotransform = ds.GetGeoTransform()

        # (originX, originY) is the upper left corner of the grid
        #   Note: this is *not* the center of the UL gridcell
        #   Note: in grid coordinates, the UL corner is at (-0.5, -0.5)
        originX = geotransform[0]
        originY = geotransform[3]
        pixelWidth = geotransform[1]
        pixelHeight = geotransform[5]  # Note this is negative for cru

        # (zeroX, zeroY) is the center of the UL gridcell
        zeroX = originX + 0.5*pixelWidth
        zeroY = originY + 0.5*pixelHeight

        # x and y are floating points
        x = (xm - zeroX)/pixelWidth
        y = (ym - zeroY)/pixelHeight

        # i and j are the rounded index values
        i = int(round(x))
        j = int(round(y))

        # Ensure that point is on the grid
        assert(x>=-0.501)
        assert(y>=-0.501)
        assert(x<=xdim-1+0.501)
        assert(y<=ydim-1+0.501)

        # Ensure that the indexes are valid
        i = max(i, 0)
        j = max(j, 0)
        i = min(i, xdim-1)
        j = min(j, ydim-1)

        return (i, j)


    def get_lon_lat_from_cru_indexes(self, i, j, month, year):
        # Based on:
        #  http://gis.stackexchange.com/questions/122335/using-gdals-getprojection-information-to-make-a-coordinate-conversion-in-pyproj

        # Note: i, j can be floating point

        temp_filename = self.get_temperature_tiff_filename(year, month)

        ds = gdal.Open(temp_filename, GA_ReadOnly)
        tiff_proj_wkt = ds.GetProjection()
        proj_converter = osr.SpatialReference()
        proj_converter.ImportFromWkt(tiff_proj_wkt)
        tiff_Proj4_string = proj_converter.ExportToProj4()
        p1 = Proj(tiff_Proj4_string)

        # Following:
        #   http://geoinformaticstutorial.blogspot.com/2012/09/
        #          reading-raster-data-with-python-and-gdal.html

        xdim = ds.RasterXSize
        ydim = ds.RasterYSize
        geotransform = ds.GetGeoTransform()

        # (originX, originY) is the upper left corner of the grid
        #   Note: this is *not* the center of the UL gridcell
        originX = geotransform[0]
        originY = geotransform[3]
        pixelWidth = geotransform[1]
        pixelHeight = geotransform[5]  # Note this is negative for cru

        # (zeroX, zeroY) is the center of the UL gridcell
        zeroX = originX + 0.5*pixelWidth
        zeroY = originY + 0.5*pixelHeight

        # (xm, ym) is the point on the projected grid in meters
        xm = zeroX + i * pixelWidth
        ym = zeroY + j * pixelHeight
        (lon, lat) = p1(xm, ym, inverse=True)

        return (lon, lat)


    # ----------------------------------------
    # Read precipitation routines from geotiff
    # ----------------------------------------

    def get_temperature_tiff_filename(self, year, month, datadir=None):
        # Generate the name of the CRU tiff file

        datadir = datadir or os.path.join(os.environ.get('PERMAMODEL_DATADIR', 
                                                         '/data'), 'tas')

        filename = "%s/tas_mean_C_cru_TS31_%02d_%4d.tif" % \
                (datadir, month, year)
        if not os.path.isfile(filename):
            print("Warning: this temperature tiff file does not exist: %s" %\
                  filename)
            exit(-1)

        return filename

    #   get_temperature_tiff_filename(year, month, [datadir])

    def get_temperature_from_cru_indexes(self, i, j, m, y):
        # Inputs are:
        #    (i, j):  the location on the grid
        #    m:       the month
        #    y:       the year
        temp_filename = self.get_temperature_tiff_filename(y, m)
        ds = gdal.Open(temp_filename, GA_ReadOnly)

        # Examine the tiff metadata, but it's just the Area or Point info
        # print("Tiff Metadata:\n%s" % ds.GetMetadata())

        # Verify that we are checking a point in the grid
        xdim = ds.RasterXSize
        ydim = ds.RasterYSize

        assert(i<xdim)
        assert(j<ydim)

        band = ds.GetRasterBand(1)

        temperatures = band.ReadAsArray(0, 0, xdim,
                           ydim).astype(gdal.GetDataTypeName(band.DataType))

        return temperatures[j][i]


    def get_temperature_from_cru(self, lon, lat, month, year):
        (i, j) = self.get_cru_indexes_from_lon_lat(lon, lat, month, year)
        temperature = self.get_temperature_from_cru_indexes(i, j, month, year)
        return temperature


    # ----------------------------------------
    # Read precipitation routines from geotiff
    # ----------------------------------------

    def get_precipitation_tiff_filename(self, year, month, datadir=None):
        datadir = datadir or os.path.join(os.environ.get('PERMAMODEL_DATADIR', 
                                                         '/data'), 'pr_AK_CRU')
        # Generate the name of the CRU precipitation tiff file

        filename = "%s/pr_total_mm_cru_TS31_%02d_%4d.tif" % \
                (datadir, month, year)
        if not os.path.isfile(filename):
            print("Warning: this temperature tiff file does not exist: %s" %\
                  filename)
            exit(-1)

        return filename

    #   get_temperature_tiff_filename(year, month, [datadir])

    def get_precipitation_from_cru_indexes(self, i, j, m, y):
        # Inputs are:
        #    (i, j):  the location on the grid
        #    m:       the month
        #    y:       the year
        prec_filename = self.get_precipitation_tiff_filename(y, m)
        ds = gdal.Open(prec_filename, GA_ReadOnly)

        # Examine the tiff metadata, but it's just the Area or Point info
        # print("Tiff Metadata:\n%s" % ds.GetMetadata())

        # Verify that we are checking a point in the grid
        xdim = ds.RasterXSize
        ydim = ds.RasterYSize

        assert(i<xdim)
        assert(j<ydim)

        band = ds.GetRasterBand(1)

        precipitations = band.ReadAsArray(0, 0, xdim,
                           ydim).astype(gdal.GetDataTypeName(band.DataType))

        return precipitations[j][i]


    def get_precipitation_from_cru(self, lon, lat, month, year):
        (i, j) = self.get_cru_indexes_from_lon_lat(lon, lat, month, year)
        temperature = self.get_precipitation_from_cru_indexes(i, j, month, year)
        return temperature

    """

    def get_permafrost_data_directory(self):
        """Note: this can now be:

        from .. import permamodel_directory
        """
        datadir = os.environ.get("PERMAMODEL_DATADIR")
        if datadir is None:
            print("Please set the environment variable PERMAMODEL_DATADIR")
            print("  to the location of the Permafrost Model data")
            print("  e.g. in ~/.bashrc:")
            print("     export PERMAMODEL_DATADIR=/data/permafrost_data/")
            assert False
        return datadir
