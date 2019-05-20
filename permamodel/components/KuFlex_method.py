# -*- coding: utf-8 -*-
"""  Kudryavtsev Model code adapted for the BMI version developed for the Topoflow model

     Author: Kang Wang, 03/29/2016
     Modified: Elchin Jafarov, 03/29/2016

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

Input:
    (1) Location:
        input_lat: Latitude
        input_lon: Longitude

    (2) Climate :
        Ta  : Mean annual air temperature (C)
        Aa  : Amplitude of air temperature (C)
        Hsn : Winter-Averaged Snow Depth (m)
        Rsn : Snow Density (kg/m3)
        vwc : Volumetric Water Content (m3 / m3)

    (3) Vegetation:
        Hvgf: Height of vegetation in frozen period (m)
        Hvgt: Height of vegetation in thawed period (m)
        Dvf : Thermal diffusivity of vegetation in frozen period (m2 s)
        Dvt : Thermal diffusivity of vegetation in thawed period (m2 s)

Output:
        1) Mean annual temperature on the top of permafrost (C)
        2) Active Layer Thickness (m)

References:

    Anisimov, O. A., Shiklomanov, N. I., & Nelson, F. E. (1997).
        Global warming and active-layer thickness: results from transient general circulation models.
        Global and Planetary Change, 15(3), 61-77.
    Romanovsky, V. E., & Osterkamp, T. E. (1997).
        Thawing of the active layer on the coastal plain of the Alaskan Arctic.
        Permafrost and Periglacial processes, 8(1), 1-22.
    Sazonova, T. S., & Romanovsky, V. E. (2003).
        A model for regional‐scale estimation of temporal and spatial variability of active layer thickness and mean annual ground temperatures.
        Permafrost and Periglacial Processes, 14(2), 125-139.
    Sturm, M., Holmgren, J., König, M., & Morris, K. (1997).
        The thermal conductivity of seasonal snow. Journal of Glaciology, 43(143), 26-41.
    Ling, F., & Zhang, T. (2004).
        A numerical model for surface energy balance and thermal regime of the active layer and permafrost containing unfrozen water.
        Cold Regions Science and Technology, 38(1), 1-15.
    Wieder, W.R., J. Boehnert, G.B. Bonan, and M. Langseth. (2014).
        Regridded Harmonized World Soil Database v1.2. Data set.
        Available on-line [http://daac.ornl.gov] from Oak Ridge National Laboratory Distributed Active Archive Center, Oak Ridge, Tennessee, USA.  http://dx.doi.org/10.3334/ORNLDAAC/1247  .

"""
from __future__ import print_function

import os
import numpy as np
from permamodel.utils import model_input
from permamodel.components import perma_base
from netCDF4 import Dataset
from .. import data_directory
# from permamodel.tests import examples_directory


class KuFlex_method( perma_base.PermafrostComponent ):

    #   get_var_units()
    #-------------------------------------------------------------------
    def check_input_types(self):

        #--------------------------------------------------
        # Notes: rho_H2O, Cp_snow, rho_air and Cp_air are
        #        currently always scalars.
        #--------------------------------------------------
        are_scalars = np.array([
                          self.is_scalar('T_air'),
                          self.is_scalar('A_air'),
                          
                          self.is_scalar('h_snow'),
                          self.is_scalar('k_snow'),
                          self.is_scalar('c_snow'),
                          self.is_scalar('rho_snow'),
                          
                          self.is_scalar('lh_soil'),
                          self.is_scalar('kt_soil'),
                          self.is_scalar('kf_soil'),
                          self.is_scalar('ct_soil'),
                          self.is_scalar('cf_soil'),
                          
                          self.is_scalar('Hvgf'),
                          self.is_scalar('Hvgt'),
                          self.is_scalar('Dvf'),
                          self.is_scalar('Dvt') ])

        self.ALL_SCALARS = np.all(are_scalars)

    #   check_input_types()
    #-------------------------------------------------------------------
    def open_input_files(self):
        # this function will work only if filename is not empty
        
        # Variables related to air:

        self.T_air_file       = self.in_directory + self.T_air_file
        self.A_air_file       = self.in_directory + self.A_air_file
        
        # Variables related to snow:
        
        self.h_snow_file      = self.in_directory + self.h_snow_file
        self.k_snow_file      = self.in_directory + self.k_snow_file
        self.c_snow_file      = self.in_directory + self.c_snow_file
        self.rho_snow_file    = self.in_directory + self.rho_snow_file
        
        # Variables related to vegetation:
        
        self.Hvgf_file        = self.in_directory + self.Hvgf_file
        self.Hvgt_file        = self.in_directory + self.Hvgt_file
        self.Dvf_file         = self.in_directory + self.Dvf_file
        self.Dvt_file         = self.in_directory + self.Dvt_file
        
        # Variables related to soil thermal:
        
        self.lh_soil_file     = self.in_directory + self.lh_soil_file
        self.kt_soil_file     = self.in_directory + self.kt_soil_file
        self.kf_soil_file     = self.in_directory + self.kf_soil_file
        self.ct_soil_file     = self.in_directory + self.ct_soil_file
        self.cf_soil_file     = self.in_directory + self.cf_soil_file
        
        # Variables for outputs:

        self.Outfile         = self.out_directory + self.outfile
        
        # File Units for each:
         
        self.T_air_unit       = self.open_file_KU(self.T_air_type,  self.T_air_file)
        self.A_air_unit       = self.open_file_KU(self.A_air_type,  self.A_air_file)
        
        self.h_snow_unit      = self.open_file_KU(self.h_snow_type,  self.h_snow_file)
        self.rho_snow_unit    = self.open_file_KU(self.rho_snow_type,  self.rho_snow_file)
        self.k_snow_unit      = self.open_file_KU(self.k_snow_type,  self.k_snow_file)
        self.c_snow_unit      = self.open_file_KU(self.c_snow_type,  self.c_snow_file)
        
        self.Hvgf_unit        = self.open_file_KU(self.Hvgf_type,  self.Hvgf_file)
        self.Hvgt_unit        = self.open_file_KU(self.Hvgt_type,  self.Hvgt_file)
        self.Dvf_unit         = self.open_file_KU(self.Dvf_type,  self.Dvf_file)
        self.Dvt_unit         = self.open_file_KU(self.Dvt_type,  self.Dvt_file)
        
        self.lh_soil_unit     = self.open_file_KU(self.lh_soil_type,  self.lh_soil_file)
        self.kt_soil_unit     = self.open_file_KU(self.kt_soil_type,  self.kt_soil_file)
        self.kf_soil_unit     = self.open_file_KU(self.kf_soil_type,  self.kf_soil_file)
        self.ct_soil_unit     = self.open_file_KU(self.ct_soil_type,  self.ct_soil_file)
        self.cf_soil_unit     = self.open_file_KU(self.cf_soil_type,  self.cf_soil_file)
        
    def open_output_files(self):
        
        n_lat = np.size(self.lat)
        n_lon = np.size(self.lon)

        self.out_fid = Dataset(self.Outfile+'.nc', 'w', format='NETCDF4')
        
        # ==== Latitude ====

        self.out_fid.createDimension('lat', n_lat) # Create Dimension
        lats = self.out_fid.createVariable('lat',np.dtype('float32').char,('lat',))
        lats.units = 'degrees_north'
        lats.standard_name = 'latitude'
        lats.long_name = 'latitude'
        lats.axis = 'Y'
        lats[:] = self.lat
        
        # ==== Longitude ====

        self.out_fid.createDimension('lon', n_lon) # Create Dimension
        lons = self.out_fid.createVariable('lon',np.dtype('float32').char,('lon',))
        lons.units = 'degrees_east'
        lons.standard_name = 'longitude'
        lons.long_name = 'longitude'
        lons.axis = 'X'
        lons[:] = self.lon
        
        # ==== Time ====

        self.out_fid.createDimension('time', None) # Create Dimension
        self.out_timeid = self.out_fid.createVariable('time',np.dtype('float32').char,('time',))
        self.out_timeid.units = 'Year'
        self.out_timeid.axis = 'Z'
             
        if self.SAVE_ALT_GRIDS:
            
            units = 'm'
            long_name = 'Active Layer Thickness'
            
            # ==== Data ====
            self.alt_out_varid = self.out_fid.createVariable('ALT',np.dtype('float32').char,
                                                   ('time','lat','lon'),
                                                   fill_value = -999)
            self.alt_out_varid.units = units
            self.alt_out_varid.long_name = long_name     

        if self.SAVE_TPS_GRIDS:
            
            units = 'deg_C'
            long_name = 'Mean Annual Temperature at Permafrost Surface'
                               
            # ==== Data ====
            self.tps_out_varid = self.out_fid.createVariable('Tps',np.dtype('float32').char,
                                                   ('time','lat','lon'),
                                                   fill_value = -999)
            self.tps_out_varid.units = units
            self.tps_out_varid.long_name = long_name 
  
    #   open_input_files()
    #-------------------------------------------------------------------
    def read_input_files(self):

        #rti = self.rti # has a problem with loading rti: do not know where its been initialized
        
        size_of_inputs = np.ones(15)
        size_of_x      = np.ones(15)

        #-------------------------------------------------------
        # All grids are assumed to have a data type of Float32.
        #-------------------------------------------------------    
 
        #--------------  AIR              
        T_air = self.read_next_modified_KU(self.T_air_unit,self.T_air_type)   
        
        if (self.T_air_type.lower() == 'grid'):
            self.T_air = T_air
            n_T_air  = np.size(T_air)
            size_of_inputs[0] = n_T_air
            size_of_x[0]      = np.shape(T_air)[0]
        elif (self.T_air_type.lower() == 'time_series'):
            self.T_air = T_air
            n_T_air = 1
        else:
            n_T_air = 1
            
        A_air = self.read_next_modified_KU(self.A_air_unit, self.A_air_type)
        if (self.A_air_type.lower() == 'grid'): 
            self.A_air = A_air
            n_A_air  = np.size(A_air)
            size_of_inputs[1] = n_A_air
            size_of_x[1]      = np.shape(A_air)[0]
        elif (self.A_air_type.lower() == 'time_series'):
            self.A_air = A_air
            n_A_air = 1
        else:
            n_A_air = 1
            
        #--------------  SNOW                          
        h_snow = self.read_next_modified_KU(self.h_snow_unit, self.h_snow_type)
        if (h_snow is not None): 
            self.h_snow = h_snow
            n_h_snow  = np.size(h_snow)
            size_of_inputs[2] = n_h_snow
            size_of_x[2]      = np.shape(h_snow)[0]
        elif (self.h_snow_type.lower() == 'time_series'):
            self.h_snow = h_snow
            n_h_snow = 1
        else:
            n_h_snow = 1
            
        k_snow = self.read_next_modified_KU(self.k_snow_unit, self.k_snow_type)
        if (self.k_snow_type.lower() == 'grid'): 
            self.k_snow = k_snow
            n_k_snow  = np.size(k_snow)
            size_of_inputs[3] = n_k_snow
            size_of_x[3]      = np.shape(k_snow)[0]
        elif (self.k_snow_type.lower() == 'time_series'):
            self.k_snow = k_snow
            n_k_snow = 1
        else:
            n_k_snow = 1
            
        c_snow = self.read_next_modified_KU(self.c_snow_unit, self.c_snow_type)
        if (self.c_snow_type.lower() == 'grid'): 
            self.c_snow = c_snow    
            n_c_snow  = np.size(c_snow)
            size_of_inputs[4] = n_c_snow
            size_of_x[4]      = np.shape(c_snow)[0]
        elif (self.c_snow_type.lower() == 'time_series'):
            self.c_snow = c_snow
            n_c_snow = 1
        else:
            n_c_snow = 1
            
        rho_snow = self.read_next_modified_KU(self.rho_snow_unit, self.rho_snow_type)
        if (self.rho_snow_type.lower() == 'grid'): 
            self.rho_snow = rho_snow
            n_rho_snow  = np.size(rho_snow)
            size_of_inputs[5] = n_rho_snow
            size_of_x[5]      = np.shape(rho_snow)[0]
        elif (self.rho_snow_type.lower() == 'time_series'):
            self.rho_snow = rho_snow
            n_rho_snow = 1
        else:
            n_rho_snow = 1
            
        #--------------   VEG                         
        Hvgf = self.read_next_modified_KU(self.Hvgf_unit, self.Hvgf_type)
        if (self.Hvgf_type.lower() == 'grid'): 
            self.Hvgf = Hvgf  
            n_Hvgf  = np.size(Hvgf)
            size_of_inputs[6] = n_Hvgf
            size_of_x[6]      = np.shape(Hvgf)[0]
        elif (self.Hvgf_type.lower() == 'time_series'): 
            self.Hvgf = Hvgf  
            n_Hvgf = 1
        else:
            n_Hvgf = 1
            
        Hvgt = self.read_next_modified_KU(self.Hvgt_unit, self.Hvgt_type)
        if (self.Hvgt_type.lower() == 'grid'): 
            self.Hvgt = Hvgt
            n_Hvgt  = np.size(Hvgt)
            size_of_inputs[7] = n_Hvgt
            size_of_x[7]      = np.shape(Hvgt)[0]
        elif (self.Hvgt_type.lower() == 'time_series'): 
            self.Hvgt = Hvgt
            n_Hvgt = 1
        else:
            n_Hvgt = 1
            
        Dvt = self.read_next_modified_KU(self.Dvt_unit, self.Dvt_type)
        if (self.Dvt_type.lower() == 'grid'): 
            self.Dvt = Dvt
            n_Dvt  = np.size(Dvt)
            size_of_inputs[8] = n_Dvt
            size_of_x[8]      = np.shape(Dvt)[0]
        elif (self.Dvt_type.lower() == 'time_series'): 
            self.Dvt = Dvt
            n_Dvt = 1
        else:
            n_Dvt = 1
            
        Dvf = self.read_next_modified_KU(self.Dvf_unit, self.Dvf_type)
        if (self.Dvf_type.lower() == 'grid'): 
            self.Dvf = Dvf
            n_Dvf  = np.size(Dvf)
            size_of_inputs[9] = n_Dvf
            size_of_x[9]      = np.shape(Dvf)[0]
        elif (self.Dvf_type.lower() == 'time_series'): 
            self.Dvf = Dvf
            n_Dvf = 1
        else:
            n_Dvf = 1
            
        #--------------    SOIL
        lh_soil = self.read_next_modified_KU(self.lh_soil_unit, self.lh_soil_type)
        if (self.lh_soil_type.lower() == 'grid'): 
            self.lh_soil = lh_soil
            n_lh_soil  = np.size(lh_soil)
            size_of_inputs[10] = n_lh_soil
            size_of_x[10]      = np.shape(lh_soil)[0]
        elif (self.lh_soil_type.lower() == 'time_series'): 
            self.lh_soil = lh_soil   
            n_lh_soil = 1
        else:
            n_lh_soil = 1            
            
        kt_soil = self.read_next_modified_KU(self.kt_soil_unit, self.kt_soil_type)
        if (self.kt_soil_type.lower() == 'grid'): 
            self.kt_soil = kt_soil
            n_kt_soil  = np.size(kt_soil)
            size_of_inputs[11] = n_kt_soil
            size_of_x[11]      = np.shape(kt_soil)[0]
        if (self.kt_soil_type.lower() == 'time_series'): 
            self.kt_soil = kt_soil
            n_kt_soil = 1
        else:
            n_kt_soil = 1 
            
        kf_soil = self.read_next_modified_KU(self.kf_soil_unit, self.kf_soil_type)
        if (self.kf_soil_type.lower() == 'grid'): 
            self.kf_soil = kf_soil
            n_kf_soil  = np.size(kf_soil)
            size_of_inputs[12] = n_kf_soil
            size_of_x[12]      = np.shape(kf_soil)[0]
        elif (self.kf_soil_type.lower() == 'time_series'): 
            self.kf_soil = kf_soil
            n_kf_soil = 1
        else:
            n_kf_soil = 1             

        ct_soil = self.read_next_modified_KU(self.ct_soil_unit, self.ct_soil_type)
        if (self.ct_soil_type.lower() == 'grid'): 
            self.ct_soil = ct_soil
            n_ct_soil  = np.size(ct_soil)
            size_of_inputs[13] = n_ct_soil
            size_of_x[13]      = np.shape(ct_soil)[0]
        if (self.ct_soil_type.lower() == 'time_series'): 
            self.ct_soil = ct_soil
            n_ct_soil = 1 
        else:
            n_ct_soil = 1 
            
        cf_soil = self.read_next_modified_KU(self.cf_soil_unit, self.cf_soil_type)
        if (self.cf_soil_type.lower() == 'grid'): 
            self.cf_soil = cf_soil
            n_cf_soil  = np.size(cf_soil)
            size_of_inputs[14] = n_cf_soil
            size_of_x[14]      = np.shape(cf_soil)[0]
        elif (self.cf_soil_type.lower() == 'time_series'): 
            self.cf_soil = cf_soil
            n_cf_soil = 1
        else:
            n_cf_soil = 1    
            
        ## check input shapes:
        
        self.n_grids = np.int(np.max(size_of_inputs))  
        
        uniq_grids_size = np.unique(size_of_inputs)
                
        assert np.size(uniq_grids_size)<=2, 'Error in inputs shapes'
        assert np.size(np.unique(size_of_x)), 'Error in inputs shapes'
        
        self.grid_shape = [np.int(np.max(size_of_x)), np.int(self.n_grids / np.max(size_of_x))]
                
        ## repeat scalers to same grid shape:

        if n_T_air ==1:
            self.T_air = self.T_air + np.zeros(self.grid_shape)
        if n_A_air ==1:
            self.A_air = self.A_air + np.zeros(self.grid_shape)

        if n_h_snow ==1:
            self.h_snow = self.h_snow + np.zeros(self.grid_shape)            
        if n_k_snow ==1:
            self.k_snow = self.k_snow + np.zeros(self.grid_shape) 
        if n_c_snow ==1:
            self.c_snow = self.c_snow + np.zeros(self.grid_shape)         
        if n_rho_snow ==1:
            self.rho_snow = self.rho_snow + np.zeros(self.grid_shape)
            
        if n_Hvgf ==1:
            self.Hvgf = self.Hvgf  + np.zeros(self.grid_shape)
        if n_Hvgt ==1:
            self.Hvgt = self.Hvgt  + np.zeros(self.grid_shape)
        if n_Dvf ==1:
            self.Dvf = self.Dvf   + np.zeros(self.grid_shape)
        if n_Dvt ==1:
            self.Dvt = self.Dvt   + np.zeros(self.grid_shape)
            
        if n_kf_soil ==1:
            self.kf_soil = self.kf_soil  + np.zeros(self.grid_shape)
        if n_kt_soil ==1:
            self.kt_soil = self.kt_soil  + np.zeros(self.grid_shape)
        if n_lh_soil ==1:
            self.lh_soil = self.lh_soil  + np.zeros(self.grid_shape)
        if n_cf_soil ==1:
            self.cf_soil = self.cf_soil   + np.zeros(self.grid_shape)
        if n_ct_soil ==1:
            self.ct_soil = self.ct_soil  + np.zeros(self.grid_shape)

    #   read_input_files()
    #-------------------------------------------------------------------

    def update_soil_heat_capacity(self):

        self.Ct = self.ct_soil + 0. 
        self.Cf = self.cf_soil + 0.
        
    #   update_soil_heat_capacity()
    #-------------------------------------------------------------------

    def update_soil_thermal_conductivity(self):
        
        self.Kt = self.kt_soil + 0.;
        self.Kf = self.kf_soil + 0.;
        
    #   update_soil_thermal_conductivity()
    #-------------------------------------------------------------------
    def update_snow_thermal_properties(self):

        self.Ksn = self.k_snow + 0.; # Unit: (W m-1 C-1)
        self.Csn = self.c_snow + 0.; # Unit: J m-3 C-1

    #   update_ssnow_thermal_properties()
    #-------------------------------------------------------------------
    def update_TOP_temperatures(self):

        #---------------------------------------------------------
        #   1.  Estimating vegetation effect
        #       deta_Tsn -- eq-7, Anisimov et al. 1997
        #       deta_Asn -- eq-2, Sazonova et al., 2003
        #       Tvg -- mean annual temperature Page-129, Sazonova et al., 2003
        #       Avg -- amplitude bellow snow OR top of vegetation
        #--------------------------------------------------

        tao = np.zeros(self.grid_shape) + self.sec_per_year;        
        
        K_diffusivity = self.Ksn/(self.rho_snow*self.Csn)
        
        temp = np.exp(-1.0*self.h_snow*np.sqrt(np.pi/(tao*K_diffusivity)))
        deta_Tsn = self.A_air*(1.0 - temp);
        deta_Asn = deta_Tsn*2.0/np.pi;

        Tvg = self.T_air + deta_Tsn;
        Avg = self.A_air - deta_Asn;
        
        self.deta_Tsn = deta_Tsn + 0.;
        self.deta_Asn = deta_Asn + 0.;

        #---------------------------------------------------------
        #   2.  Estimating Snow Effects
        #       deta_A1 -- winter vegetation thermal effects: eq-10, Anisimov et al. 1997
        #       deta_A2 -- summer vegetation thermal effects: eq-11, Anisimov et al. 1997
        #       deta_Av -- Effects of vegetation on seasonal amplitude of temperature, eq-8
        #       deta_Tv -- Effects of vegetation on an annual mean temperature, eq-9
        #       Tgs, Ags -- mean annual gs temperature and amplitude eq-13,14 Sazonova et al., 2003
        #--------------------------------------------------
        temp = 1.- np.exp(-1.*self.Hvgf*np.sqrt(np.pi/(self.Dvf*2.*self.tao1)))
        deta_A1 = (Avg - Tvg) * temp;

        temp = 1.- np.exp(-1.*self.Hvgt*np.sqrt(np.pi/(self.Dvt*2.*self.tao2)))
        deta_A2 = (Avg  + Tvg) * temp;

        deta_Av = (deta_A1*self.tao1+deta_A2*self.tao2) / tao;

        deta_Tv = (deta_A1*self.tao1-deta_A2*self.tao2) / tao * (2. / np.pi)

        Tgs = Tvg + deta_Tv;
        Ags = Avg - deta_Av;
        
        self.deta_Tv = deta_Tv + 0.;
        self.deta_Av = deta_Av + 0.;

        #---------------------------------------------------------
        #   3.  Calculates Tps_Numerator;
        #       eq-14, Anisimov et al. 1997
        #--------------------------------------------------

        Tps_numerator = 0.5*Tgs*(self.Kf+self.Kt)\
                            +(Ags*(self.Kt-self.Kf)/np.pi\
                            *(Tgs/Ags*np.arcsin(Tgs/Ags)\
                            +np.sqrt(1.-(np.pi**2.0/Ags**2.0))));

        #---------------------------------------------------------
        #   4.  Calculates temperature at the top of permafrost
        #       Tps -- eq-14 cont., Anisimov et al. 1997
        #--------------------------------------------------

        K_star = self.Kf
            
        K_star[np.where(Tps_numerator>0.0)] = self.Kt[np.where(Tps_numerator>0.0)];
       
        self.Tgs=Tgs
        self.Ags=Ags
        self.Tps_numerator=Tps_numerator

        self.Tps = self.Tps_numerator/K_star

    #   update_TOP_temperatures()
    #-------------------------------------------------------------------
    def update_ALT(self):

        #---------------------------------------------------------
        #       Calculates active layer thickness
        #       Aps  -- eq-4, Romanovsky et al. 1997
        #       Zs -- eq-5, Romanovsky et al. 1997
        #       Zal -- eq-3, Romanovsky et al. 1997
        #--------------------------------------------------

        tao = np.zeros(self.grid_shape) + self.sec_per_year;

        K = self.Kt
        C = self.Ct       
            
        K[np.where(self.Tps_numerator>0.0)] = self.Kf[np.where(self.Tps_numerator>0.0)]
        C[np.where(self.Tps_numerator>0.0)] = self.Cf[np.where(self.Tps_numerator>0.0)]
        
        Aps = np.zeros(self.grid_shape) + np.nan
        Zc  = np.zeros(self.grid_shape) + np.nan
        Zal = np.zeros(self.grid_shape) + np.nan

        idx_need = np.where((self.Tgs + self.Ags >0) & (self.Tgs - self.Ags <0))

        Aps[idx_need] = (self.Ags[idx_need] - abs(self.Tps[idx_need]))/np.log((self.Ags[idx_need]+self.L[idx_need]/(2.*C[idx_need])) / \
                    (abs(self.Tps[idx_need])+self.L[idx_need]/(2.*C[idx_need]))) - self.L[idx_need]/(2.*C[idx_need]);

        Zc[idx_need] = (2.*(self.Ags[idx_need] - abs(self.Tps[idx_need]))*np.sqrt((K[idx_need]*tao[idx_need]*C[idx_need])/np.pi)) / \
                    (2.*Aps[idx_need]*C[idx_need] + self.L[idx_need]);
                    
        Zal[idx_need] = (2.*(self.Ags[idx_need] - abs(self.Tps[idx_need]))*np.sqrt(K[idx_need]*tao[idx_need]*C[idx_need]/np.pi) \
                +(2.*Aps[idx_need]*C[idx_need]*Zc[idx_need]+self.L[idx_need]*Zc[idx_need])*self.L[idx_need]*np.sqrt(K[idx_need]*tao[idx_need]/(np.pi*C[idx_need])) \
                /(2.*Aps[idx_need]*C[idx_need]*Zc[idx_need] + self.L[idx_need]*Zc[idx_need] +\
                  (2.*Aps[idx_need]*C[idx_need]+self.L[idx_need])*\
                  np.sqrt(K[idx_need]*tao[idx_need]/(np.pi*C[idx_need])))) \
                /(2.*Aps[idx_need]*C[idx_need]+ self.L[idx_need]);

        Zal[np.where(self.Tps_numerator>0.0)] = np.nan # Seasonal Frozen Ground
        Zal[Zal <=0 ]                = 0

        self.Aps = Aps;
        self.Zc  = Zc;  
        self.Zal = Zal;                
        
    #   update_ALT()
    #-------------------------------------------------------------------
    def update_ground_temperatures(self):
        # in this method there is only one output the temperature at the top of permafrost
        # TTOP
        self.update_soil_heat_capacity()
        self.update_soil_thermal_conductivity()
        self.update_snow_thermal_properties()
        
        tao = np.zeros(self.grid_shape) + self.sec_per_year;
        self.tao1 = tao * 0.
        self.tao2 = tao * 0.

        # Update mean temperatures for warmes and coldest seasons similar to Nelson & Outcalt 87
        # Cold and Warm Season, Page-129, Sazonova, 2003
        
        idx_unthaw   = np.where(self.T_air + self.A_air <=0)
        idx_unfrozen = np.where(self.T_air - self.A_air >=0)
        
        idx_others   = np.where((self.T_air + self.A_air >0) & (self.T_air - self.A_air <0))

        self.tao1[idx_unthaw] = self.sec_per_year +0.
        self.tao2[idx_unthaw] = 0.
        
        self.tao1[idx_unfrozen] = 0.
        self.tao2[idx_unfrozen] = self.sec_per_year +0.

        self.tao1[idx_others] = tao[idx_others]*(0.5 - 1./np.pi*np.arcsin(self.T_air[idx_others]/self.A_air[idx_others]));
        self.tao2[idx_others] = tao[idx_others] - self.tao1[idx_others];

        self.update_TOP_temperatures()

    #   update_ground_temperatures()
    #-------------------------------------------------------------------
    def close_input_files(self):

        if (self.T_air_type     != 'Scalar'): self.T_air_unit.close()
        if (self.A_air_type     != 'Scalar'): self.A_air_unit.close()
        
        if (self.h_snow_type    != 'Scalar'): self.h_snow_unit.close()
        if (self.rho_snow_type  != 'Scalar'): self.rho_snow_unit.close()
        if (self.k_snow_type    != 'Scalar'): self.k_snow_unit.close()
        if (self.c_snow_type    != 'Scalar'): self.c_snow_unit.close()
        
        if (self.lh_soil_type   != 'Scalar'): self.lh_soil_unit.close()
        if (self.kt_soil_type   != 'Scalar'): self.kt_soil_unit.close()
        if (self.kf_soil_type   != 'Scalar'): self.kf_soil_unit.close()
        if (self.ct_soil_type   != 'Scalar'): self.ct_soil_unit.close()
        if (self.cf_soil_type   != 'Scalar'): self.cf_soil_unit.close()
        
        if (self.Hvgf_type      != 'Scalar'): self.Hvgf_unit.close()
        if (self.Hvgt_type      != 'Scalar'): self.Hvgt_unit.close()
        if (self.Dvf_type       != 'Scalar'): self.Dvf_unit.close()
        if (self.Dvt_type       != 'Scalar'): self.Dvt_unit.close()
        
        
    def close_output_files(self):
        if self.SAVE_ALT_GRIDS or self.SAVE_TPS_GRIDS:
            self.out_fid.close()

    #   close_input_files()
    #-------------------------------------------------------------------
        
    def initialize(self, cfg_file=None, mode="nondriver",
                   SILENT=False):
        
        #---------------------------------------------------------
        # Notes:  Need to make sure than h_swe matches h_snow ?
        #         User may have entered incompatible values.
        #---------------------------------------------------------
        # (3/14/07) If the Energy Balance method is used for ET,
        # then we must initialize and track snow depth even if
        # there is no snowmelt method because the snow depth
        # affects the ET rate.  Otherwise, return to caller.
        #---------------------------------------------------------
        if not(SILENT):
            print(' ')
            print('KuFlex model component: Initializing...')

        self.status     = 'initializing'  # (OpenMI 2.0 convention)
        self.mode       = mode

        # Set the cfg file if it exists, otherwise, provide a default.
        if cfg_file is None:
            cfg_file = "permamodel/examples/KuFlex_method.cfg"
        self.cfg_file = cfg_file

#            if os.path.isfile(cfg_file):
#                print("Default config file exists: %s" % cfg_file)
#            else:
#                print("Default config file does not exist: ")
#                print("  %s" % cfg_file)
#                raise(ValueError(
#                    "Default frostnumber config file %s does not exist" %\
#                    cfg_file))
        # Initialize the output variables (internal names)
        self.Tgs = np.float32(-999.99)
        self.Ags = np.float32(-999.99)
        self.Tps = np.float32(-999.99)
        self.Zal = np.float32(-999.99)
        self.time = 0

        #-----------------------------------------------
        # Load component parameters from a config file
        #-----------------------------------------------
        self.set_constants()
        self.initialize_config_vars()
        # At this stage we are going to ignore read_grid_info b/c
        # we do not have rti file associated with our model
        # we also skipping the basin_vars which calls the outlets
        #self.read_grid_info()
        #self.initialize_basin_vars()
        
        self.initialize_time_vars()

        # Initialize the year to the start year
        #  or to zero if it doesn't exist
        try:
            self.year = self.start_year
        except AttributeError:
            self.year = 0
            self.start_year = 0
            self.end_year = 0

        # Ensure that the end_year is not before the start_year
        # If no end_year is given,
        #   it is assumed that this will run for one year
        #   so the end_year is the same as the start_year
        try:
            assert(self.end_year >= self.start_year)
        except AttributeError:
            self.end_year = self.start_year

        #---------------------------------------------
        # Open input files needed to initialize vars
        #---------------------------------------------
        self.open_input_files()
        self.read_input_files()
        
        self.L=self.lh_soil + 0.        
        
        self.lat = np.arange(self.grid_shape[0])
        self.lon = np.arange(self.grid_shape[1])
        
        self.open_output_files()
        
        self.status = 'initialized'
        
    def read_next_modified_KU(self, file_unit, var_type, \
                  dtype='Float32', factor=1.0):
    
        #-------------------------------------------------------
        # (5/7/09) Allow "dtype" to be given using RTI types.
        # (4/21/16) Elchin Jafarov introduced this function b/c
        # he was not sure how to deal with rti in the original function
        # this version foes not have grid choice
        #-------------------------------------------------------
        rti_types = ['BYTE','INTEGER','LONG','FLOAT','DOUBLE']
        if (dtype.upper() in rti_types):
            dtype_map = {'BYTE':'uint8', 'INTEGER':'int16',
                         'LONG':'int32',
                         'FLOAT':'float32', 'DOUBLE':'float64'}
            dtype = dtype_map[dtype]
    
    
        if (var_type.lower() == 'scalar'):
            #-------------------------------------------
            # Scalar value was entered by user already
            #-------------------------------------------
            data = None
            
        elif (var_type.lower() == 'time_series'):
            #----------------------------------------------
            # Time series: Read scalar value from file.
            # File is ASCII text with one value per line.
            #----------------------------------------------
            data = model_input.read_scalar(file_unit, dtype)
            
        elif (var_type.lower() == 'grid'):
            #----------------------------------------------
            # Time series: Read scalar value from file.
            # File is ASCII text with one value per line.
            #----------------------------------------------
#            data = np.loadtxt(file_name)
#            print self.cont
            for var in file_unit.variables.keys():
                if (var != 'time' and var[0:3] !='lat' and var[0:3] != 'lon'):
                    data  = file_unit.variables[var][self.time,:,:]
                                   
        else:
            raise RuntimeError('No match found for "var_type".')
            return None
    
        #---------------------------------------------
        # Multiply by a conversion or scale factor ?
        #---------------------------------------------
        if (factor != 1) and (data is not None):
            data = (data * factor)            
    
        #-----------------------------------------------------
        # Values must usually be read from file as FLOAT32
        # but then need to be returned as FLOAT64. (5/17/12)
        # But numpy.float64( None ) = NaN. (5/18/12)
        #-----------------------------------------------------
        if (data is None):
            return
        else:
            return np.float32( data )
        
     ## def update(self, dt=-1.0, time_seconds=None):
    def update(self, dt=-1.0):

        #----------------------------------------------------------
        # Note: The read_input_files() method is first called by
        #       the initialize() method.  Then, the update()
        #       method is called one or more times, and it calls
        #       other update_*() methods to compute additional
        #       variables using input data that was last read.
        #       Based on this pattern, read_input_files() should
        #       be called at end of update() method as done here.
        #       If the input files don't contain any additional
        #       data, the last data read persists by default.
        #----------------------------------------------------------

        #-------------------------------------------------
        # Note: self.SM already set to 0 by initialize()
        #-------------------------------------------------
        if (self.comp_status == 'Disabled'): return
        self.status = 'updating'  # (OpenMI)
        #-------------------------
        # Update computed values
        #-------------------------
        self.update_ground_temperatures()
        self.update_ALT()

        #-----------------------------------------
        # Read next perm vars from input files ? NOTE: does not work see the read_input_files()
        #-------------------------------------------
        # Note that read_input_files() is called
        # by initialize() and these values must be
        # used for "update" calls before reading
        # new ones.
        #-------------------------------------------
        #----------------------------------------------
        # Write user-specified data to output files ?
        #----------------------------------------------
        # Components use own self.time_sec by default.
        #-----------------------------------------------
        
        self.save_grids()

        #-----------------------------
        # Update internal clock
        # after write_output_files()
        #-----------------------------
        
        self.read_input_files()
        self.year += self.dt
        
        self.update_time(dt)
        
        self.status = 'updated'  # (OpenMI)
        
        
    def save_grids(self):
        
        if self.SAVE_ALT_GRIDS or self.SAVE_TPS_GRIDS:
            self.out_timeid[self.time]    = self.year
                
        if (self.SAVE_ALT_GRIDS):
            self.alt_out_varid[self.time,:,:] = self.Zal

        if (self.SAVE_TPS_GRIDS):
            self.tps_out_varid[self.time,:,:] = self.Tps
    
    def open_file_KU(self, var_type, input_file):
    
        #-----------------------------------------------------
        # Note:  This method's name cannot be "open" because
        #        it calls Python's built-in "open()" method.
        #-----------------------------------------------------
        # print 'var_type   =', var_type
        # print 'input_file =', input_file
    
        #--------------------------------------------
        # A scalar input value was provided already
        #--------------------------------------------
        file_unit = None
        if (var_type.lower() == 'scalar'):
            return file_unit
        if (input_file == ''):
            print('ERROR in model_input.open_file():')
            print('    Input file is null string.')
            # print '    variable type =' + var_type
            return file_unit
    
        #----------------------------------
        # Does input file exist locally ?
        #----------------------------------
        if not(os.path.exists(input_file)):
            print('ERROR in model_input.open_file():')
            print('    Could not find input file =')
            print('    ' + input_file)
            # print '   ' + input_file
            return file_unit
    
        if (var_type.lower() == 'time_series'):
            #-----------------------------------------
            # Input file contains a time series and
            # is ASCII text with one value per line.
            #-----------------------------------------
            file_unit = open(input_file, 'r')
        else:
            #--------------------------------------------
            # Input file contains a grid or grid stack
            # as row-major, binary file with no header.
            #--------------------------------------------
            from netCDF4 import Dataset
            file_unit = Dataset(input_file, 'r')
    
        return file_unit

