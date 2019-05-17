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

        self.ALT_file         = self.out_directory + self.ALT_file
        self.TPS_file         = self.out_directory + self.TPS_file
        
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
        
    #   open_input_files()
    #-------------------------------------------------------------------
    def read_input_files(self):

        #rti = self.rti # has a problem with loading rti: do not know where its been initialized

        #-------------------------------------------------------
        # All grids are assumed to have a data type of Float32.
        #-------------------------------------------------------    
 
        #--------------  AIR              
        T_air = self.read_next_modified_KU(self.T_air_unit,self.T_air_type)        
        if (T_air is not None):
            self.T_air = T_air            
        A_air = self.read_next_modified_KU(self.A_air_unit, self.A_air_type)
        if (A_air is not None): 
            self.A_air = A_air
        #--------------  SNOW                          
        h_snow = self.read_next_modified_KU(self.h_snow_unit, self.h_snow_type)
        if (h_snow is not None): 
            self.h_snow = h_snow
        k_snow = self.read_next_modified_KU(self.k_snow_unit, self.k_snow_type)
        if (k_snow is not None): 
            self.k_snow = k_snow
        c_snow = self.read_next_modified_KU(self.c_snow_unit, self.c_snow_type)
        if (c_snow is not None): 
            self.c_snow = c_snow            
        rho_snow = self.read_next_modified_KU(self.rho_snow_unit, self.rho_snow_type)
        if (rho_snow is not None): 
            self.rho_snow = rho_snow
        #--------------   VEG                         
        Hvgf = self.read_next_modified_KU(self.Hvgf_unit, self.Hvgf_type)
        if (Hvgf is not None): 
            self.Hvgf = Hvgf           
        Hvgt = self.read_next_modified_KU(self.Hvgt_unit, self.Hvgt_type)
        if (Hvgt is not None): 
            self.Hvgt = Hvgt
        Dvt = self.read_next_modified_KU(self.Dvt_unit, self.Dvt_type)
        if (Dvt is not None): 
            self.Dvt = Dvt
        Dvf = self.read_next_modified_KU(self.Dvf_unit, self.Dvf_type)
        if (Dvf is not None): 
            self.Dvf = Dvf
        #--------------    SOIL
        lh_soil = self.read_next_modified_KU(self.lh_soil_unit, self.lh_soil_type)
        if (lh_soil is not None): 
            self.lh_soil = lh_soil
        kt_soil = self.read_next_modified_KU(self.kt_soil_unit, self.kt_soil_type)
        if (kt_soil is not None): 
            self.kt_soil = kt_soil
        kf_soil = self.read_next_modified_KU(self.kf_soil_unit, self.kf_soil_type)
        if (kf_soil is not None): 
            self.kf_soil = kf_soil
        ct_soil = self.read_next_modified_KU(self.ct_soil_unit, self.ct_soil_type)
        if (ct_soil is not None): 
            self.ct_soil = ct_soil
        cf_soil = self.read_next_modified_KU(self.cf_soil_unit, self.cf_soil_type)
        if (cf_soil is not None): 
            self.cf_soil = cf_soil            
        
    #   read_input_files()
    #-------------------------------------------------------------------

    def update_soil_heat_capacity(self):

        self.Ct = self.ct_soil + 0. 
        self.Cf = self.cf_soil + 0.
        
    #   update_soil_heat_capacity()
    #-------------------------------------------------------------------

    def update_soil_thermal_conductivity(self):
        
        self.Kt = self.kt_Soil + 0.;
        self.Kf = self.kf_Soil + 0.;
        
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

        tao = self.T_air*0.0 + self.sec_per_year;        
        
        K_diffusivity = self.Ksn/(self.rho_snow*self.Csn)
        
        temp = np.exp(-1.0*self.h_snow*np.sqrt(np.pi/(tao*K_diffusivity)))
        deta_Tsn = self.A_air*(1.0 - temp);
        deta_Asn = deta_Tsn*2.0/np.pi;

        Tvg = self.T_air + deta_Tsn;
        Avg = self.A_air - deta_Asn;
        
        self.deta_Tsn = deta_Tsn;
        self.deta_Asn = deta_Asn;

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
        
        self.deta_Tv = deta_Tv;
        self.deta_Av = deta_Av;

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

        n_grid = np.size(self.T_air)
        
        if n_grid >1 :               
        
            K_star = self.Kf
            
            if np.size(self.Kf)>1:
            	K_star[np.where(Tps_numerator>0.0)] = self.Kt[np.where(Tps_numerator>0.0)];
            
        else:
            if Tps_numerator<=0.0:
                K_star = self.Kf
            else:
                K_star = self.Kt
            

                
        self.Tgs=Tgs
        self.Ags=Ags
        self.Tps_numerator=Tps_numerator

        self.Tps = self.Tps_numerator/K_star
        
        if n_grid > 1:        
        
            self.Tps[np.where(self.Tps_numerator>0.0)] = np.nan # Seasonal Frozen Ground
            
        else:
            
            if self.Tps_numerator>0.0:
                self.Tps = np.nan


    #   update_TOP_temperatures()
    #-------------------------------------------------------------------
    def update_ALT(self):

        #---------------------------------------------------------
        #       Calculates active layer thickness
        #       Aps  -- eq-4, Romanovsky et al. 1997
        #       Zs -- eq-5, Romanovsky et al. 1997
        #       Zal -- eq-3, Romanovsky et al. 1997
        #--------------------------------------------------

        tao = self.T_air*0.0 + self.sec_per_year;

        n_grid = np.size(self.T_air) 

        if n_grid > 1:        

            K = self.Kt
            C = self.Ct       
            if np.size(self.Kf)>1:        
            	K[np.where(self.Tps_numerator>0.0)] = self.Kf[np.where(self.Tps_numerator>0.0)]
            	C[np.where(self.Tps_numerator>0.0)] = self.Cf[np.where(self.Tps_numerator>0.0)]
            
        else:

            if self.Tps_numerator<=0.0:
                K = self.Kt
                C = self.Ct
            else:
                K = self.Kf
                C = self.Cf

        Aps = (self.Ags - abs(self.Tps))/np.log((self.Ags+self.L/(2.*C)) / \
                    (abs(self.Tps)+self.L/(2.*C))) - self.L/(2.*C);

        Zc = (2.*(self.Ags - abs(self.Tps))*np.sqrt((K*tao*C)/np.pi)) / \
                    (2.*Aps*C + self.L);
                    
        Zal = (2.*(self.Ags - abs(self.Tps))*np.sqrt(K*tao*C/np.pi) \
                +(2.*Aps*C*Zc+self.L*Zc)*self.L*np.sqrt(K*tao/(np.pi*C)) \
                /(2.*Aps*C*Zc + self.L*Zc +(2.*Aps*C+self.L)*np.sqrt(K*tao/(np.pi*C)))) \
                /(2.*Aps*C+ self.L);

        if n_grid > 1:        
        
            Zal[np.where(Zal<=0.01)] = np.nan
            Zal[np.where(self.Tps_numerator>0.0)] = np.nan # Seasonal Frozen Ground
            Zal[np.where(np.isnan(Zal))] = np.nan
            
        else:
            
            if self.Tps_numerator>0.0 or Zal<=0.01 or np.isnan(Zal):
                Zal = np.nan

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
        
        tao = self.T_air*0.0 + self.sec_per_year;

        # Update mean temperatures for warmes and coldest seasons similar to Nelson & Outcalt 87
        # Cold and Warm Season, Page-129, Sazonova, 2003
        self.tao1 = tao*(0.5 - 1./np.pi*np.arcsin(self.T_air/self.A_air));
        self.tao2 = tao - self.tao1;
        self.L=334000.*1000.*self.vwc_H2O

        self.update_TOP_temperatures()

    #   update_ground_temperatures()
    #-------------------------------------------------------------------
    def close_input_files(self):

        if (self.T_air_type     != 'Scalar'): self.T_air_unit.close()
        if (self.A_air_type     != 'Scalar'): self.A_air_unit.close()
        if (self.h_snow_type    != 'Scalar'): self.h_snow_unit.close()
        if (self.rho_snow_type  != 'Scalar'): self.rho_snow_unit.close()
        if (self.vwc_H2O_type   != 'Scalar'): self.vwc_H2O_unit.close()
        if (self.Hvgf_type      != 'Scalar'): self.Hvgf_unit.close()
        if (self.Hvgt_type      != 'Scalar'): self.Hvgt_unit.close()
        if (self.Dvf_type       != 'Scalar'): self.Dvf_unit.close()
        if (self.Dvt_type       != 'Scalar'): self.Dvt_unit.close()
#        if (self.lat_type       != 'Scalar'): self.lat_unit.close()
#        if (self.lon_type       != 'Scalar'): self.lon_unit.close()    

    #   close_input_files()
    #-------------------------------------------------------------------
    
    def Extract_Soil_Texture_Loops(self):
        
        n_lat = np.size(self.lat)
        n_lon = np.size(self.lon)
        
        n_grid = n_lat*n_lon     
        
        if n_grid > 1:
            
            p_clay_list = np.zeros((n_lat,n_lon));
            p_sand_list = np.zeros((n_lat,n_lon));
            p_silt_list = np.zeros((n_lat,n_lon));
            p_peat_list = np.zeros((n_lat,n_lon));
            
#            lon = np.reshape(self.lon, (n_grid,1))
#            lat = np.reshape(self.lat, (n_grid,1))
                    
            for i in range(n_lon):
                for j in range(n_lat):
                
                    input_lat   = self.lat[j]
                    input_lon   = self.lon[i]
                    
                    [p_clay0, p_sand0, p_silt0, p_peat0] = self.Extract_Soil_Texture(input_lat, input_lon);
                
                    p_clay_list[j,i] = p_clay0            
                    p_sand_list[j,i] = p_sand0        
                    p_silt_list[j,i] = p_silt0        
                    p_peat_list[j,i] = p_peat0
        else:
            
            input_lat   = self.lat
            input_lon   = self.lon
                
            [p_clay0, p_sand0, p_silt0, p_peat0] = self.Extract_Soil_Texture(input_lat, input_lon);
            
            p_clay_list = p_clay0            
            p_sand_list = p_sand0        
            p_silt_list = p_silt0        
            p_peat_list = p_peat0*0.
                         
                    
        self.p_clay = p_clay_list
        self.p_sand = p_sand_list
        self.p_silt = p_silt_list
        self.p_peat = p_peat_list
        
    def Extract_Soil_Texture_Loops_New(self):
        
        [p_clay_list, p_sand_list, p_silt_list, p_peat_list] = self.Extract_Soil_Texture2();
        
        self.p_clay = p_clay_list
        self.p_sand = p_sand_list
        self.p_silt = p_silt_list
        self.p_peat = p_peat_list*0.0

    def Extract_Soil_Texture(self, input_lat, input_lon): 
    
        """ 
        The function is to extract the grid value from matrix,
        according to input of latitude and longitude;
        
        INPUTs:
                input_lat: Latitude;
                input_lon: Longitude;
                lon_grid : Array of longitude
                lat_grid : Array of latitude
                p_data   : Matrix of data (from NetCDF file)
                
        OUTPUTs:
                q_data: grid value (SINGLE)   
                        
        DEPENDENTs:
                None 
        """
            
        import numpy as np
        
        lon_grid_scale = 0.05;
        lat_grid_scale = 0.05;
        
        lon_grid_top = self.lon_grid + lon_grid_scale / 2.0;
        lat_grid_top = self.lat_grid + lat_grid_scale / 2.0;
        
        lon_grid_bot = self.lon_grid - lon_grid_scale / 2.0;
        lat_grid_bot = self.lat_grid - lat_grid_scale / 2.0;
        
        # Get the index of input location acccording to lat and lon inputed
        
        idx_lon = np.where((input_lon <= lon_grid_top) & (input_lon >= lon_grid_bot))          
        idx_lat = np.where((input_lat <= lat_grid_top) & (input_lat >= lat_grid_bot))
        
        idx_lon = np.array(idx_lon)
        idx_lat = np.array(idx_lat)
        
        if np.size(idx_lon) >= 1 and np.size(idx_lat) >= 1:
            clay_perc  = self.Clay_percent[idx_lat[0,0], idx_lon[0,0]]
            sand_perc  = self.Sand_percent[idx_lat[0,0], idx_lon[0,0]]
            silt_perc  = self.Silt_percent[idx_lat[0,0], idx_lon[0,0]]
            peat_perc  = self.Peat_percent[idx_lat[0,0], idx_lon[0,0]]
        else:
            clay_perc  = np.nan;
            sand_perc  = np.nan;
            silt_perc  = np.nan;
            peat_perc  = np.nan;            
    
        return clay_perc, sand_perc, silt_perc, peat_perc

    def Extract_Soil_Texture2(self): 
        
        import numpy as np
        from affine import Affine
    
        lon_cell_size = abs(self.lon_grid[0] - self.lon_grid[1])
        lat_cell_size = abs(self.lat_grid[0] - self.lat_grid[1])
    
        min_lon = min(self.lon_grid) - lon_cell_size/2.0*0.
        min_lat = min(self.lat_grid) - lat_cell_size/2.0*0. 
        
        n_lat = np.size(self.lat)
        n_lon = np.size(self.lon)
    
        aff = Affine.from_gdal(min_lon, lon_cell_size, 0.0, min_lat, 0.0, lat_cell_size)
        
        lon = np.reshape(np.repeat(self.lon,n_lat), (n_lon,n_lat));
        lat = np.transpose(np.reshape(np.repeat(self.lat,n_lon), (n_lat,n_lon)));
    
        x_coords, y_coords = ~aff * (lon, lat)
    
        x_coords = np.round(x_coords).astype(np.int)
        y_coords = np.round(y_coords).astype(np.int)
    
        if np.size(x_coords) >= 1 and np.size(y_coords) >= 1:
            
            clay_perc0  = self.Clay_percent[y_coords, x_coords]
            sand_perc0  = self.Sand_percent[y_coords, x_coords]
            silt_perc0  = self.Silt_percent[y_coords, x_coords]
            peat_perc0  = self.Peat_percent[y_coords, x_coords]

            clay_perc  =  np.transpose(np.reshape(clay_perc0, (n_lon, n_lat)))
            sand_perc  =  np.transpose(np.reshape(sand_perc0, (n_lon, n_lat)))           
            silt_perc  =  np.transpose(np.reshape(silt_perc0, (n_lon, n_lat)))
            peat_perc  =  np.transpose(np.reshape(peat_perc0, (n_lon, n_lat)))
            
            
        else:
            clay_perc  = np.nan;
            sand_perc  = np.nan;
            silt_perc  = np.nan;
            peat_perc  = np.nan;
            
        return clay_perc, sand_perc, silt_perc, peat_perc
        
    def read_whole_soil_texture_from_GSD(self):
        
        Clay_file = self.get_param_nc4_filename("T_CLAY")
        Sand_file = self.get_param_nc4_filename("T_SAND")
        Silt_file = self.get_param_nc4_filename("T_SILT")
        Peat_file = self.get_param_nc4_filename("T_OC")

        lonname    = 'lon';
        latname    = 'lat';
        
        varname    = 'T_CLAY';                                        
        [lat_grid, lon_grid, Clay_percent] = self.import_ncfile(Clay_file, 
                                                lonname, latname, varname)
        varname    = 'T_SAND';                                        
        [lat_grid, lon_grid, Sand_percent] = self.import_ncfile(Sand_file, 
                                                lonname, latname, varname)
         
        varname    = 'T_SILT';                                        
        [lat_grid, lon_grid, Silt_percent] = self.import_ncfile(Silt_file, 
                                                lonname, latname, varname)
                                                
        varname    = 'T_OC';                                        
        [lat_grid, lon_grid, Peat_percent] = self.import_ncfile(Peat_file, 
                                                lonname, latname, varname)
        
        self.Clay_percent = Clay_percent;
        self.Sand_percent = Sand_percent;
        self.Silt_percent = Silt_percent;
        self.Peat_percent = Peat_percent;
        self.lon_grid     = lon_grid;
        self.lat_grid     = lat_grid;
    
    def import_ncfile(self, input_file, lonname,  latname,  varname): 
                                           
        from netCDF4 import Dataset
        
        # Read the nc file 
        
        fh = Dataset(input_file, mode='r')
        
        # Get the lat and lon
        
        lon_grid = fh.variables[lonname][:]; 
        lat_grid = fh.variables[latname][:];
        
        p_data  = fh.variables[varname][:];
        
        return lat_grid,lon_grid,p_data
        
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
        self.Tps = np.float32(-999.99)
        self.Zal = np.float32(-999.99)
        self.cont = 0.0

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

        if (self.comp_status == 'Disabled'):
            #########################################
            #  DOUBLE CHECK THIS; SEE NOTES ABOVE
            #########################################
               ####### and (ep.method != 2):  ??????
            if not(SILENT):
                print('Permafrost component: Disabled.')
            self.lat    = self.initialize_scalar(0, dtype='float64')
            self.lon    = self.initialize_scalar(0, dtype='float64')
            self.start_year    = self.initialize_scalar(0, dtype='float64')
            self.end_year    = self.initialize_scalar(0, dtype='float64')
            self.T_air  = self.initialize_scalar(0, dtype='float64')
            self.h_snow = self.initialize_scalar(0, dtype='float64')
            self.vwc_H2O= self.initialize_scalar(0, dtype='float64')
            self.Hvgf   = self.initialize_scalar(0, dtype='float64')
            self.Hvgt   = self.initialize_scalar(0, dtype='float64')
            self.Dvf    = self.initialize_scalar(0, dtype='float64')
            self.Dvt    = self.initialize_scalar(0, dtype='float64')
            self.DONE   = True
            self.status = 'initialized'
            return

        #---------------------------------------------
        # Open input files needed to initialize vars
        #---------------------------------------------
        self.open_input_files()
        self.read_input_files()
        
        #        self.read_nc_lat_lon(self, file_name, var_type)

        self.status = 'initialized'
    
    def read_nc_lat_lon(self, file_unit, var_type):
        
        if (var_type.lower() == 'scalar'):
            #-------------------------------------------
            # Scalar value was entered by user already
            #-------------------------------------------
            lat = None
            lon = None
            
        elif (var_type.lower() == 'time_series'):
            #----------------------------------------------
            # Time series: Read scalar value from file.
            # File is ASCII text with one value per line.
            #----------------------------------------------
            lat = None
            lon = None
            
        elif (var_type.lower() == 'grid'):
            for var in file_unit.variables.keys():
                if var[0:3] =='lat':
                    lat = file_unit.variables[var][:]
                if var[0:3] =='lon':
                    lon = file_unit.variables[var][:]
#            lon = self.ncread(file_unit, 'longitude')
            
            if (np.min(lon) > 0.):
            	lon = np.mod((lon+180.),360.) -180.
            
#            lat = np.float(lat)
#            lon = np.float(lon)
            
        else:
            raise RuntimeError('No match found for "var_type".')
            return None            
        
        if (lat is None):
            return
        else:
            return lat,lon
    
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
            from netCDF4 import Dataset
            for var in file_unit.variables.keys():
                if (var != 'time' and var[0:3] !='lat' and var[0:3] != 'lon'):
                    data  = file_unit.variables[var][self.cont,:,:]
                                   
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
            return np.float64( data )

    def ncread(self, input_file, varname):

        from netCDF4 import Dataset
        
        f = Dataset(input_file, mode='r') # Open the nc file -> handle
    
        data  = f.variables[varname][:]
    
        f.close()
    
        return data
        
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
        if (self.time_index > 0):
            self.read_input_files()

        #----------------------------------------------
        # Write user-specified data to output files ?
        #----------------------------------------------
        # Components use own self.time_sec by default.
        #-----------------------------------------------
        if (self.SAVE_ALT_GRIDS):       
        	self.save_grids()

        #-----------------------------
        # Update internal clock
        # after write_output_files()
        #-----------------------------
        self.update_time( dt )
        self.status = 'updated'  # (OpenMI)

    #   update()
    
#     def finalize(self):
#         SILENT = True
# 
#         # Finish with the run
#         self._model.status = 'finalizing'  # (OpenMI)
# 
#         # Close the input files
#         self._model.close_input_files()   # Close any input files
# 
#         # Write output last output
#         # self._model.write_output_to_file(SILENT=True)
# 
#         if (self.SAVE_ALT_GRIDS):       
#             self.save_grids()
# 
#         # Close the output files
#         self._model.close_output_files()
# 
#         # Done finalizing  
#         self._model.status = 'finalized'  # (OpenMI)
# 
#         # Print final report, as desired
#         if not SILENT:
#             self._model.print_final_report(\
#                     comp_name='Permamodel Ku component')
    # finalize()

   
#    def save_grids(self):
#        # Saves the grid values based on the prescribed ones in cfg file
#
#        #if (self.SAVE_MR_GRIDS):
#        #    model_output.add_grid( self, self.T_air, 'T_air', self.time_min )
#        self.ALT_file  = self.out_directory + self.ALT_file
#        
#        if (self.SAVE_ALT_GRIDS):
#            self.write_out_ncfile(self.ALT_file,self.output_alt)
##            self.write_out_ncfile(self.ALT_file,self.Zal)
#            
#        self.TPS_file  = self.out_directory + self.TPS_file
#        
#        if (self.SAVE_TPS_GRIDS):
#            self.write_out_ncfile(self.TPS_file,self.output_tps)
#            self.write_out_ncfile(self.TPS_file,self.Tps)
            
        #if (self.SAVE_SW_GRIDS):
        #    model_output.add_grid( self, self.Tps, 'Tps', self.time_min )

        #if (self.SAVE_CC_GRIDS):
        #    model_output.add_grid( self, self.Zal, 'Zal', self.time_min )

    def close_output_files(self):
        
        tst = 'in progressing';

        #if (self.SAVE_MR_GRIDS): model_output.close_gs_file( self, 'mr')
#        if (self.SAVE_HS_GRIDS): model_output.close_gs_file( self, 'hs')
        #if (self.SAVE_SW_GRIDS): model_output.close_gs_file( self, 'sw')
        #if (self.SAVE_CC_GRIDS): model_output.close_gs_file( self, 'cc')
        #-----------------------------------------------------------------
        #if (self.SAVE_MR_PIXELS): model_output.close_ts_file( self, 'mr')
        #if (self.SAVE_HS_PIXELS): model_output.close_ts_file( self, 'hs')
        #if (self.SAVE_SW_PIXELS): model_output.close_ts_file( self, 'sw')
        #if (self.SAVE_CC_PIXELS): model_output.close_ts_file( self, 'cc')
        
    def write_out_ncfile(self, output_file, varname):

        from netCDF4 import Dataset
        import numpy as np
        
        n_lat = np.size(self.lat)
        n_lon = np.size(self.lon)
        
#        print np.shape(varname)
        
        ALT = varname + 0.0 #self.mask;
        idx = np.where(np.isnan(ALT))
        ALT[idx] = -999.99;
        
        #print output_file[-1-2]
        
        if (output_file[-1-2] == 'T'):
        	units = 'degree C';
        	long_name = 'Temperature at top of permafrost';
        else:
        	units = 'm'; 
        	long_name = 'Active Layer Thickness'
        
        # Open a file to save the final result
        w_nc_fid = Dataset(output_file+'.nc', 'w', format='NETCDF4');
        
        # ==== Latitude ====

        w_nc_fid.createDimension('lat', n_lat) # Create Dimension
        lats = w_nc_fid.createVariable('lat',np.dtype('float32').char,('lat',))
        lats.units = 'degrees_north'
        lats.standard_name = 'latitude'
        lats.long_name = 'latitude'
        lats.axis = 'Y'
        lats[:] = self.lat
        
        # ==== Longitude ====

        w_nc_fid.createDimension('lon', n_lon) # Create Dimension
        lons = w_nc_fid.createVariable('lon',np.dtype('float32').char,('lon',))
        lons.units = 'degrees_east'
        lons.standard_name = 'longitude'
        lons.long_name = 'longitude'
        lons.axis = 'X'
        lons[:] = self.lon
        
        # ==== Time ====

        w_nc_fid.createDimension('time', self.end_year-self.start_year+1.0) # Create Dimension
        time = w_nc_fid.createVariable('time',np.dtype('float32').char,('time',))
        time.units = 'Year'
#         time.standard_name = 'longitude'
#         time.long_name = 'longitude'
        time.axis = 'Z'
        time[:] = np.linspace(self.start_year, self.end_year, self.end_year-self.start_year+1.0)
               
        # ==== Data ====
        temp = w_nc_fid.createVariable('data',np.dtype('float32').char,('time','lat','lon'))
        temp.units = units
        temp.missing_value = -999.99
        temp.long_name = long_name
        temp[:] = ALT;
#        
        w_nc_fid.close()  # close the new file
        
    def write_out_txtfile(self, output_file, varname):

        import numpy as np
        
        ALT = self.Zal+0.;
        
        # Open a file to save the final result
        np.savetxt(output_file+'.txt',self.lat);
    
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

