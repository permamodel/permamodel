# -*- coding: utf-8 -*-
"""  Kudryavtsev Model code adapted for the BMI version developed for the Topoflow model

     Author: Kang Wang, 03/29/2016
     Modified: Elchin Jafarov, 03/29/2016

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

import numpy as np
from permamodel.utils import model_input
from permamodel.components import perma_base

class Ku_method( perma_base.permafrost_component ):

    #-------------------------------------------------------------------
    _att_map = {
    # NOTE: this will change in the future
        'model_name':         'PermaModel_Kudryavtsev_method',
        'version':            '0.1',
        'author_name':        'Kang Wang and Elchin Jafarov',
        'grid_type':          'none',
        'time_step_type':     'fixed',
        'step_method':        'explicit',
        #-------------------------------------------------------------
        'comp_name':          'Ku_model',
        'model_family':       'PermaModel',
        'cfg_extension':      '_ku_model.cfg',
        'cmt_var_prefix':     '/input/',
        'gui_yaml_file':      '/input/ku_model.yaml',
        'time_units':         'years' }

    _input_var_names = [
        'latitude',
        'longitude',
        'atmosphere_bottom_air__temperature',
        'atmosphere_bottom_air__temperature_amplitude',
        'snowpack__depth',
        'snowpack__density',
        'water-liquid__volumetric-water-content-soil',
        'vegetation__Hvgf',
        'vegetation__Hvgt',
        'vegetation__Dvf',
        'vegetation__Dvt' ]

    _output_var_names = [
        'soil__temperature',                                  # Tps
        'soil__active_layer_thickness' ]                      # Zal

    _var_name_map = {
    # NOTE: we need to look up for the corresponding standard names
        'latitude':                                           'lat',
        'longitude':                                          'lon',
        'atmosphere_bottom_air__temperature':                 'T_air',
        'atmosphere_bottom_air__temperature_amplitude':       'A_air',
        'snowpack__depth':                                    'h_snow',
        'snowpack__density':                                  'rho_snow',
        'water-liquid__volumetric-water-content-soil':        'vwc_H2O',
        'vegetation__Hvgf':                                   'Hvgf',
        'vegetation__Hvgt':                                   'Hvgt',
        'vegetation__Dvf':                                    'Dvf',
        'vegetation__Dvt':                                    'Dvt' }

    _var_units_map = {
    # NOTE: Kang please complete the vegetation info both on var names and units
        'latitude':                                           'lat',
        'longitude':                                          'lon',
        'atmosphere_bottom_air__temperature':                 'deg_C',
        'atmosphere_bottom_air__temperature_amplitude':       'deg_C',
        'snowpack__depth':                                    'm',
        'snowpack__density':                                  'kg m-3',
        'water-liquid__volumetric-water-content-soil':        'm3 m-3',
        'vegetation__Hvgf':                                   'm',
        'vegetation__Hvgt':                                   'm',
        'vegetation__Dvf':                                    'm2 s',
        'vegetation__Dvt':                                    'm2 s' }

    #-------------------------------------------------------------------
    def get_attribute(self, att_name):

        try:
            return self._att_map[ att_name.lower() ]
        except:
            print '###################################################'
            print ' ERROR: Could not find attribute: ' + att_name
            print '###################################################'
            print ' '

    #   get_attribute()
    #-------------------------------------------------------------------
    def get_input_var_names(self):

        #--------------------------------------------------------
        # Note: These are currently variables needed from other
        #       components vs. those read from files or GUI.
        #--------------------------------------------------------
        return self._input_var_names

    #   get_input_var_names()
    #-------------------------------------------------------------------
    def get_output_var_names(self):

        return self._output_var_names

    #   get_output_var_names()
    #-------------------------------------------------------------------
    def get_var_name(self, long_var_name):

        return self._var_name_map[ long_var_name ]

    #   get_var_name()
    #-------------------------------------------------------------------
    def get_var_units(self, long_var_name):

        return self._var_units_map[ long_var_name ]

    #   get_var_units()
    #-------------------------------------------------------------------
    def check_input_types(self):

        #--------------------------------------------------
        # Notes: rho_H2O, Cp_snow, rho_air and Cp_air are
        #        currently always scalars.
        #--------------------------------------------------
        are_scalars = np.array([
                          self.is_scalar('lat'),
                          self.is_scalar('lon'),
                          self.is_scalar('T_air'),
                          self.is_scalar('A_air'),
                          self.is_scalar('h_snow'),
                          self.is_scalar('rho_snow'),
                          self.is_scalar('vwc_H2O'),
                          self.is_scalar('Hvgf'),
                          self.is_scalar('Hvgt'),
                          self.is_scalar('Dvf'),
                          self.is_scalar('Dvt') ])

        self.ALL_SCALARS = np.all(are_scalars)

    #   check_input_types()
    #-------------------------------------------------------------------
    def open_input_files(self):
        # this function will work only if filename is not empty
        self.T_air_file       = self.in_directory + self.T_air_file
        self.A_air_file       = self.in_directory + self.A_air_file
        self.h_snow_file      = self.in_directory + self.h_snow_file
        self.rho_snow_file    = self.in_directory + self.rho_snow_file
        self.vwc_H2O_file     = self.in_directory + self.vwc_H2O_file
        self.Hvgf_file        = self.in_directory + self.Hvgf_file
        self.Hvgt_file        = self.in_directory + self.Hvgt_file
        self.Dvf_file         = self.in_directory + self.Dvf_file
        self.Dvt_file         = self.in_directory + self.Dvt_file

        self.T_air_unit       = model_input.open_file(self.T_air_type,  self.T_air_file)
        self.A_air_unit       = model_input.open_file(self.A_air_type,  self.A_air_file)
        self.h_snow_unit      = model_input.open_file(self.h_snow_type,  self.h_snow_file)
        self.rho_snow_unit    = model_input.open_file(self.rho_snow_type,  self.rho_snow_file)
        self.vwc_H2O_unit     = model_input.open_file(self.vwc_H2O_type,  self.vwc_H2O_file)
        self.Hvgf_unit        = model_input.open_file(self.Hvgf_type,  self.Hvgf_file)
        self.Hvgt_unit        = model_input.open_file(self.Hvgt_type,  self.Hvgt_file)
        self.Dvf_unit         = model_input.open_file(self.Dvf_type,  self.Dvf_file)
        self.Dvt_unit         = model_input.open_file(self.Dvt_type,  self.Dvt_file)

    #   open_input_files()
    #-------------------------------------------------------------------
    #def read_input_files(self):

        #rti = self.rti # has a problem with loading rti: do not know where its been initialized

        #-------------------------------------------------------
        # All grids are assumed to have a data type of Float32.
        #-------------------------------------------------------
        #T_air = model_input.read_next(self.T_air_unit, self.T_air_type, rti)
        #if (T_air != None): self.T_air = T_air
        #print T_air

        #T0 = model_input.read_next(self.T0_unit, self.T0_type, rti)
        #if (T0 != None): self.T0 = T0

        #rho_snow = model_input.read_next(self.rho_snow_unit, self.rho_snow_type, rti)
        #if (rho_snow != None): self.rho_snow = rho_snow

        #h0_snow = model_input.read_next(self.h0_snow_unit, self.h0_snow_type, rti)
        #if (h0_snow != None): self.h0_snow = h0_snow

       # h0_swe = model_input.read_next(self.h0_swe_unit, self.h0_swe_type, rti)
        #if (h0_swe != None): self.h0_swe = h0_swe

    #   read_input_files()
    #-------------------------------------------------------------------


    def update_soil_heat_capacity(self):

        #---------------------------------------------------------
        # Notes: we need a better documentation of this subroutine here
        #
        #
        #---------------------------------------------------------
        # Note: need to update frozen and thawed (Cf,Ct)
        #       heat capacities yearly
        #       this methods overriddes the method in the perma_base
        #
        #
        #--------------------------------------------------
        # I do not like this input file here need fix later
        #input_file = 'Parameters/Typical_Thermal_Parameters.csv'
        input_file = self.permafrost_dir + '/permamodel/components/Parameters/Typical_Thermal_Parameters.csv'
        s_data = np.genfromtxt(input_file, names = True, delimiter=',', dtype=None)

        Bulk_Density_Texture = s_data['Bulk_Density']
        Heat_Capacity_Texture = s_data['Heat_Capacity']

        # Adjusting percent of sand, silt, clay and peat ==
        tot_percent = self.p_sand+self.p_clay+self.p_silt+self.p_peat

        percent_sand = self.p_sand / tot_percent
        percent_clay = self.p_clay / tot_percent
        percent_silt = self.p_silt / tot_percent
        percent_peat = self.p_peat / tot_percent

        # Calculate heat capacity and bulk density of soil using exponential weighted.
        Heat_Capacity =  Heat_Capacity_Texture[2]**percent_clay * \
                         Heat_Capacity_Texture[1]**percent_sand * \
                         Heat_Capacity_Texture[0]**percent_silt * \
                         Heat_Capacity_Texture[3]**percent_peat       # Unit: J kg-1 C-1

        Bulk_Density  =  Bulk_Density_Texture[2]**percent_clay * \
                         Bulk_Density_Texture[1]**percent_sand * \
                         Bulk_Density_Texture[0]**percent_silt * \
                         Bulk_Density_Texture[3]**percent_peat        # Unit: kg m-3
        # Estimate heat capacity for composed soil
        # based on the empirical approaches suggested by Anisimov et al. (1997)
        self.Ct = Heat_Capacity*Bulk_Density + 4190.*self.vwc_H2O # eq-15, Anisimov et al. 1997; Unit: J m-3 C-1
        self.Cf = Heat_Capacity*Bulk_Density + 2025.*self.vwc_H2O # eq-15, Anisimov et al. 1997; Unit: J m-3 C-1

    #   update_soil_heat_capacity()
    #-------------------------------------------------------------------

    def update_soil_thermal_conductivity(self):

        #---------------------------------------------------------
        # Notes: we need a better documentation of this subroutine here
        #
        #
        #---------------------------------------------------------
        # Note: need to update frozen and thawed (kf,kt)
        #       thermal conductivities yearly
        #       this methods overriddes the method in the perma_base
        #
        #
        #--------------------------------------------------
        #input_file = 'Parameters/Typical_Thermal_Parameters.csv'
        input_file = self.permafrost_dir + '/permamodel/components/Parameters/Typical_Thermal_Parameters.csv'
        s_data = np.genfromtxt(input_file, names = True, delimiter=',', dtype=None)

        Thermal_Conductivity_Thawed_Texture = s_data['Thermal_Conductivity_Thawed']
        Thermal_Conductivity_Frozen_Texture = s_data['Thermal_Conductivity_Frozen']

        # Adjusting percent of sand, silt, clay and peat ==
        tot_percent = self.p_sand+self.p_clay+self.p_silt+self.p_peat

        percent_sand = self.p_sand / tot_percent
        percent_clay = self.p_clay / tot_percent
        percent_silt = self.p_silt / tot_percent
        percent_peat = self.p_peat / tot_percent
        # Estimate thermal conductivity for composed soil
        Kt_Soil =  Thermal_Conductivity_Thawed_Texture[0]**percent_silt * \
               Thermal_Conductivity_Thawed_Texture[2]**percent_clay * \
               Thermal_Conductivity_Thawed_Texture[1]**percent_sand * \
               Thermal_Conductivity_Thawed_Texture[3]**percent_peat

        Kf_Soil =  Thermal_Conductivity_Frozen_Texture[0]**percent_silt * \
               Thermal_Conductivity_Frozen_Texture[2]**percent_clay * \
               Thermal_Conductivity_Frozen_Texture[1]**percent_sand * \
               Thermal_Conductivity_Frozen_Texture[3]**percent_peat

        # Consider the effect of water content on thermal conductivity
        vwc=self.vwc_H2O
        self.Kt = Kt_Soil**(1.0-vwc)*0.54**vwc #   Unit: (W m-1 C-1)
        self.Kf = Kf_Soil**(1.0-vwc)*2.35**vwc #   Unit: (W m-1 C-1)

    #   update_soil_thermal_conductivity()
    #-------------------------------------------------------------------
    def update_snow_thermal_properties(self):

        #---------------------------------------------------------
        # Notes: we need a better documentation of this subroutine here
        # Conductivity of snow:  eq-4, Sturm et al., 1997:
        # Capacity of snow:
        #   eq-30, Ling et al., 2004; OR Table-1, Goodrich, 1982.
        #---------------------------------------------------------
        # Note: need to update frozen and thawed (kf,kt)
        #       thermal conductivities yearly
        #       this methods overriddes the method in the perma_base
        #
        #
        #--------------------------------------------------
        rho_sn=self.rho_snow

        self.Ksn = (rho_sn/1000.)*(rho_sn/1000.)*3.233-1.01*(rho_sn/1000.)+0.138; # Unit: (W m-1 C-1)

        self.Csn = rho_sn*2.09E3;                                                # Unit: J m-3 C-1

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
        temp = np.exp(-1.0*self.h_snow*np.sqrt(np.pi/(self.sec_per_year*self.Ksn)))
        deta_Tsn = self.A_air*(1.0 - temp);
        deta_Asn = 2.0/np.pi*deta_Tsn;

        Tvg = self.T_air + deta_Tsn;
        Avg = self.A_air - deta_Asn;

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

        deta_Av = (deta_A1*self.tao1+deta_A2*self.tao2) / self.sec_per_year;

        deta_Tv = (deta_A1*self.tao1-deta_A2*self.tao2) / self.sec_per_year * (2. / np.pi)

        Tgs = Tvg + deta_Tv;
        Ags = Avg - deta_Av;

        #---------------------------------------------------------
        #   3.  Calculates Tps_Numerator;
        #       eq-14, Anisimov et al. 1997
        #--------------------------------------------------
        Tps_numerator = 0.5*Tgs*(self.Kf+self.Kt)\
                            +(Ags*(self.Kf-self.Kt)/np.pi\
                            *(Tgs/Ags*np.arcsin(Tgs/Ags)\
                            +np.sqrt(1.-(np.pi**2.0/Ags**2.0))));

        #---------------------------------------------------------
        #   4.  Calculates temperature at the top of permafrost
        #       Tps -- eq-14 cont., Anisimov et al. 1997
        #--------------------------------------------------
        if Tps_numerator<=0.0: # PERMAFROST
            K_star = self.Kf;
        else:                  # SEASONAL FROZEN GROUND
            K_star = self.Kt;

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

        if self.Tps_numerator<=0.0:
            print 'PERMAFROST'
            K = self.Kf;
            C = self.Kf;
        else:
            print 'SEASONAL FROZEN GROUND'
            K = self.Kt;
            C = self.Kt;
            self.Zal = 999.9
            self.Tps = 999.9
            return

        Aps = (self.Ags - abs(self.Tps))/np.log((self.Ags+self.L/(2.*C)) / \
                    (abs(self.Tps)+self.L/(2.*C))) - self.L/(2.*C);

        Zc = (2.*(self.Ags - abs(self.Tps))*np.sqrt((K*self.sec_per_year*C)/np.pi)) / \
                    (2.*Aps*C + self.L);

        self.Zal = (2.*(self.Ags - abs(self.Tps))*np.sqrt(K*self.sec_per_year*C/np.pi) \
                +(2.*Aps*C*Zc+self.L*Zc)*self.L*np.sqrt(K*self.sec_per_year/(np.pi*C)) \
                /(2.*self.Ags*C*Zc + self.L*Zc +(2.*Aps*C+self.L)*np.sqrt(K*self.sec_per_year/(np.pi*C)))) \
                /(2.*Aps*C+ self.L);

    #   update_ALT()
    #-------------------------------------------------------------------
    def update_ground_temperatures(self):
        # in this method there is only one output the temperature at the top of permafrost
        # TTOP
        self.update_soil_heat_capacity()
        self.update_soil_thermal_conductivity()
        self.update_snow_thermal_properties()

        # Update mean temperatures for warmes and coldest seasons similar to Nelson & Outcalt 87
        # Cold and Warm Season, Page-129, Sazonova, 2003
        self.tao1 = self.sec_per_year*(0.5 - 1./np.pi*np.arcsin(self.T_air/self.A_air));
        self.tao2 = self.sec_per_year - self.tao1;
        self.L=self.Lf*self.vwc_H2O

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

    #   close_input_files()
    #-------------------------------------------------------------------
