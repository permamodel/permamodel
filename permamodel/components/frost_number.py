# -*- coding: utf-8 -*-
"""  Frost Number by Nelson and Outcalt 1983. DOI: 10.2307/1551363. http://www.jstor.org/stable/1551363
"""

import numpy as np
from permamodel.utils import model_input
from permamodel.components import perma_base
import os
import gdal
from gdalconst import *  # Import standard constants, such as GA_ReadOnly
import osr
from pyproj import Proj, transform

class frostnumber_method( perma_base.permafrost_component ):

    #-------------------------------------------------------------------
    _att_map = {
    # NOTE: this will change in the future
        'model_name':         'PermaModel_frostnumber_method',
        'version':            '0.1',
        'author_name':        'Scott Stewart and Elchin Jafarov',
        'grid_type':          'none',
        'time_step_type':     'fixed',
        'step_method':        'explicit',
        #-------------------------------------------------------------
        'comp_name':          'frostnumber',
        'model_family':       'PermaModel',
        'cfg_extension':      '_frostnumber_model.cfg',
        'cmt_var_prefix':     '/input/',
        'gui_yaml_file':      '/input/frostnumber_model.yaml',
        'time_units':         'years' }

    _input_var_names = [
        'latitude',
        'longitude',
        'atmosphere_bottom_air__temperature_min',
        'atmosphere_bottom_air__temperature_max',
        'datetime__start',
        'datetime__end']

    _output_var_names = [
        'frostnumber__air',            # Air Frost number
        'frostnumber__surface',        # Surface Frost number
        'frostnumber__stefan' ]        # Stefan Frost number

    _var_name_map = {
    # NOTE: we need to look up for the corresponding standard names
        'latitude':                                  'lat',
        'longitude':                                 'lon',
        'atmosphere_bottom_air__temperature_min':    'T_air_min',
        'atmosphere_bottom_air__temperature_max':    'T_air_max',
        'datetime__start':                           'start_year',
        'datetime__end':                             'end_year',
        'frostnumber__air':                          'frostnumber_air',
        'frostnumber__surface':                      'frostnumber_surface',
        'frostnumber__stefan':                       'frostnumber_stefan'}

    _var_units_map = {
        'latitude':                                           'deg',
        'longitude':                                          'deg',
        'atmosphere_bottom_air__temperature_min':             'deg_C',
        'atmosphere_bottom_air__temperature_max':             'deg_C',
        'datetime__start':                                    'year',
        'datetime__end':                                      'end',
        'frostnumber__air':                                   '',
        'frostnumber__surface':                               '',
        'frostnumber__stefan':                                '' }

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
                          self.is_scalar('T_air_min'),
                          self.is_scalar('T_air_max'),
                          self.is_scalar('start_year'),
                          self.is_scalar('end_year') ])

        self.ALL_SCALARS = np.all(are_scalars)

    #   check_input_types()
    #-------------------------------------------------------------------
    def open_input_files(self):
        # this function will work only if filename is not empty
        self.T_air_min_file   = './permamodel/examples/fn_t_air_min.dat'
        self.T_air_min_unit = open(self.T_air_min_file, "r")

        self.T_air_max_file   = './permamodel/examples/fn_t_air_min.dat'
        self.T_air_max_unit = open(self.T_air_max_file, "r")

        self.start_year_file   = './permamodel/examples/fn_start_year.dat'
        self.start_year_unit = open(self.start_year_file, "r")

        self.end_year_file   = './permamodel/examples/fn_end_year.dat'
        self.end_year_unit = open(self.end_year_file, "r")

        # lat and lon not implemented yet

        #self.lat_file   = './permamodel/examples/fn_lat.dat'
        #self.lat_unit = open(self.lat_file, "r")

        #self.lon_file   = './permamodel/examples/fn_lon.dat'
        #self.lon_unit = open(self.lon_file, "r")

        # This isn't exactly "opening an input file", but it is an init
        #self.year = self.start_year

    #   open_input_files()
    #-------------------------------------------------------------------

    def read_input_files(self):
        #print("in read_input_files for year %d" % self.year)

        #rti = self.rti # has a problem with loading rti: do not know where its been initialized

        #-------------------------------------------------------
        # All grids are assumed to have a data type of Float32.
        #-------------------------------------------------------
        T_air_min = model_input.read_next_modified(self.T_air_min_unit,
                                                   'scalar')
        if (T_air_min != None): self.T_air_min = T_air_min

        T_air_max = model_input.read_next_modified(self.T_air_max_unit,
                                                   'scalar')
        if (T_air_max != None): self.T_air_max = T_air_max

        start_year = model_input.read_next_modified(self.start_year_unit,
                                                    'scalar')
        if (start_year != None): self.start_year = start_year

        end_year = model_input.read_next_modified(self.end_year_unit,
                                                  'scalar')
        if (end_year != None): self.end_year = end_year

        # Initialize the year to the start year
        self.year = self.start_year

    #   read_input_files()
    #-------------------------------------------------------------------


    def update_dd(self):

        # Input: T_hot (avg temp of warmest month)
        #        T_cold (avg temp of coldest month)

        # Output: ddf (degree freezing days)
        #         ddt (degree thawing days)

        # In the first test case, we used T_air_max and T_air_min
        T_cold=self.T_air_min
        T_hot=self.T_air_max

        '''  this section is for later when reading temp values from CRU data
        # Now, we use the values from the temperature CRU tiff files
        # Assume that warmest month is July and coldest is the following Jan
        print("Lon: %f" % self.lon)
        print("Lat: %f" % self.lat)
        T_hot = self.get_temperature_from_cru(self.lon, self.lat, 7, self.year)
        T_cold = self.get_temperature_from_cru(self.lon, self.lat, 1, self.year+1)
        '''

        print("In update_dd, year=%d" % self.year)
        assert(T_hot > T_cold)
        T_avg = (T_hot + T_cold) / 2.0
        # Note that these conditions should cover T_hot == T_cold
        if T_hot <= 0:
            # Always freezing
            # Negative sign because ddf is + and T_avg (here) is -
            ddf = -365.0 * T_avg
            ddt = 0
        elif T_cold>0:
            # Never freezing
            ddf = 0
            ddt = 365.0 * T_avg
        elif (self.T_air_type != 'Scalar'):
            #wk = np.loadtxt('examples/temp_copy.txt', skiprows=1,unpack=False)
            temperature_filename = self.permafrost_dir +\
                "permamodel/examples/temp_copy.txt"
            wk = np.loadtxt(temperature_filename, skiprows=1,unpack=False)
            t_month = wk[:,0]
            T_month = wk[:,1]
            Th=max(T_month)
            Tc=min(T_month)
            T=(Th+Tc)/2                                             #(eqn. 2.1)
            A=(Th-Tc)/2                                             #(eqn. 2.2)
            beta=np.arccos(-T/A)                                    #(eqn. 2.3)
            Ts=T+A*np.sin(beta/beta)                                #(eqn. 2.4)
            Tw=T-A*np.sin(beta/(np.pi-beta))                        #(eqn. 2.5)
            Ls=365*(beta/np.pi)                                     #(eqn. 2.6)
            Lw = 365-Ls                                             #(eqn. 2.7)
            print('winter length:',Lw,'summer length:',Ls)
            ddt = Ts*Ls                                             #(eqn. 2.8)
            ddf = -Tw*Lw                                            #(eqn. 2.9)
            print Th,Tc
        else:
            # Assume cosine fit for temp series
            A = (T_hot - T_cold) / 2.0
            Beta = np.arccos(-T_avg / A)
            T_summer = T_avg + A * np.sin(Beta) / Beta
            T_winter = T_avg - A * np.sin(Beta) / (np.pi - Beta)
            L_summer = 365.0 * Beta / np.pi
            L_winter = 365.0 - L_summer
            ddt = T_summer * L_summer
            ddf = -T_winter * L_winter
        self.Lw=Lw
        self.beta=beta
        self.T_air=T
        self.A_air=A
        self.Tw=Tw
        self.ta_month=T_month
        self.ddt=ddt
        self.ddf=ddf
    #   update_dd_from_annual_minmax_temp()
    #-------------------------------------------------------------------
    def update_air_frost_number(self):
        # Calculating Reduced Air Frost Number (pages 280-281).
        # The reduced frost number is close 0 for long summers and close to 1 for long winters.
        self.air_frost_number = np.sqrt(self.ddf) / ( np.sqrt( self.ddf) + np.sqrt( self.ddt) )

    #   update_air_frost_number()
    #-------------------------------------------------------------------

    def update_snow_prop(self):
        # find indexes for which temp > 0 and make precip = 0
        if (self.T_air_type != 'Scalar'): # if not should stop
            #wk = np.loadtxt('examples/prec.txt', skiprows=1,unpack=False)
            precipitation_filename = self.permafrost_dir +\
                "permamodel/examples/prec.txt"
            wk = np.loadtxt(precipitation_filename, skiprows=1,unpack=False)
            t_month = wk[:,0]
            prec_month = wk[:,1]

        pos_temp_ind=np.array(np.where(self.ta_month>0))
        prec_month[pos_temp_ind]=0
        neg_temp_ind=np.array(np.where(self.ta_month<=0))

        if not pos_temp_ind.any():
        # monthly temp is always below zero
        # i.e. it constantly snows over whole year
        # the point associated with glaciaer and needs to excluded
            print 'snows constatly: stop!'

        m=np.size(neg_temp_ind)
        pp=0.5; # assume only 50% of precip change to at the beg and end of the snow season

        # this is portions of the code assumes a perfect winter season
        # needs to be used with care when there is a warm month during snow season
        if (m==1):
            s_idx=neg_temp_ind[:,0]
            e_idx=neg_temp_ind[:,m-1]
            prec_month[s_idx]=prec_month[s_idx]*pp
        else:
            s_idx=neg_temp_ind[:,0]
            e_idx=neg_temp_ind[:,m-1]
            prec_month[s_idx]=prec_month[s_idx]*pp
            prec_month[e_idx]=prec_month[e_idx]*pp

        # sum up precip to get SWE
        j=0; s=0; swe=np.zeros(m);
        for i in range(s_idx,e_idx+1):
            s=s+prec_month[i]
            swe[j]=s
            j=j+1

        #calculating snow density, depth and thermal counductivity
        r_snow=np.zeros(m); # snow density in kg/m3
        h_snow=np.zeros(m); # snow depth in m
        c_snow=np.zeros(m); # snow depth in W/mK

        rho_sn_min=200; rho_sn_max=300 # allowed min and max snow density
        tauf=0.24 # e-folding value (see Verseghy, 1991)

        s=rho_sn_min
        s=((s - rho_sn_max)*np.exp(-tauf)) + rho_sn_max
        r_snow[0] = s
        for i in range(1,m):
        # starting from month 2 tauf should be multpled by the 30 days
        # otherwise snow thermal conductivity can be low and insulate ground well enough over the season
        # usually we assume constant max snow thermal conductivity over snow season
            s=((s - rho_sn_max)*np.exp(-tauf)) + rho_sn_max
            r_snow[i] = s

        h_snow  = (swe/(r_snow*0.001))
        # snow thermal conductivity according to M. Sturm, 1997.
        c_snow = (0.138-1.01*r_snow + 3.233*(r_snow**2))*1e-6

        self.r_snow=r_snow
        self.h_snow=h_snow
        self.c_snow=c_snow

    #   update_snow_prop()
    #-------------------------------------------------------------------
    def update_surface_frost_number(self):
        # phi [scalar]: sites latitude
        # Zs [scalar]: an average winter snow thickness
        # Zss [scalar]: a damping depth in snow
        # P [scalar]: length of an annual temperature cycle
        # k [scalar]: number of winter months
        # rho_s [scalar]: density of snow [kg m-3]
        # lambda_s [scalar]: snow thermal conductivity [W m-1 C-1]
        # c_s [scalar]: snow specific heat capacity [J kg-1 C-1]
        # alpha_s [scalar]: thermal diffusivity of snow
        # Uw [scalar]: mean winter wind speed [m s-1]
        # Aplus [scalar]: temperature amplitude at the surface with snow
        # Twplus [scalar]: the mean winter surface temperature
        # DDFplus [scalar]: freezing index at the surface
        # Tplus [scalar]: mean annual tempratures at the surface
        # Fplus [scalar] : surface frost number

        rho_s=np.mean(self.r_snow)
        lambda_s=np.mean(self.c_snow)
        Zs=np.mean(self.h_snow)
        P=2*np.pi/365; # i am not sure what they mean by length of the annual temprature cycle
        # Something worthwhile discussing

        c_s=7.79*self.Tw+2115                                               #(eqn. 7)
        alpha_s=lambda_s/(c_s*rho_s)                                        #(eqn. 8)
        Zss=np.sqrt(alpha_s*P/np.pi)                                        #(eqn. 10)
        Aplus=self.A_air*np.exp(-Zs/Zss)                                    #(eqn. 9)
        Twplus=self.T_air-Aplus*np.sin(self.beta/(np.pi-self.beta))         #(eqn. 11)
        # Twplus is a mean winter surface temprature, I think, should be warmer than air temperature?
        # Here is another problem. DDFplus degree days winter should be positive.
        # The way it is written in the paper is wrong. I added a minus sign to fix it (see eqn. 2.9)
        DDFplus=-Twplus*self.Lw                                             #(eqn. 12)
        Tplus=(self.ddt-DDFplus)/365                                        #(eqn. 13)
        #Nevertheless the surface frost number is smaller than air which looks resonable to me.
        self.Fplus=np.sqrt(DDFplus)/(np.sqrt(self.ddt)+np.sqrt(DDFplus))    #(eqn. 14)
        self.Twplus=Twplus

    #   update_surface_frost_number()
    #-------------------------------------------------------------------
    def update_stefan_frost_number(self):
        # Zfplus [scalar]: the depth [m] to which forst extends
        # lambda_f [scalar]: frozen soil thermal conductivity [W m-1 C-1]
        # S [scalar]: is a const scalar factor [s d-1]
        # rho_d [scalar]: dry density of soil [kg m-3]
        # wf [scalar] : soil water content (proportion of dry weight)
        # L [scalar] : is a latent heat of fusion of water [J kg-1]

        lambda_f=1.67 # some dummy thermal conductivity
        # https://shop.bgs.ac.uk/GeoReports/examples/modules/C012.pdf
        sec_per_day=86400
        rho_d=2.798  # dry density of silt
        wf=0.4       # tipical for silty soils
        denominator=rho_d*wf*self.Lf
        self.Zfplus=np.sqrt(2*lambda_f*sec_per_day*np.abs(self.Twplus)*self.Lw/denominator)               #(eqn. 15)
        print 'Zfplus=',self.Zfplus

        # assuming 3 soil layers with thickness 0.25, 0.5 and 1.75 m
        # and thermal conductivities 0.08, 1.5, and 2.1
        soil_thick=np.array([0.25, 0.5 , 1.75])
        lambda_b=np.array([0.08, 1.5, 2.1])
        # resistivity R
        R=soil_thick/lambda_b
        QL=self.Lf/1000 # volumetric latent heat, 1000 s a density of water
        #partial freezing thawing index DD
        DD=np.zeros(3)
        Z=np.zeros(3)
        DD[0]=0.5*QL*soil_thick[0]*R[0]/sec_per_day
        S=0;
        for i in range(1,3):
            S=R[i]+S
            DD[i]=(S+0.5*R[i-1])*QL*soil_thick[i]/sec_per_day
            #The depth of the frost thaw penetration
        S=0; Z_tot=0
        for i in range(0,3):
            #The depth of the frost thaw penetration
            Z[i]=np.sqrt(2*lambda_b[i]*sec_per_day*DD[i]/QL + lambda_b[i]**2*S**2) \
                - lambda_b[i]*S
            S=R[i]+S
            Z_tot=Z_tot + Z[i]

        self.Z_tot=Z_tot
        self.stefan_number = np.sqrt(self.Fplus) / ( np.sqrt( self.Fplus) + np.sqrt( self.Z_tot) )

    #   update_stefan_frost_number()
    #-------------------------------------------------------------------

    def update_ALT(self):

        #---------------------------------------------------------
        #       coming up
        #--------------------------------------------------
        print 'OK'

    #   update_ALT()
    #-------------------------------------------------------------------
    def update_ground_temperatures(self):
        # This method does not update temps instead it does frost numbers
        self.update_dd()
        self.update_air_frost_number()
        self.update_snow_prop()
        self.update_surface_frost_number()
        self.update_stefan_frost_number()

    #   update_ground_temperatures()
    #-------------------------------------------------------------------
    def close_input_files(self):

        if (self.T_air_type     != 'Scalar'): self.T_air_unit.close()
        if (self.A_air_type     != 'Scalar'): self.A_air_unit.close()

    #   close_input_files()
    #-------------------------------------------------------------------

