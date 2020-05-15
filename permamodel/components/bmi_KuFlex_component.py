# -*- coding: utf-8 -*-
"""  Kudryavtsev model by Anisimov et al. (1997). 

Anisimov, O. A., Shiklomanov, N. I., & Nelson, F. E. (1997).
Global warming and active-layer thickness: results from transient general circulation models. 
Global and Planetary Change, 15(3-4), 61-77. DOI:10.1016/S0921-8181(97)00009-X

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

import warnings
warnings.filterwarnings("ignore",category =RuntimeWarning) 

import numpy as np
from permamodel.utils import model_input
from permamodel.components import perma_base
from permamodel.components import KuFlex_method
#from permamodel.components.perma_base import *
#from permamodel.tests import examples_directory
import os


"""
class BmiKuMethod( perma_base.PermafrostComponent ):
    _thisname = 'this name'
"""

class BmiKuFlexMethod( perma_base.PermafrostComponent ):

    """ Implement the Ku model """

    # Set up the name of this permafrost module
    _name = 'KuFlex module'

    """ Note: all of these are defined below instead!
    # Indicate the CSDMS standard names of input and output variables
    _input_var_names = ('land_surface_air__temperature',
                        'land_surface__latitude',
                        'land_surface__longitude',
                       )
                       # other standard names that might apply?
                       # land_surface__temperature
                       # model__initial_time_step
                       # model__max_allowed_time_step
                       # model__min_allowed_time_step
                       # model__run_time
                       # model__spinup_time
                       # model__start_time
                       # model__stop_time
                       # model__time
                       # model__time_step
                       # model__time_step_count
    _output_var_names = ('frost_number_air',
                         'frost_number_surface',
                         'frost_number_stefan',
                        )
                        # other standard names that might apply?
                        # soil_permafrost__thickness
                        # soil_permafrost_top__depth
                        # soil_permafrost_bottom__depth
    """

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

    # This used to be [...] instead of (...)
    _input_var_names = (
            
        'datetime__start',
        'datetime__end',
        
        'atmosphere_bottom_air__temperature',
        'atmosphere_bottom_air__temperature_amplitude',
        
        'snowpack__depth',
        'snowpack__density',
        'snow__thermal_conductivity',
        'snow__volume-specific_isochoric_heat_capacity',
        
        'water__volume_latent_fusion_heat',
        'soil-frozen__volume-specific_isochoric_heat_capacity',
        'soil-thaw__volume-specific_isochoric_heat_capacity',
        'soil-frozen__thermal_conductivity',
        'soil-thaw__thermal_conductivity',
        
        'vegetation__Hvgf',
        'vegetation__Hvgt',
        'vegetation__Dvf',
        'vegetation__Dvt' )

    _output_var_names = (
        'soil_surface__temperature',
        'soil_surface__temperature_amplitude',
        'soil__temperature',                                  # Tps
        'soil__active_layer_thickness' )                      # Zal

    _var_name_map = {
        'datetime__start':                                    'start_year',
        'datetime__end':                                      'end_year',
        
        'atmosphere_bottom_air__temperature':                 'T_air',
        'atmosphere_bottom_air__temperature_amplitude':       'A_air',
        
        'snowpack__depth':                                    'h_snow',
        'snowpack__density':                                  'rho_snow',
        'snow__thermal_conductivity':                         'k_snow',
        'snow__volume-specific_isochoric_heat_capacity':      'c_snow',
        
        'water__volume_latent_fusion_heat':                   'L',
        'soil-frozen__volume-specific_isochoric_heat_capacity':'Cf',
        'soil-thaw__volume-specific_isochoric_heat_capacity': 'Ct',
        'soil-frozen__thermal_conductivity':                  'Kf',
        'soil-thaw__thermal_conductivity':                    'Kt',
        
        'vegetation__Hvgf':                                   'Hvgf',
        'vegetation__Hvgt':                                   'Hvgt',
        'vegetation__Dvf':                                    'Dvf',
        'vegetation__Dvt':                                    'Dvt' ,
        

        'soil_surface__temperature':                          'Tgs',
        'soil_surface__temperature_amplitude':                'Ags',        
        'soil__temperature':                                  'Tps',
        'soil__active_layer_thickness':                       'Zal'}


    _var_units_map = {
        # These are the links to the model's variables' units
        'datetime__start':                                    'year',
        'datetime__end':                                      'year',
        
        'atmosphere_bottom_air__temperature':                 'deg_C',
        'atmosphere_bottom_air__temperature_amplitude':       'deg_C',
        
        'snowpack__depth':                                    'm',
        'snowpack__density':                                  'kg m-3',
        'snow__thermal_conductivity':                         'W m-1 K-1',
        'snow__volume-specific_isochoric_heat_capacity':      'J m-3 K-1',
        
        'water__volume_latent_fusion_heat':                   'J m-3',
        'soil-frozen__volume-specific_isochoric_heat_capacity':'J m-3 K-1',
        'soil-thaw__volume-specific_isochoric_heat_capacity': 'J m-3 K-1',
        'soil-frozen__thermal_conductivity':                  'W m-1 K-1',
        'soil-thaw__thermal_conductivity':                    'W m-1 K-1',
        
        'vegetation__Hvgf':                                   'm',
        'vegetation__Hvgt':                                   'm',
        'vegetation__Dvf':                                    'm2 s-1',
        'vegetation__Dvt':                                    'm2 s-1'  ,

        'soil_surface__temperature':                          'deg_C',
        'soil_surface__temperature_amplitude':                'deg_C',         
        'soil__temperature':                                  'deg_C',
        'soil__active_layer_thickness':                       'm'}

    #-------------------------------------------------------------------
    def __init__(self):
        self._model = None
        self._values = {}
        self._var_units = {}
        self._grids = {}
        self._grid_type = {}

    def initialize(self, cfg_file=None):
        
        self._model = KuFlex_method.KuFlex_method()
        
        self._name = "Permamodel Ku Component"
        self._model.initialize(cfg_file=cfg_file)
        
        # make 2 vars to store each results and used for write out.
        n_lat  = self._model.grid_shape[0]
        n_lon  = self._model.grid_shape[1]
        n_time = self._model.end_year-self._model.start_year+1
                
        self.output_alt = np.zeros((n_time,n_lat,n_lon))*np.nan;
        self.output_tps = np.zeros((n_time,n_lat,n_lon))*np.nan;
        self.output_tgs = np.zeros((n_time,n_lat,n_lon))*np.nan;
        self.output_ags = np.zeros((n_time,n_lat,n_lon))*np.nan;
        
        # Verify that all input and output variable names are in the
        # variable name and the units map
        for varname in self._input_var_names:
            assert(varname in self._var_name_map)
            assert(varname in self._var_units_map)
            #print("Input var %s is in the name map and the units map"\
            #      % varname)
        for varname in self._output_var_names:
            assert(varname in self._var_name_map)
            assert(varname in self._var_units_map)
            
#        self._model.cont = -1

        gridnumber = 0
        for varname in self._input_var_names:
            self._grids[gridnumber] = varname
            self._grid_type[gridnumber] = 'uniform_rectilinear'
            #self._grid_type[gridnumber] = 'scalar'
            gridnumber += 1
        for varname in self._output_var_names:
            self._grids[gridnumber] = varname
            self._grid_type[gridnumber] = 'uniform_rectilinear'
            #self._grid_type[gridnumber] = 'scalar'
            gridnumber += 1        

        self._values = {
                
        # These are the links to the model's variables and
        # should be consistent with _var_name_map 

            'datetime__start':                                             self._model.start_year,
            'datetime__end':                                               self._model.end_year,
            
            'atmosphere_bottom_air__temperature':                          self._model.T_air,
            'atmosphere_bottom_air__temperature_amplitude':                self._model.A_air,
            
            'snowpack__depth':                                             self._model.h_snow,
            'snowpack__density':                                           self._model.rho_snow,
            'snow__thermal_conductivity':                                  self._model.k_snow,
            'snow__volume-specific_isochoric_heat_capacity':               self._model.c_snow,
            
            'water__volume_latent_fusion_heat':                            self._model.L,
            'soil-frozen__volume-specific_isochoric_heat_capacity':        self._model.cf_soil,
            'soil-thaw__volume-specific_isochoric_heat_capacity':          self._model.ct_soil,
            'soil-frozen__thermal_conductivity':                           self._model.kf_soil,
            'soil-thaw__thermal_conductivity':                             self._model.kt_soil,
            
            'vegetation__Hvgf':                                            self._model.Hvgf,
            'vegetation__Hvgt':                                            self._model.Hvgt,
            'vegetation__Dvf':                                             self._model.Dvf,
            'vegetation__Dvt':                                             self._model.Dvt,

            'soil_surface__temperature':                                   self._model.Tgs,
            'soil_surface__temperature_amplitude':                         self._model.Ags,              
            'soil__temperature':                                           self._model.Tps,
            'soil__active_layer_thickness':                                self._model.Zal}
        
    def get_attribute(self, att_name):

        try:
            return self._att_map[ att_name.lower() ]
        except:
            print('###################################################')
            print(' ERROR: Could not find attribute: ' + att_name)
            print('###################################################')
            print(' ')

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

    def get_var_location(self, long_var_name):
        return "node"

    def update(self):
#        self._model.update(self._model.dt)
            # Ensure that we've already initialized the run
#        assert(self._model.status == 'initialized')
        
        self._model.update()
       
        self._values['soil__active_layer_thickness']        = self._model.Zal
        self._values['soil__temperature']                   = self._model.Tps
        self._values['soil_surface__temperature']           = self._model.Tgs
        self._values['soil_surface__temperature_amplitude'] = self._model.Ags
        
    def update_frac(self, time_fraction):
        time_step = self.get_time_step()
        self._model.dt = time_fraction * time_step
        self.update()
        self._model.dt = time_step
    
    def update_until(self, then):
        
        assert then <= self.get_end_time(), 'Exceed time end'
            
        n_steps = (then - self.get_current_time()) / self.get_time_step()
        for _ in range(int(n_steps)):
            self.update()
        self.update_frac(n_steps - int(n_steps))

    def finalize(self):

        # Finish with the run
        self._model.status = 'finalizing'  # (OpenMI)

        # Close the input files
        self._model.close_input_files()   # Close any input files
        self._model.close_output_files()

        # Done finalizing  
        self._model.status = 'finalized'  # (OpenMI)
        
    def get_start_time(self):
        return 0.0

    def get_current_time(self):
        return float(self._model.year - self._model.start_year)
    
    def get_time(self):
        return self.get_current_time()

    def get_end_time(self):
        return self._model.end_year - self._model.start_year + 1.0

    # ----------------------------------
    # Functions added to pass bmi-tester
    # ----------------------------------
    def get_grid_type(self, grid_number):
        return self._grid_type[grid_number]

    def get_time_step(self):
        return self._model.dt

    # Note: get_value_ref() copied from bmi_heat.py
    def get_value_ref(self, var_name):
        """Reference to values.

        Parameters
        ----------
        var_name : str
            Name of variable as CSDMS Standard Name.

        Returns
        -------
        array_like
            Value array.
        """
        return getattr(self._model, self._var_name_map[var_name]) 

    def set_value(self, var_name, new_var_values):
        assert np.size(new_var_values) == 1 or np.shape(new_var_values) == self._model.grid_shape, 'inconsistent array shape' 
        if np.size(new_var_values) ==1:
            new_var_values = np.zeros(self._model.grid_shape) + new_var_values;
        setattr(self._model, self._var_name_map[var_name], new_var_values)
#        self._values[var_name] = new_var_values

    def set_value_at_indices(self, var_name, new_var_values, indices):
        self.get_value_ref(var_name).flat[indices] = new_var_values

    def get_var_itemsize(self, var_name):
        return np.asarray(self.get_value_ref(var_name)).flatten()[0].nbytes

    def get_value_at_indices(self, var_name, indices):
        return self.get_value_ref(var_name).take(indices)

    def get_var_nbytes(self, var_name):
        return np.asarray(self.get_value_ref(var_name)).nbytes

    def get_value(self, var_name, out=None):
        """Copy of values.

        Parameters
        ----------
        var_name : str
            Name of variable as CSDMS Standard Name.

        Returns
        -------
        array_like
            Copy of values.
        """
        # Original version: from bmi_heat.py
        #return self.get_value_ref(var_name).copy()

        # Version to convert to numpy array for bmi-tester compliance
        # Note: converting to np arrays on the fly here
        # Note: float values don't have a copy() function
        #try:
        #    return np.array(self.get_value_ref(var_name).copy())
        #except AttributeError:
        #    return np.array(self.get_value_ref(var_name))
        if out is None:
            out = self.get_value_ref(var_name).copy()
        else:
            out[...] = self.get_value_ref(var_name)
        return out

    def get_var_type(self, var_name):
        """Data type of variable.

        Parameters
        ----------
        var_name : str
            Name of variable as CSDMS Standard Name.

        Returns
        -------
        str
            Data type.
        """
        return str(self.get_value_ref(var_name).dtype)

    def get_component_name(self):
        return self._name

    # Copied from bmi_heat.py
    def get_var_grid(self, var_name):
        """Grid id for a variable.

        Parameters
        ----------
        var_name : str
            Name of variable as CSDMS Standard Name.

        Returns
        -------
        int
            Grid id.
        """
        for grid_id, var_name_list in self._grids.items():
            if var_name in var_name_list:
                return grid_id
            
    def get_grid_shape(self, grid_id):
        
        return self._model.grid_shape
    
    def get_grid_spacing(self, grid_id):
        """Get distance between nodes of the computational grid.

        Parameters
        ----------
        grid_id : int
          A grid identifier.

        Returns
        -------
        array_like
          The grid spacing.

        See Also
        --------
        bmi.vars.BmiVars.get_var_grid : Obtain a `grid_id`.

        """
        return [1,1]

    def get_grid_origin(self, grid_id):
        """Get coordinates for the origin of the computational grid.

        Parameters
        ----------
        grid_id : int
          A grid identifier.

        Returns
        -------
        array_like
          The coordinates of the lower left corner of the grid.

        See Also
        --------
        bmi.vars.BmiVars.get_var_grid : Obtain a `grid_id`.

        """
        return [0, 0]

    def get_grid_size(self, grid_id):
        """Size of grid.

        Parameters
        ----------
        grid_id : int
            Identifier of a grid.

        Returns
        -------
        int
            Size of grid.

        """
        
        return self._model.n_grids

    def get_grid_rank(self, var_id):
        """Rank of grid.

        Parameters
        ----------
        grid_id : int
            Identifier of a grid.

        Returns
        -------
        int
            Rank of grid.
        """
        return 2