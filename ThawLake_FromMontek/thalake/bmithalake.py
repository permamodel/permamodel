"""Basic Model Interface (BMI) for the Diffusion model."""

import numpy as np
from basic_modeling_interface import Bmi
from thalake import thaLakeModel

class bmithaLakeModel(Bmi):
	"""BMI for the Diffusion model."""

	_name = 'thalake permafrost Component"'
	_input_var_names = ('depth_in_meters',)
	_output_var_names = ('depth_in_meters',)

	def __init__(self):
		"""Create a Diffusion model that's ready for initialization."""
		self._model = None
		self._values = {}
		self._var_units = {}
		self._grids = {}
		self._grid_type = {}
		self.ngrids = 0

	
	def initialize(self, filename=None):
		self._model = thaLakeModel()
		#self._model.initialize(filename=cfg_file)	
		#self._model = thaLakeModel(config_file=filename)
		#self._values = {
		#	'depth_in_meters': self._model.depth,
		#}
		self._var_units = {
			'depth_in_meters': 'm'
		}
		self._grids = {
			0: ['plate_surface__temperature']
		}
		self._grid_type = {
			0: 'uniform_rectilinear_grid'
		}
		print ("opening file")
	
	def update(self):
		"""Advance model by one time step."""
		self._model.advance()

	def update_frac(self, time_frac):
		time_step = self.get_time_step()
		self._model.dt = time_frac * time_step
		self.update()
		self._model.dt = time_step

	def update_until(self, then):
		n_steps = (then - self.get_current_time()) / self.get_time_step()

		for _ in range(int(n_steps)):
			self.update()

		if (n_steps - int(n_steps)) > 0.0:
			self.update_frac(n_steps - int(n_steps))

	def finalize(self):
		"""Finalize model."""
		self._model = None
		#assert(self._model.status == 'initialized')

	def get_var_type(self, var_name):
		return str(self.get_value(var_name).dtype)

	def get_var_units(self, var_name):
		return self._var_units[var_name]

	def get_var_nbytes(self, var_name):
		return self.get_value(var_name).nbytes

	def get_var_grid(self, var_name):
		for grid_id, var_name_list in self._grids.items():
			if var_name in var_name_list:
				return grid_id

	def get_grid_rank(self, grid_id):
		return len(self.get_grid_shape(grid_id))

	def get_grid_size(self, grid_id):
		size = self.bmi.get_grid_size(grid)
		return str(size)
	
	def get_value_ref(self, var_name):	
		return self._values[var_name].reshape(-1)

	def get_value(self, var_name):
		return self.get_value_ref(var_name).copy()
	
	
	def set_value(self, var_name, src):
		val = self.get_value_ref(var_name)
		val[:] = src
	
	def get_component_name(self):
		"""Name of the component."""
		return self._name

	def get_input_var_names(self):
		"""Get names of input variables."""
		return self._input_var_names

	def get_output_var_names(self):
		"""Get names of output variables."""
		return self._output_var_names

	def get_grid_shape(self, grid_id):
		"""Number of columns and rows of uniform rectilinear grid."""
		return (self._model.ny, self._model.nx)

	def get_grid_spacing(self, grid_id):
		"""Spacing of columns and rows of uniform rectilinear grid."""
		assert_true(grid_id < self.ngrids)
		return np.array([1, 1], dtype='float32')
		#return (self._model.dy, self._model.dx)

	def get_grid_origin(self, grid_id):
		"""Origin of uniform rectilinear grid."""
		return (0.0, 0.0)

	def get_grid_type(self, grid_id):
		"""Type of grid."""
		return self._grid_type[grid_id]

	def get_start_time(self):
		"""Start time of model."""
		return 0.0

	def get_end_time(self):
		"""End time of model."""
		return np.finfo('d').max

	def get_current_time(self):
		"""Current time of model."""
		return self._model.time

	def get_time_step(self):
		"""Time step of model."""
		return self._model.dt

	def get_time_units(self):
		"""Time units of model."""
		return '-'
