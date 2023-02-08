## Copyright (c) 2009-2010, Scott D. Peckham
## January 9, 2009
## April 29, 2009
## May 2010, Changed var_types from {0,1,2,3} to
#            {'Scalar', 'Time_Series', 'Grid'}, etc.
# -------------------------------------------------------------------

#  open_file()
#  read_next()
#  read_scalar()
#  read_grid()
#  close_file()
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
# -------------------------------------------------------------------
from __future__ import print_function

import os.path

import numpy
from numpy import *


# -------------------------------------------------------------------
def open_file(var_type, input_file):
    # -----------------------------------------------------
    # Note:  This method's name cannot be "open" because
    #        it calls Python's built-in "open()" method.
    # -----------------------------------------------------
    # print 'var_type   =', var_type
    # print 'input_file =', input_file

    # --------------------------------------------
    # A scalar input value was provided already
    # --------------------------------------------
    file_unit = None
    if var_type.lower() == "scalar":
        return file_unit
    if input_file == "":
        print("ERROR in model_input.open_file():")
        print("    Input file is null string.")
        # print '    variable type =' + var_type
        return file_unit

    # ----------------------------------
    # Does input file exist locally ?
    # ----------------------------------
    if not (os.path.exists(input_file)):
        print("ERROR in model_input.open_file():")
        print("    Could not find input file =")
        print("    " + input_file)
        # print '   ' + input_file
        return file_unit

    if var_type.lower() == "time_series":
        # -----------------------------------------
        # Input file contains a time series and
        # is ASCII text with one value per line.
        # -----------------------------------------
        file_unit = open(input_file, "r")
    else:
        # --------------------------------------------
        # Input file contains a grid or grid stack
        # as row-major, binary file with no header.
        # --------------------------------------------
        file_unit = open(input_file, "rb")

    return file_unit


#   open_file()
# -------------------------------------------------------------------
def read_next(file_unit, var_type, rti, dtype="Float32", factor=1.0):
    # -------------------------------------------------------
    # (5/7/09) Allow "dtype" to be given using RTI types.
    # -------------------------------------------------------
    rti_types = ["BYTE", "INTEGER", "LONG", "FLOAT", "DOUBLE"]
    if dtype.upper() in rti_types:
        dtype_map = {
            "BYTE": "uint8",
            "INTEGER": "int16",
            "LONG": "int32",
            "FLOAT": "float32",
            "DOUBLE": "float64",
        }
        dtype = dtype_map[dtype]

    if var_type.lower() == "scalar":
        # -------------------------------------------
        # Scalar value was entered by user already
        # -------------------------------------------
        data = None
    elif var_type.lower() == "time_series":
        # ----------------------------------------------
        # Time series: Read scalar value from file.
        # File is ASCII text with one value per line.
        # ----------------------------------------------
        data = read_scalar(file_unit, dtype)
    elif var_type.lower() in ["grid", "grid_sequence"]:
        # --------------------------------------
        # Single grid: Read grid from file
        # read_grid() accounts for byte order
        # ----------------------------------------------
        # NB!  grid_type argument allows DEM to be
        # read for GW vars, which might not be FLOAT'
        # ----------------------------------------------
        data = read_grid(file_unit, rti, dtype)
    else:
        raise RuntimeError('No match found for "var_type".')
        return None

    # ---------------------------------------------
    # Multiply by a conversion or scale factor ?
    # ---------------------------------------------
    if (factor != 1) and (data != None):
        data = data * factor

    # -----------------------------------------------------
    # Values must usually be read from file as FLOAT32
    # but then need to be returned as FLOAT64. (5/17/12)
    # But numpy.float64( None ) = NaN. (5/18/12)
    # -----------------------------------------------------
    if data == None:
        return
    else:
        return numpy.float64(data)
    ## return data


#   read_next()
# -------------------------------------------------------------------
def read_next_modified(file_unit, var_type, dtype="Float32", factor=1.0):
    # -------------------------------------------------------
    # (5/7/09) Allow "dtype" to be given using RTI types.
    # (4/21/16) Elchin Jafarov introduced this function b/c
    # he was not sure how to deal with rti in the original function
    # this version foes not have grid choice
    # -------------------------------------------------------
    rti_types = ["BYTE", "INTEGER", "LONG", "FLOAT", "DOUBLE"]
    if dtype.upper() in rti_types:
        dtype_map = {
            "BYTE": "uint8",
            "INTEGER": "int16",
            "LONG": "int32",
            "FLOAT": "float32",
            "DOUBLE": "float64",
        }
        dtype = dtype_map[dtype]

    if var_type.lower() == "scalar":
        # -------------------------------------------
        # Scalar value was entered by user already
        # -------------------------------------------
        data = None
    elif var_type.lower() == "time_series":
        # ----------------------------------------------
        # Time series: Read scalar value from file.
        # File is ASCII text with one value per line.
        # ----------------------------------------------
        data = read_scalar(file_unit, dtype)

    else:
        raise RuntimeError('No match found for "var_type".')
        return None

    # ---------------------------------------------
    # Multiply by a conversion or scale factor ?
    # ---------------------------------------------
    if (factor != 1) and (data != None):
        data = data * factor

    # -----------------------------------------------------
    # Values must usually be read from file as FLOAT32
    # but then need to be returned as FLOAT64. (5/17/12)
    # But numpy.float64( None ) = NaN. (5/18/12)
    # -----------------------------------------------------
    if data == None:
        return
    else:
        return numpy.float64(data)
    ## return data


#   read_next_modified()
# -------------------------------------------------------------------
##def open_file(var_type, input_file):
##
##    #-----------------------------------------------------
##    # Note:  This method's name cannot be "open" because
##    #        it calls Python's built-in "open()" method.
##    #-----------------------------------------------------
##    # print 'var_type   =', var_type
##    # print 'input_file =', input_file
##
##    #--------------------------------------------
##    # A scalar input value was provided already
##    #--------------------------------------------
##    file_unit = None
##    if (var_type == 0):  return file_unit
##    if (input_file == ''):
##        print 'ERROR in model_input.open_file():'
##        print '    Input file is null string.'
##        # print '    variable type =' + var_type
##        return file_unit
##
##    #----------------------------------
##    # Does input file exist locally ?
##    #----------------------------------
##    if not(os.path.exists(input_file)):
##        print 'ERROR in model_input.open_file():'
##        print '    Could not find input file ='
##        print '    ' + input_file
##        # print '   ' + input_file
##        return file_unit
##
##    if (var_type == 1):
##        #-----------------------------------------
##        # Input file contains a time series and
##        # is ASCII text with one value per line.
##        #-----------------------------------------
##        file_unit = open(input_file, 'r')
##    elif (var_type > 1):
##        #--------------------------------------------
##        # Input file contains a grid or grid stack
##        # as row-major, binary file with no header.
##        #--------------------------------------------
##        file_unit = open(input_file, 'rb')
##
##    return file_unit
##
###   open_file()
###-------------------------------------------------------------------
##def read_next(file_unit, var_type, rti, \
##              dtype='Float32', factor=1.0):
##
##    #-------------------------------------------------------
##    # (5/7/09) Allow "dtype" to be given using RTI types.
##    #-------------------------------------------------------
##    rti_types = ['BYTE','INTEGER','LONG','FLOAT','DOUBLE']
##    if (dtype.upper() in rti_types):
##        dtype_map = {'BYTE':'uint8', 'INTEGER':'int16',
##                     'LONG':'int32',
##                     'FLOAT':'float32', 'DOUBLE':'float64'}
##        dtype = dtype_map[dtype]
##
##
##    if (var_type == 0):
##        #-------------------------------------------
##        # Scalar value was entered by user already
##        #-------------------------------------------
##        data = None
##    elif (var_type == 1):
##        #----------------------------------------------
##        # Time series: Read scalar value from file.
##        # File is ASCII text with one value per line.
##        #----------------------------------------------
##        data = read_scalar(file_unit, dtype)
##    elif (var_type in [2,3]):
##        #--------------------------------------
##        # Single grid: Read grid from file
##        # read_grid() accounts for byte order
##        #----------------------------------------------
##        # NB!  grid_type argument allows DEM to be
##        # read for GW vars, which might not be FLOAT'
##        #----------------------------------------------
##        data = read_grid(file_unit, rti, dtype)
##    else:
##        raise RuntimeError('No match found for "var_type".')
##        return None
##
##    if (factor != 1) and (data != None):
##        data = (data * factor)
##    return data
##
###   read_next()
# -------------------------------------------------------------------
def read_scalar(file_unit, dtype="Float32"):
    # -------------------------------------------------
    # Note:  Scalar values are read from text files.
    # -------------------------------------------------

    # -------------------------------
    # Recheck if any of these work
    # -------------------------------
    # scalar = fromfile(file_unit, count=1, dtype=dtype)
    # scalar = fromfile(file_unit, count=1, dtype=dtype, sep="\n")
    # scalar = fromfile(file_unit, count=1, dtype=dtype, sep=" ")
    # scalar = loadtxt(file_unit, dtype=dtype)

    line = file_unit.readline()
    line = line.strip()
    ##    print 'line =', line
    if line != "":
        words = line.split()
        ##        print 'type(words)   =', type(words)
        ##        print 'size(words)   =', size(words)
        ##        print 'words[0]      =', words[0]
        ##        print 'type(words[0] =', type(words[0])

        # ----------------------------------------
        # This worked in Python version at home
        # ----------------------------------------
        ## scalar = numpy.float32(words[0])
        scalar = numpy.float32(eval(words[0]))
    else:
        scalar = None

    return scalar


#   read_scalar()
# -------------------------------------------------------------------
def read_grid(file_unit, rti, dtype="Float32"):
    # -------------------------------------------
    # Note:  Grids are read from binary files.
    #        Return "None" if end of file.
    # -------------------------------------------
    file_size = os.path.getsize(file_unit.name)
    file_pos = file_unit.tell()
    END_OF_FILE = file_pos == file_size
    if END_OF_FILE:
        grid = None
        return grid

    grid = fromfile(file_unit, count=rti.n_pixels, dtype=dtype)
    grid = reshape(grid, (rti.nrows, rti.ncols))
    if rti.SWAP_ENDIAN:
        grid.byteswap(True)
    return grid


#   read_grid()
# -------------------------------------------------------------------
def close_file(file_unit):
    # ----------------------------------------------------------
    # Notes:  Each process module closes its own input files
    #         with a method called "close_input_files()."
    #         This function is included for completeness and
    #         may be used instead.
    # ----------------------------------------------------------
    file_unit.close()


#   close_file()
# -------------------------------------------------------------------
