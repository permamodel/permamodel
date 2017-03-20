# bmi-gipl_v0.1
Fortran bindings for the Basic Modeling Interface.
bmi/bmi_impl.f90 contains all required subroutines. 
Currently only irf_test.f90 works. All input files are in examples folder. The gipl_cfg.cfg contains an information about input files location and parameters setup. If there no cfg file then it load the default assuming 'examples' folder exist and has all the necessary input files. 

Build
-----
To build the BMI Fortran-bindings and tests,

    $ mkdir _build && cd _build
    $ cmake ../ 
    $ make

Run some simple tests with,

    $ make test

or run IRF file
    
    $cd testing
    $./irf_test

When finished the result files will be stored in 'examples/out'.  