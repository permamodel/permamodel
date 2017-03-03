#!/bin/bash

# The script installs the permamodel repository in the directory 'pm_env'
# or in a directory supplied as a command line argument

# Get the directory name, either default or from the only command line argument
if [ "$#" -eq 0 ]; then
	dirname=pm_env
elif [ "$#" -eq 1 ]; then
	dirname=$1
else 
	echo "Too many arguments"
	echo "Usage: install_pm [directory_name_to_create_and_install_to]"
	exit
fi

echo "dirname: $dirname "

# Verify that the directory doesn't exist, otherwise quit
if [ -e $dirname ]; then
	echo "$dirname already exists, quitting"
	exit
fi

# Create and enter the directory
mkdir $dirname && cd $dirname
echo "$dirname successfully created"

# Get the latest miniconda repository
# First need to know whether on a Max or a PC
echo "Downloading latest miniconda repository..."
this_os=`uname`
echo "this_os: ...$this_os..."
if [ $this_os = 'Linux' ]; then
  curl https://repo.continuum.io/miniconda/Miniconda2-latest-Linux-x86_64.sh -o miniconda.sh
elif [ $this_os = 'Darwin' ]; then
  curl https://repo.continuum.io/miniconda/Miniconda2-latest-MacOSX-x86_64.sh -o miniconda.sh
else
	echo "Operating system not recognized as mac or linux: $this_os"
	exit
fi

# Install the latest miniconda repository
echo "Installing miniconda in $dirname..."
#bash ./miniconda.sh -b -p $(pwd)/conda
bash ./miniconda.sh -b -p $PWD/conda

# Set up python environment
#PATH=$(pwd)/conda/bin:$PATH
PATH=$PWD/conda/bin:$PATH

# Install necessary and desired conda repository
# Required, although these would be accomplished via pip from setup.py install
declare -a module_list=("nose" "numpy" "scipy" "netCDF4" "affine" "pyyaml")
# Potentially helpful packages: ipython
# Eventually required packages: gdal, pyproj, osr
for module in "${module_list[@]}"; do
  echo "Installing $module via conda..."
  conda install -y $module
done

# Create "use me" script
echo "Creating script to allow user to set up this environment for use..."
scriptname="$(pwd)/use_this_env"
echo "scriptname is: $scriptname"
echo "PATH=\$(pwd)/conda/bin:\$PATH" >> $scriptname
#echo "echo \"PATH is now \$PATH\"" >> $scriptname
echo "echo \$(pwd)/conda/bin/ has been prepended to the PATH" >> $scriptname
echo "From $(pwd), run \". ./$scriptname\" to set as the python environment"

# Note: This is too clever right now, because the install directory may 
# not be easily discoverable
#declare -a repository_list=(\
#  "https://github.com/permamodel/permamodel.git"
#	"https://github.com/csdms/bmi-tester.git"\
#	)
#
#for repo in "${repository_list[@]}"; do
#	echo "Installing repository from: $repo"
#	git clone $repo
# cd $the_repo_install_directory
# python setup.py develop
# cd ..
#done

# Get and install the permamodel repository
echo "Get and install permamodel..."
git clone https://github.com/permamodel/permamodel.git
cd permamodel
python setup.py develop
cd ..

# Get and install the bmi-tester repository
echo "Get and install bmi-tester..."
git clone https://github.com/csdms/bmi-tester.git
cd bmi-tester
python setup.py develop
cd ..

# Get and install the cruAKtemp repository
echo "Get and install cruAKtemp..."
git clone https://github.com/permamodel/cruAKtemp.git
cd cruAKtemp
python setup.py develop
cd ..

# Print final message and exist
echo " "
echo "Miniconda environment for permamodel successfully installed in $(pwd)"
echo " "
echo "Next steps"
echo "  cd $dirname"
echo "  . ./use_this_env"
echo "  cd permamodel"
echo "    nosetests -x"
echo "    bmi-tester permamodel.components.bmi_frost_number.BmiFrostnumberMethod"
echo "If you want to add a non-standard conda package--e.g. rednose for nosetests:"
echo "  conda install -y anaconda-client"
echo "  conda install -y --channel https://conda.anaconda.org/blaze rednose"

exit

# Other things to add?
#  Run nose tests to verify installation
#  Run bmi-tester to verify BMI compliance
