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
  curl http://repo.continuum.io/miniconda/Miniconda-latest-Linux-x86_64.sh -o miniconda.sh
elif [ $this_os = 'Darwin' ]; then
  curl http://repo.continuum.io/miniconda/Miniconda-latest-MacOSX-x86_64.sh -o miniconda.sh
else
	echo "Operating system not recognized as mac or linux: $this_os"
	exit
fi

# Install the latest miniconda repository
echo "Installing miniconda in $dirname..."
#bash ./miniconda.sh -b -p $(pwd)/conda  # This fails on some shells
bash ./miniconda.sh -b -p $PWD/conda

# Set up python environment
PATH=$PWD/conda/bin:$PATH

