
You can use the bmi-tester package to test BMI compliance

To get this code:
  git clone https://github.com/csdms/bmi-tester.git

Then:
  python setup.py install
	or
	python setup.py develop

	from within the resulting bmi-tester/ directory

To use, go to the permamodel root directory (which gets installed when
the permamodel repository is cloned to your machine):

  bmi-tester --infile=./permamodel/examples/Frostnumber_example_singlesite_singleyear.cfg permamodel.components.bmi_frost_number.BmiFrostnumberMethod


  [OLD] bmi-tester --infile=./permamodel/examples/Frostnumber_example_singlesite_singleyear.cfg permamodel.components.frost_number.FrostnumberMethod

