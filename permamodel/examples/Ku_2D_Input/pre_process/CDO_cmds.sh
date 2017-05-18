rm -rf ../*.nc

# Calculate annual mean of air temperature:

cdo -b 64 chname,t2m,ta -timselmean,12 -selname,t2m -subc,273.15 Raw_ERA_Dataset.nc ../ta.nc

# Calculate annual amplitude of air temperature:

cdo -b 64 chname,t2m,aa -sub -timselmax,12 -selname,t2m Raw_ERA_Dataset.nc -timselmin,12 -selname,t2m Raw_ERA_Dataset.nc ../aa.nc

# Calculate annual mean of snow depth:

cdo -b 64 chname,sd,snd -timselmean,12 -selname,sd Raw_ERA_Dataset.nc ../snd.nc

# Calculate annual mean of snow density:

cdo -b 64 chname,rsn,rsn -timselmean,12 -selname,rsn Raw_ERA_Dataset.nc ../rsn.nc

# Calculate annual max of soil water (summer):

cdo -b 64 chname,swvl1,vwc -timselmax,12 -selname,swvl1 Raw_ERA_Dataset.nc ../vwc.nc

# Make a dump vegetation file (all value are set to ZERO):

#cdo -b 64 chname,swvl1,hvgf -mulc,0 -timmax -selname,swvl1 Raw_ERA_Dataset.nc ../hvgf.nc
#cdo -b 64 chname,swvl1,hvgt -mulc,0 -timmax -selname,swvl1 Raw_ERA_Dataset.nc ../hvgt.nc

#cdo -b 64 chname,swvl1,dvf -addc,1.39E-8 -mulc,0 -timmax -selname,swvl1 Raw_ERA_Dataset.nc ../dvf.nc
#cdo -b 64 chname,swvl1,dvt -addc,2.39E-8 -mulc,0 -timmax -selname,swvl1 Raw_ERA_Dataset.nc ../dvt.nc
