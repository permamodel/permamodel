# Thaw Lake Model-1D

# This model a 1-D numerical model of permafrost and subsidence processes. 
# It aims to investigate the subsurface thermal impact of thaw lakes of various depths, 
# and to evaluate how this impact might change in a warming climate. 

# Key paper: Matell, N., Anderson, R.S., Overeem, I., Wobus, C., Urban, F.,
# Clow, G., in review 2011. Modeling the subsurface thermal impact of Arctic thaw lakes in a warming climate. Computers and Geosciences. 
 
# Copyright (C) <2011> <Nora Matell, Irina Overeem, Robert Anderson, Cameron Wobus>

# Developer can be contacted by irina.overeem@colorado.edu

# Dr. Irina Overeem
# CSDMS Community Surface Dynamics Modeling System
# INSTAAR, University of Colorado at Boulder
# PO Box 450, 80309-0450
# Boulder, CO, USA


# This program is free software; you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation; either version 2 of the License, or (at your option) any later version.
# This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
# You should have received a copy of the GNU General Public License along with this program; if not, write to the Free Software Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA. 

# This model a 1-D numerical model of permafrost and subsidence processes. 
# It aims to investigate the subsurface thermal impact of thaw lakes of various depths, 
# and to evaluate how this impact might change in a warming climate.  

# The model was designed for the Alaskan Arctic Coastal Plain
# The model uses average and amplitude of temperature from observed permafrost temperatures at 5cm in the subsurface at Drew Point data (USGS)
# Clow, G.D., 2008a. Continued permafrost warming in northern Alaska, 2008 update. NOAA/ESRL Global Monitoring Annual Conference, Boulder.  
# LINE 32 CLOW___________ 24 p. error

# The model uses observed heat flow for thermal gradient from:
# Lachenbruch, A.H., Sass, J.H., Marshall, B.V., Moses, T.H., Jr., 1982. Temperatures, heat flow, and the geothermal regime at Prudhoe By, Alaska. Journal of Geophysical Research 87, 9301-9316.  


# in 1-D, subsurface temperature changes through time and space of a semi-infinite
# half-space, then plots result.  The model code incorporates phase change by using
# the apparent heat capacity scheme for temperatures within the
# "phase-change envelope".  

# lake-permafrost model - lake freezes and thaws, permafrost changes temperature,
# when permafrost thaw it subsides, assuming that if excess ice melts all 
# the water leaves the permafrost and thus the subsurface volume decreases.  

# boundary condition of ice bottom and lake water top = Tf
# if Ts<Tf, top boundary condition for calculation is Ts.  if Ts>Tf, top
# boundary condition --> Tf and excess energy is used to melt ice.  When
# ice melted from the top, all other boundary just moved up by equivelant
# amount, if moved up a full control volume then bottom control volume
# added that is same temp as (end-1) control volume.  Assumes that when ice
# free, lake is completely mixed.  (ie Liston & Hall)

import math
import numpy as np
import csv
import operator
import copy
import scipy.io as sio
import matplotlib.pyplot as plt

# set up 1-d grid
dz0 = 0.05;                             # cell size in top 10 m
dz2 = 1;								# cell size below 10 m
deptht = 100;                           # model depth [m]

z_range1 = np.arange(0,10.05,dz0);
z_range2 = np.arange(11,deptht+1,dz2);  
z = np.append(z_range1, z_range2);		# vertical grid (distance below surface) [m]

dz = np.diff(z);
zbetween = z[1:] + np.diff(z)/2;
dzbetween = np.diff(zbetween);
numcells = len(z)                   # number of grid cells
iceinit = 0.0001;                       # initial ice depth [m]
celltype = np.int16(1*np.ones(numcells));   # substrate type - 1=regular soil with excess ice, 2=compacted soil without excess ice, 8=water, 9=ice
celltype[np.where(z<6.7)] =2;			 # depth of subsided cells, this is a function of the initial lake depth (seed lake depth/excess ice content)
celltype[np.where(z<2.0)] = 8;              # seed lake depth
celltype[np.where(z<=iceinit)] = 9;         # ice thickness
celltype[np.where(z>=50)] = 2;              # underlying cells with no excess ice
thawed = np.zeros(numcells);
thawedspec = np.zeros(numcells);

deptht = max(z);
depthsubside = 2.0;			# seed lake depth, has to match with ln 79, because the lake contains water
depthtalik = 0;				# usually unkown, so the model will trend toward equilibrium in the first 100-200 years of a simulation
depthi = iceinit;
print (depthi)
 # set time steps and simulation duration
periodyear = 3600*24*365.;       # period (1 year) [s]
dt = 3600*24;                   # length of timestep [s]
years = 1000;					# simulation duration
nt = int(((periodyear/dt)*years)+1); # number of timesteps
t = np.arange(0,dt*nt+1,dt);			# time at each timestep [s]
tday = [x/(3600*24) for x in t]; 
tyear = [x/(3600*24*365) for x in t];
tyearsonly = np.arange(0,years);
yearnum = 0;

# declare zarrays

depthin = np.zeros(nt);
depthwn = np.zeros(nt);
depthsubsiden = np.zeros(nt+1);
depthtalikn = np.zeros(nt+1);
icethickness = np.zeros(nt+1);
icegrowth = np.zeros(nt+1);
soilthickness = np.zeros(nt+1);
waterthickness = np.zeros(nt+1);
Tbot = np.zeros(nt+1);
T3m = np.zeros(nt+1);
T5m = np.zeros(nt+1);
T10m = np.zeros(nt+1);
T25m = np.zeros(nt+1);
T50m = np.zeros(nt+1);
T100m = np.zeros(nt+1);
Tbotavgann = np.zeros(years);
T3mavgann = np.zeros(years);
T5mavgann = np.zeros(years);
T10mavgann = np.zeros(years);
T25mavgann = np.zeros(years);
T50mavgann = np.zeros(years);
T100mavgann = np.zeros(years);
Trecord = np.zeros((365,numcells));
watts = np.zeros(np.size(t));

# define constants
# Model input and boundary conditions (compiled from field and modeling studies by West & Plug, 2008; Ling and Zhang, 2003;.....
# Ling & Zhang, 2004; Zhou and Huang, 2004; Romanovsky and Osterkamp, 2000; Hinzman, 1998; and French, 2007).

kelvin = 273.15;

rhoi = 917;         # ice density [kg/m3]
ci = 2108;          # ice specific heat capacity [J/kg/K]
ki = 2.18;          # ice thermal conductivity [J/s/m/K]
kappai = ki/(ci*rhoi);
rhow = 1000;        # water density [kg/m3]
cw = 4210;          # water specific heat capacity [J/kg/K]
kw = 0.58;          # water thermal conductivity [J/s/m/K]
kappaw = kw/(cw*rhow);
rhorock = 1200;     # mineral soil density [kg/m3]
crock = 1000;       # mineral soil specific heat capacity [J/kg/degree C]
krock = 1.5;        # rock thermal conductivity [J/s/m/degree C]

totalice = 0.65;    # ice content of original substrate, from bulk samples in top 5m at Drew Point
excessice = 0.3;	# excess ice content - gone once melted
porewater = totalice-excessice;     # porewater ice - still there if refreezes
W = totalice;       # total water content percent of material by mass
Wu = totalice*0.05;  # unfrozen percent water content of material by mass at 
                    #temperature T = Tpc-pcenv
Ws = porewater;
Wus = porewater*0.05;
totalice2 = totalice-Wu;
porewater2 = porewater-Wus;

cf = (ci*totalice2)+(cw*Wu)+(crock*(1-totalice));        
cu = (cw*totalice)+(crock*(1-totalice)); 
cfs = (ci*porewater2)+(cw*Wu)+(crock*(1-porewater));        
cus = (cw*porewater)+(crock*(1-porewater));
rhof = (rhoi*totalice2)+(rhow*Wu)+(rhorock*(1-totalice));
rhou = (rhow*totalice)+(rhorock*(1-totalice));
rhofs = (rhoi*porewater2)+(rhow*Wus)+(rhorock*(1-porewater));
rhous = (rhow*porewater)+(rhorock*(1-porewater));
Cf = cf*rhof;       # frozen volumetric heat capacity [J/m3/degree C]
Cu = cu*rhou;       # thawed volumetric heat capacity [J/m3/degree C]
Cfs = cfs*rhofs;
Cus = cus*rhous;
kf = pow(ki,totalice2)*pow(kw,Wu)*pow(krock,(1-totalice));
ku = pow(kw,totalice)*pow(krock,(1-totalice));
kfs = pow(ki,porewater2)*pow(kw,Wus)*pow(krock,(1-porewater));
kus = pow(kw,porewater)*pow(krock,(1-porewater));

L = 334*1000;       # latent heat of fusion for water [J/kg]
hw = 0.56;          # convective transfer coefficient =[J/s/m2/K]
Kh = 0.0;           # turbulent diffusivity for heat
Tpc = kelvin+0;     # freezing temperature [K]
pcenv = 1;          # width of phase change envelope [m]

# set up surface and initial temperatures, including geothermal gradient
# theoretical warming scenario is set here
Tbar = kelvin-11;    # MAAT [K]
warming = 5;
Tbarplus = Tbar+(float(warming)/36500)*np.arange(1,36501);
Tamp = 17.5;          # amplitude of air temperature fluctuations [K]

Tsurface = [Tbar-Tamp*math.sin(2*math.pi*x/periodyear) for x in t];
Tsurface[182499:218999] = Tbarplus[0:]-Tamp*np.sin(2.*np.pi)*t[182499:218999]/periodyear;
Tsurface[218999:] = Tbar+warming-Tamp*np.sin(2*np.pi)*t[218999:]/periodyear;

q = 0.056;          # mantle heat flow [J/m2/s]
dTdzbase = q/kf;
dTdzbase2 = q/kfs;
Tgrad_first_half = [np.multiply(dTdzbase,z[np.where(z<=10)])]
Tgrad_second_half= [np.multiply(dTdzbase2,z[np.where(z>10)])];
Tgrad = np.append(Tgrad_first_half, Tgrad_second_half);

Twi = Tsurface[0];  # start up water temperature [K]
Tsi = Tbar;         # start up permafrost temperature [K]
Tinit = np.ones((np.size(z)));     # initial temperature grid   
Tinitlake_first_half = [Tpc*Tinit[np.where(celltype == 9)]]
Tinitlake_second_half= [Twi*Tinit[np.where(celltype == 8)]];
Tinitlake = np.append(Tinitlake_first_half, Tinitlake_second_half);  
water = len(Tinitlake);
for n in range(0,numcells):
    if (celltype is 9):
        Tinit[n] = Tpc;
    elif (celltype is 8):
        Tinit[n] = Twi;
    else:
        Tinit[n] = Tsi+Tgrad[n];

# load daily average radiation, function from Drew Point meteorological data
mat_file = np.genfromtxt('radin_dailyavg.csv', delimiter=','); 

#Tinitperm = [Tsi*Tinit(find(celltype==1))+Tgrad(find(celltype==1)),Tsi*Tinit(find(celltype==2))+Tgrad(find(celltype==2))];
#Tinit = [Tinitlake, Tinitperm];

T = Tinit;
zstarthaw = math.sqrt((ku/Cu)*periodyear/(math.pi));
zstarfrozen = math.sqrt((kfs/Cf)*periodyear/(math.pi));
Tbase = Tbar+(dTdzbase2*deptht)+np.multiply(Tamp,np.exp(-deptht/zstarthaw))*np.sin((2*math.pi*np.divide(t,periodyear))-(deptht/zstarthaw));
Tright = -kelvin+Tbar+np.multiply(dTdzbase2,z)+np.multiply(Tamp,np.exp(-z/zstarthaw)); #outer edges of the funnel
Trightf = -kelvin+Tbar+np.multiply(dTdzbase2,z)+np.multiply(Tamp,np.exp(-z/zstarfrozen)); #outer edges of the funnel
Tleftf = -kelvin+Tbar+np.multiply(dTdzbase2,z)-np.multiply(Tamp,np.exp(-z/zstarfrozen));
Tleft = -kelvin+Tbar+np.multiply(dTdzbase2,z)-np.multiply(Tamp,np.exp(-z/zstarthaw));
#Ts0 = Tsurface[1];

count = 1;
daynum = 0;
numplots = 50;
countprint = dt*nt/numplots;
nplot=0;
fig_one = plt.figure(1)
zero = np.zeros((np.size(z)));
plt.plot(zero,-z, linewidth = 2);
plt.plot(Tright,-z,Tleft,-z,Trightf,-z,Tleftf,-z,)
plt.plot(Tinit-kelvin,-z)
plt.axis([Tbar-Tamp-kelvin, Tbar+Tamp-kelvin, -deptht, 0]);
plt.xlabel('Temperature (C)',fontsize = 18)
plt.ylabel('Depth (m)',fontsize = 18)
#fig_one = plt.show();


# model equations
print ('starting now')
loopcount  = 0;
print (nt);
for n in range(0,nt):
	loopcount = loopcount+1;
	print ('____LOOPCOUNT ' , loopcount); 
	thawedspec = np.zeros(numcells);
	thawed = np.zeros(numcells);	
	count = count+1;
	time = (n+1)*dt;    # time into run [s]
	error = 1;
	Told = T;       # remember temperatures from last timestep... 
	Ts0 = Tsurface[n];
	T[0]= Ts0;
	
	if depthi>0 and Tsurface[n]>Tpc:
		Ts0 = Tpc;
		
	water = -1;
	ice = -1
	subsided = -1;
	regsoil = -1;
	
	for i in range(0,numcells):		
		if celltype[i] == 8:
			water = water+1; 
		elif (celltype[i] == 9):
			ice = ice+1;
			#print 'updated ice is ', ice;
		elif (celltype[i] == 2):
			subsided = subsided+1;               
		elif (celltype[i] == 1):
			regsoil = regsoil+1;
	
	print ('ice is =', ice, ' ', 'Water is = ', water);

	day = (n % 365);
	solarrad = mat_file[day];
	#print "water on top _ ", water;
	#compute temperatures in water layer if no ice exists
	#water = 0;
	#ice = 38;
	if celltype[0] == 8:
		print ('in if loop_________')
		Twater = T[0:water];  # T(find(celltype==8));
		Tmix = np.mean(Twater); #mix again
		sextinc = 0.6;       # solar extinction constant
		albedo = 0.06;
		qrad = (1-albedo)*solarrad*np.exp(-sextinc*zbetween[0:water]);
		watts[count] = solarrad;
		error = 1;
		count = 0
		while error>0.0001:
			Tcalc = [0]
			Tit = Twater;
			for i in range (1,water-1):
				coeff1 = kw/dz[i];
				coeff2 = kw/dz[i-1];
				coeff3 = cw*rhow*dzbetween[i-1]/dt;
				coeff4 = coeff1+coeff2+coeff3;
				Tcalc = np.append(Tcalc, (coeff1*T[i+1]+coeff2*T[i-1]+coeff3*Told[i]+qrad[i-1]-qrad[i])/coeff4);
			Twater_Partial = np.append(Ts0,Tcalc[1:]);
			Twater = np.append(Twater_Partial, Tcalc[-1]);
			error = max(abs(Twater-Tit));	
		Twater2 = Twater[0:-1];
		Tmix = np.mean(Twater2);  # thoroughly mix lake
		T[1:water] = [Tmix*np.ones(np.size(x))for x in Twater2];
	# compute temperatures in ice layer if completely frozen to bottom
	
	elif water == 0:
		print ('___________________in elif loop___________________')
		Tice = (T[np.where(celltype==9)]);
		#print 'Starting in the loop ', Tice
		while error>0.0001:
			Tcalc = [0]
			#Tcalc= np.zeros(ice);
			Tit = copy.deepcopy(Tice);
			for i in range(1,ice):
				coeff1 = ki/dz[i];
				coeff2 = ki/dz[i-1];
				coeff3 = ci*rhoi*dzbetween[i-1]/dt;
				coeff4 = coeff1+coeff2+coeff3;
				Tcalc = np.append(Tcalc, (coeff1*T[i+1]+coeff2*T[i-1]+coeff3*Told[i])/coeff4);			
			#print ' Tcalc is ', Tcalc;
			Tice = np.append(Ts0,Tcalc[0:]);
			#print 'Ts0 is ', Ts0;
			#print "The final___ TICE ", Tice;
			#print 'Tice length is ', len(Tice); #(Tice is 38)
			#print 'Tit length is ', len(Tit);	#(Tit is 39)		
			error = max(abs(Tice-Tit));
			#print error;
           
		T[np.where(celltype==9)] = Tice; 
	#compute temperatures in ice and water layers if both exist
	
	else:
		print ('__in else loop')
		print (ice);
		if Tsurface[n]>=Tpc:
			T[0:ice] = Tpc;
		elif ice<3:
			T[0] = Tsurface[n];
			T[ice] = Tpc; 
			print ("whatss elif");
		else:
			Ti = T[0:ice];
			error = 1;
			while error>0.0001:
				Tcalci = [0]
				#Tcalci= np.zeros(ice);
				Titi = copy.deepcopy(Ti);
				for i in range(1,ice-1):
					coeff1 = ki/dz[i];
					coeff2 = ki/dz[i-1];
					coeff3 = ci*rhoi*dzbetween[i-1]/dt;
					coeff4 = coeff1+coeff2+coeff3;
					Tcalci = np.append(Tcalci, (coeff1*T[i+1]+coeff2*T[i-1]+coeff3*Told[i])/coeff4);
				Ti_partial=np.append(Ts0,Tcalci[1:]);
				Ti = np.append(Ti_partial,Tpc);
				error = max(abs(Ti-Titi));
			T[0:ice] = Ti;
						
		error = 1;
		if water==1:
			T[np.where(celltype==8)] = Tpc;
		else:
			print ("this loop maybe");
			error = 1;
			while error>0.0001:
				#Tcalcw = [0]
				Tit = copy.deepcopy(T);
				Tcalcw= np.zeros(ice);
				for i in range (ice,ice+water+1):
					coeff1 = kw/dz[i];
					coeff2 = kw/dz[i-1];
					coeff3 = cw*rhow*dzbetween[i-1]/dt;
					coeff4 = coeff1+coeff2+coeff3;
					Tcalcw = np.append(Tcalcw, (coeff1*T[i+1]+coeff2*T[i-1]+coeff3*Told[i])/coeff4);
				Tw = np.append(Tpc, Tcalcw[ice:]); #problem appears to be in Tcalcw[ice:]
				T[ice:ice+water+2] = Tw;
				error = max(abs(T-Tit));
			
    # Calculate temperatures in underlying permafrost
	#clear 'z2' 'zbetween2' 'Told2' 'T2' 'celltype2'
	z2 = z[ice+water-1:];       
	zbetween2 = zbetween[ice+water-1:];
	Told2 = Told[ice+water-1:];
	T2 = T[ice+water-1:];
	numcells2 = len(T2);
	celltype2 = celltype[ice+water-1:];    
	C = np.zeros(numcells2);
	k = np.zeros(numcells2);
	Tcalc = 0;
	for i in range (0,numcells2):
		if celltype2[i] == 8:
			C[i] = cw*rhow;
			k[i] = kw;
		elif celltype2[i] == 9:
			C[i] = ci*rhoi;
			k[i] = ki;
		elif celltype2[i]==1 and T2[i] < (Tpc-pcenv):
			C[i] = Cf;
			k[i] = kf;
			#print('real permafrost, frozen')
		elif celltype2[i]==1 and T2[i]>=(Tpc-pcenv) and T2[i]<=Tpc:
			C[i] = Cf+L*rhof*((W-Wu)/pcenv);
			k[i] = kf+((ku-kf)/pcenv)*(T2(i)-(Tpc-pcenv));
			print('real permafrost, thawing')
		elif celltype2[i]==1 and T2[i]>Tpc:
			C[i] = Cu;
			k[i] = ku;
			print('real permafrost, thawed')
		elif celltype2[i]==2 and T2[i]<(Tpc-pcenv):
			C[i] = Cfs;
			k[i] = kfs;
		elif celltype2[i]==2 and T2[i]>=(Tpc-pcenv) and T2[i]<=Tpc:
			C[i] = Cfs+L*rhofs*((Ws-Wus)/pcenv);
			k[i] = kfs+((kus-kfs)/pcenv)*(T2[i]-(Tpc-pcenv));
		elif celltype2[i]==2 and T2[i]>Tpc:
			C[i] = Cus;
			k[i] = kus;
	
	kbetween = k[1:]-(np.diff(k)/2);
	dz2 = np.diff(z2);
	dzbetween2 = np.diff(zbetween2);
	error = 1;
	#Tcalcperm= np.zeros(251)
	while error>0.0001:
		Tcalcperm = [0]
		Tit = copy.deepcopy(T2);
		for i in range(1,numcells2-1):
			coeff1 = kbetween[i]/dz2[i];
			coeff2 = kbetween[i-1]/dz2[i-1];
			coeff3 = C[i-1]*dzbetween2[i-1]/dt;
			coeff4 = coeff1+coeff2+coeff3;
			Tcalcperm = np.append(Tcalcperm, ((coeff1*T2[i+1]+coeff2*T2[i-1]+coeff3*Told2[i])/coeff4));
			#Tcalcperm[i] = (coeff1*T2[i+1]+coeff2*T2[i-1]+coeff3*Told2[i])/coeff4;
		T2_partial = np.append(T[ice+water],Tcalcperm[1:])
		T2 = np.append(Tcalcperm,Tbase[n]);
		error = max(abs(T2-Tit));
	T[ice+water-1:] = T2;
	
	# subside permafrost
	print (len(range(water+ice,numcells)));
	for i in range (water+ice,numcells):
		if T[i]<Tpc:
			#print 'before break;'
			break 
		elif T[i]>Tpc and celltype[i] == 1:
			thawedspec[i] = 1;
			thawed[i] = 1;
			celltype[i] = 2;
			#print 'is it here'
		elif T[i]>Tpc and celltype[i] == 2:
			thawed[i] = 1;
			#print 'where is it'
				
##[Tbar-Tamp*math.sin(2*math.pi*x/periodyear) for x in t];
	numthawedspec = np.sum(thawedspec)+1;
	depthsubsidenew = z[numthawedspec]*excessice;
	depthsubside = depthsubside+depthsubsidenew;
	depthsubsiden[n] = depthsubside;
	cellsubside = np.int16(math.floor((depthsubside)/dz0));
	maxcellsubside = cellsubside; 
	#if depthsubside>0:	#depthsubsidenew>0
	#	print celltype;
	#	for i in range(ice,maxcellsubside):
	#		celltype[i] = 8;
			
	if (Tsurface[n]<Tpc):		
		icecheckdepth = np.int16(0.33/dz0)+ice;    #z location of ~0.33 m into water column
		if icecheckdepth>ice+water:
			icecheckdepth = ice+water-1;
		Twater = T[icecheckdepth];
		if water==0:
			Twater = T[np.where(celltype==8)];
			Twater = Twater[0];
			print ("here");
		if Twater<Tpc and ice==-1:
			depthi = 0.01;
			
		print
		
		#_____DOUBLE CHECK HERE_______
		depthmix = z[icecheckdepth]-depthi;
		if depthi> 0:
			ddepthidt_partial = ((Tpc-Ts0)*pow((depthi/ki),(-1))-hw*(Twater-Tpc));
			ddepthidt = np.divide(ddepthidt_partial, rhoi*L)
			ddepthi = ddepthidt*dt;
		#### ADJUSTMENT HERE
		if ddepthi>0.05:
			ddepthi = 1.05;		
			print('__________________had to adjust________________');
		icegrowth[n] = ddepthi;	
		depthi = depthi+ddepthi;
		
		if depthi<0.00001:
			depthi = 0;
			
		depthw = depthsubside-depthi;		
		if depthi>depthsubside:
			depthi = depthsubside;
			depthw = 0;
			
	#     check if Ts>Tpc (if there is ice there already)
	elif Tsurface[n]>Tpc and depthi>0:
		Qm = (Tsurface[n]-Tpc)*(pow((dz0/ki),(-1)));
		dMdt = Qm/(L*rhoi);
		dM = dMdt*dt;
		icegrowth[n] = -dM;
		depthi = depthi-dM;
		if depthi<0:
			depthi = 0;
	
	maxcelli = np.int16(round(depthi/dz0));
	celltype[0:ice+water] = 8;

	celltype[0:maxcelli] = 9;
	icethickness[n] = depthi;
	soilthickness[n] = deptht-depthsubside;
	waterthickness[n] = deptht-icethickness[n]-soilthickness[n];
	Tbot[n] = T[ice+water]-kelvin;
	T3m[n] = T[59]-kelvin;
	T5m[n] = T[99]-kelvin;
	T10m[n] = T[199]-kelvin;
	T25m[n] = T[214]-kelvin;
	T50m[n] = T[239]-kelvin;
	T100m[n] = T[289]-kelvin;
	if (tday[n]%365)==0 and tday[n]>364:
		yearnum = yearnum+1;
		Tbotavgann[yearnum] = np.mean(Tbot[n-364:n]);
		T3mavgann[yearnum] = np.mean(T3m[n-364:n]);
		T5mavgann[yearnum] = np.mean(T5m[n-364:n]);
		T10mavgann[yearnum] = np.mean(T10m[n-364:n]);
		T25mavgann[yearnum] = np.mean(T25m[n-364:n]);
		T50mavgann[yearnum] = np.mean(T50m[n-364:n]);
		T100mavgann[yearnum] = np.mean(T100m[n-364:n]);

	if n>=nt-364:
		daynum = daynum+1;
		Trecord[daynum, 0:] = T;
			
	if((n % (365*50))==0):
		#    nplot = nplot+1
		fig_one = plt.figure(1)
		#    ice
		icethickness[n]
		#    water
		waterthickness[n]

		lakedepthx = np.arange(-20,25,1);
		lakedepth = icethickness[n]+waterthickness[n]*np.ones(np.size(lakedepthx));
		Tplot = T-kelvin;
		plt.plot(Tplot,-z,'m')
		timeplot = tday[n]/365 #will stamp time in years on your screen
		plt.xlabel('temperature',fontsize = 18)
		plt.ylabel('depth', fontsize = 18)
		#
		#figure(2)
		#clf(2)
		icedepth = icethickness[n]*np.ones(np.size(lakedepthx));
		plt.plot(zero,-z, linewidth = 2);
		plt.plot(lakedepthx,-lakedepth, lakedepthx,-icedepth,Tplot,-z)
		plt.xlabel('temperature', fontsize = 18)
		plt.ylabel('depth', fontsize = 18)
		#axis([-30 25 -5 0])

# Plot output parameters after run is completed
print ('after the loop')
'''# plot surface temperature fluctuations
fig_two = plt.figure(2)
TsK = 100;#Tsurface-kelvin;
plt.plot(tday,TsK)
plt.xlabel('time(days)',fontsize = 18)
plt.ylabel('surface temperature(degrees C)',fontsize = 18)    

# plot thickness
fig_three = figure(3)
plt.plot(tyear[1:end-1],icethickness[1:end-1],tyear[1:end-1],waterthickness[1:end-1],tyear[1:end-1],depthsubsiden[1:end-1])
plt.xlabel('time(years)','fontname','arial','fontsize',18)
plt.ylabel('thickness (m)','fontname','arial','fontsize',18)
plt.legend('ice thickness','water thickness','depth subside')

# plot temps at various depths
fig_four = figure(4)
T0C = np.zeros(size(tyear));
plt.plot(tyearsonly,Tbotavgann,tyearsonly,T3mavgann,tyearsonly,T5mavgann,tyearsonly,T10mavgann,tyearsonly,T25mavgann, tyearsonly,T50mavgann,tyearsonly,T100mavgann,tyear[1:end-1],T0C[1:end-1])
plt.xlabel('time(years)','fontname','arial','fontsize',18)
plt.ylabel('temperature (C)','fontname','arial','fontsize',18)
plt.legend('at lake bottom','at 3m depth','at 5m depth','at 10m depth','at 25m depth','at 50m depth','at 100m depth','0C')

# plot temps at lake bottom
fig_five = figure(5)
T0C = np.zeros(size(tyear));
plt.plot(tyear[1:end-1],Tbot[1:end-1],tyearsonly,Tbotavgann,tyear[1:end-1],T0C[1:end-1])
plt.xlabel('time(years)','fontname','arial','fontsize',18)
plt.ylabel('temperature (C)','fontname','arial','fontsize',18)
plt.legend('at lake bottom','lake bottom average annual temp','0C')

fig_two = plt.show();
fig_three = plt.show();
fig_four = plt.show();
fig_five = plt.show();

input()
'''
