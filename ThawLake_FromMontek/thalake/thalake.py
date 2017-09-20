# Thaw Lake Model-1D

# This model a 1-D numerical model of permafrost and subsidence processes, it is a Python version. 
# It aims to investigate the subsurface thermal impact of thaw lakes of various depths, 
# and to evaluate how this impact might change in a warming climate. 

# Key paper: Matell, N., Anderson, R.S., Overeem, I., Wobus, C., Urban, F.,
# Clow, G., 2013 Modeling the subsurface thermal impact of Arctic thaw lakes in a warming climate. Computers and Geosciences. 
 
# Originates from earlier code by Nora Matell and co-authors.
# Copyright to Python version  (C) <2017> <Irina Overeem, Montek Singh>

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

# The model uses observed heat flow for thermal gradient from:
# Lachenbruch, A.H., Sass, J.H., Marshall, B.V., Moses, T.H., Jr., 1982. Temperatures, heat flow, and the geothermal regime at Prudhoe Bay, Alaska. Journal of Geophysical Research 87, 9301-9316.  


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
#from permamodel.utils import model_input
#from permamodel.components import perma_base

class thaLakeModel():
    print ("hello maybe this is working")
    def __init__(self):
        # set up 1-d grid
        self.dz0 = 0.05;                             # cell size in top 10 m
        self.dz2 = 1;                                # cell size below 10 m
        self.deptht = 100;                           # model depth [m]
        self.z_range1 = np.arange(0,10.05,self.dz0);
        self.z_range2 = np.arange(11,self.deptht+1,self.dz2); 
        self.iceinit = 0.0001;                       # initial ice depth [m]
        self.periodyear = 3600*24*365.;       # period (1 year) [s]
        self.dt = 3600*24;                   # length of timestep [s]
        self.years = 1000;                    # simulation duration
        self.depthsubside = 2.0;            # seed lake depth, has to match with ln 79, because the lake contains water
        self.depthtalik = 0;                # usually unkown, so the model will trend toward equilibrium in the first 100-200 years of a simulation
        self.yearnum = 0;
        self.kelvin = 273.15;
        self.rhoi = 917;         # ice density [kg/m3]
        self.ci = 2108;          # ice specific heat capacity [J/kg/K]
        self.ki = 2.18;          # ice thermal conductivity [J/s/m/K]
        self.rhow = 1000;        # water density [kg/m3]
        self.cw = 4210;          # water specific heat capacity [J/kg/K]
        self.kw = 0.58;          # water thermal conductivity [J/s/m/K]
        self.rhorock = 1200;     # mineral soil density [kg/m3]
        self.crock = 1000;       # mineral soil specific heat capacity [J/kg/degree C]
        self.krock = 1.5;        # rock thermal conductivity [J/s/m/degree C]
        self.totalice = 0.65;    # ice content of original substrate, from bulk samples in top 5m at Drew Point
        self.excessice = 0.3;    # excess ice content - gone once melted    
        self.L = 334*1000;       # latent heat of fusion for water [J/kg]
        self.hw = 0.56;          # convective transfer coefficient =[J/s/m2/K]
        self.Kh = 0.0;           # turbulent diffusivity for heat
        self.pcenv = 1;          # width of phase change envelope [m]
        self.warming = 5;
        self.Tamp = 17.5;          # amplitude of air temperature fluctuations [K]
        self.q = 0.056;          # mantle heat flow [J/m2/s]
        self.count = 1;
        self.daynum = 0;
        self.numplots = 50;    
    
    def initialize(self): 
        self.z = np.append(self.z_range1, self.z_range2);        # vertical grid (distance below surface) [m]
        self.dz = np.diff(self.z);
        self.zbetween = self.z[1:] + np.diff(self.z)/2;
        self.dzbetween = np.diff(self.zbetween);
        self.numcells = len(self.z)                   # number of grid cells
        self.celltype = np.int16(1*np.ones(self.numcells));   # substrate type - 1=regular soil with excess ice, 2=compacted soil without excess ice, 8=water, 9=ice
        self.celltype[np.where(self.z<6.7)] =2;             # depth of subsided cells, this is a function of the initial lake depth (seed lake depth/excess ice content)
        self.celltype[np.where(self.z<2.0)] = 8;              # seed lake depth
        self.celltype[np.where(self.z<=self.iceinit)] = 9;         # ice thickness
        self.celltype[np.where(self.z>=50)] = 2;              # underlying cells with no excess ice
        self.thawed = np.zeros(self.numcells);
        self.thawedspec = np.zeros(self.numcells);

        self.deptht = max(self.z);
        self.depthi = self.iceinit;

         # set time steps and simulation duration

        self.nt = int(((self.periodyear/self.dt)*self.years)+1); # number of timesteps
        self.t = np.arange(0,self.dt*self.nt+1,self.dt);            # time at each timestep [s]
        self.tday = [x/(3600*24) for x in self.t]; 
        self.tyear = [x/(3600*24*365) for x in self.t];
        self.tyearsonly = np.arange(0,self.years);

        # declare zarrays

        self.depthin = np.zeros(self.nt);
        self.depthwn = np.zeros(self.nt);
        self.depthsubsiden = np.zeros(self.nt+1);
        self.depthtalikn = np.zeros(self.nt+1);
        self.icethickness = np.zeros(self.nt+1);
        self.icegrowth = np.zeros(self.nt+1);
        self.soilthickness = np.zeros(self.nt+1);
        self.waterthickness = np.zeros(self.nt+1);
        self.Tbot = np.zeros(self.nt+1);
        self.T3m = np.zeros(self.nt+1);
        self.T5m = np.zeros(self.nt+1);
        self.T10m = np.zeros(self.nt+1);
        self.T25m = np.zeros(self.nt+1);
        self.T50m = np.zeros(self.nt+1);
        self.T100m = np.zeros(self.nt+1);
        self.Tbotavgann = np.zeros(self.years);
        self.T3mavgann = np.zeros(self.years);
        self.T5mavgann = np.zeros(self.years);
        self.T10mavgann = np.zeros(self.years);
        self.T25mavgann = np.zeros(self.years);
        self.T50mavgann = np.zeros(self.years);
        self.T100mavgann = np.zeros(self.years);
        self.Trecord = np.zeros((365,self.numcells));
        self.watts = np.zeros(np.size(self.t));


        # define constants
        # Model input and boundary conditions (compiled from field and modeling studies by West & Plug, 2008; Ling and Zhang, 2003;.....
        # Ling & Zhang, 2004; Zhou and Huang, 2004; Romanovsky and Osterkamp, 2000; Hinzman, 1998; and French, 2007).

        self.kappai = self.ki/(self.ci*self.rhoi);
        self.kappaw = self.kw/(self.cw*self.rhow);

        self.porewater = self.totalice-self.excessice;     # porewater ice - still there if refreezes
        self.W = self.totalice;       # total water content percent of material by mass
        self.Wu = self.totalice*0.05;  # unfrozen percent water content of material by mass at 
                            #temperature T = Tpc-pcenv
        self.Ws = self.porewater;
        self.Wus = self.porewater*0.05;
        self.totalice2 = self.totalice-self.Wu;
        self.porewater2 = self.porewater-self.Wus;

        self.cf = (self.ci*self.totalice2)+(self.cw*self.Wu)+(self.crock*(1-self.totalice));        
        self.cu = (self.cw*self.totalice)+(self.crock*(1-self.totalice)); 
        self.cfs = (self.ci*self.porewater2)+(self.cw*self.Wu)+(self.crock*(1-self.porewater));        
        self.cus = (self.cw*self.porewater)+(self.crock*(1-self.porewater));
        self.rhof = (self.rhoi*self.totalice2)+(self.rhow*self.Wu)+(self.rhorock*(1-self.totalice));
        self.rhou = (self.rhow*self.totalice)+(self.rhorock*(1-self.totalice));
        self.rhofs = (self.rhoi*self.porewater2)+(self.rhow*self.Wus)+(self.rhorock*(1-self.porewater));
        self.rhous = (self.rhow*self.porewater)+(self.rhorock*(1-self.porewater));
        self.Cf = self.cf*self.rhof;       # frozen volumetric heat capacity [J/m3/degree C]
        self.Cu =  self.cu*self.rhou;       # thawed volumetric heat capacity [J/m3/degree C]
        self.Cfs = self.cfs*self.rhofs;
        self.Cus = self.cus*self.rhous;
        self.kf = pow(self.ki,self.totalice2)*pow(self.kw,self.Wu)*pow(self.krock,(1-self.totalice));
        self.ku = pow(self.kw,self.totalice)*pow(self.krock,(1-self.totalice));
        self.kfs = pow(self.ki,self.porewater2)*pow(self.kw,self.Wus)*pow(self.krock,(1-self.porewater));
        self.kus = pow(self.kw,self.porewater)*pow(self.krock,(1-self.porewater));

        self.Tpc = self.kelvin+0;     # freezing temperature [K]

        # set up surface and initial temperatures, including geothermal gradient
        # theoretical warming scenario is set here
        self.Tbar = self.kelvin-11;    # MAAT [K]

        self.Tbarplus = self.Tbar+(float(self.warming)/36500)*np.arange(1,36501);


        self.Tsurface = [self.Tbar-self.Tamp*math.sin(2*math.pi*x/self.periodyear) for x in self.t];
        self.Tsurface[182499:218999] = self.Tbarplus[0:]-self.Tamp*np.sin(2.*np.pi)*self.t[182499:218999]/self.periodyear;
        self.Tsurface[218999:] = self.Tbar+self.warming-self.Tamp*np.sin(2*np.pi)*self.t[218999:]/self.periodyear;

        self.dTdzbase = self.q/self.kf;
        self.dTdzbase2 = self.q/self.kfs;
        self.Tgrad_first_half = [np.multiply(self.dTdzbase,self.z[np.where(self.z<=10)])]
        self.Tgrad_second_half= [np.multiply(self.dTdzbase2,self.z[np.where(self.z>10)])];
        self.Tgrad = np.append(self.Tgrad_first_half, self.Tgrad_second_half);

        self.Twi = self.Tsurface[0];  # start up water temperature [K]
        self.Tsi = self.Tbar;         # start up permafrost temperature [K]
        self.Tinit = np.ones((np.size(self.z)));     # initial temperature grid   
        self.Tinitlake_first_half = [self.Tpc*self.Tinit[np.where(self.celltype == 9)]]
        self.Tinitlake_second_half= [self.Twi*self.Tinit[np.where(self.celltype == 8)]];
        self.Tinitlake = np.append(self.Tinitlake_first_half, self.Tinitlake_second_half);  
        self.water = len(self.Tinitlake);
        for n in range(0,self.numcells):
            if (self.celltype is 9):
                self.Tinit[n] = self.Tpc;
            elif (self.celltype is 8):
                self.Tinit[n] = self.Twi;
            else:
                self.Tinit[n] = self.Tsi+self.Tgrad[n];

        # load daily average radiation, function from Drew Point meteorological data
        self.mat_file = np.genfromtxt('radin_dailyavg.csv', delimiter=','); 

        #Tinitperm = [Tsi*Tinit(find(celltype==1))+Tgrad(find(celltype==1)),Tsi*Tinit(find(celltype==2))+Tgrad(find(celltype==2))];
        #Tinit = [Tinitlake, Tinitperm];

        self.T = self.Tinit;
        self.zstarthaw = math.sqrt((self.ku/self.Cu)*self.periodyear/(math.pi));
        self.zstarfrozen = math.sqrt((self.kfs/self.Cf)*self.periodyear/(math.pi));
        self.Tbase = self.Tbar+(self.dTdzbase2*self.deptht)+np.multiply(self.Tamp,np.exp(-self.deptht/self.zstarthaw))*np.sin((2*math.pi*np.divide(self.t,self.periodyear))-(self.deptht/self.zstarthaw));
        self.Tright = -self.kelvin+self.Tbar+np.multiply(self.dTdzbase2,self.z)+np.multiply(self.Tamp,np.exp(-self.z/self.zstarthaw)); #outer edges of the funnel
        self.Trightf = -self.kelvin+self.Tbar+np.multiply(self.dTdzbase2,self.z)+np.multiply(self.Tamp,np.exp(-self.z/self.zstarfrozen)); #outer edges of the funnel
        self.Tleftf = -self.kelvin+self.Tbar+np.multiply(self.dTdzbase2,self.z)-np.multiply(self.Tamp,np.exp(-self.z/self.zstarfrozen));
        self.Tleft = -self.kelvin+self.Tbar+np.multiply(self.dTdzbase2,self.z)-np.multiply(self.Tamp,np.exp(-self.z/self.zstarthaw));
        #Ts0 = Tsurface[1];


        self.countprint = self.dt*self.nt/self.numplots;
        self.nplot=0;
        self.Loopcount = 0;
        
    def updateModel(self):
        for n in range(0,5):
            self.Loopcount = self.Loopcount+1; # this is a total loop count for debugging a.o.
            print (self.Loopcount);
            self.thawedspec = np.zeros(self.numcells);
            self.thawed = np.zeros(self.numcells);    
            self.count = self.count+1;
            time = (n+1)*self.dt;    # time into run [s]
            error = 1;
            Told = self.T;       # remember temperatures from last timestep... 
            self.Ts0 = self.Tsurface[n];
            self.T[0]= self.Ts0;
                
            if self.depthi>0 and self.Tsurface[n]>self.Tpc:
                self.Ts0 = self.Tpc;
                    
            self.water = -1;
            self.ice = -1;
            self.subsided = -1;
            self.regsoil = -1;

            for i in range(0,self.numcells):        
                if self.celltype[i] == 8:
                    self.water = self.water+1; 
                elif (self.celltype[i] == 9):
                    self.ice = self.ice+1;
                elif (self.celltype[i] == 2):
                    self.subsided = self.subsided+1;               
                elif (self.celltype[i] == 1):
                    self.regsoil = self.regsoil+1;
                
            self.day = (n % 365);
            self.solarrad = self.mat_file[self.day];                
                
            if self.celltype[0] == 8:
                print ('in if loop')
                Twater = T[0:water];  # T(find(celltype==8));
                Tmix = np.mean(Twater); #mix again
                sextinc = 0.6;       # solar extinction constant
                albedo = 0.06;
                qrad = (1-albedo)*solarrad*np.exp(-sextinc*zbetween[0:water]);
                watts[count] = solarrad;
                error = 1;
                count = 0
                while error>0.0001:
                    self.Tcalc = [0]
                    self.Tit = Twater;
                    for i in range (1,self.water-1):
                        self.coeff1 = self.kw/self.dz[i];
                        self.coeff2 = self.kw/self.dz[i-1];
                        self.coeff3 = self.cw*self.rhow*self.dzbetween[i-1]/self.dt;
                        self.coeff4 = self.coeff1+self.coeff2+self.coeff3;
                        self.Tcalc = np.append(self.Tcalc, ((self.coeff1*self.T[i+1]+self.coeff2*T[i-1]+coeff3*Told[i]+qrad[i-1]-qrad[i])/coeff4));
                    self.Twater_Partial = np.append(self.Ts0,self.Tcalc[1:]);
                    self.Twater = np.append(self.Twater_Partial, self.Tcalc[-1]);
                    error = max(abs(self.Twater-self.Tit));    
                self.Twater2 = self.Twater[0:-1];
                self.Tmix = np.mean(self.Twater2);  # thoroughly mix lake
                self.T[1:self.water] = [self.Tmix*np.ones(np.size(x))for x in self.Twater2];
            # compute temperatures in ice layer if completely frozen to bottom            
                
            elif self.water == 0:
                print ('in elif loop')
                self.Tice = (self.T[np.where(self.celltype==9)]);
                #print 'Starting in the loop ', Tice
                while error>0.0001:
                    self.Tcalc = [0]
                    Tit = copy.deepcopy(Tice);
                    for i in range(1,self.ice):
                        self.coeff1 = self.ki/self.dz[i];
                        self.coeff2 = self.ki/self.dz[i-1];
                        self.coeff3 = self.ci*self.rhoi*self.dzbetween[i-1]/self.dt;
                        self.coeff4 = self.coeff1+self.coeff2+self.coeff3;
                        self.Tcalc = np.append(self.Tcalc, (self.coeff1*self.T[i+1]+self.coeff2*self.T[i-1]+self.coeff3*self.Told[i])/self.coeff4);            
                    self.Tice = np.append(self.Ts0,self.Tcalc[1:]);
                    #print "The final___ TICE ", Tice;            
                    error = max(abs(self.Tice-self.Tit));
                   
                self.T[np.where(self.celltype==9)] = self.Tice; 
            #compute temperatures in ice and water layers if both exist    
            
            else:
                print ("else loop")
                if self.Tsurface[n]>=self.Tpc:
                    self.T[0:self.ice] = self.Tpc;
                    print ("Else -If")
                elif self.ice<3:
                    self.T[0] = self.Tsurface[n];
                    self.T[self.ice-1] = self.Tpc; 
                    print ("Maybe Elif it is") 
                else:
                    print ('Maybe in this loop')
                    self.Ti = self.T[0:self.ice];
                    error = 1;
                    while error>0.0001:
                        self.Tcalci = [0]
                        #Tcalci= np.zeros(ice);
                        self.Titi = copy.deepcopy(Ti);
                        for i in range(1,self.ice-1):
                            coeff1 = self.ki/dz[i];
                            coeff2 = self.ki/self.dz[i-1];
                            self.coeff3 = self.ci*self.rhoi*self.dzbetween[i-1]/self.dt;
                            self.coeff4 = self.coeff1+self.coeff2+self.coeff3;
                            self.Tcalci = np.append(self.Tcalci, (self.coeff1*self.T[i+1]+self.coeff2*self.T[i-1]+self.coeff3*self.Told[i])/self.coeff4);
                        self.Ti_partial=np.append(self.Ts0,self.Tcalci[1:]);
                        self.Ti = np.append(self.Ti_partial,self.Tpc);
                        error = max(abs(self.Ti-self.Titi));
                    self.T[0:self.ice] = self.Ti;        
                error = 1;
                if self.water==1:
                    self.T[np.where(self.celltype==8)] = self.Tpc;
                else:
                    error = 1;
                    while error>0.0001:
                        Tcalcw= np.zeros(self.ice);
                        self.Tit = copy.deepcopy(self.T);
                        for i in range (self.ice,self.ice+self.water+1):
                            coeff1 = self.kw/self.dz[i];
                            coeff2 = self.kw/self.dz[i-1];
                            coeff3 = self.cw*self.rhow*self.dzbetween[i-1]/self.dt;
                            coeff4 = coeff1+coeff2+coeff3;
                            Tcalcw = np.append(Tcalcw, (coeff1*self.T[i+1]+coeff2*self.T[i-1]+coeff3*Told[i])/coeff4);
                        Tw = np.append(self.Tpc, Tcalcw[self.ice:]);
                        self.T[self.ice:self.ice+self.water+2] = Tw;
                        error = max(abs(self.T-self.Tit));
                        
            # Calculate temperatures in underlying permafrost
            #clear 'z2' 'zbetween2' 'Told2' 'T2' 'celltype2'
            #clear 'z2' 'zbetween2' 'Told2' 'T2' 'celltype2'
            z2 = self.z[self.ice+self.water-1:];       
            zbetween2 = self.zbetween[self.ice+self.water-1:];
            Told2 = Told[self.ice+self.water-1:];
            T2 = self.T[self.ice+self.water-1:];
            numcells2 = len(T2);
            celltype2 = self.celltype[self.ice+self.water-1:];    
            C = np.zeros(numcells2);
            k = np.zeros(numcells2);
            Tcalc = 0;
            for i in range (0,numcells2):
                if celltype2[i] == 8:
                    C[i] = self.cw*self.rhow;
                    k[i] = self.kw;
                elif celltype2[i] == 9:
                    C[i] = self.ci*self.rhoi;
                    k[i] = self.ki;
                elif celltype2[i]==1 and T2[i] < (self.Tpc-self.pcenv):
                    C[i] = self.Cf;
                    k[i] = self.kf;
                    #print('real permafrost, frozen')
                elif celltype2[i]==1 and T2[i]>=(self.Tpc-self.pcenv) and T2[i]<=self.Tpc:
                    C[i] = self.Cf+self.L*self.rhof*((self.W-self.Wu)/self.pcenv);
                    k[i] = self.kf+((self.ku-self.kf)/self.pcenv)*(T2[i]-(self.Tpc-self.pcenv));
                    print('real permafrost, thawing')
                elif celltype2[i]==1 and T2[i]>self.Tpc:
                    C[i] = self.Cu;
                    k[i] = self.ku;
                    print('real permafrost, thawed')
                elif celltype2[i]==2 and T2[i]<(self.Tpc-self.pcenv):
                    C[i] = self.Cfs;
                    k[i] = self.kfs;
                elif celltype2[i]==2 and T2[i]>=(self.Tpc-self.pcenv) and T2[i]<=self.Tpc:
                    C[i] = self.Cfs+self.L*self.rhofs*((self.Ws-self.Wus)/self.pcenv);
                    k[i] = self.kfs+((self.kus-self.kfs)/self.pcenv)*(T2[i]-(self.Tpc-self.pcenv));
                elif celltype2[i]==2 and T2[i]>self.Tpc:
                    C[i] = self.Cus;
                    k[i] = self.kus;
            
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
                    coeff3 = C[i-1]*dzbetween2[i-1]/self.dt;
                    coeff4 = coeff1+coeff2+coeff3;
                    Tcalcperm = np.append(Tcalcperm, ((coeff1*T2[i+1]+coeff2*T2[i-1]+coeff3*Told2[i])/coeff4));
                    #Tcalcperm[i] = (coeff1*T2[i+1]+coeff2*T2[i-1]+coeff3*Told2[i])/coeff4;
                T2_partial = np.append(self.T[self.ice+self.water],Tcalcperm[1:])
                T2 = np.append(Tcalcperm,self.Tbase[n]);
                error = max(abs(T2-Tit));
            self.T[self.ice+self.water-1:] = T2;
                        
            # subside permafrost
            print (len(range(self.water+self.ice,self.numcells)));
            
            for i in range (self.water+self.ice, self.numcells):
                if self.T[i]<self.Tpc:
                    break 
                elif T[i]>self.Tpc and celltype[i] == 1:
                    thawedspec[i] = 1;
                    thawed[i] = 1;
                    celltype[i] = 2;
                    #print 'is it here'
                elif T[i]>Tpc and celltype[i] == 2:
                    thawed[i] = 1;
                    #print 'where is it'
                        
        ##[Tbar-Tamp*math.sin(2*math.pi*x/periodyear) for x in t];
            numthawedspec = np.sum(self.thawedspec)+1;
            depthsubsidenew = self.z[numthawedspec]*self.excessice;
            #depthsubsidenew = self.z*self.excessice;
            depthsubside = self.depthsubside+depthsubsidenew;
            self.depthsubsiden[n] = self.depthsubside;
            cellsubside = np.int16(math.floor((self.depthsubside)/self.dz0));
            maxcellsubside = cellsubside; 
            #print range(self.ice,maxcellsubside)
            if self.depthsubside>0:    #depthsubsidenew>0
                for i in range(self.ice,maxcellsubside):
                    self.celltype[i] = 8;
                    #print 'p'
                    
            if (self.Tsurface[n]<self.Tpc):
                icecheckdepth = np.int16(0.33/self.dz0)+self.ice;    #z location of ~0.33 m into water column
                if icecheckdepth>self.ice+self.water:
                    icecheckdepth = self.ice+self.water-1;
                Twater = self.T[icecheckdepth];
                if self.water==1:
                    Twater = T[np.where(self.celltype==8)];
                    Twater = Twater[0];
                if Twater<self.Tpc and self.ice==0:
                    depthi = 0.01;
                
                #_____DOUBLE CHECK HERE_______
                depthmix = self.z[icecheckdepth]-depthi;
                if depthi>0:
                    ddepthidt_partial = ((self.Tpc-self.Ts0)*pow((depthi/self.ki),(-1))-self.hw*(Twater-self.Tpc));
                    ddepthidt = np.divide(ddepthidt_partial, self.rhoi*self.L)
                    ddepthi = ddepthidt*self.dt;
                #### ADJUSTMENT HERE
                if ddepthi>0.05:
                    ddepthi = 0.05;        
                    print('had to adjust');
                self.icegrowth[n] = ddepthi;
                depthi = depthi+ddepthi;
                if depthi<0.00001:
                    depthi = 0;
                depthw = self.depthsubside-depthi;
                
                if depthi>self.depthsubside:
                    depthi = self.depthsubside;
                    depthw = 0;
                    
            #     check if Ts>Tpc (if there is ice there already)
            elif self.Tsurface[n]>Tpc and depthi>0:
                Qm = (self.Tsurface[n]-Tpc)*(pow((self.dz0/self.ki),(-1)));
                dMdt = Qm/(L*rhoi);
                dM = dMdt*self.dt;
                icegrowth[n] = -dM;
                depthi = depthi-dM;
                if depthi<0:
                    depthi = 0;
            maxcelli = np.int16(round(depthi/self.dz0));
            self.celltype[0:self.ice+self.water] = 8;
            self.celltype[0:maxcelli] = 9;
            self.icethickness[n] = depthi;
            self.soilthickness[n] = self.deptht-self.depthsubside;
            self.waterthickness[n] = self.deptht-self.icethickness[n]-self.soilthickness[n];
            
            self.Tbot[n] = self.T[self.ice+self.water]-self.kelvin;
            self.T3m[n] = self.T[59]-self.kelvin;
            self.T5m[n] = self.T[99]-self.kelvin;
            self.T10m[n] = self.T[199]-self.kelvin;
            self.T25m[n] = self.T[214]-self.kelvin;
            self.T50m[n] = self.T[239]-self.kelvin;
            self.T100m[n] = self.T[289]-self.kelvin;
            if (self.tday[n]%365)==0 and self.tday[n]>364:
                yearnum = yearnum+1;
                Tbotavgann[yearnum] = np.mean(self.Tbot[n-364:n]);
                T3mavgann[yearnum] = np.mean(T3m[n-364:n]);
                T5mavgann[yearnum] = np.mean(T5m[n-364:n]);
                T10mavgann[yearnum] = np.mean(T10m[n-364:n]);
                T25mavgann[yearnum] = np.mean(T25m[n-364:n]);
                T50mavgann[yearnum] = np.mean(T50m[n-364:n]);
                T100mavgann[yearnum] = np.mean(T100m[n-364:n]);

            if n>=self.nt-364:
                daynum = daynum+1;
                Trecord[daynum, 0:] = T;
                    
            if((n % (365*50))==0):
                #    nplot = nplot+1
                fig_one = plt.figure(1)
                #    ice
                self.icethickness[n]
                #    water
                self.waterthickness[n]

                lakedepthx = np.arange(-20,25,1);
                lakedepth = self.icethickness[n]+self.waterthickness[n]*np.ones(np.size(lakedepthx));
                Tplot = self.T-self.kelvin;
                timeplot = self.tday[n]/365 #will stamp time in years on your screen)

                icedepth = self.icethickness[n]*np.ones(np.size(lakedepthx));
            
            #compute temperatures in water layer if no ice exist            

def printSomething():
    print ('This is working');

if __name__ == '__main__':
    classFunc = thaLakeModel();
    classFunc.initialize();
    classFunc.updateModel();
    print ("Loop ended")
