# -*- coding: utf-8 -*-
"""  Frost Number by Nelson and Outcalt 1983. DOI: 10.2307/1551363. http://www.jstor.org/stable/1551363
"""

import numpy as np
from permamodel.utils import model_input
from permamodel.components import perma_base

class gipl_model( perma_base.permafrost_component ):

    #-------------------------------------------------------------------
    _att_map = {
    # NOTE: this will change in the future
        'model_name':         'Geophysical_Intitute_Permafrost_Laboratory_model',
        'version':            '0.1',
        'author_name':        'Elchin Jafarov',
        'grid_type':          'none',
        'time_step_type':     'fixed',
        'step_method':        'implicit',
        #-------------------------------------------------------------
        'comp_name':          'gipl',
        'model_family':       'PermaModel',
        'cfg_extension':      '_gipl_model.cfg',
        'cmt_var_prefix':     '/input/',
        'gui_yaml_file':      '/input/gipl_model.yaml',
        'time_units':         'years' }

    _input_var_names = [
        'latitude',
        'longitude',
        'atmosphere_bottom_air__temperature',
        'snowpack__depth',
        'snowpack__density',
        'water-liquid__volumetric-water-content-soil',   
        'soil__thermal_conductivity_thawed',
        'soil__thermal_conductivity_frozen',
        'soil__heat_capacity_thawed',
        'soil__heat_capacity_frozen',
        'soil__unfrozen_water_parameter1',
        'soil__unfrozen_water_parameter2' ]
        
    _output_var_names = [
        'soil__temperature',                                  # Tps 
        'soil__active_layer_thickness' ]       				  # Zal      
      
    _var_name_map = {
    # NOTE: we need to look up for the corresponding standard names
        'latitude':                 						  'lat',
        'longitude':                 						  'lon',
        'atmosphere_bottom_air__temperature':             	  'T_air',
        'snowpack__depth':									  'sn_depth',
        'soil__thermal_conductivity_thawed':				  'k_th',
        'soil__thermal_conductivity_frozen':				  'k_fr',
        'soil__heat_capacity_thawed':						  'C_th',
        'soil__heat_capacity_frozen':						  'C_fr',
        'soil__unfrozen_water_parameter1':					  'aclv',
		'soil__unfrozen_water_parameter2':					  'bclv',
        'snowpack__density':                                  'rho_snow',
        'water-liquid__volumetric-water-content-soil':        'Wvol'}       

    _var_units_map = {
    # NOTE: Kang please complete the vegetation info both on var names and units
        'latitude':                 						  'lat',
        'longitude':                 						  'lon',
        'atmosphere_bottom_air__temperature':                 'deg_C',
        'atmosphere_bottom_air__temperature_amplitude':       'deg_C',
        'snowpack__depth':                                    'm',
        'snowpack__density':                                  'kg m-3',
        'water-liquid__volumetric-water-content-soil':        'm3 m-3' }      
    
    #-------------------------------------------------------------------
    def get_attribute(self, att_name):

        try:
            return self._att_map[ att_name.lower() ]
        except:
            print '###################################################'
            print ' ERROR: Could not find attribute: ' + att_name
            print '###################################################'
            print ' '

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
   
    #   get_var_units()
    #------------------------------------------------------------------- 
    def check_input_types(self):

        #--------------------------------------------------
        # Notes: rho_H2O, Cp_snow, rho_air and Cp_air are
        #        currently always scalars.
        #--------------------------------------------------        
        are_scalars = np.array([
                          self.is_scalar('lat'),
                          self.is_scalar('lon'),
                          self.is_scalar('A_air'),
                          self.is_scalar('sn_depth'),  
                          self.is_scalar('rho_snow'),
                          self.is_scalar('Wvol') ])

        self.ALL_SCALARS = np.all(are_scalars)
        
    #   check_input_types()
    #-------------------------------------------------------------------
    def open_input_files(self):
		# this function will work only if filename is not empty
        # I am temporarily comment T_air since in this case it is T_air_minmax
        self.T_air_file       = self.in_directory + self.T_air_file
      	self.A_air_file       = self.in_directory + self.A_air_file
        self.sn_depth_file 	  = self.in_directory + self.sn_depth_file
       	self.rho_snow_file    = self.in_directory + self.rho_snow_file
        self.Wvol_file     	  = self.in_directory + self.Wvol_file


        self.T_air_unit       = model_input.open_file(self.T_air_type,  self.T_air_file)
        self.A_air_unit       = model_input.open_file(self.A_air_type,  self.A_air_file)
        self.sn_depth_unit    = model_input.open_file(self.sn_depth_type,  self.sn_depth_file)
        self.rho_snow_unit    = model_input.open_file(self.rho_snow_type,  self.rho_snow_file)
        self.Wvol_unit        = model_input.open_file(self.Wvol_type,  self.Wvol_file)

		
    #   open_input_files()
    #-------------------------------------------------------------------  
    def read_input_files(self):

        #rti = self.rti # has a problem with loading rti: do not know where its been initialized

        #-------------------------------------------------------
        # All grids are assumed to have a data type of Float32.
        #-------------------------------------------------------
        T_air = model_input.read_next_modified(self.T_air_unit, self.T_air_type)
        if (T_air != None): self.T_air = T_air
        print T_air
        print 'OKOK'

        #T0 = model_input.read_next(self.T0_unit, self.T0_type, rti)
        #if (T0 != None): self.T0 = T0

        #rho_snow = model_input.read_next(self.rho_snow_unit, self.rho_snow_type, rti)
        #if (rho_snow != None): self.rho_snow = rho_snow

        #h0_snow = model_input.read_next(self.h0_snow_unit, self.h0_snow_type, rti)
        #if (h0_snow != None): self.h0_snow = h0_snow
        
       # h0_swe = model_input.read_next(self.h0_swe_unit, self.h0_swe_type, rti)
        #if (h0_swe != None): self.h0_swe = h0_swe
        
    #   read_input_files()       
    #-------------------------------------------------------------------  



    #-------------------------------------------------------------------
    def initialize(self, cfg_file=None, mode="nondriver",
                   SILENT=False):

        if not(SILENT):
            print ' '
            print 'GIPL component: Initializing...'
            
        self.status   = 'initializing'  # (OpenMI 2.0 convention)
        self.mode     = mode
        self.cfg_file = cfg_file    
        if cfg_file:
            fid = open(cfg_file,'r')
            std_dummy = fid.readline()
            snow_depth_file = fid.readline().rstrip()
            upper_bnd_file = fid.readline().rstrip()
            init_file = fid.readline().rstrip()
            grid_file = fid.readline().rstrip()
            grid_out_file = fid.readline().rstrip()
            therm_prop_file = fid.readline().rstrip()
            fid.close()
            self.snow_depth_file=snow_depth_file
            self.upper_bnd_file=upper_bnd_file
            self.init_file=init_file
            self.grid_file=grid_file
            self.grid_out_file=grid_out_file
            self.therm_prop_file=therm_prop_file 
        
        self.initilize_consts()
        self.load_initial_conditions()
        self.load_upper_bnd_condition()
        self.load_z_grid()
        self.load_z_grid_out()
        self.load_therm_prop()
        self.load_snow_depth()
        # starting from the calling oreder does matter
        self.interpolate_init_temp()
        self.layer_indexing()
        self.nondimensionalization()
    	
    #   initialize()  
    #-------------------------------------------------------------------  

    def initilize_consts(self):
        self.taum=0.1 
        self.tmin=0.001   #???, smallest time step division
   
        self.E1=0.01 
        self.UWK=0.01    #Iterate stopage constants
        #E=0.0;			#Smoothing paramenter
        
        self.perl=86400 # number of seconds per day hour*min*sec=24*60*60

        self.nmean=1    #average eqaul to one day, 365 average over a year 
        self.itfin=0  
        self.itfinl=1   #beginig and end of the year
        self.step=1;    #time step=1 day
        self.tini=1;    #initial time step (starting day)
        self.lbound=2;  # heat flux at the lower boundary
        self.itmax=5    #maximum number of allowed interations
        
    #   initilize_consts  
    #-------------------------------------------------------------------         
    def load_initial_conditions(self):
        data = np.loadtxt(self.init_file, skiprows=1,unpack=False)
        # initial depth and temprature
        i_depth = data[:,0] # dinit
        i_temp = data[:,1]  # tinit
        self.i_depth=i_depth 
        self.i_temp=i_temp
        self.i_time=0
        return 
        
    def load_upper_bnd_condition(self):
        data = np.loadtxt(self.upper_bnd_file, skiprows=1,unpack=False)
        # upper boundary time and temperature
        u_time = data[:,0] # tmpdayNo
        u_temp = data[:,1] # airTemp[i]
        self.T_curr=u_temp[0]
        self.u_time=u_time
        self.u_temp=u_temp
        return 
    
    def load_z_grid(self):
        data = np.loadtxt(self.grid_file, skiprows=1,unpack=False)
        # vertical grid index and depth
        z = data[:,1] # ZDEPTH[N]
        self.N=z.size
        self.z=z
        return 
    
    def load_z_grid_out(self):
        data = np.loadtxt(self.grid_out_file, skiprows=2,unpack=False)
        # vertical grid index that needs to be saved 
        z_idx = data[:,1] # iout[dout]
        self.dout=z_idx.size
        self.z_idx=z_idx.astype(int)
        return 

    def load_therm_prop(self):
        data5 = np.loadtxt(self.therm_prop_file, skiprows=8,unpack=False)
        # thermal properties
        Thick=data5[:,0]
        Wvol=data5[:,1]; aclv=data5[:,2]; bclv=data5[:,3]
        #cclv[nLayers]={0,0,0,0,0,0};
        Cond_Th=data5[:,4]; Cond_Fr=data5[:,5]; C_th=data5[:,6]; C_fr=data5[:,7]
        self.Thick=Thick
        self.Wvol=Wvol
        self.aclv=aclv
        self.bclv=bclv
        self.k_th=Cond_Th
        self.k_fr=Cond_Fr
        self.C_th=C_th
        self.C_fr=C_fr
        return 
    
    def load_snow_depth(self):
        data = np.loadtxt(self.snow_depth_file, skiprows=1,unpack=False)
        # snow day and depth
        sn_time = data[:,0]  # tmpdayNo
        sn_depth = data[:,1] # snowDepth
        
        self.hsnow=sn_depth[0]
        self.sn_time=sn_time
        self.sn_depth=sn_depth
        
        return     
        
    def interpolate_init_temp(self):
        M=self.i_temp.size
        d=self.i_depth
        temp=self.i_temp
        # need to make surethat upper and lower bound of the grid match
        if temp[M-1]!=self.z[self.N-1]:
            d=np.append(d, self.z[self.N-1])
            temp=np.append(temp, temp[M-1]) 
        
        from scipy import interpolate
        f_interp = interpolate.interp1d(d,temp,'linear')
        T=f_interp(self.z)
        self.T=T
        self.U1=T
    
        return 

    def layer_indexing(self):        
        nLayers=self.Thick.size
        numl=np.ones(self.N)*(nLayers-1)

        for i in range(nLayers-1,0,-1):
            numl[np.where(self.z<self.Thick[i-1])]=i-1
        numl[np.where(self.z<0)]=-1

        self.numl=numl
        self.nLayers=nLayers         
        
        return 
    
    def soilThermalConductivity(self,j,temp):
        ln=self.numl[j]; 
        if (self.z[j] <= -self.hsnow):
            k_eff=1.e4
        elif (self.z[j] <= 0):  # assuming ground surface starts at 0 meters
            k_eff=0.18
        else:
            theta=self.unfrWater(temp,ln)/self.Wvol[ln]
            k_eff=self.k_th[ln]**theta*self.k_fr[ln]**(1-theta)

        return k_eff
    
    def unfrWater(self,temp,i):
        E=0 # smoothing factor
        if self.Wvol[i]==0.0:
            self.Wvol[i]==0.01
        
        tfp=self.tfpw[i]
        tt=abs(temp)
        #print tfp ,tt ,temp
        if (temp <= tfp-E):
            unWater=self.aclv[i]*tt**(self.bclv[i])
        elif (temp > tfp):
            unWater=self.Wvol[i]
        else:
            tt = abs(tfp-E)
            unWater=self.aclv[i]*tt**(self.bclv[i])
            unWater=unWater+(Wvol[i]-unWater)*(temper+E-tfp)/E;
        
        return unWater
        
    def dunfrWater(self,temp,i):
        E=0 # smoothing factor
        if self.Wvol[i]==0.0:
            self.Wvol[i]==0.01
        
        tfp=self.tfpw[i]
        tt=abs(temp)
        if (temp <= tfp-E):
            dunWater=-self.bclv[i]*self.aclv[i]*tt**(self.bclv[i]-1)
        elif (temp > tfp):
            dunWater=0.0
        else:
            tt = abs(tfp-E)
            dunWater=self.aclv[i]*tt**(self.bclv[i])
            dunWater=(Wvol[i]-dunWater)/E;
        
        return dunWater
        
    def C_Q(self,temp1,temp2,i1,i2):
        h=1/(temp1-temp2)
        if (abs(temp1-temp2)<1.e-6):
            Cphase=0.5*(self.dunfrWater(temp1,i1)+self.dunfrWater(temp2,i2))
        else:
            if (i1!=i2):
                diff1=self.unfrWater(temp1,i1)-self.unfrWater(temp2,i1)
                diff2=self.unfrWater(temp1,i2)-self.unfrWater(temp2,i2)
                Cphase=0.5*(h*diff1+h*diff2)
            else:
                Cphase=h*(self.unfrWater(temp1,i1)-self.unfrWater(temp2,i2))
                
        Cphase=self.Qphase*abs(Cphase)
        
        return Cphase
    
    def Capp(self, j, V):
        ln=self.numl[j]; 
        if (self.z[j] <= -self.sn_depth[j]):
            CAP=self.C_air
        elif (self.z[j] <= 0):  # assuming ground surface starts at 0 meters
            CAP=self.C_snow
        else:
            wc=self.unfrWater(V[j],ln)/self.Wvol[ln]
            #print V[j], ln, self.Wvol[ln], wc
            CAP=self.C_th[ln]*wc+self.C_fr[ln]*(1-wc);
            if((j>0) and (j<self.N-1)):
                ww1=(V[j-1]+V[j])/2; 	nn1=self.numl[j-1];
                ww2=V[j]; 	  		nn2=self.numl[j];
                CAP=CAP+self.C_Q(ww1,ww2,nn1,nn2)*self.hy[j]/(self.hy[j+1]+self.hy[j]);
                ww1=V[j]; 	  		nn1=self.numl[j];
                ww2=(V[j+1]+V[j])/2; 	nn2=self.numl[j+1];
                CAP=CAP+self.C_Q(ww1,ww2,nn1,nn2)*self.hy[j+1]/(self.hy[j+1]+self.hy[j]);
            elif (j==0):
                ww1=V[j]; 	  		nn1=self.numl[j];
                ww2=(V[j+1]+V[j])/2; 	nn2=self.numl[j+1];
                CAP=CAP+self.C_Q(ww1,ww2,nn1,nn2);
            elif (j==self.N-1):
                ww1=(V[j-1]+V[j])/2; 	nn1=self.numl[j-1];
                ww2=V[j]; 	  		nn2=self.numl[j];
                CAP=CAP+self.C_Q(ww1,ww2,nn1,nn2);
            
            #print wc, CAP, ww1,ww2,nn1,nn2
                
        return CAP
    

            
    def nondimensionalization(self):
        
        self.C_snow=840000.0
        self.C_air=840000.0   #Heat capacity snow, air
        self.grad=0; # assuming 0 heat flux at the lower boundary
        
        y=np.zeros(self.N)
        hy=np.zeros(self.N)
        
        y=self.z/self.z[self.N-1]
        hy[1:self.N-1:1]=y[1:self.N-1:1]-y[0:self.N-2:1]
        hy[0]=hy[1];	hy[self.N-1]=hy[self.N-2]
        
        hcscale=self.z[self.N-1]*self.z[self.N-1]/self.perl
        
        self.C_fr=self.C_fr*hcscale
        self.C_th=self.C_th*hcscale
        self.tfpw=-(self.Wvol/self.aclv)**(1/self.bclv)
        self.C_snow=self.C_snow*hcscale
        self.C_air=self.C_air*hcscale
        self.Qphase=hcscale*333200012.20703;
        self.grad=self.grad*self.z[self.N-1]
        self.y=y
        self.hy=hy
        
        return 

    def GaussianElimination(self):
        #tau=self.taum# temporary 
        hz=self.hy
        #print self.N-1     
        #UU=np.zeros(self.N)#temporary
        #self.U2=self.T#temporary
        #j=1
        #d=self.Capp(j,self.T1)/self.tau
        #print j, self.T1[j], self.Capp(j,self.T1), self.tau, d

        for j in range(1,self.N-1):
            d=self.Capp(j,self.T1)/self.tau
            #print j, self.T1[j], self.Capp(j,self.T1), tau, d
            k_eff=self.soilThermalConductivity(j,self.T1[j])
            #print j,k_eff,self.Capp(j,self.T)
            a=2*k_eff/(hz[j]*(hz[j]+hz[j+1]))
            #print j, k_eff, self.soilThermalConductivity(j+1,self.T1[j+1])
            k_eff=self.soilThermalConductivity(j+1,self.T1[j+1])
            b=2*k_eff/(hz[j+1]*(hz[j]+hz[j+1]));
            c=a+b+d;
            self.alf[j+1]=b/(c-a*self.alf[j]);
            self.bet[j+1]=(a*self.bet[j]+d*self.UU[j])/(c-a*self.alf[j]);
            #print j, self.soilThermalConductivity(j,self.T[j]), self.Capp(j,self.T)
            #print j,d,a,b,c,alf[j],bet[j]
            
        rab1=self.soilThermalConductivity(self.N-1,self.T1[self.N-1])
        rab2=self.Capp(self.N-1,self.T1);
        dhz=hz[self.N-1]*hz[self.N-1];
        akapa2=2*rab1/(((rab2*dhz)/self.tau+2*rab1));
        q2=rab1*self.grad;
        amu2=(self.UU[self.N-1]*rab2/self.tau+2*q2/hz[self.N-1])/(rab2/self.tau+2*rab1/dhz);
        #print rab1,rab2,akapa2,q2,amu2 
        
        if (abs(akapa2)>1.0):
            print "YOU CAN NOT APPLY PROGONKA ON U1 - CHANGE STEPS"
        
        if (self.lbound==2):
            self.U2[self.N-1]=(amu2+akapa2*self.bet[self.N-1])/(1-self.alf[self.N-1]*akapa2)
        else:
            self.U2[self.N-1]=self.grad
        #print  U2[self.N-1]

        for j in range(2,self.N+1):
        	self.U2[self.N-j]=self.alf[self.N-j+1]*self.U2[self.N-j+1]+self.bet[self.N-j+1]
        	#print self.alf[self.N-j+1], self.U2[self.N-j+1], self.bet[self.N-j+1],  self.U2[j]
            
        
        return
    
    def Iterate(self):
        #t1=0 #temporary
        #taum=0.0 #temporary
        #tau=taum #temporary
        
        self.t=self.t1+self.tau;
        tlet=self.t;
        it=1; 
        
        self.T1=self.UU
        hz=self.hy; d=0
        
        self.alf=np.zeros(self.N); self.bet=np.zeros(self.N)        

        if(it>self.itmax):
            self.tau=self.tau/2; self.tau1=-1.0;
            
        it=1; self.alf[1]=0; self.bet[1]=self.T_curr;
        
        #print self.t, tlet, self.tau, self.alf[1], self.bet[1]
        
        while (it <= self.itmax and self.tau>self.tmin ): # will add it later
       		#print self.T1[j],self.U2[j]
       		#print 'T43=',np.around(self.T1[43], decimals=4)
        
        	self.GaussianElimination()
        	#print 'T43=',np.around(self.T1[43], decimals=4), np.around(self.U2[43], decimals=4)
        	
        	if (self.tau>self.tmin):
        		#print 'it=', it, self.tau, self.tmin
        		
        		for j in range(0,self.N): #print "skipping for for now" 
        			#j=0
        			eey=self.unfrWater(self.U2[j],self.numl[j])/self.Wvol[self.numl[j]];
        			eey1=self.unfrWater(self.T1[j],self.numl[j])/self.Wvol[self.numl[j]];
        			abs1=abs(eey-eey1);
        			abs2=abs(self.T1[j]-self.U2[j]); 
        			#print 'j=',j, np.around(self.T1[43], decimals=3)
        			#print 'j=',j, np.around(self.T1[j], decimals=3), np.around(self.U2[j], decimals=3)
        			#print j,eey,eey1,abs1,abs2
        			
        			#if (j>=43):
        				#print "OK1"
        				#d=self.Capp(j,self.T1)/self.tau
        				#print j,self.T1[j],self.Capp(j,self.T1),self.tau,d
        				#print self.T[j-1],self.T[j],self.T[j+1]		 	 
        				#k_eff=self.soilThermalConductivity(j,self.T1[j])
        				#a=2*k_eff/(hz[j]*(hz[j]+hz[j+1]))
        				#k_eff=self.soilThermalConductivity(j+1,self.T1[j+1])
        				#b=2*k_eff/(hz[j+1]*(hz[j]+hz[j+1]));
        				#c=a+b+d;
        				#self.alf[j+1]=b/(c-a*self.alf[j]);
        				#self.bet[j+1]=(a*self.bet[j]+d*self.UU[j])/(c-a*self.alf[j]);
        				#print j,d,a,b,c,self.alf[j],self.bet[j]
        				#return
            			
        			if((abs1>self.UWK) or (abs2>self.E1)):
        				#print "OK2"
        				self.T1=np.copy(self.U2)
        				it=it+1
        		
        				if (it>self.itmax):
        					self.tau=self.tau/2; self.tau1=-1.0;
        					self.t=self.t1+self.tau;
        					tlet=self.t;
        					self.T1=np.copy(self.UU)
        					it=1; self.alf[1]=0; self.bet[1]=self.T_curr;
        				break
        		#if (j>=25):
        		#	break
        		
        		if (j>=self.N-1):
        			#print j,"OKOKOK"
        			break
        
        return 
        
    def update(self):
    # EJ 05/16/16 Note: 
    # I am completely overriding this function fro this model
    # This needs to be in the future versions of the model
        
        sw=0; ttt=0; self.taum=0.1; self.tmin=0.001; 
        self.UU=np.zeros(self.N)
        self.U2=np.zeros(self.N)
        #print out current ground tempratures
        res=np.zeros(self.z_idx.size)
        res=self.U1[self.z_idx]
        #print np.around(res, decimals=3)
        #open result file and save the results
        
        self.t1=ttt; self.tau1=-1.0; self.tau=self.taum;
        #print t1, tau1, tau
        
        self.UU=np.copy(self.U1)
        count=1
        while sw<0.1: # skpping it for now
        	self.Iterate()
        	if (self.t<ttt+self.step-1.e-12):
        		#print "ok1"
        		self.t1=self.t;
        		self.UU=np.copy(self.U2)
        		#print "OK", self.t1,self.t,ttt
        	
        		if (self.tau1>0):
        			if(self.tau<self.taum):
        				self.tau=2*self.tau
        				self.tau1=-1.0
        		else: 
        			self.tau1=1.0
        		sw=0
        
        	elif(self.t>ttt+self.step+1.e-12): 
        		#print "ok2"
        		#print "tau=", np.around(self.t1, decimals=5), np.around(self.tau, decimals=5)
        		self.tau=ttt+self.step-self.t1
        		#print "tau=", np.around(self.t1, decimals=5), np.around(self.tau, decimals=5)
        		sw=0
        
        	else:
        		#print "ok3"
        		self.UU=np.copy(self.U2)
        		sw=1
        	#print count, sw, np.around(self.t, decimals=5)
        	count=count+1
        	#if count>=43:
        	#	break
        	
        self.U1=np.copy(self.UU)
        
        #for i in range(0,self.N):
        #	print i,np.around(self.U1[i], decimals=5)
        res=self.U1[self.z_idx]
        print "U1"
        print np.around(res, decimals=5)
        
        # update internal time 
        self.i_time=self.i_time+1
        self.T_curr=self.u_temp[self.i_time]
        self.hsnow=self.sn_depth[self.i_time]

    #   update()
    #-------------------------------------------------------------------
    
    def update_ALT(self):
 
        #---------------------------------------------------------
        # 		coming up
        #--------------------------------------------------
		print 'In Progress'        
        		
    #   update_ALT()
    #-------------------------------------------------------------------
    def close_input_files(self):

        if (self.T_air_type     != 'Scalar'): self.T_air_unit.close() 
      	if (self.A_air_type     != 'Scalar'): self.A_air_unit.close() 
        if (self.h_snow_type    != 'Scalar'): self.h_snow_unit.close() 
       	if (self.rho_snow_type  != 'Scalar'): self.rho_snow_unit.close() 
        if (self.vwc_H2O_type   != 'Scalar'): self.vwc_H2O_unit.close() 


    #   close_input_files()    
    #-------------------------------------------------------------------