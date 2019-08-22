import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import ode


class FineStructureHam():
    
    def __init__(self):
        
        self.mu_b = 1.3e-3 # GHz/G
        
        # initialize Hamiltonian parameters
        self.HamDict= {}
        
        self.HamDict['Dgs'] = 0.0
        self.HamDict['Des'] = 0.0

        self.HamDict['g_gs'] = [2.0,2.0]
        self.HamDict['g_es'] = [2.0,2.0]

        self.HamDict['lambdaz'] = 0.0
        self.HamDict['Delta1'] = 0.0
        self.HamDict['Delta2'] = 0.0

        self.HamDict['Bvec'] = [0,0,0]

        self.HamDict['strain_delta'] = 0.0
        self.HamDict['strain_alpha'] = 0.0
        
        #
        self.SpecParams = {}
        self.SpecParams['spec_resolution'] = 0.05
        self.SpecParams['fwhm'] = 0.5
        
        self.freqoffset = 475.0e3
        
        self.setparam = False
        self.calculatedevals = False
        self.hambuilt = False
        
        # Now define some operators
        
        self.Sx = 1/np.sqrt(2)*np.asarray([[0,1,0],[1,0,1],[0,1,0]])
        self.Sy = 1/(np.sqrt(2)*1j)*np.asarray([[0,1,0],[-1,0,1],[0,-1,0]])
        self.Sz = np.asarray([[1,0,0],[0,0,0],[0,0,-1]])
        self.Sp = self.Sx + 1j*self.Sy
        self.Sm = self.Sx - 1j*self.Sy

        self.sigma_x = np.asarray([[0,1],[1,0]])
        self.sigma_y = 1j*np.asarray([[0,1],[-1,0]])
        self.sigma_z = np.asarray([[-1,0],[0,1]])
        self.sigma_p = self.sigma_z + 1j*self.sigma_x
        self.sigma_m = self.sigma_z - 1j*self.sigma_x


        self.I2 = np.identity(2)
        self.I3 = np.identity(3)
        
    def GetHelp(self):
        print('------------------------------------------------- \n'+
            '                      Hello!       \n'+
            '------------------------------------------------- \n'+
            'I am the helper function for the FineStructureHam \n'+
            'class. The purpose of this class is to generate \n'+
            'functions required to simulate the Hamiltonian of \n'+
            'the NV center, and other isoelectronic systems. \n'+
             '-------------------------------------------------  \n'+
             'For help with the model and model parameters, \n'+
             'please use GetHamHelp().  \n'+
             '------------------------------------------------- ')
        
    def GetHamHelp(self):
        print('------------------------------------------------- \n'+
            '          Fine-Structure Hamiltonian Help      \n'+
            '------------------------------------------------- \n'+
            'This model uses two dictionaries - one to calculate\n'+
            'the Hamiltonian, one to calculate the optical spectrum\n'+
            '------------------------------------------------- \n'+
            '                     HamDict \n'
            '------------------------------------------------- \n'+
            'Dgs, Des are the ground (excited) state ZFS in GHz\n'+
            'lambdaz is the axial spin-orbit strength in GHz\n'+
            'Delta1, Delta2 are the axial and non-axial spin-spin interactions in GHz \n'+
            'Bvec is the vector magnetic field in Gauss \n'+
            'strain_alpha, strain_delta are the strain angle (deg) and strain \n'+
            'splitting (GHz) \n'
            'The function MakeNVHam autopopulates this dictionary with literature values\n' +
            'for the NV center \n' +
            '------------------------------------------------- \n'+
            '                     SpecParams \n'
            '------------------------------------------------- \n'+
            'spec_resolution is output resolution of the calculated spectrum \n' +
             'fwhm is the FWHM of a lorentzian lineshape for the peaks')
        
    def MakeNVHam(self):
        self.HamDict['Dgs'] = 2.88
        self.HamDict['Des'] = 1.42

        self.HamDict['g_gs'] = [2.03,2.01]
        self.HamDict['g_es'] = [2.03,2.01]

        self.HamDict['lambdaz'] = 5.5
        self.HamDict['Delta1'] = 3.0
        self.HamDict['Delta2'] = 0.0

        self.HamDict['Bvec'] = [0,0,0]

        self.HamDict['strain_delta'] = 0.0
        self.HamDict['strain_alpha'] = 0.0
        
        self.setparam = True
        self.hambuilt = False
        
    def SetParams(self,newdict):
        self.HamDict.update(newdict)
        self.setparam = True
        self.hambuilt = False
        
    def BuildHam(self):
        # ground state Hamiltonian
        self.Hss_g = self.HamDict['Dgs']*(self.Sz@self.Sz-2/3*self.I3)
        
        self.Hz_g = (self.HamDict['g_gs'][1]*self.mu_b*(self.HamDict['Bvec'][0]*self.Sx + self.HamDict['Bvec'][1]*self.Sy) 
                     + self.HamDict['g_gs'][0]*self.mu_b*self.HamDict['Bvec'][2]*self.Sz)
            
        self.Htot_g = self.Hss_g + self.Hz_g
        
        # excited state Hamiltonian
        self.Hso = -1*self.HamDict['lambdaz']*np.kron(self.sigma_y,self.I3)@np.kron(self.I2,self.Sz)

        self.Hss = (self.HamDict['Des']*np.kron(self.I2,(self.Sz@self.Sz-2/3*self.I3))
                    -1*self.HamDict['Delta1']/4*(np.kron(self.I2,self.Sp@self.Sp)@np.kron(self.sigma_m,self.I3)+np.kron(self.I2,self.Sm@self.Sm)@np.kron(self.sigma_p,self.I3))
                    +self.HamDict['Delta2']/(2*np.sqrt(2))*(np.kron(self.I2,self.Sp@self.Sz+self.Sz@self.Sp)@np.kron(self.sigma_p,self.I3) + np.kron(self.I2,self.Sm@self.Sz+self.Sz@self.Sm)@np.kron(self.sigma_m,self.I3)))


        self.Hz = np.kron(self.I2,
                          self.HamDict['g_es'][1]*self.mu_b*(self.HamDict['Bvec'][0]*self.Sx + self.HamDict['Bvec'][1]*self.Sy) 
                          + self.HamDict['g_es'][0]*self.mu_b*self.HamDict['Bvec'][2]*self.Sz)

        self.Hs = np.kron(self.HamDict['strain_delta']*np.cos(self.HamDict['strain_alpha']/180*np.pi)*self.sigma_z
                     +self.HamDict['strain_delta']*np.sin(self.HamDict['strain_alpha']/180*np.pi)*self.sigma_x,
                     self.I3)
        
        self.Htot = self.Hso + self.Hss + self.Hz + self.Hs
        
        self.hambuilt = True
        self.calculatedevals = False
        
    def CalcEvals(self):
        if self.setparam:
            if self.hambuilt:
                self.eigvals,self.eigvecs = np.linalg.eigh(self.Htot)
                self.eigvals_g,self.eigvecs_g = np.linalg.eigh(self.Htot_g)
            else:
                print("You haven't rebuilt the Hamiltonian since you updated the parameters!")
        else:
            print("You haven't set the Hamiltonian parameters! Set them first!")
        self.calculatedevals = True
        return self.eigvals,self.eigvecs
    
    def LorentzianFunction(self,x,scale,offset,fwhm):    
        return scale*(0.5*fwhm)**2/((x-offset)**2 + (0.5*fwhm)**2)
    
    def GenFreqAx(self,transfreq):    
        minval = transfreq.min() - 5*self.SpecParams['fwhm']
        maxval = transfreq.max() + 5*self.SpecParams['fwhm']
        npts = np.ceil((maxval-minval)/self.SpecParams['spec_resolution'])
        return np.linspace(minval,maxval,npts)

    def CalcSpectrum(self):
        if self.calculatedevals:
            self.Hamfull = np.zeros((9,9)) +0*1j
            self.Hamfull[0:3,0:3] = self.Htot_g
            self.Hamfull[3:,3:] = self.Htot + self.freqoffset*np.identity(6)

            self.TDmat = np.zeros((9,9))
            self.TDmat[(0,0,1,1,2,2),(3,6,4,7,5,8)] = 1.0
            self.TDmat = self.TDmat + (self.TDmat).T

            d,v = np.linalg.eigh(self.Hamfull)
            self.transformedTDmat = (np.conj(v)).T@self.TDmat@v
            
            transition_freq = np.asarray([(d[3:]-dval)-self.freqoffset for dval in d[0:3]])
            self.freqax = self.GenFreqAx(transition_freq)
            transition_inten = np.abs(self.transformedTDmat[0:3,:])**2
            transition_inten[transition_inten < 1e-16] = 0
            
            self.spec = []
            for freqvals,intenvals in zip(transition_freq,transition_inten[0:3,3:]):
                self.spec.append(np.sum(np.asarray([self.LorentzianFunction(self.freqax,inten,freq,self.SpecParams['fwhm']) for freq,inten in zip(freqvals,intenvals)]),0))
            self.spec = np.asarray(self.spec)
            return self.freqax,self.spec

        else:
            print("You need to recalculate the eigenstates before calculating the spectrum!")
            return np.nan,np.nan

## build model class
## order of states - GA, GB, EA, EB, Sh

class FiveLevelModel():
    
    def __init__(self):
        self.ratedict = {}
        
        # optical excitation rates from GA to EA, EB respectively (spin conseving/non-conserving)
        self.ratedict['kupA'] = 0.0
        self.ratedict['kupAB'] = 0.0
        
        self.ratedict['kupB'] = 0.0 # same but for GB
        self.ratedict['kupBA'] = 0.0
        
        self.ratedict['rateT1'] = 0.0 # spin lattice relaxation RATE
        
        # branching ratio for the radiative relaxation; this should actually be related to kupAB, 
        # but it's easier to implement the spectral stuff without it. Make the ratio of the peaks
        # in the spectrum  phi to be consistent
        self.ratedict['phi'] = 0.0 
        
        self.ratedict['phi_s'] = 1.0 # branching ratio to the ground state. If 1, no off-resonant polarization
        
        self.ratedict['kradA'] = 0.0 # radiative rate for A
        self.ratedict['kradB'] = 0.0 # radiative rate for B
        
        self.ratedict['ISC'] = 0.0 # intersystem crossing rate from E to Sh
        self.ratedict['ShG'] = 0.0 # relaxation from shelving to ground state
        
        self.ratemat = np.zeros((5,5))
        
        # now move on to the ODE stuff
        self.UpdatedFlag = False
        
        self.ODEdict = {}
        self.ODEdict['t0'] = 0.0
        self.ODEdict['t1'] = 1.0
        self.ODEdict['dt'] = 0.1
        self.ODEdict['initial_pops'] = np.asarray([0.0,0.0,0.0,0.0,1.0])
        self.ODEdict['timeax'] = np.arange(self.ODEdict['t0'],self.ODEdict['t1'],self.ODEdict['dt'])
        
        self.Results = {}
        self.Results['Labels'] = ['GA','GB','EA','EB','Sh']
        self.Results['pops_out'] = []
        self.Results['timeax'] = []
    
    def GetHelp(self):
        print('------------------------------------- \n'+
            '       Hello!       \n'+
            '------------------------------------- \n'+
            'I am the helper function for the FiveLevelModel \n'+
            'class. The purpose of this class is to generate \n'+
            'functions required to simulate the dynamcis of \n'+
            'a Five-level system commonly used to describe the\n'+
            'behavior of NV and other color centers.'+
             '------------------------------------- \n'+
             'For help with the model and model parameters, \n'+
             'please use GetModelHelp().  \n'+
             '------------------------------------- \n'+
             'For help with the ODE setup and parameters, \n'+
             'please use GetODEHelp().')
	
    def GetModelHelp(self):
        print('------------------------------------- \n'+
            '       Five-Level Model Help       \n'+
            '------------------------------------- \n'+
            'Ordering of states is GA,GB,EA,EB,Sh \n'+
            'Parameters are: \n'+
            'kupA, kupAB; kupB, kupBA - excitation rates for GA->EA,EB; GB->EB,EA \n'+
            'kradA; kradB - radiative relaxation rates for A;B \n'
            'phi - branching ratio of radiative rate between A,B \n'+
            'rateT1 - the spin-lattice relaxation rate \n'+
            'ISC - intersystem crossing rate from the excited state \n'+
            'ShG - rate from Sh to GA,GB \n'+
            'phi_s - branching ratio of ShG to A,B. 0 is all A, 1.0 is equal A,B \n'+
             '------------------------------------- \n')
        
    def GetODEHelp(self):
        print('------------------------------------- \n'+
            '     Five-Level Model ODE Help       \n'+
            '------------------------------------- \n'+
            't0,t1,dt - Start, End, Step times \n'+
            'initial_pops - initial populations \n'+
			'outputs are stored in the .Results of the class\n'+
            '------------------------------------- \n')
        
    def SetRateDict(self,newdict):
        if 'initial_pops' in newdict:
            newdict['initial_pops'] = np.asarray(newdict['initial_pops'])
        self.ratedict.update(newdict)

    def BuildRateMatrix(self):
        self.ratemat = np.zeros((5,5))

        # Eq for GA
        self.ratemat[0,0] = -self.ratedict['kupA']-self.ratedict['rateT1'] - self.ratedict['phi']*self.ratedict['kradB'] - self.ratedict['kupAB']
        self.ratemat[0,1] = self.ratedict['rateT1']
        self.ratemat[0,2] = self.ratedict['kupA'] + self.ratedict['kradA']
        self.ratemat[0,3] = (self.ratedict['kupAB'] + self.ratedict['phi']*self.ratedict['kradB'])
        self.ratemat[0,4] = self.ratedict['ShG']

        # Eq for GB
        self.ratemat[1,0] = self.ratedict['rateT1']
        self.ratemat[1,1] = -self.ratedict['kupB']-self.ratedict['rateT1'] - self.ratedict['phi']*self.ratedict['kradA'] - self.ratedict['kupBA']
        self.ratemat[1,2] = (self.ratedict['kupBA'] + self.ratedict['phi']*self.ratedict['kradA'])
        self.ratemat[1,3] = self.ratedict['kupB'] + self.ratedict['kradB']
        self.ratemat[1,4] = self.ratedict['phi_s']*self.ratedict['ShG']

        # Eq for EA
        self.ratemat[2,0] = self.ratedict['kupA']
        self.ratemat[2,1] = self.ratedict['phi']*self.ratedict['kradA'] + self.ratedict['kupBA']
        self.ratemat[2,2] = -self.ratedict['kupA'] - self.ratedict['kradA'] - (self.ratedict['kupBA'] + self.ratedict['phi']*self.ratedict['kradA']) -self.ratedict['ISC']

        # Eq for EB
        self.ratemat[3,0] = self.ratedict['phi']*self.ratedict['kradB'] + self.ratedict['kupAB']
        self.ratemat[3,1] = self.ratedict['kupB']
        self.ratemat[3,3] = -self.ratedict['kupB'] - self.ratedict['kradB'] - (self.ratedict['kupAB'] + self.ratedict['phi']*self.ratedict['kradB']) -self.ratedict['ISC']

        # Eq for Sh
        self.ratemat[4,2] = self.ratedict['ISC']
        self.ratemat[4,3] = self.ratedict['ISC']
        self.ratemat[4,4] = -(1+self.ratedict['phi_s'])*self.ratedict['ShG']
        
    def SetSolverParamaters(self,newdict):
        self.UpdatedFlag = True
        self.ODEdict.update(newdict)
    
    def SolveRateEquations(self):
        if self.UpdatedFlag:
            self.ODEdict['timeax'] = np.arange(self.ODEdict['t0'],self.ODEdict['t1'],self.ODEdict['dt'])
            self.Results['timeax'] = self.ODEdict['timeax']

            r = ode(lambda t,pops,ratemat: np.matmul(ratemat,pops)).set_integrator('vode', method='bdf')
            r.set_initial_value(self.ODEdict['initial_pops'], 0).set_f_params(self.ratemat)

            self.Results['pops_out'] = np.zeros((self.ODEdict['timeax'].size,self.ODEdict['initial_pops'].size))
            kindex = 0

            while r.successful() and r.t < (self.ODEdict['t1']-self.ODEdict['dt']):
                r.integrate(r.t+self.ODEdict['dt'])
                self.Results['pops_out'][kindex,:] = r.y
                kindex += 1
        else:
            print(r"You haven't set the ODE parameters! Do that first!")
