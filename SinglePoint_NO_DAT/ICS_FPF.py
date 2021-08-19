#-------------------------------------------------------------#
# Python Code for RMS bandwidth + Collimated Flux Calculation #
#                        Joe Crone 16/04/2021                 #
#-------------------------------------------------------------#
import csv
import scipy.constants as scipyconst
import scipy.integrate as scipyint
import numpy as np
import time

tic = time.perf_counter()

# SETUP
# Units
MeV = scipyconst.mega
keV = scipyconst.kilo
nm = scipyconst.nano
mm = scipyconst.milli
mum = scipyconst.micro
mrad = scipyconst.milli
nrad = scipyconst.nano
pC = scipyconst.pico
muJ = scipyconst.micro
MHz = scipyconst.mega
ps = scipyconst.pico
mmmrad = mm*mrad
deg2rad = scipyconst.degree

# Constants (energies in eV)
fivedeg = (5*scipyconst.pi)/180
me = scipyconst.m_e*(scipyconst.speed_of_light**2)/scipyconst.eV
h = scipyconst.h/scipyconst.eV
re = scipyconst.physical_constants['classical electron radius'][0]
echarge = scipyconst.e


# Input from file
inputdata = []

with open("indata.txt") as csvfile:
    datreader = csv.reader(csvfile, delimiter=' ', quotechar='|')
    
    for row in datreader:
        inputdata.append(row)

dataform = [[1,2,3,4,5,6,7,8,9,10,11,12,13,15,16,17,18],["Ee","wvl","normemitx","normemity","phi","sigL","Q","Epulse","f","tpulse","sigze","DeltaEe","Deltawvl","thetacol","betaIPx","betaIPy","tarrmsBW"], \
["MeV","nm","mmmrad","mmmrad","deg","mum","pC","muJ","MHz","ps","mm","none","none","mrad","m","m","none"],[MeV,nm,mmmrad,mmmrad,deg2rad,mum,pC,muJ,MHz,ps,mm,1,1,mrad,1,1,1]]

for i in range(len(dataform[0])):
            exec("%s=%.12f" % (dataform[1][i],float(inputdata[dataform[0][i]][1])*dataform[3][i]))

print("CBETA Parameters","\nEe:",Ee,"\nwvl:",wvl,"\nnormemitx:",normemitx,"\nnormemity:",normemity,"\nphi:",phi,"\nsigL:",sigL,"\nQ:",Q,"\nEpulse:",Epulse,"\nf:",f,"\ntpulse:",tpulse,"\nsigze:",sigze,"\nDeltaEe:",DeltaEe,"\nDeltawvl:",Deltawvl,"\n")
print("CBETA Variables","\nthetacol:",thetacol,"\nbetaIPx:",betaIPx,"\nbetaIPy:",betaIPy,"\ntarrmsBW:",tarrmsBW,"\n")

# Scattered Photon Energy
def gamma(Ee):
    return (Ee+me)/(me)

def speedfac(Ee):
    return np.sqrt(1-1/gamma(Ee)**2)
    
def EL(wvl):
    return (h*scipyconst.speed_of_light)/wvl

def Egamma(Ee,phi,wvl,theta):
    return (EL(wvl)*(1+speedfac(Ee)*np.cos(phi)))/(1-speedfac(Ee)*np.cos(theta)+(EL(wvl)/Ee)*(1+np.cos(phi+theta)))

# Acceptance Angle
def accang(Ee,thetacol):
    return gamma(Ee)*thetacol    

# Mandelstam Variables
def X(Ee,phi,wvl):
    return (2*gamma(Ee)*EL(wvl)*(1+speedfac(Ee)*np.cos(phi)))/me

def Y(Ee,phi,wvl,theta):
    return (2*gamma(Ee)*Egamma(Ee,phi,wvl,theta)*(1-speedfac(Ee)*np.cos(theta)))/me

# Cross Section (Berestetskii)
def dsigdY(Ee,phi,wvl,theta):
    return (((8*scipyconst.pi*re**2)/(X(Ee,phi,wvl)**2))*(((1/X(Ee,phi,wvl))-(1/Y(Ee,phi,wvl,theta)))**2 \
    +(1/X(Ee,phi,wvl))-(1/Y(Ee,phi,wvl,theta))+(1/4)*((X(Ee,phi,wvl)/Y(Ee,phi,wvl,theta))+(Y(Ee,phi,wvl,theta)/X(Ee,phi,wvl)))))

# Lorentz Invariant Y Differential
def dYdtheta(Ee,theta,phi,wvl):
    return ((X(Ee,phi,wvl)*speedfac(Ee)*np.sin(theta)*(1-speedfac(Ee)*np.cos(theta)+(1+np.cos(phi+theta))*(EL(wvl)/Ee))) \
    - (X(Ee,phi,wvl)*(1-speedfac(Ee)*np.cos(theta))*(speedfac(Ee)*np.sin(theta)-(EL(wvl)/Ee)*np.sin(phi+theta)))) \
    / (1-speedfac(Ee)*np.cos(theta)+(1+np.cos(phi+theta))*(EL(wvl)/Ee))**2
    
# Scattering Angle Cross Section
def dsigdtheta(Ee,phi,wvl,theta):
    return dsigdY(Ee,phi,wvl,theta)*dYdtheta(Ee,theta,phi,wvl)

# Integrate Cross Section
def sig(Ee,phi,wvl,thetacol):
    return scipyint.quad(lambda theta: dsigdtheta(Ee,phi,wvl,theta),0,thetacol)[0]

# Collimated Flux Intermediaries
def Ne(Q):
    return Q/echarge

def NL(Epulse,wvl):
    return Epulse/(echarge*EL(wvl))
    
def sige(betaIP,normemit,Ee):
    return np.sqrt((betaIP*normemit)/gamma(Ee))
    
def sigzl(tpulse):
    return scipyconst.speed_of_light*tpulse
    
def convxy(betaIP,normemit,Ee,sigL):
    return np.sqrt(sige(betaIP,normemit,Ee)**2 + sigL**2)
    
def convz(sigze,tpulse):
    return np.sqrt(sigze**2 + sigzl(tpulse)**2)
    
def zr(sigL,wvl):
    return (4*scipyconst.pi*(sigL**2))/wvl

# FROM HERE ON betax betay differences matter BE CAREFUL
# Head-on Luminosity
def L(Q,Epulse,wvl,betaIPx,normemitx,betaIPy,normemity,Ee,sigL):
    return (Ne(Q)*NL(Epulse,wvl))/(2*scipyconst.pi*(convxy(betaIPx,normemitx,Ee,sigL)*convxy(betaIPy,normemity,Ee,sigL)))

# Angular Crossing + Hourglass Effect
def H(phi,betaIPx,normemitx,betaIPy,normemity,Ee,sigL,sigze,tpulse):
    return np.cos(phi)*np.sqrt((convxy(betaIPx,normemitx,Ee,sigL)**2)*(convxy(betaIPy,normemity,Ee,sigL)**2) \
    / (scipyconst.pi*(convz(sigze,tpulse)**2)))

def Ux(betaIPx,normemitx,Ee,sigL,wvl):
    return np.sqrt((1/2)*((((sige(betaIPx,normemitx,Ee))**2)/(betaIPx**2))+(((sigL)**2)/((zr(sigL,wvl))**2))))

def Uy(betaIPy,normemity,Ee,sigL,wvl):
    return np.sqrt((1/2)*((((sige(betaIPy,normemity,Ee))**2)/(betaIPy**2))+((sigL**2)/((zr(sigL,wvl))**2))))
    
def hmiy(phi,betaIPx,normemitx,Ee,sigL,wvl,Zc):
    return ((np.sin(phi)**2)/((convxy(betaIPx,normemitx,Ee,sigL)**2)+(Ux(betaIPx,normemitx,Ee,sigL,wvl)**2)*(Zc**2)))+((np.cos(phi)**2)/(convz(sigze,tpulse)**2))
    
# RACHG luminosity reduction factors    
def dRACHGdZc(phi,betaIPx,normemitx,betaIPy,normemity,Ee,sigL,sigze,tpulse,wvl,Zc):
    return (H(phi,betaIPx,normemitx,betaIPy,normemity,Ee,sigL,sigze,tpulse)*np.exp(-hmiy(phi,betaIPx,normemitx,Ee,sigL,wvl,Zc)*(Zc**2))) \
    /  (np.sqrt((convxy(betaIPx,normemitx,Ee,sigL)**2)+((Ux(betaIPx,normemitx,Ee,sigL,wvl)**2)*(Zc**2)))*np.sqrt((convxy(betaIPy,normemity,Ee,sigL)**2)+((Uy(betaIPy,normemity,Ee,sigL,wvl)**2)*(Zc**2))))
    
def RACHG(phi,betaIPx,normemitx,betaIPy,normemity,Ee,sigL,sigze,tpulse,wvl):
    return (scipyint.quad(lambda Zc: dRACHGdZc(phi,betaIPx,normemitx,betaIPy,normemity,Ee,sigL,sigze,tpulse,wvl,Zc),-0.1,0.1))[0]    

# Collimated Flux
def Fcol(Ee,phi,wvl,thetacol,betaIPx,normemitx,betaIPy,normemity,sigL,sigze,tpulse,Q,Epulse,f):
    return sig(Ee,phi,wvl,thetacol)*RACHG(phi,betaIPx,normemitx,betaIPy,normemity,Ee,sigL,sigze,tpulse,wvl)*L(Q,Epulse,wvl,betaIPx,normemitx,betaIPy,normemity,Ee,sigL)*f

# Bandwidth Terms
def collterm(Ee,thetacol,phi,wvl):
    return (1/np.sqrt(12))*((accang(Ee,thetacol)**2)/(1+X(Ee,phi,wvl)+((accang(Ee,thetacol)**2)/2)))
    
def beamspreadterm(Ee,phi,wvl,thetacol,DeltaEe):
    return ((2+X(Ee,phi,wvl))/(1+X(Ee,phi,wvl)+accang(Ee,thetacol)**2))*DeltaEe

def laserspreadterm(Ee,thetacol,phi,Deltawvl):
    return ((1+accang(Ee,thetacol)**2)/(1+X(Ee,phi,wvl)+accang(Ee,thetacol)**2))*Deltawvl

def emittanceterm(Ee,betaIPx,normemitx,betaIPy,normemity,phi,wvl):
    return (np.sqrt(2)*gamma(Ee))/(1+X(Ee,phi,wvl))*np.sqrt((normemitx**2/betaIPx**2)+(normemity**2/betaIPy**2))
   
# rms Bandwidth Calc
# non-covariant case used, as dewaling with non-skewed Gaussian Beams
def RMSbandwidth(Ee,thetacol,phi,wvl,DeltaEe,Deltawvl,normemitx,betaIPx,normemity,betaIPy):
    return np.sqrt((collterm(Ee,thetacol,phi,wvl)**2)+(beamspreadterm(Ee,phi,wvl,thetacol,DeltaEe)**2)+(laserspreadterm(Ee,thetacol,phi,Deltawvl)**2) \
    +(emittanceterm(Ee,betaIPx,normemitx,betaIPy,normemity,phi,wvl)**2))

# rms BW penalty function

def penaltyfunc(tarrmsBW,Ee,thetacol,phi,wvl,DeltaEe,Deltawvl,normemitx,betaIPx,normemity,betaIPy):
    return abs(tarrmsBW-RMSbandwidth(Ee,thetacol,phi,wvl,DeltaEe,Deltawvl,normemitx,betaIPx,normemity,betaIPy))
    
# Print out calculations
print("Collimated Flux:",Fcol(Ee,phi,wvl,thetacol,betaIPx,normemitx,betaIPy,normemity,sigL,sigze,tpulse,Q,Epulse,f),"ph/s")
print("Penalty Function:",penaltyfunc(tarrmsBW,Ee,thetacol,phi,wvl,DeltaEe,Deltawvl,normemitx,betaIPx,normemity,betaIPy),"")
toc = time.perf_counter()
print("Calculation Time:",toc-tic,"s")

# Saving Calculations
# Save the optimised settings
print("\nSaving...")
outfile = open("ColFlux_penfunc_Calc.txt","w")
outfile.writelines(["Penalty Function: ",str(penaltyfunc(tarrmsBW,Ee,thetacol,phi,wvl,DeltaEe,Deltawvl,normemitx,betaIPx,normemity,betaIPy))," \n"])
outfile.writelines(["Collimated Flux: ",str(Fcol(Ee,phi,wvl,thetacol,betaIPx,normemitx,betaIPy,normemity,sigL,sigze,tpulse,Q,Epulse,f))," ph/s\n"])
outfile.close()
