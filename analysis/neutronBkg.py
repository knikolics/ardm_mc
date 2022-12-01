from ROOT import TFile, TH1, TCanvas, TLegend, gROOT
import math

gROOT.SetStyle('Plain')

hList1 = []
hList2 = []
hList3 = []

#tankVol = 1.50896e5 # cm^3

totalFlux_U238_sf = [1.09e-10,1.2805e-11,3.03767e-11,2.45196e-11,1.09e-10,1.2805e-11,5.38088e-11,2.04337e-12,1.36903e-11] # SS, HiPoly, PMT glass, PMT base, PMT electrode, PMT polyethylene, resistors, Perlite, Borotron
totalFlux_U238_sf_err = [1.83039e-11,2.15069e-12,5.10192e-12,4.11818e-12,1.83039e-11,2.15069e-12,9.03746e-12,3.43196e-13,2.29935e-12]

totalFlux_U238_an = [3.864e-11,1.90029e-11,1.45416e-10,5.36462e-11,3.864e-11,1.90029e-11,3.80755e-10,3.72831e-12,9.47187e-11] # SS, HiPoly, PMT glass, resistors, Perlite, Borotron
totalFlux_U238_an_err = [7.74001e-12,2.37842e-12,2.14661e-11,6.7845e-12,7.74001e-12,2.37842e-12,6.79417e-11,6.21923e-13,1.35946e-11]

totalFlux_Th232_an = [5.167e-11,8.11115e-12,6.28144e-11,2.46847e-11,5.167e-11,8.11115e-12,2.12966e-10,2.03839e-12,3.60525e-11] # SS, HiPoly, PMT glass, resistors, Perlite, Borotron
totalFlux_Th232_an_err = [9.33026e-12,9.65659e-13,8.91172e-12,3.19111e-12,9.33026e-12,9.65659e-13,3.64278e-11,3.35011e-13,4.86334e-12]

# contamination [Bq/kg]
material = ['SS','HiPoly','PMT glass','PMT base','PMT electrode','PMT polyethylene','resistors','Perlite','Borotron']
activityU238_combi = [0.32412,0.14674,0.67934,11.19209,1.9743,2.9976,1.85843,52.80022,0.69469] # conservative values: Pa-234m values
activityU238_low = [0.006214,0.0077,0.6432,10.122,0.1951,0.0362,1.3524,46.3819,0.01502] # average without Pa & Ra contribution
activityU238_average = [0.0472,0.0251,0.6793,11.1921,0.4175,0.0625,1.8584,52.8002,0.09996] # average with Pa & Ra
activityU238_high = [0.32412,0.14674,0.72517,18.683,1.9743,2.9976,4.1219,71.445,0.69469] # only Pa-234m values
activityU238_err = [0.0006,0.0015,0.0429,0.4808,0.0381,0.0118,0.1000,2.685,0.0023]

activityTh232 = [0.017844,0.01284,0.12920,13.246,0.11722,0.15482,1.91256,54.58463,0.0370965] # average values
activityTh232_err = [0.0002,0.0148,0.0297,0.9968,0.0408,0.0668,0.1782,4.629,0.0145]

# grams of element producing 1 Bq
# U238_g = 8.040833e-05
# Th232_g = 2.464642e-04
# --> 1 ppb: SF = 1/1e-09
U238ppb = 80.40833
Th232ppb = 246.4642

contaminU238_low = []
contaminU238_average = []
contaminU238_high = []
contaminU238_combi = []
contaminTh232 = []
contaminU238 = []

contaminU238_err = []
contaminTh232_err = []

for iAct in range(len(activityU238_low)):
    contaminU238_low.append(activityU238_low[iAct]*U238ppb)
    contaminU238_average.append(activityU238_average[iAct]*U238ppb)
    contaminU238_high.append(activityU238_high[iAct]*U238ppb)
    contaminU238_combi.append(activityU238_combi[iAct]*U238ppb)
    contaminTh232.append(activityTh232[iAct]*Th232ppb)
    contaminU238_err.append(activityU238_err[iAct]*U238ppb)
    contaminTh232_err.append(activityTh232_err[iAct]*Th232ppb)


whichContLevel = raw_input('low, average, high or combi?\n')
if whichContLevel=='low':
    contaminU238 = contaminU238_low
elif whichContLevel=='average':
    contaminU238 = contaminU238_average
elif whichContLevel=='high':
    contaminU238 = contaminU238_high
elif whichContLevel=='combi':
    contaminU238 = contaminU238_combi
else:
    print 'so what??\n'


U238_an_err = []
U238_sf_err = []
Th232_an_err = []

#compute RELATIVE error
for iErr in range(len(contaminU238_err)):
    #print 'contaminU238_err ',contaminU238_err[iAct],'\n'
    err1 = math.sqrt(pow((contaminU238_err[iErr]/contaminU238[iErr]),2)+pow((totalFlux_U238_an_err[iErr]/totalFlux_U238_an[iErr]),2))
    U238_an_err.append(err1)
    err2 = math.sqrt(pow((contaminU238_err[iErr]/contaminU238[iErr]),2)+pow((totalFlux_U238_sf_err[iErr]/totalFlux_U238_sf[iErr]),2))
    U238_sf_err.append(err2)
    err3 = math.sqrt(pow((contaminTh232_err[iErr]/contaminTh232[iErr]),2)+pow((totalFlux_Th232_an_err[iErr]/totalFlux_Th232_an[iErr]),2))
    Th232_an_err.append(err3)


infile = raw_input('choose between SS, HiPoly, PMTglass, PMTbase, PMTelectrode, PMTpoly, resistors, Perlite or Borotron\n')

if infile=='SS':
    iVol = 0
    #name = ['histos_SS_232Th_alphaN.root','histos_SS_238U_alphaN.root','histos_SS_238U_spontaneousFission.root']
    legend = 'Stainless steel'
elif infile=='HiPoly':
    iVol = 1
    #name = ['histos_Poly_232Th_an.root','histos_Poly_238U_an.root','histos_Poly_238U_sf.root']
    legend = 'HiPoly neutron shield'
elif infile=='PMTglass':
    iVol = 2
    legend = 'PMT glass'
elif infile=='PMTbase':
    iVol = 3
    legend = 'PMT base'
elif infile=='PMTelectrode':
    iVol = 4
    legend = 'PMT electrode'
elif infile=='PMTpoly':
    iVol = 5
    legend = 'PMT polyethylene'
elif infile=='resistors':
    iVol = 6
    legend = 'HV resistors'
elif infile=='Perlite':
    iVol = 7
    legend = 'Perlite'
elif infile=='Borotron':
    iVol = 8
    legend = 'Borotron'
#elif infile=='Ar':
#name = ['histos_Ar_232Th_an.root','histos_Ar_238U_an.root','histos_Ar_238U_sf.root']
#legend = 'Ar'
else:
   print 'hey, i said: choose the volume!\n'

#lowEthreshold = raw_input('lower energy threshold (keV)?\n')
#highEthreshold = raw_input('upper energy threshold (keV)?\n')
exposure = raw_input('how many days exposure?\n')
"""
# fraction w/ Borotron neutron shield - S1: 30 < E < 100 keV, S2: 10keV
fractionSS_U238_an = [0.02619,0.008976,0.005904]
fractionSS_U238_sf = [0.025192,0.008954,0.005898]
fractionSS_Th232_an = [0.025584,0.009378,0.006084]

fractionsTopFl_U238_an = [0.057002,0.050662,0.049536]
fractionsTopFl_U238_sf = [0.053182,0.04776,0.046788]
fractionsTopFl_Th232_an = [0.052834,0.047334,0.046414]

fractionsPoly_U238_an = [0,0.027466,0.015918] # Borotron!
fractionsPoly_U238_sf = [0,0.020528,0.010984] # Borotron!
fractionsPoly_Th232_an = [0,0.026944,0.016028] # Borotron!

fractionsPMTGlass_U238_an = [0.049748,0.049454,0.049214]
fractionsPMTGlass_U238_sf = [0.05774,0.057416,0.05798]
fractionsPMTGlass_Th232_an = [0.050754,0.050062,0.050516]

fractionsHVr_U238_an = [0.08319,0.080352,0.080198]
fractionsHVr_U238_sf = [0.084936,0.081556,0.081474]
fractionsHVr_Th232_an = [0.08433,0.0818,0.081202]

fractionsPMTbase_U238_an = [0.0557204,0.0581719,0.058894]
fractionsPMTbase_U238_sf = [0.0631997,0.0659863,0.0660524]
fractionsPMTbase_Th232_an = [0.0561635,0.0595445,0.0598551]

fractionsPMTelectrode_U238_an = [0.053792,0.052578,0.052516]
fractionsPMTelectrode_U238_sf = [0.049474,0.049226,0.04859]
fractionsPMTelectrode_Th232_an = [0.049378,0.048534,0.048144]

fractionsPerlite_U238_an = [0.015388,0.012082,0.011496]
fractionsPerlite_U238_sf = [0.015256,0.011772,0.010968]
fractionsPerlite_Th232_an = [0.015414,0.01209,0.011544]

#------------------------------------------------------------

#fractions Borotron neutron shield - S1: 30 < E < 100 keV, S2: 1keV
fractionSS_238U_an = [0.018714,0.006442,0.00436]
fractionSS_238U_sf = [0.018166,0.006574,0.004316]
fractionSS_232Th_an = [0.018456,0.006698,0.004552]

fractionTopFlange_238U_an = [0.042624,0.037802,0.036922]
fractionTopFlange_238U_sf = [0.0391,0.035474,0.034482]
fractionTopFlange_232Th_an = [0.039136,0.035236,0.034398]

fractionPoly_238U_an = [0,0.018782,0.010888]
fractionPoly_238U_sf = [0,0.013918,0.007532]
fractionPoly_232Th_an = [0,0.018482,0.010972]

fractionHVr_238U_an = [0.057928,0.056842,0.055614]
fractionHVr_238U_sf = [0.05912,0.058048,0.056986]
fractionHVr_232Th_an = [0.058878,0.057164,0.057594]

fractionPMTGlass_238U_an = [0.035822,0.035564,0.035428]
fractionPMTGlass_238U_sf = [0.042516,0.042278,0.042252]
fractionPMTGlass_232Th_an = [0.037078,0.036324,0.036374]

fractionPerlite_238U_an = [0.00987,0.00787,0.007998]
fractionPerlite_238U_sf = [0.010076,0.007872,0.007448]
fractionPerlite_232Th_an = [0.01019,0.008134,0.007724]

fractionPMTbase_238U_an = [0.038212,0.040772,0.040706]
fractionPMTbase_238U_sf = [0.043432,0.04537,0.046022]
fractionPMTbase_232Th_an = [0.038508,0.040986,0.04092]

fractionPMTelectrode_238U_an = [0.053792,0.052578,0.052516] #unchanged!
fractionPMTelectrode_238U_sf = [0.049474,0.049226,0.04859] #unchanged!
fractionPMTelectrode_232Th_an = [0.049378,0.048534,0.048144] #unchanged!
#-----------------------------------------------------------------

#fractions Borotron neutron shield - S1: E - [50,100], S2: 2keV
fractionSS_238U_an = [0.010808,0.00358,0.002442]
fractionSS_238U_sf = [0.010436,0.003616,0.002438]
fractionSS_232Th_an = [0.010612,0.0038,0.002582]

fractionTopFlange_238U_an = [0.026186,0.023316,0.022602]
fractionTopFlange_238U_sf = [0.023276,0.020922,0.020472]
fractionTopFlange_232Th_an = [0.02342,0.02111,0.020642]

fractionPoly_238U_an = [0,0.010876,0.006448]
fractionPoly_238U_sf = [0,0.00759,0.00425]
fractionPoly_232Th_an = [0,0.01059,0.00637]

fractionHVr_238U_an = [0.034702,0.03417,0.033604]
fractionHVr_238U_sf = [0.035552,0.034832,0.033982]
fractionHVr_232Th_an = [0.034868,0.033566,0.03363]

fractionPMTGlass_238U_an = [0.021412,0.021476,0.021374]
fractionPMTGlass_238U_sf = [0.025714,0.02562,0.025566]
fractionPMTGlass_232Th_an = [0.022658,0.022224,0.021828]

fractionPMTbase_238U_an = [0.022626,0.024478,0.024464]
fractionPMTbase_238U_sf = [0.026624,0.027434,0.027788]
fractionPMTbase_232Th_an = [0.023042,0.02426,0.024216]

fractionPerlite_238U_an = [0.005492,0.004444,0.004448]
fractionPerlite_238U_sf = [0.00562,0.004308,0.004138]
fractionPerlite_232Th_an = [0.005746,0.004642,0.004382]

fractionPMTelectrode_238U_an = [0.053792,0.052578,0.052516] #unchanged!
fractionPMTelectrode_238U_sf = [0.049474,0.049226,0.04859] #unchanged!
fractionPMTelectrode_232Th_an = [0.049378,0.048534,0.048144] #unchanged!
"""
#---------------------------------------------------------------------
# fractions Borotron neutron shield - S1: E - [50,100], S2: 2keV, FV cuts: r = 345mm && -500mm < z < 450mm
fractionSS_238U_an = [0.00893,0.00304,0.00193]
fractionSS_238U_sf = [0.00876,0.00316,0.00219]
fractionSS_232Th_an = [0.00818,0.00338,0.00218]

fractionTopFlange_238U_an = [0.00566,0.00491,0.00443]
fractionTopFlange_238U_sf = [0.00554,0.00454,0.00422]
fractionTopFlange_232Th_an = [0.00542,0.0045,0.00445]

fractionPoly_238U_an = [0,0.01079,0.0065]
fractionPoly_238U_sf = [0,0.00941,0.00559]
fractionPoly_232Th_an = [0,0.0106,0.00712]

fractionHVr_238U_an = [0.03321,0.03253,0.03215]
fractionHVr_238U_sf = [0.03428,0.03229,0.03237]
fractionHVr_232Th_an = [0.03271,0.03152,0.03171]

fractionPMTGlass_238U_an = [0.01059,0.00772,0.00859]
fractionPMTGlass_238U_sf = [0.0108,0.00877,0.00965]
fractionPMTGlass_232Th_an = [0.0114,0.00909,0.0087]

fractionPMTbase_238U_an = [0.00698,0.00581,0.00572]
fractionPMTbase_238U_sf = [0.00709,0.00637,0.00572]
fractionPMTbase_232Th_an = [0.00656,0.00592,0.00584]

fractionPerlite_238U_an = [0,0,0]
fractionPerlite_238U_sf = [0,0,0]
fractionPerlite_232Th_an = [0,0,0]

fractionPMTelectrode_238U_an = [0.053792,0.052578,0.052516] #unchanged!
fractionPMTelectrode_238U_sf = [0.049474,0.049226,0.04859] #unchanged!
fractionPMTelectrode_232Th_an = [0.049378,0.048534,0.048144] #unchanged!


fraction_U238_an = [0,fractionSS_238U_an,fractionTopFlange_238U_an,fractionPoly_238U_an,fractionPoly_238U_an,fractionHVr_238U_an,fractionPMTGlass_238U_an,fractionPMTbase_238U_an,fractionPMTelectrode_238U_an,fractionPerlite_238U_an]
fraction_U238_sf = [0,fractionSS_238U_sf,fractionTopFlange_238U_sf,fractionPoly_238U_sf,fractionPoly_238U_an,fractionHVr_238U_sf,fractionPMTGlass_238U_sf,fractionPMTbase_238U_sf,fractionPMTelectrode_238U_sf,fractionPerlite_238U_sf]
fraction_Th232_an = [0,fractionSS_232Th_an,fractionTopFlange_232Th_an,fractionPoly_232Th_an,fractionPoly_232Th_an,fractionHVr_232Th_an,fractionPMTGlass_232Th_an,fractionPMTbase_232Th_an,fractionPMTelectrode_232Th_an,fractionPerlite_232Th_an]

thickness = [0,5,10]
volumes = [0,144501,37946.9,160070.273,338846.687,113.7532,9148.008,217.14696,591,776958] # Dewar, Top Flange, neutron shield 5cm, neutron shield 10cm, resistors, PMT Glass, PMT base

if infile=='SS':
    vol = raw_input('Dewar (1) or Top Flange (2)?\n')
elif infile=='HiPoly' or infile=='Borotron':
    vol = 3
elif infile=='resistors':
    vol = 5
elif infile=='PMTglass':
    vol = 6
elif infile=='PMTbase':
    vol = 7
elif infile=='PMTelectrode':
    vol = 8
elif infile=='Perlite':
    vol = 9
else:
    vol = 0

whichVol = int(vol)

#neutronShieldVol = [0,160070.273,338846.687]

for iThickness in range(len(thickness)):
    totFluxIntegral_U238_an = totalFlux_U238_an[iVol] * volumes[whichVol] * contaminU238[iVol] * (fraction_U238_an[whichVol][iThickness]) * 24*3600
    totFluxIntegral_U238_sf = totalFlux_U238_sf[iVol] * volumes[whichVol] * contaminU238[iVol] * (fraction_U238_sf[whichVol][iThickness]) * 24*3600
    totFluxIntegral_Th232_an = totalFlux_Th232_an[iVol] * volumes[whichVol] * contaminTh232[iVol] * (fraction_Th232_an[whichVol][iThickness]) * 24*3600
    totErr = math.sqrt(pow(U238_an_err[iVol],2)+pow(U238_sf_err[iVol],2)+pow(Th232_an_err[iVol],2))
    print legend,', neutron Shield ', thickness[iThickness], 'cm: neutrons per day U238(a,n) ', round(totFluxIntegral_U238_an,3),' +/- ',round(totFluxIntegral_U238_an*U238_an_err[iVol],3),'\n'
    print legend,', neutron Shield ', thickness[iThickness], 'cm: neutrons per day U238 sf ', round(totFluxIntegral_U238_sf,3),' +/- ',round(totFluxIntegral_U238_sf*U238_sf_err[iVol],3),'\n'
    print legend,', neutron Shield ', thickness[iThickness], 'cm: neutrons per day Th232(a,n) ', round(totFluxIntegral_Th232_an,3),' +/- ',round(totFluxIntegral_Th232_an*Th232_an_err[iVol],3),'\n'
    totalNeutrons = (totFluxIntegral_U238_an + totFluxIntegral_U238_sf + totFluxIntegral_Th232_an)*int(exposure)
    print 'total number of neutrons for ',thickness[iThickness], 'cm neutron shield, ' ,exposure,' days of exposure = ',round(totalNeutrons,2),' +/- ',round(totalNeutrons*totErr,2),'\n'

raw_input()
