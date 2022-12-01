from ROOT import TFile, TTree, TH1D, TH2D, TCanvas, TLegend, TColor, gROOT, gStyle, TLine, TLatex, TPaveText
import sys
import math
import numpy as np

gROOT.SetStyle('Plain')
gStyle.SetOptStat(10)
S = np.array((0.00, 0.10, 0.25, 0.45, 0.60, 0.75, 1.00),dtype=np.float64)
R = np.array((1.00, 0.00, 0.00, 0.00, 0.97, 0.97, 0.10),dtype=np.float64)
G = np.array((1.00, 0.97, 0.30, 0.40, 0.97, 0.00, 0.00),dtype=np.float64)
B = np.array((1.00, 0.97, 0.97, 0.00, 0.00, 0.00, 0.00),dtype=np.float64)
TColor.CreateGradientColorTable(7,S,R,G,B,999)
gStyle.SetNumberContours(999)


f = TFile(sys.argv[1])

energy =  f.Get('singleErecoil')
#position = f.Get('')
scatterN_0 = f.Get('NumScatter_moreThan0')
scatterN_1 = f.Get('NumScatter_moreThan1')
zEnergy = f.Get('zFid')
rEnergy = f.Get('radialFid')

moreThan0 = scatterN_0.GetEntries()
moreThan1 = scatterN_1.GetEntries()
scatter = round(100*moreThan1/moreThan0,2)
print 'Multiple Scatter vs. single scatter', round(100*moreThan1/moreThan0,2),'% \n' 

c1 = TCanvas('c1','',800,600)
c1.cd()
energy.SetStats(0)
energy.SetTitle('E_{recoil} > 10 keV')
energy.SetXTitle('E_{recoil} (keV)')
energy.GetXaxis().SetTitleOffset(1.28)
energy.SetYTitle('Entries/(keV)')
energy.SetLineColor(2)
energy.SetLineWidth(2)
energy.Draw()
c1.Update()

c2 = TCanvas('c2','',800,600)
c2.cd()
myText = TLatex()
myText.SetNDC()
scatterN_0.SetStats(0)
scatterN_0.SetTitle('Number of elastic scatterings > 2 keV')
scatterN_0.SetXTitle('# scatterings')
scatterN_0.SetLineWidth(2)
scatterN_0.Draw()
#myText.AddText('%g percent scatter more than once'%(scatter))
myText.DrawLatex(0.4,0.85,'%g percent scatter more than once'%(scatter))
c2.Update()

c3 = TCanvas('c3','',800,600)
c3.cd()
zEnergy.SetStats(0)
zEnergy.SetXTitle('E_{recoil} (keV)')
zEnergy.GetXaxis().SetTitleOffset(1.2)
zEnergy.SetYTitle('z position (mm)')
zEnergy.GetYaxis().SetTitleOffset(1.3)
zEnergy.Draw('colz')
c3.Update()

c4 = TCanvas('c4','',800,600)
c4.cd()
rEnergy.SetStats(0)
rEnergy.SetXTitle('E_{recoil} (keV)')
rEnergy.GetXaxis().SetTitleOffset(1.2)
rEnergy.SetYTitle('radius (mm)')
rEnergy.GetYaxis().SetTitleOffset(1.3)
rEnergy.Draw('colz')
c4.Update()

raw_input()
