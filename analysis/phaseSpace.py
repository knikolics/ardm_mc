from ROOT import TFile, TH2F,TEllipse, TLine, gROOT,gStyle,TCanvas,TColor
import numpy as np

gROOT.SetStyle('Plain')
gStyle.SetOptStat(10)
S = np.array((0.00, 0.10, 0.25, 0.45, 0.60, 0.75, 1.00),dtype=np.float64)
R = np.array((1.00, 0.00, 0.00, 0.00, 0.97, 0.97, 0.10),dtype=np.float64)
G = np.array((1.00, 0.97, 0.30, 0.40, 0.97, 0.00, 0.00),dtype=np.float64)
B = np.array((1.00, 0.97, 0.97, 0.00, 0.00, 0.00, 0.00),dtype=np.float64)
TColor.CreateGradientColorTable(7,S,R,G,B,999)
gStyle.SetNumberContours(999)


f = TFile('rootFilesHistos/allSS_histos.root','read')
h = f.Get('phaseSpace')

c = TCanvas('c','',800,640)
h.SetStats(0)
h.SetXTitle('Er (keV)')
h.SetYTitle('radius (mm)')
h.GetXaxis().SetRangeUser(0,250)
h.Draw('colz')

ellipse = TEllipse(250,0,220,450,90,0,-180)
ellipse.SetLineWidth(2)
ellipse.SetLineColor(15)
ellipse.SetFillColor(0)
#ellipse.SetNoEdges()
ellipse.Draw()

Ecut = TLine(30,0,30,500)
Ecut.SetLineWidth(2)
Ecut.SetLineStyle(2)
Ecut.Draw()
Rcut = TLine(0,450,250,450)
Rcut.SetLineWidth(2)
Rcut.SetLineStyle(2)
Rcut.Draw()

raw_input()
