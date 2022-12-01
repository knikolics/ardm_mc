from ROOT import TFile,TH2,TCanvas,gROOT,gStyle,TCanvas,TColor
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
h_noCuts = f.Get('fidVol')
#h_Cuts = f.Get('signalEnergy_Cuts')

c = TCanvas('c','',800,640)
c.cd()

h_noCuts.SetStats(0)
h_noCuts.GetYaxis().SetRangeUser(-1050,530)
h_noCuts.SetXTitle('radius^{2} (cm^{2})')
h_noCuts.SetYTitle('z (cm)')
h_noCuts.Draw('colz')
#h_Cuts.Draw('same')

raw_input()
