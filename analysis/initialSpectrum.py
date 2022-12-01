from ROOT import TH1F, TCanvas, TLegend

an_232Th = []
for line in reversed(open("../data/SS_232Th.dat").readlines()):
    entry = (line.rstrip()).split()         
    #energy.append(float(entry[1]))
    an_232Th.append(float(entry[2]))

an_238U = []
sf_238U = []
for line in reversed(open("../data/SS_238U.dat").readlines()):
    entry = (line.rstrip()).split()
    an_238U.append(float(entry[2]))
    sf_238U.append(float(entry[3]))

h_232Th = TH1F('h_232Th','',100,0,10)
h_238U_an = TH1F('h_238U_an','',100,0,10)
h_238U_sf = TH1F('h_238U_sf','',100,0,10)

for iBin in range(100):
   h_232Th.SetBinContent(iBin+1,an_232Th[iBin])
   h_238U_an.SetBinContent(iBin+1,an_238U[iBin])
   h_238U_sf.SetBinContent(iBin+1,sf_238U[iBin])


c = TCanvas('c','',800,640)
h_238U_sf.SetTitle('Stainless Steel')
h_238U_sf.SetXTitle('E (MeV)')
h_238U_sf.SetYTitle('Flux (n/s/cm^{3}/MeV)')
h_238U_sf.SetStats(0)
h_238U_sf.SetLineColor(4)
h_238U_sf.SetLineWidth(2)
h_238U_sf.Draw()
h_232Th.SetStats(0)
h_232Th.SetLineColor(2)
h_232Th.SetLineWidth(2)
h_232Th.Draw('same')
h_238U_an.SetStats(0)
h_238U_an.SetLineColor(3)
h_238U_an.SetLineWidth(2)
h_238U_an.Draw('same')

leg = TLegend(0.7,0.72,0.9,0.9)
leg.SetFillColor(0)
leg.SetBorderSize(0)
leg.AddEntry(h_232Th,'232Th (#alpha,n)','l')
leg.AddEntry(h_238U_an,'238U (#alpha,n)','l')
leg.AddEntry(h_238U_sf,'238U spont. fission','l')
leg.Draw()

raw_input()
