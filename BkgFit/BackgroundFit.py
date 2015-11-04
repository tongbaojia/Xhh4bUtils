import ROOT as R

import numpy as np
from array import array

h_qcd = None
h_top = None
h_data = None


x=5

def BackgroundFit():
    global h_qcd
    global h_top
    global h_data

    ################### Setup Minuit  ###################
    minuit = R.TMinuit(2) # 2 parameter fit
    minuit.SetPrintLevel(1)
    minuit.SetErrorDef(0.5)
    minuit.SetFCN(NegLogL)
    #########################################################

    ################### Distribution Info  ###################
    dist_name   = "LeadCaloJetM"
    num_trkjet  = "4"
    num_btag    = "3"
    btag_WP     = "80"

    n_rebin     = 1
    #########################################################


    ################### Get Histograms  ###################
    channel_s = num_trkjet+"GoodTrackJetPass"+num_btag+"b"+btag_WP
    channel_c  = num_trkjet+"GoodTrackJetPass2b"+btag_WP

    folder_s  = "GoodEvent_Pass"+channel_s+"PassSBMass/"
    folder_c   = "GoodEvent_Pass"+channel_c+"PassSBMass/"

    datafile = R.TFile("hist_data.root","READ")
    topfile  = R.TFile("hist_ttbar.root","READ")

    histo_s = { "data" : datafile.Get(folder_s+dist_name).Clone("data_s"),
                "top" : topfile.Get(folder_s+dist_name).Clone("top_s")   }

    histo_c = { "data" : datafile.Get(folder_c+dist_name).Clone("data_c"),
                "top" : topfile.Get(folder_c+dist_name).Clone("top_c")   }

    h_data = histo_s["data"].Clone("data")
    h_top = histo_s["top"].Clone("top")
    h_qcd = histo_c["data"].Clone("qcd")
    h_qcd.Add( histo_c["top"], -1.0)

    
    #########################################################
  
    results = Fit( minuit )

    print "Fit Results:"
    print "mu_qcd = ", results["muqcd"], "+/-", results["muqcd_e"]
    print "top_scale = ", results["topscale"], "+/-", results["topscale_e"]
    print "correlation=", results["corr_muqcd_topscale"]

    ComputeBasicMuQCD( histo_s, histo_c )

    c=MakePlot( results["muqcd"], results["topscale"])

    return







def ComputeBasicMuQCD( histo_s, histo_c ):
    #####
    N_s_data = histo_s["data"].Integral()
    E_s_data = N_s_data**(0.5)

    E_s_top = R.Double(0)
    N_s_top = histo_s["top"].IntegralAndError(0, histo_s["top"].GetNbinsX()+1, E_s_top)

    N_s_qcd = N_s_data - N_s_top
    E_s_qcd = ( E_s_data**2 + E_s_top**2 )**(0.5)

    ####
    N_c_data = histo_c["data"].Integral()
    E_c_data = N_c_data**(0.5)

    E_c_top = R.Double(0)
    N_c_top = histo_c["top"].IntegralAndError(0, histo_c["top"].GetNbinsX()+1, E_c_top)

    N_c_qcd = N_c_data - N_c_top
    E_c_qcd = ( E_c_data**2 + E_c_top**2 )**(0.5)

    ####
    mu_qcd = N_s_qcd / N_c_qcd
    mu_qcd_err = mu_qcd * ( (E_s_qcd/N_s_qcd)**2 + (E_c_qcd/N_c_qcd) )**(0.5)

    print "No Fit: mu_qcd = ", mu_qcd, "+/-", mu_qcd_err
    print "No Fit: Ndata_s=", N_s_data, "+/-", E_s_data
    print "No Fit: Ndata_c=", N_c_data, "+/-", E_c_data
    print "No Fit: Ntop_s=", N_s_top, "+/-", E_s_top
    print "No Fit: Ntop_c=", N_c_top, "+/-", E_c_top
    print "No Fit: Nqcd_s=", N_s_qcd, "+/-", E_s_qcd
    print "No Fit: Nqcd_c=", N_c_qcd, "+/-", E_c_qcd
    
    return [mu_qcd, mu_qcd_err]




def MakePlot( muqcd, topscale):
    c=R.TCanvas()
    c.SetFillColor(0)
    c.SetFrameFillColor(0)

    h_data2 = h_data.Clone("data2")
    h_data2.SetFillColor(0)
    h_data2.SetLineColor(R.kBlack)
    h_data2.SetLineWidth(2)
    #h_data2.Rebin(nrebin)
        
    h_top2 = h_top.Clone("top2")
    h_top2.Scale( topscale )
    h_top2.SetFillColor(R.kGreen)
    h_top2.SetLineColor(R.kBlack)
    h_top2.SetLineWidth(1)
    #h_top2.Rebin(nrebin)

    h_qcd2 = h_qcd.Clone("qcd2")
    h_qcd2.Scale( muqcd )
    h_qcd2.SetLineColor(R.kRed)
    h_qcd2.SetFillColor(0)
    h_qcd2.SetLineWidth(1)
    #h_qcd2.Rebin(nrebin)

    h_pred = h_top2.Clone("pred")
    h_pred.Add( h_qcd2, 1.0)
    h_pred.SetLineColor(R.kBlue)
    h_pred.SetFillColor(0)
    h_pred.SetLineWidth(1)

    h_data2.Draw("E")
    h_top2.Draw("sameHIST")
    h_qcd2.Draw("sameHIST")
    h_pred.Draw("sameHIST")
    h_data2.Draw("sameE")

    leg = R.TLegend(0.1,0.7,0.48,0.9)
    leg.AddEntry(h_data2,"Data, 1.4 fb^{-1}","EL")
    leg.AddEntry(h_top2,"ttbar MC","F")
    leg.AddEntry(h_qcd2,"QCD model","L")
    leg.AddEntry(h_pred,"ttbar MC + QCD model","L")
    leg.SetFillColor(0)
    leg.Draw()

    #raw_input()

    print "Final Numbers after fit:"
    print "Ndata = ", h_data2.Integral()
    print "Npred = ", h_pred.Integral()
    print "Nqcd = ", h_qcd2.Integral()
    print "Ntop = ", h_top2.Integral()
    

    c.SaveAs("c.pdf")

    return c


def Fit( minuit ):
    ClearMinuit( minuit )
    minuit.Migrad()
    minuit.Command("MINOS")

    eparab = R.Double(0) #dummy
    gcc = R.Double(0)    #dummy
    
    muqcd         = R.Double(0)
    muqcd_e       = R.Double(0)
    muqcd_e_up    = R.Double(0)
    muqcd_e_dw    = R.Double(0)
    topscale      = R.Double(0)
    topscale_e    = R.Double(0)
    topscale_e_up = R.Double(0)
    topscale_e_dw = R.Double(0)

  
    minuit.GetParameter( 0, muqcd, muqcd_e )
    minuit.mnerrs( 0, muqcd_e_up, muqcd_e_dw, eparab, gcc )

    minuit.GetParameter( 1, topscale, topscale_e );
    minuit.mnerrs( 1, topscale_e_up, topscale_e_dw, eparab, gcc );

    cov= array('d', [0,0,0,0])
    minuit.mnemat(cov,2) # stored as if [ cov[0,0], cov[1,0], cov[0,1], cov[1,1] ]
    corr_muqcd_topscale = cov[1]  / ( cov[0] * cov[3] )**(0.5)

    out = { "muqcd" : muqcd,
            "muqcd_e" : muqcd_e,
            "muqcd_e_up" : muqcd_e_up,
            "muqcd_e_dw" : muqcd_e_dw,
            "topscale" : topscale,
            "topscale_e" : topscale_e,
            "topscale_e_up" : topscale_e_up,
            "topscale_e_dw" : topscale_e_dw,
            "corr_muqcd_topscale": corr_muqcd_topscale }

    
        
    return out

def ClearMinuit( minuit ):
    minuit.Command("CLEAR");
    minuit.DefineParameter(0,"muqcd", 0.01, 0.01, 0.00001, 1);
    minuit.DefineParameter(1,"topscale",1.5, 0.1, 0.00001, 5);
    return
    

def NegLogL(npar, gin, f, par, ifag):
    L = 0.0

    muqcd = par[0];
    topscale = par[1];

    Nbins = h_data.GetNbinsX()

    for ibin in range(1,Nbins+1):
        expected_i = muqcd * h_qcd.GetBinContent(ibin) + topscale * h_top.GetBinContent(ibin)
        
        if expected_i > 0:
            L += expected_i - (h_data.GetBinContent(ibin)) * np.log( expected_i );

    f[0] = L
    return



if __name__=="__main__":
	BackgroundFit()
