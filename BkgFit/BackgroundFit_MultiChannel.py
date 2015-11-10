import ROOT as R

import numpy as np
from array import array
import sys
from copy import deepcopy

from GetEigenVariations import GetEigenVariations


regions = {}
h_qcd = {}
h_top = {}
h_top_2b = {}
h_data = {}

useOneTopNuis = None
scaleTop2b = None


def BackgroundFit(datafileName="hist_data.root",
                  topfileName="hist_ttbar.root",
                  distributionName= "LeadCaloJetM",
                  n_trkjet  = ["4"],
                  n_btag    = ["4"],
                  btag_WP     = "77",
                  NRebin = 1,
                  use_one_top_nuis = False,
                  use_scale_top_2b = False,
                  nbtag_top_shape = None,
                  makePlots = False,
                  verbose = True ):
    
    global h_qcd
    global h_top
    global h_top_2b
    global h_data
    global useOneTopNuis
    global scaleTop2b
    global regions

    ################### Parse  #########################
    dist_name   = distributionName
    
    num_trkjet  = np.asarray(n_trkjet)
    if num_trkjet.shape==():
        num_trkjet = np.asarray([n_trkjet])

    num_btag    = np.asarray(n_btag)
    if num_btag.shape==():
        num_btag = np.asarray([n_btag])
    if num_btag.shape!=num_trkjet.shape:
        print "Must have same number of track jet and b-tag regions specified"
        sys.exit(0)
  
    btag_WP     = btag_WP
    
    n_rebin     = NRebin

    topShape_nbtag = nbtag_top_shape
    if nbtag_top_shape == None:
        topShape_nbtag = num_btag

    useOneTopNuis = use_one_top_nuis

    scaleTop2b = use_scale_top_2b
    
    regions = [ num_trkjet[i]+num_btag[i] for i in range(num_trkjet.shape[0]) ]

    datafile = R.TFile(datafileName,"READ")
    topfile  = R.TFile(topfileName,"READ")
    #########################################################

    

    ################### Setup Minuit  ###################
    n_param = len(regions) + ( 1 if useOneTopNuis else len(regions) )
    
    minuit = R.TMinuit(n_param) # 2 parameter fit
    minuit.SetPrintLevel( (1 if verbose else -1) )
    minuit.SetErrorDef(0.5)
    minuit.SetFCN(NegLogL)
    #########################################################


    ################### Get Histograms  ###################
    folder = lambda nt, nb, wp: "GoodEvent_Pass" + nt + "GoodTrackJetPass" + nb + "b" + wp +"PassSBMass/"

    histos = { }

    # collect all histograms
    for r in ["44","43","42","33","32"]:
        folder_r = folder( r[0], r[1], btag_WP)
        data_r   = datafile.Get(folder_r+dist_name).Clone("data_"+r)
        top_r    = topfile.Get(folder_r+dist_name).Clone("top_"+r)
        histo_r  = {"data": data_r, "top": top_r}

        histos[r] = histo_r


    # put relevant histograms in global lists for use in LogLk
    for r in regions:
        hd = histos[r]["data"].Clone("h_data_"+r)
        hq = histos[r[0]+"2"]["data"].Clone("h_qcd_"+r)
        ht2 = histos[r[0]+"2"]["top"].Clone("h_top_2b_"+r)

        if nbtag_top_shape != None:
            ht = histos[r[0]+nbtag_top_shape]["top"].Clone("h_top_"+r)
            ht.Scale( histos[r]["top"].Integral() / ht.Integral() ) #scale to correct norm for region
        else:
            ht = histos[r]["top"].Clone("h_top_"+r)

       
        hq.Add( ht2, -1.0)
        
        h_data[r] = hd 
        h_qcd[r] = hq 
        h_top[r] = ht 
        h_top_2b[r] = ht2

        if NRebin > 1:
            h_data[r].Rebin(NRebin)
            h_qcd[r].Rebin(NRebin)
            h_top[r].Rebin(NRebin)
            h_top_2b[r].Rebin(NRebin)
        
       
    
    #########################################################
  
    results = Fit( minuit )

    evars = GetEigenVariations(results["cov_m"])
    pnom  = np.asarray( results["muqcd"] + results["topscale"] )
    pvars = [ [pnom+evars[i], pnom-evars[i]] for i in range(len(evars)) ]

    results["pnom"]  = pnom
    results["pvars"] = pvars

    print pnom
    print pvars

    #print "Fit Results:"
    #print "mu_qcd = ", results["muqcd"], "+/-", results["muqcd_e"]
    #print "top_scale = ", results["topscale"], "+/-", results["topscale_e"]
    #print "correlation=", results["corr_muqcd_topscale"]

    #ComputeBasicMuQCD( histo_s, histo_c )

    if makePlots:
        for i in range(len(regions)):
            c=MakePlot(regions[i], results["muqcd"][i], results["topscale"][0 if useOneTopNuis else i])


    datafile.Close()
    topfile.Close()
    
    return results







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




def MakePlot(region, muqcd, topscale):
    c=R.TCanvas()
    c.SetFillColor(0)
    c.SetFrameFillColor(0)

    h_data2 = h_data[region].Clone("data2_"+region)
    h_data2.SetFillColor(0)
    h_data2.SetLineColor(R.kBlack)
    h_data2.SetLineWidth(2)
    #h_data2.Rebin(nrebin)
        
    h_top2 = h_top[region].Clone("top2_"+region)
    h_top2.Scale( topscale )
    h_top2.SetFillColor(R.kGreen)
    h_top2.SetLineColor(R.kBlack)
    h_top2.SetLineWidth(1)
    #h_top2.Rebin(nrebin)

    h_qcd2 = h_qcd[region].Clone("qcd2_"+region)
    h_qcd2.Scale( muqcd )
    h_qcd2.SetLineColor(R.kRed)
    h_qcd2.SetFillColor(0)
    h_qcd2.SetLineWidth(1)
    #h_qcd2.Rebin(nrebin)

    h_pred = h_top2.Clone("pred_"+region)
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
    leg.AddEntry(h_data2,"Data ("+region+"), 1.4 fb^{-1}","EL")
    leg.AddEntry(h_top2,"ttbar MC","F")
    leg.AddEntry(h_qcd2,"QCD model","L")
    leg.AddEntry(h_pred,"ttbar MC + QCD model","L")
    leg.SetFillColor(0)
    leg.Draw()

    #raw_input()

    print " "
    print "Region "+region+": Final Numbers after fit:"
    print "Ndata = ", h_data2.Integral()
    print "Npred = ", h_pred.Integral()
    print "Nqcd = ", h_qcd2.Integral()
    print "Ntop = ", h_top2.Integral()
    print " "

    c.SaveAs("fit_"+region+".pdf")

    return c


def Fit( minuit ):
    ClearMinuit( minuit )
    minuit.Migrad()
    minuit.Command("MINOS")

    eparab = R.Double(0) #dummy
    gcc = R.Double(0)    #dummy
    tmp1 = R.Double(0)
    tmp2 = R.Double(0)

    n_reg = len(regions)
    
    muqcd         = [R.Double(0)] * n_reg
    muqcd_e       = [R.Double(0)] * n_reg
    muqcd_e_up    = [R.Double(0)] * n_reg
    muqcd_e_dw    = [R.Double(0)] * n_reg
    topscale      = [R.Double(0)] * (1 if useOneTopNuis else n_reg)
    topscale_e    = [R.Double(0)] * (1 if useOneTopNuis else n_reg)
    topscale_e_up = [R.Double(0)] * (1 if useOneTopNuis else n_reg)
    topscale_e_dw = [R.Double(0)] * (1 if useOneTopNuis else n_reg)


    for i in range(n_reg):
        minuit.GetParameter( i, tmp1, tmp2  )
        muqcd[i]   = deepcopy(tmp1)
        muqcd_e[i] = deepcopy(tmp2)
        
        minuit.mnerrs( i, tmp1, tmp2, eparab, gcc )
        muqcd_e_up[i] = deepcopy(tmp1)
        muqcd_e_dw[i] = deepcopy(tmp2)

        if useOneTopNuis and i!=0:
            continue
    
        minuit.GetParameter( i+n_reg, tmp1, tmp2 )
        topscale[i]   = deepcopy(tmp1)
        topscale_e[i] = deepcopy(tmp2)
        
        minuit.mnerrs( i+n_reg, tmp1, tmp2, eparab, gcc )
        topscale_e_up[i] = deepcopy(tmp1)
        topscale_e_dw[i] = deepcopy(tmp2)
        
    npars = minuit.GetNumPars()
    cov= array('d', [0]* (npars**2) )
    minuit.mnemat(cov, npars) # stored as if [ cov[0,0], cov[1,0], cov[0,1], cov[1,1] ]

    cov_m = np.reshape( cov, (npars, npars) )
    corr_m = np.zeros( (npars, npars) )
    for i in range(npars):
        for j in range(npars):
            corr_m[i,j] = cov_m[i,j] / np.sqrt( cov_m[i,i]*cov_m[j,j] )

    out = { "muqcd" : muqcd,
            "muqcd_e" : muqcd_e,
            "muqcd_e_up" : muqcd_e_up,
            "muqcd_e_dw" : muqcd_e_dw,
            "topscale" : topscale,
            "topscale_e" : topscale_e,
            "topscale_e_up" : topscale_e_up,
            "topscale_e_dw" : topscale_e_dw,
            "cov_m": cov_m,
            "corr_m": corr_m }

    
        
    return out

def ClearMinuit( minuit ):
    minuit.Command("CLEAR")

    for i in range(len(regions)):
        minuit.DefineParameter(i,"muqcd_"+regions[i], 0.01, 0.01, 0.00001, 1)

        if useOneTopNuis and i!=0:
            continue

        topscaleName = "topscale"+("_"+regions[i] if not useOneTopNuis else '')
        minuit.DefineParameter(i+len(regions), topscaleName, 1.3, 0.01, 0.00001, 5)
        
    return
    

def NegLogL(npar, gin, f, par, ifag):
    L = 0.0

    for i in range(len(regions)):
        
        muqcd = par[i]

        if useOneTopNuis:
            topscale = par[len(regions)]
        else:
            topscale = par[i+len(regions)]

        r = regions[i]
        data_r = h_data[r]
        qcd_r = h_qcd[r]
        top_r = h_top[r]
        top2b_r = h_top_2b[r]

        Nbins = data_r.GetNbinsX()

        for ibin in range(1,Nbins+1):
            expected_i = muqcd * qcd_r.GetBinContent(ibin) + topscale * top_r.GetBinContent(ibin)
            
            if scaleTop2b:
                expected_i = expected_i + (muqcd * (1.0 - topscale) * top2b_r.GetBinContent(ibin))
                # use (1.0 - topscale) since top2b is already subtracted from data to make qcd
                # so we need to first add it back, and then subtract the newly scaled top2b
        
            if expected_i > 0:
                L += expected_i - (data_r.GetBinContent(ibin)) * np.log( expected_i );

    f[0] = L
    return



if __name__=="__main__":
	BackgroundFit()
