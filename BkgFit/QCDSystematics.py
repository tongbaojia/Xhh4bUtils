import ROOT as R

import numpy as np
from array import array
import sys
from copy import deepcopy

from GetEigenVariations import GetEigenVariations

from HistoTools import HistLocationString as HistLocStr

def QCDSystematics(datafileName="hist_data.root",
                    topfileName="hist_ttbar.root",
                    distributionName= "LeadCaloJetM",
                    n_trkjet  = ["4","4"],
                    n_btag    = ["4","3"],
                    btag_WP     = "77",
                    mu_qcd_vals = [1.0, 1.0],
                    topscale_vals = [1.0, 1.0],
                    NRebin = 1,
                    use_one_top_nuis = False,
                    use_scale_top_2b = False,
                    makePlots = False,
                    verbose = False,
                    outfileNameBase="QCDSysfit.root"):

    ##### Parse Inputs ############################################
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

    useOneTopNuis = use_one_top_nuis

    scaleTop2b = use_scale_top_2b

    n_channels = num_trkjet.shape[0]

    regions = [ num_trkjet[i]+num_btag[i] for i in range(n_channels) ]
    ##################################################################


    ##### Get Signal Region Histograms ################################
    datafile = R.TFile(datafileName,"READ")
    topfile  = R.TFile(topfileName,"READ")


    histos = {}
    
    # collect all histograms
    for r in ["44","43","42","33","32"]:
        folder_r = HistLocStr(dist_name, r[0], r[1], btag_WP, "CR")  #folder( r[0], r[1], btag_WP)
        
        data_r   = datafile.Get(folder_r).Clone("data_"+r)
        data_r.SetDirectory(0)
        
        top_r    = topfile.Get(folder_r).Clone("top_"+r)
        top_r.SetDirectory(0)

        for ibin in range(1, top_r.GetNbinsX()+1):
            if top_r.GetBinContent(ibin) < 0:
                top_r.SetBinContent(ibin, 0)
                top_r.SetBinError(ibin, 0)
                

        histos[r]     = {"data": data_r,            "top": top_r}

    datafile.Close()
    topfile.Close()
    ##################################################################

    ####### outpue object ###################
    QCDSyst_Dict = {}
    

    ##### scaling and subtractions #################################
    for ir in range(len(regions)):
        r = regions[ir]
                
        r_2b = r[0]+"2"
        r_3b = r[0]+"3"

        top_2b = histos[r_2b]["top"].Clone("top_2b__"+r)
        if scaleTop2b:
            top_2b.Scale( (topscale_vals[0] if use_one_top_nuis else topscale_vals[ir]) )

        qcd_r = histos[r_2b]["data"].Clone("qcd__"+r)
        qcd_r.Add( top_2b, -1)      # added by Qi --- we still want top to be subtracted, given that their fraction is increasing in Run 2.
        qcd_int = qcd_r.Integral()


        top_r = histos[r]["top"].Clone("top__"+r)
        top_int = top_r.Integral()


        mu_qcd = mu_qcd_vals[ir]
        top_scale = (topscale_vals[0] if use_one_top_nuis else topscale_vals[ir])
        
        qcd_r.Scale( mu_qcd )
        top_r.Scale( top_scale )


        #now do ratio
        bkg_r = qcd_r.Clone("bkg__"+r)
        bkg_r.Add( top_r )

        ## c=R.TCanvas()
        ## bkg_r_c = bkg_r.Clone("bkg_clone__"+r)
        ## bkg_r_c.SetDirectory(0)
        ## bkg_r_c.Draw("HISTs")
        ## data_r_c = histos[r]["data"].Clone("data_clone__"+r)
        ## data_r_c.SetDirectory(0)
        ## data_r_c.Draw("sames")
        ## c.SaveAs(dist_name+"_CR_Quick_"+r+".root")
        

        #bkg_r.Divide( histos[r]["data"] )
        #ratio = bkg_r
        
        histos[r]["data"].Divide(  bkg_r )
        ratio = histos[r]["data"]

        #search for last bin with data, will be used for upper fit range
        lastbin = 0
        for ibin in reversed(range(ratio.GetNbinsX()+1)):
            if ratio.GetBinContent(ibin) != 0:
                lastbin = ibin
                break

        lastbin_Xval = ratio.GetBinLowEdge(lastbin) + ratio.GetBinWidth(lastbin)

        #fitting
        fitName = "fit_"+outfileNameBase[:-5]
    
        npar = 2
        fitChoice = LinearFunc
        #fitRange = [bkg_r.GetBinLowEdge(1),  bkg_r.GetBinLowEdge(bkg_r.GetNbinsX()) + bkg_r.GetBinWidth(bkg_r.GetNbinsX()) ]
        fitRange = [500, lastbin_Xval]
        
        func = R.TF1(fitName, fitChoice, fitRange[0], fitRange[1], 2)
        func.SetParameters(1.0, 0.0)

        Vmode = ("Q" if not verbose else "")
        fitResult = ratio.Fit(fitName, "S0"+Vmode, "", fitRange[0], fitRange[1]) # Qi Question: Why we fit on data/TotalBkgEst ? Shouldn't we fit on (data-ttbar)/(QCD Est.) ? Afterare, this is what we applied in SR

        
        if fitResult.Status() != 0:
            print "Error in QCD Sys fit: did not terminate properly. Exiting"
            sys.exit(0)

        cov_TMatrix = fitResult.GetCovarianceMatrix()
        cov = np.zeros( (npar, npar) )
        for i in range(npar):
            for j in range(npar):
                cov[i,j] = cov_TMatrix[i][j]
        

        fitFunc = ratio.GetFunction(fitName)

        params_array = array('d',[0]*npar)
        fitFunc.GetParameters( params_array )

        params = np.asarray( params_array )


        outRange = [500, 3500]

        fcen = R.TF1("QCDShape_f_"+r, fitChoice, outRange[0], outRange[1], npar)
        fcen.SetParameters( params )

        fup = R.TF1("QCDShape_fup_"+r, fitChoice, outRange[0], outRange[1], npar)
        fup.SetParameters( params[0] - np.sqrt(cov[1,1])*(fitRange[1]+fitRange[0])/2.0, params[1] + np.sqrt(cov[1,1]) )
        fup.SetLineColor(R.kBlue)

        fdw = R.TF1("QCDShape_fdw_"+r, fitChoice, outRange[0], outRange[1], npar)
        fdw.SetParameters( params[0] + np.sqrt(cov[1,1])*(fitRange[1]+fitRange[0])/2.0, params[1] - np.sqrt(cov[1,1]) )
        fdw.SetLineColor(R.kBlue)

        QCDSyst_Dict["Shape_"+r] = {"f":fcen, "fup":fup, "fdw":fdw}
        QCDSyst_Dict["Scale_"+r] = np.max( np.abs( [ (1.0-params[0]),  (1.0 / np.sqrt(histos[r]["data"].Integral())) ] ) ) #scale is max of ratio non-unity and CR stat error # Qi Question. histos[r]["data"] has become the ratio before?!

        if makePlots:
            c=R.TCanvas()
            #R.SetOwnership(c,False)
            ratio.SetLineColor(R.kBlack)
            ratio.Draw()
            fitFunc.Draw("same")
            fup.Draw("same")
            fdw.Draw("same")
            c.SaveAs(outfileNameBase.split(".root")[0] + "_" + r + ".root")


    datafile.Close()
    topfile.Close()
    
    return QCDSyst_Dict


def LinearFunc(x, par):
    return (par[0] + par[1]*x[0])




'''
derivation of constant term under variation of slope,
which should keep overall normalization fixed 
by ensuring that the scaling function has the same integral
before and after slope change

int_l^h a + bx = ax + bx^2 / 2  (l...h) = a*(h-l) + b/2 * (h^2 - l^2)

a*(h-l) + b/2 * (h^2 - l^2) = z*(h-l) + (b+d)/2 * (h^2 - l^2)

z*(h-l) = a*(h-l)  + b/2 * (h^2 - l^2) - (b)/2 * (h^2 - l^2) - (d)/2 * (h^2 - l^2)

z*(h-l) = a*(h-l) - (d)/2 * (h^2 - l^2) = a*(h-l) - (d)/2 * (h - l)*(h+l)

z= a - d/2*(h+l)

'''
