import ROOT as R

import numpy as np
from array import array
import sys
from copy import deepcopy

from GetEigenVariations import GetEigenVariations

from HistoTools import HistLocationString as HistLocStr
from HistoTools import CheckAndGet


# hard-coded in!!!!
_extraNormCRSysDict = {
    # "44": 0.258,
    # "43": 0.126,
}

def QCDSystematics(datafileName="hist_data.root",
                    topfileName="hist_ttbar.root",
                    zjetfileName="hist_Zjets.root",
                    distributionName= "mHH_l",
                    n_trkjet  = ["4","3","2"],
                    n_btag    = ["4","3","2"],
                    btag_WP     = "77",
                    mu_qcd_vals = [1.0, 1.0],
                    topscale_vals = [1.0, 1.0],
                    NRebin = 1,
                    use_one_top_nuis = False,
                    use_scale_top_0b = False,
                    nbtag_top_shape_for4b = None,
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

    scaleTop0b = use_scale_top_0b

    n_channels = num_trkjet.shape[0]

    regions = [ num_trkjet[i]+num_btag[i] for i in range(n_channels) ]
    ##################################################################


    ##### Get Signal Region Histograms ################################
    datafile = R.TFile(datafileName,"READ")
    topfile  = R.TFile(topfileName,"READ")
    zjetfile  = ( R.TFile(zjetfileName,"READ") if zjetfileName!=None else None)


    histos = {}
    
    # collect all histograms
    for r in ["44","33","22","40","30","20"]:
        folder_r = HistLocStr(dist_name, r[0], r[1], btag_WP, "CR")  #folder( r[0], r[1], btag_WP)
        
        data_r   = datafile.Get(folder_r).Clone("data_"+r)
        data_r.SetDirectory(0)
        
        top_r    = topfile.Get(folder_r).Clone("top_"+r)
        top_r.SetDirectory(0)

        zjet_r   = CheckAndGet(zjetfile, folder_r, top_r).Clone("zjet_"+r)
        zjet_r.SetDirectory(0)

        for ibin in range(1, top_r.GetNbinsX()+1):
            if top_r.GetBinContent(ibin) < 0:
                top_r.SetBinContent(ibin, 0)
                top_r.SetBinError(ibin, 0)

        data_r.Rebin(n_rebin)
        top_r.Rebin(n_rebin)
        zjet_r.Rebin(n_rebin)


        histos[r]     = {"data": data_r,  "top": top_r,  "zjet":zjet_r}

    datafile.Close()
    topfile.Close()
    if zjetfile != None:
        zjetfile.Close()
    ##################################################################

    ####### outpue object ###################
    QCDSyst_Dict = {}
    

    ##### scaling and subtractions #################################
    for ir in range(len(regions)):
        r = regions[ir]
                
        r_0b = r[0]+"0"
        r_3b = r[0]+"3"

        top_0b = histos[r_0b]["top"].Clone("top_0b__"+r)
        if scaleTop0b:
            top_0b.Scale( (topscale_vals[0] if use_one_top_nuis else topscale_vals[ir]) )

        zjet_0b = histos[r_0b]["zjet"].Clone("zjet_0b__"+r)


        qcd_r = histos[r_0b]["data"].Clone("qcd__"+r)
        qcd_r.Add( top_0b, -1)      # added by Qi --- we still want top to be subtracted, given that their fraction is increasing in Run 2.
        qcd_r.Add( zjet_0b, -1)
        qcd_int = qcd_r.Integral()

        if nbtag_top_shape_for4b != None:
            top_r = histos[nbtag_top_shape_for4b]["top"].Clone("top__"+r)
            top_r.Scale( histos[r]["top"].Integral() / top_r.Integral() ) #scale to correct norm for region
        else:
            top_r = histos[r]["top"].Clone("top__"+r)

        top_int = top_r.Integral()

        zjet_r = histos[r]["zjet"].Clone("zjet__"+r)



        mu_qcd = mu_qcd_vals[ir]
        top_scale = (topscale_vals[0] if use_one_top_nuis else topscale_vals[ir])
        
        qcd_r.Scale( mu_qcd )
        top_r.Scale( top_scale )

        N_qcd_r = qcd_r.Integral()

        #now do ratio
        bkg_r = qcd_r.Clone("bkg__"+r)
        bkg_r.Add( top_r )
        bkg_r.Add( zjet_r )


        N_bkg_r = bkg_r.Integral()

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

        ratio = histos[r]["data"].Clone("ratio__"+r)
        ratio.SetDirectory(0)
        #ratio.Add( top_r, -1)

        #store integral and error
        Err_N_data_CR_r = R.Double(0)
        N_data_CR_r = ratio.IntegralAndError(0, ratio.GetNbinsX()+1, Err_N_data_CR_r)
        
        #do division
        ratio.Divide(  bkg_r )

        #search for last bin with data, will be used for upper fit range
        lastbin = 0
        for ibin in reversed(range(ratio.GetNbinsX()+1)):
            if ratio.GetBinContent(ibin) != 0:
                lastbin = ibin
                break

        lastbin_Xval = ratio.GetBinLowEdge(lastbin) + ratio.GetBinWidth(lastbin)
        fitRange = [1000, lastbin_Xval]

        

        ## fitting
        if verbose:
            print "QCD Ratio fit in SR=",r
        fitName = "fit_"+outfileNameBase[:-5]
        fitFunc, fitChoice, npar, params, cov = LinearFit(ratio, fitName, fitRange, verbose)




        
        ## Make fit variations, for shape uncertainties
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

        #scale is max of ratio non-unity and CR stat error 
        QCDSyst_Dict["Scale_"+r] = np.max( np.abs( [ (N_bkg_r - N_data_CR_r)/N_bkg_r,  (Err_N_data_CR_r / N_data_CR_r), _extraNormCRSysDict.get(r, 0.) ] ) )
        print "Scale_"+r, QCDSyst_Dict["Scale_"+r], N_bkg_r, N_data_CR_r, Err_N_data_CR_r,  (N_bkg_r - N_data_CR_r)/N_bkg_r, Err_N_data_CR_r / N_data_CR_r
        #QCDSyst_Dict["Scale_"+r] = np.max( np.abs( [ (1.0-params[0]),  (1.0 / np.sqrt(histos[r]["data"].Integral())) ] ) )  #this one was bugged a bit, keep anyway

        if makePlots:
            c=R.TCanvas()
            #R.SetOwnership(c,False)
            leg = R.TLegend(0.1,0.7,0.48,0.9)
            leg.SetFillColor(0)
            leg.AddEntry(ratio, "Ratio", "LP")
            leg.AddEntry(fcen, "Ratio Fit", "L")
            leg.AddEntry(fup, "Ratio Fit Variations", "L")
            
            ratio.SetLineColor(R.kBlack)
            ratio.SetXTitle("m_{JJ} [GeV]")
            ratio.SetYTitle("Ratio Data/Prediction")
            ratio.Draw()
            fcen.Draw("same")
            fup.Draw("same")
            fdw.Draw("same")
            leg.Draw("same")
            c.SaveAs(outfileNameBase.split(".root")[0] + "_" + r + ".root")
            c.SaveAs(outfileNameBase.split(".root")[0] + "_" + r + ".pdf")


    datafile.Close()
    topfile.Close()
    
    return QCDSyst_Dict



def ttbarShapeSysSR(topfileName="hist_ttbar.root",
                    distributionName= "mHH_l",
                    signal_region = "33",
                    compare_region = "44",
                    btag_WP     = "77",
                    makePlots = False,
                    verbose = False,
                    outfileNameBase="TopShapeSRSysfit.root"):
    
    topfile  = R.TFile(topfileName,"READ")
    
    ## get top SR shape
    folder_sig = HistLocStr(distributionName, signal_region[0], signal_region[1], btag_WP, "SR")  #folder( r[0], r[1], btag_WP)
    top_sig    = topfile.Get(folder_sig).Clone("top_sig_"+signal_region)
    top_sig.SetDirectory(0)
    top_sig.Rebin(5)


    ## get top comparison shape
    folder_comp = HistLocStr(distributionName, compare_region[0], compare_region[1], btag_WP, "SR")  #folder( r[0], r[1], btag_WP)
    top_comp    = topfile.Get(folder_comp).Clone("top_comp_"+compare_region)
    top_comp.SetDirectory(0)
    top_comp.Rebin(5)

    ## remove negative values
    ## assume same binning, else division won't work later
    for ibin in range(1, top_sig.GetNbinsX()+1):
        if top_sig.GetBinContent(ibin) < 0:
            top_sig.SetBinContent(ibin, 0)
            top_sig.SetBinError(ibin, 0)
            
        if top_comp.GetBinContent(ibin) < 0:
            top_comp.SetBinContent(ibin, 0)
            top_comp.SetBinError(ibin, 0)

    ## normalize to same area
    top_sig.Scale( 1.0 / top_sig.Integral() )
    top_comp.Scale( 1.0 / top_comp.Integral() )

    ## compute ratio
    top_ratio = top_comp.Clone("top_ratio_sig"+signal_region+"_comp"+compare_region)
    top_ratio.Divide( top_sig )

    ## search for last bin with data, will be used for upper fit range
    lastbin = 0
    for ibin in reversed(range(top_ratio.GetNbinsX()+1)):
        if top_ratio.GetBinContent(ibin) != 0:
            lastbin = ibin
            break

    lastbin_Xval = top_ratio.GetBinLowEdge(lastbin) + top_ratio.GetBinWidth(lastbin)
    #fitRange = [500, lastbin_Xval]
    #print "ttbar fit range", fitRange
    fitRange = [500, 1600]  ## reasonable range of bins with data

    ## fitting
    fitName = "fit_"+outfileNameBase[:-5]
    fitFunc, fitChoice, npar, params, cov = LinearFit(top_ratio, fitName, fitRange, verbose)


    ## Make fit variations, for shape uncertainties
    outRange = [500, 3500]

    fcen = R.TF1("ttbarShapeSR_f_sig"+signal_region+"_comp"+compare_region, fitChoice, outRange[0], outRange[1], npar)
    fcen.SetParameters( params )

    frev = R.TF1("ttbarShapeSR_frev_sig"+signal_region+"_comp"+compare_region, fitChoice, outRange[0], outRange[1], npar)
    frev.SetParameters( 2 - params[0], -params[1] )
    frev.SetLineColor(R.kOrange)

    fneg = R.TF1("ttbarShapeSR_fneg_sig"+signal_region+"_comp"+compare_region, fitChoice, outRange[0], outRange[1], npar)
    fneg.SetParameters(params[0], -params[1] )
    fneg.SetLineColor(R.kMagenta)

    fup = R.TF1("ttbarShapeSR_fup_sig"+signal_region+"_comp"+compare_region, fitChoice, outRange[0], outRange[1], npar)
    fup.SetParameters( params[0] - np.sqrt(cov[1,1])*(fitRange[1]+fitRange[0])/2.0, params[1] + np.sqrt(cov[1,1]) )
    fup.SetLineColor(R.kBlue)

    fdw = R.TF1("ttbarShapeSR_fdw_sig"+signal_region+"_comp"+compare_region, fitChoice, outRange[0], outRange[1], npar)
    fdw.SetParameters( params[0] + np.sqrt(cov[1,1])*(fitRange[1]+fitRange[0])/2.0, params[1] - np.sqrt(cov[1,1]) )
    fdw.SetLineColor(R.kBlue)

    ttbarShapeSRSyst_Dict = {"f":fcen, "frev":frev, "fup":fup, "fdw":fdw}

    if makePlots:
        c=R.TCanvas()
        #R.SetOwnership(c,False)
        leg = R.TLegend(0.1,0.7,0.48,0.9)
        leg.SetFillColor(0)
        leg.AddEntry(top_ratio, "Ratio", "LP")
        leg.AddEntry(fcen, "Ratio Fit", "L")
        leg.AddEntry(fup, "Slope Variations", "L")
        
        top_ratio.SetLineColor(R.kBlack)
        top_ratio.SetXTitle("m_{JJ} [GeV]")
        top_ratio.SetYTitle("Ratio Data/Prediction")
        top_ratio.Draw()
        fcen.Draw("same")
        #frev.Draw("same")
        #fneg.Draw("same")
        fup.Draw("same")
        fdw.Draw("same")
        leg.Draw("same")
        c.SaveAs(outfileNameBase.split(".root")[0] + "_sig"+signal_region+"_comp"+compare_region+ ".root")


    topfile.Close()
    
    return ttbarShapeSRSyst_Dict






def LinearFit(histo, fitName, fitRange, verbose):

    npar = 2
    fitChoice = LinearFunc
    #fitRange = [bkg_r.GetBinLowEdge(1),  bkg_r.GetBinLowEdge(bkg_r.GetNbinsX()) + bkg_r.GetBinWidth(bkg_r.GetNbinsX()) ]

    func = R.TF1(fitName, fitChoice, fitRange[0], fitRange[1], 2)
    func.SetParameters(1.0, 0.0)

    Vmode = ("Q" if not verbose else "")
    fitResult = histo.Fit(fitName, "S0"+Vmode, "", fitRange[0], fitRange[1]) 


    if fitResult.Status() != 0:
        print "Error in  SystematicsTools LinearFit: did not terminate properly. Exiting"
        sys.exit(0)

    cov_TMatrix = fitResult.GetCovarianceMatrix()
    cov = np.zeros( (npar, npar) )
    for i in range(npar):
        for j in range(npar):
            cov[i,j] = cov_TMatrix[i][j]


    fitFunc = histo.GetFunction(fitName)

    params_array = array('d',[0]*npar)
    fitFunc.GetParameters( params_array )

    params = np.asarray( params_array )

    return fitFunc, fitChoice, npar, params, cov





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


'''
fine where ax + b = 1  then reverse slope around that point

ax+b = 1
x = (1-b)/a

what is new offset, slope is -a?
-ax + c = 1  when x = (1-b)/a

-a(1-b)/a + c = 1
b - 1 + c = 1
c = 2-b

check
-ax + 2 - b
-a(1-b)/a +2 - b = b-1 + 2 - b = 1  
'''
