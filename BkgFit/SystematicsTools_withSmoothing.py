import ROOT as R
import numpy as np
from array import array
import sys
from copy import deepcopy
from GetEigenVariations import GetEigenVariations
from HistoTools import HistLocationString as HistLocStr
from HistoTools import CheckAndGet
import smoothfit


# hard-coded in!!!!
_extraNormCRSysDict = {
     "44": 0.148,
     "33": 0.0708,
     "22": 0.0285
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
                    smoothing_func = "Dijet",
                    SmoothRange = (1100, 3000),# (100, 2500),
                    use_one_top_nuis = False,
                    use_scale_top_0b = False,
                    nbtag_top_shape_for4b = None,
                    makePlots = False,
                    verbose = False,
                    outfileNameBase="QCDSysfitSmooth.root"):
    ##this cannot have the norm fixed
    # global rfunc1
    # global rfunc2
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


    colorlist = [ R.kGreen, R.kOrange, R.kMagenta, R.kCyan, R.kPink, (R.kAzure+1), R.kGreen+2, R.kOrange+5]        


    ##### Get Signal Region Histograms ################################
    datafile = R.TFile(datafileName,"READ")
    topfile  = R.TFile(topfileName,"READ")
    zjetfile  = ( R.TFile(zjetfileName,"READ") if zjetfileName!=None else None)


    histos = {}
    
    # collect all histograms
    for r in ["44","33","22","40","30","20"]:
        folder_r = HistLocStr(dist_name, r[0], r[1], btag_WP, "CR")  #folder( r[0], r[1], btag_WP)
        print folder_r
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
        #bkg_r.Add( zjet_r )


        N_bkg_r = bkg_r.Integral()

        Err_N_data_CR_r = R.Double(0)
        N_data_CR_r =  histos[r]["data"].IntegralAndError(0, histos[r]["data"].GetNbinsX()+1, Err_N_data_CR_r)


        bkg_r.Scale(histos[r]["data"].Integral() / bkg_r.Integral())


        c=R.TCanvas("c1_cr_"+r,"c1_cr_"+r)
        xleg, yleg = 0.52, 0.7
        leg = R.TLegend(xleg, yleg, xleg+0.3, yleg+0.2)
        leg.SetFillColor(0)
        leg.SetBorderSize(0)
        leg.SetMargin(0.3)
        histos[r]["data"].SetXTitle("m_{JJ} [GeV]")
        histos[r]["data"].SetYTitle("Entries")
        histos[r]["data"].GetXaxis().SetRangeUser(500, 4000)
        histos[r]["data"].Draw("E1")
        leg.AddEntry(histos[r]["data"], "CR data", "LP")

        ##################################
        ## smooth bkg and data
        ##################################
        data_sm   = smoothfit.smoothfit(histos[r]["data"], fitFunction = smoothing_func, fitRange = SmoothRange, makePlots = False, verbose = False, useLikelihood=True, outfileName="data_smoothfit_CRsyst_"+r+".root")
        data_sm_h = smoothfit.MakeSmoothHisto(histos[r]["data"], data_sm["nom"], keepNorm=False)

        data_sm["nom"].SetNameTitle("data_smoothfit_CRsyst_"+r,"data_smoothfit_CRsyst_"+r)
        data_sm["nom"].SetLineColor(R.kBlack)
        data_sm["nom"].Draw("same")
        leg.AddEntry(data_sm["nom"], "CR data smoothed", "L")

        bkg_sm   = smoothfit.smoothfit(bkg_r, fitFunction = smoothing_func, fitRange = SmoothRange, makePlots = False, verbose = False, outfileName="bkg_smoothfit_CRsyst_"+r+".root")
        bkg_sm_h = smoothfit.MakeSmoothHisto(bkg_r, bkg_sm["nom"], keepNorm=False)

        bkg_sm["nom"].SetNameTitle("bkg_smoothfit_CRsyst_"+r,"bkg_smoothfit_CRsyst_"+r)
        bkg_sm["nom"].SetLineColor(R.kBlue)
        bkg_sm["nom"].Draw("same")
        leg.AddEntry(bkg_sm["nom"], "CR Prediction smoothed", "L")


  
        rfunc1 = data_sm["nom"]
        rfunc2 = bkg_sm["nom"]
        def rfunc_ratio(x):
            return rfunc1.Eval(x[0]) / rfunc2.Eval(x[0])
        xMax   = histos[r]["data"].GetXaxis().GetBinUpEdge(histos[r]["data"].GetXaxis().GetNbins())
        ratio_sm = R.TF1("ratio_crsys_sm"+r, rfunc_ratio, SmoothRange[0], xMax, 0)
        ## ratio_sm.SetLineColor(R.kGray)



        for ivar in range(len(data_sm["vars"])):
            dup = data_sm["vars"][ivar][0]
            ddw = data_sm["vars"][ivar][1]

            dup.SetLineColor(colorlist[ivar])
            ddw.SetLineColor(colorlist[ivar])

            dup.Draw("same")
            ddw.Draw("same")
            leg.AddEntry(dup, "CR data smoothed variation", "L")


            data_r_qup = smoothfit.MakeSmoothHisto(histos[r]["data"], dup, keepNorm=False)
            data_r_qdw = smoothfit.MakeSmoothHisto(histos[r]["data"], ddw, keepNorm=False)

            for ibin in range(1,  data_sm_h.GetNbinsX()+1):
                err_val = np.max( np.abs( [ data_sm_h.GetBinContent(ibin) - data_r_qup.GetBinContent(ibin), data_sm_h.GetBinContent(ibin) - data_r_qdw.GetBinContent(ibin)] ) )
                data_sm_h.SetBinError(ibin, np.sqrt( data_sm_h.GetBinError(ibin)**2 + err_val**2) )

        c.SetLogy(1)
        leg.Draw("same")
        c.SaveAs(outfileNameBase.split(".root")[0] + "_" + r + ".root")
        c.SaveAs(outfileNameBase.split(".root")[0] + "_" + r + ".pdf")
        c.Close()


        h_ratio_cr_nom = data_sm_h.Clone("data_sm_h_CRsyst_"+r)
        h_ratio_cr_nom.Divide( data_sm["nom"] )
        h_ratio_cr_nom.SetDirectory(0)
                
        h_ratio_cr = data_sm["nom"].GetHistogram()
        h_ratio_cr.Divide( bkg_sm["nom"] )
        h_ratio_cr.SetDirectory(0)


        

        QCDSyst_Dict["Shape_"+r] = ratio_sm.Clone(ratio_sm.GetName() + r)

        #scale is max of ratio non-unity and CR stat error 
        QCDSyst_Dict["Scale_"+r] = np.max( np.abs( [ (N_bkg_r - N_data_CR_r)/N_bkg_r,  (Err_N_data_CR_r / N_data_CR_r), _extraNormCRSysDict.get(r, 0.) ] ) )
        print "Scale_"+r, QCDSyst_Dict["Scale_"+r], N_bkg_r, N_data_CR_r, Err_N_data_CR_r,  (N_bkg_r - N_data_CR_r)/N_bkg_r, Err_N_data_CR_r / N_data_CR_r
        
        
        c2=R.TCanvas("c2_cr_"+r,"c2_cr_"+r)
        leg = R.TLegend(0.2,0.7,0.5,0.9)
        leg.SetFillColor(0)
        h_ratio_cr_nom.SetFillColor(R.kBlack)
        h_ratio_cr_nom.SetFillStyle(3004)
        h_ratio_cr_nom.SetMarkerSize(0)
        h_ratio_cr_nom.GetXaxis().SetRangeUser(1000, 4000)
        h_ratio_cr_nom.GetYaxis().SetRangeUser(0, 3)
        h_ratio_cr_nom.GetXaxis().SetLabelSize(0.04)
        h_ratio_cr_nom.GetYaxis().SetLabelSize(0.04)
        h_ratio_cr_nom.SetXTitle("m_{JJ} [GeV]")
        h_ratio_cr_nom.SetYTitle("Ratio")
        h_ratio_cr_nom.Draw("E2")
        leg.AddEntry(h_ratio_cr_nom, "CR data", "LF")


        h_ratio_cr.SetLineColor( R.kBlue )
        h_ratio_cr.Draw("same")
        leg.AddEntry(h_ratio_cr, "CR Prediction", "L")


        ratio_sm.Draw("same")
        leg.Draw("same")
        c2.SaveAs(outfileNameBase.split(".root")[0] + "_ratio_" + r + ".root")
        c2.SaveAs(outfileNameBase.split(".root")[0] + "_ratio_" + r + ".pdf")
        c2.Close()

        ratio_sm = None
        rfunc1 = None
        rfunc2 = None

    datafile.Close()
    topfile.Close()
    
    return QCDSyst_Dict


def ttbarShapeSysSR(topfileName="hist_ttbar.root",
                    distributionName= "mHH_l",
                    signal_region  = "33",
                    compare_region = "44",
                    btag_WP        = "77",
                    smoothing_func = "Exp",
                    SmoothRange = (1200, 3000),# (100, 2500),
                    makePlots = False,
                    verbose = False,
                    outfileNameBase="TopShapeSRSysfitSmooth.root"):
    
    topfile  = R.TFile(topfileName,"READ")

    ttbarShapeSRSyst_Dict= {}

    colorlist = [ R.kGreen, R.kOrange, R.kMagenta, R.kCyan, R.kPink, (R.kAzure+1), R.kGreen+2, R.kOrange+5]        

    
    ## get top SR shape
    folder_sig = HistLocStr(distributionName, signal_region[0], signal_region[1], btag_WP, "SB")  #folder( r[0], r[1], btag_WP)
    top_sig    = topfile.Get(folder_sig).Clone("top_sig_"+signal_region)
    top_sig.SetDirectory(0)
    top_sig.Rebin(10)


    ## get top comparison shape
    folder_comp = HistLocStr(distributionName, compare_region[0], compare_region[1], btag_WP, "SB")  #folder( r[0], r[1], btag_WP)
    top_comp    = topfile.Get(folder_comp).Clone("top_comp_"+compare_region)
    top_comp.SetDirectory(0)
    top_comp.Rebin(10)

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
    top_sig.Scale( top_comp.Integral() / top_sig.Integral() )



    c=R.TCanvas("c1_topsys","c1_topsys")
    xleg, yleg = 0.52, 0.7
    leg = R.TLegend(xleg, yleg, xleg+0.3, yleg+0.2)
    leg.SetFillColor(0)
    leg.SetBorderSize(0)
    leg.SetMargin(0.3)
    top_comp.SetXTitle("m_{JJ} [GeV]")
    top_comp.SetYTitle("Entries")
    top_comp.Draw("E1")
    leg.AddEntry(top_comp, "Top Comparison Distribution", "LP")

    #################################
    ## smooth bkg and data
    ##################################
    top_comp_sm = smoothfit.smoothfit(top_comp, fitFunction = smoothing_func, fitRange = SmoothRange, makePlots = False, verbose = False, outfileName="top_comp_smoothfit_TopShape4b.root")
    top_comp_sm_h = smoothfit.MakeSmoothHisto(top_comp, top_comp_sm["nom"])

    top_comp_sm["nom"].SetLineColor(R.kBlack)
    top_comp_sm["nom"].Draw("same")
    leg.AddEntry(top_comp_sm["nom"], "Top Comparison Distribution Smooth", "L")


    top_sig_sm = smoothfit.smoothfit(top_sig, fitFunction = smoothing_func, fitRange = SmoothRange, makePlots = False, verbose = False, outfileName="top_sig_smoothfit_TopShape4.root")
    top_sig_sm_h = smoothfit.MakeSmoothHisto(top_sig, top_sig_sm["nom"])

    top_sig_sm["nom"].SetLineColor(R.kBlue)
    top_sig_sm["nom"].Draw("same")
    leg.AddEntry(top_sig_sm["nom"], "Top Nominal Distribution Smooth", "L")


    rfunc1 = top_comp_sm["nom"]
    rfunc2 = top_sig_sm["nom"]
    def rfunc_ratio(x):
        return rfunc1.Eval(x[0]) / rfunc2.Eval(x[0])
    xMax = top_comp.GetXaxis().GetBinUpEdge(top_comp.GetXaxis().GetNbins())
    ratio_sm = R.TF1("ratio_topsys_sm", rfunc_ratio, SmoothRange[0], xMax, 0)

    for ivar in range(len(top_comp_sm["vars"])):
        dup = top_comp_sm["vars"][ivar][0]
        ddw = top_comp_sm["vars"][ivar][1]

        dup.SetLineColor(colorlist[ivar])
        ddw.SetLineColor(colorlist[ivar])

        dup.Draw("same")
        ddw.Draw("same")
        leg.AddEntry(dup, "Top Comparison Smooth Variation", "L")


        top_comp_r_qup = smoothfit.MakeSmoothHisto(top_comp, dup)
        top_comp_r_qdw = smoothfit.MakeSmoothHisto(top_comp, ddw)

        for ibin in range(1,  top_comp_sm_h.GetNbinsX()+1):
            err_val = np.max( np.abs( [ top_comp_sm_h.GetBinContent(ibin) - top_comp_r_qup.GetBinContent(ibin), top_comp_sm_h.GetBinContent(ibin) - top_comp_r_qdw.GetBinContent(ibin)] ) )
            top_comp_sm_h.SetBinError(ibin, np.sqrt( top_comp_sm_h.GetBinError(ibin)**2 + err_val**2) )

    c.SetLogy(1)
    leg.Draw("same")
    c.SaveAs(outfileNameBase.split(".root")[0] + "_sig"+signal_region+"_comp"+ compare_region + ".root")
    c.SaveAs(outfileNameBase.split(".root")[0] + "_sig"+signal_region+"_comp"+ compare_region + ".pdf")
    c.Close()


    h_ratio_cr_nom = top_comp_sm_h.Clone("top_comp_sm_h_TopShape4b")
    h_ratio_cr_nom.Divide( top_comp_sm["nom"] )
    h_ratio_cr_nom.SetDirectory(0)
            
    h_ratio_cr = top_comp_sm["nom"].GetHistogram()
    h_ratio_cr.Divide( top_sig_sm["nom"] )
    h_ratio_cr.SetDirectory(0)
        

    ttbarShapeSRSyst_Dict["Shape"] = ratio_sm
        
        
    c2=R.TCanvas("c2_topsys","c2_topsys")
    leg = R.TLegend(0.2,0.7,0.5,0.9)
    leg.SetFillColor(0)
    h_ratio_cr_nom.SetFillColor(R.kBlack)
    h_ratio_cr_nom.SetFillStyle(3004)
    h_ratio_cr_nom.SetMarkerSize(0)
    h_ratio_cr_nom.GetXaxis().SetRangeUser(1000, 4000)
    h_ratio_cr_nom.GetYaxis().SetRangeUser(0, 3)
    h_ratio_cr_nom.GetXaxis().SetLabelSize(0.04)
    h_ratio_cr_nom.GetYaxis().SetLabelSize(0.04)
    h_ratio_cr_nom.SetXTitle("m_{JJ} [GeV]")
    h_ratio_cr_nom.SetYTitle("Ratio")
    h_ratio_cr_nom.Draw("E2")
    leg.AddEntry(h_ratio_cr_nom, "Nominal", "LF")


    h_ratio_cr.SetLineColor( R.kBlue )
    h_ratio_cr.Draw("same")
    leg.AddEntry(h_ratio_cr, "Predicted", "L")


    #ratio_sm.Draw("same")

    leg.Draw("same")
    c2.SaveAs(outfileNameBase.split(".root")[0] + "_sig"+signal_region+"_comp"+ compare_region +"_ratio.root")
    c2.SaveAs(outfileNameBase.split(".root")[0] + "_sig"+signal_region+"_comp"+ compare_region +"_ratio.pdf")
    c2.Close()

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
