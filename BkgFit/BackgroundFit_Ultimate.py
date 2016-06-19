import ROOT as R
import numpy as np
from array import array
import sys
import os
from copy import deepcopy
from GetEigenVariations import GetEigenVariations
from HistoTools import HistLocationString as HistLocStr
from HistoTools import CheckAndGet
R.gROOT.LoadMacro("AtlasStyle.C") 
R.SetAtlasStyle()


regions = {}
h_qcd = {}#qcd histogram
h_top = {}#top histogram
h_top_model = {}#top shape modeling
h_zjet = {}#zjet histogram
h_zjet_model = {}#zjet modeling
h_data = {}#data modeling

useOneTopNuis = None
scaleTop_model = None

def BackgroundFit(datafileName        ="hist_data.root",
                  topfileName         ="hist_ttbar.root",
                  zjetfileName        ="hist_Zjets.root",
                  distributionName    = ["LeadCaloJetM"],
                  n_trkjet            = ["4", "3", "2", "1"],
                  n_btag              = ["4", "3", "2s", "2", "1"], #["4", "3", "2s", "2", "1"],
                  btag_WP             = "77", #not useful for Xhh Framework
                  NRebin              = 1,
                  BKG_model           = 0, #define the bkg models
                  use_one_top_nuis    = False, #True to fix one top parameter
                  use_scale_top_model = False,
                  nbtag_top_shape     = True, #fix 4b and 3b shape to be the same
                  makePlots           = True,
                  whichFunc           = "XhhBoosted",
                  output              = "",
                  verbose             = False ):
    
    global h_qcd
    global h_top
    global h_top_model
    global h_data
    global dist_name
    global useOneTopNuis
    global scaleTop_model
    global regions
    global Output
    ################### Parse  #########################
    # num_trkjet  = np.asarray(n_trkjet)
    # if num_trkjet.shape==():
    #     num_trkjet = np.asarray([n_trkjet])
    # num_btag    = np.asarray(n_btag)
    # if num_btag.shape==():
    #     num_btag = np.asarray([n_btag])
    # if num_btag.shape != num_trkjet.shape:
    #     print "Must have same number of track jet and b-tag regions specified"
    #     sys.exit(0)
    btag_WP     = btag_WP
    n_rebin     = NRebin
    #setup top shape constrains
    useOneTopNuis = use_one_top_nuis
    scaleTop_model = use_scale_top_model
    dist_name = distributionName
    ########################################################
    #setup regions to fit
    regions = [ "i" + n_btag[i] for i in range(len(n_btag)) ]
    #print regions
    ########################################################
    #load the histogram files
    datafile = R.TFile(datafileName,"READ")
    topfile  = R.TFile(topfileName,"READ")
    zjetfile = R.TFile(zjetfileName,"READ")
    Output   = output
    ########################################################
    ################### Setup Minuit  ###################
    n_param = len(regions) + ( 1 if useOneTopNuis else len(regions) )
    minuit = R.TMinuit(n_param) # 2 parameter fit
    minuit.SetPrintLevel( (1 if verbose else -1) )
    minuit.SetErrorDef(0.5)
    minuit.SetFCN(NegLogL)
    #########################################################
    ################### Get Histograms  ###################
    histos = { }
    # collect all histograms; ntrkjets, nbtags
    hist_region_lst = ["i" + x for x in n_btag]
    #load the specific trackjet regions
    hist_region_lst.append("i0")
    hist_region_lst.append("2" + str(BKG_model))
    hist_region_lst.append("3" + str(BKG_model))
    hist_region_lst.append("4" + str(BKG_model))
    #print hist_region_lst
    #load the histograms
    for r in hist_region_lst:
        data_r = {}
        top_r = {}
        zjet_r = {}
        for h in dist_name:
            hist_fullpath = HistLocStr(h, r[0], r[1:], btag_WP, "Sideband", whichFunc)  #folder( r[0], r[1], btag_WP)
            #print hist_fullpath
            data_r[h] = datafile.Get(hist_fullpath).Clone("data_"+r+h)
            top_r[h]  = topfile.Get(hist_fullpath).Clone("top_"+r+h)
            zjet_r[h] = CheckAndGet(zjetfile, hist_fullpath, top_r).Clone("zjet_"+r+h)
            # do rebin if necessary
            if NRebin > 1:
                data_r[h].Rebin(NRebin)
                top_r[h].Rebin(NRebin)
                zjet_r[h].Rebin(NRebin)
        
        histo_r  = {"data": data_r, "top": top_r, "zjet": zjet_r}
        histos[r] = histo_r
    # put relevant histograms in global lists for use in LogLk
    for r in regions:
        h_data[r] = {} #this is the data to fit
        h_qcd[r] = {} #this is the qcd fit
        h_top[r] = {} #this is the top fit
        h_top_model[r] = {}
        h_zjet[r] = {}
        h_zjet_model[r] = {}
        for h in dist_name:
            hd = histos[r]["data"][h].Clone("h_data_"+r+h)
            if nbtag_top_shape != None:
                ht = histos[r]["top"][h].Clone("h_top_"+r+h)
                #change for top selection, from 4b to 3b
                if r[1:] == "4":
                    ht = histos[r[0] + "3"]["top"][h].Clone("h_top_"+r+h)
                ht.Scale( histos[r]["top"][h].Integral() / ht.Integral() ) #scale to correct norm for region
            else:
                ht = histos[r]["top"][h].Clone("h_top_"+r+h)
            hz = histos[r]["zjet"][h].Clone("h_zjet_"+r+h)

            #start background modeling
            bkg_model = r[0]+str(BKG_model)
            if r[1:] == "2s":
                bkg_model = "2"+str(BKG_model)
            elif r[1:] == "3":
                bkg_model = "3"+str(BKG_model)
            elif r[1:] == "4":
                bkg_model = "4"+str(BKG_model)
            #bkg_model = "42"
            #load the histograms
            hq = histos[bkg_model]["data"][h].Clone("h_qcd_"+r+h) #use 0 tag to fit now
            ht2 = histos[bkg_model]["top"][h].Clone("h_top_model_"+r+h) #use 0 tag to fit now
            hz2 = histos[bkg_model]["zjet"][h].Clone("h_zjet_"+r+h)
            #substract top and Zjet contributions from data       
            hq.Add( ht2, -1.0)
            hq.Add( hz2, -1.0)
            #link the dictionaries, now as a dictionary again
            h_data[r][h] = hd #this is the data to fit
            h_qcd[r][h] = hq #this is the qcd fit
            h_top[r][h] = ht #this is the top fit
            h_top_model[r][h] = ht2
            h_zjet[r][h] = hz 
            h_zjet_model[r][h] = hz2

    #########################################################
    #Start the fit
    results = Fit( minuit )
    #########################################################
    #Gather the results
    evars = GetEigenVariations(results["cov_m"])
    pnom  = np.asarray( results["muqcd"] + results["muttbar"] )
    pvars = [ [pnom+evars[i], pnom-evars[i]] for i in range(len(evars)) ]

    results["pnom"]  = pnom
    results["pvars"] = pvars

    # store the input histograms for fitting
    # h_store_2b_data = histos[r[0]+"2"]["data"].Clone()
    # h_store_2b_data.SetDirectory(0)
    # results["inputhist_2b_data"] = h_store_2b_data

    # h_store_2b_ttbar = histos[r[0]+"2"]["top"].Clone()
    # h_store_2b_ttbar.SetDirectory(0)
    # results["inputhist_2b_ttbar"] = h_store_2b_ttbar

    # h_store_4b_data = histos[r[0]+"4"]["data"].Clone()
    # h_store_4b_data.SetDirectory(0)
    # results["inputhist_4b_data"] = h_store_4b_data

    # h_store_4b_ttbar = histos[r[0]+"4"]["top"].Clone()
    # h_store_4b_ttbar.SetDirectory(0)
    # results["inputhist_4b_ttbar"] = h_store_4b_ttbar

    #print pnom
    #print pvars

    # print "Fit Results:"
    # print "mu_qcd = ", results["muqcd"], "+/-", results["muqcd_e"]
    # print "top_scale = ", results["muttbar"], "+/-", results["muttbar_e"]
    # print "correlation=", results["corr_m"]

    texoutpath = Output + "Tables/"
    if not os.path.exists(texoutpath):
        os.makedirs(texoutpath)
    fit_outtex = open( texoutpath + "normfit.tex", "w")
    #print len(n_btag)
    WriteFitResult(results, fit_outtex, nfit=len(n_btag))

    #ComputeBasicMuQCD( histo_s, histo_c )
    outroot = R.TFile.Open(Output + "fitNorm.root", "recreate")
    if makePlots:
        for i in range(len(regions)):
            MakePlot(regions[i], results["muqcd"][i], results["muttbar"][0 if useOneTopNuis else i])

    #finish and clean up
    outroot.Close()
    datafile.Close()
    topfile.Close()
    zjetfile.Close()

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


def MakePlot(region, muqcd, muttbar):
    for h in dist_name:

        c=R.TCanvas()
        c.SetName("fit" + region)
        c.SetFillColor(0)
        c.SetFrameFillColor(0)

        h_data2 = h_data[region][h].Clone("data2_"+region)
        h_data2.SetMarkerSize(1)
        h_data2.SetMarkerColor(R.kBlack)
        h_data2.SetFillColor(0)
        h_data2.SetLineColor(R.kBlack)
        h_data2.SetLineWidth(2)
        h_data2.SetMaximum(h_data2.GetMaximum() * 1.3)
        h_data2.SetXTitle( "Jet mass [GeV]")
        h_data2.SetYTitle( "Entries" )
        #h_data2.Rebin(nrebin)
            
        h_top2 = h_top[region][h].Clone("top2_"+region)
        h_top2.Scale( muttbar )
        h_top2.SetFillColor(R.kGreen)
        h_top2.SetLineColor(R.kBlack)
        h_top2.SetLineWidth(1)
        #h_top2.Rebin(nrebin)

        h_zjet2 = h_zjet[region][h].Clone("zjet2_"+region)
        h_zjet2.SetFillColor(R.kOrange)
        h_zjet2.SetLineColor(R.kBlack)
        h_zjet2.SetLineWidth(1)
        #h_top2.Rebin(nrebin)

        h_qcd2 = h_qcd[region][h].Clone("qcd2_"+region)
        h_qcd2.Scale( muqcd )
        h_qcd2.SetLineColor(R.kRed)
        h_qcd2.SetFillColor(0)
        h_qcd2.SetLineWidth(1)
        #h_qcd2.Rebin(nrebin)

        h_pred = h_top2.Clone("pred_"+region)
        h_pred.Add( h_qcd2, 1.0)
        h_pred.Add( h_zjet2, 1.0)
        h_pred.SetLineColor(R.kBlue)
        h_pred.SetFillColor(0)
        h_pred.SetLineWidth(1)


        h_data2.Draw("E")
        h_top2.Draw("sameHIST")
        h_zjet2.Draw("sameHIST")
        h_qcd2.Draw("sameHIST")
        h_pred.Draw("sameHIST")
        h_data2.Draw("sameE")

        leg = R.TLegend(0.65,0.7,0.9,0.9)
        leg.AddEntry(h_data2,"Data ("+region+"), 3.2 fb^{-1}","EL")
        leg.AddEntry(h_top2,"ttbar MC","F")
        leg.AddEntry(h_zjet2,"Z+jets MC","F")
        leg.AddEntry(h_qcd2,"QCD model","L")
        leg.AddEntry(h_pred,"ttbar MC + QCD model","L")
        leg.SetFillColor(0)
        leg.SetBorderSize(0)
        leg.Draw()

        l = R.TLatex(0.65, 0.65, h.replace("_", " "))
        l.SetNDC()
        l.SetTextSize(0.03)
        l.Draw("same")

        #raw_input()
        # print " "
        # print "Region "+region+": Final Numbers after fit:"
        # print "Ndata = ", h_data2.Integral()
        # print "Npred = ", h_pred.Integral()
        # print "Nqcd = ", h_qcd2.Integral()
        # print "Ntop = ", h_top2.Integral()
        # print "Nzjet = ", h_zjet2.Integral()
        # print " "

        c.Write()
        if not os.path.exists(Output + "Fit"):
            os.makedirs(Output + "Fit")
        c.SaveAs(Output + "Fit/" + "fitNorm_"+region+"_"+h+".pdf")
        c.Close()
        del(h_data2)
        del(h_top2)
        del(h_zjet2)
        del(h_qcd2)
        del(h_pred)
    return


def Fit( minuit ):
    ClearMinuit( minuit )
    
    migradStat = minuit.Migrad()
    if migradStat != 0:
        print "Error in background fit: Migrad did not terminate properly. Exiting"
        sys.exit(0)
    
    minosStat  = minuit.Command("MINOS")
    if minosStat != 0:
        print "Error in background fit: Minos did not terminate properly. Exiting"
        sys.exit(0)


    eparab = R.Double(0) #dummy
    gcc = R.Double(0)    #dummy
    tmp1 = R.Double(0)
    tmp2 = R.Double(0)

    n_reg = len(regions)
    
    muqcd         = [R.Double(0)] * n_reg
    muqcd_e       = [R.Double(0)] * n_reg
    muqcd_e_up    = [R.Double(0)] * n_reg
    muqcd_e_dw    = [R.Double(0)] * n_reg
    muttbar      = [R.Double(0)] * (1 if useOneTopNuis else n_reg)
    muttbar_e    = [R.Double(0)] * (1 if useOneTopNuis else n_reg)
    muttbar_e_up = [R.Double(0)] * (1 if useOneTopNuis else n_reg)
    muttbar_e_dw = [R.Double(0)] * (1 if useOneTopNuis else n_reg)


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
        muttbar[i]   = deepcopy(tmp1)
        muttbar_e[i] = deepcopy(tmp2)
        
        minuit.mnerrs( i+n_reg, tmp1, tmp2, eparab, gcc )
        muttbar_e_up[i] = deepcopy(tmp1)
        muttbar_e_dw[i] = deepcopy(tmp2)
        
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
            "muttbar" : muttbar,
            "muttbar_e" : muttbar_e,
            "muttbar_e_up" : muttbar_e_up,
            "muttbar_e_dw" : muttbar_e_dw,
            "cov_m": cov_m,
            "corr_m": corr_m }
    return out


def ClearMinuit( minuit ):
    minuit.Command("CLEAR")
    #initialize the parameters
    for i, reg in enumerate(regions):
        intial_muqcd = 0.1
        intial_top   = 1.2
        steps_top    = 20.0
        steps_muqcd  = 50.0
        #needs to trick the fit to offset it a bit?
        if "4" in reg:
            intial_muqcd = 0.0004
            intial_top   = 1.3
            steps_muqcd  = 200.0
            steps_top    = 200.0
        elif "3" in reg:
            intial_muqcd = 0.0055
            intial_top   = 1.3
            steps_muqcd  = 100.0
            steps_top    = 50.0
        elif "2s" in reg:
            intial_muqcd = 0.025
            intial_top   = 1.3
        elif "2" in reg:
            intial_muqcd = 0.032
            intial_top   = 3.1
        elif "1" in reg:
            intial_muqcd = 0.28
            intial_top   = 2.9
        #DefineParameter(parNo, name, initVal, initSTEP!, lowerLimit, upperLimit)
        minuit.DefineParameter(i, "muqcd_"+regions[i], intial_muqcd, intial_muqcd * 1/steps_muqcd, 0.000001, 100)
        if useOneTopNuis and i!=0:
            continue
        muttbarName = "muttbar"+("_"+regions[i] if not useOneTopNuis else '')
        # minuit.DefineParameter(i+len(regions), muttbarName, 1.3, 0.01, 0.00001, 5)
        minuit.DefineParameter(i+len(regions), muttbarName, intial_top, intial_top * 1/steps_top, 0.0001, 200)
    return
    

def NegLogL(npar, gin, f, par, ifag):
    L = 0.0

    for i in range(len(regions)):
        #setup the muqcd and ttbar scaling factors
        muqcd = par[i]
        if useOneTopNuis:
            muttbar = par[len(regions)]
        else:
            muttbar = par[i+len(regions)]
        #get the region
        r = regions[i]
        #loop over all the distributions
        for h in dist_name:
            data_r = h_data[r][h]
            qcd_r = h_qcd[r][h]
            top_r = h_top[r][h]
            zjet_r = h_zjet[r][h]
            top_model_r = h_top_model[r][h]
            Nbins = data_r.GetNbinsX()
            #add the distribution in, with the same region normalization
            for ibin in range(1,Nbins+1):
                expected_i = muqcd * qcd_r.GetBinContent(ibin) + muttbar * top_r.GetBinContent(ibin) + zjet_r.GetBinContent(ibin)
                if scaleTop_model:
                    expected_i = expected_i + (muqcd * (1.0 - muttbar) * top_model_r.GetBinContent(ibin))
                    # use (1.0 - muttbar) since top_model is already subtracted from data to make qcd
                    # so we need to first add it back, and then subtract the newly scaled top_model
                if expected_i > 0:
                    L += expected_i - (data_r.GetBinContent(ibin)) * np.log( expected_i )

    f[0] = L
    return

### 
def WriteFitResult(inputdic, outFile, nfit=3):
    ### 
    tableList = []
    ###
    tableList.append("\\begin{footnotesize}")
    tableList.append("\\begin{tabular}{c|c|c|c}")
    tableList.append("Sample & $\mu_{qcd}$ & $\\alpha_{t\\bar{t}}$ & $\\rho(\\mu_{qcd}, \\alpha_{t\\bar{t}})$ \\\\")
    tableList.append("\\hline\\hline")
    tableList.append("& & & \\\\")

    for i, cut in enumerate(regions):
    #get the mass plot
        outstr = ""
        outstr += cut.replace("i", "SB, Nb=")
        outstr += " & "
        outstr += str(round_sig(inputdic["muqcd"][i], 2))
        outstr += " $\\pm$ "
        outstr += str(round_sig(inputdic["muqcd_e"][i], 2))
        outstr += " & "
        outstr += str(round_sig(inputdic["muttbar"][i], 2))
        outstr += " $\\pm$ "
        outstr += str(round_sig(inputdic["muttbar_e"][i], 2))
        outstr += " & "
        outstr += str(round_sig(inputdic["corr_m"][i][i + nfit], 2))
        outstr+="\\\\"
        tableList.append(outstr)

    tableList.append("& & & \\\\")
    tableList.append("\\hline\\hline")
    tableList.append("\\end{tabular}")
    tableList.append("\\end{footnotesize}")
    tableList.append("\\newline")

    #return the table
    for line in tableList:
        print line
        outFile.write(line+" \n")

#round the significant numbers
def round_sig(x, sig=2):
    if x == 0:
        return 0
    if x > 1:
        return round(x, sig)
    else:
        return round(x, sig-int(R.TMath.Log10(abs(x))))

if __name__=="__main__":
    BackgroundFit()
