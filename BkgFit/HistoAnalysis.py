import ROOT as R

import numpy as np
from copy import deepcopy

import smoothfit
import BackgroundFit_MultiChannel as BkgFit


func1 = None
func2 = None

def HistoAnalysis(datafileName="hist_data.root",
                  topfileName="hist_ttbar.root",
                  distributionName= "DiJetMass",
                  n_trkjet  = ["4","4"],
                  n_btag    = ["4","3"],
                  btag_WP     = "77",
                  NRebin = 1,
                  use_one_top_nuis = False,
                  use_scale_top_2b = False,
                  nbtag_top_shape = "3"):

    global func1
    global func2
    
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

    topShape_nbtag = nbtag_top_shape
    if nbtag_top_shape == None:
        topShape_nbtag = num_btag

    useOneTopNuis = use_one_top_nuis

    scaleTop2b = use_scale_top_2b

    n_channels = num_trkjet.shape[0]

    regions = [ num_trkjet[i]+num_btag[i] for i in range(n_channels) ]
    ##################################################################


    ##### Do Background Fits ############################################
    bkgFitResults = BkgFit. BackgroundFit(datafileName=datafileName,
                                      topfileName=topfileName,
                                      distributionName= "LeadCaloJetM",
                                      n_trkjet  = n_trkjet,
                                      n_btag    = n_btag,
                                      btag_WP     = btag_WP,
                                      NRebin = NRebin,
                                      use_one_top_nuis = use_one_top_nuis,
                                      use_scale_top_2b = use_scale_top_2b,
                                      nbtag_top_shape = nbtag_top_shape,
                                      makePlots = False,
                                      verbose = False )

    ##################################################################


    ##### Get Signal Region Histograms ################################
    datafile = R.TFile(datafileName,"READ")
    topfile  = R.TFile(topfileName,"READ")

    folder = lambda nt, nb, wp: "GoodEvent_Pass" + nt + "GoodTrackJetPass" + nb + "b" + wp +"PassSRMass/"

    histos = {}
    histos_int = {}
    # collect all histograms
    for r in ["44","43","42","33","32"]:
        folder_r = folder( r[0], r[1], btag_WP)
        data_r   = datafile.Get(folder_r+dist_name).Clone("data_"+r)
        top_r    = topfile.Get(folder_r+dist_name).Clone("top_"+r)        

        histos[r]     = {"data": data_r,            "top": top_r}
        histos_int[r] = {"data": data_r.Integral(), "top": top_r.Integral()}
    

    # scaling and subtractions
    for ir in range(len(regions)):
        r = regions[ir]
        
        r_2b = r[0]+"2"
        r_3b = r[0]+"3"

        top_2b = histos[r_2b]["top"].Clone("top_2b__"+r)
        if scaleTop2b:
            top_2b.Scale( (bkgFitResults["topscale"][0] if use_one_top_nuis else bkgFitResults["topscale"][ir]) )

        qcd_r = histos[r_2b]["data"].Clone("qcd__"+r)
        #qcd_r.Add( top_2b, -1)

        top_r = histos[r]["top"].Clone("top__"+r)
        if nbtag_top_shape =="3":
            temp_scaler = top_r.Integral() / histos[r_3b]["top"].Integral()
            top_r = histos[r_3b]["top"].Clone("top__"+r)
            top_r.Scale( temp_scaler )

        qcd_r.Scale( bkgFitResults["muqcd"][ir] )
        top_r.Scale( (bkgFitResults["topscale"][0] if use_one_top_nuis else bkgFitResults["topscale"][ir]) )

        qcd_sm = smoothfit.smoothfit(qcd_r, fitFunction = "Exp", fitRange = (900, 2000), makePlots = True, outfileName="qcd_fit.root")
        top_sm = smoothfit.smoothfit(top_r, fitFunction = "Exp", fitRange = (850, 1200), makePlots = True, outfileName="top_fit.root")

        qcd_final = smoothfit.MakeSmoothHisto(qcd_r, qcd_sm["nom"])
        top_final = smoothfit.MakeSmoothHisto(top_r, top_sm["nom"])
        
        ## low=R.Double(0.0)
        ## high=R.Double(0.0)
        ## qcd_sm["nom"].GetRange(low, high)
        
        ## qcd_final = qcd_r.Clone("qcd_final__"+r)
        ## for ibin in range(1, qcd_final.GetNbinsX()+1):
        ##     if qcd_final.GetBinCenter(ibin) >= low:
        ##         qcd_final.SetBinContent(ibin, 0)
        ##         qcd_final.SetBinError(ibin, 0)
        ## qcd_final.Add(qcd_sm["nom"], 1.0)


        ## top_sm["nom"].GetRange(low, high)
        ## top_final = top_r.Clone("top_final__"+r)
        ## for ibin in range(1, top_final.GetNbinsX()+1):
        ##     if top_final.GetBinCenter(ibin) >= low:
        ##         top_final.SetBinContent(ibin, 0)
        ##         top_final.SetBinError(ibin, 0)
        ## top_final.Add(top_sm["nom"], 1.0)



        if True:
            pred_final = qcd_final.Clone("pred_final__"+r)
            pred_final.Add( top_final )


            func1 = qcd_sm["nom"]
            func2 = top_sm["nom"]

            pred_sm = R.TF1("pred_sm", FuncSum, 900, 3000)

            pred_sm.Draw("same")
            top_sm["nom"].Draw("same")

            pred_final_raw = qcd_r.Clone("qcd_final_raw__"+r)
            pred_final_raw.Add(top_r)

            outfile = R.TFile("outfile_"+r+".root","RECREATE")

            c=R.TCanvas()
            pred_final_raw.Draw("HIST")
            top_r.SetLineColor(R.kBlack)
            top_r.SetFillColor(R.kGreen)
            top_r.Draw("sameHIST")

            pred_sm.Draw("same")
            top_sm["nom"].Draw("same")

            c.Write()

            c=R.TCanvas()

            pred_final.Draw("HIST")

            top_final.SetLineColor(R.kBlack)
            top_final.SetFillColor(R.kGreen)

            top_final.Draw("sameHIST")

            c.Write()

            outfile.Close()

    


    return


def FuncSum(x):
    return ( func1.Eval(x[0]) + func2.Eval(x[0]))



if __name__=="__main__":
    HistoAnalysis()
