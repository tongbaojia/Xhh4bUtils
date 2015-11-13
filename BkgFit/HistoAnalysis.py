import ROOT as R

import numpy as np
from copy import deepcopy

import smoothfit
import BackgroundFit_MultiChannel as BkgFit

from HistoTools import HistLocationString as HistLocStr


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
                  nbtag_top_shape = None):

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

    pvars = bkgFitResults["pvars"]
    ##################################################################



    

    ##### Get Signal Region Histograms ################################
    datafile = R.TFile(datafileName,"READ")
    topfile  = R.TFile(topfileName,"READ")


    histos = {}
    
    # collect all histograms
    for r in ["44","43","42","33","32"]:
        folder_r = HistLocStr(dist_name, r[0], r[1], btag_WP, "SR")  #folder( r[0], r[1], btag_WP)
        
        data_r   = datafile.Get(folder_r).Clone("data_"+r)
        data_r.SetDirectory(0)
        
        top_r    = topfile.Get(folder_r).Clone("top_"+r)
        top_r.SetDirectory(0)     

        histos[r]     = {"data": data_r,            "top": top_r}

    datafile.Close()
    topfile.Close()
    ##################################################################



    
    ##### scaling and subtractions #################################
    for ir in range(len(regions)):
        r = regions[ir]
        
        outfileStat = R.TFile("outfile_boosted_"+r+".root","RECREATE")
        
        r_2b = r[0]+"2"
        r_3b = r[0]+"3"

        top_2b = histos[r_2b]["top"].Clone("top_2b__"+r)
        if scaleTop2b:
            top_2b.Scale( (bkgFitResults["topscale"][0] if use_one_top_nuis else bkgFitResults["topscale"][ir]) )

        qcd_r = histos[r_2b]["data"].Clone("qcd__"+r)
        qcd_int = qcd_r.Integral()
        #qcd_r.Add( top_2b, -1)

        top_r = histos[r]["top"].Clone("top__"+r)
        if nbtag_top_shape =="3":
            temp_scaler = top_r.Integral() / histos[r_3b]["top"].Integral()
            top_r = histos[r_3b]["top"].Clone("top__"+r)
            top_r.Scale( temp_scaler )
        top_int = top_r.Integral()


        mu_qcd = bkgFitResults["muqcd"][ir]
        top_scale = (bkgFitResults["topscale"][0] if use_one_top_nuis else bkgFitResults["topscale"][ir])
        
        qcd_r.Scale( mu_qcd )
        top_r.Scale( top_scale )



        ## Now do smoothing

        qcd_sm = smoothfit.smoothfit(qcd_r, fitFunction = "Exp", fitRange = (900, 2000), makePlots = True, outfileName="qcd_fit.root")
        top_sm = smoothfit.smoothfit(top_r, fitFunction = "Exp", fitRange = (850, 1200), makePlots = True, outfileName="top_fit.root")

        qcd_final = smoothfit.MakeSmoothHisto(qcd_r, qcd_sm["nom"])
        top_final = smoothfit.MakeSmoothHisto(top_r, top_sm["nom"])

        outfileStat.WriteTObject(qcd_final, "qcd_hh_nominal","Overwrite")
        outfileStat.WriteTObject(top_final, "top_hh_nominal","Overwrite")

        

        ### propagate correlated systematics from the smoothing procedure---> these "replace" the stat error on the bins  #############
        for ivar in range(len(qcd_sm["vars"])):
            qup = qcd_sm["vars"][ivar][0]
            qdw = qcd_sm["vars"][ivar][1]

            outfileStat.WriteTObject(smoothfit.MakeSmoothHisto(qcd_r, qup), "qcd_hh_smoothQ"+str(ivar)+"Up","Overwrite")
            outfileStat.WriteTObject(smoothfit.MakeSmoothHisto(qcd_r, qdw), "qcd_hh_smoothQ"+str(ivar)+"Down","Overwrite")

        for ivar in range(len(top_sm["vars"])):
            tup = top_sm["vars"][ivar][0]
            tdw = top_sm["vars"][ivar][1]

            outfileStat.WriteTObject(smoothfit.MakeSmoothHisto(top_r, tup), "top_hh_smoothT"+str(ivar)+"Up","Overwrite")
            outfileStat.WriteTObject(smoothfit.MakeSmoothHisto(top_r, tdw), "top_hh_smoothT"+str(ivar)+"Down","Overwrite")


            

        ### propagate correlated systematics from normalization fits for mu_qcd and top_scale ###############
        for ivar in range(len(pvars)):
            for iUD in range(2):
                mu_qcd_var = pvars[ivar][iUD][ir]
                top_scale_var = pvars[ivar][iUD][n_channels + (0 if use_one_top_nuis else ir) ]

                qvar = qcd_r.Clone("qvar")
                qvar.Scale( mu_qcd_var * qcd_int / qvar.Integral() )

                tvar = top_r.Clone("tvar")
                tvar.Scale( top_scale_var * top_int / tvar.Integral() )

                ## Now do smoothing

                qvar_sm = smoothfit.smoothfit(qvar, fitFunction = "Exp", fitRange = (900, 2000), makePlots = False)
                tvar_sm = smoothfit.smoothfit(tvar, fitFunction = "Exp", fitRange = (850, 1200), makePlots = False)
    
                qvar_final = smoothfit.MakeSmoothHisto(qvar, qvar_sm["nom"])
                tvar_final = smoothfit.MakeSmoothHisto(tvar, tvar_sm["nom"])

                UpDw = ("Up" if iUD ==0 else "Down")
                outfileStat.WriteTObject(qvar_final, "qcd_hh_normY"+str(ivar)+UpDw,"Overwrite")
                outfileStat.WriteTObject(tvar_final, "top_hh_normY"+str(ivar)+UpDw,"Overwrite")
            
            




        
        
        outfileStat.Close()
        

        

        ## if False:
        ##     pred_final = qcd_final.Clone("pred_final__"+r)
        ##     pred_final.Add( top_final )


        ##     func1 = qcd_sm["nom"]
        ##     func2 = top_sm["nom"]

        ##     pred_sm = R.TF1("pred_sm", FuncSum, 900, 3000)

        ##     pred_sm.Draw("same")
        ##     top_sm["nom"].Draw("same")

        ##     pred_final_raw = qcd_r.Clone("qcd_final_raw__"+r)
        ##     pred_final_raw.Add(top_r)

        ##     outfile = R.TFile("outfile_"+r+".root","RECREATE")

        ##     c=R.TCanvas()
        ##     pred_final_raw.Draw("HIST")
        ##     top_r.SetLineColor(R.kBlack)
        ##     top_r.SetFillColor(R.kGreen)
        ##     top_r.Draw("sameHIST")

        ##     pred_sm.Draw("same")
        ##     top_sm["nom"].Draw("same")

        ##     c.Write()

        ##     c=R.TCanvas()

        ##     pred_final.Draw("HIST")

        ##     top_final.SetLineColor(R.kBlack)
        ##     top_final.SetFillColor(R.kGreen)

        ##     top_final.Draw("sameHIST")

        ##     c.Write()

        ##     outfile.Close()

    


    return


def FuncSum(x):
    return ( func1.Eval(x[0]) + func2.Eval(x[0]))



if __name__=="__main__":
    HistoAnalysis()
