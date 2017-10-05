import ROOT as R
import numpy as np
from copy import deepcopy
import sys, time, argparse, os, glob
#import smoothfit_Ultimate as smoothfit
import smoothfit_Ultimate as smoothfit
import BackgroundFit_Ultimate as BkgFit
import SystematicsTools_withSmoothing as SystToolsSmooth
import ExpModGaussSmoothingSystematics as EMGSmoothSyst
from HistoTools import HistLocationString as HistLocStr
from HistoTools import CheckAndGet
#R.gROOT.LoadMacro("AtlasStyle.C") 
#R.SetAtlasStyle()

R.gROOT.SetBatch(True)
func1 = None
func2 = None
# rebinFinal -- added by Qi. should be array object. Do the rebinning before writing into output files
# nbtag_top_shape_normFit --- what top shape to be used in NORMALIZATION FIT?
# nbtag_top_shape_SRPred --- what top shape to be used in SR prediction?

#define functions
def options():
    parser = argparse.ArgumentParser()
    parser.add_argument("--hist", default="mHH_l")
    parser.add_argument("--massregion", default="SR")
    return parser.parse_args()

def HistoAnalysis(datafileName="/afs/cern.ch/user/b/btong/work/bbbb/MoriondAnalysis/Output/Moriond_bkg_9/data_test/hist-MiniNTuple.root",
                  topfileName="/afs/cern.ch/user/b/btong/work/bbbb/MoriondAnalysis/Output/Moriond_bkg_9/ttbar_comb_test/hist-MiniNTuple.root",
                  zjetfileName="/afs/cern.ch/user/b/btong/work/bbbb/MoriondAnalysis/Output/Moriond_bkg_9/zjets_test/hist-MiniNTuple.root",
                  distributionName= "mHH_l",
                  n_trkjet  = ["4","3","2"],
                  n_btag    = ["4","3","2"],
                  btag_WP   = "70",
                  NRebin    = 10,
                  use_one_top_nuis = False,
                  use_scale_top_0b = False,
                  nbtag_top_shape_SRPred_for4b = "22",
                  rebinFinal = None,
                  smoothing_func = "MJ8",
                  top_smoothing_func = "MJ8",
                  inputFitResult = None,
                  inputQCDSyst_Dict = None,
                  doSmoothing = True,
                  addSmoothErrorBin = False,
                  qcdSmoothRange = (1200, 3000), #(1200, 3000),
                  topSmoothRange = (1200, 3000), #(1200, 3000),
                  isSystematicVariation = False,
                  verbose = False,
                  makeOutputFiles = True,
                  MassRegionName = "SR",
                  do_variable_rebin = False
                  ):
    ##### Parse Inputs ############################################
    fitzjets    = True
    dist_name   = distributionName
    print "the chosen hist is: ", dist_name
    
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

    nbtag_top_shape_for4b = nbtag_top_shape_SRPred_for4b

    useOneTopNuis = use_one_top_nuis

    scaleTop0b = use_scale_top_0b

    n_channels = num_trkjet.shape[0]

    regions = [ num_trkjet[i]+num_btag[i] for i in range(n_channels) ]

    ##for outputing
    isMhhDistribution = (distributionName=="mHH_l" or distributionName=="mHH_pole")
    do_smoothing      = (doSmoothing if isMhhDistribution else False)   # qi
    ##################################################################
    ##### Storage Variables ############################################
    output_Dict = { }
    Nbkg_dict    = {  }
    Nbkg_SysList = {  }
    for ir in regions:
        Nbkg_dict[ir]    = { "qcd":0,  "top":0,  "zjet":0,  "bkg":0, "data":0 }
        Nbkg_SysList[ir] = { "qcd":[], "top":[], "zjet":[], "bkg":[], "data":[] }
    vartxt = ''
    ##################################################################    
    ##### Do Background Fits ############################################
    ##### This is just the same fitting procedure
    if inputFitResult == None:
        # bkgFitResults = BkgFit.BackgroundFit(datafileName=datafileName,
        #                                       topfileName=topfileName,
        #                                       zjetfileName=zjetfileName,
        #                                       distributionName = ["leadHCand_Mass"],
        #                                       whichFunc = "XhhBoosted",
        #                                       n_trkjet  = n_trkjet,
        #                                       n_btag    = n_btag,
        #                                       btag_WP   = btag_WP,
        #                                       NRebin    = 2,#NRebin, #this is reset to be fine binned
        #                                       use_one_top_nuis = use_one_top_nuis,
        #                                       makePlots = True,
        #                                       BKG_lst   = ["FourTag", "ThreeTag", "TwoTag_split"],
        #                                       BKG_dic   = {"FourTag":"NoTag_4Trk",  "ThreeTag":"NoTag_3Trk", "TwoTag_split":"NoTag_2Trk_split",  "TwoTag":"OneTag",  "OneTag":"NoTag"},
        #                                       fitzjets  = fitzjets)

        bkgest_lst    = ["FourTag", "ThreeTag", "TwoTag_split"] 
        bkgest_dict   = {"FourTag":"NoTag_4Trk",  "ThreeTag":"NoTag_3Trk", "TwoTag_split":"NoTag_2Trk_split",  "TwoTag":"OneTag",  "OneTag":"NoTag"}
        bkgFitResults = BkgFit.BackgroundFit(datafileName=datafileName,
                                              topfileName=topfileName,
                                              zjetfileName=zjetfileName,
            distributionName = ["leadHCand_Mass"], whichFunc = "XhhBoosted", NRebin=1, \
            BKG_lst=bkgest_lst, BKG_dic=bkgest_dict, use_one_top_nuis=False, fitzjets=True, a_ttbar=1.0)
        global best_attbar
        best_attbar = 1
        ##iterative fit method
        while (abs(bkgFitResults["muttbar"][2] - best_attbar) > 0.01):##use 2bs here for the normalization estiamte
            print "Refit!!!"
            best_attbar = bkgFitResults["muttbar"][2]
            bkgFitResults = BkgFit.BackgroundFit(datafileName=datafileName,
                                              topfileName=topfileName,
                                              zjetfileName=zjetfileName,
                distributionName = ["leadHCand_Mass"], whichFunc = "XhhBoosted", NRebin=1, \
                BKG_lst=bkgest_lst, BKG_dic=bkgest_dict, use_one_top_nuis=False, fitzjets=True, a_ttbar=best_attbar) #Weight_dic = weight_dict, 
    else:
        bkgFitResults = inputFitResult

    pvars = bkgFitResults["pvars"]
    output_Dict["fitResults"] = bkgFitResults
    ##################################################################
    ##### Get QCD Shape Systematics from CR  ##############################
    print "STEP: Get QCD Shape Systematics from CR"
    ##### This is smoothing the CR region distributions
    if MassRegionName == "SR" or "CR":
        # should only affect SR
        if inputQCDSyst_Dict == None and isMhhDistribution:  # qi
            QCDSyst_Dict =  SystToolsSmooth.QCDSystematics(datafileName=datafileName,
                                        topfileName=topfileName,
                                        zjetfileName=zjetfileName,
                                        distributionName = "mHH_l",   # this has been decided to fix on DiJetMass
                                        n_trkjet      = n_trkjet,
                                        n_btag        = n_btag,
                                        btag_WP       = btag_WP,
                                        mu_qcd_vals   = bkgFitResults["muqcd"],
                                        topscale_vals = bkgFitResults["muttbar"],
                                        NRebin        = 10, #10 just like SR; 5 is a finner bin
                                        smoothing_func = smoothing_func,
                                        SmoothRange    = qcdSmoothRange,# (100, 2500), #this is fixed...
                                        use_one_top_nuis = use_one_top_nuis,
                                        use_scale_top_0b = use_scale_top_0b,
                                        nbtag_top_shape_for4b = "33",
                                        makePlots = True,
                                        verbose = False,
                                        outfileNameBase="QCDSysfitSmooth.root") 
            
        elif inputQCDSyst_Dict != None:
            QCDSyst_Dict = inputQCDSyst_Dict

        else:
            QCDSyst_Dict = None
    else:
        QCDSyst_Dict = None

    output_Dict["QCDSystCR"] = QCDSyst_Dict

    ##################################################################
    ##### Get Signal Region Histograms ################################
    print "STEP: Get Signal Region Histograms"
    ##### This is loding input file histograms
    datafile  = R.TFile(datafileName,"READ")
    topfile   = R.TFile(topfileName,"READ")
    zjetfile  = ( R.TFile(zjetfileName,"READ") if fitzjets is True else None)
    histos    = {}
    
    # collect all histograms
    for r in ["44","33","22","40","30","20"]:
        folder_r = HistLocStr(dist_name, r[0], r[1], btag_WP, MassRegionName)  #folder( r[0], r[1], btag_WP)
        
        data_r   = datafile.Get(folder_r).Clone("data_"+r)
        data_r.SetDirectory(0)

        top_r    = topfile.Get(folder_r).Clone("top_"+r)
        top_r.SetDirectory(0)

        zjet_r   = CheckAndGet(zjetfile, folder_r, top_r).Clone("zjet_"+r)
        zjet_r.SetDirectory(0)

        #DO NOT clear the negative weight bins for ttbar
        #ClearNegBin(top_r)
        #print folder_r, top_r.Integral()

        if do_variable_rebin:
            data_r = smoothfit.VariableRebin(data_r,5,2000)
            top_r  = smoothfit.VariableRebin(top_r,5,2000)
            zjet_r = smoothfit.VariableRebin(zjet_r,5,2000)
        else:
            data_r.Rebin(n_rebin)
            top_r.Rebin(n_rebin)
            zjet_r.Rebin(n_rebin)

        histos[r] = {"data": data_r,  "top": top_r,  "zjet":zjet_r}

    datafile.Close()
    topfile.Close()
    if zjetfile != None:
        zjetfile.Close()

    ##################################################################
    ##### scaling and subtractions #################################
    print "STEP: scaling and subtractions"
    ##### This is loding input file histograms
    for ir in range(len(regions)):
        # print ir
        r = regions[ir]

        output_Dict[r] = {"qcd":{}, "ttbar":{}, "zjet":{}}
        
        if makeOutputFiles:
            cut_lst = {"44":"FourTag", "33":"ThreeTag", "22":"TwoTag_split"}
            outfileStat = R.TFile("outfile_boosted_"+cut_lst[r]+".root","RECREATE")
        
        r_0b = r[0]+"0"
        #r_3b = r[0]+"3"
        top_0b = histos[r_0b]["top"].Clone("top_0b__"+r)
        if scaleTop0b:
            top_0b.Scale( (bkgFitResults["muttbar"][0] if use_one_top_nuis else bkgFitResults["muttbar"][ir]) )

        zjet_0b = histos[r_0b]["zjet"].Clone("zjet_0b__"+r)

        qcd_r = histos[r_0b]["data"].Clone("qcd__"+r)
        print "zjet total:", zjet_0b.Integral(), "top total:", top_0b.Integral(), "qcd total:", qcd_r.Integral()
        qcd_r.Add( top_0b, -1 * best_attbar)
        qcd_r.Add( zjet_0b, -1)

        #clear the negative weight bins for qcd as well
        ClearNegBin(qcd_r)
        qcd_int = qcd_r.Integral()
        print "qcd total:", qcd_int, "here! 2", bkgFitResults["muqcd"][ir], best_attbar

        #print histos[r]["top"].GetName(), r
        top_r = histos[r]["top"].Clone("top__"+r)
        #print top_r.Integral(), "here! 0"
        if (nbtag_top_shape_for4b == "22") and (r == "44") and (MassRegionName == "SR"):   # the 3b top shape is only used during the SR prediction for 44 region
            temp_scaler = top_r.Integral() / histos["22"]["top"].Integral()
            top_r = histos["22"]["top"].Clone("top__"+r)
            top_r.Scale( temp_scaler )
        if (nbtag_top_shape_for4b == "22") and (r == "33") and (MassRegionName == "SR"):   # the 3b top shape is only used during the SR prediction for 44 region
            temp_scaler = top_r.Integral() / histos["22"]["top"].Integral()
            top_r = histos["22"]["top"].Clone("top__"+r)
            top_r.Scale( temp_scaler )
        top_int = top_r.Integral()
        #print top_r.Integral(), "here! 1"

        zjet_r = histos[r]["zjet"].Clone("zjet__"+r)

        mu_qcd = bkgFitResults["muqcd"][ir]
        top_scale = (bkgFitResults["muttbar"][0] if use_one_top_nuis else bkgFitResults["muttbar"][ir])
        
        qcd_r.Scale( mu_qcd )
        top_r.Scale( top_scale )
        print "top total:", top_r.Integral(), " ; qcd total:", qcd_r.Integral(), "here! 2"

        bkg_r = qcd_r.Clone("bkg__" + r)
        bkg_r.Add( top_r, 1)
        bkg_r.Add( zjet_r, 1)
        # store some numbers for the output table later
        e_qcd = R.Double(0.0)
        e_top = R.Double(0.0)
        e_bkg = R.Double(0.0)
        e_data = R.Double(0.0)
        Nbkg_dict[r]["qcd"] = qcd_r.IntegralAndError(0, qcd_r.GetNbinsX()+1, e_qcd)
        Nbkg_dict[r]["top"] = top_r.IntegralAndError(0, top_r.GetNbinsX()+1, e_top)
        Nbkg_dict[r]["bkg"] = bkg_r.IntegralAndError(0, bkg_r.GetNbinsX()+1, e_bkg)
        Nbkg_dict[r]["data"] = histos[r]["data"].IntegralAndError(0, histos[r]["data"].GetNbinsX()+1, e_data)

        Nbkg_SysList[r]["qcd"].append( float(e_qcd) )
        Nbkg_SysList[r]["top"].append( float(e_top) )
        Nbkg_SysList[r]["bkg"].append( float(e_bkg) )   # Qi Question; Tony Question as well...
        Nbkg_SysList[r]["data"].append( float(e_data) )
        
        ## Now do smoothing ###########################################################################################
        print "start smoothing: ", ir
        if do_smoothing:
            ##for qcd, trick and add zjet into the total normalization; consistent as smoothing part
            ##renormalize to Zjet normalization
            qcd_r.Scale((zjet_r.Integral() + qcd_r.Integral())/qcd_r.Integral())
            int_pre = qcd_r.Integral()
            qcd_sm = smoothfit.smoothfit(qcd_r, fitFunction = smoothing_func, fitRange = qcdSmoothRange, makePlots = False, verbose = False, outfileName="qcd_smoothfit_"+r+".root", maxPlotRange=4000)
            top_sm = smoothfit.smoothfit(top_r, fitFunction = top_smoothing_func, fitRange = topSmoothRange, makePlots = False, verbose = False, outfileName="top_smoothfit_"+r+".root", maxPlotRange=4000)
            
            #print "top total:", top_r.Integral(), " ; qcd total:", qcd_r.Integral(), "here! 2.5"
            if addSmoothErrorBin:
                qcd_final = smoothfit.MakeSmoothHistoWithError(qcd_r, qcd_sm)
                top_final = smoothfit.MakeSmoothHistoWithError(top_r, top_sm)
            else:
                qcd_final = smoothfit.MakeSmoothHisto(qcd_r, qcd_sm["nom"])
                top_final = smoothfit.MakeSmoothHisto(top_r, top_sm["nom"])
            qcd_final.SetNameTitle("qcd_hh_"+r+"__clone",   "qcd_hh_"+r+"__clone")
            top_final.SetNameTitle("ttbar_hh_"+r+"__clone", "ttbar_hh_"+r+"__clone")
        else:
            qcd_final = qcd_r.Clone("qcd_hh_"+r+"__clone")
            top_final = top_r.Clone("ttbar_hh_"+r+"__clone")
        
        int_aft = qcd_final.Integral()
        if int_aft > 0:
            qcd_final.Scale(int_pre/int_aft) ##fix normalization hard way

        print "after smoothing: ", qcd_final.Integral()
        #print "top total:", top_final.Integral(), " ; qcd total:", qcd_final.Integral(), "here! 3"
        zjet_final = zjet_r.Clone("zjet_hh_"+r+"__clone")
        
        if rebinFinal is not None:
            qcd_final = qcd_final.Rebin(len(rebinFinal)-1, qcd_final.GetName()+"_rebinFinal", rebinFinal)
            top_final = top_final.Rebin(len(rebinFinal)-1, top_final.GetName()+"_rebinFinal", rebinFinal)
            zjet_final = zjet_final.Rebin(len(rebinFinal)-1, zjet_final.GetName()+"_rebinFinal", rebinFinal)

        if makeOutputFiles:
            outfileStat.WriteTObject(qcd_final, "qcd_hh","Overwrite")
            outfileStat.WriteTObject(top_final, "ttbar_hh","Overwrite")
            outfileStat.WriteTObject(zjet_final, "zjet_hh","Overwrite")

        qcd_final.SetDirectory(0)
        top_final.SetDirectory(0)
        zjet_final.SetDirectory(0)
        output_Dict[r]["qcd"]["nom"] = qcd_final
        output_Dict[r]["ttbar"]["nom"] = top_final
        output_Dict[r]["zjet"]["nom"] = zjet_final

        # for systematics, don't need anything after this in loop
        if isSystematicVariation:
            continue
        ##################################################################################################################################
        ### propagate correlated systematics from the smoothing procedure---> these "replace" the stat error on the bins     #############
        ##################################################################################################################################
        ##### This is adding smoothing systematics
        if do_smoothing:
            ## qcd smoothing variations#################################################################
            if not addSmoothErrorBin:
                for ivar in range(len(qcd_sm["vars"])):
                    qup = qcd_sm["vars"][ivar][0]
                    qdw = qcd_sm["vars"][ivar][1]

                    qcd_r_qup = smoothfit.MakeSmoothHisto(qcd_r, qup, keepNorm=False)
                    qcd_r_qdw = smoothfit.MakeSmoothHisto(qcd_r, qdw, keepNorm=False)

                    qcd_r_qup.SetNameTitle("qcd_hh_"+r+"_smoothQ"+str(ivar)+"up__clone", "qcd_hh_"+r+"_smoothQ"+str(ivar)+"up__clone")
                    qcd_r_qdw.SetNameTitle("qcd_hh_"+r+"_smoothQ"+str(ivar)+"down__clone", "qcd_hh_"+r+"_smoothQ"+str(ivar)+"down__clone")

                    if rebinFinal is not None:
                        qcd_r_qup = qcd_r_qup.Rebin(len(rebinFinal)-1, qcd_r_qup.GetName()+"_rebinFinal", rebinFinal)
                        qcd_r_qdw = qcd_r_qdw.Rebin(len(rebinFinal)-1, qcd_r_qdw.GetName()+"_rebinFinal", rebinFinal)

                    if makeOutputFiles:
                        outfileStat.WriteTObject(qcd_r_qup, "qcd_hh_smoothQ"+str(ivar)+"up","Overwrite")
                        outfileStat.WriteTObject(qcd_r_qdw, "qcd_hh_smoothQ"+str(ivar)+"down","Overwrite")

                    qcd_r_qup.SetDirectory(0)
                    qcd_r_qdw.SetDirectory(0)
                    output_Dict[r]["qcd"]["smoothQ"+str(ivar)+"up"] = qcd_r_qup
                    output_Dict[r]["qcd"]["smoothQ"+str(ivar)+"down"] = qcd_r_qdw
   
            ## qcd smoothing function variations #################################################################
            if smoothing_func == "ExpModGauss":
                smoothFuncCompSyst = EMGSmoothSyst.smoothFuncCompare(qcd_r, fitFunction = smoothing_func,
                                                                     fitRange = qcdSmoothRange, funcCompareRange=(900, qcdSmoothRange[1]),
                                                                     makePlots = True, verbose = False, outfileName="EMGSmoothFuncCompare_"+r+".root", plotExtra=False)  # Qi
            else:
                smoothFuncCompSyst = smoothfit.smoothFuncCompare(qcd_r, fitFunction = smoothing_func, fitRange = qcdSmoothRange,            # qi
                                                                 makePlots = True, verbose = False, outfileName="smoothFuncCompare_"+r+".root", plotExtra=False)  # Qi
                
            qcd_r_func_up = smoothFuncCompSyst["up"]
            qcd_r_func_dw = smoothFuncCompSyst["dw"]
            qcd_r_func_up_super = smoothFuncCompSyst["up_super"]
            qcd_r_func_dw_super = smoothFuncCompSyst["dw_super"]

            if rebinFinal is not None:
                qcd_r_func_up = qcd_r_func_up.Rebin(len(rebinFinal)-1, qcd_r_func_up.GetName()+"_rebinFinal", rebinFinal)
                qcd_r_func_dw = qcd_r_func_dw.Rebin(len(rebinFinal)-1, qcd_r_func_dw.GetName()+"_rebinFinal", rebinFinal)

            if makeOutputFiles:
                outfileStat.WriteTObject(qcd_r_func_up, "qcd_hh_smoothFuncup","Overwrite")
                outfileStat.WriteTObject(qcd_r_func_dw, "qcd_hh_smoothFuncdown","Overwrite")
                
                outfileStat.WriteTObject(qcd_r_func_up_super, "qcd_hh_smoothFuncSuperup","Overwrite")
                outfileStat.WriteTObject(qcd_r_func_dw_super, "qcd_hh_smoothFuncSuperdown","Overwrite")

            # treat negative bin
            ClearNegBin(qcd_r_func_up)
            ClearNegBin(qcd_r_func_dw)
            ClearNegBin(qcd_r_func_up_super)
            ClearNegBin(qcd_r_func_dw_super)

            qcd_r_func_up.SetDirectory(0)
            qcd_r_func_dw.SetDirectory(0)
            output_Dict[r]["qcd"]["smoothFuncup"] = qcd_r_func_up
            output_Dict[r]["qcd"]["smoothFuncdown"] = qcd_r_func_dw

            qcd_r_func_up_super.SetDirectory(0)
            qcd_r_func_dw_super.SetDirectory(0)
            output_Dict[r]["qcd"]["smoothFuncup_super"] = qcd_r_func_up_super
            output_Dict[r]["qcd"]["smoothFuncdown_super"] = qcd_r_func_dw_super
            
            stepped_min_vals = []
            stepped_max_vals = []
            stepped_fitting = True
            if stepped_fitting == True:
                starting_bin = qcd_r_qup.FindBin(qcdSmoothRange[0])
                stepped_min_vals.append(str(qcdSmoothRange[0]))
                for step in range(0, 3):
                    current_starting_bin = starting_bin + step*1
                    current_starting_mass = qcd_r_qup.GetBinCenter(current_starting_bin)
                    stepped_min_vals.append(str(int(current_starting_mass)))
                stepped_min_vals = ["1200", "1300", "1400"]
                stepped_max_vals = ["2800", "3000", "3200"]
            else:
                stepped_max_vals = ["1850","2000","2250","2500"]
                stepped_min_vals = [str(qcdSmoothRange[0]),"1300","1400"]    
            print "MAX AND MIN ARE:"
            print stepped_max_vals
            print stepped_min_vals

            smoothfit.smoothFuncRangeCompare(qcd_r, fitFunction = smoothing_func, fitRange = qcdSmoothRange, fitMaxVals = stepped_max_vals, fitMinVals=stepped_min_vals,
                                            makePlots = True, plotExtra = False, verbose = False, outfileName="smoothFuncRangeCompare_"+r+".root")   # Qi

            
            ## ttbar smoothing variations##############################################################################
            if not addSmoothErrorBin:
                for ivar in range(len(top_sm["vars"])):
                    tup = top_sm["vars"][ivar][0]
                    tdw = top_sm["vars"][ivar][1]

                    top_r_tup = smoothfit.MakeSmoothHisto(top_r, tup)
                    top_r_tdw = smoothfit.MakeSmoothHisto(top_r, tdw)

                    top_r_tup.SetNameTitle("ttbar_hh_"+r+"_smoothQ"+str(ivar)+"up__clone",   "ttbar_hh_"+r+"_smoothQ"+str(ivar)+"up__clone")
                    top_r_tdw.SetNameTitle("ttbar_hh_"+r+"_smoothQ"+str(ivar)+"down__clone", "ttbar_hh_"+r+"_smoothQ"+str(ivar)+"down__clone")

                    if rebinFinal is not None:
                        top_r_tup = top_r_tup.Rebin(len(rebinFinal)-1, top_r_tup.GetName()+"_rebinFinal", rebinFinal)
                        top_r_tdw = top_r_tdw.Rebin(len(rebinFinal)-1, top_r_tdw.GetName()+"_rebinFinal", rebinFinal)

                    if makeOutputFiles:
                        outfileStat.WriteTObject(top_r_tup, "ttbar_hh_smoothQ"+str(ivar)+"up","Overwrite")
                        outfileStat.WriteTObject(top_r_tdw, "ttbar_hh_smoothQ"+str(ivar)+"down","Overwrite")

                    ClearNegBin(top_r_tup)
                    ClearNegBin(top_r_tdw)
                    top_r_tup.SetDirectory(0)
                    top_r_tdw.SetDirectory(0)
                    output_Dict[r]["ttbar"]["smoothQ"+str(ivar)+"up"] = top_r_tup
                    output_Dict[r]["ttbar"]["smoothQ"+str(ivar)+"down"] = top_r_tdw

            ## ttbar smootthing function variations##############################################################################
            if False:
                topsmoothFuncCompSyst = smoothfit.smoothFuncCompare(top_r, fitFunction = smoothing_func, fitRange = topSmoothRange,            # qi
                                                                 makePlots = True, verbose = False, outfileName="topsmoothFuncCompare_"+r+".root", plotExtra=False)  # Qi


                top_r_func_up = topsmoothFuncCompSyst["up"]
                top_r_func_dw = topsmoothFuncCompSyst["dw"]

                if rebinFinal is not None:
                    top_r_func_up = top_r_func_up.Rebin(len(rebinFinal)-1, top_r_func_up.GetName()+"_rebinFinal", rebinFinal)
                    top_r_func_dw = top_r_func_dw.Rebin(len(rebinFinal)-1, top_r_func_dw.GetName()+"_rebinFinal", rebinFinal)

                if makeOutputFiles:
                    outfileStat.WriteTObject(top_r_func_up, "ttbar_hh_smoothFuncup","Overwrite")
                    outfileStat.WriteTObject(top_r_func_dw, "ttbar_hh_smoothFuncdown","Overwrite")
                    

                # treat negative bin
                ClearNegBin(top_r_func_up)
                ClearNegBin(top_r_func_dw)

                top_r_func_up.SetDirectory(0)
                top_r_func_dw.SetDirectory(0)
                output_Dict[r]["ttbar"]["smoothFuncup"] = top_r_func_up
                output_Dict[r]["ttbar"]["smoothFuncdown"] = top_r_func_dw
                
                stepped_min_vals = []
                stepped_max_vals = []
                stepped_fitting = True
                if stepped_fitting == True:
                    starting_bin = top_r_tup.FindBin(topSmoothRange[0])
                    stepped_min_vals.append(str(topSmoothRange[0]))
                    for step in range(0, 3):
                        current_starting_bin = starting_bin + step*1
                        current_starting_mass = top_r_tup.GetBinCenter(current_starting_bin)
                        stepped_min_vals.append(str(int(current_starting_mass)))
                    stepped_min_vals = ["1200", "1300", "1400"]
                    stepped_max_vals = ["2800", "3000", "3200"]
                else:
                    stepped_max_vals = ["1850","2000","2250","2500"]
                    stepped_min_vals = [str(topSmoothRange[0]),"1300","1400"]    
                print "MAX AND MIN ARE:"
                print stepped_max_vals
                print stepped_min_vals

                smoothfit.smoothFuncRangeCompare(top_r, fitFunction = smoothing_func, fitRange = topSmoothRange, fitMaxVals = stepped_max_vals, fitMinVals=stepped_min_vals,
                                                makePlots = True, plotExtra = False, verbose = False, outfileName="topsmoothFuncRangeCompare_"+r+".root")   # Qi


        ########################################################################################################
        ### propagate correlated systematics from normalization fits for mu_qcd and top_scale    ###############
        ########################################################################################################
        ##### This is adding systematics from the fit
        #print pvars
        for ivar in range(len(pvars)):
            sys_qcd = []
            sys_top = []
            sys_bkg = []
            for iUD in range(2):
                upDw = ("up" if iUD ==0 else "down")

                mu_qcd_var = pvars[ivar][iUD][ir]
                top_scale_var = pvars[ivar][iUD][n_channels + (0 if use_one_top_nuis else ir) ]

                qvar = qcd_r.Clone("qvar")
                qvar.Scale( mu_qcd_var * qcd_int / qvar.Integral() )

                ## for ibin in range(1, qvar.GetNbinsX()+1):
                ##     if qvar.GetBinError(ibin) > qvar.GetBinContent(ibin):
                ##         qvar.SetBinError(ibin, qvar.GetBinContent(ibin))
                tvar = top_r.Clone("tvar")
                tvar.Scale( top_scale_var * top_int / tvar.Integral() )

                ### store some numbers for table
                sys_qcd.append( qvar.Integral() - Nbkg_dict[r]["qcd"] )
                sys_top.append( tvar.Integral() - Nbkg_dict[r]["top"] )
                sys_bkg.append( qvar.Integral() + tvar.Integral() - Nbkg_dict[r]["bkg"] )

                #vartxt = vartxt + str(r) + ' ' + str(ivar) + ' ' + str(iUD) + ' ' + str(qvar.Integral()) + ' ' + str(tvar.Integral()) + ' ' + str( (qvar.Integral() + tvar.Integral())) + '\n'
                ## Now do smoothing #######
                if do_smoothing:
                    qvar_sm = smoothfit.smoothfit(qvar, fitFunction = smoothing_func, fitRange = qcdSmoothRange, makePlots = False, verbose = verbose,
                                                  outfileName="qcd_smoothfit_"+r+"_Norm"+str(ivar)+str(iUD)+".root")
                    tvar_sm = smoothfit.smoothfit(tvar, fitFunction = top_smoothing_func, fitRange = topSmoothRange, makePlots = False, verbose = verbose,
                                                  outfileName="top_smoothfit_"+r+"_Norm"+str(ivar)+str(iUD)+".root")

                    if addSmoothErrorBin:
                        qvar_final = smoothfit.MakeSmoothHistoWithError(qvar, qvar_sm)
                        tvar_final = smoothfit.MakeSmoothHistoWithError(tvar, tvar_sm)
                    else:
                        qvar_final = smoothfit.MakeSmoothHisto(qvar, qvar_sm["nom"])
                        tvar_final = smoothfit.MakeSmoothHisto(tvar, tvar_sm["nom"])

                    qvar_final.SetNameTitle("qcd_hh_"+r+"_normY"+str(ivar)+upDw+"__clone",   "qcd_hh_"+r+"_normY"+str(ivar)+upDw+"__clone")
                    tvar_final.SetNameTitle("ttbar_hh_"+r+"_normY"+str(ivar)+upDw+"__clone", "ttbar_hh_"+r+"_normY"+str(ivar)+upDw+"__clone")

                else:
                    qvar_final = qvar.Clone("qcd_hh_"+r+"_normY"+str(ivar)+upDw+"__clone")
                    tvar_final = tvar.Clone("ttbar_hh_"+r+"_normY"+str(ivar)+upDw+"__clone")

                if rebinFinal is not None:
                    qvar_final = qvar_final.Rebin(len(rebinFinal)-1, qvar_final.GetName()+"_rebinFinal", rebinFinal)
                    tvar_final = tvar_final.Rebin(len(rebinFinal)-1, tvar_final.GetName()+"_rebinFinal", rebinFinal)

                if makeOutputFiles:
                    outfileStat.WriteTObject(qvar_final, "qcd_hh_normY"+str(ivar)+upDw,"Overwrite")
                    outfileStat.WriteTObject(tvar_final, "ttbar_hh_normY"+str(ivar)+upDw,"Overwrite")

                qvar_final.SetDirectory(0)
                tvar_final.SetDirectory(0)
                output_Dict[r]["qcd"]["normY"+str(ivar)+upDw] = qvar_final
                output_Dict[r]["ttbar"]["normY"+str(ivar)+upDw] = tvar_final

                
            # store some numbers for table later
            e_qcd_i = np.max( np.abs(sys_qcd) )
            e_top_i = np.max( np.abs(sys_top) )
            e_bkg_i = np.max( np.abs(sys_bkg) )
            Nbkg_SysList[r]["qcd"].append( e_qcd_i )
            Nbkg_SysList[r]["top"].append( e_top_i )
            Nbkg_SysList[r]["bkg"].append( e_bkg_i )


        ########################################################################################################
        ####### QCD Shape and Norm estimated from CR            ################################################
        ########################################################################################################
        ##### This is adding systematics from the CR
        if QCDSyst_Dict!=None and isMhhDistribution:  # qi
            original_norm = qcd_final.Integral(0, 4000) ##restricted range
            qvar_shape_up = qcd_r.Clone("qvar_QCDshape_up")
            qvar_shape_dw = qcd_r.Clone("qvar_QCDshape_dw")
            tvar_shape_up_final = top_final.Clone("tvar_QCDshape_up")
            tvar_shape_dw_final = top_final.Clone("tvar_QCDshape_dw")

            ClearNegBin(qvar_shape_up)
            ClearNegBin(qvar_shape_dw)
            ClearNegBin(tvar_shape_up_final)
            ClearNegBin(tvar_shape_dw_final)
        
            ## Now do smoothing
            if do_smoothing:
                qvar_shape_up_sm = smoothfit.smoothfit(qvar_shape_up, fitFunction = smoothing_func, fitRange = qcdSmoothRange, makePlots = False, verbose = verbose,
                                                        outfileName="qcd_smoothfit_"+r+"_QCDShapeup.root")
                qvar_shape_dw_sm = smoothfit.smoothfit(qvar_shape_dw, fitFunction = smoothing_func, fitRange = qcdSmoothRange, makePlots = False, verbose = verbose,
                                                        outfileName="qcd_smoothfit_"+r+"_QCDShapedown.root")

                if addSmoothErrorBin:
                    qvar_shape_up_final = smoothfit.MakeSmoothHistoWithError(qvar_shape_up, qvar_shape_up_sm)
                    qvar_shape_dw_final = smoothfit.MakeSmoothHistoWithError(qvar_shape_dw, qvar_shape_dw_sm)
                else:
                    qvar_shape_up_final = smoothfit.MakeSmoothHisto(qvar_shape_up, qvar_shape_up_sm["nom"], keepNorm=True)
                    qvar_shape_dw_final = smoothfit.MakeSmoothHisto(qvar_shape_dw, qvar_shape_dw_sm["nom"], keepNorm=True)

                qvar_shape_up_final.Multiply( QCDSyst_Dict["Shape_"+r] )
                qvar_shape_dw_final.Divide( QCDSyst_Dict["Shape_"+r] )
                ##sacale the top as well
                tvar_shape_up_final.Multiply( QCDSyst_Dict["Shape_"+r] )
                tvar_shape_dw_final.Divide( QCDSyst_Dict["Shape_"+r] )

                # for i in range(1, qvar_shape_up_final.GetNbinsX()):
                #     print r, i, qvar_shape_up_final.GetBinContent(i), qvar_shape_up_final.GetBinCenter(i), qcd_r.GetBinContent(i), QCDSyst_Dict["Shape_"+r].Integral(qcd_r.GetBinLowEdge(i), qcd_r.GetBinLowEdge(i+1))

                print "HERE!!!!!!!!", QCDSyst_Dict["Shape_"+r], qvar_shape_up_final.GetName(), original_norm, qvar_shape_up_final.Integral(0, 4000), qvar_shape_dw_final.Integral(0, 4000)

                qvar_shape_up_final.SetNameTitle("qcd_hh_"+r+"_QCDShapeCRup__clone",     "qcd_hh_"+r+"_QCDShapeCRup__clone")
                qvar_shape_dw_final.SetNameTitle("qcd_hh_"+r+"_QCDShapeCRdown__clone",   "qcd_hh_"+r+"_QCDShapeCRdown__clone")
            else:
                qvar_shape_up_final = qvar_shape_up.Clone("qcd_hh_"+r+"_QCDShapeCRup__clone")
                qvar_shape_dw_final = qvar_shape_dw.Clone("qcd_hh_"+r+"_QCDShapeCRdown__clone")

            #make sure normalization is correct; not necessary for now
            #qvar_shape_up_final.Scale( original_norm/qvar_shape_up_final.Integral(0, 4000) )
            #qvar_shape_dw_final.Scale( original_norm/qvar_shape_dw_final.Integral(0, 4000) )
            if rebinFinal is not None:
                qvar_shape_up_final = qvar_shape_up_final.Rebin(len(rebinFinal)-1, qvar_shape_up_final.GetName()+"_rebinFinal", rebinFinal)
                qvar_shape_dw_final = qvar_shape_dw_final.Rebin(len(rebinFinal)-1, qvar_shape_dw_final.GetName()+"_rebinFinal", rebinFinal)

            # treat increasing bin, flat it out
            for i in range(1, qvar_shape_up_final.GetNbinsX()):
                if qvar_shape_up_final.GetBinLowEdge(i) > qcdSmoothRange[0]:
                   if qvar_shape_up_final.GetBinContent(i) < qvar_shape_up_final.GetBinContent(i + 1):
                        qvar_shape_up_final.SetBinContent(i + 1, qvar_shape_up_final.GetBinContent(i))
                        qvar_shape_up_final.SetBinError(i + 1, qvar_shape_up_final.GetBinError(i))
            for i in range(1, qvar_shape_dw_final.GetNbinsX()):
                if qvar_shape_dw_final.GetBinLowEdge(i) > qcdSmoothRange[0]:
                   if qvar_shape_dw_final.GetBinContent(i) < qvar_shape_dw_final.GetBinContent(i + 1):
                        qvar_shape_dw_final.SetBinContent(i + 1, qvar_shape_dw_final.GetBinContent(i))
                        qvar_shape_dw_final.SetBinError(i + 1, qvar_shape_dw_final.GetBinError(i))
            for i in range(1, tvar_shape_dw_final.GetNbinsX()):
                if tvar_shape_dw_final.GetBinLowEdge(i) > qcdSmoothRange[0]:
                   if tvar_shape_dw_final.GetBinContent(i) < tvar_shape_dw_final.GetBinContent(i + 1):
                        tvar_shape_dw_final.SetBinContent(i + 1, tvar_shape_dw_final.GetBinContent(i))
                        tvar_shape_dw_final.SetBinError(i + 1, tvar_shape_dw_final.GetBinError(i))
            for i in range(1, tvar_shape_dw_final.GetNbinsX()):
                if tvar_shape_dw_final.GetBinLowEdge(i) > qcdSmoothRange[0]:
                   if tvar_shape_dw_final.GetBinContent(i) < tvar_shape_dw_final.GetBinContent(i + 1):
                        tvar_shape_dw_final.SetBinContent(i + 1, tvar_shape_dw_final.GetBinContent(i))
                        tvar_shape_dw_final.SetBinError(i + 1, tvar_shape_dw_final.GetBinError(i))

            ###split the QCD shape systematic into two parts
            qvar_shape_up_low  = qvar_shape_up_final.Clone(qvar_shape_up_final.GetName() + "_low")
            qvar_shape_up_low.SetDirectory(0)
            qvar_shape_up_high = qvar_shape_up_final.Clone(qvar_shape_up_final.GetName() + "_high")
            qvar_shape_up_high.SetDirectory(0)
            qvar_shape_dw_low  = qvar_shape_dw_final.Clone(qvar_shape_dw_final.GetName() + "_low")
            qvar_shape_dw_low.SetDirectory(0)
            qvar_shape_dw_high = qvar_shape_dw_final.Clone(qvar_shape_dw_final.GetName() + "_high")
            qvar_shape_dw_high.SetDirectory(0)
            tvar_shape_up_low  = tvar_shape_up_final.Clone(tvar_shape_up_final.GetName() + "_low")
            tvar_shape_up_low.SetDirectory(0)
            tvar_shape_up_high = tvar_shape_up_final.Clone(tvar_shape_up_final.GetName() + "_high")
            tvar_shape_up_high.SetDirectory(0)
            tvar_shape_dw_low  = tvar_shape_dw_final.Clone(tvar_shape_dw_final.GetName() + "_low")
            tvar_shape_dw_low.SetDirectory(0)
            tvar_shape_dw_high = tvar_shape_dw_final.Clone(tvar_shape_dw_final.GetName() + "_high")
            tvar_shape_dw_high.SetDirectory(0)
            
            #print "FUCK THIS", qcd_r.GetNbinsX(), qvar_shape_up_final.GetNbinsX(), qvar_shape_up_high.GetNbinsX()
            for k in range(1, qcd_final.GetNbinsX()):
                if qcd_final.GetBinLowEdge(k) >= 2000:
                    qvar_shape_up_low.SetBinContent(k, qcd_final.GetBinContent(k))
                    qvar_shape_up_low.SetBinError(k, 0)
                    qvar_shape_dw_low.SetBinContent(k, qcd_final.GetBinContent(k))
                    qvar_shape_dw_low.SetBinError(k, 0)
                else:
                    qvar_shape_up_high.SetBinContent(k, qcd_final.GetBinContent(k))
                    qvar_shape_up_high.SetBinError(k, 0)
                    qvar_shape_dw_high.SetBinContent(k, qcd_final.GetBinContent(k))
                    qvar_shape_dw_high.SetBinError(k, 0)
            for k in range(1, top_final.GetNbinsX()):
                if top_final.GetBinLowEdge(k) >= 2000:
                    tvar_shape_up_low.SetBinContent(k, top_final.GetBinContent(k))
                    tvar_shape_up_low.SetBinError(k, 0)
                    tvar_shape_dw_low.SetBinContent(k, top_final.GetBinContent(k))
                    tvar_shape_dw_low.SetBinError(k, 0)
                else:
                    tvar_shape_up_high.SetBinContent(k, top_final.GetBinContent(k))
                    tvar_shape_up_high.SetBinError(k, 0)
                    tvar_shape_dw_high.SetBinContent(k, top_final.GetBinContent(k))
                    tvar_shape_dw_high.SetBinError(k, 0)

            if makeOutputFiles:
                outfileStat.WriteTObject(qvar_shape_up_final, "qcd_hh_QCDShapeCRup")
                outfileStat.WriteTObject(qvar_shape_dw_final, "qcd_hh_QCDShapeCRdown")
                outfileStat.WriteTObject(qvar_shape_up_low,   "qcd_hh_QCDShapeCRLowup")
                outfileStat.WriteTObject(qvar_shape_dw_low,   "qcd_hh_QCDShapeCRLowdown")
                outfileStat.WriteTObject(qvar_shape_up_high,  "qcd_hh_QCDShapeCRHighup")
                outfileStat.WriteTObject(qvar_shape_dw_high,  "qcd_hh_QCDShapeCRHighdown")
                outfileStat.WriteTObject(tvar_shape_up_final, "top_hh_QCDShapeCRup")
                outfileStat.WriteTObject(tvar_shape_dw_final, "top_hh_QCDShapeCRdown")
                outfileStat.WriteTObject(tvar_shape_up_low,   "top_hh_QCDShapeCRLowup")
                outfileStat.WriteTObject(tvar_shape_dw_low,   "top_hh_QCDShapeCRLowdown")
                outfileStat.WriteTObject(tvar_shape_up_high,  "top_hh_QCDShapeCRHighup")
                outfileStat.WriteTObject(tvar_shape_dw_high,  "top_hh_QCDShapeCRHighdown")

            qvar_shape_up_final.SetDirectory(0)
            qvar_shape_dw_final.SetDirectory(0)
            qvar_shape_up_low.SetDirectory(0)
            qvar_shape_dw_low.SetDirectory(0)
            qvar_shape_up_high.SetDirectory(0)
            qvar_shape_dw_high.SetDirectory(0)

            output_Dict[r]["qcd"]["QCDShapeCRup"]      = qvar_shape_up_final
            output_Dict[r]["qcd"]["QCDShapeCRdown"]    = qvar_shape_dw_final
            output_Dict[r]["qcd"]["QCDShapeCRLowup"]   = qvar_shape_up_low
            output_Dict[r]["qcd"]["QCDShapeCRLowdown"] = qvar_shape_dw_low
            output_Dict[r]["qcd"]["QCDShapeCRHighup"]  = qvar_shape_up_high
            output_Dict[r]["qcd"]["QCDShapeCRHighdown"]= qvar_shape_dw_high

            tvar_shape_up_final.SetDirectory(0)
            tvar_shape_dw_final.SetDirectory(0)
            tvar_shape_up_low.SetDirectory(0)
            tvar_shape_dw_low.SetDirectory(0)
            tvar_shape_up_high.SetDirectory(0)
            tvar_shape_dw_high.SetDirectory(0)
            
            output_Dict[r]["ttbar"]["QCDShapeCRup"]      = tvar_shape_up_final
            output_Dict[r]["ttbar"]["QCDShapeCRdown"]    = tvar_shape_dw_final
            output_Dict[r]["ttbar"]["QCDShapeCRLowup"]   = tvar_shape_up_low
            output_Dict[r]["ttbar"]["QCDShapeCRLowdown"] = tvar_shape_dw_low
            output_Dict[r]["ttbar"]["QCDShapeCRHighup"]  = tvar_shape_up_high
            output_Dict[r]["ttbar"]["QCDShapeCRHighdown"]= tvar_shape_dw_high
        ###########################################################################################
        ### Norm comparison in CR      ############################################################
        ###########################################################################################
        if QCDSyst_Dict != None:
            
            qvar_normCR_up =  qcd_final.Clone("qcd_hh_"+r+"_QCDnormCRup__clone")
            qvar_normCR_up.Scale( 1.0 + QCDSyst_Dict["Scale_"+r] )

            qvar_normCR_dw =  qcd_final.Clone("qcd_hh_"+r+"_QCDnormCRdown__clone")
            qvar_normCR_dw.Scale( 1.0 - QCDSyst_Dict["Scale_"+r] )

            if rebinFinal is not None:
                qvar_normCR_up = qvar_normCR_up.Rebin(len(rebinFinal)-1, qvar_normCR_up.GetName()+"_rebinFinal", rebinFinal)
                qvar_normCR_dw = qvar_normCR_dw.Rebin(len(rebinFinal)-1, qvar_normCR_dw.GetName()+"_rebinFinal", rebinFinal)

            if makeOutputFiles:
                outfileStat.WriteTObject(qvar_normCR_up, "qcd_hh_QCDNormCRup")
                outfileStat.WriteTObject(qvar_normCR_dw, "qcd_hh_QCDNormCRdown")

            qvar_normCR_up.SetDirectory(0)
            qvar_normCR_dw.SetDirectory(0)
            output_Dict[r]["qcd"]["QCDNormCRup"] = qvar_normCR_up
            output_Dict[r]["qcd"]["QCDNormCRdown"] = qvar_normCR_dw


        #####################################################################################################################
        ### top shape systematics in 4b region, if using 3b shape ###########################################################
        #####################################################################################################################
        if r == "44" and nbtag_top_shape_SRPred_for4b == "22" and MassRegionName == "SR"  and isMhhDistribution:   # qi
            ## ttbarShapeSRSyst_Dict = SystTools.ttbarShapeSysSR(topfileName,
            ##                                                     distributionName,
            ##                                                     signal_region = "22",
            ##                                                     compare_region = "33",
            ##                                                     btag_WP     = btag_WP,
            ##                                                     makePlots = True,
            ##                                                     verbose = False,
            ##                                                     outfileNameBase="TopShapeSRSysfit.root")
            ttbarShapeSRSyst_Dict = SystToolsSmooth.ttbarShapeSysSR(topfileName,
                                                                distributionName,
                                                                signal_region  = "33",
                                                                compare_region = "22",
                                                                btag_WP        = btag_WP,
                                                                smoothing_func = top_smoothing_func,
                                                                SmoothRange    = topSmoothRange,# (100, 2500),
                                                                makePlots      = True,
                                                                verbose        = False,
                                                                outfileNameBase="TopShapeSRSysfitSmooth.root")

            tvar_shape_up = top_r.Clone("tvar_ttbarShapeSR_up")
            #tvar_shape_up.Multiply( ttbarShapeSRSyst_Dict["fup"] )

            tvar_shape_dw = top_r.Clone("tvar_ttbarShapeSR_dw")
            #tvar_shape_dw.Multiply( ttbarShapeSRSyst_Dict["fdw"] )

            ClearNegBin(tvar_shape_up)
            ClearNegBin(tvar_shape_dw)
                    
            tvar_shape_up.Scale( top_r.Integral() / tvar_shape_up.Integral() )
            tvar_shape_dw.Scale( top_r.Integral() / tvar_shape_dw.Integral() )


            ## Now do smoothing ##########################
            if do_smoothing:
                tvar_shape_up_sm = smoothfit.smoothfit(tvar_shape_up, fitFunction = top_smoothing_func, fitRange = topSmoothRange, makePlots = False, verbose = verbose,
                                                        outfileName="top_smoothfit_"+r+"_ttbarShapeSRup.root")

                tvar_shape_dw_sm = smoothfit.smoothfit(tvar_shape_dw, fitFunction = top_smoothing_func, fitRange = topSmoothRange, makePlots = False, verbose = verbose,
                                                        outfileName="top_smoothfit_"+r+"_ttbarShapeSRedown.root")

                if addSmoothErrorBin:
                    tvar_shape_up_final = smoothfit.MakeSmoothHistoWithError(tvar_shape_up, tvar_shape_up_sm)
                    tvar_shape_dw_final = smoothfit.MakeSmoothHistoWithError(tvar_shape_dw, tvar_shape_dw_sm)
                else:
                    tvar_shape_up_final = smoothfit.MakeSmoothHisto(tvar_shape_up, tvar_shape_up_sm["nom"])
                    tvar_shape_dw_final = smoothfit.MakeSmoothHisto(tvar_shape_dw, tvar_shape_dw_sm["nom"])

                tvar_shape_up_final.Multiply( ttbarShapeSRSyst_Dict["Shape"] )
                tvar_shape_dw_final.Divide( ttbarShapeSRSyst_Dict["Shape"] )

                tvar_shape_up_final.SetNameTitle("ttbar_hh_"+r+"_ttbarShapeSRup__clone",     "ttbar_hh_"+r+"_ttbarShapeSRup__clone")
                tvar_shape_dw_final.SetNameTitle("ttbar_hh_"+r+"_ttbarShapeSRdown__clone",   "ttbar_hh_"+r+"_ttbarShapeSRdown__clone")


            else:
                tvar_shape_up_final = tvar_shape_up.Clone("ttbar_hh_"+r+"_ttbarShapeSRup__clone")
                tvar_shape_dw_final = tvar_shape_dw.Clone("ttbar_hh_"+r+"_ttbarShapeSRdown__clone")


            if rebinFinal is not None:
                tvar_shape_up_final = tvar_shape_up_final.Rebin(len(rebinFinal)-1, tvar_shape_up_final.GetName()+"_rebinFinal", rebinFinal)
                tvar_shape_dw_final = tvar_shape_dw_final.Rebin(len(rebinFinal)-1, tvar_shape_dw_final.GetName()+"_rebinFinal", rebinFinal)

            if makeOutputFiles:
                outfileStat.WriteTObject(tvar_shape_up_final, "ttbar_hh_ttbarShapeSRup")
                outfileStat.WriteTObject(tvar_shape_dw_final, "ttbar_hh_ttbarShapeSRdown")

            tvar_shape_up_final.SetDirectory(0)
            tvar_shape_dw_final.SetDirectory(0)
            output_Dict[r]["ttbar"]["ttbarShapeSRup"] = tvar_shape_up_final
            output_Dict[r]["ttbar"]["ttbarShapeSRdown"] = tvar_shape_dw_final
        
        ### close outfiles, if used ###
        if makeOutputFiles:
            outfileStat.Close()

    ### Print tables ###
    PrintTable( Nbkg_dict, Nbkg_SysList, regions)
    #print vartxt

    output_Dict['regions'] = regions
    #print output_Dict
    return 
    #disable the dictionary output for now
    #return output_Dict

def ClearNegBin(hist):
    for ibin in range(0, hist.GetNbinsX()+1):
        if hist.GetBinContent(ibin) < 0:
            hist.SetBinContent(ibin, 0)
            hist.SetBinError(ibin, 0)
    return


def FuncSum(x):
    return ( func1.Eval(x[0]) + func2.Eval(x[0]))

def PrintTable( Nbkg_dict, Nbkg_SysList, Regions):
    print Nbkg_dict
    e_qcd_tot = {}
    e_top_tot = {}
    e_bkg_tot = {}

    for iR in Regions:
        e_qcd_tot[iR] = 0
        e_top_tot[iR] = 0
        e_bkg_tot[iR] = 0

        for ierr in range(len(Nbkg_SysList[iR]["qcd"])):
            e_qcd_tot[iR] = e_qcd_tot[iR] + Nbkg_SysList[iR]["qcd"][ierr]**2
            e_top_tot[iR] = e_top_tot[iR] + Nbkg_SysList[iR]["top"][ierr]**2
            e_bkg_tot[iR] = e_bkg_tot[iR] + Nbkg_SysList[iR]["bkg"][ierr]**2

        e_qcd_tot[iR] = np.sqrt(e_qcd_tot[iR])
        e_top_tot[iR] = np.sqrt(e_top_tot[iR])
        e_bkg_tot[iR] = np.sqrt(e_bkg_tot[iR])


    columnStructure = '| l |'
    for ic in range(len(Regions)):
        columnStructure = columnStructure + ' c |'

        
    outtext = ''
    outtext  = outtext + ' \n'
    outtext  = outtext + ' \n'
    outtext  = outtext + '\\begin{table}[htbp!] \n'
    outtext  = outtext + '\\begin{center} \n'
    outtext  = outtext + '\\begin{tabular}{' + columnStructure + ' } \n'
    outtext  = outtext + '\\hline \n'
    outtext  = outtext + ' Sample '
    for iR in Regions:
        outtext  = outtext + ' & ' + iR[1] + 'b SR Prediction '
    outtext  = outtext + ' \\\\ \n'
    
    outtext  = outtext + '\\hline \n'
    outtext  = outtext + '\\hline \n'
    
    outtext  = outtext + 'QCD '
    for iR in Regions:
        outtext  = outtext + ' & ' + str(float('%.3g' % Nbkg_dict[iR]["qcd"] )) + ' $\pm$ ' + str(float('%.3g' % e_qcd_tot[iR] ))
    outtext  = outtext + ' \\\\ \n'
    
    outtext  = outtext + '$ t \\bar{t}$ '
    for iR in Regions:
        outtext  = outtext + ' & ' + str(float('%.3g' % Nbkg_dict[iR]["top"] )) + ' $\pm$ ' + str(float('%.3g' % e_top_tot[iR] ))
    outtext  = outtext + ' \\\\ \n'

    
    outtext  = outtext + '\\hline \n'
    outtext  = outtext + 'Total '
    for iR in Regions:
        outtext  = outtext + ' & ' + str(float('%.3g' % Nbkg_dict[iR]["bkg"] )) + ' $\pm$ ' + str(float('%.3g' % e_bkg_tot[iR] ))
    outtext  = outtext + ' \\\\ \n'
    
    outtext  = outtext + '\\hline \n'
    outtext  = outtext + '\\hline \n'
    outtext  = outtext + 'Data '
    for iR in Regions:
        outtext  = outtext + ' & ' + str(float('%.3g' % Nbkg_dict[iR]["data"] )) + ' $\pm$ ' + str(float('%.3g' % np.sqrt(Nbkg_dict[iR]["data"]) ))
    outtext  = outtext + ' \\\\ \n'

    outtext  = outtext + '\\hline \n'
    outtext  = outtext + '\\end{tabular}  \n'
    outtext  = outtext + '\\caption{Background predictions in SR}  \n'
    outtext  = outtext + '\\label{tab:boosted-SR-yields-wsys}  \n'
    outtext  = outtext + '\\end{center}  \n'
    outtext  = outtext + '\\end{table}  \n'
    outtext  = outtext + '  \n'
    outtext  = outtext + '  \n'

    print outtext


if __name__=="__main__":
    start_time = time.time()
    global ops
    ops = options()
    HistoAnalysis(distributionName= ops.hist, MassRegionName=ops.massregion)
    print("--- %s seconds ---" % (time.time() - start_time))

