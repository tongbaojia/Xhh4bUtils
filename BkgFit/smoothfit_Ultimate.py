import ROOT as R
import os
import numpy as np
import scipy.special
from array import array
import sys

from GetEigenVariations import GetEigenVariations

import cPickle as pickle

import time

def smoothfit(histo, fitFunction = "Exp", fitRange = (900, 3000), outrange_start = None, makePlots = False, verbose = False, useLikelihood=False, outfileName="fit", ouutfilepath="", initpar=[]):
    npar = None
    func = None
    fitChoice = None
    colorlist = [R.kBlue, R.kGreen, R.kOrange, R.kMagenta, R.kCyan, R.kPink]

    fitName = "fit_"+outfileName+"_%s" % (time.time())  # need to make this name unique, otherwise it will give wrong answer when two consecutive fitting is performed
    
    if fitFunction == "Exp":
        npar = 2
        fitChoice = ExpoFunc
        func = R.TF1(fitName, fitChoice, fitRange[0], fitRange[1], npar)
        # func.SetParameters(0.006, 5.0)
        func.SetParameters(0.005, 5.0)

    elif fitFunction == "Dijet":
        npar = 3
        fitChoice = DijetFunc
        func = R.TF1(fitName, fitChoice, fitRange[0], fitRange[1], npar)
        func.SetParameters(60, 100, 10)
        if len(initpar) == 3:
            func.SetParameters(initpar[0], initpar[1], initpar[2])

    elif fitFunction == "MJ2":
        npar = 3
        fitChoice = MJ2Func
        func = R.TF1(fitName, fitChoice, fitRange[0], fitRange[1], npar)
        func.SetParameters(10, 10, 2)

    elif fitFunction == "MJ3":
        npar = 3
        fitChoice = MJ3Func
        func = R.TF1(fitName, fitChoice, fitRange[0], fitRange[1], npar)
        func.SetParameters(10, 10, 2)

    elif fitFunction == "MJ4":
        npar = 3
        fitChoice = MJ4Func
        func = R.TF1(fitName, fitChoice, fitRange[0], fitRange[1], npar)
        func.SetParameters(1, 10, 1)

    elif fitFunction == "MJ5":
        npar = 3
        fitChoice = MJ5Func
        func = R.TF1(fitName, fitChoice, fitRange[0], fitRange[1], npar)
        func.SetParameters(1, 10, 1)

    elif fitFunction == "MJ6":
        npar = 3
        fitChoice = MJ6Func
        func = R.TF1(fitName, fitChoice, fitRange[0], fitRange[1], npar)
        func.SetParameters(1, 10, 1)

    elif fitFunction == "MJ7":
        npar = 3
        fitChoice = MJ7Func
        func = R.TF1(fitName, fitChoice, fitRange[0], fitRange[1], npar)
        func.SetParameters(1, 10, 1)

    elif fitFunction == "MJ8": ##maybe new default
        npar = 3
        fitChoice = MJ8Func
        func = R.TF1(fitName, fitChoice, fitRange[0], fitRange[1], npar)
        func.SetParameters(1, 10, 1)

    elif fitFunction == "GaussExp":
        npar = 5
        fitChoice = GaussExp
        func = R.TF1(fitName, fitChoice, 500, fitRange[1], npar)
        func.SetParameters(histo.Integral(), histo.GetMean() / 3000., histo.GetRMS() / 3000., 1.0, 900./3000.)
        #func.SetParLimits(0, histo.Integral()*0.1, histo.Integral()*10.0)
        func.SetParLimits(1, 0.0001, 5.)
        func.SetParLimits(2, 0.0001, 5.)
        #func.SetParLimits(3, 0, 1000.)
        #func.FixParameter(4, 1000./3000.)

    elif fitFunction == "ExpModGauss":
        npar = 4
        fitChoice = ExpModGauss
        func = R.TF1(fitName, fitChoice, fitRange[0], fitRange[1], npar)
        func.SetParameters(histo.Integral(), histo.GetMean() / 3000., histo.GetRMS() / 3000., 10.0)
        #func.SetParLimits(0, histo.Integral()*0.1, histo.Integral()*10.0)
        func.SetParLimits(1, 0.0001, 5.)
        func.SetParLimits(2, 0.0001, 5.)
        #func.SetParLimits(3, 0.001, 100.)
    
    elif fitFunction == "Dijet4Param":
         npar = 4
         fitChoice = Dijet4ParamFunc
         func = R.TF1(fitName, fitChoice, fitRange[0], fitRange[1], npar)
         func.SetParameters(-1, 10, -4, 0.01)


    Vmode = ("Q" if not verbose else "")
    Lmode = ("L" if useLikelihood else "")
    fitResult = histo.Fit(fitName, "S0"+Vmode+Lmode, "", fitRange[0], fitRange[1])

    retry = 0
    while retry < 10 and fitResult.Status() != 0:
        print "Retry...", retry
        if retry == 0:
            func.SetParameters(1, 10, 1)
        else:
            func.SetParameters(1 + retry * 2, 10, 1 - retry * 2)

        fitResult = histo.Fit(fitName, "S0"+Vmode+Lmode, "", fitRange[0], fitRange[1])
        retry += 1

    if fitResult.Status() != 0:
        print "\x1b[1;33;41m Error!!! \x1b[0m", "in smoothing fit: did not terminate properly. Exiting", ouutfilepath
        c=R.TCanvas()
        #R.SetOwnership(c,False)
        histo.Draw()
        c.SaveAs("Xhh4bUtils/failed__"+outfileName+".pdf")
        sys.exit(0)
            
    
    cov_TMatrix = fitResult.GetCovarianceMatrix()
    cov = np.zeros( (npar, npar) )
    corr = np.zeros( (npar, npar) )
    for i in range(npar):
        for j in range(npar):
            cov[i,j] = cov_TMatrix[i][j]
            corr[i,j] = cov_TMatrix[i][j] / np.sqrt(cov_TMatrix[i][i] * cov_TMatrix[j][j])

    if verbose:
        print "covariance matrix"
        print cov
        print "correlation matrix"
        print corr
    

    S, U= np.linalg.eigh( cov )
    Sd = np.diag( S )
    sigma = np.sqrt(S)

    evars = GetEigenVariations( cov )

    
    fitFunc = histo.GetFunction(fitName)
    fitProb = fitFunc.GetProb()

    params = array('d',[0]*npar)
    fitFunc.GetParameters( params )


    fit = []
    fiterr = []
    for i in range(npar):
        fit.append(fitFunc.GetParameter(i))
        fiterr.append(fitFunc.GetParError(i))

    z = np.asarray( params )
    z_variations = []

    for i in range(len(evars)): 
        zu = z + evars[i]
        zd = z - evars[i]
        z_variations.append( [zu, zd] )


    drawFunc = R.TF1("drawfit_"+outfileName, fitChoice, fitRange[0], 5000, npar)
    drawFunc.SetParameters( params )

    if makePlots:
        c=R.TCanvas()
        #R.SetOwnership(c,False)
        leg = R.TLegend(0.65,0.68,0.8,0.88)
        leg.SetFillColor(0)
        leg.AddEntry(histo, "Distribution", "LP")
        leg.AddEntry(drawFunc, "Fit", "L")

        histo.SetMarkerStyle(20)
        histo.SetMarkerColor(1)
        histo.SetMarkerSize(1)
        histo.SetXTitle("m_{JJ} [GeV]")
        histo.SetYTitle("Entries")
        histo.Draw()
        drawFunc.Draw("same")


    fvar = []
    for ivar in range(len(z_variations)):
        fup = R.TF1("fup_"+str(ivar)+"_"+outfileName, fitChoice, fitRange[0], 5000, npar)
        fup.SetParameters( z_variations[ivar][0] )
        fup.SetLineColor(colorlist[ivar])
        if makePlots:
            fup.Draw("same")
            # leg.AddEntry(drawFunc, "Fit Variation "+str(ivar), "L")
            leg.AddEntry(fup, "Fit Variation "+str(ivar), "L")          # Qi

        fdw = R.TF1("fdw_"+str(ivar)+"_"+outfileName, fitChoice, fitRange[0], 5000, npar)
        fdw.SetParameters( z_variations[ivar][1] )
        fdw.SetLineColor( colorlist[ivar] )
        if makePlots:
            fdw.Draw("same")

        fvar.append([fup, fdw])
    

    #raw_input()

    if makePlots:
        outplotpath = ouutfilepath
        

        leg.SetBorderSize(0)
        leg.SetMargin(0.3)
        leg.SetTextSize(0.04)
        leg.Draw("same")

        #print fit information
        l = R.TLatex(0.68, 0.9, "Prob: %s; Chi2/NDOF: %s" % \
            (str('%.2g' % func.GetProb()), str('%.2g' % (func.GetChisquare()/(func.GetNDF() * 1.0)))))
        l.SetNDC()
        l.SetTextSize(0.03)
        l.Draw("same")

        #for fitting debug
        #print outfileName, "Prob: %s; Chi2/NDOF: %s" % \
        #    (str('%.2g' % func.GetProb()), str('%.2g' % (func.GetChisquare()/(func.GetNDF() * 1.0))))

        if not os.path.exists(outplotpath):
            os.makedirs(outplotpath)
        c.SetLogy(0)
        c.SaveAs(outplotpath + outfileName + ".pdf")
        c.SetLogy(1)
        c.SaveAs(outplotpath + outfileName + "_l.pdf")

        #fTxT = open(ouutfilepath + outfileName+".pkl", "w")
        # fitResultDict = {
        #   "params": params,
        #   "cov": cov,
        # }
        # pickle.dump(fitResultDict, fTxT)
        # fTxT.close()

    fitResultDict = {
          "params": fit,
          "paramerrs": fiterr,
          "corr": corr,
        }
    return {"nom": drawFunc, "vars":fvar, "res":fitResultDict, "prob":fitProb}


def smoothFuncCompare(histo, fitFunction = "Dijet", fitRange = (900, 3000),  minProb = 0.01, integralMaxRatio = 2.0, makePlots = False, plotExtra = False, verbose = False, outfileName="smoothFuncCompare.root"):

    colorlist = [R.kBlue, R.kGreen, R.kOrange, R.kMagenta, R.kCyan, R.kPink, (R.kAzure+1), R.kGreen+2, R.kOrange+5]

    funclist = ["Dijet","MJ2","MJ3","MJ4","MJ5","MJ6","MJ7","MJ8"]

    funclist_pass = {}    

    namestr = outfileName.split(".root")[0]

    h_clone = histo.Clone()
    h_clone.GetXaxis().SetRangeUser(500, 4000)
    h_clone.GetYaxis().SetRangeUser(1e-2, h_clone.GetMaximum() * 100)
    h_clone.SetDirectory(0)
    
    results = {}
    results_hist = {}
    results_hist_ud = {}
    
    for theFunc in funclist:
        curr_result = smoothfit(h_clone, fitFunction = theFunc, fitRange = fitRange, makePlots = False, verbose = verbose, outfileName = theFunc+"_"+outfileName)
        if curr_result["prob"] <  minProb:
            funclist_pass[theFunc] = False
            print "failed prob", h_clone.GetName(), theFunc
            continue

        if integralMaxRatio is not None and not PassIntegralCondition(h_clone, curr_result["nom"], integralMaxRatio):
            funclist_pass[theFunc] = False
            print "failed norm", h_clone.GetName(), theFunc
            continue
        
        funclist_pass[theFunc] = True

        
        results[theFunc] = curr_result
        results_hist[theFunc] = MakeSmoothHisto(h_clone, results[theFunc]["nom"])
        
        #print results_hist[theFunc].Integral(), results[theFunc]["nom"].Integral(1000, 4000), theFunc

        results_hist_ud[theFunc] = {}
        for ivar in range(len(results[theFunc]["vars"])):
            results_hist_ud[theFunc]["up"+str(ivar)] = MakeSmoothHisto(h_clone, results[theFunc]["vars"][ivar][0])
            results_hist_ud[theFunc]["dw"+str(ivar)] = MakeSmoothHisto(h_clone, results[theFunc]["vars"][ivar][1])


    #print results_hist
    histo_up = results_hist[fitFunction].Clone(histo.GetName() + "_" + namestr + "_up")
    histo_up.SetDirectory(0)
    histo_dw = results_hist[fitFunction].Clone(histo.GetName() + "_" + namestr + "_dw")
    histo_dw.SetDirectory(0)

    histo_up_super = results_hist[fitFunction].Clone(histo.GetName() + "_" + namestr + "_up_super")
    histo_up_super.SetDirectory(0)
    histo_dw_super = results_hist[fitFunction].Clone(histo.GetName() + "_" + namestr + "_dw_super")
    histo_dw_super.SetDirectory(0)

    for ibin in range(1, histo.GetNbinsX()+1):
        deltas = []
        deltas_super = []
        for theFunc in funclist:
            if funclist_pass[theFunc] == False:
                continue
            
            deltas.append( np.abs( results_hist[fitFunction].GetBinContent(ibin) - results_hist[theFunc].GetBinContent(ibin) ) )
            #print results_hist[theFunc].GetBinContent(ibin), results_hist[theFunc].GetBinCenter(ibin), deltas[-1], theFunc
            for ivarh in results_hist_ud[theFunc]:
                deltas_super.append( np.abs( results_hist[fitFunction].GetBinContent(ibin) - results_hist_ud[theFunc][ivarh].GetBinContent(ibin) ) )
            

        theDelta = np.max( deltas )
        theDelta_super = np.max( deltas_super )
        #print theDelta, histo_up.GetBinContent(ibin), histo_up.GetBinCenter(ibin)
        histo_up.SetBinContent(ibin, histo_up.GetBinContent(ibin) + theDelta)
        histo_dw.SetBinContent(ibin, histo_dw.GetBinContent(ibin) - theDelta)
        
        histo_up_super.SetBinContent(ibin, histo_up.GetBinContent(ibin) + theDelta_super)
        histo_dw_super.SetBinContent(ibin, histo_dw.GetBinContent(ibin) - theDelta_super)

    smoothFuncCompSyst = {"up":histo_up, "dw":histo_dw, "up_super":histo_up_super, "dw_super":histo_dw_super}
        
        

    if makePlots:
        f = R.TFile(outfileName, "RECREATE")
        
        c=R.TCanvas("c1_func","c1_func")
        #R.SetOwnership(c,False)
        leg1 = R.TLegend(0.6,0.6,0.9,0.9)
        leg1.SetFillColor(0)
        leg1.SetBorderSize(0)
        leg1.SetMargin(0.3)

        leg2 = R.TLegend(0.2,0.6,0.55,0.9)
        leg2.SetFillColor(0)
        leg2.SetBorderSize(0)
        leg2.SetMargin(0.3)
    
        h_clone.SetLineColor(R.kBlack)
        h_clone.Draw()
        leg1.AddEntry(histo, "Histogram", "L")
        leg2.AddEntry(histo, "Histogram", "P")

        icol = 0
        ivar0 = True
        err_hist_ratio = None
        for theFunc in funclist:

            if funclist_pass[theFunc] == False:
                continue

            if theFunc==fitFunction:
                err_hist = MakeSmoothHisto(histo, results[theFunc]["nom"])
                err_hist.SetDirectory(0)

                for ivar in range(len(results[theFunc]["vars"])):
                    err_hist_up = MakeSmoothHisto(histo, results[theFunc]["vars"][ivar][0])
                    err_hist_dw = MakeSmoothHisto(histo, results[theFunc]["vars"][ivar][1])

                    for ibin in range(1, err_hist.GetNbinsX()+1):
                        err_val = np.max( np.abs( [ err_hist.GetBinContent(ibin) - err_hist_up.GetBinContent(ibin), err_hist.GetBinContent(ibin) - err_hist_dw.GetBinContent(ibin)] ) )
                        err_hist.SetBinError(ibin, np.sqrt( err_hist.GetBinError(ibin)**2 + err_val**2) )

                err_hist_ratio = err_hist.Clone("err_hist_ratio__"+namestr)
                err_hist_ratio.SetDirectory(0)
                err_hist_ratio.Divide( results[theFunc]["nom"] )


                err_hist.SetFillColor(R.kBlack)
                err_hist.SetMarkerSize(0)
                err_hist.SetFillStyle(3001)
                err_hist.Draw("sameE3")
                leg1.AddEntry(err_hist, "smoothing error", "F")
                leg2.AddEntry(err_hist, "smoothing error", "F")



            #print results[theFunc]["nom"], results[theFunc]["nom"].Eval(1000), results[theFunc]["nom"].Eval(2000), results[theFunc]["nom"].Eval(3000)

            results[theFunc]["nom"].SetLineColor( colorlist[icol] )
            results[theFunc]["nom"].Draw("same")
            leg1.AddEntry(results[theFunc]["nom"], theFunc, "L")
            leg2.AddEntry(results[theFunc]["nom"], theFunc, "L")


            if plotExtra:
                for ivar in range(len(results[theFunc]["vars"])):
                    results[theFunc]["vars"][ivar][0].SetLineColor( R.kGray+2 )
                    results[theFunc]["vars"][ivar][0].Draw("same")
                    results[theFunc]["vars"][ivar][1].SetLineColor( R.kGray+2 )
                    results[theFunc]["vars"][ivar][1].Draw("same")
                    
            

            icol += 1
            
        if plotExtra:
            results_hist_ud["MJ2"]["up0"].SetLineColor( R.kGray+2 )
            leg1.AddEntry(results_hist_ud["MJ2"]["up0"], "Param Variations", "L")
            leg2.AddEntry(results_hist_ud["MJ2"]["up0"], "Param Variations", "L")
        leg1.Draw()
        


        c2=R.TCanvas("c2_func","c2_func")
        #R.SetOwnership(c,False)
        print "err_hist_ratio",err_hist_ratio
        err_hist_ratio.SetFillColor(R.kBlack)
        err_hist_ratio.SetFillStyle(3004)
        err_hist_ratio.SetMarkerSize(0)
        err_hist_ratio.GetXaxis().SetRangeUser(1000, 4000)
        err_hist_ratio.GetYaxis().SetRangeUser(0.4, 2)
        err_hist_ratio.GetXaxis().SetLabelSize(0.04)
        err_hist_ratio.GetYaxis().SetLabelSize(0.04)
        err_hist_ratio.Draw("E2")

        icol = 0
        f_ratio = {}

        delta_ratio_super = {}
        h_ratio_lst  = []
        for theFunc in funclist:
            if funclist_pass[theFunc] == False:
                continue
            
            ## func_ratio[theFunc] = lambda x: (results[theFunc]["nom"].Eval(x[0]) / results["Exp"]["nom"].Eval(x[0]))
            ## f_ratio[theFunc] = R.TF1(theFunc+"_ratio_"+namestr, func_ratio[theFunc], fitRange[0], 3000, 0)
            ## f_ratio[theFunc].SetLineColor( colorlist[icol] )
            ## f_copy = f_ratio[theFunc].DrawCopy("same")
            h_ratio_lst.append(results[theFunc]["nom"].GetHistogram().Clone())
            h_ratio_lst[-1].SetName("err_hist_ratio__" + theFunc + "_copy")
            h_ratio_lst[-1].SetDirectory(0)
            h_ratio_lst[-1].Divide( results[fitFunction]["nom"] )
            h_ratio_lst[-1].SetLineColor( colorlist[icol] )
            h_ratio_lst[-1].Draw("hist same")

            if plotExtra:
                for ivar in range(len(results[theFunc]["vars"])):
                    h_ratio_ud = results[theFunc]["vars"][ivar][0].GetHistogram()
                    h_ratio_ud.Divide( results[fitFunction]["nom"] )
                    h_ratio_ud.SetDirectory(0)
                    h_ratio_ud.SetLineColor(R.kGray+2)
                    h_ratio_ud.Draw("same")

                    delta_ratio_super[theFunc+"_"+str(ivar)+"_up"] = h_ratio_ud.GetBinContent( h_ratio_ud.FindBin(3000) )


                    h_ratio_ud = results[theFunc]["vars"][ivar][1].GetHistogram()
                    h_ratio_ud.Divide( results[fitFunction]["nom"] )
                    h_ratio_ud.SetDirectory(0)
                    h_ratio_ud.SetLineColor(R.kGray+2)
                    h_ratio_ud.Draw("same")

                    delta_ratio_super[theFunc+"_"+str(ivar)+"_dw"] = h_ratio_ud.GetBinContent( h_ratio_ud.FindBin(3000) )
            
            #print f_copy, f_ratio[theFunc], f_ratio[theFunc].Eval(1000), f_ratio[theFunc].Eval(2000), f_ratio[theFunc].Eval(3000)
            icol += 1

        leg2.Draw()
        #for drs in delta_ratio_super:
        #    print drs, delta_ratio_super[drs]

        c.SetLogy()

        f.WriteTObject(c)
        f.WriteTObject(c2)

        c.SaveAs(outfileName.split(".root")[0] + "_comp" + ".pdf")
        c2.SaveAs(outfileName.split(".root")[0] + "_comp" + "_ratio.pdf")
        c.Close()
        c2.Close()
        f.Close()
        
    return smoothFuncCompSyst


def smoothFuncRangeCompare(histo, fitFunction = "Dijet", fitRange = (900, 3000), fitMaxVals = ["2000", "1500", "1750"], fitMinVals=["900","1000","1100"], minProb = 0.01, integralMaxRatio = 2.0, makePlots = False, plotExtra = False, verbose = False, outfileName="smoothFuncRangeCompare.root"):

    colorlist = [R.kBlue, R.kGreen, R.kOrange, R.kMagenta, R.kCyan, R.kPink, (R.kAzure+1), R.kGreen+2, R.kOrange+5]        

    namestr = outfileName.split(".root")[0]

    strInMax = str(fitRange[1])

    strNom = str(fitRange[0]) + "_" + str(fitRange[1])

    h_clone = histo.Clone()
    h_clone.SetDirectory(0)
    
    results = {}
    results_hist = {}
    results_hist_ud = {}

    fitPairs={}
    for maxRange in fitMaxVals:
        fitPairs[str(fitRange[0])+"_"+maxRange] =   (fitRange[0], float(maxRange)) 
    for minRange in fitMinVals:
        fitPairs[minRange + "_"+str(fitRange[1])] = (float(minRange),fitRange[1])

    fitPairs_pass={}
    

    for fpair in fitPairs:
        curr_result = smoothfit(h_clone, fitFunction = fitFunction, fitRange = fitPairs[fpair], outrange_start = fitRange[0], makePlots = False, verbose = verbose, outfileName =maxRange+"_"+outfileName)
        if curr_result["prob"] < minProb:
            fitPairs_pass[fpair] = False
            continue

        if integralMaxRatio is not None and not PassIntegralCondition(h_clone, curr_result["nom"], integralMaxRatio):
            fitPairs_pass[fpair] = False
            continue

        fitPairs_pass[fpair] = True
        
        results[fpair] = curr_result
        results_hist[fpair] = MakeSmoothHisto(h_clone, results[fpair]["nom"])

        results_hist_ud[fpair] = {}
        for ivar in range(len(results[fpair]["vars"])):
            results_hist_ud[fpair]["up"+str(ivar)] = MakeSmoothHisto(h_clone, results[fpair]["vars"][ivar][0])
            results_hist_ud[fpair]["dw"+str(ivar)] = MakeSmoothHisto(h_clone, results[fpair]["vars"][ivar][1])


    if verbose:
        for maxRange in fitMaxVals:
            print ("MaxRange=",maxRange, "Integral(",maxRange,",3000)=", results[maxRange]["nom"].Integral(float(maxRange), 3000),
                   "nominal Integral(",maxRange,",3000)=", results["3000"]["nom"].Integral(float(maxRange), 3000) )

    histo_up = results_hist[strNom].Clone(histo.GetName() + "_" + namestr + "_up")
    histo_up.SetDirectory(0)
    histo_dw = results_hist[strNom].Clone(histo.GetName() + "_" + namestr + "_dw")
    histo_dw.SetDirectory(0)

    histo_up_super = results_hist[strNom].Clone(histo.GetName() + "_" + namestr + "_up_super")
    histo_up_super.SetDirectory(0)
    histo_dw_super = results_hist[strNom].Clone(histo.GetName() + "_" + namestr + "_dw_super")
    histo_dw_super.SetDirectory(0)

    for ibin in range(1, histo.GetNbinsX()+1):
        deltas = []
        deltas_super = []
        for fpair in fitPairs:
            if fitPairs_pass[fpair] == False:
                continue
            
            deltas.append( np.abs( results_hist[strNom].GetBinContent(ibin) - results_hist[fpair].GetBinContent(ibin) ) )
            #print results_hist[strNom].GetBinContent(ibin), results_hist[strNom].GetBinCenter(ibin), deltas[-1]

            for ivarh in results_hist_ud[fpair]:
                deltas_super.append( np.abs( results_hist[strNom].GetBinContent(ibin) - results_hist_ud[fpair][ivarh].GetBinContent(ibin) ) )
            

        theDelta = np.max( deltas )
        theDelta_super = np.max( deltas_super )
        #print theDelta, histo_up.GetBinContent(ibin), histo_up.GetBinCenter(ibin)
        histo_up.SetBinContent(ibin, histo_up.GetBinContent(ibin) + theDelta)
        histo_dw.SetBinContent(ibin, histo_dw.GetBinContent(ibin) - theDelta)
        
        histo_up_super.SetBinContent(ibin, histo_up.GetBinContent(ibin) + theDelta_super)
        histo_dw_super.SetBinContent(ibin, histo_dw.GetBinContent(ibin) - theDelta_super)

    smoothFuncCompSyst = {"up":histo_up, "dw":histo_dw, "up_super":histo_up_super, "dw_super":histo_dw_super}

    if makePlots:
        f = R.TFile(outfileName, "RECREATE")
        
        c=R.TCanvas("c1_range","c1_range")
        #R.SetOwnership(c,False)
        leg1 = R.TLegend(0.6,0.6,0.9,0.9)
        leg1.SetFillColor(0)
        leg1.SetBorderSize(0)
        leg1.SetMargin(0.3)

        leg2 = R.TLegend(0.3,0.6,0.65,0.9)
        leg2.SetFillColor(0)
        leg2.SetBorderSize(0)
        leg2.SetMargin(0.3)
    
        h_clone.SetLineColor(R.kBlack)
        h_clone.GetXaxis().SetRangeUser(500, 4000)
        h_clone.GetYaxis().SetRangeUser(1e-2, h_clone.GetMaximum() * 100)
        h_clone.SetTitle("")
        h_clone.Draw()
        
        leg1.AddEntry(histo, "Histogram", "P")
        leg2.AddEntry(histo, "Histogram", "L")

        icol = 0
        ivar0 = True
        err_hist_ratio = None
        for fpair in fitPairs:
            if fitPairs_pass[fpair] == False:
                continue

            if fpair==strNom:
                err_hist = MakeSmoothHisto(histo, results[fpair]["nom"])
                err_hist.SetDirectory(0)

                for ivar in range(len(results[fpair]["vars"])):
                    err_hist_up = MakeSmoothHisto(histo, results[fpair]["vars"][ivar][0])
                    err_hist_dw = MakeSmoothHisto(histo, results[fpair]["vars"][ivar][1])

                    for ibin in range(1, err_hist.GetNbinsX()+1):
                        err_val = np.max( np.abs( [ err_hist.GetBinContent(ibin) - err_hist_up.GetBinContent(ibin), err_hist.GetBinContent(ibin) - err_hist_dw.GetBinContent(ibin)] ) )
                        err_hist.SetBinError(ibin, np.sqrt( err_hist.GetBinError(ibin)**2 + err_val**2) )

                err_hist_ratio = err_hist.Clone("err_hist_ratio__"+namestr)
                err_hist_ratio.SetDirectory(0)
                err_hist_ratio.Divide( results[fpair]["nom"] )

                err_hist.SetFillColor(R.kBlack)
                err_hist.SetMarkerSize(0)
                err_hist.SetFillStyle(3001)
                err_hist.Draw("sameE3")

                leg1.AddEntry(err_hist, "smoothing error", "F")
                leg2.AddEntry(err_hist, "smoothing error", "F")


            #print results[theFunc]["nom"], results[theFunc]["nom"].Eval(1000), results[theFunc]["nom"].Eval(2000), results[theFunc]["nom"].Eval(3000)

            results[fpair]["nom"].SetLineColor( colorlist[icol] )
            results[fpair]["nom"].Draw("same")

            leg1.AddEntry(results[fpair]["nom"], fpair, "L")
            leg2.AddEntry(results[fpair]["nom"], fpair, "L")


            if plotExtra:
                for ivar in range(len(results[fpair]["vars"])):
                    results[fpair]["vars"][ivar][0].SetLineColor( R.kGray+2 )
                    results[fpair]["vars"][ivar][0].Draw("same")
                    results[fpair]["vars"][ivar][1].SetLineColor( R.kGray+2 )
                    results[fpair]["vars"][ivar][1].Draw("same")
                    
            icol += 1

        results_hist_ud[strNom]["up0"].SetLineColor( R.kGray+2 )
        if plotExtra:
            leg1.AddEntry(results_hist_ud[strNom]["up0"], "Param Variations", "L")
            leg2.AddEntry(results_hist_ud[strNom]["up0"], "Param Variations", "L")
        
        leg1.Draw()
        c.SetLogy(1)

        c2=R.TCanvas("c2_range","c2_range")
        c2.cd()
        #R.SetOwnership(c,False)
        print "err_hist_ratio",err_hist_ratio
        err_hist_ratio.SetFillColor(R.kBlack)
        err_hist_ratio.SetFillStyle(3004)
        err_hist_ratio.SetMarkerSize(0)
        err_hist_ratio.GetXaxis().SetRangeUser(1000, 4000)
        err_hist_ratio.GetYaxis().SetRangeUser(0.4, 2)
        err_hist_ratio.GetXaxis().SetLabelSize(0.04)
        err_hist_ratio.GetYaxis().SetLabelSize(0.04)
        err_hist_ratio.SetTitle("")
        err_hist_ratio.Draw("E2")

        icol = 0
        f_ratio = {}

        delta_ratio_super = {}
        h_ratio_lst  = []
        for fpair in fitPairs:
            if fitPairs_pass[fpair] == False:
                continue
            
            ## func_ratio[theFunc] = lambda x: (results[theFunc]["nom"].Eval(x[0]) / results["Exp"]["nom"].Eval(x[0]))
            ## f_ratio[theFunc] = R.TF1(theFunc+"_ratio_"+namestr, func_ratio[theFunc], fitRange[0], 3000, 0)
            ## f_ratio[theFunc].SetLineColor( colorlist[icol] )
            ## f_copy = f_ratio[theFunc].DrawCopy("same")
            h_ratio_lst.append(results[fpair]["nom"].GetHistogram().Clone())
            h_ratio_lst[-1].SetName("err_hist_ratio__" + fpair + "_copy")
            h_ratio_lst[-1].SetDirectory(0)
            h_ratio_lst[-1].Divide( results[strNom]["nom"] )
            h_ratio_lst[-1].SetLineColor( colorlist[icol] )
            h_ratio_lst[-1].Draw("hist same")

            if plotExtra:
                for ivar in range(len(results[fpair]["vars"])):
                    h_ratio_ud = results[fpair]["vars"][ivar][0].GetHistogram()
                    h_ratio_ud.Divide( results[strNom]["nom"] )
                    h_ratio_ud.SetDirectory(0)
                    h_ratio_ud.SetLineColor(R.kGray+2)
                    h_ratio_ud.Draw("same")

                    delta_ratio_super[fpair+"_"+str(ivar)+"_up"] = h_ratio_ud.GetBinContent( h_ratio_ud.FindBin(3000) )


                    h_ratio_ud = results[fpair]["vars"][ivar][1].GetHistogram()
                    h_ratio_ud.Divide( results[strNom]["nom"] )
                    h_ratio_ud.SetDirectory(0)
                    h_ratio_ud.SetLineColor(R.kGray+2)
                    h_ratio_ud.Draw("same")

                    delta_ratio_super[fpair+"_"+str(ivar)+"_dw"] = h_ratio_ud.GetBinContent( h_ratio_ud.FindBin(3000) )
            
            #print h_ratio, h_ratio.Eval(1000), h_ratio.Eval(2000), h_ratio.Eval(3000)
            icol += 1

        leg2.Draw()

        if verbose:
            for drs in delta_ratio_super:
                print drs, delta_ratio_super[drs]

        c.SaveAs(outfileName.split(".root")[0] + "_comp" + ".pdf")
        c2.SaveAs(outfileName.split(".root")[0] + "_comp" + "_ratio.pdf")
        f.WriteTObject(c)
        f.WriteTObject(c2)
        c.Close()
        c2.Close()
        f.Close()
        
    return smoothFuncCompSyst


#Pass in a histo and it will return a histogram with scaled bins above a min_mass_to_rebin.
def VariableRebin(histo, bin_multiplier, min_mass_to_rebin):
    bin_list = [0]
    cur_bin = 1
    first_transition = True
    current_bin_width = histo.GetBinWidth(cur_bin)
    while cur_bin < histo.GetNbinsX():
        if histo.GetBinLowEdge(cur_bin) <= 0:
            cur_bin += 1
        elif histo.GetBinLowEdge(cur_bin) < min_mass_to_rebin:
            bin_list.append(histo.GetBinLowEdge(cur_bin))
            cur_bin += 1
        else:
            if first_transition:
                cur_bin -=1
                first_transition = False
            varied_edge = histo.GetBinLowEdge(cur_bin) + histo.GetBinWidth(cur_bin)*bin_multiplier
            cur_bin = histo.FindBin(varied_edge)
            if cur_bin > histo.GetNbinsX():
                break
            bin_list.append(varied_edge)
            
    bin_array = np.array(bin_list)
    print bin_array            
    return histo.Rebin(len(bin_list)-1,"variable rebinned histo", bin_array)


def MakeSmoothHisto(hist, fitCurve, lowFillVal = 500, keepNorm=True):   # qi
    low=R.Double(0.0)
    high=R.Double(0.0)
    fitCurve.GetRange(low, high)

    oldIntegral = 0
    newIntegral = 0

    hist_smooth = hist.Clone(hist.GetName()+"__smooth")
    for ibin in range(0, hist_smooth.GetNbinsX()+1):

        if hist_smooth.GetBinLowEdge(ibin) >= high: ##don't do anything above 4 TeV for now!!
            hist_smooth.SetBinContent(ibin, 0)
            hist_smooth.SetBinError(ibin, 0)
            continue

        if hist_smooth.GetBinLowEdge(ibin) >= low:
            oldIntegral += hist_smooth.GetBinContent(ibin)
            newIntegral += fitCurve.Integral(hist_smooth.GetBinLowEdge(ibin), hist_smooth.GetBinLowEdge(ibin+1))/hist_smooth.GetBinWidth(ibin)
            #print oldIntegral, newIntegral, ibin, hist_smooth.GetBinCenter(ibin)

            hist_smooth.SetBinContent(ibin, 0)
            hist_smooth.SetBinError(ibin, 0)
        

    if keepNorm:
        #print oldIntegral, newIntegral
        hist_smooth.Add(fitCurve, oldIntegral/newIntegral)
    else:
        hist_smooth.Add(fitCurve, 1.0)

    ##clear head
    for ibin in range(1, hist_smooth.FindBin(lowFillVal)):
        hist_smooth.SetBinContent(ibin, 0)
        hist_smooth.SetBinError(ibin, 0)
    ##clear tail
    for ibin in range(hist_smooth.FindBin(high), hist_smooth.GetNbinsX() + 1):
        hist_smooth.SetBinContent(ibin, 0)
        hist_smooth.SetBinError(ibin, 0)

    #print fitCurve.GetName(), oldIntegral, newIntegral, low, high, hist.Integral(), hist_smooth.Integral(), fitCurve.Integral(1200, 4000)
    return hist_smooth

def MakeSmoothHistoWithError(hist, smoothResult, lowFillVal=500, keepNorm=True):
    h_nominal = MakeSmoothHisto(hist, smoothResult["nom"], lowFillVal, keepNorm)

    sysList = []
    for ivar in range(len(smoothResult["vars"])):
        h_var_up = MakeSmoothHisto(hist, smoothResult["vars"][ivar][0], lowFillVal, keepNorm)
        h_var_down = MakeSmoothHisto(hist, smoothResult["vars"][ivar][1], lowFillVal, keepNorm)

        sysList.append( [h_var_up, h_var_down] )

    h_sysMerged = AddSysErrorToHist(h_nominal, sysList)
    return h_sysMerged["h_sysSum"]


# add systematics on histogram as error bar
# each item in sysList should be a set of either length 2 or 1; Length 2 means up/down; Length 1 means we symmetrize the variation
# two histogram will be returned: one is with systematic variation only, and another one is the sum of systematic variation and whatever original error on nominal histogram
# the error on normalization will also be returned
def AddSysErrorToHist(h, sysList):
    h_original = h.Clone()
    h_original.SetDirectory(0)

    nbins = h.GetNbinsX()

    h_sysOnly = h.Clone(h.GetName()+"_sysOnly")
    h_sysOnly.SetDirectory(0)
    for ibin in range(0, h_sysOnly.GetNbinsX()+2):
        h_sysOnly.SetBinError(ibin, 0.)

    h_sysSum = h.Clone(h.GetName()+"_sysSum")
    h_sysSum.SetDirectory(0)

    normSysError = 0.

    for oneSys in sysList:
        # distribution
        for ibin in range(0, h.GetNbinsX()+2):
            if len(oneSys) == 2:
                sysError1 = abs(oneSys[0].GetBinContent(ibin) - h.GetBinContent(ibin))
                sysError2 = abs(oneSys[1].GetBinContent(ibin) - h.GetBinContent(ibin))
            elif len(oneSys) == 1:
                sysError1 = abs(oneSys[0].GetBinContent(ibin) - h.GetBinContent(ibin))
                sysError2 = sysError1
            else:
                print "Invalid systematic item:",oneSys,"Zero systematic will be set!"
                sysError1 = 0.
                sysError2 = 0.

            sysError = max(sysError1, sysError2)

            h_sysOnly.SetBinError(ibin, R.TMath.Sqrt((h_sysOnly.GetBinError(ibin))**2 + sysError**2) )
            h_sysSum.SetBinError(ibin, R.TMath.Sqrt((h_sysSum.GetBinError(ibin))**2 + sysError**2) )

        # normalization
        if len(oneSys) == 2:
            normSysError1 = abs(oneSys[0].Integral(0, nbins+1) - h.Integral(0, nbins+1))
            normSysError2 = abs(oneSys[1].Integral(0, nbins+1) - h.Integral(0, nbins+1))
        elif len(oneSys) == 1:
            normSysError1 = abs(oneSys[0].Integral(0, nbins+1) - h.Integral(0, nbins+1))
            normSysError2 = normSysError1

        # debug
        # print h_original.Integral(0, nbins+1), oneSys[0].Integral(0, nbins+1) - h.Integral(0, nbins+1)

        currentSysError = max(normSysError1, normSysError2)
        normSysError = R.TMath.Sqrt(normSysError**2 + currentSysError**2)

    normStatsError = R.Double(0.)
    norm = h_original.IntegralAndError(0, nbins+1, normStatsError)

    normSumError = R.TMath.Sqrt( normStatsError**2 + normSysError**2 )

    outputDict = {
                   'h_sysOnly': h_sysOnly,
                   'h_sysSum':  h_sysSum,
                   'h_original': h_original,

                   'normSysError': normSysError,
                   "normStatsError": normStatsError,
                   "normSumError": normSumError,

                   "norm": norm,
                 }

    return outputDict

def PassIntegralCondition(hist, func, integralMaxRatio, testRanges = [[1200, 1500], [1500, 2000], [2000, 2500]] ):
    for itest in testRanges:
        h_int = hist.Integral( hist.FindBin( itest[0] ), hist.FindBin( itest[1] ), "width")
        f_int = func.Integral( hist.GetXaxis().GetBinLowEdge(hist.FindBin( itest[0] )), hist.GetXaxis().GetBinLowEdge(hist.FindBin( itest[1]) + 1))
        
        if f_int == 0:
            print "PassIntegralCondition:: function integral is Zero!!!"
            sys.exit(0)

        if h_int / f_int < 1/(integralMaxRatio) or  h_int / f_int > integralMaxRatio:
            print itest, h_int, f_int, hist.GetXaxis().GetBinCenter(hist.FindBin( itest[0] ))
            return False

    return True

############################################################################################
### functions
############################################################################################
def ExpoFunc(x, par):
    return np.exp(-par[0]*x[0] + par[1])

def DijetFunc(x, par):
    #x = [np.array(i, dtype="float128") for i in x]
    #par = [np.array(i,dtype="float128") for i in par]
    z = x[0] / 13000.0
    #return np.exp( par[0] ) * np.power((1.0 - z), par[1]) * np.power(z, par[2])
    return np.exp(par[0] + par[1]* np.log((1.0 - z)) + par[2]* np.log(z)  )
    #return np.exp(par[0]) * np.power(1.0-z, par[1]) * np.power(z, par[2])

def MJ2Func(x, par):
    #https://cds.cern.ch/record/2026080
    z = x[0] / 13000.0
    return np.exp(par[0] + par[1]* np.log((1.0 - x[0] / 13000.0)) + par[2]*(x[0] / 13000.0)**2  )

def MJ3Func(x, par):
    #https://cds.cern.ch/record/2026080
    z = x[0] / 13000.0
    return np.exp(par[0] + par[1]* np.log((1.0 - z)) + par[2]*z*np.log(z)  )

def MJ4Func(x, par):
    #https://cds.cern.ch/record/2026080
    z = x[0] / 13000.0
    return np.exp( par[0] ) * np.power((1.0 - z), par[1]) * np.power(z, par[2]*np.log(z))

def MJ5Func(x, par):
    #https://cds.cern.ch/record/2026080
    z = x[0] / 13000.0
    return np.exp( par[0] ) * np.power((1.0 - z), par[1]) * np.power((1.0 + z), par[2]*z)

def MJ6Func(x, par):
    #https://cds.cern.ch/record/2026080
    z = x[0] / 13000.0
    return np.exp( par[0] ) * np.power((1.0 - z), par[1]) * np.power((1.0 + z), par[2]*np.log(z))

def MJ7Func(x, par):
    #https://cds.cern.ch/record/2026080
    z = x[0] / 13000.0
    return np.exp( par[0] ) * (1.0 / z) * np.power((1.0 - z), par[1] - par[2]*np.log(z))

def MJ8Func(x, par):
    #https://cds.cern.ch/record/2026080
    z = x[0] / 13000.0
    return np.exp( par[0] ) * (1.0 / z**2) * np.power((1.0 - z), par[1] - par[2]*np.log(z))

def Dijet4ParamFunc(x, par):
    z = x[0] / 13000.0
    return np.exp(par[0] + par[1]* np.log((1.0 - z)) + par[2]* np.log(z) ) * np.power(z, par[3]*np.log(z))

def GaussExp(x, par):
    z = x[0] / 3000.
    norm = par[0]
    mu = par[1]
    sigma = par[2]
    beta = par[3]
    crossover = par[4]

    ## nu = beta - (crossover - mu) / (sigma**2)
    ## #crossover = (beta- nu) * (sigma**2) + mu
    ## normExp = np.exp( -0.5 * ( (crossover - mu) / sigma )**2 + beta * crossover )
    ## f = norm * ( np.exp( -0.5 * ( (z - mu) / sigma )**2 ) * (1.0 / (1 + np.exp(-1.0 * nu * (crossover - z)) ) ) +
    ##              normExp * np.exp( -1.0 * beta * z )     * (1.0 / (1 + np.exp(-1.0 * nu * (z - crossover )) ) ) )
    

    if z < crossover:
        f = norm * np.exp( -0.5 * ( (z - mu) / sigma )**2 )
    else:
        f = norm * np.exp( -0.5 * ( (crossover - mu) / sigma )**2 + beta * crossover ) * np.exp( -1.0 * beta * z )

    return f


def ExpModGauss(x, par):
    z = x[0] / 3000.
    norm = par[0]
    mu = par[1]
    sigma = par[2]
    beta = par[3]

    #return norm  *(beta/2.0) * (np.power( np.exp(  (2*mu + beta*np.power(sigma,2) - 2*z)), (beta/2.0)) *
    #                            scipy.special.erfc( (mu + beta*np.power(sigma,2) - z)/(sigma * np.sqrt(2)) ) )

    #return (np.exp( np.log(norm) +
    #                np.log(beta/2.0) +
    #                (beta/2.0)*(2*mu + beta*np.power(sigma,2) - 2*z)) *
    #        scipy.special.erfc( (mu + beta*np.power(sigma,2) - z)/(sigma * np.sqrt(2)) ) )

    return (R.TMath.Exp( R.TMath.Log(norm) + R.TMath.Log(beta/2.0) +
                        (beta/2.0)*(2*mu + beta*R.TMath.Power(sigma,2) - 2*z)) *
            R.TMath.Erfc( (mu + beta*R.TMath.Power(sigma,2) - z)/(sigma * np.sqrt(2)) ) )
    




############################################################################################

if __name__=="__main__":
    datafile = R.TFile("hist_data.root ","READ")
    h = datafile.Get("GoodEvent_Pass4GoodTrackJetPass2b77PassSRMass/DiJetMass").Clone()
    smoothfit( h, fitFunction = "Dijet", fitRange = (900, 3000), makePlots = True, verbose = True )
    
