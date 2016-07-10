import ROOT as R

import numpy as np
import scipy.special
from array import array
import sys

from GetEigenVariations import GetEigenVariations

import smoothfit

import cPickle as pickle

import time



def smoothFuncCompare(histo, fitFunction = "ExpModGauss", fitRange = (100, 3000), funcCompareRange=(900,3000), makePlots = False, plotExtra = True, verbose = False, outfileName="ExpModGaussSmoothFuncCompare.root"):


    
    colorlist = [R.kBlue, R.kGreen, R.kOrange, R.kMagenta, R.kCyan, R.kPink, (R.kAzure+1), R.kGreen+2]        

    namestr = outfileName.split(".root")[0]

    h_clone = histo.Clone()
    h_clone.SetDirectory(0)


    nominal_result = smoothfit.smoothfit(h_clone, fitFunction = fitFunction, fitRange = fitRange, makePlots = False, verbose = verbose, outfileName = fitFunction+"_"+outfileName)

    nominal_hist = smoothfit.MakeSmoothHisto(h_clone, nominal_result["nom"])
    
    results = {}
    results_hist = {}
    results_hist_ud = {}

    for theFunc in ["Exp","MJ2","MJ3","MJ4","MJ5","MJ6","MJ7","MJ8"]:
        results[theFunc] = smoothfit.smoothfit(h_clone, fitFunction = theFunc, fitRange = funcCompareRange, makePlots = False, verbose = verbose, outfileName = theFunc+"_"+outfileName)
        results_hist[theFunc] = smoothfit.MakeSmoothHisto(h_clone, results[theFunc]["nom"])

        results_hist_ud[theFunc] = {}
        for ivar in range(len(results[theFunc]["vars"])):
            results_hist_ud[theFunc]["up"+str(ivar)] = smoothfit.MakeSmoothHisto(h_clone, results[theFunc]["vars"][ivar][0])
            results_hist_ud[theFunc]["dw"+str(ivar)] = smoothfit.MakeSmoothHisto(h_clone, results[theFunc]["vars"][ivar][1])



    histo_up = nominal_hist.Clone(histo.GetName() + "_" + namestr + "_up")
    histo_up.SetDirectory(0)
    histo_dw = nominal_hist.Clone(histo.GetName() + "_" + namestr + "_dw")
    histo_dw.SetDirectory(0)

    histo_up_super = nominal_hist.Clone(histo.GetName() + "_" + namestr + "_up_super")
    histo_up_super.SetDirectory(0)
    histo_dw_super = nominal_hist.Clone(histo.GetName() + "_" + namestr + "_dw_super")
    histo_dw_super.SetDirectory(0)

    for ibin in range(1, histo.GetNbinsX()+1):
        if histo.GetBinLowEdge(ibin) + histo.GetBinWidth(ibin) < funcCompareRange[0]:
            continue
        
        deltas = []
        deltas_super = []
        for theFunc in ["Exp","MJ2","MJ3","MJ4","MJ5","MJ6","MJ7","MJ8"]:
            deltas.append( np.abs( nominal_hist.GetBinContent(ibin) - results_hist[theFunc].GetBinContent(ibin) ) )

            for ivarh in results_hist_ud[theFunc]:
                deltas_super.append( np.abs( nominal_hist.GetBinContent(ibin) - results_hist_ud[theFunc][ivarh].GetBinContent(ibin) ) )
            

        theDelta = np.max( deltas )
        theDelta_super = np.max( deltas_super )
        histo_up.SetBinContent(ibin, histo_up.GetBinContent(ibin) + theDelta)
        histo_dw.SetBinContent(ibin, histo_dw.GetBinContent(ibin) - theDelta)
        
        histo_up_super.SetBinContent(ibin, histo_up.GetBinContent(ibin) + theDelta_super)
        histo_dw_super.SetBinContent(ibin, histo_dw.GetBinContent(ibin) - theDelta_super)

    smoothFuncCompSyst = {"up":histo_up, "dw":histo_dw, "up_super":histo_up_super, "dw_super":histo_dw_super}
        
        

    if makePlots:
        f = R.TFile(outfileName, "RECREATE")
        
        c=R.TCanvas("c1_err_hist","c1")
        #R.SetOwnership(c,False)
        xleg, yleg = 0.52, 0.7
        leg = R.TLegend(xleg, yleg, xleg+0.3, yleg+0.2)
        leg.SetFillColor(0)
        leg.SetBorderSize(0)
        leg.SetMargin(0.3)
    
        h_clone.SetLineColor(R.kBlack)
        h_clone.Draw()
        leg.AddEntry(histo, "Histogram", "L")

        icol = 0
        ivar0 = True
        err_hist_ratio = None


        err_hist = smoothfit.MakeSmoothHisto(histo, nominal_result["nom"])
        err_hist.SetDirectory(0)

        for ivar in range(len(nominal_result["vars"])):
            err_hist_up = smoothfit.MakeSmoothHisto(histo, nominal_result["vars"][ivar][0])
            err_hist_dw = smoothfit.MakeSmoothHisto(histo, nominal_result["vars"][ivar][1])

            for ibin in range(1, err_hist.GetNbinsX()+1):
                err_val = np.max( np.abs( [ err_hist.GetBinContent(ibin) - err_hist_up.GetBinContent(ibin), err_hist.GetBinContent(ibin) - err_hist_dw.GetBinContent(ibin)] ) )
                err_hist.SetBinError(ibin, np.sqrt( err_hist.GetBinError(ibin)**2 + err_val**2) )

        err_hist_ratio = err_hist.Clone("err_hist_ratio__"+namestr)
        err_hist_ratio.SetDirectory(0)
        err_hist_ratio.Divide( nominal_result["nom"] )


        err_hist.SetFillColor(R.kBlack)
        err_hist.SetFillStyle(3001)
        err_hist.Draw("sameE3")
        leg.AddEntry(err_hist, "smoothing error", "F")

        
        for theFunc in ["Exp","MJ2","MJ3","MJ4","MJ5","MJ6","MJ7","MJ8"]:

            #print results[theFunc]["nom"], results[theFunc]["nom"].Eval(1000), results[theFunc]["nom"].Eval(2000), results[theFunc]["nom"].Eval(3000)

            results[theFunc]["nom"].SetLineColor( colorlist[icol] )
            results[theFunc]["nom"].Draw("same")
            leg.AddEntry(results[theFunc]["nom"], theFunc, "L")


            if plotExtra:
                for ivar in range(len(results[theFunc]["vars"])):
                    results[theFunc]["vars"][ivar][0].SetLineColor( R.kGray+2 )
                    results[theFunc]["vars"][ivar][0].Draw("same")
                    results[theFunc]["vars"][ivar][1].SetLineColor( R.kGray+2 )
                    results[theFunc]["vars"][ivar][1].Draw("same")
                    
            

            icol += 1
            
        if plotExtra:
            results_hist_ud["MJ2"]["up0"].SetLineColor( R.kGray+2 )
            leg.AddEntry(results_hist_ud["MJ2"]["up0"], "Param Variations", "L")
        leg.Draw()
        


        c2=R.TCanvas("c2_err_hist_ratio","c2")
        #R.SetOwnership(c,False)
        print "err_hist_ratio",err_hist_ratio
        err_hist_ratio.SetFillColor(R.kBlack)
        err_hist_ratio.SetFillStyle(3001)
        err_hist_ratio.Draw("E2")

        icol = 0
        f_ratio = {}


        delta_ratio_super = {}
        for theFunc in ["Exp","MJ2","MJ3","MJ4","MJ5","MJ6","MJ7","MJ8"]:

            h_ratio  = results[theFunc]["nom"].GetHistogram()
            h_ratio.Divide( nominal_result["nom"] )
            h_ratio.SetDirectory(0)

            h_ratio.SetLineColor( colorlist[icol] )
            h_ratio.Draw("same")

            if plotExtra:
                for ivar in range(len(results[theFunc]["vars"])):
                    h_ratio_ud = results[theFunc]["vars"][ivar][0].GetHistogram()
                    h_ratio_ud.Divide( nominal_result["nom"] )
                    h_ratio_ud.SetDirectory(0)
                    h_ratio_ud.SetLineColor(R.kGray+2)
                    h_ratio_ud.Draw("same")

                    delta_ratio_super[theFunc+"_"+str(ivar)+"_up"] = h_ratio_ud.GetBinContent( h_ratio_ud.FindBin(3000) )


                    h_ratio_ud = results[theFunc]["vars"][ivar][1].GetHistogram()
                    h_ratio_ud.Divide( nominal_result["nom"] )
                    h_ratio_ud.SetDirectory(0)
                    h_ratio_ud.SetLineColor(R.kGray+2)
                    h_ratio_ud.Draw("same")

                    delta_ratio_super[theFunc+"_"+str(ivar)+"_dw"] = h_ratio_ud.GetBinContent( h_ratio_ud.FindBin(3000) )
            
            #print f_copy, f_ratio[theFunc], f_ratio[theFunc].Eval(1000), f_ratio[theFunc].Eval(2000), f_ratio[theFunc].Eval(3000)
            icol += 1

        leg.Draw()
        #for drs in delta_ratio_super:
        #    print drs, delta_ratio_super[drs]

        c.SaveAs(namestr + "_" + c.GetName() +  ".pdf")
        c2.SaveAs(namestr + "_" + c2.GetName() +  ".pdf")
        f.WriteTObject(c)
        f.WriteTObject(c2)
        f.Close()

        
    return smoothFuncCompSyst


    
