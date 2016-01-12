import ROOT as R

import numpy as np
from array import array
import sys

from GetEigenVariations import GetEigenVariations

import cPickle as pickle

def smoothfit(histo, fitFunction = "Exp", fitRange = (900, 3000), makePlots = False, verbose = False, outfileName="fit.root"):
    npar = None
    func = None
    fitChoice = None
    colorlist = [R.kBlue, R.kGreen, R.kOrange, R.kMagenta, R.kCyan, R.kPink]

    fitName = "fit_"+outfileName[:-5]
    
    if fitFunction == "Exp":
        npar = 2
        fitChoice = ExpoFunc
        func = R.TF1(fitName, fitChoice, fitRange[0], fitRange[1], npar)
        func.SetParameters(0.006, 5.0)

    elif fitFunction == "Dijet":
        npar = 3
        fitChoice = DijetFunc
        func = R.TF1(fitName, fitChoice, fitRange[0], fitRange[1], npar)
        func.SetParameters(1, 10, 1)

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

    elif fitFunction == "MJ8":
        npar = 3
        fitChoice = MJ8Func
        func = R.TF1(fitName, fitChoice, fitRange[0], fitRange[1], npar)
        func.SetParameters(1, 10, 1)

    Vmode = ("Q" if not verbose else "")
    fitResult = histo.Fit(fitName, "S0"+Vmode, "", fitRange[0], fitRange[1])

    if fitResult.Status() != 0:
        print "Error in smoothing fit: did not terminate properly. Exiting"
        c=R.TCanvas()
        #R.SetOwnership(c,False)
        histo.Draw()
        c.SaveAs("failed__"+outfileName)
        sys.exit(0)
            
    
    cov_TMatrix = fitResult.GetCovarianceMatrix()
    cov = np.zeros( (npar, npar) )
    corr = np.zeros( (npar, npar) )
    for i in range(npar):
        for j in range(npar):
            cov[i,j] = cov_TMatrix[i][j]
            corr[i,j] = cov_TMatrix[i][j] / np.sqrt(cov_TMatrix[i][i] * cov_TMatrix[j][j])

    #print "correlation matrix"
    #print corr
    

    S, U= np.linalg.eigh( cov )
    Sd = np.diag( S )
    sigma = np.sqrt(S)

    evars = GetEigenVariations( cov )

    
    fitFunc = histo.GetFunction(fitName)

    params = array('d',[0]*npar)
    fitFunc.GetParameters( params )

    z = np.asarray( params )
    z_variations = []

    for i in range(len(evars)): 
        zu = z + evars[i]
        zd = z - evars[i]
        z_variations.append( [zu, zd] )


    namestr = outfileName.split(".root")[0]

    drawFunc = R.TF1("drawfit_"+namestr, fitChoice, fitRange[0], 5000, npar)
    drawFunc.SetParameters( params )

    if makePlots:
        c=R.TCanvas()
        #R.SetOwnership(c,False)
        histo.Draw()
        drawFunc.Draw("same")

    fvar = []
    for ivar in range(len(z_variations)):
        fup = R.TF1("fup_"+str(ivar)+"_"+namestr, fitChoice, fitRange[0], 5000, npar)
        fup.SetParameters( z_variations[ivar][0] )
        fup.SetLineColor(colorlist[ivar])
        if makePlots:
            fup.Draw("same")

        fdw = R.TF1("fdw_"+str(ivar)+"_"+namestr, fitChoice, fitRange[0], 5000, npar)
        fdw.SetParameters( z_variations[ivar][1] )
        fdw.SetLineColor( colorlist[ivar] )
        if makePlots:
            fdw.Draw("same")

        fvar.append([fup, fdw])
    

    #raw_input()

    if makePlots:
        c.SaveAs(outfileName)

        fTxT = open(outfileName[:-5]+".pkl", "w")
        fitResultDict = {
          "params": params,
          "cov": cov,
        }
        pickle.dump(fitResultDict, fTxT)
        fTxT.close()

    
    return {"nom": drawFunc, "vars":fvar}




def smoothFuncCompare(histo, fitFunction = "Exp", fitRange = (900, 3000), makePlots = False, verbose = False, outfileName="smoothFuncCompare.root"):

    colorlist = [R.kBlue, R.kGreen, R.kOrange, R.kMagenta, R.kCyan, R.kPink, (R.kAzure+1), R.kGray]        

    namestr = outfileName.split(".root")[0]
    
    results = {}
    results_hist = {}
    for theFunc in ["Exp","MJ2","MJ3","MJ4","MJ5","MJ6","MJ7","MJ8"]:
        h_clone = histo.Clone()
        h_clone.SetDirectory(0)
        results[theFunc] = smoothfit(h_clone, fitFunction = theFunc, fitRange = fitRange, makePlots = False, verbose = verbose, outfileName = theFunc+"_"+outfileName)
        results_hist[theFunc] = MakeSmoothHisto(h_clone, results[theFunc]["nom"])


    histo_up = results_hist[fitFunction].Clone(histo.GetName() + "_" + namestr + "_up")
    histo_up.SetDirectory(0)
    histo_dw = results_hist[fitFunction].Clone(histo.GetName() + "_" + namestr + "_dw")
    histo_dw.SetDirectory(0)

    for ibin in range(1, histo.GetNbinsX()+1):
        deltas = []
        for theFunc in ["Exp","MJ2","MJ3","MJ4","MJ5","MJ6","MJ7","MJ8"]:
            deltas.append( np.abs( results_hist[fitFunction].GetBinContent(ibin) - results_hist[theFunc].GetBinContent(ibin) ) )

        theDelta = np.max( deltas )
        histo_up.SetBinContent(ibin, histo_up.GetBinContent(ibin) + theDelta)
        histo_dw.SetBinContent(ibin, histo_dw.GetBinContent(ibin) - theDelta)

    smoothFuncCompSyst = {"up":histo_up, "dw":histo_dw}
        
        

    if makePlots:
        f = R.TFile(outfileName, "RECREATE")
        
        c=R.TCanvas("c1","c1")
        #R.SetOwnership(c,False)
        leg = R.TLegend(0.1,0.7,0.48,0.9)
        leg.SetFillColor(0)
    
        histo.SetLineColor(R.kBlack)
        histo.Draw()
        leg.AddEntry(histo, "Histogram", "L")

        icol = 0
        err_hist_ratio = None
        for theFunc in ["Exp","MJ2","MJ3","MJ4","MJ5","MJ6","MJ7","MJ8"]:

            if theFunc==fitFunction:
                err_hist = MakeSmoothHisto(histo, results[theFunc]["nom"])
                err_hist.SetDirectory(0)

                for ivar in range(len(results[theFunc]["vars"])):
                    err_hist_up = MakeSmoothHisto(histo, results[theFunc]["vars"][ivar][0])
                    err_hist_dw = MakeSmoothHisto(histo, results[theFunc]["vars"][ivar][1])

                    for ibin in range(1, err_hist.GetNbinsX()+1):
                        err_val = np.max( np.abs( [ err_hist.GetBinContent(ibin) - err_hist_up.GetBinContent(ibin), err_hist.GetBinContent(ibin) - err_hist_dw.GetBinContent(ibin)] ) )
                        err_hist.SetBinError(ibin, np.sqrt( err_hist.GetBinError(ibin)**2 + err_val**2) )

                err_hist_ratio = err_hist.Clone()
                err_hist_ratio.SetDirectory(0)
                err_hist_ratio.Divide( results[theFunc]["nom"] )



                err_hist.SetFillColor(R.kBlack)
                err_hist.SetFillStyle(3001)
                err_hist.Draw("sameE3")
                leg.AddEntry(err_hist, "smoothing error", "F")



            #print results[theFunc]["nom"], results[theFunc]["nom"].Eval(1000), results[theFunc]["nom"].Eval(2000), results[theFunc]["nom"].Eval(3000)

            results[theFunc]["nom"].SetLineColor( colorlist[icol] )
            results[theFunc]["nom"].Draw("same")
            leg.AddEntry(results[theFunc]["nom"], theFunc, "L")

            icol += 1

        leg.Draw()



        c2=R.TCanvas("c2","c2")
        #R.SetOwnership(c,False)
        err_hist_ratio.Draw("E2")

        icol = 0
        f_ratio = {}

        for theFunc in ["Exp","MJ2","MJ3","MJ4","MJ5","MJ6","MJ7","MJ8"]:
            ## func_ratio[theFunc] = lambda x: (results[theFunc]["nom"].Eval(x[0]) / results["Exp"]["nom"].Eval(x[0]))
            ## f_ratio[theFunc] = R.TF1(theFunc+"_ratio_"+namestr, func_ratio[theFunc], fitRange[0], 3000, 0)
            ## f_ratio[theFunc].SetLineColor( colorlist[icol] )
            ## f_copy = f_ratio[theFunc].DrawCopy("same")

            h_ratio  = results[theFunc]["nom"].GetHistogram()
            h_ratio.Divide( results["Exp"]["nom"] )
            h_ratio.SetDirectory(0)

            h_ratio.SetLineColor( colorlist[icol] )
            h_ratio.Draw("same")

            #print f_copy, f_ratio[theFunc], f_ratio[theFunc].Eval(1000), f_ratio[theFunc].Eval(2000), f_ratio[theFunc].Eval(3000)
            icol += 1

        leg.Draw()

        f.WriteTObject(c)
        f.WriteTObject(c2)
        f.Close()
        
    return smoothFuncCompSyst



def MakeSmoothHisto(hist, fitCurve):
    low=R.Double(0.0)
    high=R.Double(0.0)
    fitCurve.GetRange(low, high)

    hist_smooth = hist.Clone(hist.GetName()+"__smooth")
    for ibin in range(1, hist_smooth.GetNbinsX()+1):
        if hist_smooth.GetBinCenter(ibin) >= low:
            hist_smooth.SetBinContent(ibin, 0)
            hist_smooth.SetBinError(ibin, 0)
    hist_smooth.Add(fitCurve, 1.0)

    return hist_smooth





############################################################################################
### functions
############################################################################################
def ExpoFunc(x, par):
    return np.exp(-par[0]*x[0] + par[1])

def DijetFunc(x, par):
    z = x[0] / 13000.0
    #return np.exp( par[0] ) * np.power((1.0 - z), par[1]) * np.power(z, par[2])
    return np.exp(par[0] + par[1]* np.log((1.0 - z)) + par[2]* np.log(z)  )

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

############################################################################################

if __name__=="__main__":
    datafile = R.TFile("hist_data.root ","READ")
    h = datafile.Get("GoodEvent_Pass4GoodTrackJetPass2b77PassSRMass/DiJetMass").Clone()
    smoothfit( h, fitFunction = "Dijet", fitRange = (900, 3000), makePlots = True, verbose = True )
    
