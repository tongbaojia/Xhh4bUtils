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
        func.SetParameters(1, 10, 1)

    elif fitFunction == "MJ3":
        npar = 3
        fitChoice = MJ3Func
        func = R.TF1(fitName, fitChoice, fitRange[0], fitRange[1], npar)
        func.SetParameters(1, 10, 1)

    elif fitFunction == "MJ4":
        npar = 3
        fitChoice = MJ4Func
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

    

    drawFunc = R.TF1("drawfit", fitChoice, fitRange[0], 5000, npar)
    drawFunc.SetParameters( params )

    if makePlots:
        c=R.TCanvas()
        #R.SetOwnership(c,False)
        histo.Draw()
        drawFunc.Draw("same")

    fvar = []
    for ivar in range(len(z_variations)):
        fup = R.TF1("fup_"+str(ivar), fitChoice, fitRange[0], 5000, npar)
        fup.SetParameters( z_variations[ivar][0] )
        fup.SetLineColor(colorlist[ivar])
        if makePlots:
            fup.Draw("same")

        fdw = R.TF1("fdw_"+str(ivar), fitChoice, fitRange[0], 5000, npar)
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


if __name__=="__main__":
    datafile = R.TFile("hist_data.root ","READ")
    h = datafile.Get("GoodEvent_Pass4GoodTrackJetPass2b77PassSRMass/DiJetMass").Clone()
    smoothfit( h, fitFunction = "Dijet", fitRange = (900, 3000), makePlots = True, verbose = True )
    
