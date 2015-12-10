import ROOT as R

import numpy as np
from array import array

from GetEigenVariations import GetEigenVariations


def smoothfit(histo, fitFunction = "Exp", fitRange = (900, 3000), makePlots = False, verbose = False, outfileName="fit.root"):
    npar = None
    func = None
    fitChoice = None
    colorlist = [R.kBlue, R.kGreen, R.kOrange, R.kMagenta, R.kCyan, R.kPink]

    fitName = "fit_"+outfileName[:-5]
    
    if fitFunction == "Exp":
        npar = 2
        fitChoice = ExpoFunc
        func = R.TF1(fitName, fitChoice, fitRange[0], fitRange[1], 2)
        func.SetParameters(0.006, 5.0)

    elif fitFunction == "Dijet":
        npar = 3
        fitChoice = DijetFunc
        func = R.TF1(fitName, fitChoice, fitRange[0], fitRange[1], 3)
        func.SetParameters(0.3, 30, -3)

    Vmode = ("Q" if not verbose else "")
    fitResult = histo.Fit(fitName, "S0"+Vmode, "", fitRange[0], fitRange[1])

    if fitResult.Status() != 0:
        print "Error in smoothing fit: did not terminate properly. Exiting"
        sys.exit(0)
    
    cov_TMatrix = fitResult.GetCovarianceMatrix()
    cov = np.zeros( (npar, npar) )
    for i in range(npar):
        for j in range(npar):
            cov[i,j] = cov_TMatrix[i][j]

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
    
    return {"nom": drawFunc, "vars":fvar}


def ExpoFunc(x, par):
    return np.exp(-par[0]*x[0] + par[1])

def DijetFunc(x, par):
    return par[0] * np.power((1.0 - x[0] / 13000.0), par[1]) * np.power((x[0] / 13000.0), par[2])



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
    h = datafile.Get("GoodEvent_Pass4GoodTrackJetPass2b80PassSRMass/DiJetMass").Clone()
    smoothfit( h )
    
