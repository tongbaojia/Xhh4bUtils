import ROOT as R

import numpy as np

def makeSmoothedMJJPlots( infileName, outfileName ):

    f = R.TFile(infileName,"READ")

    
    qcd = f.Get("qcd_hh").Clone()
    top = f.Get("ttbar_hh").Clone()
    bkg = qcd.Clone("bkg_hh")

    #stack qcd on top of top
    bkg.Add(top)
    

    AddErrors( qcd, top, bkg, f)
    

    #make canvas
    c=R.TCanvas()
    R.gPad.SetLogy()
       
    bkg.SetLineColor(R.kBlack)
    bkg.SetFillColor(R.kAzure-9)
    bkg.GetXaxis().SetRangeUser(500, 3500)
    bkg.SetXTitle("m_{JJ} [GeV]")
    bkg.SetYTitle("Events")
    bkg.Draw("HIST")

    bkg_err = bkg.Clone("bkg_err")
    bkg_err.SetFillColor(R.kBlack)
    bkg_err.SetFillStyle(3001)
    bkg_err.Draw("sameE2")
    
    top.SetLineColor(R.kBlack)
    top.SetFillColor(R.kRed)
    top.Draw("HISTsame")

    top_err = top.Clone("top_err")
    top_err.SetFillColor(R.kBlack)
    top_err.SetFillStyle(3001)
    top_err.Draw("sameE2")

    leg = R.TLegend(0.1,0.7,0.48,0.9)
    leg.SetFillColor(0)
    leg.AddEntry(bkg, "QCD", "F")
    leg.AddEntry(top, "t #bar{t}", "F")
    leg.Draw()

    c.SaveAs(outfileName)


    return

def AddErrors( qcd, top, bkg, inFile):
    sys_List = ["smoothQ0","smoothQ1","smoothT0","smoothT1", "normY0","normY1","normY2","normY3"]

    NotFound = inFile.Get("NotAKey")

    qcd_errors = [0] + [qcd.GetBinError(ibin) for ibin in range(1, qcd.GetNbinsX()+1)]
    top_errors = [0] + [top.GetBinError(ibin) for ibin in range(1, qcd.GetNbinsX()+1)]
    bkg_errors = [0] + [bkg.GetBinError(ibin) for ibin in range(1, qcd.GetNbinsX()+1)]

    for sys in sys_List:
        qcd_sys_up = inFile.Get("qcd_hh_"+sys+"Up")
        qcd_sys_dw = inFile.Get("qcd_hh_"+sys+"Down")
        top_sys_up = inFile.Get("ttbar_hh_"+sys+"Up")
        top_sys_dw = inFile.Get("ttbar_hh_"+sys+"Down")

        if qcd_sys_up == NotFound or qcd_sys_dw == NotFound:
            qcd_sys_up = qcd.Clone("qcd_hh_"+sys+"__Up")
            qcd_sys_dw = qcd.Clone("qcd_hh_"+sys+"__Down")

        if top_sys_up == NotFound or top_sys_dw == NotFound:
            top_sys_up = top.Clone("ttbar_hh_"+sys+"__Up")
            top_sys_dw = top.Clone("ttbar_hh_"+sys+"__Down")

            
        for ibin in range(1, qcd.GetNbinsX()+1):
            qcd_v = qcd.GetBinContent(ibin)
            top_v = top.GetBinContent(ibin)
            bkg_v = qcd.GetBinContent(ibin)

            
            e_up = qcd_sys_up.GetBinContent(ibin) - qcd_v
            e_dw = qcd_sys_dw.GetBinContent(ibin) - qcd_v
        
            qcd_sys_e = np.max(np.abs( [e_up, e_dw] ))
            qcd_errors[ibin] = np.sqrt(qcd_errors[ibin]**2 + qcd_sys_e**2)

            

            e_up = top_sys_up.GetBinContent(ibin) - top_v
            e_dw = top_sys_dw.GetBinContent(ibin) - top_v

            top_sys_e = np.max(np.abs( [e_up, e_dw] ))
            top_errors[ibin] = np.sqrt(top_errors[ibin]**2 + top_sys_e**2)

            

            e_up =  (qcd_sys_up.GetBinContent(ibin) - qcd_v) + (top_sys_up.GetBinContent(ibin) - top_v)
            e_dw =  (qcd_sys_dw.GetBinContent(ibin) - qcd_v) + (top_sys_dw.GetBinContent(ibin) - top_v)

            bkg_sys_e = np.max(np.abs( [e_up, e_dw] ))
            bkg_errors[ibin] = np.sqrt(bkg_errors[ibin]**2 + bkg_sys_e**2)


            
    for ibin in range(1, qcd.GetNbinsX()+1):
        qcd.SetBinError(ibin, qcd_errors[ibin])
        top.SetBinError(ibin, top_errors[ibin])
        bkg.SetBinError(ibin, bkg_errors[ibin])


    return qcd, top, bkg







if __name__=="__main__":
    makeSmoothedMJJPlots("../BkgFit/outfile_boosted_43.root", "./smooth_result_43.root")
    
