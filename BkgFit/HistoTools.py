
def HistLocationString(distName="DiJetMass", ntrackjet = "4", nbtag="4", btagWP="77", massRegion="SB", whichFunc="SLAC"):
    output = None
    
    if whichFunc == "SLAC":
        output = "GoodEvent_Pass" + ntrackjet + "GoodTrackJetPass" + nbtag + "b" + btagWP +"Pass"+massRegion+"Mass/"+distName
        #folder = lambda nt, nb, wp: "GoodEvent_Pass" + nt + "GoodTrackJetPass" + nb + "b" + wp +"PassSRMass/"
    elif whichFunc == "XhhBoosted":
        word_dict = {"4": "Four", "3" : "Three", "2" : "Two"}
        output = word_dict[nbtag] + "Tag_" + massRegion + "/" + distName

    return output


def CollectHistos(datafile, topfile, distName="DiJetMass", massRegion="SB", btagWP="77", tagRegions=["44","43","42","33","32"]):
    histos = {}
    
    for r in tagRegions:
        folder_r = HistLocStr(distName, r[0], r[1], btagWP, massRegion)  #folder( r[0], r[1], btag_WP)
        
        data_r   = datafile.Get(folder_r).Clone("data_"+r)
        data_r.SetDirectory(0)
        
        top_r    = topfile.Get(folder_r).Clone("top_"+r)
        top_r.SetDirectory(0)     

        histos[r]     = {"data": data_r,            "top": top_r}

    return histos

    
