
def HistLocationString(distName="DiJetMass", ntrackjet = "4", nbtag="4", btagWP="77", massRegion="SB", whichFunc="2016"):
    output = None
    
    if whichFunc == "SLAC":
        output = "GoodEvent_Pass" + ntrackjet + "GoodTrackJetPass" + nbtag + "b" + btagWP +"Pass"+massRegion+"Mass/"+distName
        #folder = lambda nt, nb, wp: "GoodEvent_Pass" + nt + "GoodTrackJetPass" + nb + "b" + wp +"PassSRMass/"
    
    #this is the setting!!!
    elif whichFunc == "XhhBoosted":
        word_dict = {"4": "FourTag", "3" : "ThreeTag", "2" : \
        "TwoTag", "2s" : "TwoTag_split", "1": "OneTag", "0": "NoTag"}
        bkg_dict = {"0":"", "1":"", "2" : "_2Trk_split_", "3" : "_3Trk_", "4" : "_4Trk_"} ##if put in specific numbers, then 
        if (ntrackjet != "i") and (nbtag == "s"): #this is for the specific bcg modelingl without "i"; will look for NoTag_blah, thus "s"
          output = word_dict["0"] + bkg_dict[ntrackjet] + massRegion + "/" + distName
        else:
          output = word_dict[nbtag] + "_" + massRegion + "/" + distName

    #this is not my setting!
    elif whichFunc =="2016":
        tag_dict = {"4": "Four", "3" : "Three", "2" : "Two", "1":"One", "0":"No"}
        reg_dict = {"SB":"Sideband", "CR":"Control", "SR":"Signal"}
        if nbtag == "3" or nbtag == "4" or nbtag == "1" :
            output = tag_dict[nbtag] + "Tag_" +  reg_dict[massRegion] + "/" + distName
        elif nbtag == "2":
            output = tag_dict[nbtag] + "Tag_split_" +  reg_dict[massRegion] + "/" + distName
        elif nbtag == "0":
            output = tag_dict[nbtag] + "Tag_" + ntrackjet + "Trk_"+ ("split_" if ntrackjet=="2" else "") + reg_dict[massRegion] + "/" + distName

    return output


def CollectHistos(datafile, topfile, distName="DiJetMass", massRegion="SB", btagWP="77", tagRegions=["44","43","42","33","32"]):
    histos = {}
    
    for r in tagRegions:
        folder_r = HistLocStr(distName, r[0], r[1], btagWP, massRegion)  #folder( r[0], r[1], btag_WP)
        
        data_r   = datafile.Get(folder_r).Clone("data_"+r)
        data_r.SetDirectory(0)
        if (r == "42") and (massRegion == "SR"):
          data_r = BlindData2bSR(data_r)
        
        top_r    = topfile.Get(folder_r).Clone("top_"+r)
        top_r.SetDirectory(0)     

        histos[r]     = {"data": data_r,            "top": top_r}

    return histos

def BlindData2bSR(h_data, blindThreshold=2000.):
  h_data_blinded = h_data.Clone()
  h_data_blinded.SetDirectory(0)

  blind_bin = h_data_blinded.GetXaxis().FindFixBin(blindThreshold)
  nbins = h_data_blinded.GetNbinsX()

  for ibin in range(blind_bin, nbins+2):
    h_data_blinded.SetBinContent(ibin, 0)
    h_data_blinded.SetBinError(ibin, 0)

  return h_data_blinded

def CheckAndGet(infile, folder, alternative_histo):
    if infile == None:
        out_hist = alternative_histo.Clone(folder+"__clone")
        out_hist.Reset("ICE")
        out_hist.SetDirectory(0)
        return out_hist
    
    NotFound = infile.Get("NotAKey")
    if infile.Get(folder) == NotFound:
        out_hist = alternative_histo.Clone(folder+"__clone")
        out_hist.Reset("ICE")
    else:
        out_hist = infile.Get(folder)

    out_hist.SetDirectory(0)
    return out_hist
        
        
    
def GetSignalHistos(infileName="hist_RSG_c10.root",
                    distName="ChannelNumber_DiJetMass",
                    tagRegions=["44","33","22","40","30","20"],
                    massRegion="SR", btagWP="77",
                    massPoints=["all"]):

    tagRegions  = np.asarray(tagRegions)
    if tagRegions.shape==():
        tagRegions = np.asarray([tagRegions])

    massPoints  = np.asarray(massPoints)
    if massPoints.shape==():
        massPoints = np.asarray([massPoints])

    if massPoints[0]=="all":
        massPoints = ["300","500","600","700","800","900","1000",
                      "1100","1200","1300","1400","1500","1600","1800",
                      "2000","2250","2500","2750","3000"]


    sfile = TFile(infileName,"READ")

    sigDict={}

    for tr in tagRegions:
        dist_tr = sfile.Get(HistLocationString(distName, tr[0], tr[1], btagWP, massRegion)).Clone(distName+"_"+tr)

        proj_x = dist_tr.ProjectionX()

        mass_dict = {}
        for im in massPoints:
            chNum =  dict_RSG_Mass_channelNumber_c10[ int(im) ]

            xbin = proj_x.FindBin( chNum )
            proj_y = dist_tr.ProjectionY(distName+"_"+tr+"_mass"+im, xbin, xbin)

            mass_dict[im] = proj_y

        sigDict[tr]=mass_dict
        
        
    return sigDict, massPoints


dict_RSG_channelNumber_Mass = {301488: ('c10', 300),
                               301490: ('c10', 500),
                               301491: ('c10', 600),
                               301492: ('c10', 700),
                               301493: ('c10', 800),
                               301494: ('c10', 900),
                               301495: ('c10', 1000),
                               301496: ('c10', 1100),
                               301497: ('c10', 1200),
                               301498: ('c10', 1300),
                               301499: ('c10', 1400),
                               301500: ('c10', 1500),
                               301501: ('c10', 1600),
                               301502: ('c10', 1800),
                               301503: ('c10', 2000),
                               301504: ('c10', 2250),
                               301505: ('c10', 2500),
                               301506: ('c10', 2750),
                               301507: ('c10', 3000) }

dict_RSG_Mass_channelNumber_c10 = { 300  : 301488,
                                    500  : 301490,
                                    600  : 301491,
                                    700  : 301492,
                                    800  : 301493,
                                    900  : 301494,
                                    1000 : 301495,
                                    1100 : 301496,
                                    1200 : 301497,
                                    1300 : 301498,
                                    1400 : 301499,
                                    1500 : 301500,
                                    1600 : 301501,
                                    1800 : 301502,
                                    2000 : 301503,
                                    2250 : 301504,
                                    2500 : 301505,
                                    2750 : 301506,
                                    3000 : 301507}



    

    

    
