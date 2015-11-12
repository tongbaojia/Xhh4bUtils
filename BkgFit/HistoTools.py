
def HistLocationString(distName="DiJetMass", ntrackjet = "4", nbtag="4", btagWP="77", massRegion="SB", whichFunc="SLAC"):
    output = None
    
    if whichFunc=="SLAC":
        output = "GoodEvent_Pass" + ntrackjet + "GoodTrackJetPass" + nbtag + "b" + btagWP +"Pass"+massRegion+"Mass/"+distName

    return output
