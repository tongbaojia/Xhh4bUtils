import numpy as np

def GetEigenVariations( cov ):
    S, U= np.linalg.eigh( cov )
    sigma = np.sqrt(S)

    variations = []

    for i in range(sigma.shape[0]):
        var = np.zeros( sigma.shape )
        var[i] = sigma[i]
        
        var = np.dot( U, var)

        variations.append( var )

    return variations
