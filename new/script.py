import numpy as np
import os
import matplotlib.pyplot as plt

import mahotas as mh

myPath = "/Volumes/GoogleDrive/My Drive/Janes Lab/Projects/Mouse glioma/Data/images/Ki76-pRB/"

f = sorted(os.listdir(myPath))

for iFile in f:
    I = plt.imread(myPath + iFile)
    I = I / 2**12
    plt.imsave(myPath + "python_" + iFile[:-3] + "png", I, cmap = "gray")
