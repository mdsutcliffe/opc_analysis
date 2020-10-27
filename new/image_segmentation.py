import os
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
import numpy as np
import skimage
from skimage import data, io, filters, morphology
import scipy as sp
from scipy import ndimage as ndi

myPath = "/Volumes/GoogleDrive/My Drive/Janes Lab/Projects/Mouse glioma/Data/images/Ki76-pRB/"

f = sorted(os.listdir(myPath))

tdt = mpimg.imread(myPath + f[3])
tdt = tdt / 3096
io.imshow(tdt)

tdt2 = filters.gaussian(tdt,sigma = 4)
io.imshow(tdt2 > 0.0525)

tdt3 = morphology.remove_small_objects(tdt2 > 0.0525,1000)
io.imshow(tdt3)

tdt4 = sp.ndimage.morphology.binary_fill_holes(tdt3)
io.imshow(tdt4)

tdt5 = morphology.dilation(tdt4)
io.imshow(tdt5)


# Ki67
ki67 = mpimg.imread(myPath + f[0])
ki67 = np.divide(ki67,np.max(ki67))

io.imshow(ki67)

plt.hist(np.ndarray.flatten(ki67),bins = 100)

tdt6 = np.subtract(np.matrix(tdt5),np.matrix(tdt4))
io.imshow(np.logical_and(tdt5,np.invert(tdt4)))

tdt6, _ = ndi.label(tdt5)
tdt6 / 32
plt.imshow(tdt6 / 32)

t = morphology.skeletonize(tdt5)
plt.imshow(t)

b = mpimg.imread(myPath + f[0]) / 1000
g = mpimg.imread(myPath + f[1]) / 1000
r = mpimg.imread(myPath + f[3]) / 500

rgb = np.dstack((r,g,b))

plt.imshow(rgb)
