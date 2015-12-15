# -*- coding: utf-8 -*-
"""
XY scatterplot with trend curves
--------------------------------
Do XY scatterplots with trend curves and optionally plot errorbars + sigma curves.
If desired, you can skip the plot and just return coordinates of the trend curve.

The main function is **trendplot()**. See its help for examples.

Other helper functions calculate bootstrapping, mean, std, etc. It is easy to 
define new functions for custom statistics.

Dependencies
------------
1. astropy : to provide some statistics
2. matplotlib : to provide basic graphics

Todo
----
- [ ] Pass extra keyword parameters to ``plt.plot`` routine

Credits
-------
Emilio Donoso

@author: edonoso
"""

import numpy as np
import matplotlib.pyplot as plt
from astropy.stats import bootstrap,biweight_location,biweight_midvariance

# Available y-axis statistics =================================================
def fmed(y):
    """Helper function of trendplot. Returns median of input 1d-array"""
    return np.median(y)
    
def fmean(y):
    """Helper function of trendplot. Returns mean of input 1d-array"""
    return np.mean(y)

def fbimean(y):
    """Helper function of trendplot. Returns biweight mean of input 1d-array (see astropy.stats)"""
    return biweight_location(y)


# Available y-axis error bars =================================================
def fstd(y):
    """Helper function of trendplot. Returns std of input 1d-array"""
    return np.std(y)

def fpoisson(y):
    """Helper function of trendplot. Returns poisson error (srtq(N)) of input 1d-array"""
    return np.sqrt(len(y))

def fboot(y):
    """
    Helper function of trendplot. Returns bootstrap error of input 1d-array.

    Calculates the biweight mean of each bootstrap sample, and then gets the
    biweight standard deviation of all samples.
    """
    if len(y)>5:
        bb=bootstrap(y,bootnum=250,bootfunc=biweight_location)
        res=biweight_midvariance(bb)
    else:
        res=np.float('nan')
    return res

def fbistd(y):
    """Helper function of trendplot. Returns biweight std of input 1d-array (see astropy.stats)"""
    return biweight_midvariance(y)


# Main function ===============================================================
def trendplot(xin,yin,crange=None,nbin=None,xycut=None,prange=None,
                  ystat=None,estat=None,psize=None,pcolor=None,
                  cmarker=None,ccolor=None,csize=None,cls=None,clw=None,
                  sigfact=None,sigcolor=None,sigls=None,
                  noscatter=False,noplot=False):
    """
    XY scatter plot with tendency curves and optional errorbars + sigma curves.
    
    Parameters
    ----------
    **xin,yin** : 1d-arrays.
       Input data. Required.

    **crange** : list of form [xlim1,xlim2].
       x-axis range that the tendency curve will span. Defaults to data [min,max].
       
    **nbin** : int.
       Nr of bins of tendency curve.
       
    **xycut** : list of form [xlim1,xlim2,ylim1,ylim2].
       Rectangular window to cut input data by. Defaults to data [min,max].

    **prange** : list of form [xlim1,xlim2,ylim1,ylim2].
       Rectangular window for plotting. Defaults to data [min,max].
       
    **ystat** : y-axis statistic. Default: "median".
       * "median" : median of data in bin.

       * "mean" :  mean of data in bin.
       
       * "bimean" : biweight mean of data in bin (from astropy.stats).

    **estat** : y-axis error bar statistic.
       * "std" : standard deviation.
       
       * "bistd" : biweight standard deviation(from astropy.stats).
       
       * "poisson" : poisson error.
       
       * "boot" : boostrap sampling error, i.e. the error of each bin is bistd(bimean(bootsamples)). Defaults to 250 samples.
    
    **psize** : point size of scatter plot. Default: 1
    
    **pcolor** : point color of scatter plot. Default: "0.5" (gray)
    
    **cmarker** : marker of tendency curve. Default: "o"
    
    **ccolor** : color of tendency curve. Default: "red"
    
    **csize** : size of tendency curve markers. Default: 4
    
    **cls** : linestyle of tendency curve. Default: "-"
    
    **clw** : line width of tendency curve. Default: 1
    
    **sigfact** : float. Plot 1,2,3 or x sigma lines above and below tendency curve. Uses bistd().
    
    **sigcolor** : color of sigma lines. Default: "k"
    
    **sigls** :  linestyle of sigma lines. Default: "--"
    
    **noscatter** : bool. Do not plot scatter points. Default: False
    
    **noplot** : bool. Do not plot anything and return curves instead. Default: False

    Returns
    -------
    (only when noplot=True)
    
    * **(px,py)** : coordinates of curve, if estat is not set
    
    * **(px,py,ey)**: coordinates of curve + errors, if estat is set 
    
    Examples
    --------
    ::
    
        import trend as t
        x = np.arange(0,100,0.1)
        y = 1.2*x + np.random.randn(len(x))*15
        y[300:600] = 55. + np.random.randn(300)*6

        # Simple plot
        t.trendplot(x,y)

        # More bins, 1-sigma error lines
        t.trendplot(x,y, nbin=20, sigfact=1, estat='std')

        # Custom curve range, custom ystat, boostrap error bars
        t.trendplot(x,y,crange=[30,80],ystat='mean',estat='boot',ccolor='blue')

        # Don't plot but get the biweight mean tendency curve with boostrap errors
        px,py,ey = t.trendplot(x,y,ystat='bimean',estat='boot',noplot=True)
    """

    # Default plotting styles
    if psize is None:    psize    = 1       # scatter point size
    if pcolor is None:   pcolor   = '0.5'   # scatter point color
    if cmarker is None:  cmarker  = 'o'     # t. curve point style
    if csize is None :   csize    = 4       # t. curve point size
    if ccolor is None:   ccolor   = 'red'   # t. curve color
    if cls is None:      cls      = '-'     # t. curve line style
    if clw is None:      clw      = 1       # t. curve line width
    if sigcolor is None: sigcolor = 'k'     # sigma line color
    if sigls is None:    sigls    = '--'    # sigma line style
    
    # Default number of bins
    if nbin is None:  nbin=10

    # Default y-axis statistic
    if ystat is None:  ystat='median'

    # Default xrange that the curve will spawn
    if crange is None:  crange=[np.min(xin),np.max(xin)]
    
    # Cut input values. Add some tolerance
    if xycut is None:
        tol=0.01
        xycut=[0.,0.,0.,0.]
        xycut[0]=np.min(xin)-tol
        xycut[1]=np.max(xin)+tol
        xycut[2]=np.min(yin)-tol
        xycut[3]=np.max(yin)+tol
    idxin=np.where((xin>xycut[0]) & (xin<xycut[1]) & (yin>xycut[2]) & (yin<xycut[3]))
    xin, yin = xin[idxin], yin[idxin]

    # Default plot range for entire plot. Do this after cutting input?
    if prange is None:
        prange=[np.min(xin),np.max(xin),np.min(yin),np.max(yin)]
       
    # Bin x-coordinates
    xcounts,xedges=np.histogram(xin,bins=nbin,range=crange)

    # Find curve x-axis (mid-point of bins)
    hbsiz=abs(xedges[0]-xedges[1])*0.5
    px=xedges[:-1]+hbsiz   #drop last

    # Find curve y-axis (mean, median, etc.)
    py=np.zeros(nbin)
    ey=np.zeros(nbin)
    sigmav=np.zeros(nbin)
    ystatoptions={'median':fmed,'mean':fmean,'bimean':fbimean}
    estatoptions={'std':fstd,'poisson':fpoisson,'boot':fboot,'bistd':fbistd}
    for l,r,i in zip(xedges[:-1],xedges[1:],range(nbin)):
        ydat=yin[(xin>l) & (xin<r)]
        py[i]=ystatoptions[ystat](ydat)
        if estat is not None:
            ey[i]=estatoptions[estat](ydat)
        if sigfact is not None:
            sigmav[i]=biweight_midvariance(ydat)
        #print(i,'--',l,'--',r,'--',elm,'--',py[i])
    
    # Do plot or get results  =================================================
    if noplot==False:
        if noscatter==False:
            plt.scatter(xin,yin,marker='.',s=psize,color=pcolor)  # scatter points
        print cls
        plt.plot(px,py,color=ccolor,linestyle=cls,linewidth=clw)  # curve line
        plt.plot(px,py,marker=cmarker,markersize=csize,color=ccolor,markeredgecolor=ccolor,linestyle='none')  #curve points
        if estat is not None:
            plt.errorbar(px,py,yerr=ey,fmt="none",ecolor=ccolor)  # curve errorbars
        if sigfact is not None:
            plt.plot(px,py+sigfact*sigmav,color=sigcolor,linestyle=sigls)  #sigma lines
            plt.plot(px,py-sigfact*sigmav,color=sigcolor,linestyle=sigls)
        # Reset plot limits to desired range
        plt.xlim(prange[0],prange[1])
        plt.ylim(prange[2],prange[3])
    else:
        if estat is not None:
            res=(px,py,ey)
        else:
            res=(px,py)
        return res
