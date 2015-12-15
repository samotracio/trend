Trendplot : XY scatterplot with trend curves
============================================
Do XY scatterplots with trend curves and optionally plot errorbars + sigma curves.
If desired, you can skip the plot and just return coordinates of the trend curve.

The main function is **trendplot()**. See its help for examples.

Other helper functions calculate bootstrapping, mean, std, etc. It is easy to 
define new functions for custom statistics.

Basic Usage
--------
    import trend as t
    x = np.arange(0,100,0.1)
    y = 1.2*x + np.random.randn(len(x))*15
    y[300:600] = 55. + np.random.randn(300)*6

    # Simple plot
    t.trendplot(x,y)

<img src="https://github.com/samotracio/trend/raw/master/ex1.png?raw=true" alt="example1" width="400px">
    
    # More bins, 1-sigma error lines, stdev error bar
    t.trendplot(x,y, nbin=20, sigfact=1, estat='std')

<img src="https://github.com/samotracio/trend/raw/master/ex2.png?raw=true" alt="example2" width="400px">
    
    # Boostrap error bars, custom curve range, custom ystat
    t.trendplot(x,y,prange=[30,90,30,100],crange=[30,80],ystat='mean',estat='boot',ccolor='blue',pcolor='green')

<img src="https://github.com/samotracio/trend/raw/master/ex3.png?raw=true" alt="example3" width="400px">

    # Don't plot but get the biweight mean tendency curve with boostrap errors
    px,py,ey = t.trendplot(x,y,ystat='bimean',estat='boot',noplot=True)

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