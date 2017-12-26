import pandas as pd
from astropy.coordinates import SkyCoord
from astropy import units as u
import numpy as np
import matplotlib.pyplot as plt
from astropy.cosmology import Planck15
from astropy import units as u
from scipy.odr import Model, Data, ODR
from scipy.stats import linregress
import scipy.stats.distributions as dist
import scipy

def latex_float(f):
    string =[]
    for i in f:
        float_str = "{0:.2g}".format(i)
        if "e" in float_str:
            base, exponent = float_str.split("e")
            string = string + [r"{0} \times 10^{{{1}}}".format(base, int(exponent))]
        else:
            string = string + [float_str]
    return string

def latex_float_single(f):
    string =[]
    float_str = "{0:.2g}".format(f)
    if "e" in float_str:
        base, exponent = float_str.split("e")
        string = r"{0} \times 10^{{{1}}}".format(base, int(exponent))
    else:
        string = float_str
    return string

def latex_float_single_e1(f):
    string =[]
    float_str = "{0:.1g}".format(f)
    if "e" in float_str:
        base, exponent = float_str.split("e")
        string = r"{0} \times 10^{{{1}}}".format(base, int(exponent))
    else:
        string = float_str
    return string

def latex_float_err(f,f_err):
    string =[]
    for i in range(len(f)):
        float_str = "{0:.2g}".format(f[i])
        float_err = "{0:.2g}".format(f_err[i])
        if "e" in float_str:
            base, exponent = float_str.split("e")
            base_err, exponent_err = float_err.split("e")
            if int(exponent) == int(exponent_err):
                string = string + [r"$(%s\pm%s) \times 10^{{%s}}$" % (base,base_err, int(exponent))]
            else: 
                x = int(exponent)-int(exponent_err)
                string = string + [r"$(%s\pm %.1f) \times 10^{{%s}}$" % (base,float(base_err)/(10**float(x)), int(exponent))]
        else:
            string = string + [float_str]
    return string

def latex_float_err_res(f,f_err,res):
    string =[]
    for i in range(len(f)):
        float_str = "{0:.2g}".format(f[i])
        float_err = "{0:.2g}".format(f_err[i])
        if "e" in float_str:
            base, exponent = float_str.split("e")
            base_err, exponent_err = float_err.split("e")
            if int(exponent) == int(exponent_err):
                string = string + [r"$(%s\pm%s) \times 10^{{%s}}$" % (base,base_err, int(exponent))]
            else: 
                x = int(exponent)-int(exponent_err)
                string = string + [r"$(%s\pm %.1f) \times 10^{{%s}}$" % (base,float(base_err)/(10**float(x)), int(exponent))]
        else:
            string = string + [float_str]
    for i in range(len(string)):
        if res[i] == 1:
            string[i] = '$>$'+string[i]
            
    return string


def match_catalogues(cat1, cat2, name1, name2, RA1, Dec1, unit1, RA2, Dec2, unit2, columns1, columns2, distance,keep_cat1,keep_cat2):
    ##Check for duplicates
    where = np.where(pd.concat([pd.DataFrame(cat1.columns),pd.DataFrame(cat2.columns)]).reset_index(drop=True).duplicated()==True)[0]-len(pd.DataFrame(cat1.columns))
    print columns2
    for i in where:
        x = cat2.columns[i]
        cat2 = cat2.rename(index=str, columns={x: x+'_1'})
        if x in columns2:
            columns2[np.where(columns2==x)[0][0]] = x+'_1' ## Replace duplicates in columns 2
            if x.startswith(RA2):
                RA2 = x+'_1'
            if x.startswith(Dec2):
                Dec2 = x+'_1'
        #print 'Duplicates:%s' % cat2.columns[i]
    coo_cat1 = SkyCoord(cat1[RA1], cat1[Dec1],unit=unit1)
    coo_cat2 = SkyCoord(cat2[RA2], cat2[Dec2],unit=unit2)
    idxc, idxcatalog, d2d, d3d = coo_cat1.search_around_sky(coo_cat2, distance)
    frames = []
    for i in columns1:
        frames = frames + [cat1[i].iloc[idxcatalog].reset_index(drop=True)]
    for i in columns2:
        frames = frames + [cat2[i].iloc[idxc].reset_index(drop=True)]
    plt.clf()
    plt.figure(1)
    plt.hist(d2d.arcsec, histtype='step', range=(0,distance.value))
    plt.xlabel('separation [%s]'% distance.unit)
    plt.title('On-sky separation')
    plt.savefig('separation_%s_%s' % (name1,name2))
    plt.show()
    #plt.figure(2)
    #plt.hist(d2d.arcsec, histtype='step')
    #plt.hist(d2d_ran.arcsec, histtype='step')
    #plt.xlabel('separation [%s]'% distance.unit)
    #plt.title('On-sky separation')
    #plt.tight_layout()
    x = pd.concat(frames,axis=1)
    if keep_cat1 == True:
        x = cat1.merge(x, how='left')
    if keep_cat2 == True:
        x = cat1.merge(x, how='left')
    return x

def nearest_match_to(coordinates,catalog,RA1,DEC1,unit1,RA2,DEC2,unit2,distance):
    coords = SkyCoord(ra=coordinates[RA1], dec=coordinates[DEC1], unit=unit1)  
    catalog_coords = SkyCoord(ra=catalog[RA2], dec=catalog[DEC2], unit=unit2)  
    idx, d2d, d3d = coords.match_to_catalog_sky(catalog_coords)
    separation = pd.DataFrame({'Separation':d2d.arcsec})
    RA, Dec = coords.spherical_offsets_to(catalog_coords[idx])
    DeltaRA = pd.DataFrame({'DeltaRA':RA.arcsec})
    DeltaDEC = pd.DataFrame({'DeltaDec':Dec.arcsec})
    x = pd.concat([catalog.iloc[idx].reset_index(drop=True),coordinates,separation,DeltaRA,DeltaDEC],axis=1)
    return x[x['Separation']<distance].reset_index(drop=True)

def Monte_carlo_d_L_errs(z,z_upp,z_low,plot_path):
    x = []
    f, (ax1, ax2) = plt.subplots(1, 2, sharey=True)
    for i in range(len(z)):
        z_err = np.random.normal(z[i],(z_upp[i]+z_low[i])/2,100000)
        z_err =z_err[z_err > 0]
        d_L = Planck15.luminosity_distance(z_err).to(u.m).value
        ax1.hist(z_err, 100,alpha=0.3)
        ax2.hist(d_L, 100,alpha=0.3)
        print np.sqrt(np.var(d_L))
        x = x + [np.sqrt(np.var(d_L))]
    ax1.set_xlabel('z')
    ax2.set_xlabel(r'$d_{L}$')
    plt.show()
    plt.savefig(plot_path+'d_L_errors.pdf',bbox_inches='tight')
    return np.array(x)

def f(p, x):
    """Basic linear regression 'model' for use with ODR"""
    return (p[0] * x) + p[1]

def orthoregress(x, y):
    """Perform an Orthogonal Distance Regression on the given data,
    using the same interface as the standard scipy.stats.linregress function.
    Arguments:
    x: x data
    y: y data
    Returns:
    [m, c, nan, nan, nan]
    Uses standard ordinary least squares to estimate the starting parameters
    then uses the scipy.odr interface to the ODRPACK Fortran code to do the
    orthogonal distance calculations.
    """
    linreg = linregress(x, y)
    mod = Model(f)
    dat = Data(x, y)
    od = ODR(dat, mod, beta0=linreg[0:2])
    out = od.run()
    return list(out.beta) + [np.nan, np.nan, np.nan]

def resolver(fitted_size, beam):
    if fitted_size < beam:
        return '$<$%.1f' % beam
    else: 
        return '%.1f' % fitted_size
    
def resolver_val(fitted_size, beam):
    if fitted_size < beam:
        return beam
    else: 
        return fitted_size
    
def minimum_resolvable_size(SNR,weighting,axes):
    return (2.**(2.-(weighting/2.)))*axes*np.sqrt((np.log(2)/np.pi)*np.log(SNR/(SNR-1)))

def variability(epoch1flux,epoch2flux,epoch1_err,epoch2_err):
    ## First pick sources with modulation index > 4.3
    m = []
    variable = (epoch1flux-epoch2flux)/np.sqrt(epoch1_err**2 + epoch2_err**2)
    for i in range(len(variable)):
        if np.abs(variable[i]) >= 4.3:
            m = m+[(2*(epoch1flux[i]/epoch2flux[i])-1)/((epoch1flux[i]/epoch2flux[i])+1)]
        else:
            m = m +[np.nan]
    return np.array(variable), np.array(m)

def Bayesian_binomial_confidence_interval(n,p,sigma):
    ## From Cameron+2011
    c = scipy.special.erf(float(sigma)/np.sqrt(2)) ## create confidence interval based upon sigma
    k = p
    n = n
    p_lower = dist.beta.ppf((1-c)/2.,k+1,n-k+1) 
    p_upper = dist.beta.ppf(1-(1-c)/2.,k+1,n-k+1)
    return p_lower, p_upper,k/n, (k/n)-p_lower, p_upper-(k/n)
        
def convert_columns(catalogue,suffix):
    for i in catalogue.columns:
        catalogue = catalogue.rename(index=str, columns={i: '%s_%s' % (i,suffix)})
    return catalogue

def convert_degrees_to_hmsdms(ra, dec, roundingRA, roundingDec, form):
    c = SkyCoord(ra, dec, unit=('deg', 'deg'))
    if form == 'hmsdms':
        Coords = ['%s:%s:%s' % (str(int(c.ra.hms[0])).rjust(2, '0'),str(int(c.ra.hms[1])).rjust(2, '0'),\
                                '.'.join([str(c.ra.hms[2]).split('.')[0].rjust(2, '0'),\
     str(np.round(float(str(c.ra.hms[2])),roundingRA)).split('.')[1].ljust(roundingRA,'0')])),\
          '+%s:%s:%s' % (str(int(c.dec.dms[0])).rjust(2,'0'),str(int(c.dec.dms[1])).rjust(2, '0'),\
                        '.'.join([str(c.dec.dms[2]).split('.')[0].rjust(2, '0'),\
                         str(np.round(float(str(c.dec.dms[2])),roundingDec)).split('.')[1].ljust(roundingDec,'0')]))]
    elif form == 'dmsdms':
        Coords = ['%s:%s:%s' % (c.ra.dms[0],c.ra.dms[1],\
                                str(np.round(c.ra.dms[2],roundingRA)).rjust(3+roundingRA, '0').ljust(3+roundingRA, '0')),\
          '%s:%s:%s' % (c.dec.dms[0],c.dec.dms[1],str(np.round(c.dec.dms[2],roundingDec)).rjust(3+roundingDec, '0').ljust(3+roundingDec, '0'))]
    else:
        print 'Idiot, this converts from deg to \'hmsdms\' or \'dmsdms\' formats'
    return Coords