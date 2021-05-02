from veusz.plugins import *
import string
import os.path
import zipfile
import xml.dom.minidom as xmlminidom
from pylab import *
from scipy.integrate import trapz
from scipy.interpolate import interp1d

versionstring='ProcSpec plugin v140417'

####################
# Veusz Import Plugin for reading Ocean Optics ProcSpec files
# Note: the ProcSpec files should contain a 'Processed' spectrum
# (e.g. an absorbance spectrum, or a background-corrected emission spectrum)
# Files containing only a raw spectrum will not be read, resulting
# in an error message.
# 
# by Martinus Werts, CNRS, ENS Rennes, France
#
# v140416   first release, 

#####################
#
# mwave class, enables storing abscissa and ordinate data simultaneously
#    in a single object (such as wavelength and intensity for a spectrum)
#    allowing for certain operations to be more easily achieved
#    this module has been copied from the original source file
#    in order to have this plugin as a monolithic file
#
#####################

# new mwave object working!
# to do: thorough test baseline correct

class mwave(object):
    def __init__(self,xdata,ydata):
        if (len(xdata)!=len(ydata)):
            raise Exception('xdata and ydata have different lengths')
        self.x=array(xdata)
        xdif12=self.x[1]-self.x[0]
        self.x_ascending=(xdif12>0)
        self.x_descending=(xdif12<0)
        if self.x_ascending:
            if not (all(self.x[:-1] < self.x[1:])):
                raise Exception('xdata not monotonously ascending')
        else:
            if not (all(self.x[:-1] > self.x[1:])):
                raise Exception('xdata not monotonously descending')
        self.y=array(ydata)
        self.yoriginal=self.y # yoriginal is just a shallow copy initially
        self.yinterp = interp1d(self.x, self.y) # linear interpolate, TO DO: add more options
        # however, can be replaced from main code

    def __call__(self, x):
        # if single argument then return interpolated value
        # TO DO: if two arguments then return copy with region of interest (range)
        return self.yinterp(x)

    def identicalx(self,other):
        if len(self.x)!=len(other.x):
            raise Exception('abscissae incompatible')
            #return False
        elif not (all(self.x == other.x)):
            raise Exception('abscissae incompatible')
            #return False
        else:
            return True

    # simple maths
    # TO DO: invoke interpolator if not 'identicalx'
    # (in which case the Exception in identicalx should not be raised, only a False returned) 
        
    def __add__(self,other):
        if isinstance(other,(int,float)):
            return mwave(self.x,(self.y+other))
        elif self.identicalx(other):
            return mwave(self.x,(self.y+other.y))
        
    def __sub__(self,other):
        if isinstance(other,(int,float)):
            return mwave(self.x,(self.y-other))
        elif self.identicalx(other):
            return mwave(self.x,(self.y-other.y))

    def __mul__(self,other):
        if isinstance(other,(int,float)):
            return mwave(self.x,(self.y*other))
        elif self.identicalx(other):
            return mwave(self.x,(self.y*other.y))

    def __div__(self,other):
        if isinstance(other,(int,float)):
            return mwave(self.x,(self.y/other))
        if self.identicalx(other):
            return mwave(self.x,(self.y/other.y))
            
    def __pow__(self,other):
        if isinstance(other,(int,float)):
            return mwave(self.x,pow(self.y,other))
        elif self.identicalx(other):
            return mwave(self.x,pow(self.y,other.y))

    def __radd__(self,other):
        return self.__add__(other)
    
    def __rsub__(self,other):
        if isinstance(other,(int,float)):
            return mwave(self.x,(other-self.y))
        else:
            raise Exception('mwave.__rsub__ exception')

    def __rmul__(self,other):
        return self.__mul__(other)
        
    def __rdiv__(self,other):
        if isinstance(other,(int,float)):
            return mwave(self.x,(other/self.y))
        else:
            raise Exception('mwave.__rdiv__ exception')
            
    def __rpow__(self,other):
        if isinstance(other,(int,float)):
            return mwave(self.x,pow(other,self.y))
        else:
            raise Exception('mwave.__rpow__ exception')
            
    def __neg__(self):
        return mwave(self.x,-self.y)

    #__add__(self, other)
    #__sub__(self, other)
    #__mul__(self, other)
    #__div__(self, other)
    #__mod__(self, other)
    #__divmod__(self, other)
    #__pow__(self, other[, modulo])
    #__lshift__(self, other)
    #__rshift__(self, other)
    #__and__(self, other)
    #__xor__(self, other)
    #__or__(self, other) 
             
    def closestidx(self,xcoord):
        dist=abs(xcoord-self.x)
        return argmin(dist)

    def idxrange(self,xlim1,xlim2):
        '''return as a 2-tuple the indices that define a range 
that lies between x=xlim1 and x=xlim2'''
        tlim1idx = self.closestidx(xlim1)
        tlim2idx = self.closestidx(xlim2)
        # should always obtain a range idx1...idx2 with idx1<idx2
        # (symmetric behav.)
        # therefore if idx1>idx2 ...
        if (tlim1idx>tlim2idx):
            lim1idx=tlim2idx
            lim2idx=tlim1idx+1 # remember: slice is [a,b>
        else:
            lim1idx=tlim1idx
            lim2idx=tlim2idx+1 # remember: slice is [a,b>
        # TO DO: shift boundaries to ensure that range is always "in between"
        # xlim1,xlim2
        # currently, simply takes smallest distance
        return (lim1idx,lim2idx)
        
    def integral(self,xlim1,xlim2):
        (idx1,idx2)=self.idxrange(xlim1,xlim2)
        roi_x=self.x[idx1:idx2]
        roi_y=self.y[idx1:idx2]
        result=trapz(roi_y,roi_x)
        return result
        
    def average(self,xlim1,xlim2):
        return (1.0/(xlim2-xlim1))*self.integral(xlim1,xlim2)
    
    # these would later be part of a 'spectrum' child clase of 'mwave'   
    def polyfit_coeffs(self,fit_sections,polyorder=0):
        '''returns the coefficients of a polynomial fit of order polyorder
based on the x-intervals (2-tuples) given in the list fit_sections'''
        fitdatax = None
        fitdatay = None
        for interval in fit_sections:
            (idx1,idx2)=self.idxrange(interval[0],interval[1])
            fitdatax_ext = self.x[idx1:idx2]
            fitdatay_ext = self.y[idx1:idx2]
            if fitdatax == None:
                fitdatax = fitdatax_ext
                fitdatay = fitdatay_ext
            else:
                fitdatax = concatenate((fitdatax,fitdatax_ext))
                fitdatay = concatenate((fitdatay,fitdatay_ext))
        fitdata_pol = polyfit(fitdatax,fitdatay,polyorder)
        return fitdata_pol

    def baseline_subtract(self, fit_sections,polyorder=0):
        '''subtract a polynomial baseline from the mwave'''
        bl_coeffs = self.polyfit_coeffs(fit_sections,polyorder)
        baseliney = polyval(bl_coeffs,self.x)
        self.yoriginal=(self.y).copy()
        self.y = self.yoriginal - baseliney
        # may return some result... to be implemented
        
    def findmax(self, windowA, windowB): 
        '''find maximum using parabolic fit
    window1: tuple that indicates absolute x values defining the working window
             this is used to find numerical maximum (step 1)
    window2: tuple that indicates relative x values for parabolic fit
             around numerical maximum (step 2)
    returns a tuple containing maximum-x maximum-y'''
        maxsearchidx = self.idxrange(windowA[0],windowA[1])
        fitcenter = self.x[(self.y[maxsearchidx[0]:maxsearchidx[1]].argmax() + maxsearchidx[0])]
        (a,b,c)=self.polyfit_coeffs([((fitcenter+windowB[0]),(fitcenter+windowB[1]))],polyorder=2)
        maxpos = -b/(2*a)
        (a,b,c)=self.polyfit_coeffs([((maxpos+windowB[0]),(maxpos+windowB[1]))],polyorder=2)
        maxpos = -b/(2*a)
        (a,b,c)=self.polyfit_coeffs([((maxpos+windowB[0]),(maxpos+windowB[1]))],polyorder=2)
        maxpos = -b/(2*a)
        maxval = polyval((a,b,c),maxpos)
        return (maxpos,maxval)


#####################
#
# ProcSpec file reader code and mwave wrapper
#    this module has been copied from the original source file
#    in order to have this plugin as a monolithic file
#
#####################

def _readProcSpec(pathname,getAcqTime=False):
    '''read OceanOptics .ProcSpec processed spectrum file
    returns a 2-tuple of lists of wavelengths and value (as lists of floats)'''
    # v06: now returns lists of float
    # v05: added new option: getAcqTime
    # if set True, then the return tuple contains the 'millitime' as third element
    OOzip = zipfile.ZipFile(pathname,mode='r')
    OOzipfiles = OOzip.namelist()
    if len(OOzipfiles) != 3:
        raise Exception('ProcSpec container with more than 3 files')
    OOspecfile = None
    # the processed spectrum file is probably always called ps_<xxx>.xml
    for zf in OOzipfiles:
        if zf.startswith('ps_'):
            OOspecfile=zf
    OOspecstr = OOzip.read(OOspecfile)
    OOzip.close()
    # replace all non-standard ASCII chars in the string
    # (else minidom parser chokes)
    # this problem has to do with the OO XML being encoded
    # in cp1252 and xmlminidom not accepting this
    # xmlminidom seems to be UNICODE! and we should be able to codec this 
    # in a later version
    # currently this is a work-around, which simply replaces
    # each byte in the string with an X for values > 127
    ttblA = string.maketrans('','')
    ttblB = 'X'*128
    ttbl = ttblA[0:128]+ttblB
    OOspecstr_cnv = OOspecstr.translate(ttbl)
    OOdoc = xmlminidom.parseString(OOspecstr_cnv)
    ndprocspec = OOdoc.childNodes[0]
    if ndprocspec.nodeName != 'com.oceanoptics.spectrasuite.processors.ProcessedSpectrum':
        raise Exception('Unexpected Node')
    procpix = ndprocspec.getElementsByTagName('processedPixels')
    # in case of missing 'processedPixels', try to get 'pixelValues'
    # uncomment following line, comment first 'procpix = ...' line
    # procpix = ndprocspec.getElementsByTagName('pixelValues')
    valnd=None
    for nd in procpix:
        if nd.hasChildNodes():
            valnd=nd
    chanwl = ndprocspec.getElementsByTagName('channelWavelengths')
    wvlnnd=None
    for nd in chanwl:
        if nd.hasChildNodes():
            wvlnnd=nd
    values=[]
    for nd in valnd.getElementsByTagName('double'):
        values.append(float(nd.firstChild.nodeValue))
    wavelengths=[]
    for nd in wvlnnd.getElementsByTagName('double'):
        wavelengths.append(float(nd.firstChild.nodeValue))
    acqtime_elems = ndprocspec.getElementsByTagName('acquisitionTime')
    milliel = acqtime_elems[0].getElementsByTagName('milliTime') 
    # apparently the first acquistiontime in the file contains the actual spectrum acq time
    # the others are for dark, and blank apparently
    acqmillitime = milliel[0].firstChild.nodeValue
    # print acqmillitime
    if getAcqTime:
        # print >> sys.stderr, 'Warning! I hope I extracted the right aquisition timestamp'
        # we should aim to find a more reliable way of extracting the right acquisitionTime
        returntuple = (wavelengths,values,acqmillitime)
    else:
        returntuple = (wavelengths,values)
    return returntuple

def mwave_loadProcSpec(pathname,xcoordinates=None):
    # TO DO:xmwave can contain abscissa positions, will load and interpolate   
#    # check if pkl exists, load from pkl
#    pklpathname=pathname+'.pkl'
#    if os.path.exists(pklpathname):
#        # load from pkl
#        # print 'fast load from Pickle'
#        picklefile=open(pklpathname,'rb')
#        tuple_readProcSpec_result=cPickle.load(picklefile)
#        picklefile.close()
#    else:
#        # if not exists, load original, create pkl
#        # print 'slow load from original'
        # read ProcSpec file
    tuple_readProcSpec_result=_readProcSpec(pathname,getAcqTime=True)
#        # pickle this tuple containing xdata ydata and acquisitiontime
#        picklefile=open(pklpathname,'wb')
#        cPickle.dump(tuple_readProcSpec_result,picklefile,2) # pickle using most efficient protocol 2
#        picklefile.close()
    # convert to a mwave
    (xdata,ydata,acqmillitime)=tuple_readProcSpec_result
    specmwave=mwave(xdata,ydata)
    if xcoordinates!=None:
        specmwave=mwave(xcoordinates,specmwave(xcoordinates))
    specmwave.ProcSpec_time=acqmillitime
    return specmwave

#####################
#
# PLUGIN CODE
#
#####################

class ImportPluginExample(ImportPlugin):
    """An example plugin for reading a set of unformatted numbers
    from a file."""

    name = "ProcSpec"
    author = "Martinus Werts"
    description = "OceanOptics ProcSpec file reader"

    # Uncomment this line for the plugin to get its own tab
    promote_tab='ProcSpec'

    file_extensions = set(['.ProcSpec'])

    def __init__(self):
        ImportPlugin.__init__(self)
        self.fields = [
            ImportFieldText("name", descr="Dataset name", default="name"),
            ImportFieldCheck("bl_flag", 
                descr="Baseline correction (order 0)",default=False),
            ImportFieldFloat("bl_start", descr="Baseline start", default=800.),
            ImportFieldFloat("bl_end", descr="Baseline end", default=900.),
            ImportFieldCheck("slice_flag", 
                descr="Slice",default=False),
            ImportFieldFloat("slice_start", descr="Slice start", default=360.),
            ImportFieldFloat("slice_end", descr="Slice end", default=1010.)
#            ImportFieldCombo("subtract", items=("0", "1", "2"),
#                            editable=False, default="0")
            ]

    def doImport(self, params):
        """Actually import data
        params is a ImportPluginParams object.
        Return a list of ImportDataset1D, ImportDataset2D objects
        """
        spec_mwave=mwave_loadProcSpec(params.filename)       
#        f = params.openFileWithEncoding()
#        data = []
#        mult = params.field_results["mult"]
#        sub = float(params.field_results["subtract"])
        if params.field_results["bl_flag"]:
            bl_start=params.field_results["bl_start"]
            bl_end=params.field_results["bl_end"]
            spec_mwave.baseline_subtract([(bl_start,bl_end)])
        if params.field_results["slice_flag"]:
            slice_start=params.field_results["slice_start"]
            slice_end=params.field_results["slice_end"]
            (idx1,idx2)=spec_mwave.idxrange(slice_start,slice_end)
            slice_xdata=spec_mwave.x[idx1:idx2]
            slice_ydata=spec_mwave.y[idx1:idx2]
            spec_mwave=mwave(slice_xdata,slice_ydata)            
#        for line in f:
#            data += [float(x)*mult-sub for x in line.split()]
        return [ImportDataset1D(params.field_results["name"]+'_wl',
                spec_mwave.x),ImportDataset1D(params.field_results["name"]+'_val',
                spec_mwave.y)]

    def getPreview(self, params):
        """Generate preview
        In the case of ProcSpec, no preview is available
        """
        return ('No preview available.\n'+versionstring,True)
   
# add the class to the registry. An instance also works, but is deprecated
importpluginregistry.append(ImportPluginExample)

