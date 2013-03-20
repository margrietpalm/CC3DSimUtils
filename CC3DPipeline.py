import time,string
import numpy as np
from .Readers import *
from .AnalysisUtils import *
from .ImageUtils import *

#----- Pre-processing ----#
def createPBSScripts(runid,joblist,command,time,ncores=8,ppn=8,path='batchScripts/'):
    """ Create a set of PBS scripts to run a simulation on a cluster. Each script starts with something like:
    
        #PBS -S /bin/bash
        #PBS -lnodes=1:cores12:ppn=11
        #PBS -lwalltime=12:00:00
    
    If these commands are not correct or complete for the cluster you use, edit :func:`~CC3DPipeline.createPBS`. 
    
    For each job in joblist a single line command is added to the script:
        
        python command jobid > log/jobid.out 2> log/jobid.err &
            
    :param runid: identifier for the scripts
    :type runid: str
    :param joblist: list of job identifiers
    :param command: command that runs the simulation
    :type command: str
    :param time: requested walltime on the cluster (hh:mm:ss)
    :param ncores: numbor of cores in the requested node
    :type ncores: int
    :param ppn: number of processers per node that will be used
    :type ppn: int
    :param path: location where pbs scripts are saved
    :type path: str
    
    .. seealso:: :func:`~CC3DPipeline.createPBS`, :func:`~CC3DPipeline.addCommandToPBS`, :func:`~CC3DPipeline.finishPBS`        
    """
    nscripts = 0
    if ppn > ncores:
        ppn = ncores
    for i,jobid in enumerate(joblist):
        if (i%ppn == 0):
            if nscripts > 0:
                finishPBS(fn)
            fn = path+'/run_'+runid+'_part_'+str(nscripts)
            createPBS(fn,time,ncores,ppn,xserver,python27)
            nscripts += 1
        addCommandToPBS(fn,command+' '+jobid,'log/'+jobid)
    finishPBS(fn)

def createPBS(filename,time,ncores=None,ppn=None):
    """ Create a new pbs script and add initial commands and settings.
    
    :param filename: filename of the new pbs script
    :type filename: str
    :param time: requested walltime on the cluster (hh:mm:ss)
    :param ncores: numbor of cores in the requested node
    :type ncores: int
    :param ppn: number of processers per node that will be used
    :type ppn: int
    """
    f = open(filename,'w')
    f.write('#PBS -S /bin/bash\n')
    f.write('#PBS -lnodes=1')
    if (ncores == 8) or (ncores == 12) or (ncores == 16):
        f.write(':cores'+str(ncores))
    if ppn:
        f.write(':ppn='+str(ppn))
    f.write('\n')
    f.write('#PBS -lwalltime='+str(time)+'\n')
    f.write('cd $HOME\n')
    f.close()

def addCommandToPBS(filename,command,log):
    """ Add single line command to existing PBS script:
    
    :param filename: filename of the new pbs script
    :type filename: str    
    :param command: command that runs the simulation
    :type command: str   
    :param log: name (with path) of the log files (without extension)
    :type log: str
    """
    f = open(filename,'a')
    f.write(command+' > '+log+'.out 2> '+log+'.err &\n')
    f.close()

def finishPBS(filename):
    """ Finish pbs file 
    
    :param filename: filename of the new pbs script
    :type filename: str      
    """
    
    f = open(filename,'a')
    f.write('wait\n')
    f.close()


#----- Post-processing ----#
def makeImages(id,trange,inpath,outpath,cm='default.ctb',gzipped=False,timestamp=False,label=False,scale=1,bc=None,fontsize=6,fieldname=None,border=True):
    """ Make images for a single simulation simulation
    
    :param id: simulation identifier
    :type id: str
    :param trange: list of time steps for which images are created
    :param inpath: path to data 
    :type inpath: str 
    :param outpath: path to save images to
    :type outpath: str
    :param cm: file containing the colormap
    :type cm: str
    :param gzipped: data is gzipped
    :type gzipped: bool
    :param timestamp: add time stamp to the image
    :type timestamp: bool
    :param label: add id as label to the image
    :type label: bool
    :param scale: scaling of the image
    :type scale: number
    :param bc: color of cell boundaries (r,g,b)
    :param fontsize: size of the fonts used for label and time stamp; font size will be multiplied by scale.
    :type fontsize: int
    :param fieldname: name of chemical field
    :type fieldname: str
    :param border: cut of border pixels
    :type border: bool
    
    .. seealso:: :func:`~ImageUtils.makeImage`
    """
    t0 = time.time()
    tlen = len(str(trange[-1]))
    colormap = readColorMap(cm)
    for t in trange:
        outfile = id+'_'+string.zfill(str(t),tlen)+'.png'    
        im = makeImage(id,inpath,t,colormap,timestamp,label,scale,bc,fontsize,border,gzipped,fieldname)
        im.save(outpath+outfile)
    print str(len(trange)) + ' images drawn in '+str(time.time()-t0)+' seconds'

def getOrderParameterForSim(id,trange,inpath,radii,gzipped=False,border=True,outpath=None):
    """ Calculate orderparameters for one simulation. All order parameters are collected and save in a file inpath/id_orderparameter.data
    
    :param id: simulation identifier
    :type id: str
    :param trange: list of time steps for which the order parameter is calculated
    :param inpath: path to data 
    :type inpath: str
    :param radii: list of radii for wich the order parameter is calculates
    :param gzipped: if True, data is expected to be gzipped, and stored in inpath/id/, and the output file will be gzipped and stored in outpath/id/
    :type gzipped: bool
    :param border: remove border pixels from data        
    :type border: bool
    :param outpath: path where order parameter data will be saved, if omitted outpath = inpath
    :type outpath: str
    
    .. seealso:: :func:`~AnalysisUtils.getOrderParameterFromGrid`
    """    
    if outpath is None:
        outpath = inpath
    t0 = time.time()
    if gzipped:
        f = gzip.open(outpath+'/'+id+'/'+id+'_orderparameter.data.gz','w')
    else:
        f = open(outpath+'/'+id+'_orderparameter.data','w')
    data = np.zeros((len(trange),len(radii)+1))
    f.write('#MCS')
    for r in radii:
        f.write('\t'+str(r))
    f.write('\n')
    for i,t in enumerate(trange):
        sigma = readSigma(id,t,inpath,gzipped,border)
        (nx,ny) = sigma.shape
        data[i,0] = t
        angles = getAngleField(sigma)
        for j,r in enumerate(radii):
            data[i,j+1] = getOrderParameterFromGrid(sigma,angles,r)          
    np.savetxt(f,data)
    f.close()
    print 'Order parameter calculated for '+str(len(trange))+' simulations and '+str(len(radii))+' radii in '+str(time.time()-t0)+' seconds'

def getCompactnessForSim(id,trange,inpath,gzipped=False,border=True,outpath=None):
    """ Calculate compactness for one simulation
    
    :param id: simulation identifier
    :type id: stre
    :param trange: list of time steps for which the compactness is calculated
    :param inpath: path to data 
    :type inpath: str 
    :param gzipped: if True, data is expected to be gzipped, and stored in inpath/id/, and the output file will be gzipped and stored in outpath/id/
    :type gzipped: bool
    :param border: remove border pixels from data        
    :type border: bool
    :param outpath: path where order parameter data will be saved, if omitted outpath = inpath
    :type outpath: str
    
    .. seealso:: :func:`~AnalysisUtils.getCompactness`
    """
    t0 = time.time()
    if outpath is None:
        outpath = inpath
    if gzipped:
        f = gzip.open(outpath+'/'+id+'/'+id+'_compactness.data.gz','w')
    else:
        f = open(outpath+'/'+id+'_compactness.data','w')    
    data = np.zeros((len(trange),2))            
    for i,t in enumerate(trange):            
        sigma = readSigma(id,t,inpath,gzipped,border)          
        data[i,:] = [t,getCompactness(getLCC(sigma))]
    f.write('#MCS\tcompactness\n')
    np.savetxt(f,data)
    f.close()
    print 'compactness calculated for '+str(len(trange))+' morphologies in '+str(time.time()-t0)+' seconds'

def getClustersForSim(id,trange,inpath,r,th,minlabsize,opendisk,mincellsize,gzipped=False,border=False,outpath=None):
    """ Calculate clusters and mean squared displacement and rotation for each cell in a simulation. For more details on clustering see the documentation of :func:`~AnalysisUtils.getCellClusters`.

    :param id: simulation identifier
    :type id: str
    :param trange: list of time steps for which the clusters are calculated
    :param inpath: path to data
    :type inpath: str
    :param r: radius for relative director field
    :type r: number
    :param th: threshold value for step 1
    :type th: number
    :param minlabsize: labelled areas smaller than this value are ignored (2b)
    :type minlabelsize: int
    :param opendisk: disk size for opening operation (2a)
    :type opendisk: int
    :param mincellsize: minimal fraction of the cell that must be on the labelled area to be added to the  cluster    
    :type mincellsize: int
    :param gzipped: if True, data is expected to be gzipped, and stored in inpath/id/, and the output file will be gzipped and stored in outpath/id/
    :type gzipped: bool
    :param border: remove border pixels from data        
    :type border: bool
    :param outpath: path where order parameter data will be saved, if omitted outpath = inpath
    :type outpath: str   
    
    .. seealso:: :func:`~AnalysisUtils.getCellClusters`, :func:`~AnalysisUtils.getRelativeDirField`, :func:`~AnalysisUtils.calcMSDTransForCellTC`, :func:`~AnalysisUtils.calcMSDRotForCellTC`, :class:`~AnalysisUtils.ClusterCellTC`
    """    
    if outpath is None:
        outpath = inpath
    cellTCdict = {}
    t0 = time.time()
    for t in trange:
        sigma = readSigma(id,t,inpath,gzipped=gzipped,border=border)
        diffield = getRelativeDirField(sigma,r)
        clusters = getCellClusters(diffield,sigma,th,minlabsize,opendisk,mincellsize)            
        cellids = np.unique(sigma)
        # remove ECM (sigma=0) from cellids
        for cid,cluster in clusters.iteritems():
            csz = cluster.getClusterSize()
            # save data per cells
            for cellid in cluster.cells:
                p = np.where(sigma==cellid)
                cellTCdict.setdefault(cellid,ClusterCellTC(cellid)).addTimeStep(t,p,cid,csz) 
                cellids = np.delete(cellids,np.where(cellids==cellid))                
        # save data for cells not in a cluster
        for cellid in cellids:
            if cellid == 0:
                continue
            p = np.where(sigma==cellid)
            cellTCdict.setdefault(cellid,ClusterCellTC(cellid)).addTimeStep(t,p,np.nan,np.nan)
    # write data to file per cell
    for cellid,celltc in cellTCdict.iteritems():
        data = np.column_stack((celltc.time,calcMSDTransForCellTC(celltc.com),
                            calcMSDRotForCellTC(celltc.laxis),celltc.clusterId,celltc.clusterSize))
        f = open(outpath+'/'+id+'_cell_'+str(int(cellid))+'_cluster+MSD.data','w')
        f.write('#MCS\tMSDTrans\tMSDRot\tclusterId\tclusterSize\n')
        np.savetxt(f,data)
        f.close() 
    print 'clusters calculated for '+str(len(trange))+' morphologies in '+str(time.time()-t0)+' seconds'

