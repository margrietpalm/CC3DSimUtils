#-------------------------------------------------------------------#
# Copyright 2012-2013 Margriet Palm                                 #
#                                                                   #
# CC3DSimUtils is free software; you can redistribute               #
# it and/or modify it under the terms of the GNU General Public     #    
# License as published by the Free Software Foundation; either      #
# version 3 of the License, or (at your option) any later version.  #
#                                                                   #
# CC3DSimUtils is distributed in the hope that it will be useful,   #
# but WITHOUT ANY WARRANTY; without even the implied warranty of    #
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU  #
# General Public License for more details.                          #
#                                                                   #
# You should have received a copy of the GNU General Public License #
# along with CC3DSimUtils; if not, write to the Free Software       #
# Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA         #
# 02110-1301 USA                                                    #
#                                                                   #
#-------------------------------------------------------------------#

import copy
import numpy as np
from mahotas import labeled
from mahotas import _convex,polygon
import pymorph as m
from Readers import *

norm = lambda v: np.array(v)/np.linalg.norm(np.array(v))
angle = lambda x,y: np.arccos(np.dot(x,y)/(np.linalg.norm(x)*np.linalg.norm(y)))

def getCoM(pix):
    """ Calculate center of mass of a cell 
    
    :param pix: cell coordinates ([x1,...,xn],[y1,...,yn])
    :return: center of mass (x,y)
    """
    try:
        a = len(pix[0])
        cx = 1/float(a)*(np.sum(pix[0]))
        cy = 1/float(a)*(np.sum(pix[1]))
    except:
        print pix
    return (cx,cy)

def getCellInertiaTensor(pix):
    """ Get inertia tensor for a cell 
    
    :param pix: cell coordinates ([x1,...,xn],[y1,...,yn])
    :return: inertia tensor [[Ixx,Ixy],[Ixy,Iyy]]
    """
    a = len(pix[0])
    cx = 1/float(a)*(np.sum(pix[0]))
    cy = 1/float(a)*(np.sum(pix[1]))
    Ixx = np.sum(np.power(pix[0]-cx,2))
    Iyy = np.sum(np.power(pix[1]-cy,2))
    Ixy = -np.sum((pix[0]-cx)*(pix[1]-cy))
    return np.array([[Ixx,Ixy],[Ixy,Iyy]])

def getCellOrientation(pix):
    """ Calculate orientation of a cell. The orientation is the eigenvector corresponding to the largest eigenvalue of the cells' inertia tensor.
    
    :param pix: cell coordinates ([x1,...,xn],[y1,...,yn])
    :return: unit vector of the cell orientation

    .. seealso:: :func:`~AnalysisUtils.getCellInertiaTensor`        
    """    
    [V,D] = np.linalg.eig(getCellInertiaTensor(pix))
    return [D[V==max(V)][0][0],D[V==max(V)][0][1]]

def getCellAngle(pix):
    """ Calculate angle of a cell on interval :math:`[0,\pi]`
    
    :param pix: cell coordinates ([x1,...,xn],[y1,...,yn])
    :return: angle in radians on interval :math:`[0,\pi]`

    .. seealso:: :func:`~AnalysisUtils.getCellOrientation`            
    """
    l = getCellOrientation(pix)
    a = np.arctan2(l[1],l[0])
    if a < 0:
        a += np.pi
    return a

def getAngleField(sigma):
    """ Get field with the cell angles 
    
    :param sigma: numpy array with cell id's
    :return: numpy array with cell angles in radians
    """
    angles = np.zeros(sigma.shape)    
    for id in np.unique(sigma):
        if id == 0:
            continue
        angles[sigma==id] = getCellAngle(np.where(sigma==id))
    return angles

def getDoubledOrientationField(sigma):
    """ Get field with the doubled cell orientations. By doubling the orientations we map the orientations on the domain :math:`[0,2\pi)`.
    
    :param sigma: CPM grid
    :return: field with doubled orientations
    """    
    orients = np.zeros((sigma.shape[0],sigma.shape[1],2))    
    for id in np.unique(sigma):
        if id == 0:
            continue
        vec = np.array(getCellOrientation(np.where(sigma==id)))
        # remap orientation vectors on 0-pi
        if vec[1] < 0:
            vec = -1*vec
        orients[sigma==id] = [vec[0]*vec[1]*2,vec[0]*vec[0]-vec[1]*vec[1]]
    return orients
    
def getOrientationField(sigma):
    """ Get field with the cell orientations. 
    
    :param sigma: CPM grid
    :return: field with orientations
    """      
    orients = np.zeros((sigma.shape[0],sigma.shape[1],2))    
    for id in np.unique(sigma):
        if id == 0:
            continue
        vec = np.array(getCellOrientation(np.where(sigma==id)))
        # remap orientation vectors on 0-pi
        if vec[1] < 0:
            vec = -1*vec
        orients[sigma==id] = vec
    return orients

def getDirector(pos,r,sigma,orientations):
    """ Calculate the director at a given position. The director is the mean direction of all cells within a radius r. Because we represent the direction of cells as vectors, the mean direction is the resulting vector.  
        
    :param com: position (x,y) for which the order parameter is calculated
    :param r: radius of neighborhood
    :param sigma: CPM grid
    :param orientations: array with doubled cell orientations (as returned by :func:`~CC3DSimUtils.AnalysisUtils.getDoubledOrientationField`)
    :return: vector representing the director director
    """    
    # cut out part of grid we are interested in
    xmin = 0 if (com[0]-r) < 0 else int(np.floor(com[0]-r))
    xmax = sigma.shape[0] if (com[0]+r+1) > sigma.shape[0] else int(np.ceil(com[0]+r+1))
    ymin = 0 if (com[1]-r) < 0 else int(np.floor(com[1]-r))
    ymax = sigma.shape[1] if (com[1]+r+1) > sigma.shape[1] else int(np.ceil(com[1]+r+1))
    sigma = sigma[xmin:xmax,ymin:ymax]
    orientations = orientations[xmin:xmax,ymin:ymax,:]
    (x,y) = np.mgrid[xmin:xmax,ymin:ymax]
    # calculate distances on the grid to find pixels within the circle
    d = np.sqrt(np.power(x-com[0],2)+np.power(y-com[1],2))
    # select angles within r
    sa = sigma[d<=r].astype(np.int)
    cells = np.unique(sa[sa>0])
    o = np.array([orientations[np.where(sigma==id)][0] for id in cells])
    bins = np.bincount(sa[sa>0])
    w = bins[bins>0]
    sx = np.sum(w*o[:,0])
    sy = np.sum(w*o[:,1])    
    if np.linalg.norm([sx,sy]) > 0:
        return [sx,sy]/np.linalg.norm([sx,sy])
    else:
        return [sx,sy]

def getOrderParameter(sigma,orientations,r):
    """ Calculate order parameter for a morphology. 
    
    :param sigma: CPM grid
    :param orientations: array with doubled cell orientations (as returned by :func:`~CC3DSimUtils.AnalysisUtils.getDoubledOrientationField`)    
    :param r: radius of neighborhood
    :return: order parameter    
    
    .. seealso:: :func:`~CC3DSimUtils.AnalysisUtils.getLocalOrderParameter2D`, :func:`~CC3DSimUtils.AnalysisUtils.getGlobalOrderParameter2D`
    """       
    (nx,ny) = sigma.shape
    rmax = np.sqrt(nx**2+ny**2)
    if r < rmax:
        return getLocalOrderParameter2D(sigma,orientations,r)
    else:
        return getGlobalOrderParameter2D(sigma,orientations)    

def getGlobalOrderParameter(sigma,orientations):
    """ Calculate order parameter for a morphology using a radius larger than the diagonal of the CPM grid.
    
    :param sigma: CPM grid
    :param orientations: array with doubled cell orientations (as returned by :func:`~CC3DSimUtils.AnalysisUtils.getDoubledOrientationField`)    
    :return: order parameter    
    
    .. seealso:: :func:`~CC3DSimUtils.AnalysisUtils.getOrderParameter2D`
    """     
    cells = np.unique(sigma[sigma>0])
    bins = np.bincount(sigma[sigma>0].astype(np.int))
    w = bins[bins>0]  
    o = np.array([orientations[np.where(sigma==id)][0] for id in cells])
    sx = np.sum(w*o[:,0])
    sy = np.sum(w*o[:,1])    
    A = [sx,sy]/np.linalg.norm([sx,sy])
    s = 0.
    n = 0.
    for id in cells:
        a = orientations[np.where(sigma==id)][0]
        dA = getAngle(a,A)
        s += np.cos(dA)
        n += 1
    return s/n
    
def getLocalOrderParameter(sigma,orientations,r):
    """ Calculate order parameter for a morphology. 
    
    :param sigma: CPM grid
    :param orientations: array with doubled cell orientations (as returned by :func:`~CC3DSimUtils.AnalysisUtils.getDoubledOrientationField`)    
    :param r: radius of neighborhood
    :return: order parameter    
    
    .. seealso:: :func:`~CC3DSimUtils.AnalysisUtils.getOrderParameter2D`, :func:`~CC3DSimUtils.AnalysisUtils.getDirector2D`

    """          
    s = 0
    n = 0
    for id in np.unique(sigma):
        if id == 0:
            continue
        a = orientations[np.where(sigma==id)][0]
        A = getDirector2D(getCoM(np.where(sigma==id)),r,sigma,orientations)        
        dA = getAngle(a,A)
        s += np.cos(dA)
        n += 1
    return s/n

def getRelativeDirField(sigma,r):
    """ Calculate field with relative director for each pixel. The relative director is the difference to the angle of the cell at that pixel and the relative director on the pixel. Pixels with high values represent unordered areas, such as branchpoints.
    
    :param sigma: CPM grid
    :param r: radius used for calculating the director
    :return: field with relative director values
    
    .. seealso:: :func:`~CC3DSimUtils.AnalysisUtils.getDoubledOrientationField`
    """
    orientations = getDoubledOrientationField(sigma)  
    da = np.zeros_like(sigma)
    # calculate director at every pixel
    (nx,ny) = sigma.shape    
    for i in range(nx):
        for j in range(ny):
            if sigma[i,j] > 0:
                dir = getDirector2D((i,j),r,sigma,orientations)
                da[i,j] = getAngle(dir,orientations[i,j])/2.
    return da    

def getCellAngle(pix):
    """ Calculate angle of a cell on interval :math:`[0,\pi]`
    
    :param pix: cell coordinates ([x1,...,xn],[y1,...,yn])
    :return: angle in radians on interval :math:`[0,\pi]`
    """
    l = np.array(getCellOrientation(pix))
    a = np.arctan2(l[1],l[0])
    return a%np.pi
    
def realconvexhull(bwimg):
    """ 
    Calcuate the real convex hull of an image. The default convex hull algorithms, part of mahotas, returns the coordinates of the pixels on the hull. Therefore, that convex hull intersects with the outer pixels. Here we calculate the hull over a set containing the 4 border points of each pixel.
    
    :param bwimg: black-and-white image with 1 for occupied pixels and 0 at empty pixels
    :type bwimg: numpy array
    :return: numpy array of points in the hull    
    """
    # based on convexhull from mahotas.polygon
    # extra points are added to prevent intersection between the hull and pixel centers
    Y,X = np.where(bwimg)
    newim = np.zeros((2*bwimg.shape[0],2*bwimg.shape[1]))
    for i in range(len(Y)):
        y = Y[i]
        x = X[i]
        newim[2*x-1:2*x+2,2*y-1:2*y+2] = 1
    return np.array(polygon.convexhull(newim))/2.0

def getCompactness(sigma,minval=0):
    """ 
    Calculate compactness of a morphology: :math:`\\frac{A_{\\text{area}}}{A_{\\text{convex hull}}}` .

    :param sigma: numpy array with cell id's
    :param minval: minimum cell id for non medium pixels
    :type minval: int
    :return: compactness
    """
    # 1: create negative image
    (nx,ny) = sigma.shape    
    imneg = np.uint8(sigma>minval)
    # 2: find borders of negative image, these are all pixels where a neighbor is not cell
    #    thus, all borders are 2 pixels thick! 
    bim = labeled.borders(np.asarray(imneg))    
    # 3: get convex hull with intersection between bim and imneg -> the boundary pixels of the cells
    #    ph is an array of points on the hull
    # Note that the input are only the pixels on the border of the image, thus speeds up the 
    # calculations A LOT.
    ph = realconvexhull(imneg&bim)
    # 4: calculate area of hull according to http://en.wikipedia.org/wiki/Polygon#Area_and_centroid
    #    note that first point on the hull is added to the end for this calculation
    n = len(ph)
    ph = np.vstack((ph,ph[0,:]))
    A = -0.5*np.sum(ph[0:n,0]*ph[1:n+1,1]-ph[1:n+1,0]*ph[0:n,1])
    # return compactness
    return np.sum(imneg)/float(A)

def getLCC(sigma):
    """ Find largest connected component of an image 
    
    :param sigma: numpy array with cell id's
    :return: image with largest connected component
    """
    lab,nr = labeled.label(sigma)
    sizes = np.bincount(np.reshape(lab,(lab.shape[0]*lab.shape[1])))
    return lab==np.where(sizes==max(sizes[1:len(sizes)]))[0][0]
    
def getCellClusters(field,sigma,th=15,minlabsize=50,opendisk=1,mincellsize=.25):    
    """ Get clusters for a single morphology. 
    
    :param field: numpy array with values on which data is seperated
    :param cells: dict with cell identifiers as keys and :class:`~Cell.ClusterCell` instances as values
    :param sigma: CPM grid
    :param th: threshold value for step 1
    :type th: number
    :param minlabsize: labelled areas smaller than this value are ignored (2b)
    :type minlabsize: int
    :param opendisk: disk size for opening operation (2a)
    :type opendisk: int
    :param mincellsize: minimal fraction of the cell that must be on the labelled area to be added to the cluster
    :type mincellsize: number
    :return: dictionary with cluster id as key and :class:`~AnalysisUtils.Cluster` instances
    
    .. seealso:: :class:`~AnalysisUtils.Cluster` 
    """
    cfield = np.zeros_like(field)
    cfield[field<th] = 120
    cfield[sigma==0] = 0
    celldict = dict((cellid,[]) for cellid in np.unique(sigma[sigma>0]))
    labnum,lab = _getBlobs(cfield,minsize=minlabsize,opendisk=opendisk)
    clusters = {}
    for n in labnum:
        c = Cluster(n)
        if n == 0:
            continue
        cellpos = list(sigma[np.where(lab==n)])
        for cellid in np.unique(sigma[sigma>0]):
            if cellpos.count(cellid) >= mincellsize*np.sum(sigma==cellid):                
                celldict[cellid].append(n)
                c.addCell(cellid)
        if c.getClusterSize() > 0:
            clusters[n] = c
    nolab = 0
    for cellid,labels in celldict.iteritems():
        if len(labels) == 0:
            nolab += 1           
        elif len(labels) > 1:
            lc = [clusters[lab].getClusterSize() for lab in labels]
            cmax = labels[lc.index(max(lc))]
            for lab in labels:
                if lab is not cmax:
                    empty = clusters[lab].removeCell(cellid)
                    if empty:
                        del clusters[lab]
    return clusters

def _getBlobs(im,minsize=10,opendisk=0):
    """ Find blobs in an image. 
    
    :param im: input image
    :type im: numpy array
    :param minsize: minimal size of detected blobs
    :type minsize: int
    :param opendisk: disk size for opening operation (ommitted when opendisk <= 0)
    :type opendisk: int
    :return: identifiers for areas with blobs > minsize and labeled image (lac,lab)
    """
    if opendisk > 0:
        im = m.open(im.astype(np.int),m.sedisk(opendisk))  
    n = np.array([[0,1,0],[1,1,1],[0,1,0]])
    lab,nr = labeled.label(im,n)
    sizes = np.bincount(np.reshape(lab,(lab.shape[0]*lab.shape[1])))
    lac = np.nonzero(sizes>minsize)[0]
    return lac,lab

def calcMSDTransForCellTC(com):
    """ Calculate the translational MSD for a single object.
    
    :param com: list of centers of mass at each time step
    :return: list with MSD for each time step
    """
    return np.power(com-com[0,:],2).sum(axis=1)

def calcMSDRotForCellTC(vecset):
    """ Calculate the rotational MSD for a single object.
    
    :param vecset: list of orientation vectors for each time step
    :return: list with MSD for each time step
    """
    
    # 1) Calculate angles between two vectors of consecutive time steps: $a$, this is used to find the best fitting orientation of each cell:
    a = np.array([angle(vecset[i],vecset[i-1]) for i in range(1,len(vecset))])
    # Find all indices in $a$ where the $a > \frac{1}{2}\pi$; here the opposite orientation should be used
    # Rotate the vector corresponding to that angle 180 degrees.
    # Update angles that are calculated with that vector.
    # Repeat these steps untill all angles are correct.
    # Now we have a set of vectors which in which all cells have the most likely orientation: vecset
    while len(np.where(a>.5*np.pi)[0]):
        i = np.where(a>.5*np.pi)[0][0]
        vecset[i+1] = -1*vecset[i+1]
        a[i] = angle(vecset[i+1],vecset[i])
        if i+2 < len(vecset):
            a[i+1] = angle(vecset[i+2],vecset[i+1])
    
    # 2) Calculate the angle between the possitive x-axis and each vector
    xa = np.arctan2(vecset[:,1],vecset[:,0])
    # Remap these angles from the [0,-pi] domain to the [pi,2pi] domain.
    xa[np.where(xa<0)] = 2*np.pi+xa[np.where(xa<0)]
    
    # 3) Calculate differences between consecutive angles
    # calculate differences between angles
    da = xa[1:len(xa)] - xa[0:-1]
    # detect transition over 2*np.pi and correct the angular differences
    da[np.where(da<-1*np.pi)] = 2*np.pi+da[np.where(da<-1*np.pi)]
    da[np.where(da>1*np.pi)] = -2*np.pi+da[np.where(da>1*np.pi)]
    # summarize all angular differences to get a(t)-a(0) = sum(a(t)-a(t-1)) for eacht t
    
    # 4) Calculate cumulative sum of the angles and square
    sda = np.cumsum(da)   
    return np.power(np.concatenate(([0],sda)),2)

class Cluster:
    """ Container for a cell cluster
    
    :param id: cluster id
    :ivar cells: list of ids of the cells in the clusters
    """
    def __init__(self,id):
        self.id = id
        self.cells = []
        self.cellinterface = -1
        self.ecminterface = -1
        self.perimeter = -1
    
    def addCell(self,cellid):
        """ Add cell to cluster 
        
        :param cellid: id of cell
        """
        if cellid not in self.cells:
            self.cells.append(cellid)
    
    def removeCell(self,cellid):
        """ Remove cell from cluster
        
        :param cellid: id of cell
        """
        if cellid in self.cells:
            self.cells.remove(cellid)
        if len(self.cells) == 0:
            return True
        else:
            return False
    
    def getClusterSize(self):
        """ Calculate number of cell in cluster
        
        :return: number of cells in cluster
        """
        return len(self.cells)

class ClusterCellTC():
    """ A class that holds properties related to a cell at each measured time step. These properties are:

        * cluster id and size at each time step
        * long axis at each time step
        * center of mass at each time step
    
    :param id: cell id
    :ivar id: cell id
    :ivar clusterId: list of cluster id's
    :ivar clusterSize: list of cluster sizes
    :ivar time: list of time steps
    :ivar laxis: 2D array with long axes of the cell
    :ivar com: 2D array with centers of mass of the cell
    """    
    def __init__(self,id):
        self.id = id
        self.clusterId = np.array([],dtype=np.int)
        self.clusterSize = np.array([],dtype=np.int)
        self.time = np.array([],dtype=np.int)
        self.laxis = np.empty((0,2),dtype=np.float)  
        self.com = np.empty((0,2),dtype=np.float)
        
    def addTimeStep(self,t,pix,cid,csz):   
        """ Add time step
        
        :param t: time step
        :type t: int
        :param pix: cell coordinates ([x1,...,xn],[y1,...,yn])        
        :param cid: cluster id
        :type cid: int
        :param csz: cluster size
        :type csz: int
        """
        self.time = np.append(self.time,(int(t)))
        self.com = np.vstack((self.com,getCoM(pix)))
        self.laxis = np.vstack((self.laxis,getCellOrientation(pix)))
        self.clusterId = np.append(self.clusterId,(cid))
        self.clusterSize = np.append(self.clusterSize,(csz))  