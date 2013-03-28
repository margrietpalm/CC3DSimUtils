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

import numpy as np

def readColorMap(filename):
    """ Read colormap from a file, formatted like: celltype r g b
    
    :param filename: file with the colormap
    :type filename: str
    :return: dictionary with cell type as keys and colors (r,g,b) as values.
    """
    f = open(filename,'r')
    lines = f.readlines()
    colormap = {}
    for line in lines:
        try:
            sline = line.split()
            colormap[int(sline[0])]=(int(sline[1]),int(sline[2]),int(sline[3]))
        except:
            print 'bla'            
    f.close()
    return colormap

def _readGrid(simid,t,suffix,indir,gzipped,border):
    if gzipped:
        grid = np.loadtxt(indir+'/'+simid+'/'+simid+suffix+str(t)+'.data.gz')
    else:
        grid = np.loadtxt(indir+'/'+simid+suffix+str(t)+'.data')
    (nx,ny) = grid.shape
    if border:
        return _removeBorder(grid)
    else:
        return grid    

def _removeBorder(grid):
    # cut off border cell
    (nx,ny) = grid.shape
    grid[0:nx,0] = 0
    grid[0:nx,ny-1] = 0
    grid[0,0:ny] = 0
    grid[nx-1,0:ny] = 0  
    return grid 

def readSigma(simid,t,indir,gzipped=True,border=True):
    """ Read cell field (sigma) from file.
    
    :param simid: simulation identifier
    :type simid: str
    :param t: time step
    :type t: int
    :param indir: path to data files
    :type indir: str
    :param gzipped: data is gzipped (gzipped data is expected to be in indir/simid/)
    :type gzipped: bool
    :param border: cut of border pixels
    :type border: bool
    :return: numpy array with cell id's
    """
    return _readGrid(simid,t,'_CF_',indir,gzipped,border)
    
def readTau(simid,t,indir,gzipped=True,border=True):
    """ Read type field (tau) from file.
    
    :param simid: simulation identifier
    :type simid: str
    :param t: time step
    :type t: int
    :param indir: path to data files
    :type indir: str
    :param gzipped: data is gzipped (gzipped data is expected to be in indir/simid/)
    :type gzipped: bool
    :param border: cut of border pixels
    :type border: bool
    :return: numpy array with cell types
    """
    return _readGrid(simid,t,'_TF_',indir,gzipped,border)

def readChemField(simid,t,indir,fieldname,gzipped=True,border=True):
    """ Read chemical field from file.
    
    :param simid: simulation identifier
    :type simid: str
    :param t: time step
    :type t: int
    :param indir: path to data files
    :type indir: str
    :param fieldname: name of the chemical field
    :type fieldname: str
    :param gzipped: data is gzipped (gzipped data is expected to be in indir/simid/)
    :type gzipped: bool
    :param border: cut of border pixels
    :type border: bool
    :return: numpy array with the levels of the chemical field at each position
    """
    return _readGrid(simid,t,'_'+fieldname+'_',indir,gzipped,border)    
    