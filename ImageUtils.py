#*******************************************************************#
#                                                                   #
#-------------------------------------------------------------------#
#                                                                   #
#*******************************************************************#
import copy,os
import numpy as np
from PIL import Image,ImageDraw,ImageFont
from mahotas import labeled
from .Readers import *

def _getFont(fontpath,fontsize):
    try:
        if os.path.isfile(fontpath+"/times.ttf"):
            return ImageFont.truetype(fontpath+"/times.ttf", fontsize)        
        else:
            #~ return ImageFont.load("ariel.pil",
            return ImageFont.truetype("times.ttf", fontsize) 
    except:
        print 'I cannot use a nice font (freetype support may be missing)'
        return ImageFont.load_default()
        
def addTimeStamp(im,stamp,fs=6,fc=None,fontpath='/usr/share/fonts/msttcore/'):
    """ Draw a timestamp at the right bottom of the image. """
    if fc is None:
        fc = (0,0,0)
    draw = ImageDraw.Draw(im)
    fontsize = int(fs)
    font = _getFont(fontpath,fontsize)
    (w,h) = draw.textsize(str(stamp),font=font)
    (nx,ny) = im.size
    y = ny-h
    x = nx-w
    draw.text((x,y),str(stamp),fill=fc,font=font)

def addLabel(im,label,fs=8,fc=None,fontpath='/usr/share/fonts/msttcore/'):
    """ Draw label at the top center of the image. """
    if fc is None:
        fc = (0,0,0)    
    draw = ImageDraw.Draw(im)
    fontsize = int(fs)
    font = _getFont(fontpath,fontsize)    
    (w,h) = draw.textsize(str(label),font=font)
    (nx,ny) = im.size
    y = h
    x = int((0.5*nx)-w)
    #~ print 'add label ',x,y
    draw.text((x,y),str(label),fill=fc,font=font) 

def _getImAsArray(sigma,types,field,colormap,scale=1,nx=None,ny=None,transparant=False,gv=None,bc=None):
    if bc is None:
        bc = (0,0,0)
    (nx,ny) = sigma.shape 
    sigma = np.transpose((sigma.astype(np.float)))
    types = np.transpose((types.astype(np.float)))
    field = np.transpose((field.astype(np.float)))
    if np.sum(sigma) == 0:
        im = Image.fromarray(255*np.uint8(np.ones_like(sigma)))
        im = im.resize((int(scale*nx),int(scale*ny)))
        return np.asarray(im),sigma,np.zeros_like(sigma)
    if gv is None:
        gv = [100,255]
    if ((np.max(field)-np.min(field)) > 0):
        field = (((field-np.max(field))/(np.min(field)-np.max(field)))*(gv[1]-gv[0]))+gv[0]
    fim = Image.fromarray(field)
    fim = fim.resize((int(scale*nx),int(scale*ny)))
    fim.convert("RGB")
    im = Image.fromarray(sigma)
    im = im.resize((int(scale*nx),int(scale*ny)))
    sigma = np.asarray(im)
    tim = Image.fromarray(np.uint8(types))
    tim = tim.resize((int(scale*nx),int(scale*ny)))            
    types = np.asarray(tim)
    sigMin = 0
    bim = None
    while sigMin < np.max(sigma):
        imtemp = copy.deepcopy(np.asarray(im))
        imtemp[np.where(imtemp < sigMin)] = 0
        imtemp[np.where(imtemp > sigMin+254)] = 0
        if bim is None:
            bim = labeled.borders(imtemp)
        else:
            bim += labeled.borders(imtemp)
        sigMin += 254
    if np.sum(field) == 0:
        imnew = colormap[0]*np.ones((scale*ny,scale*nx,3))
    else:
        imnew = np.dstack((np.asarray(fim),np.asarray(fim),np.asarray(fim)))
    if not transparant:
        for tp in np.unique(types):
            if tp == 0:
                continue
            imnew[types==tp] = colormap[tp]
    #~ imnew[types==0] = colormap[0]
    # set border black
    imnew[bim] = bc
    return imnew,sigma,np.asarray(tim)

def makeImage(id,inpath,t,colormap,timestamp=False,label=False,scale=1,bc=None,fontsize=6,border=True,gzipped=True,fieldname=None):
    """ Draw morphology for one timestep 
    
    :param id: simulation identifier
    :param inpath: path containing data files
    :param t: time step
    :param outpath: path to save images to
    :param colormap: dictionary with cell types as keys and colors (r,g,b) as values
    :param timestamp: add time stamp to the image
    :type timestamp: bool
    :param label: add id as label to the image
    :type label: bool
    :param scale: scaling of the image
    :type scale: number    
    :param bc: color of cell boundaries (r,g,b)
    :param fontsize: size of the fonts used for label and time stamp; font size will be multiplied by scale.
    :type fontsize: int    
    :param border: cut of border pixels
    :type border: bool
    :param gzipped: data is gzipped
    :type gzipped: bool    
    :param fieldname: name of chemical field
    :type fieldname: str    
    :return: image object
    
    .. seealso:: :func:`~ImageUtils.drawCells`, :func:`~ImageUtils.addTimeStamp`, :func:`~ImageUtils.addLabel`    
    """          
    sigma = readSigma(id,t,inpath,gzipped,border)
    types = readTau(id,t,inpath,gzipped,border)
    if fieldname is not None:
        field = readChemField(id,t,inpath,fieldname,gzipped,border)
    else:
        field = np.zeros_like(sigma) 
    (nx,ny) = sigma.shape
    (imnew,sigma,types) = _getImAsArray(sigma,types,field,colormap,scale=scale,bc=bc) 
    im = Image.fromarray(np.uint8(imnew))
    im = im.rotate(90)
    if timestamp:
        fs = int(fontsize*(nx/200.0)*scale)
        addTimeStamp(im,str(t),fs,fc=bc)
    if label:
        fs = int(1.25*fontsize*(nx/200.0)*scale)
        addLabel(im,str(id),fs,fc=bc)
    return im

def drawRelDirField(field,sigma,scale=1):
    """ Draw gray-scale image of a field of the relative director
    
    :param field: numpy array with field
    :param scale: scaling factor
    :param gv: minimum and maximum color value
    :param raw: use raw values of field, otherwise values are mapped on the range of gv
    :return: image object
    """
    (nx,ny) = field.shape
    gv = [0,220]
    field = (((field-np.max(field))/(np.min(field)-np.max(field)))*(gv[1]-gv[0]))+gv[0]
    field[np.where(sigma==0)] = 255
    fim = Image.fromarray(field)
    fim = fim.resize((int(scale*nx),int(scale*ny)))
    fim.convert("RGB")
    im = Image.fromarray(np.uint8(np.dstack((np.asarray(fim),np.asarray(fim),np.asarray(fim)))))
    im = im.transpose(Image.FLIP_TOP_BOTTOM)
    #~ im = im.rotate(90)    
    return im
    
def stackImages(images,geometry,filename,label=False,title=None,fontsize=20,border=False,scale=1,fontpath="/usr/share/fonts/msttcore/"):
    """ Stack a set of images together in one image.
    
    :param images: dictionary with labels as keys and image filenames as values
    :param geometry: number of rows and columns (x,y)
    :param filename: target of the stacked image
    :param label: add labels to the subimages
    :param title: overall title for image
    :param fontsize: fontsize for labels and title
    :param border: add border to subimages
    :type border: bool
    :param scale: scaling factor of the created picture
    """
    if (label) or (title is not None):
        font = _getFont(fontpath,fontsize)            
    fontsize=int(fontsize*scale)    
    labels = images.keys()
    nx = 0
    ny = 0
    for im in images.values():
        sz = Image.open(im).size
        if sz[0] > nx:
            nx = sz[0]
            ny = sz[1]
    imsize = (scale*nx,scale*ny)
    if label:
        im = Image.new('RGB',(10,10))
        draw = ImageDraw.Draw(im)
        lsize = draw.textsize(str(labels[0]),font=font)
        imsize = (int(scale*imsize[0]),int(scale*imsize[1])+lsize[1])
        #~ imsize[1] += lsize[1]
    #---- Create new image ----#
    offsetX = 5
    offsetY = 0
    imw = int(np.ceil(imsize[0]*geometry[0]+2*offsetX))
    imh = int(np.ceil(imsize[1]*geometry[1]+offsetY))
    if title is not None:
        im = Image.new('RGB',(10,10))
        draw = ImageDraw.Draw(im)
        tsize = draw.textsize(str(title),font=font)
        im = Image.new('RGBA',(imw,imh+tsize[1]),(255,255,255))
        offsetY = tsize[1]
    else:
        im = Image.new('RGBA',(imw,imh),(255,255,255))
        offsetY = 0
    draw = ImageDraw.Draw(im)
    
    #----- Put plots on canvas ----#
    labels.sort()
    for i,l in enumerate(labels):
        x0 = int(i%geometry[0]*imsize[0]+offsetX)
        y0 = int((i/geometry[0])*imsize[1]+offsetY)
        newim = Image.open(images[l])
        im.paste(newim.resize((int(nx*scale),int(ny*scale))),(x0,y0))
        if border:
            draw.rectangle([(x0,y0),(x0+int(nx*scale),y0+int(ny*scale))],outline=(0,0,0))        
    if label:
        for i,l in enumerate(labels):
            tsize = draw.textsize(str(l),font=font)
            x0 = (i%geometry[0])*imsize[0]+0.5*imsize[0]-0.5*tsize[0]+offsetX
            y0 = (i/geometry[0])*imsize[1]+1*imsize[1]-tsize[1]+offsetY
            draw.text((x0,y0),str(l),fill=(0,0,0),font=font)
    if not (title is None):
        tsize = draw.textsize(str(title),font=font)        
        draw.text((im.size[0]/2.0-0.5*tsize[0]+offsetX,0),str(title),fill=(0,0,0),font=font)
    
    #----- Save image -----#
    im.save(filename)
    
def morphImages(images,filename,xlabel=None,ylabel=None,xtics=None,ytics=None,fontsize=20,scale=1,border=False,title=None,bcolor=(255,255,255),fcolor=(0,0,0),fontpath='/usr/share/fonts/msttcore/',delta=0):
    """ Stack a set of images together in one morphospace.
    
    :param images: 2D array with image filenames
    :param filename: target of the stacked image
    :param xlabel: label to be plotted on x-axis
    :param ylabel: label to be plotted on y-axis
    :param xtics: items on x-axis
    :param ytics: items on y-axis
    :param fontsize: fontsize for labels and title
    :param scale: scaling factor of the created picture
    :param border: add border to subimages    
    :param title: overall title for image    
    """
    if (xlabel is not None) or (ylabel is not None) or (ytics is not None) or (xtics is not None):
        font = _getFont(fontpath,fontsize)
        tfont = _getFont(fontpath,int(1.1*fontsize))
    orgsize = Image.open(images[0][0]).size
    subimsize = (int(scale*orgsize[0]+delta),int(scale*orgsize[1]+delta))
    imsize = [subimsize[0]*len(images[0])-delta,subimsize[1]*len(images)-delta]
    oX = 0
    oY = 0
    toX = 0
    toY = 0
    if xlabel is not None:
        im = Image.new('RGB',(10,10))
        draw = ImageDraw.Draw(im)
        lsize = draw.textsize(str(xlabel),font=font)
        imsize[1] = imsize[1]+lsize[1]
        oY += lsize[1]
        toY += lsize[1]
    if ylabel is not None:
        im = Image.new('RGB',(10,10))
        draw = ImageDraw.Draw(im)
        lsize = draw.textsize(str(ylabel),font=font)
        imsize[0] += lsize[1]
        oX += lsize[1]
        toX += lsize[1]
    if xtics is not None:
        im = Image.new('RGB',(10,10))
        draw = ImageDraw.Draw(im)
        lsize = draw.textsize(str(xtics[0]),font=font)
        imsize[1] += lsize[1]
        oY += lsize[1]
    if ytics is not None:
        im = Image.new('RGB',(10,10))
        draw = ImageDraw.Draw(im)
        lsize = draw.textsize(str(ytics[0]),font=font)
        imsize[0] += lsize[1]
        oX += lsize[1]
    loY = 0
    if title is not None:
        im = Image.new('RGB',(10,10))
        draw = ImageDraw.Draw(im)
        lsize = draw.textsize(str(title),tfont)
        offset = lsize[1]+20*scale        
        imsize[1] += offset
        oY += offset
        if xlabel is not None:
            loY += offset
            toY += offset
    #---- Create new image ----#    
    im = Image.new('RGBA',imsize,bcolor)
    draw = ImageDraw.Draw(im)        
    #----- Put plots on canvas ----#
    for i in range(len(images)):
        for j in range(len(images[0])):
            #~ print i,j,len(images),len(images[0])
            newim = Image.open(images[i][j])
            x0 = j*subimsize[0]+oX
            y0 = i*subimsize[0]+oY
            im.paste(newim.resize((int(scale*orgsize[0]),int(scale*orgsize[1]))),(x0,y0))
            #~ im.paste(newim.resize(subimsize),(x0,y0))
            if border:
                draw.rectangle([(x0,y0),(x0+newim.size[0],y0+newim.size[1])],outline=(0,0,0))        
    #----- Draw axis and tics -----#
    if (ytics is not None) or (xtics is not None):
        ticsfont = _getFont(fontpath,int(.75*fontsize))
    if xtics is not None:
        for i in range(len(xtics)):
            tsize = draw.textsize(str(xtics[i]),font=ticsfont)
            x0 = (i+0.5)*subimsize[0]+oX-int(tsize[0]/2.0)
            y0 = toY
            draw.text((x0,y0),str(xtics[i]),fill=fcolor,font=ticsfont)
    if ytics is not None:
        for i in range(len(ytics)):
            x0 = toX
            y0 = int((i+.5)*subimsize[1]+oY)
            #~ draw.text((x0,y0),str(ytics[i]),fill=(0,0,0),font=ticsfont)    
            tsize = draw.textsize(str(ytics[i]),font=ticsfont)
            tim = Image.new('RGBA',(tsize[0],tsize[1]),bcolor)
            tdraw = ImageDraw.Draw(tim)        
            tdraw.text((0,0),str(ytics[i]),fill=fcolor,font=ticsfont)
            tim = tim.rotate(90)
            y0 -= int(tsize[0]/2.0)
            im.paste(tim,(x0,y0))
    if xlabel is not None:
        tsize = draw.textsize(str(xlabel),font=font)        
        draw.text((im.size[0]/2.0-0.5*tsize[0],loY),str(xlabel),fill=fcolor,font=font)        
    if ylabel is not None:
        tsize = draw.textsize(str(ylabel),font=font)
        tim = Image.new('RGBA',(tsize[0],tsize[1]),bcolor)
        tdraw = ImageDraw.Draw(tim)        
        tdraw.text((0,0),str(ylabel),fill=fcolor,font=font)
        tim = tim.rotate(90)
        im.paste(tim,(0,int(im.size[1]/2.0-tsize[0])))
    if title is not None:
        tsize = draw.textsize(str(title),font=tfont)
        draw.text((im.size[0]/2.0-0.5*tsize[0],0),str(title),fill=fcolor,font=tfont)
    #----- Save image -----#
    im.save(filename)