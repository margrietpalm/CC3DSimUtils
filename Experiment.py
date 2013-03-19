#*******************************************************************#
# written by: Margriet Palm (08-03-10)		                        #
#-------------------------------------------------------------------#
#                                                                   #
# tools for setting up batch jobs for CC3D                          #
#*******************************************************************#
import xml.dom.minidom as xdm

class Experiment:
    """ The class Experiments holds the data for a CC3D experiment as an xml object. This objects is created from a template xml file. Elements that are already present in the template can be changed or removed and new elements can be added. When the xml is written to file all unchanged content of the template, including commented xml, is written to the new file.
    """       
    def __init__(self,templatefile):
        """ Instantiate experiment based on a given template 
        
        :param templatefile: path to xml template
        """
        self.template = xdm.parse(templatefile)
        self.cc3d = self.template.getElementsByTagName("CompuCell3D")[0]
        self.potts = self.template.getElementsByTagName("Potts")[0]
        self.plugins = self.template.getElementsByTagName("Plugin")
        self.steppables = self.template.getElementsByTagName("Steppable")
    
    #-- Private functions --#
    def _setNodeValue(self,node,value):
        """ Set value of an element.
        
        :param node: xml element
        :param value: value of xml element
        """        
        if len(node.childNodes) == 0:
            child = self.template.createTextNode(str(value))
            node.appendChild(child)
        else:
            node.childNodes[0].nodeValue = str(value)
    
    def _getElement(self,name,parent):
        """ Get element
            
        :param name: element name
        :param parent: parent node
        :return: xml element
        """        
        if (len(parent.getElementsByTagName(name))==0):
            element = self.template.createElement(name)
            parent.appendChild(element)
        else:
            element = parent.getElementsByTagName(name)[0]
        return element

    def _getPluginByName(self,name):
        """ Get CC3D plugin element by plugin name
           
        :param name: name of the plugin
        :return: xml element with plugin   
        """
        for plugin in self.plugins:
            if plugin.getAttribute("Name") == name:
                return plugin
        plugin = self.template.createElement(name)
        self.cc3d.appendChild(plugin)
        return plugin

    def _getSteppableByType(self,type):
        """ Get CC3D steppable element by steppable type
            
        :param type: type of the steppable
        :return: xml element with steppable specification  
        """
        for steppable in self.steppables:
            if steppable.getAttribute("Type") == type:
                return steppable
        steppable = self.template.createElement("Steppable")
        self._setAttribute("Type",type,steppable)
        self.cc3d.appendChild(steppable)
        return steppable

    def _setAttribute(self,name,value,parent):
        """ Set attribute value of an xml element
            
        :param name: attribute name
        :param value: attribute value
        :param parent: parent node
        """        
        if not parent.hasAttribute(name):
            parent.setAttributeNode(self.template.createAttribute(str(name)))
        parent.setAttribute(name,str(value))

    def setPottsProperty(self,tagname,attributes={},value=None):
        """ General function to set a parameter in the Potts element
            
        :param tagname: name of the element
        :param attributes: dictionary with attribute names as key and attribute values as values.
        :param value: element value
        """
        tag = self._getElement(tagname,self.potts)
        for key,val in attributes.iteritems():
            self._setAttribute(key, str(val), tag)
        if value:
            self._setNodeValue(tag,str(value))
    
    def deletePottsProperty(self,tagname):
        """ Remove element from the Potts element
            
        :param tagname: name of element
        """
        self.potts.removeChild(tagname)
    
    def setGenericPlugin(self,name,elements=[]):
        """ Set generic plugin. If the plugin is already in the template, it is updated. If not, the plugin is added to the model. This function does not support multi-level xml elements.
            
        :param name: plugin name
        :param elements: list of elements described by a dictionary: {'name':name,'value':val,'attributes':{}}
        """
        # element = {'name':name,'value':val,'attributes':{}}
        plugin = self._getPluginByName(name)
        for e in elements:
            element = self._getElement(e['name'], plugin)
            if 'value' in e:            
                self._setNodeValue(element, e['value'])
            if 'attributes' in e:
                for key,val in e['attributes'].iteritems():
                    self._setAttribute(key, val, element)
    
    def deletePlugin(self,name):
        """ Delete plugin
            
        :param name: plugin name        
        """
        self.cc3d.removeChild(self._getPluginByName(name))
    
    def setGenericSteppable(self,type,freq=None,elements=[]):     
        """ Set generic steppable. If the steppable is already in the template, it is updated. If not, the plugin is added to the model. This function does not support multi-level xml elements.
        
        :param type: steppable type
        :param freq: frequency
        :param elements: list of elements described by a dictionary with keys name, value and attributes        
        """
        # element = {'name':name,'value':val,'attributes':{}}
        steppable = self._getSteppableByType(type)
        if freq is not None:
            steppable.setAttribute("Frequency",str(freq))
        for e in elements:
            element = self._getElement(e['name'], steppable)
            if 'value' in e:
            #~ if len(e['value'] > 0):
                self._setNodeValue(element, str(e['value']))
            if 'attributes' in e:
                for key,val in e['attributes'].iteritems():
                    self._setAttribute(key, val, element)
    
    def deleteSteppable(self,type):
        """ Delete steppable
        
        :param type: steppable type        
        """
        self.cc3d.removeChild(self._getSteppableByType(type))    

    #-- Functions for specific settings, plugins, etc --#
    def setMultiCore(self,threads=1,cores=1):
        """ Set number of threads and cores for a simulation
            
        :param threads: number of threads per core
        :param cores: number of cores
        """        
        vpu = self._getElement("VirtualProcessingUnits",self._getElement("Metadata",self.cc3d))
        vpu.setAttribute("ThreadsPerVPU",str(threads))
        self._setNodeValue(vpu,cores)
    
    def setDebugFrequencyInMeta(self,freq):
        """ Set debug frequency
            
        :param frequency: frequency
        """        
        debug = self._getElement("DebugOutputFrequency",self._getElement("Metadata",self.cc3d))
        self._setNodeValue(debug,freq)
    
    def setSeed(self,seed):
        """ Set simulation seed
            
        :param seed: random seed
        """
        self.setPottsProperty("RandomSeed",value=seed)

    def setMCS(self,mcs):
        """ Set number of Monte Carlo steps
            
        :param mcs: number of MCS (note that mcs+1 appears in the xml)
        """
        self.setPottsProperty("Steps",value=mcs+1)
       
    def setTemp(self,T):
        """ Set temperature tag
            
        :param T: temperature 
        """
        self.template.getElementsByTagName("Temperature")[0].childNodes[0].nodeValue = str(T)

    def setMotility(self,celltype,T):
        """ Edit motility parameters per cell type
            
        :param celltype: name of the cell type
        :param T: motility
        """ 
        mot = self._getElement("CellMotility", self.potts)
        mottypes = mot.getElementsByTagName("MotilityParameters")
        typeFound = False
        for e in mottypes:
            if e.getAttribute("CellType") == celltype:
                e.setAttribute("Motility",str(T))
                typeFound = True
                break
        if not typeFound:
            element = self.template.createElement("MotilityParameters")
            element.setAttribute("CellType",str(celltype))
            element.setAttribute("Motility",str(T))
            mot.appendChild(element)            
    
    def addStatistics(self,freq,basename):
        """ Edit statistics save options
        
        :param freq: save frequency
        :param basename: basename for files (full path)
        """
        func = self._getElement("EnergyFunctionCalculator",self.potts)
        self._setAttribute("Type","Statistics",func)
        out = self._getElement("OutputFileName",func)
        self._setAttribute("Frequency", freq, out)
        self._setNodeValue(out, basename)

    def setVolume(self,celltype,vol,lam):
        """ Set volume per cell type
            
        :param celltype: cell type name
        :param vol: target volume
        :param lam: lambda volume
        """
        steppable = self._getPluginByName("VolumeFlex")
        elements = steppable.getElementsByTagName("VolumeEnergyParameters")
        typeFound = False
        for e in elements:
            if celltype == e.getAttribute("CellType"):
                self._setAttribute("TargetVolume",vol,e)
                self._setAttribute("LambdaVolume",str(lam),e)
                typeFound = True
                break
        if not typeFound:
            element = self.template.createElement("VolumeEnergyParameters")
            self._setAttribute("CellType",str(celltype),element)
            self._setAttribute("TargetVolume",str(vol),element)
            self._setAttribute("LambdaVolume",str(lam),element)
            steppable.appendChild(element)
    
    def setContact(self,_type1,_type2,J):
        """ Edit contact energy
        
        :param type1: name of the first cell type
        :param type2: name of the second cell type
        :param J: contact energy between type1 and type2        
        """
        plugin = self._getPluginByName("Contact")
        energy = plugin.getElementsByTagName("Energy")
        eFound = False
        for rule in energy: 
            type1 = rule.getAttribute("Type1")
            type2 = rule.getAttribute("Type2")
            if ((type1 in _type1) and (type2 in _type2)) or ((type1 in _type2) and (type2 in _type1)) :
                rule.childNodes[0].nodeValue = str(J)
                eFound = True
                break
        if not eFound:
            element = self.template.createElement("Energy")
            self._setAttribute("Type1",_type1,element)
            self._setAttribute("Type2",_type2,element)
            self._setNodeValue(element,J)
            plugin.appendChild(element)            
    
    def setSecretion(self,celltype,s,solver="FastDiffusionSolver2DFE"):
        """ Set secretion coefficient for specific cell type
            
        :param celltype: name of the cell type
        :param s: secretion coefficient
        :param solver: solver name
        """
        plugin = self._getSteppableByType("FastDiffusionSolver2DFE")
        elements = plugin.getElementsByTagName("Secretion")
        typeFound = False
        for e in elements:
            if e.getAttribute("Type") == celltype:
                self._setNodeValue(e,str(s))
                typeFound = True
                break
        if not typeFound:
            sdata = plugin.getElementsByTagName("SecretionData")[0]
            element = self.template.createElement("Secretion")
            element.setAttribute("Type",str(celltype))
            self._setNodeValue(element,str(s))
            sdata.appendChild(element)
    
    def setChemotaxis(self,celltype,towards,lam):
        """ Set chemotaxis for cell type
        
        :param celltype: name of the cell type
        :param towards: cell type towards chemotaxis occurs
        :param lam: chemotactic strength
        """
        plugin = self._getPluginByName("Chemotaxis")
        elements = plugin.getElementsByTagName("ChemotaxisByType")
        typeFound = False
        for e in elements:
            if e.getAttribute("Type") == celltype:
                e.setAttribute("ChemotactTowards",str(towards))
                e.setAttribute("Lambda",str(lam))
                typeFound = True
                break
        if not typeFound:
            element = self.template.createElement("ChemotaxisByType")
            element.setAttribute("Type",str(celltype))
            e.setAttribute("ChemotactTowards",str(towards))
            e.setAttribute("Lambda",str(lam))
            plugin.getElementsByTagName("ChemicalField")[0].appendChild(element)
    
    def write(self,filename):
        """ Save xml to file 
        
        :param filename: filename of new xml
        """       
        f = open(filename,'w')
        self.template.writexml(f)
        f.close()     