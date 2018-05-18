###############################################################################################
#    Varna SVG color script
#
#      Author: Gregg Rice
#              gmr@unc.edu
#               Steve modified to plot T1 cleavage data
#               Steve modified to plot additional custom nuc colors and lines between nucs
#               Tony modified it to simplify class structure (reduce redundnacy) and fix bugs
#
#
# Affiliation: Weeks Lab, UNC Chemistry
#
#        Date: 11.07.12+
#     Version: 0.92+
#
# released under GPL 2.0
###############################################################################################


import math,sys
import numpy as np

import RNAtools


###############################################################################################
# End of ssDiagram class
#    Combines both varna and xrna classes into one, since they replicate the same functionality
###############################################################################################

# known bugs:
# - 23S rRNA structure fails to load. There are two RNAs present in this fileIn
#   and that causes an issue

class ssDiagram:

    def __init__(self, fIN, filetype = None):
        """
        Init the secondary structure diagram object
        Will infer object type from file extension, or can specified manually using filetype
        """

        if filetype is None:
            ext = fIN.split('.')[-1].upper()
            if ext == 'VARNA':
                filetype = 'VARNA'
            elif ext == 'XRNA':
                filetype = 'XRNA'
            else:
                sys.exit("Extension %s of file %s is not recognized" % (ext, fIN))

        self.filetype = filetype
        self.name = fIN

        if filetype == 'VARNA':
            self.num, self.pairMap, self.baseMap,self.center,self.shape,self.ct,self.scale,self.period= self.readVarna(fIN)
        else:
            self.num, self.pairMap, self.baseMap,self.shape,self.ct,self.scale, self.period = self.parseXRNA(fIN)
            self.center = self.pairMap


        self.numMap = self.baseMap
        self.circleColors = False
        self.colorIDs = []
        self.extraLines = False
        self.extraLineData = []
        self.threshA = 0.025
        self.threshB = 0.055
        self.diff = False
        self.colorFunction = self.getSHAPEcolor


    def __str__(self):
        a = '{ Name= %s, Length = %s}' % (self.name, str(len(self.num)))


    def setdiff(self):
        self.diff = True
        self.colorFunction = self.getDiffSHAPEcolor


    def readVarna(self,x):
        import xml.etree.ElementTree as ET
        tree = ET.parse(x)
        root = tree.getroot()

        #initialize some arrays
        offset = 15/2.
        num,bases,x_pos,y_pos,shape,x_cen,y_cen = [],[],[],[],[],[],[]

        #read the varna xml file, nucleotides are under bases
        for nt in root.findall('./RNA/bases/nt'):
            num.append(int(nt.get('num')))
            shape.append(float(nt.get('val')))
            base = nt.find('base').text
            for i in nt:
                if i.get('r')=='pos':
                    x_pos.append(float(i.get('x')))
                    y_pos.append(float(i.get('y')))
                    bases.append(base)
                if i.get('r')=='center':
                    x_cen.append(float(i.get('x')))
                    y_cen.append(float(i.get('y')))
        #determine offset
        x_pos,y_pos = np.array(x_pos),np.array(y_pos)

        vec = ((x_pos[0]-x_pos[1])**2+(y_pos[0]-y_pos[1])**2)**0.5
        offset = vec/4

        #transpose coords to positive space
        xoff = abs(min(x_pos))+offset
        yoff = abs(min(y_pos))+offset
        x_pos = x_pos + xoff
        y_pos = y_pos + yoff

        x_cen += xoff
        y_cen += yoff
        center = list(zip(x_cen,y_cen))

        #make expected arrays
        coord = list(zip(x_pos,y_pos))
        basemap = list(zip(bases,zip(x_pos,y_pos+offset)))

        #read varna xml file for pairing information
        ct = np.zeros(len(num))
        for pair in root.findall('./RNA/BPs/bp'):
            p5,p3 = int(pair.get('part5')),int(pair.get('part3'))
            #print p5+1,p3+1
            ct[p5] = p3+1
            ct[p3] = p5+1
        ct = list(map(int,ct))

        #get the number period
        period = int(root.find('config').get('numperiod'))
        return num,coord,basemap,center,shape,ct,offset,period




    def parseXRNA(self,x):
        import xml.etree.ElementTree as ET
        tree = ET.parse(x)
        root = tree.getroot()

        nucList = root.findall('./Complex/RNAMolecule/')
        nucLists = []
        for i in nucList:
            #print i.tag
            if i.tag == 'NucListData':nucLists.append(i)

        #print nucLists
        startNT = int(nucLists[0].get('StartNucID'))
        #print startNT
        num = []
        bases,x_pos,y_pos,x_cen,y_cen = [],[],[],[],[]
        for nt in nucLists[0].text.split('\n'):
            if nt == '':continue
            line = nt.split()
            num.append(startNT)
            startNT+=1
            bases.append(line[0]),x_pos.append(float(line[1])),y_pos.append(float(line[2]))

        #print x_pos
        #determine offset
        x_pos,y_pos = np.array(x_pos),-1*np.array(y_pos)

        vec = ((x_pos[0]-x_pos[1])**2+(y_pos[0]-y_pos[1])**2)**0.5
        offset = vec/1.5

        #transpose coords to positive space
        xoff = abs(min(x_pos))+offset
        yoff = abs(min(y_pos))+offset
        x_pos = x_pos + xoff
        y_pos = y_pos + yoff
        y_pos=2.3*y_pos
        x_pos=2.3*x_pos

        x_cen += xoff
        y_cen += yoff
        center = list(zip(x_cen,y_cen))

        #make expected arrays
        coord = list(zip(x_pos,y_pos))
        basemap = list(zip(bases,list(zip(x_pos,y_pos+offset))))
        #print basemap

        shape = np.zeros(len(num))
        clist = {'bbbbbb':-999,'999999':-999,'ff0000':1.0, '0':0.0, '000000':0.0, '1c75bc':-0.45,'00ff00':0.45, 'ff9900':0.45, 'f57e1f':0.45}
        for shapeLine in root.findall('./Complex/RNAMolecule/Nuc'):
            nucRange = shapeLine.get('RefIDs')
            preColor = shapeLine.get('Color')
            if not preColor:continue
            try:nucColor = clist[shapeLine.get('Color')]
            except:nucColor = 0.0
            if not nucRange:continue
            for i in nucRange.split(','):
                if len(i.split('-'))==1:
                    try:shape[int(i)-1]=nucColor
                    except:pass
                else:
                    line = i.split('-')
                    line = list(map(int,line))
                    for j in range(line[0],line[1]+1):
                        try:shape[j-1]=nucColor
                        except:pass

        shape = list(map(float,shape))
        period = 20
        #get pairing informationo
        ct = np.zeros(len(num))
        for pair in root.findall('./Complex/RNAMolecule/BasePairs'):
            pairStart,pairLength,pairEnd = pair.get('nucID'),pair.get('length'),pair.get('bpNucID')
            for nt in range(int(pairLength)):
                p5 = int(pairStart)+nt
                p3 = int(pairEnd)-nt
                try:
                    ct[p5-1] = p3
                    ct[p3-1] = p5
                except:pass
        ct = list(map(int,ct))

        return num,coord,basemap,shape,ct,offset,period



    def readSHAPE(self,z):
        self.shape = RNAtools.readSHAPE(z)


    def parseCircleColorLine(self, line):
        """
        first column is nuc number, second and third columns are color specifications
        available color formats are as follows:
        0/1/2          : the values will be assigned to 'white', 'red, nofill', 'red, fill' respectively
        float [0/1]    : the values will be converted to colors based on SHAPE or differential SHAPE color scheme
                         optional 0/1 indicates whether to fill (defaults to 0 -- nofill)
        color [color]  : Colorname must be SVG color. Optional second color specifies the fill (defaults to white)
        r,g,b [r,g,b]  : Input r,g,b values. Optional second color specifies the fill (defaults to white)
        """

        defaultmap = ['white', 'red', 'red']
        color = ['white', 'white']

        spl = line.split()

        try:
            cindex = int(spl[1])
            color = [defaultmap[cindex], defaultmap[cindex-1]]

        except ValueError:
            try:
                color[0] = float(spl[1])
                if len(spl) == 3 and spl[2]=='1':
                    color[1] = color[0]

            except ValueError:
                for i in range(1,len(spl)):
                    if ',' in spl[i]:
                        color[i-1] = 'rgb('+spl[i]+')'
                    else:
                        color[i-1] = spl[i]

        except:
            sys.exit("Improperly formatted circle color :: %s" % line)

        return tuple(color)


    def readCircleColors(self, z, ):
        self.circleColors = True

        with open(z, 'rU') as f:
            for line in f:
                # parse out any headers
                if line[0] == '#':
                    continue

                self.colorIDs.append(self.parseCircleColorLine(line))


    def readExtraLines(self,z):
        self.extraLines = True
        splitZ = z.split(',')
        filename=splitZ[0]

        if len(splitZ)>1:
            self.threshA = float(splitZ[1])
            self.threshB = float(splitZ[2])

        extra = []
        with open(filename, 'rU') as f:
            for line in f:
                spl = line.split()
                try:
                    extra.append((int(spl[0]), int(spl[1]), float(spl[2])))
                except (ValueError, IndexError):
                    pass

        self.extraLineData = extra


    def getSHAPEcolor(self, x, suppress=False):
        """
        return shape color if x is expected float, 'black' if x is None, and the
          same input string otherwise
          if suppress is True, will suppress black colors by setting them to white
        """

        if isinstance(x, str):
            try:
                x = float(x)
            except ValueError:
                return x

        elif x is None:
            return 'rgb(1,0,0)'

        if x < -4: col = '160,160,160'  # grey
        elif x > 0.85: col = '255,0,0'  # red
        elif 0.85 >= x >0.4: col ='255,164,26' # orange
        elif not suppress and 0.4 >= x > -4: col = '1,0,0' # black
        else: col = '255,255,255'

        return 'rgb(%s)' % col

    def getDiffSHAPEcolor(self, x, suppress=False):
        """
        return diff shape color if x is expected float, 'white' if x is None, and the
          and return the same input string otherwise
        """

        if isinstance(x, str):
            try:
                x = float(x)
            except ValueError:
                return x

        elif x is None:
            return 'rgb(1,0,0)'

        if x < -500: col = '160,160,160'
        elif x <= -0.2: col = '41,171,226'
        elif x >= 0.2: col = '0,210,0'
        elif not suppress: col = '1,0,0'
        else: col = '255,255,255'

        return 'rgb(%s)' % col



###############################################################################################
# Start Enzyme class
###############################################################################################


class Enzyme:

    def __init__(self, filePath):
        self.arrowSizes = []
        self.nums = []
        if filePath != "":
            self.nums, self.arrowSizes = self.parseEnzymeFile(filePath)

    def parseEnzymeFile(self,filePath):
        fileIn = open(filePath, "rU")
        fileIn.readline() # skip header line
        nums = []
        arrowSizes = []
        for line in fileIn:
            splitLine = line.strip().split()
            nums.append(int(splitLine[0]))
            arrowSizes.append(int(splitLine[1]))
        fileIn.close()
        return nums, arrowSizes

    def rotatePoint(self, coords, theta):
        x = coords[0]
        y = coords[1]
        xRot = x*math.cos(-theta)+ y*math.sin(-theta)
        yRot = -x*math.sin(-theta) + y*math.cos(-theta)
        #print "x,y= %f, %f\ttheta= %f\t xrot,yrot= %f, %f"%(x,y,theta,xRot,yRot)
        return [xRot,yRot]

    def translatePoint(self, coords, offset):
        return [coords[i] + offset[i] for i in range(len(coords))]

    def midpoint(self, c1, c2):
        return [c1[0]+(c2[0]-c1[0])/2, c1[1]+(c2[1]-c1[1])/2]


    def calcArrowCoords(self, origin, h, w, theta):
        # define initial wedge V
        left = [-w/2,h]
        right = [w/2,h]
        # rotate initial wedge about 0,0
        leftRot = self.rotatePoint(left,theta)
        rightRot = self.rotatePoint(right,theta)
        # translate to given origin
        left = self.translatePoint(leftRot, origin)
        right = self.translatePoint(rightRot, origin)
        # return three coords
        return origin, left, right

    def arrowToString(self, pt1, pt2, pt3):
        svgString = '<polygon fill="rgb(80,100,255)" stroke="green" stroke-width="0" points="%f,%f  %f,%f %f,%f" />'%(pt1[0],pt1[1],pt2[0],pt2[1],pt3[0],pt3[1])
        return svgString

    def drawArrows(self, varna):
        heights = [0,20,40,80]
        width = 10
        arrowString = ""
        for i in range(len(self.nums)):
            num = self.nums[i]
            arrowSize = self.arrowSizes[i]
            loc1 = [varna.baseMap[num-1][1][0],4.0+varna.baseMap[num-1][1][1]]
            loc2 = [varna.baseMap[num][1][0],4.0+varna.baseMap[num][1][1]]
            loc = self.midpoint(loc1, loc2) # put arrow 3prime of cleaved nuc
            #print loc
            # 0=no arrow, 1=lo, 2=mid, 3=high
            if arrowSize != 0:
                height = heights[arrowSize]

                #pt = findFarPoint(num-1,varna)
                #xDiff = pt[0] - loc[0]
                #yDiff = pt[1] - loc[1]
                #theta = math.atan2(yDiff, xDiff)

                # assuming clockwise nuc order, find normal angle
                diff = [loc2[n]-loc1[n] for n in [0,1]]
                theta = math.atan2(diff[1], diff[0]) + math.pi

                coords = self.calcArrowCoords(loc,height,width,theta)
                arrowString += self.arrowToString(coords[0],coords[1],coords[2])

        return arrowString


###############################################################################################
# General functions below
###############################################################################################


def evalStructures(rna1,rna2):
    # Returns shared basepairs, those only in rna1, those only in rna2
    # n number of accepted pairs x2, p number of predicted pairs x2
    n,p = 0,0
    shared,acceptedOnly,predictedOnly = [],[],[]
    #print(rna1.ct)
    rna1.ct = list(rna1.ct)
    rna2.ct = list(rna2.ct)
    for i in range(len(rna1.ct)):
        #clean out duplicates and nonbasepairs
        if rna1.ct[i] == 0 and rna2.ct[i] ==0:continue

        #count for Sens,PPV
        if rna1.ct[i] != 0: n+=1
        if rna2.ct[i] != 0: p+=1

        #shared bps
        if rna1.ct[i] == rna2.ct[i]:
            if rna1.ct[i] < rna1.num[i]:continue
            if rna2.ct[i] < rna2.num[i]:continue
            shared.append((rna1.num[i],rna1.ct[i]))
            continue

        if rna1.ct[i] != 0 and rna1.ct[i] < rna1.num[i]:
            acceptedOnly.append((rna1.num[i],rna1.ct[i]))
        if rna2.ct[i] != 0 and rna2.ct[i] < rna2.num[i]:
            predictedOnly.append((rna2.num[i],rna2.ct[i]))

    return acceptedOnly,predictedOnly,shared,len(shared)/(float(n)/2),len(shared)/(float(p)/2)


def offPoint(p,q,r):
    p,q = np.array(list(p)),np.array(list(q))
    v_u = (q-p)/(np.sum((q-p)**2)**0.5)
    return v_u*r+p

def newLines(pointSet,locMap,r,struct=False):
    a = []
    for i,j in pointSet:
        p1 = offPoint(locMap[i-1],locMap[j-1],r)
        p2 = offPoint(locMap[j-1],locMap[i-1],r)
        #check noncanonical
        #if struct:
            #canonical
            #print struct.baseMap[i-1][0],struct.baseMap[j-1][0]
        a.append((p1,p2))
    return a



def drawBases(varna, fontsize = 24):

    bases = varna.baseMap
    shape = varna.shape

    line = ''

    max_x, max_y = 0.0, 0.0

    for i in range(len(bases)):

        color = varna.colorFunction(shape[i])

        #if varna.diff and color != 'black':
        #    line += '<text x="%s" y="%s" text-anchor="middle" font-family="Sans-Serif" font-weight="bold" font-size="%d" stroke="%s" fill="%s" >%s</text>'\
        #            % (bases[i][1][0], bases[i][1][1], fontsize+2, color, color, bases[i][0])
        #else:
        if max_x < bases[i][1][0]:
            max_x = float(bases[i][1][0])
        if max_y < bases[i][1][1]:
            max_y = float(bases[i][1][1])

        line += '<text x="%s" y="%s" text-anchor="middle" font-family="Sans-Serif" font-weight="bold" font-size="%d" fill="%s" >%s</text>' \
                % (bases[i][1][0], bases[i][1][1]+4.0, fontsize, color, bases[i][0])
    max_x += 60
    max_y += 60
    return line, max_x, max_y


def findFarPoint(pt,varna):
    #finds a good direction to draw the number label
    # goes through the distance pairs, finds all nts within ~ 63pts and finds the center of mass
    x,y = [],[]
    for i,j in varna.pairMap:
        x.append(i),y.append(j)
    x,y = np.array(x),np.array(y)
    point = np.array((x[pt],y[pt]))
    #print point
    dist = np.sum((point-np.transpose((x,y)))**2,axis=1)
    #print len(dist[dist<5000]),'#'
    #print (x,y)
    cutoff=4000
    length = len(x[dist<cutoff])
    centerMass = np.sum(x[dist<cutoff])/length, np.sum(y[dist<cutoff])/length
    #print str(np.array(centerMass))
    return np.array(centerMass)


def findCWNormalPoint(pt,varna):
    x,y = [],[]
    for i,j in varna.pairMap:
        x.append(i),y.append(j)
    x,y = np.array(x),np.array(y)
    point = np.array((x[pt],y[pt]))
    #print "point: "+str(point)
    try:
        pointAfter = np.array((x[pt+1],y[pt+1]))
    except IndexError:
        pointAfter = np.array((x[pt]+1,y[pt]))
    # assuming clockwise nuc order, find normal angle
    diff = [pointAfter[n]-point[n] for n in [0,1]]
    #print "diff: "+str(diff)
    theta = math.atan2(diff[1], diff[0]) + math.pi/2
    #print "theta: "+str(theta)
    distance = 20
    newX = point[0]+math.cos(theta)*distance
    newY = point[1]+math.sin(theta)*distance
    newPoint = np.array((newX,newY))
    #print "newPoint: "+str(newPoint)+"\n"
    return newPoint

def drawNums(varna,offset):
    period = varna.period
    #draw over i
    nums = []
    lines = []
    for i in [1]+list(range(0,len(varna.num),period))[1:]+[len(varna.num)]:
        key = i-1
        if varna.filetype == 'XRNA':
            center = findFarPoint(key,varna)
            #center = findCWNormalPoint(key,varna)
        else:
            center = np.array(varna.center[key])

        a = np.array(varna.pairMap[key])
        base = np.array(varna.pairMap[key])
        #print "base: "+str(base)
        #print "center: "+str(center)
        norm = np.sum((base-center)**2)**0.5
        u_vec = (base-center)/norm*varna.scale*7 + base
        nums.append((str(i),list(map(float,u_vec))))

        p1 = offPoint(map(float,u_vec),map(float,base),varna.scale*2)
        p2 = offPoint(map(float,base),map(float,u_vec),varna.scale*2)
        lines.append((p1,p2))

    #add lines connecting base and letters
    line = processConnect(lines,(3,0,0),lineMap=True)
    #add numbering
    for i in list(range(len(nums))):
        line += '<text x="%s" y="%s" text-anchor="middle" font-family="Sans-Serif" font-weight="bold" font-size="18" fill="rgb(0,1,0)" >%s</text>' \
                % (nums[i][1][0],nums[i][1][1]+varna.scale,str(int(nums[i][0])+offset))
    return line


def drawOutline(varna):
    outlineString = '<polyline points= "'
    for nuc in varna.baseMap:
        pt = nuc[1]
        outlineString += '%f,%f ' % (pt[0],pt[1]-4)
    outlineString += '" stroke="rgb(200,200,200)" stroke-width="3.0" fill="none"/>'
    return outlineString


def drawCircles(varna):

    outlineString = ''

    for i in range(len(varna.baseMap)):
        if i < len(varna.colorIDs):
            pt = varna.baseMap[i][1]

            col = [varna.colorFunction(x, True) for x in varna.colorIDs[i]]

            outlineString += '<circle cx="%f" cy="%f" r="14" stroke="%s" stroke-width="2" fill="%s"/>'% (pt[0], pt[1]-4, col[0], col[1])

    return outlineString


def processConnect(pairlist,color,dashed=False,lineMap=False,strokewidth=4.0):
    def rbgform(x):
        return ','.join(map(str, x))

    out = ''
    for i,j in pairlist:

        line = '<line x1="%s" y1="%s" x2="%s" y2="%s" ' % (i[0],i[1],j[0],j[1])
        if lineMap == True:
            line += 'stroke="rgb(%s)" stroke-width="0.9" opacity="0.95" />' % rbgform(color)
            out+=line
            continue

        if dashed==False:line += 'stroke="rgb(%s)" stroke-width="%.1f" opacity="0.95" />' % (rbgform(color),strokewidth)
        if dashed==True:line += 'stroke="rgb(%s)" stroke-width="1.2" opacity="0.85" />' % rbgform(color)
        out+=line
    return out


def drawExtraLines(varna):

    out = ""
    if varna.extraLines == False:
        return out


    bins = [-1.0, -varna.threshB, -varna.threshA, 0, varna.threshA, varna.threshB, 1.0]
    colorList = [(44,123,182), (44,123,182), (171,217,233), (255,255,255),
                 (253,174,97), (215,25,28), (215,25,28)]
    gradient = [False, True, False, False, True, False]

    for fields in varna.extraLineData:

        fromNuc = fields[0]-1
        fromNucCoord = varna.baseMap[fromNuc][1]
        toNuc = fields[1]-1
        toNucCoord = varna.baseMap[toNuc][1]

        corrCoeff = fields[2]

        # filter out low correlations
        if abs(corrCoeff) < varna.threshA or abs(fromNuc-toNuc) < 6:
            continue


        for i in range(len(bins)-1):
            if bins[i]<corrCoeff<bins[i+1]:

                if gradient[i]:
                    col = colorGrad(corrCoeff, colorList[i], colorList[i+1], bins[i], bins[i+1])
                else:
                    col = colorList[i]

                line = processConnect([( ([fromNucCoord[0],fromNucCoord[1]]),
                                         ([toNucCoord[0],toNucCoord[1]])  )],
                                    col,strokewidth=5)
                out+=line
                break


    #print out
    return out


def colorGrad(value, colorMin, colorMax, minValue, maxValue):
    """
    returns a middle rgb color value based on the distance between max and min
    """
    c_range = abs(maxValue - minValue)
    v       = value - min(maxValue,minValue)
    v_pct = v/c_range

    #v_pct *= v_pct

    #print value, v_pct, colorMin, colorMax

    i= v_pct*(colorMax[0] - colorMin[0]) + colorMin[0]
    j= v_pct*(colorMax[1] - colorMin[1]) + colorMin[1]
    k= v_pct*(colorMax[2] - colorMin[2]) + colorMin[2]

    alphaThresh = 0.25
    if v_pct > alphaThresh: alpha = 1.0
    else:
        alpha = v_pct / alphaThresh

    if value > maxValue: return colorMax #, alpha
    elif value < minValue: return colorMin #, alpha
    else: return map(int,(i,j,k)) #, alpha



def parseArgs():
    import argparse
    prs = argparse.ArgumentParser(description='Colors and optionally compares a VARNA or XRNA file with a reference ct and .SHAPE file')
    prs.add_argument('input',action='store',type=str,help='input file')
    prs.add_argument('output',action='store',type=str,help='output .svg file')
    prs.add_argument('--ct',action='store',type=str,help='compare structure with a reference ct file')
    #prs.add_argument('-x','--xrna',action='store_true',default=False,help='changes input file type to XRNA')
    prs.add_argument('-e','--enzyme',action='store',type=str,help='draw enzymatic cleavage data from file')
    prs.add_argument('-d','--diff',action='store_true',default=False,help='changes chemical probing type to differential, coloring cutoffs +=0.3')
    prs.add_argument('-c','--colors',action='store',type=str,help='color behind nucs with custom colors')
    prs.add_argument('-l','--lines',action='store',type=str,help='draw additional lines between certain nucs')
    prs.add_argument('-s','--shape',action='store',type=str,help='overide stored chemical probing values from varna')
    prs.add_argument('--offset', action='store',type=int,default=0,help='numbering ofset, adds this to the numbering in the file')
    prs.add_argument('--switch',action='store_true',default=False,help='reverse the pairing coloring scheme')
    o=prs.parse_args()
    return o

def hasCT(correct,varna,switch=False):
    ac,pred,both,s,p = evalStructures(correct,varna)
    print('PPV: ', round(p,2),'SENS: ',round(s,2))

    setScale = varna.scale*2
    if varna.filetype == 'XRNA':
        setScale = varna.scale*2
    #make lines
    both_pt = newLines(both,varna.pairMap,setScale)
    pred_pt = newLines(pred,varna.pairMap,setScale)
    ac_pt = newLines(ac,varna.pairMap,setScale)

    #define colors
    green,red,purple = (0,50,0),(100,0,0),(39,17,56)
    if switch:
        # switch the predicted and accepted colors
        red,purple=(39,17,56),(100,0,0)

    #draw lines
    drawnLines = processConnect(pred_pt,purple) + processConnect(ac_pt,red,dashed=True) + processConnect(both_pt,green)
    return drawnLines

def noCT(varna):
    ac,pred,both,s,p = evalStructures(varna,varna)

    setScale = varna.scale*2
    if varna.filetype == 'XRNA':
        setScale = varna.scale*2
    both_pt = newLines(both,varna.pairMap,setScale,struct=varna)
    black = (0,0,0)
    drawnLines = processConnect(both_pt,black)
    return drawnLines

def main():

    arg = parseArgs()

    #read input file
    svgStruct = ssDiagram(arg.input)

    #overwrite stored chemical probing values
    if arg.shape:
        svgStruct.readSHAPE(arg.shape)

    if arg.diff:
        svgStruct.setdiff()

    # do custom colors if given
    if arg.colors:
        svgStruct.readCircleColors(arg.colors)

    # do extra lines
    if arg.lines:
        svgStruct.readExtraLines(arg.lines)

    # read enzymatic cleavage data
    if arg.enzyme:
        enzyme = Enzyme(arg.enzyme)

    # draw nucleotide outline
    outline = drawOutline(svgStruct)

    # draw circles behind nucleotides
    circles = drawCircles(svgStruct)
    #print svgStruct
    #print circles

    # extra lines
    extraLines = drawExtraLines(svgStruct)

    # draw lines and compare structures if needed
    if arg.ct:
        correct = RNAtools.CT(arg.ct)
        drawnLines = hasCT(correct,svgStruct,arg.switch)
    else:
        drawnLines = noCT(svgStruct)

    letters, max_x, max_y = drawBases(svgStruct)
    letters += drawNums(svgStruct,arg.offset)
    #construct header and letters
    header = '<svg width="{0}" height="{1}" version="1.1" xmlns="http://www.w3.org/2000/svg">'.format(max_x, max_y)
    background = '<rect x="0" y="0" width="100%" height="100%" style="fill:rgb(255,255,255)"/>'.format(max_x, max_y)
    if arg.enzyme:
        arrows = enzyme.drawArrows(svgStruct)
    else:
        arrows = ""

    #write file
    out = header + background + outline + arrows + circles + letters + drawnLines + extraLines + '</svg>'
    w = open(arg.output,"w")
    w.write(out)
    w.close()

if __name__ == "__main__":
    main()
