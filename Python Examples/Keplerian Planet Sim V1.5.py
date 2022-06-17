#https://ssd.jpl.nasa.gov/txt/aprx_pos_planets.pdf
from math import pi, sqrt
import math

def sin(x):
    return math.sin(math.radians(x))
def cos(x):
    return math.cos(math.radians(x))


#Planet Class Defintion:
########################
class Planet:
    def __init__(self,
                 a0, e0, I0, L0, W0, o0,
                 aR, eR, IR, LR, WR, oR):
        self.a0, self.e0, self.I0, self.L0, self.W0, self.o0 = a0, e0, I0, L0, W0, o0
        self.aR, self.eR, self.IR, self.LR, self.WR, self.oR = aR, eR, IR, LR, WR, oR
        self.extra = (0,0,0,0)
    def compute(self, T):
        a = self.a0 + self.aR*T
        e = self.e0 + self.eR*T
        I = self.I0 + self.IR*T
        L = self.L0 + self.LR*T
        W = self.W0 + self.WR*T
        o = self.o0 + self.oR*T

        b, c, s, f = self.extra
        
        w = W - o
        M = L - W + b*T**2 + c*cos(f*T) + s*sin(f*t)

        self.params = a,e,I,L,W,o,w

        return self.computeFromM(M)

    def computeFromM(self, M):
        a,e,I,L,W,o,w = self.params
        
        M = ((M+180)%360)-180
        En = M + (180/pi)*e*sin(M)
        iterate = True
        while iterate:
            dM = M - (En - (180/pi)*e*sin(En))
            dE = dM / (1 - e*cos(En))
            En += dE
            if dE <= 1e-6:
                iterate = False
        E = En

        oX = a*(cos(E)-e)
        oY = a*sqrt(1-e**2)*sin(E)
        oZ = 0

        eX = (cos(w)*cos(o)-sin(w)*sin(o)*cos(I))*oX + (-sin(w)*cos(o)-cos(w)*sin(o)*cos(I))*oY
        eY = (cos(w)*sin(o)+sin(w)*cos(o)*cos(I))*oX + (-sin(w)*sin(o)+cos(w)*cos(o)*cos(I))*oY
        eZ = (sin(w)*sin(I))*oX                      + (cos(w)*sin(I))*oY

        obliquity = 23.43928
        eqX = eX
        eqY = cos(obliquity)*eY - sin(obliquity)*eZ
        eqZ = sin(obliquity)*eY + cos(obliquity)*eZ

        return eqX, eqY, eqZ


#Planet Object Initialisation From Parameters:
##############################################
Mercury = Planet( 0.38709843, 0.20563661, 7.00559432,   252.25166724, 77.45771895, 48.33961819,
                  0.00000000, 0.00002123,-0.00590158,149472.67486623,  0.15940013, -0.12214182)
Venus   = Planet( 0.72332102, 0.00676399, 3.39777545,   181.97970850,131.76755713, 76.67261496,
                 -0.00000026,-0.00005107, 0.00043494, 58517.81560260,  0.05679648, -0.27274174)
Earth   = Planet( 1.00000018, 0.01673163,-0.00054346,   100.46691572,102.93005885, -5.11260389,
                 -0.00000003,-0.00003661,-0.01337178, 35999.37306329,  0.31795260, -0.24123856)
Mars    = Planet( 1.52371243, 0.09336511, 1.85181869,    -4.56813164,-23.91744784, 49.71320984,
                  0.00000097, 0.00009149,-0.00724757, 19140.29934243,  0.45223625, -0.26852431)
Jupiter = Planet( 5.20248019, 0.04853590, 1.29861416,    34.33479152, 14.27495244,100.29282654,
                 -0.00002864, 0.00018026,-0.00322699,  3034.90371757,  0.18199196,  0.13024619)
Saturn  = Planet( 9.54149883, 0.05550825, 2.49424102,    50.07571329, 92.86136063,113.63998702,
                 -0.00003065,-0.00032044, 0.00451969,  1222.11494724,  0.54179478, -0.25015002)
Uranus  = Planet(19.18797948, 0.04685740, 0.77298127,   314.20276625,172.43404441, 73.96250215,
                 -0.00020455,-0.00001550,-0.00180155,   428.49512595,  0.09266985,  0.05739699)
Neptune = Planet(30.06952752, 0.00895439, 1.77005520,   304.22289287, 46.68158724,131.78635853,
                  0.00006447, 0.00000818, 0.00022400,   218.46515314,  0.01009938, -0.00606302)
Pluto   = Planet(39.48686035, 0.24885238,17.14104260,   238.96535011,224.09702598,110.30167986,
                  0.00449751, 0.00006016, 0.00000501,   145.18042903, -0.00968827, -0.00809981)


Jupiter.extra = (-0.00012452, 0.06064060,-0.35635438,38.35125000)
Saturn.extra  = ( 0.00025899,-0.13434469, 0.87320147,38.35125000)
Uranus.extra  = ( 0.00058331,-0.97731848, 0.17689245, 7.67025000)
Neptune.extra = (-0.00041348, 0.68346318,-0.10162547, 7.67025000)
Pluto.extra   = (-0.01262724, 0         , 0         , 0         )

#Camera Rotation Control and Transformation:
############################################
cameraRot = (0,0,0)
mousePressed = False
lastPos = [0, 0]
def motion(event):
    global lastPos, cameraRot
    if mousePressed:
        x,y,z = cameraRot
        cameraRot = (x+((event.y - lastPos[1]) * 0.01), y-((event.x - lastPos[0]) * 0.01), z)
        lastPos = [event.x, event.y]
def mouseDown(event):
    global mousePressed, lastPos
    mousePressed = True
    lastPos = [event.x, event.y]
def mouseUp(event):
    global mousePressed
    mousePressed = False
def scroll(event):
    global screenR, factor
    direction = int(event.delta/abs(event.delta))
    screenR += -direction*0.05*screenR
##    screenR = max([screenR, 0.1])
    factor = (min([screenW, screenH])/2)/screenR


def mMProduct(M, N):
    if len(M[0]) != len(N):
        raise ValueError("Matrices product is undefined.", M, N)
    return [[sum([M[r][i]*N[i][c] for i in range(len(M[0]))]) for c in range(len(N[0]))] for r in range(len(M))]
def rotYmatrix(yRotation):
    return [[ math.cos(yRotation), 0, math.sin(yRotation)],
            [                   0, 1,                   0],
            [-math.sin(yRotation), 0, math.cos(yRotation)]]
def rotXmatrix(xRotation):
    return [[1,                   0,                   0],
            [0, math.cos(xRotation),-math.sin(xRotation)],
            [0, math.sin(xRotation), math.cos(xRotation)]]
def rotZmatrix(zRotation):
    return [[math.cos(zRotation),-math.sin(zRotation), 0],
            [math.sin(zRotation), math.cos(zRotation), 0],
            [                  0,                   0, 1]]
def rotXYZ(rotation):
    return mMProduct(mMProduct(rotZmatrix(rotation[2]), rotXmatrix(rotation[0])), rotYmatrix(rotation[1]))
def matrixNoteSwap(V):
    return [[c[r] for c in V] for r in range(len(V[0]))]

def rotateToCamera(coords):
    global cameraRot
    rotMatrix = rotXYZ([-j for j in list(cameraRot)])
    coords = tuple(matrixNoteSwap(mMProduct(rotMatrix, matrixNoteSwap([list(coords)])))[0])
    return coords


#Tkinter Setup and Drawing:
###########################
screenR = 6.5#AU
screenW, screenH = 800, 800
factor = (min([screenW, screenH])/2)/screenR

from tkinter import *
root = Tk()
canvas = Canvas(root, width=screenW, height=screenH, bg="#000000")
canvas.bind('<Motion>', motion)
canvas.bind('<ButtonPress-1>', mouseDown)
canvas.bind('<ButtonRelease-1>', mouseUp)
canvas.bind("<MouseWheel>", scroll)
canvas.pack()

import datetime
epoch = datetime.date(2000,1,1)

def drawPlanet(coords, r=5, colour="red", outline=""):
    X,Y = spaceToScreen(coords)
    canvas.create_oval(X-r,Y-r,X+r,Y+r, fill=colour, outline=outline)

def spaceToScreen(coords):
    x,y,z = rotateToCamera(coords)
    X = x*factor + screenW/2
    Y = y*factor + screenH/2
    return X,Y
    
def displayPlanet(planet, step=10, colour="#ff0000"):
    RGB = [colour[1:3], colour[3:5], colour[5:7]]
    orbitColour = "#"+"".join([hex(min([255, int(i,base=16)+64]))[2:].zfill(2) for i in RGB])
    
    planet.compute(t)
    points = [spaceToScreen(planet.computeFromM(i)) for i in range(0,360,step)]
    canvas.create_line(points+[points[0]], fill=orbitColour, width=1.0)
    drawPlanet(planet.computeFromM(0),   outline=orbitColour, r=3, colour=orbitColour)
    drawPlanet(planet.computeFromM(180), outline=orbitColour, r=3, colour="black")
    drawPlanet(planet.compute(t), r=5, colour=colour)

def displayScale():
##    order = int(math.log10(20/factor))
##    length = 10**order * factor
    size = 200/factor
    size = round(size, 1-int(math.log10(size)))
    length = size*factor
    canvas.create_line(screenW-10,10, screenW-10-length,10, fill="white")
    canvas.create_line(screenW-10,5, screenW-10,15, fill="white")
    canvas.create_line(screenW-10-length,5, screenW-10-length,15, fill="white")
    canvas.create_text(screenW-10,20,anchor=NE,fill="white",text=str(size)+"AU")


#Main Program Loop:
###################
t = -50
def simulate():  
    global t
    canvas.delete(ALL)
    canvas.create_oval((screenW/2)-10,(screenH/2)-10,(screenW/2)+10,(screenH/2)+10, fill="yellow")
    
    displayPlanet(Mercury, colour="#d5d2d1")
    displayPlanet(Venus,   colour="#9b870c")
    displayPlanet(Earth,   colour="#0000ff")
    displayPlanet(Mars,    colour="#b7410e")
    displayPlanet(Jupiter, colour="#ffa500")
    displayPlanet(Saturn,  colour="#32cd32")
    displayPlanet(Uranus,  colour="#00ffff")
    displayPlanet(Neptune, colour="#003366")
    displayPlanet(Pluto,   colour="#ffc0cb")

    displayScale()

    if t > -19.99:
        date = epoch + datetime.timedelta(days=(t*100*365.25636))
        canvas.create_text(10,10,anchor=NW,fill="white",text=str(date))
    else:
        canvas.create_text(10,10,anchor=NW,fill="white",text=str(abs(round((t+20)*100)))+"BC")
    
    t += 0.01
    root.after(1, simulate)

simulate()
mainloop()
