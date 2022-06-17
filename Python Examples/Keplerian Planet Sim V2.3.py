import math, time
from math import pi, sqrt


def sin(x):
    return math.sin(math.radians(x))
def cos(x):
    return math.cos(math.radians(x))
def mMProduct(M, N):
    if len(M[0]) != len(N):
        raise ValueError("Matrices product is undefined.", M, N)
    return [[sum([M[r][i]*N[i][c] for i in range(len(M[0]))]) for c in range(len(N[0]))] for r in range(len(M))]
def matrixNoteSwap(V): # Primarily for converting vectors from row to column notation, and vice versa. Performing this operation twice returns the original value.
    return [[c[r] for c in V] for r in range(len(V[0]))]

def getMatrixMinor(m,i,j):
    return [row[:j] + row[j+1:] for row in (m[:i]+m[i+1:])]
def getMatrixDeternminant(m):
    #base case for 2x2 matrix
    if len(m) == 2:
        return m[0][0]*m[1][1]-m[0][1]*m[1][0]
    determinant = 0
    for c in range(len(m)):
        determinant += ((-1)**c)*m[0][c]*getMatrixDeternminant(getMatrixMinor(m,0,c))
    return determinant
def getMatrixInverse(m):
    determinant = getMatrixDeternminant(m)
    #special case for 2x2 matrix:
    if len(m) == 2:
        return [[m[1][1]/determinant, -1*m[0][1]/determinant],
                [-1*m[1][0]/determinant, m[0][0]/determinant]]
    #find matrix of cofactors
    cofactors = []
    for r in range(len(m)):
        cofactorRow = []
        for c in range(len(m)):
            minor = getMatrixMinor(m,r,c)
            cofactorRow.append(((-1)**(r+c)) * getMatrixDeternminant(minor))
        cofactors.append(cofactorRow)
    cofactors = matrixNoteSwap(cofactors)
    for r in range(len(cofactors)):
        for c in range(len(cofactors)):
            cofactors[r][c] = cofactors[r][c]/determinant
    return cofactors

##############################################################################################################################################################################################################################

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

##############################################################################################################################################################################################################################
screenW, screenH = 1200, 800
fps = 60

from tkinter import *
root = Tk()
canvas = Canvas(root, width=screenW, height=screenH, bg="#000000")


move_state = {}
move_factor = {}
def movement_key_press(event):
    global move_state
    try:
        char = event.keysym
    except:
        char = event
    char = char.lower()
    move_state[char] = True
def movement_key_release(event):
    global move_state
    try:
        char = event.keysym
    except:
        char = event
    char = char.lower()
    move_state[char] = False

accelRate = 0.2
def move():
    if 'control_l' in move_state and move_state['control_l']:
        moveSpeed = 2
    else:
        moveSpeed = 5
    for char in move_state:
        if char not in move_factor:
            move_factor[char] = 0.0
        if move_state[char] == False:
            if move_factor[char] > 0:
                move_factor[char] -= accelRate
        elif move_state[char] == True:
            if move_factor[char] < 1:
                move_factor[char] += accelRate
        move_factor[char] = round(1000*move_factor[char])/1000
                
        if move_factor[char] > 0:
            viewVector = [0,0,0]
            if char == 'w':
                viewVector = [ 0, 0, 1]
            elif char == 's':
                viewVector = [ 0, 0,-1]
            elif char == 'a':
                viewVector = [-1, 0, 0]
            elif char == 'd':
                viewVector = [ 1, 0, 0]
            elif char == 'space':
                viewVector = [ 0, 1, 0]
            elif char == 'shift_l':
                viewVector = [ 0,-1, 0]
            worldVector = matrixNoteSwap(mMProduct(camera.invRotMatrix, matrixNoteSwap([viewVector])))[0]                
            camera.pos = tuple([camera.pos[i]+worldVector[i]*move_factor[char]*moveSpeed for i in range(3)])



mousePressed = False
lastPos = [0, 0]
def motion(event):
    global lastPos
    if mousePressed:
        x,y,z = camera.rot
        camera.rot = (x+((event.y - lastPos[1]) * 0.01), y+((event.x - lastPos[0]) * 0.01), z)
    lastPos = [event.x, event.y]
    

def mouseDown(event):
    global mousePressed, lastPos
    mousePressed = True
    lastPos = [event.x, event.y]

def mouseUp(event):
    global mousePressed
    mousePressed = False



canvas.bind("<KeyPress>", movement_key_press)
canvas.bind("<KeyRelease>", movement_key_release)
canvas.bind('<Motion>', motion)
canvas.bind('<ButtonPress-1>', mouseDown)
canvas.bind('<ButtonRelease-1>', mouseUp)
canvas.focus_set()

canvas.pack()



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



class Camera:
    def __init__(self, pos=(0,0,0), rot=(0,0,0), fov=120, f=None):
        self.pos = pos
        self.rot = rot
        self.fov = math.radians(fov)
        if f == None:
            self.f = 1/math.tan(0.5*self.fov)
        else:
            self.f = f
            
        self.oldRot = None
        self.generateMatrices()


    def generateMatrices(self):
        if self.oldRot != self.rot:
            self.oldRot = self.rot
            self.rotMatrix = rotXYZ([-j for j in list(self.rot)])
            self.invRotMatrix = getMatrixInverse(self.rotMatrix)

    def render(self, planets):
        self.rotMatrix = rotXYZ([-j for j in list(self.rot)])
        self.invRotMatrix = getMatrixInverse(self.rotMatrix)
        buffer = []
        buffer += self.renderObject((0,0,0), colour="#ffff00", r=10)
        for planet, colour, name in planets:
            buffer += self.renderObject(planet.compute(t), colour=colour)
            buffer += self.renderOrbit(planet, colour=colour)

        canvas.delete(ALL)
        buffer.sort(reverse=True, key=lambda x: x[1])
        for i in buffer:
            if i[0] == "line":
                canvas.create_line(i[2], fill=i[3], width=1.0)
            elif i[0] == "oval":
                canvas.create_oval(i[2], fill=i[3], outline="")

    def worldToView(self, coords, factor=100):
        coords = coords[0], -coords[2], coords[1]
        coords = tuple([coords[i]*factor for i in range(3)])
        coords = tuple([coords[i]-self.pos[i] for i in range(3)])
        coords = tuple(matrixNoteSwap(mMProduct(self.rotMatrix, matrixNoteSwap([list(coords)])))[0])
        return coords

    def viewToScreen(self, coords):
        x,y = (coords[0]/coords[2])*self.f, (coords[1]/coords[2])*self.f
        x = (x*screenW/2)+screenW/2
        y = screenH-((y*screenW/2)+screenH/2)
        return x,y

    def calcDot(self, A):
        a = self.worldToView(A)
        if a[2] < self.f:
            return None, None
        return self.viewToScreen(a), a[2]

    def calcLine(self, A, B):
        a, b = self.worldToView(A), self.worldToView(B)
        if min([a[2],b[2]]) < self.f:
            return None, None
        return (self.viewToScreen(a), self.viewToScreen(b)), (a[2]+b[2])/2

    def renderOrbit(self, planet, colour="#ffffff", step=10):
        RGB = [colour[1:3], colour[3:5], colour[5:7]]
        orbitColour = "#"+"".join([hex(min([255, int(i,base=16)+64]))[2:].zfill(2) for i in RGB])
        
        lines = []
        points = [planet.computeFromM(i) for i in range(0,360,step)]
        for i in range(len(points)):
            AB, z = self.calcLine(points[i], points[(i+1)%len(points)])
            if AB != None:
                lines.append(["line", z, AB, orbitColour])
        return lines

    def renderObject(self, coords, colour="#ffffff", r=5):
        A, z = self.calcDot(coords)
        if A != None:
            X,Y = A
            corners = (X-r, Y-r), (X+r, Y+r)
            return [["oval", z, corners, colour]]
        return []

##    def renderPlanet(self, planet, colour="#ffffff", r=5):
##        A, z = self.calcDot(planet.compute(t))
##        if A != None:
##            X,Y = A
##            corners = (X-r, Y-r), (X+r, Y+r)
##            return [["oval", z, corners, colour]]
##        return []
##
##    def renderSun(self, coords=(0,0,0), colour="#ffff00", r=10):
##        A, z = self.calcDot(coords)
##        if A != None:
##            X,Y = A
##            corners = (X-r, Y-r), (X+r, Y+r)
##            return [["oval", z, corners, colour]]
##        return []

##############################################################################################################################################################################################################################

hudSelection = None
def hud(planets, r=5):
    global hudSelection
    for planet, colour, name in planets:
        planetViewCoords = camera.worldToView(planet.compute(t))
        if planetViewCoords[2] > 0:
            planetPos = camera.viewToScreen(planetViewCoords)
            mousePos = lastPos[:]
            distance = sqrt((mousePos[0]-planetPos[0])**2 + (mousePos[1]-planetPos[1])**2)
            if mousePressed:
                if distance <= r:
                    hudSelection = name
            if name == hudSelection:
    ##            canvas.create_text(planetPos, text="test", fill="#ffffff", anchor=S, font=('Helvetica', 16))
                canvas.create_line(planetPos, (planetPos[0]+50,planetPos[1]-50), fill="#ffffff", width=2.0)
                canvas.create_line((planetPos[0]+50,planetPos[1]-50), (planetPos[0]+250,planetPos[1]-50), fill="#ffffff", width=2.0)
                canvas.create_text((planetPos[0]+50,planetPos[1]-50), fill="#ffffff", anchor=SW, font=('Helvetica', 16),
                                   text=str(name))
                
                height = sqrt(sum([i**2 for i in planet.compute(t)]))
                canvas.create_text((planetPos[0]+50,planetPos[1]-50), fill="#ffffff", anchor=NW, font=('Helvetica', 12),
                                   text=str(round(height,3)))
                canvas.create_text((planetPos[0]+100,planetPos[1]-50), fill="#ffffff", anchor=NW, font=('Helvetica', 12),
                                   text="AU")

                A, B = planet.compute(t), planet.compute(t+(0.01/31557600))
                velocity = [B[i]-A[i] for i in range(3)]
                speed = sqrt(sum([i**2 for i in velocity]))*1.496e+8
                canvas.create_text((planetPos[0]+150,planetPos[1]-50), fill="#ffffff", anchor=NW, font=('Helvetica', 12),
                                   text=str(round(speed,3)))
                canvas.create_text((planetPos[0]+200,planetPos[1]-50), fill="#ffffff", anchor=NW, font=('Helvetica', 12),
                                   text="km/s")

##############################################################################################################################################################################################################################

camera = Camera(pos=(0,0,-1),rot=(0,0,0))

planets = [(Mercury, "#d5d2d1", "Mercury"),
           (Venus,   "#9b870c", "Venus"  ),
           (Earth,   "#0000ff", "Earth"  ),
           (Mars,    "#b7410e", "Mars"   ),
           (Jupiter, "#ffa500", "Jupiter"),
           (Saturn,  "#32cd32", "Saturn" ),
           (Uranus,  "#00ffff", "Uranus" ),
           (Neptune, "#003366", "Neptune"),
           (Pluto,   "#ffc0cb", "Pluto"  )]

t = -50

frames = 0
lastcheck = 0
current_fps = 0
def main():
    frame_start = time.perf_counter()

    move()
    canvas.delete(ALL)
    camera.render(planets)
    hud(planets, r=10)
##    print(lastPos)

    global frames, lastcheck, current_fps
    frames += 1
    now = time.perf_counter()-lastcheck
    if now > 2:
        lastcheck += 1
        current_fps = frames + 0
        frames = 0
    canvas.create_text(20, 20, text=str(current_fps), font=10, fill="white")

    global t
    t += 0.0001
    root.after(max([1,int(1000*((1/fps) - time.perf_counter() + frame_start))]), main)


main()
mainloop()
