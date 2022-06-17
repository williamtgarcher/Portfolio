################################################################################
# Version 1.11                                                                 #            
#                                                                              #
# The intention for this version is to add classes and rendering support for   #
# hyperprisms and/or hyperspheres. The hyperprisms will be able to generate    #
# from a given 3D shape, either by entering coordinates or by loading a .obj   #
# model file.                                                                  #
#                                                                              #
# Update Status = COMPLETE                               HYPERSPHERES:CANCELED #                          
################################################################################

import math, random

class HyperPlane:
    def __init__(self, n=[0,0,0,1], d=0):
        self.n = n
        self.d = d

    def reductionMatrix(self):
        normal = matrixNoteSwap([self.n])
        
        angleXW = math.atan2(normal[0][0], normal[3][0])
        rotMatXW = [
            [ math.cos(angleXW), 0, 0,-math.sin(angleXW)],
            [ 0                , 1, 0, 0                ],
            [ 0                , 0, 1, 0                ],
            [ math.sin(angleXW), 0, 0, math.cos(angleXW)]]
        normal = mMProduct(rotMatXW, normal)
        

        angleYW = math.atan2(normal[1][0], normal[3][0])
        rotMatYW = [
            [ 1,                 0, 0, 0                ],
            [ 0, math.cos(angleYW), 0,-math.sin(angleYW)],
            [ 0,                 0, 1, 0                ],
            [ 0, math.sin(angleYW), 0, math.cos(angleYW)]]
        normal = mMProduct(rotMatYW, normal)

        angleZW = math.atan2(normal[2][0], normal[3][0])
        rotMatZW = [
            [ 1, 0,                 0, 0                ],
            [ 0, 1,                 0, 0                ],
            [ 0, 0, math.cos(angleZW),-math.sin(angleZW)],
            [ 0, 0, math.sin(angleZW), math.cos(angleZW)]]
        normal = mMProduct(rotMatZW, normal)

        rotMat = mMProduct(mMProduct(rotMatZW, rotMatYW), rotMatXW)
        redMat = rotMat[:3]
        
        return redMat

class Polychoron:
    def __init__(self, pos=[0,0,0,0], scale=[1,1,1,1], rot=[0,0,0,0,0,0]):
        self.pos = pos
        self.scale = scale
##        self.rot = rot#XY,YZ,ZX, WX,WY,WZ
##        self.rotMap = [[1,1,0,0],[0,1,1,0],[1,0,1,0], [1,0,0,1],[0,1,0,1],[0,0,1,1]]

        self.rot = rot#YZ,XZ,XY, WX,WY,WZ
        self.rotMap = [[0,1,1,0],[1,0,1,0],[1,1,0,0], [1,0,0,1],[0,1,0,1],[0,0,1,1]]

        self.init()
        self.update()

    def update(self):
        self.verts = self.init_verts[:]
        #scale
        self.verts = [[i[j]*self.scale[j] for j in range(4)] for i in self.verts]
        #rotate
        R = [[1,0,0,0],[0,1,0,0],[0,0,1,0],[0,0,0,1]]
        for i in range(6):
            M = []
            for r in range(4):
                row = []
                for c in range(4):
                    if self.rotMap[i][r] == 1 and self.rotMap[i][c] == 1:
                        if self.rotMap[i][:r+1].count(1) == self.rotMap[i][:c+1].count(1):
                            element = math.cos(self.rot[i])
                        elif self.rotMap[i][:r+1].count(1) == 1 and self.rotMap[i][:c+1].count(1) == 2:
                            element = math.sin(self.rot[i])
                        elif self.rotMap[i][:r+1].count(1) == 2 and self.rotMap[i][:c+1].count(1) == 1:
                            element = -math.sin(self.rot[i])
                    elif r == c:
                        element = 1
                    else:
                        element = 0
                    row.append(element)
                M.append(row)
            R = mMProduct(M, R)
        self.verts = [matrixNoteSwap(mMProduct(R, matrixNoteSwap([i])))[0] for i in self.verts]
        #translate
        self.verts = [[i[j]+self.pos[j] for j in range(4)] for i in self.verts]

    def intersection(self, hPlane):
        self.update()
        ALLpoints = []
        sides = []
        for cell in self.cells:
            lines = []
            for face in cell:
                points = []
                for edge in self.faces[face]:
                    OA, OB = self.verts[self.edges[edge][0]], self.verts[self.edges[edge][1]]
                    AB = vectorOp(OB,'-',OA)
                    d = hPlane.d
                    n = hPlane.n

                    if vectorOp(AB,'.',n) != 0:
                        λ = ((d - vectorOp(OA,'.',n)) / vectorOp(AB,'.',n))
                        if 0 <= λ <= 1:
                            r = vectorOp(OA,'+',vectorOp(λ, '*', AB))
                            points.append(r)
                            ALLpoints.append(r)
                lines.append(points)
            sides.append(lines)
        ALLpointsTWO = []
        for i in ALLpoints:
            if i not in ALLpointsTWO:
                ALLpointsTWO.append(i)

        faces = []
        sides = [[edge for edge in side if edge != []] for side in sides]
        sides = [[edge for edge in side if len(edge)==2] for side in sides]
        for side in sides:
            face = []
            if len(side) != 0:
                face.append(side[0][0])
                face.append(side[0][1])
                loops = 0
                while len(face) < len(side):
                    loops += 1
                    for edge in side[1:]:
                        if edge[0] == face[-1]:
                            face.append(edge[1])
                        elif edge[1] == face[-1]:
                            face.append(edge[0])
                    if len(face) < len(side):
                        for edge in side[1:]:
                            if edge[0] == face[0]:
                                face.insert(0,edge[1])
                            elif edge[1] == face[0]:
                                face.insert(0,edge[0])
                    if loops > 100:
##                        raise StopIteration("While loop exceeded 100 loops (therefore it is stuck)\n    Data:\n"+str([len(face),len(side), face, side]))
                        break
            faces.append(face)
                        
        return faces

class HyperCube(Polychoron):
    def init(self):
        self.verts = []#in coords
        self.edges = []#in verts ids
        self.faces = []#in edges ids
        self.cells = []#in faces ids
        for x in range(2):
            for y in range(2):
                for z in range(2):
                    for w in range(2):
                        self.verts.append([x-0.5, y-0.5, z-0.5, w-0.5])
        for v in range(len(self.verts)):
            for i in range(4):
                points = [self.verts[v][:] for x in range(2)]
                for x in range(2):
                    points[x][i] *= 1-2*x
                edge = [self.verts.index(points[x]) for x in range(2)]
                edge.sort()
                if edge not in self.edges:
                    self.edges.append(edge)
        for v in range(len(self.verts)):
            for i in range(4):
                for j in range(i+1,4):
                    points = [[self.verts[v][:] for x in range(2)] for y in range(2)]
                    for x in range(2):
                        for y in range(2):
                            points[x][y][i] *= 1-2*x
                            points[x][y][j] *= 1-2*y
                    faceVerts = [self.verts.index(points[x][y]) for x in range(2) for y in range(2)]
                    faceEdges = []
                    for e in range(len(self.edges)):
                        if len([1 for q in self.edges[e] if q in faceVerts]) == 2:
                            faceEdges.append(e)
                    faceEdges.sort()
                    if faceEdges not in self.faces:
                        self.faces.append(faceEdges)
        for v in range(len(self.verts)):
            for i in range(4):
                for j in range(i+1,4):
                    for k in range(j+1,4):
                        points = [[[self.verts[v][:] for x in range(2)] for y in range(2)] for z in range(2)]
                        for x in range(2):
                            for y in range(2):
                                for z in range(2):
                                    points[x][y][z][i] *= 1-2*x
                                    points[x][y][z][j] *= 1-2*y
                                    points[x][y][z][k] *= 1-2*z
                        cellVerts = [self.verts.index(points[x][y][z]) for x in range(2) for y in range(2) for z in range(2)]
                        cellEdges = []
                        for e in range(len(self.edges)):
                            if len([1 for q in self.edges[e] if q in cellVerts]) == 2:
                                cellEdges.append(e)
                        cellFaces = []
                        for f in range(len(self.faces)):
                            if len([1 for q in self.faces[f] if q in cellEdges]) == 4:
                                cellFaces.append(f)
                        cellFaces.sort()
                        if cellFaces not in self.cells:
                            self.cells.append(cellFaces)
                            
        self.init_verts = self.verts[:]

        

    

class FiveCell(Polychoron):
    def init(self):
        self.edges = []#in verts ids
        self.faces = []#in edges ids
        self.cells = []#in faces ids

        q = math.sqrt(5)
        self.verts = [[1,1,1,-1/q],#in coords
                      [1,-1,-1,-1/q],
                      [-1,1,-1,-1/q],
                      [-1,-1,1,-1/q],
                      [0,0,0,q-1/q]]
        self.verts = [[j/(2*math.sqrt(2)) for j in i] for i in self.verts]
        for i in range(5):
            for j in range(i+1, 5):
                self.edges.append([i, j])
        for i in range(5):
            for j in range(i+1, 5):
                for k in range(j+1, 5):
                    if [i, k] in self.edges and [j, k] in self.edges:
                        self.faces.append([self.edges.index([i,k]),
                                           self.edges.index([j,k]),
                                           self.edges.index([i,j])])
##        [i.sort() for i in self.faces]
##        print(self.faces)
        for i in range(10):
            for j in range(i+1, 10):
                for k in range(j+1, 10):
                    for m in range(k+1, 10):
                        if sum([1 for f in [m,j,k] if [1 for e in range(3) if self.faces[i][e] in self.faces[f]]]) == 3:
                            if sum([1 for f in [i,m,k] if [1 for e in range(3) if self.faces[j][e] in self.faces[f]]]) == 3:
                                if sum([1 for f in [i,j,m] if [1 for e in range(3) if self.faces[k][e] in self.faces[f]]]) == 3:
                                    self.cells.append([i,j,k,m])

        self.init_verts = self.verts[:]

        


class SixteenCell(Polychoron):
    def init(self):
        self.edges = []#in verts ids
        self.faces = []#in edges ids
        self.cells = []#in faces ids

        self.verts = [[0.5,0,0,0], [-0.5,0,0,0],
                      [0,0.5,0,0], [0,-0.5,0,0],
                      [0,0,0.5,0], [0,0,-0.5,0],
                      [0,0,0,0.5], [0,0,0,-0.5]]
        for i in range(8):
            for j in range(i+1, 8):
                if [-k for k in self.verts[j]] != self.verts[i]:
                    self.edges.append([i,j])
        for i in range(24):
            for j in range(i+1, 24):
                for k in range(j+1, 24):
                    if [i, k] in self.edges and [j, k] in self.edges and [i, j] in self.edges:
                        self.faces.append([self.edges.index([i,k]),
                                           self.edges.index([j,k]),
                                           self.edges.index([i,j])])
        for i in range(32):
            for j in range(i+1, 32):
                for k in range(j+1, 32):
                    for m in range(k+1, 32):
##                        if sum([1 for f in [m,j,k] if [1 for e in range(3) if self.faces[i][e] in self.faces[f]]]) == 3:
##                            if sum([1 for f in [i,m,k] if [1 for e in range(3) if self.faces[j][e] in self.faces[f]]]) == 3:
##                                if sum([1 for f in [i,j,m] if [1 for e in range(3) if self.faces[k][e] in self.faces[f]]]) == 3:
##                                    self.cells.append([i,j,k,m])
                        if len(set([v for f in [i,j,k,m] for e in self.faces[f] for v in self.edges[e]])) == 4:
                            self.cells.append([i,j,k,m])

        self.init_verts = self.verts[:]




class HyperPrism(Polychoron):
    def __init__(self, modelName, length, pos=[0,0,0,0], scale=[1,1,1,1], rot=[0,0,0,0,0,0]):
        self.modelName = modelName
        self.length = length
        Polychoron.__init__(self, pos=pos, scale=scale, rot=rot)


    def init(self):
        data = [i.split() for i in open(self.modelName+".obj", "r").read().split("\n")]
        model_vertices = [[float(j) for j in i[1:4]] for i in data if len(i)>0 if i[0]=="v"]
        model_faces = [[int(j)-1 for j in i[1:4]] for i in data if len(i)>0 if i[0]=="f"]
        
        model_edges = []
        for face in model_faces:
            for i in range(len(face)):
                edge = [face[i], face[(i+1)%len(face)]]
                edge.sort()
                if edge not in model_edges:
                    model_edges.append(edge)
        
        self.verts = []#in coords
        self.edges = []#in verts ids
        self.faces = []#in edges ids
        self.cells = []#in faces ids

        self.verts = [i+[-self.length/2] for i in model_vertices]+[i+[self.length/2] for i in model_vertices]
        self.edges = model_edges + [[j+len(model_vertices) for j in i] for i in model_edges] + [[i,i+len(model_vertices)] for i in range(len(model_vertices))]
        
        faces_in_vert_ids = model_faces + [[j+len(model_vertices) for j in i] for i in model_faces] + [[i[0], i[1], i[1]+len(model_vertices), i[0]+len(model_vertices)] for i in model_edges]
        for face in faces_in_vert_ids:
            self.faces.append([])
            for i in range(len(face)):
                edge_in_verts = [face[i], face[(i+1)%len(face)]]
                edge_in_verts.sort()
                self.faces[-1].append(self.edges.index(edge_in_verts))

        self.cells = [[i for i in range(len(model_faces))]] + [[i+len(model_faces) for i in range(len(model_faces))]] + [[i, i+len(model_faces)]+[j+2*len(model_faces) for j in self.faces[i]] for i in range(len(model_faces))]#cells between end faces
                            
        self.init_verts = self.verts[:]

        
                        

def vectorOp(a, op, b):
    if op == ".":#vector a . vector b
        return sum([a[i]*b[i] for i in range(len(a))])
    if op == "+":#vector a + vector b
        return [a[i]+b[i] for i in range(len(a))]
    if op == "-":#vector a - vector b
        return [a[i]-b[i] for i in range(len(a))]
    if op == "*":#scalar a * vector b
        return [a*b[i] for i in range(len(b))]

def mMProduct(M, N):
    if len(M[0]) != len(N):
        raise ValueError("Matrices product is undefined.", M, N)
    return [[sum([M[r][i]*N[i][c] for i in range(len(M[0]))]) for c in range(len(N[0]))] for r in range(len(M))]

def mSProduct(M, S):
    return [[j*S for j in i] for i in M]


def matrixNoteSwap(V): # Primarily for converting vectors from row to column notation, and vice versa. Performing this operation twice returns the original value.
    return [[c[r] for c in V] for r in range(len(V[0]))]


##############################################################################################################################################################################################################################
screenW, screenH = 1200, 800
fps = 60

from tkinter import *
root = Tk()
canvas = Canvas(root, width=screenW, height=screenH, bg="#EEEEEE")


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
            x,y,z = camera.pos
            if char == 'w':
                x += moveSpeed * math.sin(camera.rot[1]) * move_factor[char]
                z += moveSpeed * math.cos(camera.rot[1]) * move_factor[char]
            elif char == 's':
                x -= moveSpeed * math.sin(camera.rot[1]) * move_factor[char]
                z -= moveSpeed * math.cos(camera.rot[1]) * move_factor[char]
            elif char == 'a':
                x -= moveSpeed * math.cos(camera.rot[1]) * move_factor[char]
                z += moveSpeed * math.sin(camera.rot[1]) * move_factor[char]
            elif char == 'd':
                x += moveSpeed * math.cos(camera.rot[1]) * move_factor[char]
                z -= moveSpeed * math.sin(camera.rot[1]) * move_factor[char]
            elif char == 'space':
                y += moveSpeed * move_factor[char]
            elif char == 'shift_l':
                y -= moveSpeed * move_factor[char]
            camera.pos = (x,y,z)
            
##            x,y,z = camera.rot
##            rotSpeed = math.pi/128
##            if char == 'q':
##                y -= rotSpeed
##            elif char == 'e':
##                y += rotSpeed
##            elif char == 'r':
##                x -= rotSpeed
##            elif char == 'f':
##                x += rotSpeed
##            elif char == 't':
##                z += rotSpeed
##            elif char == 'g':
##                z -= rotSpeed
##            camera.rot = (x,y,z)


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


class Dial:
    def __init__(self, pos, attribute, increments=0.1256, r=20, colour="#FFFFFF", value=0, maxValue=6.28,
                 labelText="", labelFont="10", labelFill=""):
        self.pos = pos
        self.attribute = attribute
        self.inc = increments
        self.r = r
        self.colour = colour
        self.value = value
        self.initValue = value
        self.maxValue = maxValue
        self.updateAttribute()
        self.labelText = labelText
        self.labelFont = labelFont
        self.labelFill = labelFill
    def draw(self):
        pos, r = self.pos, self.r
        canvas.create_oval(pos[0]-r, pos[1]-r, pos[0]+r, pos[1]+r, fill=self.colour, width=3.0)
        if self.maxValue != None:
            canvas.create_line(pos[0],pos[1],
                               pos[0]+r*math.sin((self.value/self.maxValue)*2*math.pi),
                               pos[1]-r*math.cos((self.value/self.maxValue)*2*math.pi), width=3.0)
        if len(self.labelText)>0 and self.labelFill != "":
            canvas.create_text(pos[0], pos[1], fill=self.labelFill, font=self.labelFont, text=self.labelText)
    def scroll(self, mouse, direction):
        if (mouse[0]-self.pos[0])**2 + (mouse[1]-self.pos[1])**2 < self.r**2:
            self.value += self.inc * direction
            if self.maxValue != None:
                self.value %= self.maxValue
            self.updateAttribute()
            return 1
    def reset(self, mouse):
        if (mouse[0]-self.pos[0])**2 + (mouse[1]-self.pos[1])**2 < self.r**2:
            self.value = self.initValue + 0
            self.updateAttribute()
            return 1
    def updateAttribute(self):
        global selectionID
        if selectionID == len(objects):
            for obj in objects:
                exec("obj."+self.attribute+" = self.value")
        else:
            exec("objects[selectionID]."+self.attribute+" = self.value")
            

def scroll(event):
    direction = int(event.delta/abs(event.delta))
    if 1 not in [i.scroll([event.x, event.y], direction) for i in dials]:
        hPlane.d += direction*10

def right_click(event):
    if 1 not in [i.reset([event.x, event.y]) for i in dials]:
        hPlane.d = 0

def tab(event):
    global selectionID, objects
    selectionID += 1
    selectionID %= len(objects)+1

canvas.bind("<KeyPress>", movement_key_press)
canvas.bind("<KeyRelease>", movement_key_release)
canvas.bind('<Motion>', motion)
canvas.bind('<ButtonPress-1>', mouseDown)
canvas.bind('<ButtonRelease-1>', mouseUp)
canvas.bind("<MouseWheel>", scroll)
canvas.bind('<ButtonPress-3>', right_click)
canvas.bind("<Tab>", tab)
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


    def render(self, objects, plane, colourMap=["red", "orange", "yellow", "green", "blue", "pink", "purple", "grey"]):
        self.rot = (min([max([self.rot[0], -1.57]), 1.57]), self.rot[1], self.rot[2])
        rotMatrix = rotXYZ([-j for j in list(self.rot)])
        canvas.delete(ALL)

        allFaces = [i.intersection(plane) for i in objects]
        redMat = plane.reductionMatrix()
        allFaces = [[[matrixNoteSwap(mMProduct(redMat, matrixNoteSwap([j])))[0] for j in i] for i in faces] for faces in allFaces]

        
        buffer = []
        for faces in range(len(allFaces)):
##            print(len(allFaces[faces]))
            for f in range(len(allFaces[faces])):
                face = allFaces[faces][f]
                verts = [(x-self.pos[0],y-self.pos[1],z-self.pos[2]) for x,y,z in face]
                verts = [matrixNoteSwap(mMProduct(rotMatrix, matrixNoteSwap([list(i)])))[0] for i in verts]
                points = []
                depths = []
                positions = []
                for i in range(len(face)):
                    a, b = verts[i], verts[(i+1)%len(face)]
                    if max([a[2],b[2]]) >= self.f:
                        points += self.renderLine(a, b)
##                        depths.append(math.sqrt(sum([a[k]**2 for k in range(3)])))
                        positions.append(a)
                if len(points) >= 3:
##                    z = sum(depths)/len(depths)
                    z = [sum([j[i] for j in positions])/len(positions) for i in range(3)]
                    z = z[0]**2 + z[1]**2 + z[2]**2
                    buffer.append([points, z, colourMap[f%len(colourMap)], faces])
        buffer.sort(reverse=True, key=lambda x: x[1])
        for i in buffer:
            if selectionID == i[3] or selectionID == len(objects):
                canvas.create_polygon(i[0], fill=i[2], outline="#000000")
            else:
                canvas.create_polygon(i[0], fill=i[2], outline="")
            

    def renderLine(self, a, b):
        if min([a[2],b[2]]) < self.f:
            if a[2] > b[2]:
                BA = [a[i]-b[i] for i in range(3)]
                p = (self.f - b[2])/BA[2]
                BA = [i*p for i in BA]
                b = [b[i]+BA[i] for i in range(3)]
            else:
                AB = [b[i]-a[i] for i in range(3)]
                p = (self.f - a[2])/AB[2]
                AB = [i*p for i in AB]
                a = [a[i]+AB[i] for i in range(3)]

        xA,yA = (a[0]/a[2])*self.f, (a[1]/a[2])*self.f
        xA = (xA*screenW/2)+screenW/2
        yA = screenH-((yA*screenW/2)+screenH/2)

        xB,yB = (b[0]/b[2])*self.f, (b[1]/b[2])*self.f
        xB = (xB*screenW/2)+screenW/2
        yB = screenH-((yB*screenW/2)+screenH/2)

        return [[xA,yA], [xB,yB]]






hPlane = HyperPlane(n=[0,0,0,1],d=0)

hCube = HyperCube(pos=[0,0,0,0],scale=[100,100,100,100], rot = [0,0,0,0,0,0])
fcell = FiveCell(pos=[-150,0,0,0],scale=[100,100,100,100], rot = [0,0,0,0,0,0])
scell = SixteenCell(pos=[0,150,0,10],scale=[100,100,100,100], rot = [0,0,0,0,0,0])

objects = [hCube, fcell, scell]
selectionID = 0

dials = [Dial([1100,650], "rot[0]", labelText="YZ", labelFill="#FF0000"),
         Dial([1100,700], "rot[1]", labelText="XZ", labelFill="#00FF00"),
         Dial([1100,750], "rot[2]", labelText="XY", labelFill="#0000FF"),
         Dial([1150,650], "rot[3]", labelText="WX", labelFill="#FF7777"),
         Dial([1150,700], "rot[4]", labelText="WY", labelFill="#77FF77"),
         Dial([1150,750], "rot[5]", labelText="WZ", labelFill="#7777FF"),

         Dial([1150,400], "pos[0]", increments=10, r=15, colour="#FF0000", value=0, maxValue=None, labelText="X", labelFill="#FFFFFF"),
         Dial([1150,450], "pos[1]", increments=10, r=15, colour="#00FF00", value=0, maxValue=None, labelText="Y", labelFill="#FFFFFF"),
         Dial([1150,500], "pos[2]", increments=10, r=15, colour="#0000FF", value=0, maxValue=None, labelText="Z", labelFill="#FFFFFF"),
         Dial([1150,550], "pos[3]", increments=10, r=15, colour="#333333", value=0, maxValue=None, labelText="W", labelFill="#FFFFFF")
         ]

camera = Camera(pos=(0,0,-1),rot=(0,0,0))

##random.seed(1)
##colourmap = ["#"+"".join(["0123456789ABCDEF"[random.randint(0,15)] for i in range(6)]) for i in range(16)]

def mapping(angle):#in degrees
    angle %= 360
    if angle > 180:
        angle = 360 - angle
    if 0 <= angle < 60:
        return 1
    elif 120 < angle <= 180:
        return 0
    else:
        return (1 - (angle-60)/60)
angles = [i*220 for i in range(18)]
colourmap = ["#"+"".join([hex(256+abs(int(255 * mapping(angles[j]-(120*i)))))[-2:].upper() for i in range(3)]) for j in range(18)]


import time
frames = 0
lastcheck = 0
current_fps = 0


q = 0
def main():
    frame_start = time.perf_counter()

    move()
    camera.render(objects, hPlane, colourMap=colourmap)
    [i.draw() for i in dials]

    global frames, lastcheck, current_fps
    frames += 1
    now = time.perf_counter()-lastcheck
    if now > 2:
        lastcheck += 1
##        print(frames)
        current_fps = frames + 0
        frames = 0
    canvas.create_text(20, 20, text=str(current_fps), font=10)
    
    root.after(max([1,int(1000*((1/fps) - time.perf_counter() + frame_start))]), main)
##    print(int(1000*((1/fps) - time.perf_counter() + frame_start)))


main()
mainloop()
