import math

screenW, screenH = 800, 600
px_per_m = 100





class Polygon:
    def __init__(self, points, position, velocity=[0, 0], color="#ff0000", scale=1.0, mass=1.0, e=0.5, f=0.5):
        self.f = f
        self.e = e
        self.coefficient = e**2
        self.color = color
        self.points = [[j*scale for j in i] for i in points]
        self.n = len(points)
        self.position = position
        self.velocity = velocity
        self.angle = 0.0
        self.angvel = 0.0
        self.m = mass

        A = 0.5 * sum([(self.points[i][0]*self.points[(i+1)%self.n][1] - self.points[i][1]*self.points[(i+1)%self.n][0]) for i in range(self.n)])
        C_x = (1/(6*A)) * sum([(self.points[i][0]+self.points[(i+1)%self.n][0]) * (self.points[i][0]*self.points[(i+1)%self.n][1] - self.points[i][1]*self.points[(i+1)%self.n][0]) for i in range(self.n)])
        C_y = (1/(6*A)) * sum([(self.points[i][1]+self.points[(i+1)%self.n][1]) * (self.points[i][0]*self.points[(i+1)%self.n][1] - self.points[i][1]*self.points[(i+1)%self.n][0]) for i in range(self.n)])
        COM = [C_x, C_y]
        self.points = [[i[j]-COM[j] for j in range(2)] for i in self.points]

        numerator = 0.0
        denominator = 0.0
        for n in range(self.n):
            A, B = self.points[n], self.points[(n+1)%self.n]
            cross = A[0]*B[1]-B[0]*A[1]
            AA = A[0]**2 + A[1]**2
            BB = B[0]**2 + B[1]**2
            AB = A[0]*B[0] + A[1]*B[1]
            numerator += cross*(AA+AB+BB)
            denominator += cross
        self.I = (self.m/6)*(numerator/denominator)

    def vertices(self):
        c = math.cos(self.angle)
        s = math.sin(self.angle)
        return [[c*i[0]+s*i[1], -s*i[0]+c*i[1]] for i in self.points]#rotation
    
    def globalPos(self):
        return [[i[j]+self.position[j] for j in range(2)] for i in self.vertices()]

    def draw(self, canvas):
        verts = [[i[0]*px_per_m, screenH-i[1]*px_per_m] for i in self.globalPos()]
        canvas.create_polygon(verts, fill=self.color)

        max_y = max([i[1] for i in verts])
        if max_y+0.1 < screenH:
            return
        for i in verts:
            if abs(i[1] - max_y) < 0.1:
                canvas.create_oval(i[0]-5, i[1]-5, i[0]+5, i[1]+5, fill="#000000")

##    def physics(self, t, steps=1, g=9.81):
##        dt = t/steps
##        for i in range(steps):
##            self.physics_step(dt, g=g)
    def physics_step(self, t, g=9.81):
        old_state = self.position, self.velocity, self.angle
        
        self.position[0] += self.velocity[0]*t
        self.position[1] += self.velocity[1]*t - 0.5*g*(t**2)

        self.velocity[1] -= g*t

        self.angle += self.angvel*t

        new_gpos = self.globalPos()
        minmax = [min, max]
        screen = [screenW, screenH]
        for j in range(2):
            outside_screen = [lambda x: x < 0, lambda x: x > screen[j]/px_per_m]
            for p in range(2):
                for k in range(self.n):
                    if outside_screen[p](new_gpos[k][j]):
                        return False, old_state, new_gpos
        return True, old_state, new_gpos
                        

    def physics_collision(self, gpos, g=9.81):        
        m = self.m
        I = self.I

        minmax = [min, max]
        screen = [screenW, screenH]
        for j in range(2):
            outside_screen = [lambda x: x < 0, lambda x: x > screen[j]/px_per_m]
            for p in range(2):
                collision_points = []
                minmax_p = minmax[p]([i[j] for i in gpos])
                for k in range(self.n):
                    if outside_screen[p](gpos[k][j]) and abs(gpos[k][j] - minmax_p) < 0.001:
                        collision_points.append(k)

                if len(collision_points) > 0:
                    center = [sum([gpos[k][i] for k in collision_points])/len(collision_points) for i in range(2)]
                            
                    d = self.position[1-j] - center[1-j]
                    if j == 0:
                        d *= -1
                    w = self.angvel
                    v = self.velocity[j]
                    A = v+d*w
                    B = (1/m)+((d**2)/I)
                    impulse = (-A+((-1)**p)*math.sqrt(abs(A**2 - B*(1-self.coefficient)*(m*v**2+I*w**2))))/B# +/- root

                    self.velocity[j] += impulse/m
                    self.angvel += d * impulse/I

                    #Friction
                    perp_d = self.position[j] - center[j]
                    contact_v = self.velocity[1-j]+((-1)**(p+j))*self.angvel*math.sqrt(d**2+perp_d**2)
                    delta_v = self.f * abs(impulse/m)
                    if abs(contact_v) < delta_v:
                        delta_v = -contact_v
                    elif contact_v > 0:
                        delta_v *= -1
                    self.velocity[1-j] += delta_v
                    friction_impulse = delta_v * m
                    self.angvel += ((-1)**j) * perp_d * friction_impulse/I                        

                    local_center_j = sum([self.vertices()[k][j] for k in collision_points])/len(collision_points)
                    self.position[j] = p*(screen[j]/px_per_m)-local_center_j#+0.001


        

class RegularPolygon(Polygon):
    def __init__(self, n, r, position, velocity=[0, 0], color="#ff0000", mass=1.0, e=0.5, f=0.2):
        self.f = f
        self.e = e
        self.coefficient = e**2
        self.color = color
        self.n = n
        self.r = r
        self.position = position
        self.velocity = velocity
        self.angle = 0.0
        self.angvel = 0.0
        self.m = mass
        self.I = 0.5*mass*r**2 * (1-(2/3)*math.sin(math.pi/n)**2)
##        self.edge_length = 2*r*math.sin(math.pi/n)

    def vertices(self):
        return [[self.r*math.sin(-self.angle+(2*i+1)*math.pi/self.n), -self.r*math.cos(-self.angle+(2*i+1)*math.pi/self.n)] for i in range(self.n)]



class World:
    def __init__(self, polygons):
        self.polygons = polygons

    def physics(self, t, steps=1, g=9.81):
        dt = t/steps
        for i in range(steps):
            elapsed = 0
##            while elapsed < 1:
            elapsed += self.physics_step((1-elapsed)*dt, g=g)
            elapsed += self.physics_step((1-elapsed)*dt, g=g)
                
                

    def physics_step(self, t, g=9.81):
        old_poly_states = []
        new_poly_gpos = []
        for i in self.polygons:
            success, old_state, new_gpos = i.physics_step(t, g=g)
            old_poly_states.append(old_state)
            new_poly_gpos.append(new_gpos)
            if not success:
                break
        if success:
            return 1.0
        #else:
        for i in range(len(old_poly_states)):
            position, velocity, angle = old_poly_states[i]
            self.polygons[i].position = position
            self.polygons[i].velocity = velocity
            self.polygons[i].angle = angle

        interval = [0.0, 1.0]
        for j in range(4):
            midpoint = sum(interval)/2
            old_poly_states = []
            new_poly_gpos = []
            for i in self.polygons:
                success, old_state, new_gpos = i.physics_step((midpoint-interval[0])*t, g=g)
                old_poly_states.append(old_state)
                new_poly_gpos.append(new_gpos)
                if not success:
                    break
            if success:
                interval[0] = midpoint
                pass
            else:
                interval[1] = midpoint
                for i in range(len(old_poly_states)):
                    position, velocity, angle = old_poly_states[i]
                    self.polygons[i].position = position
                    self.polygons[i].velocity = velocity
                    self.polygons[i].angle = angle
                    
##        success, old_state, new_gpos = i.physics_step((interval[1]-interval[0])*t, g=g)
        new_poly_gpos = []
        for i in self.polygons:
            success, old_state, new_gpos = i.physics_step((interval[1]-interval[0])*t, g=g)
            new_poly_gpos.append(new_gpos)
        self.physics_collision(new_poly_gpos, g=g)
        return interval[1]
##        for i in self.polygons:
##            success, old_state, new_gpos = i.physics_step((1-interval[1])*t, g=g)

    def physics_collision(self, gpos, g=9.81):
        for i in range(len(self.polygons)):
            self.polygons[i].physics_collision(gpos[i], g=9.81)

    def draw(self, canvas):
        [i.draw(canvas) for i in self.polygons]




square = RegularPolygon(4, 1.41, [4.0, 3.0], [6.0, 2.0])

shapes = []
shapes = [RegularPolygon(i%4+3, 5.0/(i%4+8), [4.0, 4.0], [-4.0+i, 2.0+i], color="#"+hex(int(i%2)*200)[2:].zfill(2)+hex(int(i/4)*200)[2:].zfill(2)+hex(int((i%4)/2)*200)[2:].zfill(2)) for i in range(8)]
##shapes = [square]


shapes += [Polygon([[-8, -1], [8, -1], [8, 1], [-8, 1]], [4.0, 4.0], [8.0, 7.0], scale=0.125, e=0.5, f=0.4)]
##shapes += [Polygon([[-1, -1], [15, -1], [15, 1], [-1, 1]], [4.0, 4.0], [8.0, 7.0], scale=0.125, e=0.5, f=0.4)]

##shapes += [Polygon([[-2, -1], [2, -1], [2, 1], [1, 1], [1, 0], [-1, 0], [-1, 1], [-2, 1]], [4.0, 4.0], [3.0, 7.0], scale=0.25, e=0.5, f=0.4)]



from tkinter import *
root = Tk()
canvas = Canvas(root, width=screenW, height=screenH)
canvas.pack()

time_step = 0.01# seconds
time_step_ms = int(1000*time_step)

world = World(shapes)
def main():
    root.after(time_step_ms, main)
    canvas.delete(ALL)
    world.physics(time_step, steps=1)
    world.draw(canvas)


main()
mainloop()
