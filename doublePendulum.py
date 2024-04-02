Web VPython 3.2

#lengths of the strings
L1 = 1.5 #top string
L2 = 1 #middle string

#masses
M1 = 2 #top moving ball
M2 = 1 #middle moving ball

#g
g = 9.8

t = 0
dt = 0.001

#starting angles and velocities
theta1=80.5*pi/180 #release angle top ball
theta2= 0*pi/180 #release angle middle ball
theta1dot = 0
theta2dot = 0

pivot = sphere(pos=vector(0,L1,0), radius=0.05)

m1 = sphere(pos=pivot.pos-vector(0,L1,0),radius=0.05,color=color.white)
stick1 = cylinder(pos=pivot.pos,axis=m1.pos-pivot.pos,radius=0.015,color=color.yellow)

m2 = sphere(pos=m1.pos-vector(0,L2,0),radius=0.05,color=color.white)
stick2 = cylinder(pos=m1.pos, axis=m2.pos-m1.pos, radius=0.015,color=color.yellow)


m1.pos = pivot.pos+vector(L1*sin(theta1),-L1*cos(theta1),0)
m2.pos = m1.pos+vector(L2*sin(theta2),-L2*cos(theta2),0)

stick1.axis = m1.pos - pivot.pos
stick2.pos = m1.pos
stick2.axis = m2.pos - m1.pos


attach_trail(m2, retain=200, color=color.red)

eChart = gdisplay(x=500, y=0, width=600, height=400, title="", xtitle="", ytitle="", foreground=color.black, background=color.white)
ePlot = gdots(color=color.red)

while t < 50:  
    rate(1000)

    
    a = (M1 + M2)*L1**2
    b = M2*L1*L2
    c =- (M1 + M2)*g*L1*theta1
    d = M2*L1*L2
    k = M2*L2**2
    f =- M2*g*L2*theta2
    
    #angular acceleration
    theta1ddot = (-g*(2*M1 + M2)*sin(theta1) -M2*g*sin(theta1 -2*theta2)-2*sin(theta1 - theta2)*M2*(L2*theta2dot**2 + L1*cos(theta1 - theta2)*theta1dot**2))/(L1*(2*M1 + M2 - M2*cos(2*theta1 -2*theta2)))
    theta2ddot = (2*sin(theta1 - theta2)*((M1 + M2)*L1*theta1dot**2 + g*(M1 + M2)*cos(theta1) + L2*M2*cos(theta1 - theta2)*theta2dot**2))/(L2*(2*M1 + M2 - M2*cos(2*theta1 -2*theta2)))
    
    #angular velocity
    theta1dot = theta1dot + theta1ddot*dt
    theta2dot = theta2dot + theta2ddot*dt
    
    #angular position
    theta1 = theta1 + theta1dot*dt
    theta2 = theta2 + theta2dot*dt
    
    m1.pos = pivot.pos + vector(L1*sin(theta1),-L1*cos(theta1),0)
    m2.pos = m1.pos + vector(L2*sin(theta2),-L2*cos(theta2),0)
    stick1.axis = m1.pos - pivot.pos
    stick2.pos = m1.pos
    stick2.axis = m2.pos - m1.pos

    t = t+dt

    x1 = L1*sin(theta1)
    y1 = -L1*cos(theta1)
    x1dot = L1*theta1dot*cos(theta1)
    y1dot = L1*theta1dot*sin(theta1)
    x2 = x1 + L2*sin(theta2)
    y2 = y1 - L2*cos(theta2)
    x2dot = L1*theta1dot*cos(theta1)+L2*theta2dot*cos(theta2)
    y2dot = L1*theta1dot*sin(theta1)+L2*theta2dot*sin(theta2)
    T = 0.5*M1*(x1dot**2 + y1dot**2)+0.5*M2*(x2dot**2 + y2dot**2)
    U = M1*g*y1 + M2*g*y2
    
    E = T + U
    #ePlot.plot((t, E))
    
    
