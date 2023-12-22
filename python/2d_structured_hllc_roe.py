import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import sys#,os
# sys.path.append('../')
from exact_riemann import wetbed
from scipy.optimize import fsolve
from scipy.integrate import solve_ivp


C=0.5
n = 0.0 #manning friction
tend = 60.0
g = 9.81;
htol = 1.0e-2 #height tolerance below which hl&hr ; ul&ur is set to zero
lx1 = -2.5; lx2 = 2.5;
ly1 = -2.5; ly2 = 2.5;
Nx = 100; #cells x-axis
Ny = 100; #cells y-axis

# bslope = ['flatbed','slope_upwards','slope_downwards','obstruction',
#           'bump','flatbed+sloping_stage','2d_parabolic_dam','riemann_problem',
#           '2d_dambreak','circular_dambreak','bump_steady_flow','bump_partial_cover']

# choose the problem you need from the list of bslope below by setting ibslope
# to the correct index
bslope = ['bump','2d_parabolic_dam','riemann_problem',
          '2d_dambreak','circular_dambreak','bump_steady_flow','bump_partial_cover']
ibslope = 4

#choose the bump configuration for the bump problem mentioned in swashes
bump_config = ['noflow','emerged','subcritical']
ibump = 1

bump_config = bump_config[ibump]
bslope = bslope[ibslope]
# choose the boundary condition needed
bc_opt=['flow','noflow']
bc = bc_opt[0]
# choose solver
solver = ['hllc','roe']
solver = solver[1]

#setting up of initial problem definitions
if bslope =='2d_parabolic_dam' or bslope =='bump':
    lx1,ly1 = 0.,0.
    lx2 = 25.
    ly2 = 10.
    if bslope =='2d_parabolic_dam':
        Nx = 10;
        Ny = 10;
    elif bslope =='bump':
        
        if bump_config == 'emerged':
            C=0.01
            Nx = 200; 
            Ny = 200; 
        elif bump_config == 'noflow':
            Nx = 10; 
            Ny = 10; 
        elif bump_config == 'subcritical':
            Nx = 50; 
            Ny = 50; 
        
elif bslope =='riemann_problem':
    
    lx1,ly1 = -2.5,-2.5
    lx2 = 2.5
    ly2 = 2.5
    Nx = 100; 
    Ny = 100; 
elif bslope =='2d_dambreak':
    C = 0.1
    lx1,ly1 = -1.0,-1.0#2.5#10
    lx2,ly2 = 1.0,1.0#2.5#5
    Nx = 40; 
    Ny = 40; 

elif bslope =='circular_dambreak':
    lx1 = -2.5; lx2 = 2.5;
    ly1 = -2.5; ly2 = 2.5;
    Nx = 30; 
    Ny = 30; 
elif bslope in ['bump_steady_flow','bump_partial_cover']:
    C = 0.1
    lx1,ly1 = -0.5,-0.5#2.5#10
    lx2,ly2 = 0.5,0.5#2.5#5
    Nx = 10; 
    Ny = 10; 
    if bslope == 'bump_partial_cover':
        Nx = 30;
        Ny = 30;

Lx = lx2 - lx1;
Ly = ly2 - ly1;
dx = Lx/Nx; 
dy = Ly/Ny; 

xx = np.linspace(lx1-1.5*dx,lx2+1.5*dx,Nx+4);
yy = np.linspace(ly1-1.5*dy,ly2+1.5*dy,Ny+4);
x = np.arange(lx1-2.0*dx,lx2+2.01*dx,dx); #x coordinate of cell center
y = np.arange(ly1-2.0*dy,ly2+2.01*dy,dy); #y coordinate of cell center


#Preallocation
t = 0; #start time
h = np.zeros((Nx+4,Ny+4),float);
hlx = np.zeros((Nx+4,Ny+4),float);
hrx = np.zeros((Nx+4,Ny+4),float);
hly = np.zeros((Nx+4,Ny+4),float);
hry = np.zeros((Nx+4,Ny+4),float);
H = np.zeros((Nx+4,Ny+4),float); #total head
u = np.zeros((Nx+4,Ny+4),float);
v = np.zeros((Nx+4,Ny+4),float);
hu = np.multiply(h,u);  #discharge along x
hv = np.multiply(h,v);  #ischarge along y
f1 = np.zeros((Nx+3,Ny+2),float); #fluxes in x
f2 = np.zeros((Nx+3,Ny+2),float); #fluxes in x
f3 = np.zeros((Nx+3,Ny+2),float); #fluxes in x
g1 = np.zeros((Nx+2,Ny+3),float); #fluxes in y
g2 = np.zeros((Nx+2,Ny+3),float); #fluxes in y
g3 = np.zeros((Nx+2,Ny+3),float); #fluxes in y
z = np.zeros((Nx+4,Ny+4),float);
zbx = np.zeros((Nx+4,Ny+4),float);
zby = np.zeros((Nx+4,Ny+4),float);

# Defining classes for problem definitions
class initial_conditions():
    def __init__(self,Nx,Ny,h=0.0,q=0.0):
        self.Nx = Nx+4
        self.Ny = Ny+4
        self.h = np.full((Ny+4,Nx+4),h)
        # self.u = np.full((N,1),u)
        self.q = np.full((Ny+4,Nx+4),q)

class spatial_domain():
    def __init__(self,lx1,ly1,lx2,ly2,Nx,Ny,n=0):
        #meeting 
        # add S0
        self.lx1 = lx1 #length
        self.ly1 = ly1 #length
        self.Nx = Nx #number of nodes
        self.Ny = Ny #number of nodes
        dx = (lx2-lx1)/Nx 
        dy = (ly2-ly1)/Ny 
        self.dx = dx
        self.dy = dy
        self.n = n #Manning's roughness
        # self.x = np.linspace(-3*dx/2,l+3*dx/2,N)
        self.x = np.linspace(lx1-1.5*dx,lx2+1.5*dx,Nx+4);
        self.y = np.linspace(ly1-1.5*dy,ly2+1.5*dy,Ny+4);
        self.x1 = np.arange(lx1-2.0*dx,lx2+2.01*dx,dx); #x coordinate of cell center
        self.y1 = np.arange(ly1-2.0*dy,ly2+2.01*dy,dy); #y coordinate of cell center

        self.zbed = np.zeros((Nx+4,Ny+4),float) #channel bed
        self.S0x = np.zeros((Nx+4,Ny+4),float)
        self.S0y = np.zeros((Nx+4,Ny+4),float)

#setting up the benchmark problems
D = spatial_domain(lx1=lx1,ly1=ly1,lx2=lx2,ly2=ly2,Nx=Nx,Ny=Ny,n=n)
IC = initial_conditions(Nx=Nx,Ny=Ny)


if bslope == 'flat_bed' or bslope == 'obstruction':
    for i in range(Nx+4):
        for j in range(Ny+4):
            D.zbed[i,j] = 0.0
elif bslope == 'slope_upwards':
    # bed is sloping upwards
    for i in range(Nx+4):
        for j in range(Ny+4):
            #row - y, col -x 
            D.zbed[i,j] = 0.0*D.y[j] + 0.01*D.x[i] 
elif bslope == 'slope_downwards':
    hl1 = 4.0
    # bed is sloping downwards
    for i in range(Nx+4):
        for j in range(Ny+4):
            #row - y, col -x 
            D.zbed[i,j] = 0.0*D.y[j] - 0.1*D.x[i]  + 3.25
elif bslope == 'bump':
    
    if bump_config == 'subcritical':
        C = 0.1
        N = 500
        hl1 =2.0
        hr1 = 2.0
        ul1=0.0
        ur1=0.0
        qb = 4.42
        hb = 2.0
        tend = 300.0
    elif bump_config == 'emerged':
        C = 0.5
        N = 2000 #number of cells
        hl1 = 0.1
        hr1 = 0.1
        ul1=0.0
        ur1=0.0
        qb = 0.0
        hb = 0.1
        tend = 60.0
    elif bump_config == 'noflow':
        C = 0.5
        N = 1000 #number of cells
        hl1 = 0.5
        hr1 = 0.5
        ul1=0.0
        ur1=0.0
        qb = 0.0
        hb = 0.5
        tend = 60.0
    for i in range(D.Nx+4):
        for j in range(D.Ny+4):
            if D.x[i] > 8 and D.x[i] < 12:
                D.zbed[i,j] = 0.2 - 0.05*(D.x[i]-10)**2
                D.S0x[i,j] = -0.1*(D.x[i]-10)
                D.S0y[i,j] = 0
        #         #S0 = dz/dx or -dx/dx
    # IC.h = 2 - D.zbed
elif bslope == 'flat_bed+sloping_stage':
    for i in range(D.N):
        D.zbed[i] = 0.0
        D.S0[i] = 0.0
elif bslope == '2d_parabolic_dam':
    dam_h = 0.5 #m
    dam_w = 1.0 #m
    dam_d = 10.0 #m
    D.zbed = np.zeros((Ny+4,Nx+4))
    for i in range(D.Nx+4):
        for j in range(D.Ny+4):
            
            phi_xy = D.x[i] - ((D.y[j] - 0.5*(Ly))/(Ly))**2 - dam_d
            alpha = (np.exp(-((0.5*dam_w)**2)) - np.exp(-(2.0*dam_w)))/dam_h
            beta = ((np.exp(-((0.5*dam_w)**2)))/alpha) - dam_h
            
            term2 = max(0,((np.exp(-(phi_xy**2)))/alpha)-beta)
            
            D.zbed[i,j] = min(0.5,term2)
            if phi_xy <0:
                IC.h[j,i] = 0.5
            else:
                IC.h[j,i] = 0.0
elif bslope == '2d_dambreak':
    for i in range(Nx+4):
        for j in range(Ny+4):
            #row - y, col -x 
            # if D.x[j,i]>-0.2 and D.x[j,i]<0.2 && (D.y[j,i]<-0.2 || D.y[j,i]>0.2)
            # if D.x[i]>(lx-1.2) and D.x[i]<(lx-0.8) and ((D.y[j]<(Ly-1.2) or D.y[j]>(Ly-0.8))):
            if D.x[i]>(-0.2) and D.x[i]<(0.2) and ((D.y[j]<(-0.2) or D.y[j]>(0.2))):
                D.zbed[i,j] = 1.0
elif bslope == 'bump_steady_flow':
    hl1 = 0.5
    ul1=0.0
    for i in range(Nx+4):
        for j in range(Ny+4):
            if abs(D.y[j]) < 0.2 and abs(D.x[i]) < 0.2:
                D.zbed[i,j] = max(0,0.25 - 5*( (D.x[i]-0.0)**2 + (D.y[j]-0.0)**2))
elif bslope == 'bump_partial_cover':
    hl1 = 0.1
    ul1= 0.0
    for i in range(Nx+4):
        for j in range(Ny+4):
            if abs(D.y[j]) < 0.2 and abs(D.x[i]) < 0.2:
                D.zbed[i,j] = max(0,0.25 - 5*( (D.x[i]-0.0)**2 + (D.y[j]-0.0)**2))
X, Y = np.meshgrid(xx, yy) 
z = (D.zbed)
for i in range(0,Nx+4): # x<=0
    for k in range(0,Ny+4):
        r = np.sqrt((D.x[i]-0)**2+(D.y[k]-0)**2);
        if bslope == 'circular_dambreak':
            hl1 =2.0
            hr1 = 1.0
            ul1 = 0.0
            ur1 = 0.0
            # xdam =0.0
            if (r <= 0.5):
                h[i,k] = hl1 - z[i,k];
                u[i,k] = ul1;
            else:
                h[i,k] = hr1;
                u[i,k] = ur1;
        elif bslope == '2d_parabolic_dam':
            hl1 =1.0
            hr1 = 0.0
            ul1 = 0.0
            ur1 = 0.0
            # xdam =10.0
            tend = 40.0
            phi_xy = D.x[i] - ((D.y[k] - 0.5*(Ly))/(Ly))**2 - dam_d
            if (phi_xy < 0.0):
                h[i,k] = hl1 - z[i,k];
                u[i,k] = ul1;
            else:
                h[i,k] = hr1;
                u[i,k] = ur1;
        elif bslope in ['riemann_problem']:
            tend = 0.2
            hl1 = 3.0
            hr1 = 1.0
            ul1 = 0.0
            ur1 = 0.0
            xdam = 0.0
            if (D.x[i] < 0.0):
                h[i,k] = hl1 - z[i,k];
                u[i,k] = ul1;
            else:
                h[i,k] = hr1 - z[i,k];
                u[i,k] = ur1;
        elif bslope in ['2d_dambreak']:
            hl1 = 1.0
            hr1 = 0.5
            ul1 =0.0
            ur1 = 0.0
            if (D.x[i] < -0.2):
                h[i,k] = hl1 - z[i,k];
                u[i,k] = ul1;
            else:
                h[i,k] = hr1 - z[i,k];
                u[i,k] = ur1;
        elif (bslope == 'bump' and bump_config == 'emerged'):
            h[i,k] = max(hl1,z[i,k]) - z[i,k];
            u[i,k] = ul1;
        # elif bslope == 'bump_partial_cover':
        #     h[i,k] = max(hl1,z[i,k]) - z[i,k];
        #     u[i,k] = ul1;
        else :
            # h[i,k] = max(hl1,z[i,k]) - z[i,k];
            
            h[i,k] = hl1 - z[i,k];
            u[i,k] = ul1;
H = h + z;
def minmod(a,b):
    if (abs(a)<=abs(b)) and (a*b>0):
        m=a;
    elif (abs(b)<=abs(a)) and (a*b>0):
        m=b;
    elif (a*b<=0):
        
        m=0;
    return m
#function for boundary conditions for the benchmark Problems
def BC(h,H, hu,hv):
    
    if bslope == 'bump' and bump_config == 'subcritical':
        #boundary condition
        # x direction
        #left side
        h[0,:] = h[3,:];
        h[1,:] = h[2,:];
            
        H[0,:] = H[3,:];
        H[1,:] = H[2,:];
        hu[0,:] = qb#-hu[3,:];
        hu[1,:] = qb#-hu[2,:];
        hv[0,:] = hv[3,:];
        hv[1,:] = hv[2,:];              
        # right side
        h[Nx+3,:] = hb#h[Nx,:];
        h[Nx+2,:] = hb#h[Nx+1,:];
        H[Nx+3,:] = H[Nx,:];
        H[Nx+2,:] = H[Nx+1,:];
        hu[Nx+3,:] = hu[Nx,:];
        hu[Nx+2,:] = hu[Nx+1,:];
        hv[Nx+3,:] = hv[Nx,:];
        hv[Nx+2,:] = hv[Nx+1,:];              
        # y direction
        #bottom side
        h[:,0] = h[:,3];
        h[:,1] = h[:,2];              
        H[:,0] = H[:,3];
        H[:,1] = H[:,2];
        hu[:,0] = -hu[:,3];
        hu[:,1] = -hu[:,2];
        hv[:,0] = -hv[:,3];
        hv[:,1] = -hv[:,2];              
        # Top side
        h[:,Ny+3] = h[:,Ny];
        h[:,Ny+2] = h[:,Ny+1];
        H[:,Ny+3] = H[:,Ny];
        H[:,Ny+2] = H[:,Ny+1];
        hu[:,Ny+3] = -hu[:,Ny];
        hu[:,Ny+2] = -hu[:,Ny+1];
        hv[:,Ny+3] = -hv[:,Ny];
        hv[:,Ny+2] = -hv[:,Ny+1];
    elif (bslope == 'bump' and bump_config != 'subcritical') or bc == 'noflow':
        #left side
        if bump_config in ['noflow','emerged']:
            h[0,:] = hb
            h[1,:] = hb
            h[Nx+3,:] = hb;
            h[Nx+2,:] = hb;
        else:
            h[0,:] = h[:,3];
            h[1,:] = h[:,2];
            h[Nx+3,:] = h[Nx,:];
            h[Nx+2,:] = h[Nx+1,:];
            
        H[0,:] = H[3,:];
        H[1,:] = H[2,:];
        hu[0,:] = -hu[3,:];
        hu[1,:] = -hu[2,:];
        hv[0,:] = -hv[3,:];
        hv[1,:] = -hv[2,:];              
        # right side
        H[Nx+3,:] = H[Nx,:];
        H[Nx+2,:] = H[Nx+1,:];
        hu[Nx+3,:] = -hu[Nx,:];
        hu[Nx+2,:] = -hu[Nx+1,:];
        hv[Nx+3,:] = -hv[Nx,:];
        hv[Nx+2,:] = -hv[Nx+1,:];              
        # y direction
        #bottom side
        h[:,0] = h[:,3];
        h[:,1] = h[:,2];              
        H[:,0] = H[:,3];
        H[:,1] = H[:,2];
        hu[:,0] = -hu[:,3];
        hu[:,1] = -hu[:,2];
        hv[:,0] = -hv[:,3];
        hv[:,1] = -hv[:,2];              
        # Top side
        h[:,Ny+3] = h[:,Ny];
        h[:,Ny+2] = h[:,Ny+1];
        H[:,Ny+3] = H[:,Ny];
        H[:,Ny+2] = H[:,Ny+1];
        hu[:,Ny+3] = -hu[:,Ny];
        hu[:,Ny+2] = -hu[:,Ny+1];
        hv[:,Ny+3] = -hv[:,Ny];
        hv[:,Ny+2] = -hv[:,Ny+1];
    elif bslope in ['bump_steady_flow' or 'bump_partial_cover']:
        h[0,:] = h[3,:];
        h[1,:] = h[2,:];              
        H[0,:] = H[3,:];
        H[1,:] = H[2,:];
        hu[0,:] = -hu[3,:];
        hu[1,:] = -hu[2,:];
        hv[0,:] = -hv[3,:];
        hv[1,:] = -hv[2,:];              
        # right side
        h[Nx+3,:] = h[Nx,:];
        h[Nx+2,:] = h[Nx+1,:];
        H[Nx+3,:] = H[Nx,:];
        H[Nx+2,:] = H[Nx+1,:];
        hu[Nx+3,:] = -hu[Nx,:];
        hu[Nx+2,:] = -hu[Nx+1,:];
        hv[Nx+3,:] = -hv[Nx,:];
        hv[Nx+2,:] = -hv[Nx+1,:];              
        # y direction
        #bottom side
        h[:,0] = h[:,3];
        h[:,1] = h[:,2];              
        H[:,0] = H[:,3];
        H[:,1] = H[:,2];
        hu[:,0] = -hu[:,3];
        hu[:,1] = -hu[:,2];
        hv[:,0] = -hv[:,3];
        hv[:,1] = -hv[:,2];              
        # Top side
        h[:,Ny+3] = h[:,Ny];
        h[:,Ny+2] = h[:,Ny+1];
        H[:,Ny+3] = H[:,Ny];
        H[:,Ny+2] = H[:,Ny+1];
        hu[:,Ny+3] = -hu[:,Ny];
        hu[:,Ny+2] = -hu[:,Ny+1];
        hv[:,Ny+3] = -hv[:,Ny];
        hv[:,Ny+2] = -hv[:,Ny+1];
    else:
        #boundary condition
        # x direction
        #left side
        h[0,:] = h[3,:];
        h[1,:] = h[2,:];              
        H[0,:] = H[3,:];
        H[1,:] = H[2,:];
        hu[0,:] = -hu[3,:];
        hu[1,:] = -hu[2,:];
        hv[0,:] = -hv[3,:];
        hv[1,:] = -hv[2,:];              
        # right side
        h[Nx+3,:] = h[Nx,:];
        h[Nx+2,:] = h[Nx+1,:];
        H[Nx+3,:] = H[Nx,:];
        H[Nx+2,:] = H[Nx+1,:];
        hu[Nx+3,:] = hu[Nx,:];
        hu[Nx+2,:] = hu[Nx+1,:];
        hv[Nx+3,:] = hv[Nx,:];
        hv[Nx+2,:] = hv[Nx+1,:];              
        # y direction
        #bottom side
        h[:,0] = h[:,3];
        h[:,1] = h[:,2];              
        H[:,0] = H[:,3];
        H[:,1] = H[:,2];
        hu[:,0] = -hu[:,3];
        hu[:,1] = -hu[:,2];
        hv[:,0] = -hv[:,3];
        hv[:,1] = -hv[:,2];              
        # Top side
        h[:,Ny+3] = h[:,Ny];
        h[:,Ny+2] = h[:,Ny+1];
        H[:,Ny+3] = H[:,Ny];
        H[:,Ny+2] = H[:,Ny+1];
        hu[:,Ny+3] = -hu[:,Ny];
        hu[:,Ny+2] = -hu[:,Ny+1];
        hv[:,Ny+3] = -hv[:,Ny];
        hv[:,Ny+2] = -hv[:,Ny+1];
    return h,H, hu,hv
time = 0.0
tstep = 0
#function to plot
def plot(tstep,H):
    ttime = time
    # hh1 = np.transpose(hh[:,:,tstep])
    hh1 = np.transpose(H[:,:])
    hh2 =hh1 - z
    zz1 = np.transpose(z)
    if tstep % 1 == 0:
        if bslope == 'bump':
            # qb = 0#4.42
            # hb = 2.0
            N = 500
            qa1 = np.full(N, qb)
            ha1 = np.zeros(N)
            zb1 = np.zeros(N)
            # Analytical Solution
            distance = np.linspace(0,lx2-lx1,N)
            for i in range(N):
                x11 = round(distance[i], 1)
                if 8.0 < x11 < 12.0:
                    zb1[i] = 0.2 - 0.05*(x11 - 10)**2  # elevation
                else:
                    zb1[i] = 0.0  # elevation
    
            for j in range(N):
                # fun = lambda h: h**3 + (zb[j]-qb**2/(2*g*2**2) - 2)*h**2 + qb**2/(2*g)
                def fun(h): return h**3 + \
                    (zb1[j]-qb**2/(2*g*hb**2) - hb)*h**2 + qb**2/(2*g)
                ha1[j] = fsolve(fun, hb)
            za1 = ha1 + zb1
            xa1 = distance
            ua1 = qa1/ha1
            
            fig1 = plt.figure(figsize=(15,10))
            ax1 = fig1.add_subplot( 111)
            ax1.scatter(xx,hh1[0,:],color='k')
            ax1.plot(xa1,za1)
            ax1.set_title(f'{bslope} time = {round(time,3)} s dx = {dx}')
            ax1.set_xlabel(f'x')
            ax1.set_ylabel(f'y')
            ax1.set_ylim([0,np.max(hh1)+0.1])
            
            plt.show()
            
        elif bslope == 'riemann_problem':
            
            xa1,za1,ua1 = wetbed(h_L=hl1, h_R=hr1, u_L=ul1, u_R=ur1, t=ttime, lx1=lx1,
                                  lx2=lx2,n=Ny*Nx,dam_location=xdam,g=g)
            # ax.plot(xa1,np.full(len(xa1),ly1),zs=za1,label = 't= '+str(ttime)+'(A)',zorder=1)
            # ax.plot(xa1,np.full(len(xa1),ly2),zs=za1)
            # ax.plot(xa1,za1,label = 't= '+str(ttime)+'(A)',zorder=1)
            
            
            
            fig1 = plt.figure(figsize=(15,10))
            ax1 = fig1.add_subplot( 111)
            ax1.scatter(xx,hh1[0,:],color='k')
            ax1.plot(xa1,za1)
            ax1.set_title(f'{bslope} time = {round(time,3)} s dx = {dx}')
            ax1.set_xlabel(f'x')
            ax1.set_ylabel(f'y')
            ax1.set_ylim([0,np.max(hh1)+0.1])
            plt.show()
        else:
            fig = plt.figure(figsize=(15,10))
            ax = fig.add_subplot( 111, projection='3d')
            plotz = ax.plot_surface(X, Y,zz1,cmap='cool',alpha=0.05)
           
            plot = ax.plot_surface(X, Y, hh1,cmap='viridis')
            
            # plotz = ax.plot_surface(X[2:-2,2:-2], Y[2:-2,2:-2],zz1[2:-2,2:-2],cmap='cool',alpha=0.05)
            # plot = ax.plot_surface(X[2:-2,2:-2], Y[2:-2,2:-2], hh1[2:-2,2:-2],cmap='viridis')
            ax.set_title(f'{bslope} time = {round(time,3)} s dx = {dx}')
            ax.set_xlabel(f'x')
            ax.set_ylabel(f'y')
            # ax.view_init(0, 0) 
            fig.colorbar(plot,shrink=0.5, aspect=5)
            
        plt.show()
        
plot(tstep,H) #plots the initial condition

#for friction flux terms that are solved implicitly
def flux_implicit( t,Y_fn ):
    
    DY_imp = np.zeros(3*(Nx+4)*(Ny+4),float)
    # Y_fn = Y
    for i in range ((Nx+4)*(Ny+4)):
        y = Y_fn[i]
        
        if y < htol:
            u = 0.0
            v = 0
            # c = 0
        else:
            y = Y_fn[i]
            u = Y_fn[i+1*(Nx+4)*(Ny+4)]/y
            v = Y_fn[i+2*(Nx+4)*(Ny+4)]/y
            # c = Y_fn[i + 3*NumEle]/y
        un = np.sqrt(u*u + v*v)
        
        #original expression
        temp = 0 if y < 0 else y**(4/3)
        Rough = D.n
        Sfx = 0 if temp <= 0 else Rough**2*u*un/temp
        Sfy = 0 if temp <= 0 else Rough**2*v*un/temp
        # print('fric',i,Sfx)
        
        MODEL = 'flow'
        if MODEL == 'flow':
            E = 0
            Sd = 0
        DY_imp[i] = 0.0
        DY_imp[i+ 1*(Nx+4)*(Ny+4)] = -g*y*Sfx #+ (rho0 - rho)*E*u/rho/(1-mu)
        DY_imp[i+ 2*(Nx+4)*(Ny+4)] = -g*y*Sfy #+ (rho0 - rho)*E*u/rho/(1-mu)
        
    return DY_imp 

#roe solver
def roe(hl,hr,qxL,qxR,qyL,qyR,nx,ny):
    # Flux0 = np.zeros((Nx+4,Ny+4),float)
    # Flux1 = np.zeros((Nx+4,Ny+4),float)
    # Flux2 = np.zeros((Nx+4,Ny+4),float)
    df = np.zeros((Nx+2,Ny+2),float)
    dg = np.zeros((Nx+2,Ny+2),float)
    dh = np.zeros((Nx+2,Ny+2),float)
    Flux = np.zeros((5))#np.zeros((N+1,3))
    
    eigen_range = 3
    
    alpha = np.zeros((eigen_range))
    eig = np.zeros((eigen_range))
    eigL = np.zeros((eigen_range))
    eigR = np.zeros((eigen_range))
    eigvec = np.zeros((eigen_range*eigen_range))
    waveL = np.zeros((eigen_range))
    waveR = np.zeros((eigen_range))
    
    
    if hl <=0 and hr <=0:
        uRoe = 0
        vRoe = 0
    elif (hl <=0):
        uRoe = qxR/hr/2.0
        vRoe = qyR/hr/2.0
        
    
    elif (hr <=0):
        uRoe = qxL/hl/2.0
        vRoe = qyL/hl/2.0
        
    else:
        uRoe = (qxL/np.sqrt(hl) + qxR/np.sqrt(hr))/(np.sqrt(hl) + np.sqrt(hr))
        vRoe = (qyL/np.sqrt(hl) + qyR/np.sqrt(hr))/(np.sqrt(hl) + np.sqrt(hr))

    hRoe = (hl + hr)/2.0
    # print(hRoe)
    # hRoe = abs(hRoe)
    aRoe = np.sqrt(g*hRoe)
    # print(aRoe,hRoe,uRoe)
    
    # //Alpha coefficients based on Roe's average
    if(aRoe == 0):
        for k in range(eigen_range):
            alpha[k] = 0;

    else:
        # Ahmad, M.F., Mamat, M., Nik, W.B.W., Kartono, A., 2013. 
        # Numerical method for dam break problem by using Godunov approach

        alpha[0] = ((uRoe + aRoe)*(hr - hl) - (qxR - qxL))/(2.0*aRoe)
        alpha[1] = (qyR - qyL) - vRoe*(hr - hl)
        alpha[2] = ((-(uRoe - aRoe))*(hr - hl) + (qxR - qxL))/(2.0*aRoe)
    
    # //Eigenvectors based on Roe's average
    #1d
    # eigvec[0] = 1
    # eigvec[1] = 1
    # eigvec[2] = uRoe - aRoe
    # eigvec[3] = uRoe + aRoe
    
    #2d
    eigvec[0] = 1
    eigvec[1] = 0
    eigvec[2] = 1
    eigvec[3] = uRoe - aRoe
    eigvec[4] = 0
    eigvec[5] = uRoe + aRoe
    eigvec[6] = vRoe
    eigvec[7] = 0
    eigvec[8] = vRoe
    


    # Eigenvalues based on Roe's average
    
    eig[0] = uRoe - aRoe
    eig[1] = uRoe
    eig[2] = uRoe + aRoe
    
    
    if hl !=0 and hr !=0 : #and entropy_fix == 'on'
        # /******************Entropy fix*****************************/
        #hl and hr go to zero in some cases in the lines of code below
        
        # //Left Eigenvalues
        
        eigL[0]  = (qxL)/hl - np.sqrt(g*hl)
        eigL[1] = (qxL)/hl
        eigL[2]  = (qxL)/hl + np.sqrt(g*hl)
        # print(ttime,i,hl)
    
        # //Right Eigenvalues
        
        eigR[0]  = (qxR)/hr - np.sqrt(g*hr)
        eigL[1] = (qxR)/hr
        eigR[2]  = (qxR)/hr + np.sqrt(g*hr)
        
        # //Entropy fix
        for m in range(eigen_range):
            epsilon = max(0,max(eig[m] - eigL[m],eigR[m] - eig[m]));
            # print(epsilon)
            if(abs(eig[m]) < epsilon):
                if eig[m]>0:
                    eig[m] = epsilon
                else:
                    eig[m] = -epsilon
        
        # /******************Entropy fix*****************************/
        
    # //Waves based on Roe's average
    for m in range(eigen_range):
        waveL[m]  = 0;
        waveR[m]  = 0;

    for m in range(eigen_range):
        for k in range(eigen_range):
            # //left-going waves
            if(eig[k]<0):
                # waveL[m]  = waveL[m]  + eigvec[k+m*4]*alpha[k]*eig[k];
                waveL[m]  = waveL[m]  + eigvec[k+m*eigen_range]*alpha[k]*eig[k];
                # //right-going waves
            else:
                # waveR[m]  = waveR[m]  + eigvec[k+m*4]*alpha[k]*eig[k];
                waveR[m]  = waveR[m]  + eigvec[k+m*eigen_range]*alpha[k]*eig[k];
    # print(waveL,waveR)

    if(hl<=0 and hr<=0):
        for k in range(eigen_range):
            Flux[k] = 0;
          
    elif(hl<=0 ):
        Flux[0] = (qxR*nx + qyR*ny);
        Flux[1] = (qxR*qxR/hr + g*hr*hr/2)*nx + (qxR*qyR/hr)*ny;
        Flux[2] = (qxR*qyR/hr)*nx + (qyR*qyR/hr + g*hr*hr/2.0)*ny;
      # Flux[3] = (pR*qxR*nx + pR*qyR*ny)/hr;

    elif(hr<=0):
        Flux[0] = (qxL*nx + qyL*ny)
        Flux[1] = (qxL*qxL/hl + g*hl*hl/2.0)*nx + (qxL*qyL/hl)*ny
        Flux[2] = (qxL*qyL/hl)*nx + (qyL*qyL/hl + g*hl*hl/2.0)*ny
      # Flux[3] = (pL*qxL*nx + pL*qyL*ny)/hl;

    else:
                    
        #changing it to right flux - left flux
        Flux[0] = (qxL*nx + qyL*ny) + (qxR*nx + qyR*ny);
        Flux[1] = (qxL*qxL/hl + g*hl*hl/2.0)*nx + (qxL*qyL/hl)*ny + (qxR*qxR/hr + g*hr*hr/2.0)*nx + (qxR*qyR/hr)*ny;
        Flux[2] = (qxL*qyL/hl)*nx + (qyL*qyL/hl + g*hl*hl/2.0)*ny + (qxR*qyR/hr)*nx + (qyR*qyR/hr + g*hr*hr/2.0)*ny;

        
        # print('flux',i,hl,qxL,hr,qxR,Flux[0],Flux[1])
        # print('wave',waveL[0],waveR[0],waveL[1],waveR[1])
    for k in range(eigen_range):
        #setting jacobian terms to zero
        if jacobian_terms == 'on':
            Flux[k] = Flux[k] + waveL[k] - waveR[k];
        else:
            Flux[k] = Flux[k] #+ waveL[k] - waveR[k];
       

    for k in range(eigen_range):
        Flux[k] = Flux[k]/2.0
    return Flux

#solves the problem based on solver chosen
def solver_opt(solver,step):
    if step == 0:
    # computation along x
        for j in range(2,Ny+2): #start at interface i-1/2
            for i in range(2,Nx+3):
                
                # slope limiter
                # left side
                hldotx = minmod((h[i-1,j]-h[i-2,j])/dx,(h[i,j]-h[i-1,j])/dx);
                Hldotx = minmod((H[i-1,j]-H[i-2,j])/dx,(H[i,j]-H[i-1,j])/dx);
                qxLdot = minmod((hu[i-1,j]-hu[i-2,j])/dx,(hu[i,j]-hu[i-1,j])/dx);
                qyLdot = minmod((hv[i-1,j]-hv[i-2,j])/dx,(hv[i,j]-hv[i-1,j])/dx);
                # right side
                hrdot = minmod((h[i,j]-h[i-1,j])/dx,(h[i+1,j]-h[i,j])/dx);
                Hrdot = minmod((H[i,j]-H[i-1,j])/dx,(H[i+1,j]-H[i,j])/dx);
                qxRdot = minmod((hu[i,j]-hu[i-1,j])/dx,(hu[i+1,j]-hu[i,j])/dx);
                qyRdot = minmod((hv[i,j]-hv[i-1,j])/dx,(hv[i+1,j]-hv[i,j])/dx);
                # Data reconstruction
                # left side
                hl = h[i-1,j]+0.5*dx*hldotx ;
                Hl = H[i-1,j]+0.5*dx*Hldotx;
                qxL = hu[i-1,j]+0.5*dx*qxLdot;
                qyL = hv[i-1,j]+0.5*dx*qyLdot;
                zl = Hl-hl;
                # right side
                hr = h[i,j]-0.5*dx*hrdot;
                Hr = H[i,j]-0.5*dx*Hrdot;
                qxR = hu[i,j]-0.5*dx*qxRdot;
                qyR = hv[i,j]-0.5*dx*qyRdot;
                zr = Hr-hr;

                # Single topography value
                zbx[i,j] = max(zl,zr);
                # compute h well-balanced ;
                hl = max(0,Hl-zbx[i,j]);
                hr = max(0,Hr-zbx[i,j]);
                hlx[i,j] = hl;
                hrx[i,j] = hr;
                # calc velocity u,v
                # dry bed
                if (hl<htol):
                    ul = 0.0;
                    vl = 0.0;
                # wet bed
                else:
                    ul = qxL/hl;
                    vl = qyL/hl;
                # dry bed
                if (hr<htol):
                    ur = 0.0;
                    vr = 0.0;
                else:
                    ur = qxR/hr;
                    vr = qyR/hr;
                
                al = np.sqrt(g*hl)
                ar = np.sqrt(g*hr)
                
                if solver == 'hllc':
                    # wave speed (Toro 2001 : Pg 279)
                    um = 0.5*(ul+ur)+al-ar;
                    hm = (1/g)*((0.5*(al+ar)+0.25*(ul-ur))**2);
                    am = np.sqrt(g*hm)
                    #equation 38 (https://ascelibrary.org/doi/10.1061/%28ASCE%290733-9399%282008%29134%3A4%28277%29#core-c8)
                    if (hl>0) and (hr>0):
                        Sl = min(ul-(al),um-(am));
                        Sr = max(ur+(ar),um+(am));
                        if(hr*(ur-Sr))==(hl*(ul-Sl)):
                        	Sm = 0;
                        else:
                            Sm = ((Sl*hr*(ur-Sr))-(Sr*hl*(ul-Sl)))/((hr*(ur-Sr))-(hl*(ul-Sl)));
                    elif (hl==0) and (hr==0):
                        Sl = 0;
                        Sm = 0;
                        Sr = 0;
                    elif (hl>0) and (hr==0):
                        Sl = ul-al;
                        Sr = ul+(2*al);
                        Sm = Sr;
                    elif (hl==0) and (hr>0):
                        Sl = ur-(2*ar);
                        Sr = ur+ar;
                        Sm = Sl;
                    #compute HLLC flux
                    # left side
                    f1l = hl*ul;
                    f2l = hl*ul**2+0.5*g*hl**2;
                    f3l = hl*ul*vl;
                    # right side
                    f1r = hr*ur;
                    f2r = hr*ur**2+0.5*g*hr**2;
                    f3r = hr*ur*vr;
                    # middle side
                    if(Sr==Sl):
                        
                        f1m=0;
                        f2m=0;
                    else:
                        f1m = (Sr*f1l-Sl*f1r+Sl*Sr*(hr-hl))/(Sr-Sl);
                        f2m = (Sr*f2l-Sl*f2r+Sl*Sr*(hr*ur-hl*ul))/(Sr-Sl);
                    
                    # compute F flux
                    if (Sl==0) and (Sr==0):
                        f1[i,j] = 0;
                        f2[i,j] = 0;
                        f3[i,j] = 0;
                    elif (Sl>=0):
                        f1[i,j] = f1l;
                        f2[i,j] = f2l;
                        f3[i,j] = f3l;
                    elif (Sl<=0) and (Sm>=0):
                        f1[i,j] = f1m;
                        f2[i,j] = f2m;
                        f3[i,j] = f1m*vl;
                    elif (Sm<=0) and (Sr>=0):
                        f1[i,j] = f1m;
                        f2[i,j] = f2m;
                        f3[i,j] = f1m*vr;
                    else: #(Sr<=0)
                        f1[i,j] = f1r;
                        f2[i,j] = f2r;
                        f3[i,j] = f3r;
                elif solver == 'roe':
                    Flux = roe(hl,hr,qxL,qxR,qyL,qyR,nx=1,ny=0)
                        
                    f1[i,j] = Flux[0]
                    f2[i,j] = Flux[1]
                    f3[i,j] = Flux[2]
                    
        
        return f1,f2,f3#,g1,g2,g3
    if step == 1:
        # computation along y
        for i in range(2,Nx+2): #start at interface i-1/2
            for j in range(2,Ny+3):
                # slope limiter
                # left side
                hldoty = minmod((h[i,j-1]-h[i,j-2])/dy,(h[i,j]-h[i,j-1])/dy);
                Hldoty = minmod((H[i,j-1]-H[i,j-2])/dy,(H[i,j]-H[i,j-1])/dy);
                qxLdoty = minmod((hu[i,j-1]-hu[i,j-2])/dy,(hu[i,j]-hu[i,j-1])/dy);
                qyLdoty = minmod((hv[i,j-1]-hv[i,j-2])/dy,(hv[i,j]-hv[i,j-1])/dy);
                # right side
                hrdoty = minmod((h[i,j]-h[i,j-1])/dy,(h[i,j+1]-h[i,j])/dy);
                Hrdoty = minmod((H[i,j]-H[i,j-1])/dy,(H[i,j+1]-H[i,j])/dy);
                qxRdoty = minmod((hu[i,j]-hu[i,j-1])/dy,(hu[i,j+1]-hu[i,j])/dy);
                qyRdoty = minmod((hv[i,j]-hv[i,j-1])/dy,(hv[i,j+1]-hv[i,j])/dy);
                # Data reconstruction
                # left side
                hl = h[i,j-1]+0.5*dy*hldoty;
                Hl = H[i,j-1]+0.5*dy*Hldoty;
                qxL = hu[i,j-1]+0.5*dy*qxLdoty ;
                qyL = hv[i,j-1]+0.5*dy*qyLdoty ;
                zl = Hl-hl;
                # right side
                hr = h[i,j]-0.5*dy*hrdoty;
                Hr = H[i,j]-0.5*dy*Hrdoty;
                qxR = hu[i,j]-0.5*dy*qxRdoty;
                qyR = hv[i,j]-0.5*dy*qyRdoty;
                zr = Hr-hr;
                
                # Single topography value
                zby[i,j]=max(zl,zr);
                # compute h wellbance ;
                hl = max(0,Hl-zby[i,j]);
                hr = max(0,Hr-zby[i,j]);
                hly[i,j] = hl;
                hry[i,j] = hr;
                
                # calc velocity u,v
                # dry bed
                if (hl<htol):
                     ul = 0.0;
                     vl = 0.0;
                # wet bed
                else:
                     ul = qxL/hl;
                     vl = qyL/hl;
                # dry bed
                if (hr<htol):
                     ur = 0.0;
                     vr = 0.0;
                else:
                    ur = qxR/hr;
                    vr = qyR/hr;
                    
                al = np.sqrt(g*hl)
                ar = np.sqrt(g*hr)
                
                if solver == 'hllc':
                    um = 0.5*(ul+ur)+al-ar;
                    hm = (1/g)*((0.5*(al+ar)+0.25*(ul-ur))**2);
                    am = np.sqrt(g*hm)
                    # wave speed
                    vm=0.5*(vl+vr)+al-ar;
                    hm=(1/g)*((0.5*(al+ar)+0.25*(vl-vr))**2);
                    if (hl>0) and (hr>0):
                        Sl = min(vl-(al),vm-(am));
                        Sr = max(vr+(ar),vm+(am));
                        if (hr*(vr-Sr))==(hl*(vl-Sl)):
                            Sm = 0;
                        else:
                            Sm = ((Sl*hr*(vr-Sr))-(Sr*hl*(vl-Sl)))/((hr*(vr-Sr))-(hl*(vl-Sl)));
                        
                    elif (hl==0) and (hr==0):
                        Sl = 0;
                        Sm = 0;
                        Sr = 0;
                    elif (hl>0) and (hr==0):
                        Sl = vl-al;
                        Sr = vl+(2*al);
                        Sm = Sr;
                    elif (hl==0) and (hr>0):
                        Sl = vr-(2*ar);
                        Sr = vr+ar;
                        Sm = Sl;
                    # left side
                    g1l = hl*vl;
                    g2l = hl*ul*vl;
                    g3l = hl*vl**2+0.5*g*hl**2;
                    # right side
                    g1r = hr*vr;
                    g2r = hr*ur*vr;
                    g3r = hr*vr**2+0.5*g*hr**2;
                    # middle side
                    if(Sr==Sl):
                        g1m = 0;
                        g3m = 0;
                    else:
                        g1m = (Sr*g1l-Sl*g1r+Sl*Sr*(hr-hl))/(Sr-Sl);
                        g3m = (Sr*g3l-Sl*g3r+Sl*Sr*(hr*vr-hl*vl))/(Sr-Sl);
                    
                    # compute G flux
                    if (Sl==0) and (Sr==0):
                        g1[i,j] = 0;
                        g2[i,j] = 0;
                        g3[i,j] = 0;
                    elif (Sl>=0):
                        g1[i,j] = g1l;
                        g2[i,j] = g2l;
                        g3[i,j] = g3l;
                    elif (Sl<=0) and (Sm>=0):
                        g1[i,j] = g1m;
                        g2[i,j] = g1m*ul;
                        g3[i,j] = g3m;
                    elif (Sm<=0) and (Sr>=0):
                        g1[i,j] = g1m;
                        g2[i,j] = g1m*ur;
                        g3[i,j] = g3m;
                    else: #(Sr<=0)
                        g1[i,j] = g1r;
                        g2[i,j] = g2r;
                        g3[i,j] = g3r;
 
                elif solver == 'roe':
                    Flux = roe(hl,hr,qxL,qxR,qyL,qyR,nx=0,ny=1)
                    g1[i,j] = Flux[0]
                    g2[i,j] = Flux[1]
                    g3[i,j] = Flux[2]
                
    
        return g1,g2,g3 #f1,f2,f3,

while time < tend:
    print(time)
    h[h<0] = 0 # set any negative values of h to zero
    h,H, hu,hv = BC(h,H, hu,hv) #set the Boundary condition
    
#--------------------------------------------------------------------------      
 # calc velocity u,v
    for i in range(0,Nx+4):
        for j in range(0,Ny+4):
            if (h[i,j]<=htol):
                u[i,j]=0.0;
                v[i,j]=0.0;
            else:
                u[i,j]=hu[i,j]/h[i,j];
                v[i,j]=hv[i,j]/h[i,j];

    #computation of time step
    dt = 0.1;
       
    dt = C*min(dx/np.max(abs((u))+np.sqrt(g*h)),
                    dy/np.max(abs((v))+np.sqrt(g*h)));
    # dt = 0.001
    old_ttime = time # keept track of previous timestep for implicit
    time = time + dt
    tstep = tstep+1
    
    jacobian_terms = 'on'
    
    f1,f2,f3 = solver_opt(solver,step = 0) #,g1,g2,g3,
    # g1,g2,g3 = solver_opt(solver,step = 1) #f1,f2,f3,
    
    #Updating of solution at dt/2 using the fluxes along x
    # this is related to muscl hancock scheme
    for i in range(2,Nx+2):
        for j in range(2,Ny+2):
            #source term
            sox = ((hlx[i+1,j]+hrx[i,j])/2)*(((zbx[i+1,j]-zbx[i,j])/dx)); 
            soy = ((hly[i,j+1]+hry[i,j])/2)*(((zby[i,j+1]-zby[i,j])/dy));

            h[i,j] = h[i,j]-(0.5)*(dt/dx)*(f1[i+1,j]-f1[i,j])-(0.5)*(dt/dy)*(g1[i,j+1]-g1[i,j]);
            if solver == 'hllc':
                h[i,j] = h[i,j]-(0.5)*(dt/dx)*(f1[i+1,j]-f1[i,j])-(0.5)*(dt/dy)*(g1[i,j+1]-g1[i,j]);
                hu[i,j] = hu[i,j]-(0.5)*(dt/dx)*(f2[i+1,j]-f2[i,j])-(0.5)*(dt/dy)*(g2[i,j+1]-g2[i,j])-(0.5)*dt*g*sox;
                hv[i,j] = hv[i,j]-(0.5)*(dt/dx)*(f3[i+1,j]-f3[i,j])-(0.5)*(dt/dy)*(g3[i,j+1]-g3[i,j])-(0.5)*dt*g*soy;
            else:
                h[i,j] = h[i,j]-(0.5)*(dt/dx)*(f1[i+1,j]-f1[i,j])-(0.5)*(dt/dy)*(g1[i,j+1]-g1[i,j]);
                hu[i,j] = hu[i,j]-(0.5)*(dt/dx)*(f2[i+1,j]-f2[i,j])-(0.5)*(dt/dy)*(g2[i,j+1]-g2[i,j])-(0.5)*dt*g*sox;
                hv[i,j] = hv[i,j]-(0.5)*(dt/dx)*(f3[i+1,j]-f3[i,j])-(0.5)*(dt/dy)*(g3[i,j+1]-g3[i,j])-(0.5)*dt*g*soy;
    H = h+z;
    g1,g2,g3 = solver_opt(solver,step = 1) #f1,f2,f3,
    
    #Updating of solution at n+1 using the fluxes along y
    for i in range(2,Nx+2):
        for j in range(2,Ny+2):
            #source term
            sox = ((hlx[i+1,j]+hrx[i,j])/2)*(((zbx[i+1,j]-zbx[i,j])/dx)); 
            soy = ((hly[i,j+1]+hry[i,j])/2)*(((zby[i,j+1]-zby[i,j])/dy));

            if solver == 'hllc':
                h[i,j] = h[i,j]-(0.5)*(dt/dx)*(f1[i+1,j]-f1[i,j])-(0.5)*(dt/dy)*(g1[i,j+1]-g1[i,j]);
                hu[i,j] = hu[i,j]-(0.5)*(dt/dx)*(f2[i+1,j]-f2[i,j])-(0.5)*(dt/dy)*(g2[i,j+1]-g2[i,j])-(0.5)*dt*g*sox;
                hv[i,j] = hv[i,j]-(0.5)*(dt/dx)*(f3[i+1,j]-f3[i,j])-(0.5)*(dt/dy)*(g3[i,j+1]-g3[i,j])-(0.5)*dt*g*soy;
            else:
                h[i,j] = h[i,j]-(0.5)*(dt/dx)*(f1[i+1,j]-f1[i,j])-(0.5)*(dt/dy)*(g1[i,j+1]-g1[i,j]);
                hu[i,j] = hu[i,j]-(0.5)*(dt/dx)*(f2[i+1,j]-f2[i,j])-(0.5)*(dt/dy)*(g2[i,j+1]-g2[i,j])-(0.5)*dt*g*sox;
                hv[i,j] = hv[i,j]-(0.5)*(dt/dx)*(f3[i+1,j]-f3[i,j])-(0.5)*(dt/dy)*(g3[i,j+1]-g3[i,j])-(0.5)*dt*g*soy;
    
    # manning's roughness greater than zero requires implicit solver
    if D.n !=0:
        #implicit section
        Yinit = np.vstack([h,hu,hv])
        method_name = ['Radau','BDF','LSODA','RK45','RK23','DOP853'] #
        method_name = method_name[1]
    
    
        abs_tol=1e-4
        rel_tol=1e-2
    
        # start = time.time()
        
        y0 = Yinit.ravel()#np.reshape(Yinit,((Nx+4)*3*(Nx+4)))
        # t_span = np.array([0,t_final[-1]])#np.array([0,t_final[-1]])
        t_span = np.array([old_ttime,time])#np.array([0,t_final[-1]])
        dns_out = False
        sol_radau = solve_ivp(flux_implicit, t_span,y0,method=method_name,dense_output = dns_out,atol=abs_tol,rtol=rel_tol,max_step=0.1) # 
        
        # end = time.time()
        # print('Solving took: ' + str((end-start)/60) + ' min')
        
        it = -1 #final calculation
        # negs = np.where(sol_radau.y[:,it] < 0)[0]
        # sol_radau.y[negs,it] = 0
        NumEle = (Nx+4)*(Ny+4)
        hu_new = sol_radau.y[NumEle:2*NumEle,-1]
        hv_new = sol_radau.y[2*NumEle:3*NumEle,-1]
        h_new = sol_radau.y[0:NumEle,-1]
        
        h_new = h_new.reshape(Nx+4,Ny+4)
        hu_new = hu_new.reshape(Nx+4,Ny+4)
        hv_new = hv_new.reshape(Nx+4,Ny+4)
        
        width = 1
        # Vnew = (Qnew)/ynew
        #bc coming from explicit update
        for i in range(2,Nx+2):
            for j in range(2,Ny+2): 
                h[i,j] = h_new[i,j] 
                hu[i,j] = hu_new[i,j] 
                hv[i,j] = hv_new[i,j] 

        Y_imp = np.copy(sol_radau.y[0:3*NumEle,-1])
    H = h + z;
    #plot statement with indent like this plots at each timestep
    #useful for a quick check
    #remove indent if you only need to plot at final
    plot(tstep,H)