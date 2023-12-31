# analytical solutions from supplementary material of Leakey et al., 2020
import numpy as np
from scipy.optimize import fsolve

def depth(h, h_L, h_R, u_L, u_R, g):
    
    a_L = np.sqrt(g*h_L)
    a_R = np.sqrt(g*h_R)
    a = np.sqrt(g*h)
    du = u_R - u_L
    
    if h <= h_L:            #left rarefaction
        f_L = 2*(a-a_L)
    else:                   #left shock
        f_L = (h-h_L)*np.sqrt(0.5*g*(h + h_L)/(h*h_L))
        
    if h <= h_R:            #right rarefaction
        f_R = 2*(a-a_R)
    else:                   #right shock
        f_R = (h-h_R)*np.sqrt(0.5*g*(h + h_R)/(h*h_R))
        
    f = f_L + f_R + du
    u = 0.5*(u_L + u_R) + 0.5*(f_R - f_L)
    
    return(f,u)

# def wetbed(h_L, h_R, u_L, u_R, t, l,dx,g=9.81):
# def wetbed(h_L, h_R, u_L, u_R, t, l,n,dam_location,g=9.81):
def wetbed(h_L, h_R, u_L, u_R, t, lx1,lx2,n,dam_location,g):
    #celerity
    a_L = np.sqrt(g*h_L)
    a_R = np.sqrt(g*h_R)
    
    #find star region
    h0 = (0.5*(a_L+a_R)-0.25*(u_R-u_L))**2/g
    h_star = fsolve(lambda h: depth(h,h_L,h_R,u_L,u_R,g)[0],h0)[0]
    a_star = np.sqrt(g*h_star)
    u_star = depth(h_star,h_L,h_R,u_L,u_R,g)[1]
    
    #determine wave structure
    if h_star <= h_L:       #left rarefaction
        S_HL = u_L - a_L
        S_TL = u_star - a_star
        S_L = np.nan
    else:                   #left shock
        q_L = np.sqrt(0.5*((h_star + h_L)*h_star)/(h_L**2))
        S_L = u_L - a_L*q_L
        S_HL = np.nan
        S_TL = np.nan
        
    if h_star <= h_R:       #right rarefaction
        S_HR = u_R + a_R
        S_TR = u_star + a_star
        S_R = np.nan
    else:                   #right shock
        q_R = np.sqrt(0.5*((h_star + h_R)*h_star)/(h_R**2))
        S_R = u_R + a_R*q_R
        S_HR = np.nan
        S_TR = np.nan
    
    #sampling the solution
    # N = 4*l + 1
    # N =  int(l/dx)+1 #nodes
    N = n
    # x_array = np.linspace(start = 0, stop = l, num = N)
    # x_array = np.linspace(start = -l/2, stop = l/2, num = N)
    
    l = (lx2-lx1)
    num = dam_location/l
    num1 = 1.0 - num
   
    # x_array = np.linspace(start = -l*num, stop = l*num1, num = N)
    x_array = np.linspace(start = lx1, stop = lx2, num = N)
    h_array = np.zeros(N)
    u_array = np.zeros(N)
    for i in range(N):
        x = x_array[i]+abs(dam_location)
        if x/t < np.nanmin([S_L,S_HL]):
            #left of left wave
            h_array[i] = h_L
            u_array[i] = u_L
        elif x/t > np.nanmax([S_R,S_HR]):
            #right of right wave
            h_array[i] = h_R
            u_array[i] = u_R
        elif x/t > np.nanmin([S_L,S_TL]) and x/t < np.nanmax([S_R,S_TR]):
            #star region
            h_array[i] = h_star
            u_array[i] = u_star
        elif x/t > S_HL and x/t < S_TL:
            #left fan
            h_array[i] = ((1/3)*(u_L + 2*a_L - x/t))**2/g         
            u_array[i] = (1/3)*(u_L + 2*a_L + 2*x/t) 
        elif x/t > S_TR and x/t < S_HR:
            #right fan
            h_array[i] = ((1/3)*(-u_R + 2*a_R + x/t))**2/g         
            u_array[i] = (1/3)*(u_R - 2*a_R + 2*x/t)
    # return(x_array + l/2, h_array, u_array)
    # return(x_array + l*num, h_array, u_array)
    return(x_array, h_array, u_array)
