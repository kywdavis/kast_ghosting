import scipy.optimize as opt
import numpy as np
import matplotlib.pyplot as plt
import pickle


def tdmc_gau((x,y), amplitude1, xo1, yo1, sigma_x1, sigma_y1, theta1,
                           amplitude2, xo2, yo2, sigma_x2, sigma_y2, theta2, offset):
    
    xo1 = float(xo1)
    yo1 = float(yo1)    
    a1 = (np.cos(theta1)**2)/(2*sigma_x1**2) + (np.sin(theta1)**2)/(2*sigma_y1**2)
    b1 = -(np.sin(2*theta1))/(4*sigma_x1**2) + (np.sin(2*theta1))/(4*sigma_y1**2)
    c1 = (np.sin(theta1)**2)/(2*sigma_x1**2) + (np.cos(theta1)**2)/(2*sigma_y1**2)
     
    g1 = amplitude1*np.exp( - (a1*((x-xo1)**2) + 2*b1*(x-xo1)*(y-yo1) + c1*((y-yo1)**2)))
    
    
    xo2 = float(xo2)
    yo2 = float(yo2)    
    a2 = (np.cos(theta2)**2)/(2*sigma_x2**2) + (np.sin(theta2)**2)/(2*sigma_y2**2)
    b2 = -(np.sin(2*theta2))/(4*sigma_x2**2) + (np.sin(2*theta2))/(4*sigma_y2**2)
    c2 = (np.sin(theta2)**2)/(2*sigma_x2**2) + (np.cos(theta2)**2)/(2*sigma_y2**2)
     
    g2 = amplitude2*np.exp( - (a2*((x-xo2)**2) + 2*b2*(x-xo2)*(y-yo2) + c2*((y-yo2)**2)))
    
    
    g = g1 + g2 + offset
    
    return g.ravel()