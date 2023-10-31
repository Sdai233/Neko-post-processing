#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jul  7 12:12:36 2023

@author: sdai
"""

## fully-developed flow is assumed, i.e., d()/dx=0, d()/dz=0.
## therefore only a function of y.

import numpy as np

## the input for this module are as follows:
## all the field file here are .bin file extracted from the .fld field file
## an example is shown in mean_example.py
## mean_data: the mean velocity file
## stats_data: the statistics file
## grad_data: the derivative of the mean velocity
## derivative_data1 to derivative_data9: the derivative obtained from the postprocess_fluid_stats subroutine 
## nu, rho: defined as cases
## difference: to fix mean pressure. The value of 'difference' is the absolute value of the mean pressure at the wall
   

def comp(mean_data, stats_data, grad_data, derivative_data1, derivative_data2, derivative_data3,
         derivative_data4, derivative_data5, derivative_data6, derivative_data7, 
         derivative_data8, derivative_data9, nu, rho, difference):
    
    # Mean velocity
    U = mean_data[0,:]
    V = mean_data[1,:]
    W = mean_data[2,:]

    # Reynolds stress tensor
    uu = stats_data[0,:] - U*U
    vv = stats_data[1,:] - V*V
    ww = stats_data[2,:] - W*W
    uv = stats_data[3,:] - U*V
    uw = stats_data[5,:] - U*W
    vw = stats_data[6,:] - V*W

    # Mean, RMS, skewness and flatness of pressure
    P = mean_data[3,:]
    pp = stats_data[4,:] - P*P
    ppp = stats_data[20,:] - 3*P*pp - P*P*P
    pppp = stats_data[21,:] - 4*P*ppp - 6*P*P*pp - P*P*P*P
    #Normalize pressure
    prms = np.sqrt(pp)
    pskew = np.divide(ppp,pp**1.5)
    pflat = np.divide(pppp,pp**2)

    # Skewness&Flatness tensor
    uuu = stats_data[7,:] - 3*U*uu - U*U*U 
    vvv = stats_data[8,:] - 3*V*vv - V*V*V
    www = stats_data[9,:] - 3*W*ww - W*W*W
    uuv = stats_data[10,:] - 2*U*uv - V*uu - U*U*V
    uuw = stats_data[11,:] - 2*U*uw - W*uu - U*U*W
    uvv = stats_data[12,:] - 2*V*uv - U*vv - V*V*U
    vvw = stats_data[14,:] - 2*V*vw - W*vv - V*V*W
    uww = stats_data[15,:] - 2*W*uw - U*ww - W*W*U
    vww = stats_data[16,:] - 2*W*vw - V*ww - W*W*V
    uvw = stats_data[13,:] - U*vw - V*uw - W*uv - U*V*W

    uuuu = stats_data[17,:] - 4*U*uuu - 6*U*U*uu - U*U*U*U
    vvvv = stats_data[18,:] - 4*V*vvv - 6*V*V*vv - V*V*V*V
    wwww = stats_data[19,:] - 4*W*www - 6*W*W*ww - W*W*W*W
    
    
    
    uskew = np.divide(uuu,uu**1.5)
    uflat = np.divide(uuuu,uu**2)
    
    vskew = np.divide(vvv,vv**1.5)
    vflat = np.divide(vvvv,vv**2)
    
    wskew = np.divide(www,ww**1.5)
    wflat = np.divide(wwww,ww**2)
    
    # Velocity gradient
    # dUdx = grad_data[0,:]
    dUdy = grad_data[1,:]
    # dUdz = grad_data[2,:]
    # dVdx = grad_data[4,:]
    dVdy = grad_data[3,:]
    # dVdz = grad_data[5,:]
    # dWdx = grad_data[6,:]
    dWdy = grad_data[7,:]
    # dWdz = grad_data[8,:]
    
    dUdx = np.zeros_like(dUdy)
    dVdx = np.zeros_like(dUdy)
    dWdx = np.zeros_like(dUdy)
    
    dUdz = np.zeros_like(dUdx)
    dVdz = np.zeros_like(dUdx)
    dWdz = np.zeros_like(dUdx)

    # Production tensor
    # Pxx = -2.0*(uu*dUdx + uv*dUdy)
    # Pyy = -2.0*(uv*dVdx + vv*dVdy)
    # Pzz = -2.0*(uw*dWdx + vw*dWdy)
    # Pxy = -(uu*dVdx + uv*dVdy + uv*dUdx + vv*dUdy)
    # Pyz = -(uv*dWdx + vv*dWdy + uw*dVdx + vw*dVdy)
    # Pxz = -(uu*dWdx + uv*dWdy + uw*dUdx + vw*dUdy)
    
    # production tensor
    
    Pxx = -2.0*(uu*dUdx + uv*dUdy)
    Pyy = -2.0*(uv*dVdx + vv*dVdy)
    Pzz = -2.0*(uw*dWdx + vw*dWdy)
    Pxy = -(uu*dVdx + uv*dVdy + uv*dUdx + vv*dUdy)
    Pyz = -(uv*dWdx + vv*dWdy + uw*dVdx + vw*dVdy)
    Pxz = -(uu*dWdx + uv*dWdy + uw*dUdx + vw*dUdy)
    
    
    production = Pxx+Pyy+Pzz+Pxy+Pyz+Pxz
    production = production/2


    # Dissipation tensor
    e11_tot = stats_data[34,:]
    e22_tot = stats_data[35,:]
    e33_tot = stats_data[36,:]
    e12_tot = stats_data[37,:]
    e13_tot = stats_data[38,:] 
    e23_tot = stats_data[39,:]

    e11 = e11_tot - dUdx*dUdx - dUdy*dUdy
    e22 = e22_tot - dVdx*dVdx - dVdy*dVdy
    e33 = e33_tot - dWdx*dWdx - dWdy*dWdy
    e12 = e12_tot - dUdx*dVdx - dUdy*dVdy
    e13 = e13_tot - dUdx*dWdx - dUdy*dWdy
    e23 = e23_tot - dVdx*dWdx - dVdy*dWdy

    Dxx = -2*nu*e11
    Dyy = -2*nu*e22
    Dzz = -2*nu*e33
    Dxy = -2*nu*e12 
    Dyz = -2*nu*e23
    Dxz = -2*nu*e13
    
    dissipation= Dxx+Dyy+Dzz+Dxy+Dyz+Dxz
    dissipation = dissipation/2
    
    
    
    # mean convection tensor
    # duudx = derivative_data1[1,:] - 2*U*dUdx
    # dvvdx = derivative_data1[4,:] - 2*V*dVdx
    # dwwdx = derivative_data1[6,:] - 2*W*dWdx
    # duvdx = derivative_data2[1,:] - U*dVdx - V*dUdx
    # duwdx = derivative_data2[6,:] - U*dWdx - W*dUdx
    # dvwdx = derivative_data2[4,:] - V*dWdx - W*dVdx

    duudy = derivative_data1[0,:] - 2*U*dUdy
    dvvdy = derivative_data1[3,:] - 2*V*dVdy
    dwwdy = derivative_data1[7,:] - 2*W*dWdy
    duvdy = derivative_data2[0,:] - U*dVdy - V*dUdy
    duwdy = derivative_data2[7,:] - U*dWdy - W*dUdy
    dvwdy = derivative_data2[3,:] - V*dWdy - W*dVdy
    
    duudx = np.zeros_like(duudy)
    dvvdx = np.zeros_like(duudy)
    dwwdx = np.zeros_like(duudy)
    duvdx = np.zeros_like(duudy)
    duwdx = np.zeros_like(duudy)
    dvwdx = np.zeros_like(duudy)
    
    Cxx = U*duudx + V*duudy
    Cyy = U*dvvdx + V*dvvdy
    Czz = U*dwwdx + V*dwwdy
    Cxy = U*duvdx + V*duvdy
    Cxz = U*duwdx + V*duwdy
    Cyz = U*dvwdx + V*dvwdy
    
    convection = Cxx+Cyy+Czz+Cxy+Cxz+Cyz
    convection = convection/2
    
    
    ## turbulence transport tensor

    
    # duuudx = derivative_data3[1,:] - 3*U*U*dUdx - 3*(U*duudx+uu*dUdx)
    # dvvudx = derivative_data4[6,:] - 2*(V*duvdx+uv*dVdx) - (U*dvvdx+vv*dUdx) - (V*V*dUdx+2*U*V*dVdx)
    # dwwudx = derivative_data9[1,:] - 2*(W*duwdx+uw*dWdx) - (U*dwwdx+ww*dUdx) - (W*W*dUdx+2*U*W*dWdx)
    # duvudx = derivative_data4[1,:] - 2*(U*duvdx+uv*dUdx) - (V*duudx+uu*dVdx) - (U*U*dVdx+2*U*V*dUdx)
    # duwudx = derivative_data4[4,:] - 2*(U*duwdx+uw*dUdx) - (W*duudx+uu*dWdx) - (U*U*dWdx+2*U*W*dUdx)
    # dvwudx = derivative_data9[6,:] - (U*dvwdx+vw*dUdx) - (V*duwdx+uw*dVdx) - (W*duvdx+uv*dWdx) - (U*V*dWdx+U*W*dVdx+V*W*dUdx)
    
        
    duuvdy = derivative_data4[0,:] - 2*(U*duvdy+uv*dUdy) - (V*duudy+uu*dVdy) - (U*U*dVdy+2*U*V*dUdy)
    dvvvdy = derivative_data3[3,:] - 3*(V*dvvdy+vv*dVdy) - 3*V*V*dVdy
    dwwvdy = derivative_data9[3,:] - 2*(W*dvwdy+vw*dWdy) - (V*dwwdy+ww*dVdy) - (W*W*dVdy+2*V*W*dWdy)
    duvvdy = derivative_data4[7,:] - 2*(V*duvdy+uv*dVdy) - (U*dvvdy+vv*dUdy) - (V*V*dUdy+2*U*V*dVdy)
    duwvdy = derivative_data9[7,:] - (U*dvwdy+vw*dUdy) - (V*duwdy+uw*dVdy) - (W*duvdy+uv*dWdy) - (U*V*dWdy+U*W*dVdy+V*W*dUdy)
    dvwvdy = derivative_data9[10,:] - 2*(V*dvwdy+vw*dVdy) - (W*dvvdy+vv*dWdy) - (V*V*dWdy+2*V*W*dVdy)
    
    duuudx = np.zeros_like(dvvvdy)
    dvvudx = np.zeros_like(dvvvdy)
    dwwudx = np.zeros_like(dvvvdy)
    duvudx = np.zeros_like(dvvvdy)
    duwudx = np.zeros_like(dvvvdy)
    dvwudx = np.zeros_like(dvvvdy)
    
    
    duuwdz = np.zeros_like(duuudx)
    dvvwdz = np.zeros_like(duuudx)
    dwwwdz = np.zeros_like(duuudx)
    duvwdz = np.zeros_like(duuudx)
    duwwdz = np.zeros_like(duuudx)
    dvwwdz = np.zeros_like(duuudx)
    
    Txx = - (duuudx + duuvdy + duuwdz)
    Tyy = - (dvvudx + dvvvdy + dvvwdz)
    Tzz = - (dwwudx + dwwvdy + dwwwdz)
    Txy = - (duvudx + duvvdy + duvwdz)
    Txz = - (duwudx + duwvdy + duwwdz)
    Tyz = - (dvwudx + dvwvdy + dvwwdz)
    
    turbulence_transport = Txx+Tyy+Tzz+Txy+Txz+Tyz
    turbulence_transport = turbulence_transport/2
    
    
    ## viscous diffusion tensor
    
    # d2uudx2 = derivative_data6[1,:] - 2*(U*derivative_data5[1,:]+dUdx*dUdx)
    # d2vvdx2 = derivative_data6[4,:] - 2*(V*derivative_data5[4,:]+dVdx*dVdx)
    # d2wwdx2 = derivative_data6[6,:] - 2*(W*derivative_data5[6,:]+dWdx*dWdx)
    # d2uvdx2 = derivative_data7[1,:] - (V*derivative_data5[1,:]+U*derivative_data5[4,:]+2*dUdx*dVdx)
    # d2uwdx2 = derivative_data7[4,:] - (U*derivative_data5[6,:]+W*derivative_data5[1,:]+2*dUdx*dWdx)
    # d2vwdx2 = derivative_data7[6,:] - (V*derivative_data5[6,:]+W*derivative_data5[4,:]+2*dVdx*dWdx)


    d2uudy2 = derivative_data6[0,:] - 2*(U*derivative_data5[0,:]+dUdy*dUdy)
    d2vvdy2 = derivative_data6[3,:] - 2*(V*derivative_data5[3,:]+dVdy*dVdy)
    d2wwdy2 = derivative_data6[7,:] - 2*(W*derivative_data5[7,:]+dWdy*dWdy)
    d2uvdy2 = derivative_data7[0,:] - (V*derivative_data5[0,:]+U*derivative_data5[3,:]+2*dUdy*dVdy)
    d2uwdy2 = derivative_data7[3,:] - (U*derivative_data5[7,:]+W*derivative_data5[0,:]+2*dUdy*dWdy)
    d2vwdy2 = derivative_data7[7,:] - (V*derivative_data5[7,:]+W*derivative_data5[3,:]+2*dVdy*dWdy)

    d2uudx2 = np.zeros_like(d2uudy2)
    d2vvdx2 = np.zeros_like(d2uudy2)
    d2wwdx2 = np.zeros_like(d2uudy2)
    d2uvdx2 = np.zeros_like(d2uudy2)
    d2uwdx2 = np.zeros_like(d2uudy2)
    d2vwdx2 = np.zeros_like(d2uudy2)
    

    d2uudz2 = np.zeros_like(d2uudx2)
    d2vvdz2 = np.zeros_like(d2uudx2)
    d2wwdz2 = np.zeros_like(d2uudx2)
    d2uvdz2 = np.zeros_like(d2uudx2)
    d2uwdz2 = np.zeros_like(d2uudx2)
    d2vwdz2 = np.zeros_like(d2uudx2)

    VDxx = nu*(d2uudx2 + d2uudy2 + d2uudz2)
    VDyy = nu*(d2vvdx2 + d2vvdy2 + d2vvdz2)
    VDzz = nu*(d2wwdx2 + d2wwdy2 + d2wwdz2)
    VDxy = nu*(d2uvdx2 + d2uvdy2 + d2uvdz2)
    VDxz = nu*(d2uwdx2 + d2uwdy2 + d2uwdz2)
    VDyz = nu*(d2vwdx2 + d2vwdy2 + d2vwdz2)
    
    
    
    viscous_diffusion = VDxx+VDyy+VDzz+VDxy+VDyz+VDxz
    viscous_diffusion = viscous_diffusion/2
    
    
    
    
    
    # Pressure transport
    # there is a constant for Pressure in Neko
    # the difference is used to fix this constant
    P = P - difference
    
    # dpudx = (derivative_data8[1,:]-difference*dUdx) - P*dUdx - U*grad_data[9,:]
    # dpvdx = (derivative_data8[4,:]-difference*dVdx) - P*dVdx - V*grad_data[9,:]
    # dpwdx = (derivative_data8[6,:]-difference*dWdx) - P*dWdx - W*grad_data[9,:]
    
    
    dpudy = (derivative_data8[0,:]-difference*dUdy) - P*dUdy - U*grad_data[10,:]
    dpvdy = (derivative_data8[3,:]-difference*dVdy) - P*dVdy - V*grad_data[10,:]
    dpwdy = (derivative_data8[7,:]-difference*dWdy) - P*dWdy - W*grad_data[10,:]
    
    dpudx = np.zeros_like(dpudy)
    dpvdx = np.zeros_like(dpudy)
    dpwdx = np.zeros_like(dpudy)
    
    dpudz = np.zeros_like(dpudx)
    dpvdz = np.zeros_like(dpudx)
    dpwdz = np.zeros_like(dpudx)
    
    PTxx = 2.0/rho*dpudx
    PTyy = 2.0/rho*dpvdy
    PTzz = 2.0/rho*dpwdz
    PTxy = 1.0/rho*(dpudy+dpvdx)
    PTxz = 1.0/rho*(dpudz+dpwdx)
    PTyz = 1.0/rho*(dpvdz+dpwdy)
    
    pressure_transport = PTxx+PTyy+PTzz+PTxy+PTxz+PTyz
    pressure_transport = pressure_transport/2
    
    
    
    # Pressure strain tensor
    
    pdudx = (stats_data[25,:]-difference*dUdx) - P*dUdx
    pdudy = (stats_data[26,:]-difference*dUdy) - P*dUdy
    pdudz = (stats_data[27,:]-difference*dUdz) - P*dUdz

    pdvdx = (stats_data[28,:]-difference*dVdx) - P*dVdx
    pdvdy = (stats_data[29,:]-difference*dVdy) - P*dVdy
    pdvdz = (stats_data[30,:]-difference*dVdz) - P*dVdz

    pdwdx = (stats_data[31,:]-difference*dWdx) - P*dWdx
    pdwdy = (stats_data[32,:]-difference*dWdy) - P*dWdy
    pdwdz = (stats_data[33,:]-difference*dWdz) - P*dWdz
    

    PSxx = 2.0/rho*pdudx
    PSyy = 2.0/rho*pdvdy
    PSzz = 2.0/rho*pdwdz
    PSxy = 1.0/rho*(pdudy+pdvdx)
    PSxz = 1.0/rho*(pdudz+pdwdx)
    PSyz = 1.0/rho*(pdvdz+pdwdy)

#   Construct velocity-pressure-gradient-tensor
    Pixx = -1.0/rho*(PTxx - PSxx)
    Piyy = -1.0/rho*(PTyy - PSyy)
    Pizz = -1.0/rho*(PTzz - PSzz)
    Pixy = -1.0/rho*(PTxy - PSxy)
    Pixz = -1.0/rho*(PTxz - PSxz)
    Piyz = -1.0/rho*(PTyz - PSyz)
    
  
    
#   Budget for each component of the Reynolds stress tensor 
#   Without mean convection
    Sxx = Pxx+Dxx+Txx+VDxx+Pixx
    Syy = Pyy+Dyy+Tyy+VDyy+Piyy
    Szz = Pzz+Dzz+Tzz+VDzz+Pizz
    Sxy = Pxy+Dxy+Txy+VDxy+Pixy
    Sxz = Pxz+Dxz+Txz+VDxz+Pixz
    Syz = Pyz+Dyz+Tyz+VDyz+Piyz

#   With mean convection
    Scxx = Pxx+Dxx+Txx+VDxx+Pixx+Cxx
    Scyy = Pyy+Dyy+Tyy+VDyy+Piyy+Cyy
    Sczz = Pzz+Dzz+Tzz+VDzz+Pizz+Czz
    Scxy = Pxy+Dxy+Txy+VDxy+Pixy+Cxy
    Scxz = Pxz+Dxz+Txz+VDxz+Pixz+Cxz
    Scyz = Pyz+Dyz+Tyz+VDyz+Piyz+Cyz
    

    pressure_diffusion = Pixx+Piyy+Pizz+Pixy+Pixz+Piyz
    pressure_diffusion = pressure_diffusion/2
    
    test111 = pressure_diffusion
    
    
    # Save averaged fields in a dictionary
    data = {'U': U, 'V': V,'W': W, \
            'P': P, \
            'pp': pp, \
            'pppp': pppp, \
            'uu': uu, 'vv': vv,'ww': ww, 'uv': uv, 'uw': uw, 'vw': vw,\
            'uuu':uuu, 'vvv':vvv, 'www':www,\
            'uuv':uuv, 'wwv':vww, 'uvv':uvv, 'uwv':uvw, 'vwv':vvw, \
            'uuuu':uuuu, 'vvvv':vvvv, 'wwww':wwww,\
            
            'dUdx': dUdx, 'dUdy': dUdy, \
            'dUdz': dUdz, 'dVdx': dVdx, \
            'dVdy': dVdy, 'dVdz': dVdz, \
            'dWdx': dWdx, 'dWdy': dWdy, \
            'dWdz': dWdz, \
                
            'uskew': uskew, 'uflat': uflat, \
            'vskew': vskew, 'vflat': vflat, \
            'wskew': wskew, 'wflat': wflat, \
                
            'Cxx': Cxx, 'Cyy': Cyy, 'Czz': Czz, \
            'Cxy': Cxx, 'Cxz': Cyy, 'Cyz': Czz, \
                
            'Sxx':Sxx, 'Syy':Syy, 'Szz':Szz, 'Sxy':Sxy, 'Sxz':Sxz, 'Syz':Syz, \
            'Scxx':Scxx, 'Scyy':Scyy, 'Sczz':Sczz, 'Scxy':Scxy, 'Scxz':Scxz, 'Scyz':Scyz, \
                
            'dissipation' : dissipation, 'viscous_diffusion' : viscous_diffusion, \
            'turbulence_transport' : turbulence_transport, 'pressure_transport' : pressure_transport, \
            'convection' : convection, 'production' : production, \
            'pressure_diffusion' : pressure_diffusion, \
            'test111':test111,\
           }
    return data


























