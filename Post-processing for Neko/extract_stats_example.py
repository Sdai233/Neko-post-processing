#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul  4 15:41:53 2023

@author: sdai
"""

#!/usr/bin/python
# an example for extracting data from .fld file

import numpy as np
from pymech.neksuite import readnek


data = readnek('stats0_avg_z_x0.f00000')
nel = data.nel
order = 7
total_number = nel*(order+1)**3
data1 = []
y = []

print('complete')


# defining the size of the output array
size = 0
for e in range(0,nel):
    for nx in range(0, order+1):
        for ny in range(0, order+1):
            for nz in range(0, order+1):
                
                if abs(data.elem[e].pos[0,nx,ny,nz] - 4.93169) < 1e-4:
                    if abs(data.elem[e].pos[2,nx,ny,nz] - 1.32579) < 1e-4:
                        y.append(data.elem[e].pos[1,nx,ny,nz])
                        
                        size = size +1
                        
                        
                        
data_sum = np.zeros((40, size))

for field_number in range(0,40):
    
    ## different variable fields
    
    for e in range(0,nel):
        for nx in range(0, order+1):
            for ny in range(0, order+1):
                for nz in range(0, order+1):
                    
                    if abs(data.elem[e].pos[0,nx,ny,nz] - 4.93169) < 1e-4:
                        if abs(data.elem[e].pos[2,nx,ny,nz] - 1.32579) < 1e-4:
                            
                            if field_number == 0:
                                field_data = data.elem[e].vel[0,nx,ny,nz]
                            if field_number == 1:
                                field_data = data.elem[e].vel[1,nx,ny,nz]  
                            if field_number == 2:
                                field_data = data.elem[e].vel[2,nx,ny,nz]
                            if field_number == 3:
                                field_data = data.elem[e].temp[0,nx,ny,nz]
                            if field_number == 4:
                                field_data = data.elem[e].pres[0,nx,ny,nz]
                            if field_number == 5:
                                field_data = data.elem[e].scal[0,nx,ny,nz]
                            if field_number == 6:
                                field_data = data.elem[e].scal[1,nx,ny,nz]
                            if field_number == 7:
                                field_data = data.elem[e].scal[2,nx,ny,nz]
                            if field_number == 8:
                                field_data = data.elem[e].scal[3,nx,ny,nz]
                            if field_number == 9:
                                field_data = data.elem[e].scal[4,nx,ny,nz]
                            if field_number == 10:
                                field_data = data.elem[e].scal[5,nx,ny,nz]
                            if field_number == 11:
                                field_data = data.elem[e].scal[6,nx,ny,nz]
                            if field_number == 12:
                                field_data = data.elem[e].scal[7,nx,ny,nz]
                            if field_number == 13:
                                field_data = data.elem[e].scal[8,nx,ny,nz]
                            if field_number == 14:
                                field_data = data.elem[e].scal[9,nx,ny,nz]
                            if field_number == 15:
                                field_data = data.elem[e].scal[10,nx,ny,nz]
                            if field_number == 16:
                                field_data = data.elem[e].scal[11,nx,ny,nz]
                            if field_number == 17:
                                field_data = data.elem[e].scal[12,nx,ny,nz]
                            if field_number == 18:
                                field_data = data.elem[e].scal[13,nx,ny,nz]
                            if field_number == 19:
                                field_data = data.elem[e].scal[14,nx,ny,nz]
                            if field_number == 20:
                                field_data = data.elem[e].scal[15,nx,ny,nz]
                            if field_number == 21:
                                field_data = data.elem[e].scal[16,nx,ny,nz]
                            if field_number == 22:
                                field_data = data.elem[e].scal[17,nx,ny,nz]
                            if field_number == 23:
                                field_data = data.elem[e].scal[18,nx,ny,nz]
                            if field_number == 24:
                                field_data = data.elem[e].scal[19,nx,ny,nz]
                            if field_number == 25:
                                field_data = data.elem[e].scal[20,nx,ny,nz]
                            if field_number == 26:
                                field_data = data.elem[e].scal[21,nx,ny,nz]
                            if field_number == 27:
                                field_data = data.elem[e].scal[22,nx,ny,nz]
                            if field_number == 28:
                                field_data = data.elem[e].scal[23,nx,ny,nz]
                            if field_number == 29:
                                field_data = data.elem[e].scal[24,nx,ny,nz]
                            if field_number == 30:
                                field_data = data.elem[e].scal[25,nx,ny,nz]
                            if field_number == 31:
                                field_data = data.elem[e].scal[26,nx,ny,nz]
                            if field_number == 32:
                                field_data = data.elem[e].scal[27,nx,ny,nz]
                            if field_number == 33:
                                field_data = data.elem[e].scal[28,nx,ny,nz]
                            if field_number == 34:
                                field_data = data.elem[e].scal[29,nx,ny,nz]
                            if field_number == 35:
                                field_data = data.elem[e].scal[30,nx,ny,nz]
                            if field_number == 36:
                                field_data = data.elem[e].scal[31,nx,ny,nz]
                            if field_number == 37:
                                field_data = data.elem[e].scal[32,nx,ny,nz]
                            if field_number == 38:
                                field_data = data.elem[e].scal[33,nx,ny,nz]
                            if field_number == 39:
                                field_data = data.elem[e].scal[34,nx,ny,nz] 
                            
                                
                                
                            data1.append(field_data)
                            
                
    data1 = np.array(data1)         
     
    data_sum[field_number,:] = data1[:]
    data1 = []
                    
                    
                    
                  
                        
                        
    
    
print('complete')

y = np.array(y)
y.tofile('y')
data_sum.tofile('output')


    
   







