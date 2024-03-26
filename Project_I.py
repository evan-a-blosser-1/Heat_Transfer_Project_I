import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
###################
######### Constants
N_fin = 12
t_fin = 0.0025
L_fin = 0.05
####################
# 0deg C = 273.15 K
Delta_T = (180+273.15)-(35+273.15)
##########
# assumed heat transfer
# coefficient for air
h  = 50  # [Watts / m^2 kelvin]
# Aluminum conductivity
k  = 237 # [Watts / meter kelvin]
#####################################
######## Parameter D setting ########
D = np.linspace(0.01,0.5,1000000)



###########################
############## Longitudinal
L_c = L_fin + t_fin/2
w = np.pi * D
A_long = w*t_fin
P_long = 2*(w+ t_fin)
m = np.sqrt((2*h)/(k*t_fin))
Q_long_Fin = N_fin * np.sqrt(h*k*P_long*A_long) *Delta_T*((np.sinh(m*L_c) + (h/(m*k))*np.cosh(m*L_c))/(np.cosh(m*L_c) + (h/(m*k))*np.sinh(m*L_c)))  

# Un-finned 
A_long_fin = 2*w*L_c
Spacing_long = np.pi*D - N_fin*A_long_fin
A_long_unfin = np.pi*D*Spacing_long
Q_long_unfin = h*A_long_unfin*Delta_T
Q_long_total = Q_long_Fin + Q_long_unfin 

###########################
#################### Radial 
r_1 = D/2
r_2 = D/2 + L_fin
r_cor = r_2 + t_fin/2
A_rad_fin = 2*np.pi*(r_cor**2 - r_1**2) + 2*np.pi*r_2*t_fin
Q_rad_fin = h*A_rad_fin*Delta_T

# Un-finned 
Spacing_rad = np.pi*D - N_fin*t_fin
A_rad_unfin = np.pi * D * Spacing_rad
Q_rad_unfin = h*A_rad_unfin*Delta_T
Q_rad_total = Q_rad_fin + Q_rad_unfin




##############################
############ Equivalency check
for i in range(len(D)):
    if round(Q_rad_total[i],2) == round(Q_long_total[i],2):
        print(f'{round(Q_rad_total[i],2)} (W) is the equivalent heat X-fer @: {i}')
##################################################################
# Save Data @ Path
Data_Out = pd.DataFrame({'Longi Q': Q_long_total,
                         'Radial Q': Q_rad_total,
                         'D': D,
                        })
Data_Out_File_Name = 'Heat_Transfer_Data.csv'
Data_Out.to_csv(Data_Out_File_Name, sep=' ' ,index=False)
##################################
############## Plot ##############
fig, ax = plt.subplots(figsize=(8, 6))
ax.plot(D, Q_long_total, 
        linestyle='-', 
        label='Longitudinal')
ax.plot(D, Q_rad_total, 
        linestyle='-', 
        label='Annular')
########################
################# Labels
ax.set_xlabel('Parameter "D" (meters)',
             fontsize=14, 
             fontweight='bold',
             fontdict={'family': 'Times New Roman'})
ax.set_ylabel('\u0394 Q (Watts)',
             fontsize=14, 
             fontweight='bold',
             fontdict={'family': 'Times New Roman'})
plt.suptitle('Temperature as a function of D', 
             fontsize=24, 
             fontweight='bold',
             fontdict={'family': 'Times New Roman', 'style': 'italic'})
#########################
# Show the legend & plot
ax.legend()
plt.grid(which='both')
plt.minorticks_on()
plt.show()
##########
