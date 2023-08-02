import numpy as np
import tkinter
import matplotlib.pyplot as plt
import random
import math

gallium_density = 5.1 #g/cm3
a_cons_gallium = np.array([[1.00000e-3,1.697e3],
[1.05613e-3,1.492e3],
[1.11540e-3,1.312e3],
[1.11541e-3,3.990e3],
[1.12877e-3,4.887e3],
[1.14230e-3,5.664e3],
[1.14231e-3,7.405e3],
[1.21752e-3,7.138e3],
[1.29770e-3,6.358e3],
[1.29771e-3,7.206e3],
[1.50000e-3,5.087e3],
[2.00000e-3,2.515e3],
[3.00000e-3,8.857e2],
[4.00000e-3,4.130e2],
[5.00000e-3,2.266e2],
[6.00000e-3,1.382e2],
[8.00000e-3,6.302e1],
[1.00000e-2,3.421e1],
[1.03671e-2,3.099e1],
[1.03672e-2,2.214e2],
[1.50000e-2,8.537e1],
[2.00000e-2,3.928e1],
[3.00000e-2,1.281e1],
[4.00000e-2,5.726e0],
[5.00000e-2,3.076e0],
[6.00000e-2,1.868e0],
[8.00000e-2,8.823e-1],
[1.00000e-1,5.197e-1],
[1.50000e-1,2.387e-1],
[2.00000e-1,1.619e-1],
[3.00000e-1,1.123e-1],
[4.00000e-1,9.325e-2],
[5.00000e-1,8.236e-2],
[6.00000e-1,7.487e-2],
[8.00000e-1,6.466e-2],
[1.00000e0,5.767e-2],
[1.25000e0,5.139e-2],
[1.50000e0,4.692e-2],
[2.00000e0,4.113e-2],
[3.00000e0,3.538e-2],
[4.00000e0,3.280e-2],
[5.00000e0,3.156e-2],
[6.00000e0,3.099e-2],
[8.00000e0,3.086e-2],
[1.00000e1,3.130e-2],
[1.50000e1,3.300e-2],
[2.00000e1,3.479e-2]],float)

from tkinter import *

def mu_coeff(single_photon_energy):
    if single_photon_energy < a_cons_gallium[0,0]:
        return a_cons_gallium[0,1]
    else:
        j=1
        while single_photon_energy > a_cons_gallium[j,0]:
            j+=1
        return (a_cons_gallium[j-1,1]*abs(a_cons_gallium[j-1,0]-single_photon_energy)+a_cons_gallium[j,1]*abs(a_cons_gallium[j,0]-single_photon_energy))/abs(a_cons_gallium[j,0]-a_cons_gallium[j-1,0])

#def show_values():
#    print (distance_scale.get(), length_scale.get())

def run():
    result = np.zeros(1)
    missing_ray_counter=0
    accurate_ray_counter=-1
    distance = float(distance_scale.get()/10)                   #cm
    length = float(length_scale.get()/10)                       #cm
    radius = float(radius_scale.get()/10)                       #cm
    photoenergy = float(photoenergy_scale.get()/1000)           #MeV
    activity_per_second = float(activity_scale.get()*37000)     #Bq
    i=0
    while i <= activity_per_second:
        theta=round(random.uniform(0, 17999))/100               #max: 179.99
        if (math.tan(theta*math.pi/180)*distance>radius):
            missing_ray_counter+=1
            #result = np.append(result, np.array([0]), axis=0)
            i+=1
        else:
            accurate_ray_counter+=1
            x_inside=radius/math.tan(theta*math.pi/180)-distance
            x_squared=math.pow(x_inside , 2)
            if x_inside > length:
                x_squared = math.pow(length, 2)
                x_inside = length
            y_squared=math.pow((x_inside)*math.tan(theta*math.pi/180),2)
            alpha_const=gallium_density*math.sqrt(x_squared+y_squared)*mu_coeff(photoenergy)
            absorbed_energy=photoenergy*(1-math.exp(-alpha_const))
            result = np.append(result, np.array([absorbed_energy]), axis=0)
            i+=1   
    total_energy = 0
    for i in range(result.size):
        total_energy =+ result[i]
    print (missing_ray_counter, "photons have missed.")
    print (accurate_ray_counter, "photons have hit.")
    print ("Detector efficiency is:", total_energy/(missing_ray_counter*photoenergy))
    bins=200
    plt.hist(result, bins, log=True)
    plt.title('Counts vs Energy Histogram')
    plt.xlabel('Detected Energy of photons (MeV)')
    plt.ylabel('Number of detected photons (Counts)')
    plt.show()

master = Tk()
distance_scale = Scale(master, from_=0, to=1000, length=400,tickinterval=100, orient=HORIZONTAL, label='Distance (mm)')
distance_scale.set(1)
distance_scale.pack()
length_scale = Scale(master, from_=0, to=1000, length=400,tickinterval=100, orient=HORIZONTAL, label='Length (mm)')
length_scale.set(1)
length_scale.pack()
radius_scale = Scale(master, from_=0, to=1000, length=400,tickinterval=100, orient=HORIZONTAL, label='Radius (mm)')
radius_scale.set(1)
radius_scale.pack()
photoenergy_scale = Scale(master, from_=0, to=20000, length=400,tickinterval=2000, orient=HORIZONTAL, label='Energy of a single photon (KeV)')
photoenergy_scale.set(1)
photoenergy_scale.pack()
activity_scale = Scale(master, from_=0, to=50, length=400,tickinterval=10, orient=HORIZONTAL, label='Source activity (μCi)')
activity_scale.set(1)
activity_scale.pack()
Button(master, text='Run for 1 second', command=run).pack()
mainloop()