#import required libraries
import numpy as np
import matplotlib.pyplot as plt

#define the required parameters
length = 4
stepsize = 0.1
gamma = 1.4
R = 287
density = 1.125
initial_temp = [300,500,700,900,1100]
slope = [0,-50,-100,-150,-200]
N = int(input("Enter the harmonics value from this array [0,1,2,3,4] > "))
n = []
n.append(N)
print("Eigen frequency are for n = {}".format(n[len(n)-1]))
x = np.arange(0,length+stepsize,stepsize)
y_store5, u_store5, indexes_store = [],[],[]
for m in range(len(n)):
    y_store4, u_store4, index_store = [],[],[]
    for l in range(len(initial_temp)):
        end_temp = initial_temp[0]+slope[0]*length
        omega_store, freq_store, speed_store, temp_store = [], [], [], []
        for i in x:
            instantaneous_temp = initial_temp[l]+slope[l]*i
            instantaneous_speed = np.sqrt(gamma*R*instantaneous_temp)
            freq = (2*n[m]+1)*instantaneous_speed/(4*length)
            omega = 2*np.pi*freq
            freq_store.append(freq)
            omega_store.append(omega)
            speed_store.append(instantaneous_speed)
            temp_store.append(instantaneous_temp)
        y = 2000
        z = 0

        y_store1, y_store2, z_store1, z_store2, u_store1, u_store2 = [],[],[],[],[],[]
        for i in range(len(x)):
            y, z = 2000, 0
            for j in range(len(x)): 
                
                k = omega_store[i]/speed_store[j]

                k1 = z
                l1 = -(slope[l]/(initial_temp[l]+slope[l]*x[j]))*z - (k**2)*y

                k2 = z+(stepsize*l1/2)
                l2 = -(slope[l]/(initial_temp[l]+slope[l]*(x[j]+stepsize/2)))*(z+(stepsize*l1/2)) - (k**2*(y+(stepsize*k1/2)))

                k3 = z+(stepsize*l2/2)
                l3 = -(slope[l]/(initial_temp[l]+slope[l]*(x[j]+stepsize/2)))*(z+(stepsize*l2/2)) - (k**2*(y+(stepsize*k2/2)))

                k4 = z+(stepsize*l3)
                l4 = -(slope[l]/(initial_temp[l]+slope[l]*(x[j]+stepsize)))*(z+(stepsize*l3)) - (k**2*(y+(stepsize*k3)))
                
                avg_k = stepsize*(k1+2*k2+2*k3+k4)/6
                avg_l = stepsize*(l1+2*l2+2*l3+l4)/6
                
                u = - z * temp_store[j] / (end_temp*density*omega_store[i])
                
                y_store1.append(y)
                z_store1.append(z)
                u_store1.append(u)
                
                y = y + avg_k
                z = z + avg_l
            y_store2.append(y_store1)
            z_store2.append(z_store1)
            u_store2.append(u_store1)
            y_store1, z_store1, u_store1= [],[],[]
        y_store4.append(y_store2)
        u_store4.append(u_store2)
        y_store3 = np.array(y_store2)
        g = np.abs(y_store3[:,-1])
        index = np.where(g== min(g))[0][0]
        index_store.append(index)
        #print("min value of Pressure = {}, frequency =  initial temp = {} and slope = {} at index = for n = {}".format(min(g), initial_temp[l], slope[l], n[m]))
        print("{}".format(round(freq_store[index],5)))
        
    indexes_store.append(index_store)
    y_store5.append(y_store4)
    u_store5.append(u_store4)

plt.figure(figsize=[15,10])
plt.grid()
plt.plot(x,np.abs(y_store5[len(n)-1][0][indexes_store[len(n)-1][0]]))
plt.plot(x,np.abs(y_store5[len(n)-1][1][indexes_store[len(n)-1][1]]))
plt.plot(x,np.abs(y_store5[len(n)-1][2][indexes_store[len(n)-1][2]]))
plt.plot(x,np.abs(y_store5[len(n)-1][3][indexes_store[len(n)-1][3]]))
plt.plot(x,np.abs(y_store5[len(n)-1][4][indexes_store[len(n)-1][4]]))
plt.xlabel("Length of tube (m)")
plt.ylabel("Acoustic pressure (Pa)")
plt.xlim([0,4])
plt.ylim(0)
plt.title("Acoustic pressure (Pa) for n = {}".format(n[len(n)-1]))
plt.legend(['Initial Temp = {} at slope = {}'.format(initial_temp[0],slope[0]),'Initial Temp = {} at slope = {}'.format(initial_temp[1],slope[1]),'Initial Temp = {} at slope = {}'.format(initial_temp[2],slope[2]),'Initial Temp = {} at slope = {}'.format(initial_temp[3],slope[3]),'Initial Temp = {} at slope = {}'.format(initial_temp[4],slope[4])], loc='upper left')
plt.show()

plt.figure(figsize=[15,10])
plt.grid()
plt.plot(x,np.abs(u_store5[len(n)-1][0][indexes_store[len(n)-1][0]]))
plt.plot(x,np.abs(u_store5[len(n)-1][1][indexes_store[len(n)-1][1]]))
plt.plot(x,np.abs(u_store5[len(n)-1][2][indexes_store[len(n)-1][2]]))
plt.plot(x,np.abs(u_store5[len(n)-1][3][indexes_store[len(n)-1][3]]))
plt.plot(x,np.abs(u_store5[len(n)-1][4][indexes_store[len(n)-1][4]]))
plt.xlabel("Length of tube (m)")
plt.ylabel("Acoustic velocity(m/s)")
plt.xlim([0,4])
plt.ylim(0)
plt.title("Acoustic velocity (Pa) for n = {}".format(n[len(n)-1]))
plt.legend(['Initial Temp = {} at slope = {}'.format(initial_temp[0],slope[0]),'Initial Temp = {} at slope = {}'.format(initial_temp[1],slope[1]),'Initial Temp = {} at slope = {}'.format(initial_temp[2],slope[2]),'Initial Temp = {} at slope = {}'.format(initial_temp[3],slope[3]),'Initial Temp = {} at slope = {}'.format(initial_temp[4],slope[4])], loc='upper right')
plt.show()