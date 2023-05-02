import csv
import numpy as np
import matplotlib.pyplot as plt
#####################
savefile = "data.csv"   #Rename to required output file
#####################
#Read in the .csv data file: 
file = open(savefile)
csvreader = csv.reader(file)
rows = []
for row in csvreader:
  rows.append(row)
file.close()

#Extract header information:
dt = float(rows[0][0]) 
T = int(rows[0][1])
Tsave = int(rows[0][2])
N = int(rows[0][3])
NumSave = int(np.ceil((T+1)/Tsave)) #Number of save entries

#Initialise storage for data:
pos = np.zeros([NumSave,N,3])
vel = np.zeros([NumSave,N,3])

#Extract and reshape the data:
data = np.zeros(3*N)

for i in range(NumSave):
      
  start = i*N + 1 #Beginning of save N=3
   
  #Vectors variables:
  for j in range(N):
    for k in range(3): 
      pos[i,j,k] = float(rows[start+j][k])
      vel[i,j,k] = float(rows[start+j][k+3])
      
#-------------------------------------------
#Plot the paths of particles on 12x12x12 grid, centre (0,0,0):
ax = plt.figure(figsize=(15,10),dpi=80).add_subplot(projection='3d')

for j in range(N):
    ax.plot(pos[:,j,0], pos[:,j,1], pos[:,j,2], '.', markersize = 5, label = f'particle{j}', alpha = 0.2)

#ax.legend()
ax.set_xlim([-6,6])
ax.set_ylim([-6,6])
ax.set_zlim([-6,6])
ax.set_title('Particles in a box')
ax.set_xlabel('$x$')
ax.set_ylabel('$y$')
ax.set_zlabel('$z$')

plt.tight_layout()

figname = 'particles in a box small.png'
plt.savefig(figname, dpi=80, bbox_inches = 'tight')

plt.show()

#-------------------------------------------
#Plot the paths of particles on 30x30x30 grid, centre (0,0,0):
ax = plt.figure(figsize=(15,10),dpi=80).add_subplot(projection='3d')

for j in range(N):
    ax.plot(pos[:,j,0], pos[:,j,1], pos[:,j,2], '.', markersize = 1, label = f'particle{j}', alpha = 0.2)

#ax.legend()
ax.set_xlim([-15,15])
ax.set_ylim([-15,15])
ax.set_zlim([-15,15])
ax.set_title('Particles in a box')
ax.set_xlabel('$x$')
ax.set_ylabel('$y$')
ax.set_zlabel('$z$')

plt.tight_layout()

figname = 'particles in a box large.png'
plt.savefig(figname, dpi=80, bbox_inches = 'tight')

plt.show()