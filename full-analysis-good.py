import pytraj as pt
import numpy as np
import math
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA
import time

#TRAJECTORY_FILENAME = "/media/xavi/99e463f2-6cec-40d5-85bb-b6aa84b4b648/hrmova_waters/cov-complex-WT/wt-cov-complex/rep3-cov-wt/24312/wt-cov-complex-1000ns-rep3-nowat.dcd"
#TOPOLOGY_FILENAME = "/media/xavi/99e463f2-6cec-40d5-85bb-b6aa84b4b648/hrmova_waters/cov-complex-WT/wt-cov-complex/rep0/wt-cov-complex-nowat.prmtop"

TRAJECTORY_FILENAME = "/media/xavi/99e463f2-6cec-40d5-85bb-b6aa84b4b648/hrmova_waters/cov-complex-e220a/rep3-cov-e220a/24313/e220a-cov-complex-1000ns-rep3-nowat.dcd"
TOPOLOGY_FILENAME = "/media/xavi/99e463f2-6cec-40d5-85bb-b6aa84b4b648/hrmova_waters/cov-complex-e220a/rep0/e220a-cov-complex-nowat.prmtop"

TRAJECTORY_LENGTH = 1000 #Length of the trajectory in ns (for correct axis labels in plots)

print(TRAJECTORY_FILENAME)

trajectory = pt.iterload(TRAJECTORY_FILENAME, TOPOLOGY_FILENAME, stride=2)
trajectory = trajectory[:]
print(trajectory)
time.sleep(3)
print("I sleep...")

SELECTION = ':10-590@CA'

trajectory[SELECTION]

distances = pt.analysis.rmsd.rmsd(traj=trajectory, mask=SELECTION, ref=0, update_coordinate=False)

RMSD_PLOT_FILENAME = '/media/xavi/99e463f2-6cec-40d5-85bb-b6aa84b4b648/hrmova_waters/e220a-cov-r3-RMSD.pdf' #Output filename for the plot

RMSD_COLOR = 'forestgreen' #Color of the line
MAX_RMSD = math.ceil(distances.max()) #It can be modified if you want to change the scale of the plot

plt.figure()
plt.plot([TRAJECTORY_LENGTH*a/len(distances) for a in range(len(distances))], distances, linewidth=0.2, color=RMSD_COLOR)
plt.ylim(0, MAX_RMSD)
plt.xlabel('Time (ns)')
plt.ylabel('RMSD (Å)')
plt.savefig(RMSD_PLOT_FILENAME, bbox_inches='tight')

#Save data into file for later plotting
data_rmsd = np.column_stack((np.arange(len(distances)) * TRAJECTORY_LENGTH / len(distances), distances))
np.savetxt('rmsd_data-e220a-cov-r3.txt', data_rmsd, header='Time (ns) RMSD (Å)', comments='')




"""fig, ax = plt.subplots()
ax.grid(True)
ax.figure()
ax.plot([TRAJECTORY_LENGTH*a/len(distances) for a in range(len(distances))], distances, linewidth=0.2, color=RMSD_COLOR)
ax.ylim(0, MAX_RMSD)
ax.xlabel('Time (ns)')
ax.ylabel('RMSD (Å)')
plt.savefig(RMSD_PLOT_FILENAME, bbox_inches='tight')"""
MIN_RMSD_VALUE = 0.0 #It can be modified if you want to force the color scale of the plot to a predefined RMSD range
MAX_RMSD_VALUE = math.ceil(distances.max()) #It can be modified if you want to force the color scale of the plot to a predefined RMSD range
ALL_FILENAME = '/media/xavi/99e463f2-6cec-40d5-85bb-b6aa84b4b648/hrmova_waters/All_rmsd-e220a-cov-r3.pdf' #Output filename for the plot

# Calculate all-to-all RMSD matrix (it could take some time if you choose a small step size!)
step_all = math.ceil(float(len(trajectory)/1000)) # Step size: trajectory frames that will be used in the analysis (e.g. if step is 20, frames 0, 20, 40 ... will be selected)
traj_all = trajectory[::step_all]
#traj_all = traj_all.atom_slice(SELECTION)
distances = np.empty((math.ceil(traj_all.n_frames), math.ceil(traj_all.n_frames)))
for i, frame in enumerate(range(0, traj_all.n_frames)):
    distances[i] = pt.analysis.rmsd.rmsd(traj=traj_all, mask=SELECTION, ref=i, update_coordinate=False)
    
fig, ax = plt.subplots()
im = ax.imshow(distances, origin='lower', vmin=MIN_RMSD_VALUE, vmax=MAX_RMSD_VALUE)
cbar = ax.figure.colorbar(im, ax=ax)
cbar.set_label('RMSD (Å)', rotation=-90, va="bottom")

ticks = np.arange(0, len(distances), len(distances)//5)
plt.xticks(ticks, [math.ceil(a*TRAJECTORY_LENGTH/len(distances)) for a in ticks])
plt.yticks(ticks, [math.ceil(a*TRAJECTORY_LENGTH/len(distances)) for a in ticks])
ax.set_xlabel('Time (ns)')
ax.set_ylabel('Time (ns)')

fig.tight_layout()
plt.savefig(ALL_FILENAME, bbox_inches='tight')
plt.show()
    
"""ALL_FILENAME = '/media/xavi/99e463f2-6cec-40d5-85bb-b6aa84b4b648/hrmova_waters/All_rmsd-e220a-cov-r3.pdf' #Output filename for the plot

MIN_RMSD_VALUE = 0.0 #It can be modified if you want to force the color scale of the plot to a predefined RMSD range
MAX_RMSD_VALUE = math.ceil(distances.max()) #It can be modified if you want to force the color scale of the plot to a predefined RMSD range
"""
step_PCA = math.ceil(float(len(trajectory)/10000)) # SI VOLS QUE AGAFI MENYS, FES MES GRAN EL NUMERO ¿O NO? Step size: trajectory frames that will be used in the analysis (e.g. if step is 20, frames 0, 20, 40 ... will be selected)
traj_PCA = trajectory[::step_PCA]
data = pt.pca(traj_PCA, mask=SELECTION, n_vecs=2)
projection_data = data[0]

PCA_FILENAME = '/media/xavi/99e463f2-6cec-40d5-85bb-b6aa84b4b648/hrmova_waters/PCA_e220a-cov-r3_full.pdf' #Output filename for the plot


MIN_X_VALUE = math.floor(projection_data[0].min(axis=0)) #It can be modified if you want to force the scale of the plot
MAX_X_VALUE = math.ceil(projection_data[0].max(axis=0)) #It can be modified if you want to force the scale of the plot
MIN_Y_VALUE = math.floor(projection_data[1].min(axis=0)) #It can be modified if you want to force the scale of the plot
MAX_Y_VALUE = math.ceil(projection_data[1].max(axis=0)) #It can be modified if you want to force the scale of the plot
POINTS_SIZE = 3 #You can adjust the size of the points of the plot

plt.figure()
plt.scatter(projection_data[0], projection_data[1], marker='.', c=range(len(traj_PCA)), vmin=0, vmax=len(traj_PCA), s=POINTS_SIZE)
plt.xlabel('PC1')
plt.ylabel('PC2')
plt.xlim([MIN_X_VALUE, MAX_X_VALUE])
plt.ylim([MIN_Y_VALUE, MAX_Y_VALUE])
cbar = plt.colorbar(aspect=30)
ticks = np.arange(0, len(traj_PCA), len(traj_PCA)//10)
cbar.set_ticks(ticks)
cbar.ax.tick_params()
cbar.set_ticklabels([math.ceil(a*TRAJECTORY_LENGTH/len(traj_PCA)) for a in ticks])
cbar.set_label('Time (ns)', rotation=-90, va="bottom")
plt.savefig(PCA_FILENAME, bbox_inches='tight')

pt.superpose(trajectory, ref=0, mask=SELECTION)

# compute rmsf
rmsf_data = pt.rmsf(trajectory, mask=SELECTION, options='byres')

RMSF_PLOT_FILENAME = '/media/xavi/99e463f2-6cec-40d5-85bb-b6aa84b4b648/hrmova_waters/RMSF_e220a-cov-r3_full.pdf' #Output filename for the plot
RMSF_COLOR = 'forestgreen' #Color of the line

plt.figure()
plt.bar(rmsf_data.T[0], rmsf_data.T[1], color=RMSF_COLOR)
plt.xlabel('Residue')
plt.ylabel('RMSF (Angstrom)')

plt.savefig(RMSF_PLOT_FILENAME, bbox_inches='tight')

#data_rmsf = np.column_stack((np.aange(len(distances)) * TRAJECTORY_LENGTH / len(distances), distances))
np.savetxt('rmsf_data-e220a-cov-r3.txt', rmsf_data, header='Residue RMSF (Å)', comments='')

