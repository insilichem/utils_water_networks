import numpy as np
import matplotlib.pyplot as plt

plot = True

#Get WT data
data_file1 = np.loadtxt("../../DCD/WT/r0/r0-wt-RDP.dat", delimiter=" ")
residue1, rdf1 = data_file1[:, 0], data_file1[:, 1]

data_file2 = np.loadtxt("../../DCD/WT/r1/r1-wt-RDP.dat", delimiter=" ")
residue2, rdf2 = data_file2[:, 0], data_file2[:, 1]

data_file3 = np.loadtxt("../../DCD/WT/r2/r2-wt-RDP.dat", delimiter=" ")
residue3, rdf3 = data_file3[:, 0], data_file3[:, 1]

#data_file4 = np.loadtxt("./WT/wt-cov-r3-gofr.dat", delimiter=" ")
#residue4, rdf4 = data_file4[:, 0], data_file4[:, 1]

#Get E220A data
edata_file1 = np.loadtxt("../../DCD/E220A/r0/r0-e220a-RDP.dat", delimiter=" ")
eresidue1, erdf1 = edata_file1[:, 0], edata_file1[:, 1]

edata_file2 = np.loadtxt("../../DCD/E220A/r1/r1-e220a-RDP.dat", delimiter=" ")
eresidue2, erdf2 = edata_file2[:, 0], edata_file2[:, 1]

edata_file3 = np.loadtxt("../../DCD/E220A/r2/r2-e220a-RDP.dat", delimiter=" ")
eresidue3, erdf3 = edata_file3[:, 0], edata_file3[:, 1]

#edata_file4 = np.loadtxt("./E220A/e220a-cov-r3-gofr.dat", delimiter=" ")
#eresidue4, erdf4 = edata_file4[:, 0], edata_file4[:, 1]

#Correct the numeration
"""wt_res = [residue1, residue2, residue3, residue4]

for subarray in wt_res:
    for i in range(len(subarray)):
        subarray[i] = subarray[i]-4

e220a_res = [eresidue1, eresidue2, eresidue3, eresidue4]

for subarray in e220a_res:
    for i in range(len(subarray)):
        subarray[i] = subarray[i]-4"""

#Plotting time

cmap = plt.cm.viridis
ecmap = plt.cm.magma
dcmap = plt.cm.spring


#Prepare the WT average values
wt_val = [rdf1, rdf2, rdf3]
wt_rdf_avg = average_values = np.mean(wt_val, axis=0)

#Prepare the E220A average values
e220a_val = [erdf1, erdf2, erdf3]
e220a_rdf_avg = average_values = np.mean(e220a_val, axis=0)

if plot:
    
    fig, ax1 = plt.subplots(nrows=1, ncols=1, figsize=(10, 7))
    ax1.plot(residue1, wt_rdf_avg, color=cmap(0.1), linewidth=0.9, label="WT Average") # avg WT
    ax1.plot(residue1, e220a_rdf_avg, color=ecmap(0.8), linewidth=0.9, label="E220A Average") #avg E220A
    #ax.plot(residue1, rmsf1, color=cmap(0.1), linewidth=0.9, label="WT")
    #ax.plot(residue2, rmsf2, color=cmap(0.1), linewidth=0.9)
    #ax.plot(residue3, rmsf3, color=cmap(0.1), linewidth=0.9)
    #ax.plot(residue4, rmsf4, color=cmap(0.1), linewidth=0.9)
    #ax.plot(eresidue1, ermsf1, color=ecmap(0.8), linewidth=0.9, label="E220A")
    #ax.plot(eresidue2, ermsf2, color=ecmap(0.8), linewidth=0.9)
    #ax.plot(eresidue3, ermsf3, color=ecmap(0.8), linewidth=0.9)
    #ax.plot(eresidue4, ermsf4, color=ecmap(0.8), linewidth=0.9)
    ax1.set_xlabel('r (Å)')
    ax1.set_ylabel('g(r)')
    ax1.set_title("WT & E220A cov-complex Radial Pair Distribution")
    ax1.set_xlim(2,17.5)
    ax1.legend()
    ax1.grid(True)
    

    
    #Write annotations of residues with RMSF greater than threshold
    

    plt.tight_layout() 
    # Add grid lines
    #ax2.grid(True)
    plt.savefig("./RDF-WT-E220A-3reps-avg.png")
