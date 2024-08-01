import numpy as np
import matplotlib.pyplot as plt

plot = True
plt.rcParams['font.family'] = 'DejaVu Sans'

# Read data from File WT
data_file1 = np.loadtxt('./WT-cov-1000ns/r0/rmsf_data-wt-cov-r0.txt', skiprows=1)
residue1, rmsf1 = data_file1[:, 0], data_file1[:, 1]

data_file2 = np.loadtxt('./WT-cov-1000ns/r1/rmsf_data-wt-cov-r1.txt', skiprows=1)
residue2, rmsf2 = data_file2[:, 0], data_file2[:, 1]

data_file3 = np.loadtxt('./WT-cov-1000ns/r2/rmsf_data-wt-cov-r2.txt', skiprows=1)
residue3, rmsf3 = data_file3[:, 0], data_file3[:, 1]

#data_file4 = np.loadtxt('./WT-cov-1000ns/r3/rmsf_data-wt-cov-r3.txt', skiprows=1)
#residue4, rmsf4 = data_file4[:, 0], data_file4[:, 1]

#Read from File E220A
# Read data from File 1
edata_file1 = np.loadtxt('./E220A-cov-1000ns/r0/rmsf_data-e220a-cov-r0.txt', skiprows=1)
eresidue1, ermsf1 = edata_file1[:, 0], edata_file1[:, 1]

# Read data from File 2
edata_file2 = np.loadtxt('./E220A-cov-1000ns/r1/rmsf_data-e220a-cov-r1.txt', skiprows=1)
eresidue2, ermsf2 = edata_file2[:, 0], edata_file2[:, 1]

edata_file3 = np.loadtxt('./E220A-cov-1000ns/r2/rmsf_data-e220a-cov-r2.txt', skiprows=1)
eresidue3, ermsf3 = edata_file3[:, 0], edata_file3[:, 1]

# Read data from File 2
#edata_file4 = np.loadtxt('./E220A-cov-1000ns/r3/rmsf_data-e220a-cov-r3.txt', skiprows=1)
#eresidue4, ermsf4 = edata_file4[:, 0], edata_file4[:, 1]


#Correct the numeration

wt_res = [residue1, residue2, residue3]

for subarray in wt_res:
    for i in range(len(subarray)):
        subarray[i] = subarray[i]-4

e220a_res = [eresidue1, eresidue2, eresidue3]

for subarray in e220a_res:
    for i in range(len(subarray)):
        subarray[i] = subarray[i]-4


#Plotting time

cmap = plt.cm.viridis
ecmap = plt.cm.magma
dcmap = plt.cm.spring


#Prepare the WT average values
wt_val = [rmsf1, rmsf2, rmsf3]
wt_rmsf_avg = average_values = np.mean(wt_val, axis=0)

#Prepare the E220A average values
e220a_val = [ermsf1, ermsf2, ermsf3]
e220a_rmsf_avg = average_values = np.mean(e220a_val, axis=0)

#Delta RMSF
delta_rmsf = (np.array(e220a_rmsf_avg) - np.array(wt_rmsf_avg))


if plot:
    
    fig, (ax1, ax2) = plt.subplots(nrows=2, ncols=1, figsize=(15, 10), sharex=True, gridspec_kw={'height_ratios': [3, 1]})
    ax1.plot(residue1, wt_rmsf_avg, color=cmap(0.1), linewidth=0.9, label="WT") # avg WT
    ax1.plot(residue1, e220a_rmsf_avg, color=ecmap(0.8), linewidth=0.9, label="E220A") #avg E220A
    #ax.plot(residue1, rmsf1, color=cmap(0.1), linewidth=0.9, label="WT")
    #ax.plot(residue2, rmsf2, color=cmap(0.1), linewidth=0.9)
    #ax.plot(residue3, rmsf3, color=cmap(0.1), linewidth=0.9)
    #ax.plot(residue4, rmsf4, color=cmap(0.1), linewidth=0.9)
    #ax.plot(eresidue1, ermsf1, color=ecmap(0.8), linewidth=0.9, label="E220A")
    #ax.plot(eresidue2, ermsf2, color=ecmap(0.8), linewidth=0.9)
    #ax.plot(eresidue3, ermsf3, color=ecmap(0.8), linewidth=0.9)
    #ax.plot(eresidue4, ermsf4, color=ecmap(0.8), linewidth=0.9)
    ax1.set_xlabel('Residue number')
    ax1.set_ylabel('RMSF (Å)')
    ax1.set_title("WT & E220A cov-complex")
    ax1.set_xlim(0,580)
    ax1.set_ylim(0,4)
    ax1.legend()
    ax1.grid(True)
    
    ax2.plot(residue1, delta_rmsf, color=cmap(0.9), linewidth=0.9, label="ΔRMSF") #avg E220A
    ax2.set_xlabel('Residue number')
    ax2.set_ylabel('ΔRMSF (Å)')
    ax2.set_ylim(-1.5,1.5)
    ax2.legend()
    plt.xticks(np.arange(0, len(residue1)+1, step=20), rotation='vertical')
    
    #Write annotations of residues with RMSF greater than threshold
    annotations = []
    
    for x, y in zip(residue1, delta_rmsf):
        if abs(y) > 0.5:
            annotations.append(int(x))

    
    ax2.text(50, 1.2, "Residues with ΔRMSF > 0.5: "+str(annotations), fontsize=10, color='black')

    plt.tight_layout() 
    # Add grid lines
    #ax2.grid(True)
    plt.savefig("./RMSF-WT-E220A-3replicates-paper-quality.png")
