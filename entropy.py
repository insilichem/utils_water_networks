import functions
from rdkit import Chem
from rdkit.Chem import rdmolfiles
import os
from collections import defaultdict

"""
Whole script idea: https://stackoverflow.com/questions/35709562/how-to-calculate-clustering-entropy-a-working-example-or-software-code%E2%80%8B

Script to calculate the entropy of an array of clusters. 

by XFF
"""

# assign directories
#directory = "/media/xavi/99e463f2-6cec-40d5-85bb-b6aa84b4b648/hrmova_waters/NEW-COV/DCD/WT/r0/r0-entropy-water"  # PUT YOUR TEST FILES HERE
#ref_directory = "./wt-ref-r0"  # PUT YOUR XRAY STRUCT HERE. ALIGN PREVIOUSLY WITH TEST FILES.

directory = "/media/xavi/99e463f2-6cec-40d5-85bb-b6aa84b4b648/hrmova_waters/NEW-COV/DCD/WT/r1/r1-entropy-water"  # PUT YOUR TEST FILES HERE
ref_directory = "./wt-ref-r1"

#directory = "/media/xavi/99e463f2-6cec-40d5-85bb-b6aa84b4b648/hrmova_waters/NEW-COV/DCD/WT/r2/r2-entropy-water"  # PUT YOUR TEST FILES HERE
#ref_directory = "./wt-ref-r2"

####################################################

#directory = "/media/xavi/99e463f2-6cec-40d5-85bb-b6aa84b4b648/hrmova_waters/NEW-COV/DCD/E220A/r0/r0-entropy-water"  # PUT YOUR TEST FILES HERE
#ref_directory = "./e220a-ref-r0"

#directory = "/media/xavi/99e463f2-6cec-40d5-85bb-b6aa84b4b648/hrmova_waters/NEW-COV/DCD/E220A/r1/r1-entropy-water"  # PUT YOUR TEST FILES HERE
#ref_directory = "./e220a-ref-r1"

#directory = "/media/xavi/99e463f2-6cec-40d5-85bb-b6aa84b4b648/hrmova_waters/NEW-COV/DCD/E220A/r2/r2-entropy-water"  # PUT YOUR TEST FILES HERE
#ref_directory = "./e220a-ref-r2"

######################################################33


threshold = 2.7  # CHANGE THE THRESHOLD VALUE
total_entropy = []
N_total = 0
xray_positions = functions.get_xray_positions(ref_directory)

watlim = 50000  #CHANGE LIMIT OF WATER MOLECULES IN CHANNEL
channel_pdbs = []
# iterate over files in
# that directory
print("Loading water channels in pdb format")

for filename in os.listdir(directory):

    f = os.path.join(directory, filename)
    # checking if it is a file
    if os.path.isfile(f):

        channel_pdbs.append(
            rdmolfiles.MolFromPDBFile(
                f,
                sanitize=True,
                removeHs=False,
                proximityBonding=False,
            )
        )


print("Loaded ", len(channel_pdbs), " PDBs")
print("Begin cluster entropy calculation")

k=0
for mol in channel_pdbs:
    k +=1

    if mol == None:
        del mol
        break
    atomic_count = defaultdict(lambda: 0)
    try:
        for at in mol.GetAtoms():

            atomic_count[at.GetAtomicNum()] += 1
    except:

        pass


    if int(atomic_count[8]) <= watlim:
        #print("????")
        # Define arrays for X and O points within channel
        o_ = []
        x_ = []
        for i, at in enumerate(mol.GetAtoms()):  # for each atom in the given molecule


            position = mol.GetConformer().GetAtomPosition(
                i
            )  # calculate its position and
            for xyz in xray_positions:  # given the xray positions
                #print("new cycle\n")
                val = functions.compute_geom_distance(
                    position, xyz
                )  
                #print(val)# check the distance between all of them
                if val <= threshold:  # append 1 to O if the condition is met
                    o_.append(1)
                    #print("appending to O_")
                    break
                elif val > threshold:  # otherwise continue
                    #print("passing")
                    pass

            else:  # and if it fails to append to O, append to X
                #print("appending to X_")
                x_.append(1)

        N_w = len(x_) + len(o_)  # define total number of points in channel
        
        total_entropy.append(
            (N_w, functions.calculate_entropy_channel(N_w, len(x_), len(o_)))
        )
        #print(N_w, len(x_), len(o_))
        

    else:

        pass
print("Done.\nCalculating entropy of system")
# calculate total number of channels. Convert tuple into dict for easier iteration

for i, j in total_entropy:
    N_total += int(i)

print("Number of PDB analysed: ", len(total_entropy))
result = functions.total_entropy(N_total, total_entropy)
print("Done.\nTotal entropy of system: ", result)
zeroes = []
print("Classification of zero entropy clusters\n")
for i, j in total_entropy:
    if j == -0.0 or j == 0.0:
        zeroes.append(j)
print("len zeroes: ", len(zeroes))

"""Pseudocode:
channel_entropies = []

Get reference PDB (previously aligned with the channels)

for each xyz in reference_struct:
    store position in xray_positions_array
    
For each PDB in data folder:
    channel = Channel([],0) #create channel object to calculate entropy
    convert pdb in rdkit mol
    for at in mol.GetAtoms():
        for i in xray_positions_array:
            if distance(at, i) <= thresh:
                channel.append("O")
            else:
                channel.append("X")
    channel_entropies.append(calculate_entropy_of_channel(channel)      
    """
