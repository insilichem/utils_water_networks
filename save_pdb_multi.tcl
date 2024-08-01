# Loop through the frames and save every frame
set num_frames [molinfo top get numframes]
for {set i 0} {$i < $num_frames} {incr i 2} {
    # Define the selection for each frame
    set sel [atomselect top "protein or (resname ROH 0GB 3GB COV or (same resid as ((((type OW and (within 3.5 of resid 495 and name OE1)) or (type OW and (within 3.5 of resid 495 and name OE2))) and (type OW and (within 4.5 of resid 289 and name C1))))))"]
    $sel frame $i
    set frame_num [format "%04d" $i]
    set outfile "protein_$frame_num.pdb"
    $sel writepdb $outfile
    puts "Saved frame $i to $outfile"

}
$sel delete

