# VMD text mode:  vmd -dispdev text -e select_fixed.tcl
# Run AFTER minimisation. Residues (whole) with >=1 atom within POCKET_CUTOFF of any ion stay mobile;
# all other protein atoms get beta=1 (frozen). Water + ions: beta=0 (free).
set cut $env(POCKET_CUTOFF)
mol new ionized.psf type psf waitfor all
mol addfile min.coor type namdbin waitfor all

set all [atomselect top all]
$all set beta 0
set prot [atomselect top protein]
$prot set beta 1
set pocket [atomselect top "same residue as (protein and within $cut of (segname =~ \"I0.*\"))"]
$pocket set beta 0
$all writepdb fixed.pdb

set f [open pocket_residues.txt w]
puts $f "# segid new_resid resname   (map back to PDB numbering with resmap.tsv)"
set seen {}
foreach r [$pocket get {segname resid resname}] {
    if {[lsearch -exact $seen $r] < 0} { lappend seen $r ; puts $f $r }
}
close $f
puts "POCKET_RESIDUES [llength $seen]  FIXED_ATOMS [$fixed num]  MOBILE_PROTEIN_ATOMS [$pocket num]"
exit
