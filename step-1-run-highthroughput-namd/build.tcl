# VMD text mode:  vmd -dispdev text -e build.tcl
# Builds PSF (psfgen), solvates with PADDING A, writes cell.tcl for NAMD.
package require psfgen
package require solvate
set ff  $env(FF_DIR)
set pad $env(PADDING)

resetpsf
topology $ff/top_all36_prot.rtf
topology $ff/toppar_water_ions.str
topology $ff/top_all36_cgenff.rtf
foreach t [split $env(EXTRA_TOP) ":"] { if {$t ne ""} { topology $t } }
pdbalias atom ILE CD1 CD

set fp [open segments.txt r]
while {[gets $fp line] >= 0} {
    lassign $line type seg file extra
    if {$type eq "protein"} {
        segment $seg { first $extra ; last CTER ; pdb $file }
    } else {
        segment $seg { first NONE ; last NONE ; residue 1 $extra }
    }
    coordpdb $file $seg
}
close $fp

if {[file exists disu.txt]} {
    set fp [open disu.txt r]
    while {[gets $fp line] >= 0} { lassign $line s1 r1 s2 r2 ; patch DISU $s1:$r1 $s2:$r2 }
    close $fp
}
regenerate angles dihedrals
guesscoord                       ;# rebuilds missing heavy atoms + all hydrogens
writepsf prot_ions.psf
writepdb prot_ions.pdb

solvate prot_ions.psf prot_ions.pdb -t $pad -o solvated

# periodic cell from the solvated system -> cell.tcl (sourced by NAMD)
mol new solvated.psf type psf waitfor all
mol addfile solvated.pdb type pdb waitfor all
set all [atomselect top all]
lassign [measure minmax $all] mn mx
set f [open cell.tcl w]
foreach d {x y z} i {0 1 2} {
    set lo [lindex $mn $i] ; set hi [lindex $mx $i]
    puts $f "set cb$d [expr {$hi - $lo}]"
    puts $f "set co$d [expr {0.5*($hi + $lo)}]"
}
close $f
set q [vecsum [$all get charge]]
set f [open charge.txt w] ; puts $f $q ; close $f
puts "NET_CHARGE $q"

# ---------------------------------------------------------------- neutralise + salt
set fp [open salt.txt r] ; lassign [gets $fp] cat_res cat_atom ; close $fp
set q  [vecsum [$all get charge]]
set qi [expr {int(round($q))}]
set nw [[atomselect top "resname TIP3 and name OH2"] num]
set nsalt [expr {int(round($conc*$nw/55.5))}]
set ncat [expr {$nsalt + ($qi < 0 ? -$qi : 0)}]
set nan  [expr {$nsalt + ($qi > 0 ?  $qi : 0)}]
puts "IONIZE: net charge $q, waters $nw, salt pairs $nsalt -> $ncat $cat_res + $nan CLA"
 
set cand [atomselect top "resname TIP3 and name OH2 and not within $d_sol of (not resname TIP3)"]
set sites [pick_sites [shuffle [$cand get {segname resid x y z}]] [expr {$ncat+$nan}] $d_ion]
if {[llength $sites] < $ncat+$nan} { error "only [llength $sites] ion sites found, need [expr {$ncat+$nan}]" }
mol delete top
 
resetpsf
readpsf water_only.psf
coordpdb water_only.pdb
foreach s $sites { lassign $s wseg wres x y z ; delatom $wseg $wres }
set i 0
foreach s $sites {
    incr i
    lassign $s wseg wres x y z
    if {$i <= $ncat} {
        set seg [format C%03d $i] ; set res $cat_res ; set atom $cat_atom
    } else {
        set seg [format A%03d [expr {$i-$ncat}]] ; set res CLA ; set atom CLA
    }
    segment $seg { first NONE ; last NONE ; residue 1 $res }
    coord $seg 1 $atom [list $x $y $z]
    if {$i <= $ncat && $res in $nh4list} { add_nh4_h $seg 1 $x $y $z }
}
guesscoord
writepsf solvated.psf
writepdb solvated.pdb
 
# ---------------------------------------------------------------- report
mol new solvated.psf type psf waitfor all
mol addfile solvated.pdb type pdb waitfor all
set q [vecsum [[atomselect top all] get charge]]
set f [open charge.txt w] ; puts $f $q ; close $f
puts "NET_CHARGE_FINAL $q"


# sodiums are replaced by the chosen cation (salt.txt: NH4 / POT / SOD) -> ionized.{psf,pdb}
package require psfgen
package require autoionize

autoionize -psf solvated.psf -pdb solvated.pdb -sc 0.15 -cation SOD -anion CLA -o autoion

# ---- 3) replace them with the chosen cation
resetpsf
topology $ff/top_all36_prot.rtf
topology $ff/toppar_water_ions.str
topology $ff/top_all36_cgenff.rtf
foreach t [split $env(EXTRA_TOP) ":"] { if {$t ne ""} { topology $t } }
readpsf autoion.psf
coordpdb autoion.pdb
 set use_tetra 1
if {$cat_res ne "SOD"} {
    foreach s $sites { lassign $s sg rs x y z ; delatom $sg $rs }
    set i 0
    foreach s $sites {
        incr i
        lassign $s sg rs x y z
        set seg [format C%03d $i]
        segment $seg { first NONE ; last NONE ; residue 1 $cat_res }
        coord $seg 1 $cat_atom [list $x $y $z]
        if {$use_tetra && $cat_res in $nh4list} {
            set d [expr {$NH_DIST/sqrt(3.0)}] ; set n 0
            foreach v {{1 1 1} {1 -1 -1} {-1 1 -1} {-1 -1 1}} {
                incr n ; lassign $v a b c
                coord $seg 1 H$n [list [expr {$x+$d*$a}] [expr {$y+$d*$b}] [expr {$z+$d*$c}]]
            }
        }
    }
    guesscoord                   ;# places the missing NH4 hydrogens
}
writepsf ionized.psf
writepdb ionized.pdb
# ---- 4) report
mol new ionized.psf type psf waitfor all
mol addfile ionized.pdb type pdb waitfor all
set q [vecsum [[atomselect top all] get charge]]
set f [open charge.txt w] ; puts $f $q ; close $f
puts "NET_CHARGE_FINAL $q   ($cat_res: [llength $sites], CLA: $ncl)"

exit

