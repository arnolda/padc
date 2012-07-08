# Dies ist Teil der Vorlesung Physik auf dem Computer, SS 2012,
# Axel Arnold, Universitaet Stuttgart.
# 
# Dieses Werk ist unter einer Creative Commons-Lizenz vom Typ
# Namensnennung-Weitergabe unter gleichen Bedingungen 3.0 Deutschland
# zugänglich. Um eine Kopie dieser Lizenz einzusehen, konsultieren Sie
# http://creativecommons.org/licenses/by-sa/3.0/de/ oder wenden Sie sich
# schriftlich an Creative Commons, 444 Castro Street, Suite 900, Mountain
# View, California, 94041, USA.
#############################################################
#  Parameters                                               #
#############################################################

set n_part  [lindex $argv 0]
set box_l   [lindex $argv 1]
set epss    [lindex $argv 2]
set limit   [lindex $argv 3]

if {"$limit" == ""} { set limit 0 }

setmd time_step 0.005
setmd skin      0.4

set tcl_precision 6

#############################################################
#  Setup System                                             #
#############################################################

set density [expr $n_part/pow($box_l, 2)]
puts "box is $box_l, density is $density"

setmd box_l $box_l $box_l $box_l

inter 0 0 lennard-jones 1 1 [expr pow(2,1./6.)] auto

#############################################################
# Main loop                                                 #
#############################################################

set coulomb_params ""

foreach eps $epss {
    puts "epsilon = $eps"

#############################################################
#      Particles and Warmup                                 #
#############################################################
    
    for {set i 0} {$i < $n_part} {incr i} {
	set q [expr 2*($i < $n_part/2) - 1]
	part $i pos [expr rand()*$box_l] [expr rand()*$box_l] 0 \
	    q $q type 0 fix 0 0 1
    }

    thermostat langevin 1.0 1.0

    puts "Warmup"
    if {$coulomb_params == ""} {
	inter coulomb 4 p3m tune accuracy 1e-4 mesh 16
	set coulomb_params [lindex [inter coulomb] 0]
    } {
	eval "inter $coulomb_params"
    }
    for {set i 1} {$i < 1000} {incr i 100} {
	inter ljforcecap $i
	integrate 1000
    }
    inter ljforcecap 0

    puts "At end of warmum T=[expr [analyze energy kinetic]/$n_part]"
    puts "potential energy is [analyze energy nonbonded 0 0]"

#############################################################
#      Integration                                          #
#############################################################

    setmd time 0

    set out [open "lj_E_$eps.data" "w"]
    puts $out "# time Epot"
    set T 1.0
    while {$T > 0.001 && ($limit == 0 || [setmd time] < $limit)} {
	thermostat langevin $T 1.0

	eval "inter [lreplace $coulomb_params 1 1 [expr 1.0/$T]]"

	set E_pot [expr [analyze energy total] - [analyze energy kinetic]]
	puts "energy is [analyze energy]"
	puts $out "[setmd time] $E_pot"
	integrate 100
	set T [expr (1-$eps)*$T]
	puts "Ttgt=$T, T=[expr [analyze energy kinetic]/$n_part]"
	flush $out
    }
    close $out

    set f [open "lj_snapshot_A_$eps.coord" "w"]
    for {set i 0} {$i < $n_part/2} {incr i} {
	puts $f [lrange [part $i pr folded] 0 1]
    }
    close $f
    set f [open "lj_snapshot_B_$eps.coord" "w"]
    for {set i [expr $n_part/2]} {$i < $n_part} {incr i} {
	puts $f [lrange [part $i pr folded] 0 1]
    }
    close $f
}
