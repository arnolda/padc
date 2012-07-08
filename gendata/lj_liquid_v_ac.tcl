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
set density [lindex $argv 1]
set samples [lindex $argv 2]

setmd time_step 0.01
setmd skin      0.4

thermostat langevin 1.0 1.0
set tcl_precision 6

#############################################################
#  Setup System                                             #
#############################################################

cellsystem nsquare

set box_l [expr pow($n_part/$density, 1./3)]

puts "box is $box_l"

setmd box_l $box_l $box_l $box_l

inter 0 0 lennard-jones 1 1 2.5 auto

for {set i 0} {$i < $n_part} {incr i} {
    part $i pos [expr rand()*$box_l] [expr rand()*$box_l] [expr rand()*$box_l] type 0
}

#############################################################
#      Warmup                                               #
#############################################################

puts "warmup"
for {set i 10} {$i < 100} {incr i 10} {
    inter ljforcecap $i
    integrate 100
}
inter ljforcecap 0
puts "T=[expr [analyze energy kinetic]/$n_part/1.5]"

#############################################################
#      Integration                                          #
#############################################################

puts "integration"

set out [open "v_$density.data" "w"]
puts $out "# time v"
for {set i 0} {$i < $samples/$n_part} { incr i} {
    integrate 1
    puts -nonewline $out "[setmd time]"
    for {set p 0} {$p < $n_part} {incr p} {
	puts -nonewline $out " [part $p pr v]"
    }
    puts $out ""
}

close $out
