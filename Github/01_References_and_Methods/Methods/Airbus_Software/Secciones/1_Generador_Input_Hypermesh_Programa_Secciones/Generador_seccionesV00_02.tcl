# ===================================================================
# TCL CODE WHICH GENERATES THE .DAT INPUT FILE FOR 'SECCIONES' TOOL.
# ===================================================================
# Control of versions
# -------------------
# V00_01:
# 	Author:		Marina Cantador Tamurejo
# 	Mail:		marina.cantador.external@airbus.com
#   Description:
#		- First issue
# V00_02:
# 	Author:		Jorge Pérez Lozano
# 	Mail:		jorge.p.perez.external@airbus.com
# 	Description:
#		- Updated 19/05/2025 to write down the data by sides, maintaining the 
#		  order defined by the user by selecting the nodes individually. In addition,
#		  messages are also displayed in the console due to HyperMesh bugs not allowing
#		  dynamic updates of the GUI during *createmarkpanel procedure, and an instruction
#		  to exit the program using Escape key
#
# The code reads the ID of the nodes selected in Hypermesh screen and 
# generates .dat file with the file format requested by the tool.
#
# ===================================================================
#
# Tool for 8-character writing. Fill the gaps with spaces, for example
# "12345" turns into "   12345"
proc to_eight_char {string_to_transform} {
	set string_spaces "        $string_to_transform"
	return [string range $string_spaces [expr [string length $string_spaces]-8] [string length $string_spaces]]
}
# Stop the program if Escape key is pressed
bind . <Escape> {set stop_for 1}
# -------------------------------------------------------------------
# Define the file where the data will be stored
hm_usermessage "Select dat file. If it does not exist, create a new one by writing the name follwed by .dat"
set nombre_dat [tk_getSaveFile -title "Select a File"]
puts $nombre_dat
set dat [open $nombre_dat "w"]
# Header of .dat
puts $dat "                        PROGRAMA SECCION (nudos)"
puts $dat "             Generado con tcl Generador_seccionesV00_02"
puts $dat "========================================================================"
puts $dat "TITULO:  MALE,  PROGRAMA SECCION 2014"
puts $dat "=======*================================================================"
puts $dat "       0        NUMERO DE SECCIONES A PROCESAR"
puts $dat "=======*================================================================"
# -------------------------------
set latest_node no
set section_number 0
puts "\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n"
while {$latest_node != yes} {
set section_number [expr $section_number+1]
# -----------------BLOCK 1-----------------------------
puts "Select the node which defines the section"
*createmarkpanel nodes 1 "Select the node which defines the section"
set section_node [hm_getmark nodes 1]
hm_markclearall 1
puts $section_node
# Define the number of nodes in the section to loop over
set nodes_count_1_user [hm_getstring "Enter the number of nodes defining the section"]
set nodes_count_1 [expr $nodes_count_1_user - 1]
set side1_nodes {}
for {set i 0} {$i < $nodes_count_1} {incr i} {
	set counter [expr $i + 1]
	puts "Select the node $counter in one of the sides of the section"
	*createmarkpanel nodes 1 "Select the node $counter in one of the sides of the section"
	set selection_node_individual [hm_getmark nodes 1]
	hm_markclearall 1
	set side1_nodes [concat $side1_nodes $selection_node_individual]
	puts $selection_node_individual
	}

set node_number_1 [expr [llength $section_node]+[llength $side1_nodes]]
set block_1 [concat $section_node $side1_nodes]
#set node_number_1 [expr [llength $section_node]+[llength $side1_nodes]+[llength $side2_nodes]]
#set block_1 [concat $section_node $side1_nodes $side2_nodes]
# Write in the .dat file
puts $dat ""
puts $dat "                    DATOS DE LA SECCION  $section_number"
puts $dat ""
puts $dat "       0        NUMERO DE CASOS DE CARGA PARA ESTA SECCION"
puts $dat "-------*-----------------"
puts $dat "[to_eight_char $node_number_1]        NUMERO DE NUDOS DE LA SECCION Y LISTA"
foreach node $block_1 {
	puts $dat [to_eight_char $node]
}
set side2_nodes {}
# -----------------BLOCK 2-----------------------------
# Define the number of nodes in the section to loop over
set nodes_count_2 [hm_getstring "Enter the number of nodes in front of the section"]
for {set i 0} {$i < $nodes_count_2} {incr i} {
	set counter [expr $i + 1]
	puts "Select the node $counter in the other side of the section"
	*createmarkpanel nodes 1 "Select the node $counter in the other side of the section"
	set selection_node_individual [hm_getmark nodes 1]
	hm_markclearall 1
	set side2_nodes [concat $side2_nodes $selection_node_individual]
	puts $selection_node_individual
	}
#
set node_number_2 [expr [llength $side2_nodes]]
set block_2 [concat $side2_nodes]
# Write in the .dat file
puts $dat "     0.0        OFSET DESDE EL REVESTIMIENTO"
puts $dat "-------*-----------------"
puts $dat "[to_eight_char $node_number_2]        NUMERO DE NUDOS AL OTRO LADO DE LA SECCION "
foreach node $block_2 {
	if {$node == [lindex $block_2 0]} {
		puts $dat "[to_eight_char $node]        Y LISTA"
	} else {
		puts $dat [to_eight_char $node]}
}
# -----------------BLOCK 3-----------------------------
puts "Select the node in X-direction, in the zone to remove"
*createmarkpanel nodes 1 "Select the node in Xdirection, in the zone to remove"
set x_node [hm_getmark nodes 1]
hm_markclearall 1
puts $x_node
puts "Select the node in the XY-plane, which defines the Y-direction"
*createmarkpanel nodes 1 "Select the node in the XY-plane, which defines the Y-direction"
set y_node [hm_getmark nodes 1]
hm_markclearall 1
puts $y_node
# Write in the .dat file
puts $dat " -------*-----------------"
puts $dat "[to_eight_char $section_node]        ORIGEN     NUDOS PARA SISTEMA DE COORDENADAS"
puts $dat "[to_eight_char $x_node]        EJE  X  LOCAL (EN LA DIRECCION DE LA CUADERNA)"
puts $dat "[to_eight_char $y_node]        NUDO EN  XY  CUADRANTE X+,Y+"
puts $dat "=======*================================================================"
set latest_node [hm_getstring "Latest nodes?" "Write yes if it is the latest section defined. If not, keep it blank"]
close $dat
# Writes the number of sections processed
set in [open $nombre_dat "r+"]
set data [read $in]
regsub "[to_eight_char [expr $section_number-1]]        NUMERO DE SECCIONES A PROCESAR" $data "[to_eight_char $section_number]        NUMERO DE SECCIONES A PROCESAR" data
seek $in 0
puts -nonewline $in $data
close $in
set dat [open $nombre_dat "r+"]
seek $dat 0 end 
if {$stop_for} {
			break
		}
}
close $dat