# -------------------------------------------------------------------
# Author:	Marina Cantador Tamurejo
# Mail:		marina.cantador.external@airbus.com
# Version:  V00_01
# -------------------------------------------------------------------
# TCL CODE WHICH GENERATES THE .DAT INPUT FILE FOR 'SECCIONES' TOOL.
# -------------------------------------------------------------------
#
# The code reads the ID of the nodes selected in Hypermesh screen and 
# generates .dat file with the file format requested by the tool.
#
# V00_01: The file returned is filled with the data despite an error
# occurs during the nodes selection
# Define de file where the data will be stored
# -------------------------------------------------------------------
# Tool for 8-character writting. Fill the gaps with spaces, for example
# "12345" turns into "   12345"
proc to_eight_char {string_to_transform} {
	set string_spaces "        $string_to_transform"
	return [string range $string_spaces [expr [string length $string_spaces]-8] [string length $string_spaces]]
}
# -------------------------------------------------------------------
hm_usermessage "Select dat file. If it does not exist, create a new one by writting the name follwed by .dat"
set nombre_dat [tk_getSaveFile -title "Select a File"]
puts $nombre_dat
set dat [open $nombre_dat "w"]
# Header of .dat
puts $dat "                        PROGRAMA SECCION (nudos)"
puts $dat "             Generado con tcl Generador_seccionesV00_01"
puts $dat "========================================================================"
puts $dat "TITULO:  MALE,  PROGRAMA SECCION 2014"
puts $dat "=======*================================================================"
puts $dat "       0        NUMERO DE SECCIONES A PROCESAR"
puts $dat "=======*================================================================"
# -------------------------------
set latest_node no
set section_number 0
while {$latest_node != yes} {
set section_number [expr $section_number+1]
# -----------------BLOCK 1-----------------------------
*createmarkpanel nodes 1 "Select the node which define the section"
set section_node [hm_getmark nodes 1]
hm_markclearall 1
*createmarkpanel nodes 1 "Select the nodes in one of the side of the section"
set side1_nodes [hm_getmark nodes 1]
hm_markclearall 1
*createmarkpanel nodes 1 "Select the nodes in the other side of the section"
set side2_nodes [hm_getmark nodes 1]
hm_markclearall 1
#
set node_number_1 [expr [llength $section_node]+[llength $side1_nodes]+[llength $side2_nodes]]
set block_1 [concat $section_node $side1_nodes $side2_nodes]
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
# -----------------BLOCK 2-----------------------------
*createmarkpanel nodes 1 "Select the node in Xdirection, in the zone to remove"
set x_node [hm_getmark nodes 1]
hm_markclearall 1
*createmarkpanel nodes 1 "Select the nodes in one of the side to remove"
set side1_removenodes [hm_getmark nodes 1]
hm_markclearall 1
*createmarkpanel nodes 1 "Select the nodes in the other side to remove"
set side2_removenodes [hm_getmark nodes 1]
hm_markclearall 1
#
set node_number_2 [expr [llength $x_node]+[llength $side1_removenodes]+[llength $side2_removenodes]]
set block_2 [concat $x_node $side1_removenodes $side2_removenodes]
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
*createmarkpanel nodes 1 "Select the node which define the Y direction"
set y_node [hm_getmark nodes 1]
hm_markclearall 1
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
}
close $dat