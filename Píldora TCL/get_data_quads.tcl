# TCL macro to return several parameters of each CQUAD by its ID. 
# https://2021.help.altair.com/2021/hwdesktop/hwd/topics/reference/hm/data_names-properties.htm


# Option to exit the script
set stop_for 0

bind . <Escape> {set stop_for 1}

proc min {args} {
	# Function to find the minimum in a list
    set min [lindex $args 0]
    foreach val $args {
        if {$val < $min} {
            set min $val
        }
    }
    return $min
}
proc max {args} {
	# Function to find the maximum in a list
    set max [lindex $args 0]
    foreach val $args {
        if {$val > $max} {
            set max $val
        }
    }
    return $max
}

puts "\n\n\n\n\n\n\n\n"

# Define the input and output documents
set input_data [open "D:/Jorge/Macros TCL/Pildora/prueba_quads.txt" r]
set element_data [read $input_data]
set output_data "D:/Jorge/Macros TCL/Pildora/propiedades_quads.txt"
set file_id [open $output_data w]


# Split the input into lines and take the first element. Then, retrieve its property ID and write it down in the output file
set panels [split $element_data "\n"]
set len [llength $panels]
puts $len
for {set i 0} {$i < $len} {incr i} {
	# Iterate in i for each line in the input document
	set elms [lindex $panels $i]
	# Split the eids in the line by creating a list using , as separator
	set elm_list [split $elms ","]
	for {set j 0} {$j < [llength $elm_list]} {incr j} {
		set elm [lindex $elm_list $j]
		# Get the property ID
		set pid [hm_getvalue elems id=$elm dataname=propertyid]
		set thickness [hm_getvalue props id=$pid dataname=thickness]
		set nodes [hm_getvalue elems id=$elm dataname=nodes]
		# Check if the element is a QUAD or a TRIA by counting the number of nodes it is composed of. Then, retrieve the parameters accordingly
		if {[llength $nodes] == 4} { 
			set y1 [hm_getvalue nodes id=[lindex $nodes 0] dataname=globaly]
			set y2 [hm_getvalue nodes id=[lindex $nodes 1] dataname=globaly]
			set y3 [hm_getvalue nodes id=[lindex $nodes 2] dataname=globaly]
			set y4 [hm_getvalue nodes id=[lindex $nodes 3] dataname=globaly]
			set z1 [hm_getvalue nodes id=[lindex $nodes 0] dataname=globalz]
			set z2 [hm_getvalue nodes id=[lindex $nodes 1] dataname=globalz]
			set z3 [hm_getvalue nodes id=[lindex $nodes 2] dataname=globalz]
			set z4 [hm_getvalue nodes id=[lindex $nodes 3] dataname=globalz]
			# Get the maximum and minimum y and z node-coordinates to define the width of the element
			set miny [min $y1 $y2 $y3 $y4]
			set minz [min $z1 $z2 $z3 $z4]
			set maxy [max $y1 $y2 $y3 $y4]
			set maxz [max $z1 $z2 $z3 $z4]
			# Get the x-position of the centroid and store it 
			set x_centroid [hm_getvalue elems id=$elm dataname=centerx]

		} else {
			set y1 [hm_getvalue nodes id=[lindex $nodes 0] dataname=globaly]
			set y2 [hm_getvalue nodes id=[lindex $nodes 1] dataname=globaly]
			set y3 [hm_getvalue nodes id=[lindex $nodes 2] dataname=globaly]
			set z1 [hm_getvalue nodes id=[lindex $nodes 0] dataname=globalz]
			set z2 [hm_getvalue nodes id=[lindex $nodes 1] dataname=globalz]
			set z3 [hm_getvalue nodes id=[lindex $nodes 2] dataname=globalz]
			# Get the maximum and minimum y and z node-coordinates to define the width of the element
			set miny [min $y1 $y2 $y3]
			set minz [min $z1 $z2 $z3]	
			set maxy [max $y1 $y2 $y3]
			set maxz [max $z1 $z2 $z3]	
			# Get the x-position of the centroid and store it 
			set x_centroid [hm_getvalue elems id=$elm dataname=centerx]
		}
		puts $file_id "$elm $pid $x_centroid $miny $maxy $minz $maxz $thickness"
	}
	puts $file_id "\n"
}



close $file_id