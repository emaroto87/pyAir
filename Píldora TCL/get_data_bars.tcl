# TCL macro to return several 1D element parameters of an element by its ID.
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
set input_data [open "D:/Jorge/Macros TCL/Pildora/prueba_bars.txt" r]
set element_data [read $input_data]
set output_data "D:/Jorge/Macros TCL/Pildora/propiedades_bars.txt"
set file_id [open $output_data w]


# Split the input into lines and take the first element. Then, retrieve the parameters and write it down in the output file
set panels [split $element_data "\n"]
set len [llength $panels]
for {set i 0} {$i < $len} {incr i} {
	# Iterate in i for each line in the input document
	set elms [lindex $panels $i]
	# Split the eids in the line by creating a list using , as separator
	set elm_list [split $elms ","]
	for {set j 0} {$j < [llength $elm_list]} {incr j} {
		set elm [lindex $elm_list $j]
		# Get the property ID
		set pid [hm_getvalue elems id=$elm dataname=property]
		# Check if it is a CBAR or a CBEAM and retrieve the parameters accordingly
		set type [hm_getvalue elems id=$elm dataname=typename]
		if {$type == "CBAR"} { 
			set beamid [hm_getvalue props id=$pid dataname=BeamSec]
			set dim1 [hm_getvalue beamsects id=$beamid dataname=beamsect_dim1]
			set dim2 [hm_getvalue beamsects id=$beamid dataname=beamsect_dim2]
			set dim3 [hm_getvalue beamsects id=$beamid dataname=beamsect_dim3]
			set dim4 [hm_getvalue beamsects id=$beamid dataname=beamsect_dim4]
			set type_sec [hm_getvalue beamsects id=$beamid dataname=standard_subtype]
			set Iyy_cog [hm_getvalue elems id=$elm dataname=IYYcog]
			set Izz_cog [hm_getvalue elems id=$elm dataname=IZZcog]
			# set Iyy [hm_getvalue beamsects id=$beamid dataname=results_Icentroid1]
			# set Izz [hm_getvalue beamsects id=$beamid dataname=results_Icentroid0]
			puts $file_id "$type $elm $pid $type_sec $dim1 $dim2 $dim3 $dim4 $Iyy_cog $Izz_cog"
			# puts $file_id "$Izz $Iyy"
		} else {
			set beamid [hm_getvalue props id=$pid dataname=BeamSecA]
			set dim1 [hm_getvalue props id=$pid dataname=pbeamlDIM1A]
			set dim2 [hm_getvalue props id=$pid dataname=pbeamlDIM2A]
			set dim3 [hm_getvalue props id=$pid dataname=pbeamlDIM3A]
			set dim4 [hm_getvalue props id=$pid dataname=pbeamlDIM4A]
			set type_sec [hm_getvalue props id=$pid dataname=pbeamlcstype]
			#set Iyy_cog [hm_getvalue elems id=$elm dataname=IYYcog]
			#set Izz_cog [hm_getvalue elems id=$elm dataname=IZZcog]
			set Iyy [hm_getvalue beamsects id=$beamid dataname=results_Icentroid1]
			set Izz [hm_getvalue beamsects id=$beamid dataname=results_Icentroid0]
			puts $file_id "$type $elm $pid $type_sec $dim1 $dim2 $dim3 $dim4 $Iyy_cog $Izz_cog"
			set cogy [hm_getvalue beamsects id=$beamid dataname=results_centroid1]
			set cogz [hm_getvalue beamsects id=$beamid dataname=results_centroid0]
			# puts $file_id "$elm $cogy $cogz"
		}
	}
	puts $file_id "\n"
}

close $file_id