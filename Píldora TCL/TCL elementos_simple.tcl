set save_file "D:/Jorge/Macros TCL/Pildora/prueba_bars.txt"
set file_id [open $save_file w]

puts "\n\n\n\n\n\n\n\n"

proc update_counter {} {
	# Function to update the counter
	global fr_num
	global str_num
	incr fr_num
	set str_num 1
	puts "Next: FR$fr_num // STR$str_num"
}	

# Set a counter of stringers and frames
set str_num 1
set fr_num 1
set stop_for 0

bind . <Escape> {set stop_for 1}

# Writes in the beginning of the file
# Transform the file in a Tcl list where each line is an element
# -------------------------------
for {set i 0} {$i < 3} {incr i} {
	# Create mark in HM interface and store it
	*createmarkpanel elements 1
	set elems [hm_getmark elems 1]
	hm_markclearall 1
	# Write the element ids in the file separated by commas
	puts $file_id [join $elems ","]
	incr str_num
	puts "Next: FR$fr_num // STR$str_num"
	
	# Update counter if enter is pressed
	bind . <Return> {update_counter}
	if {$stop_for} {
		break
	}
}

close $file_id