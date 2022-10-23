#!/bin/bash

# Declear vars for number of files and directories
num_of_files=0
num_of_dirs=0

# print the the instruction a screen
echo -e "Please enter a directory : \c"

# Read the input directory
read dir

# Move to the input diretory
cd $dir

# Loop over the contains in the input directory
for f in $(ls $dir); do

    # Check for file existence
    if [ -f $f ]; then
        num_of_files=$(($num_of_files+1))

    # Check for directory existence
    elif [ -d $f ]; then
        num_of_dirs=$(($num_of_dirs+1))
    fi

done

# Print the result on the screen
echo 'Number of files = ' $num_of_files
echo 'Number of directories = ' $num_of_dirs

        


