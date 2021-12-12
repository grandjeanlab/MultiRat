# code to substitute code in job files.

filename='job_sub-0101007_ses-1.sh'
search='N4_reg --anat_autobox'
replace='no_reg'
sed -i "s/$search/$replace/" $filename