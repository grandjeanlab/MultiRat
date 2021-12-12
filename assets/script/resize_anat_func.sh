cd /project/4180000.19/multiRat/tmp_crop/ds50

ls ./ | while read line
do
echo $line
cd ${line}

ls ./ | while read ses
do
#ses='ses-1'
cd ${ses}/anat
nii=$(ls *.nii.gz)
mv $nii tmp.nii.gz
fslroi tmp.nii.gz $nii 40 170 0 -1 70 145
slicer $nii -a ../../../../$(remove_ext $nii).png

cd ../func
nii=$(ls *.nii.gz)
mv $nii tmp.nii.gz
fslroi tmp.nii.gz $nii 20 60 0 -1 27 53
fslmaths $nii -Tmean tmean
slicer tmean -a ../../../../$(remove_ext $nii).png
rm tmp.nii.gz
rm tmean.nii.gz

cd ../../
done
cd ../
done