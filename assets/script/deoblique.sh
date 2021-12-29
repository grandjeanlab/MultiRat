module load afni

cd /project/4180000.19/multiRat/to_trim

ls ./ | while read line
do
echo $line
cd ${line}

cd anat
nii=$(ls *.nii.gz)
mv $nii tmp.nii.gz
3dWarp -oblique2card -prefix $nii tmp.nii.gz
rm tmp.nii.gz

cd ../func
nii=$(ls *.nii.gz)
mv $nii tmp.nii.gz
3dWarp -oblique2card -prefix $nii tmp.nii.gz
rm tmp.nii.gz

cd ../..
done
