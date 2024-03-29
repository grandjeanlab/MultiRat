root_dir='/project/4180000.19/multiRat/datadir'

cd $root_dir/data

output_dir=$root_dir'/export'
mkdir -p $output_dir


dataset_name='dataset01'

ds_type=01
ds_id=001
id=0


ls . | while read line
do
cd $line

sub=$ds_type$ds_id'0'$id

ses='1'


output_sub_dir=$output_dir'/sub-'$sub'/ses-'$ses
mkdir -p $output_sub_dir'/anat'
mkdir -p $output_sub_dir'/func'

anat_name='sub-'$sub'_ses-'$ses'_T2w.nii.gz'
func_name='sub-'$sub'_ses-'$ses'_run-1_bold.nii.gz'


Bru2 -a -z -o tmp 3/acqp
fslchpixdim tmp.nii.gz 0.11449226 0.11449226 0.114140127
fslorient -deleteorient tmp.nii.gz
fslswapdim tmp.nii.gz x -y z tmp2.nii.gz
3dresample -inset tmp2.nii.gz -prefix  $output_sub_dir'/anat/'$anat_name
rm tmp.nii.gz
rm tmp2.nii.gz

Bru2 -a -z -o tmp 8/acqp
fslorient -deleteorient tmp.nii.gz
fslswapdim tmp.nii.gz x z -y tmp2.nii.gz
3dresample -inset tmp2.nii.gz -prefix  $output_sub_dir'/func/'$func_name
rm tmp.nii.gz
rm tmp2.nii.gz




echo $dataset_name > $output_sub_dir/'info'
echo $line >> $output_sub_dir/'info'


id=$((id + 1))


cd ..
done
cd ../export
tree

