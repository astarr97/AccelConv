for folder in run*;
do
	cd $folder
	sbatch -p hbfraser,hns $folder.sh
	cd ..
done
