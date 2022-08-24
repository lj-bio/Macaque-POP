path="/SCENIC/matrix"
cd $path
name=$(cat allsample.txt)
echo $name
for i in $name
do
{
	echo $i' start' `date`
	cd $path
	mkdir $i
	cd $i
	#grn
	pyscenic grn \
	--num_workers 20 \
	--output adj.${i}.tsv \
	--method grnboost2 \
	$path/${i}.csv \
	$path/allTFs.txt
	#ctx
	pyscenic ctx \
	adj.${i}.tsv \
	$path/hg19-tss-centered-10kb-7species.mc9nr.feather \
	--annotations_fname $path/motifs-v9-nr.hgnc-m0.001-o0.0.tbl \
	--expression_mtx_fname $path/${i}.csv \
	--mode "dask_multiprocessing" \
	--output reg_${i}.csv \
	--num_workers 20 
	#aucell
	pyscenic aucell \
	$path/${i}.csv \
	reg_${i}.csv \
	--output auc_${i}.csv \
	--num_workers 20
	echo $i' end' `date`
}
done
