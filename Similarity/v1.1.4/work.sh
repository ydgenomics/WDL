rds="~{rds}"
output_name="~{output_name}"
batch_key="~{batch_key}"
cluster_key="~{cluster_key}"
threshold_value=0.95
only_metaNeighbor="yes"  # Set to "yes" to skip Jaccard clustering

if [ "$only_metaNeighbor" != "yes" ]; then
    /opt/conda/bin/Rscript /WDL/Similarity/v1.1.4/jaccard_hclust.R \
    --input_file $rds --output_name $output_name --batch_key $batch_key --cluster_key $cluster_key
    rm Rplots.pdf
fi

/opt/conda/bin/Rscript /WDL/Similarity/v1.1.4/metaNeighbor.R \
--input_file $rds --output_name $output_name --batch_key $batch_key --cluster_key $cluster_key --threshold_value $threshold_value

path=$(find "$(pwd)" -maxdepth 1 -name '*_metaNeighbor.csv' -exec readlink -f {} \;)
path=$(echo "$path" | head -n 1)
echo "Path to celltype_NV of metaNeighbor output: $path"
seq='seq.txt'
slimit=0.95
python /WDL/Similarity/v1.1.4/sanky_plot.py \
--path $path --seq $seq --slimit $slimit