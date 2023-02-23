ls |cut -d "_" -f 1 |sort|uniq| xargs mkdir 
for dir in `ls -d */`; do mv ${dir/\/}_* $dir; done


# GSM6281320 GSM6281322 GSM6281321 GSM6281323 GSM6281324 GSM6281325 GSM6281326 GSM6281327

 

for dir in `ls -d */`; do
  cd $dir

  gunzip *gz
  mkdir -p outs/spatial
  mkdir -p outs/raw_feature_bc_matrix
  
  mv *spatial_* outs/spatial
  mv *raw_feature_bc_matrix_* outs/raw_feature_bc_matrix
  
  cd outs/spatial
  for file in ./*; do mv $file  ${file/*_spatial_/}; done
  
  cd  ../raw_feature_bc_matrix
  for file in ./*; do mv $file  ${file/*raw_feature_bc_matrix_/}; done
  
  cd ../../../
done
