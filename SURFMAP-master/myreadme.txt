This tool gets a pdb and maps a parameter to a 2D map (it is possible to change the parameter - by changing the b-factor)


ubuntu
cd /mnt/c/Users/namih/Documents/BGU/SURFMAP-master

python3.10.exe surfmap.py -pdb test .pdb -tomap "bfactor" --keep

Rscript scripts/computeCoordList.R -f output_SURFMAP_test_bfactor/test_bfactor_partlist.out -s 5 -P sinusoidal

Rscript scripts/computeMatrices.R -i output_SURFMAP_test_bfactor/coord_lists/test_bfactor_coord_list.txt -s 5 -P sinusoidal

Rscript scripts/computeMaps.R -i output_SURFMAP_test_bfactor/smoothed_matrices/test_bfactor_smoothed_matrix.txt --bfactor -s 5  -P sinusoidal
