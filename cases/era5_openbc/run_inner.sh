# Cleanup!
rm test/dom1/*00*
rm test/dom1/*.nc

base_dir=$(pwd) 
nproc=8

python era5_openbc_input.py --domain=1

cd test/dom1

mpiexec -n $nproc ./microhh init era5_openbc

find . -maxdepth 1 -type f -name '*_overwrite*' | while read -r file; do
    newname="${file/_overwrite/}"
    echo "Renaming: $file -> $newname"
    mv "$file" "$newname"
done

mpiexec -n $nproc ./microhh run era5_openbc

python cross_to_nc.py -n 12

cd $base_dir
