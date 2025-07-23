mpiexec -n 8 ./microhh init eurec4a
#./microhh init eurec4a

find . -maxdepth 1 -type f -name '*_ext*' | while read -r file; do
    newname="${file/_ext/}"
    echo "Renaming: $file -> $newname"
    mv "$file" "$newname"
done

mpiexec -n 8 ./microhh run eurec4a
#./microhh run eurec4a

python cross_to_nc.py -p single -n 16
