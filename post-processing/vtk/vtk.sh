
for i in {1..10}
do
FILES1=field_${i}.dat
FILES1a=field_${i}.vtk

for f in $FILES1
do 
    echo "Processing $f file..."
    cat vtkhead $FILES1 > $FILES1a

    
done

done
echo "Job Done!"
