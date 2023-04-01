echo
echo %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
echo

echo "8.000.000 DATAPOINTS"
echo "Number of threads: 1"
mpirun --bind-to none -n 1 ./stencil /home/maya/public/PDP_Assignment1/input8000000.txt ./output.txt 1
echo "Number of threads: 2"
mpirun --bind-to none -n 2 ./stencil /home/maya/public/PDP_Assignment1/input8000000.txt ./output.txt 1
echo "Number of threads: 4"
mpirun --bind-to none -n 4 ./stencil /home/maya/public/PDP_Assignment1/input8000000.txt ./output.txt 1
echo "Number of threads: 8"
mpirun --bind-to none -n 8 ./stencil /home/maya/public/PDP_Assignment1/input8000000.txt ./output.txt 1
