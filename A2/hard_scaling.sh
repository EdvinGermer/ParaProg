echo "1 process"
mpirun --bind-to none -n 1 ./matmul ./input1800.txt output.txt
echo
echo "2 process"
mpirun --bind-to none -n 2 ./matmul ./input1800.txt output.txt
echo
echo "4 process"
mpirun --bind-to none -n 4 ./matmul ./input1800.txt output.txt
echo
echo "8 process"
mpirun --bind-to none -n 8 ./matmul ./input1800.txt output.txt
echo
