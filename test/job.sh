#PBS -N test
#PBS -l nodes=1:ppn=12
#PBS -r n
#PBS -l walltime=252:0:0


cd /path/test

./full
