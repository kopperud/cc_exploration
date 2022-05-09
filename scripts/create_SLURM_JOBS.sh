#!/bin/bash

N_CORES=4
JOB_DIR="jobs"
LOG_DIR="logs"

## Delete old jobs
if [ ${JOB_DIR} != "" ]; then
  if [ ! -d ${JOB_DIR} ]; then
    mkdir ${JOB_DIR}
  else
    rm -f ${JOB_DIR}/*
  fi
fi

## Delete old logs
if [ ${LOG_DIR} != "" ]; then
  if [ ! -d ${LOG_DIR} ]; then
    mkdir ${LOG_DIR}
  else
    rm -f ${LOG_DIR}/*
  fi
fi


# need to find a way how to loop over all trees
for treepath in trees/*.tre;
do
    ds=`basename $treepath .tre`
    echo "${ds}"

    echo "#!/bin/bash
#SBATCH --job-name=HSMRF_${ds}
#SBATCH --output=${LOG_DIR}/HSMRF_${ds}.log
#SBATCH --error=${LOG_DIR}/HSMRF_${ds}.err
#SBATCH --ntasks=${N_CORES}
#SBATCH --nodes=1
#SBATCH --mem=${N_CORES}G
#SBATCH --qos=low
#
#SBATCH --mail-user b.kopperud@lmu.de
#SBATCH --mail-type=FAIL

mpirun -np ${N_CORES} rb-mpi-dev scripts/analysis.Rev --args ${ds} ${N_CORES} 1234 > ${LOG_DIR}/HSMRF_${ds}.out" > ${JOB_DIR}/${ds}.sh
#sbatch ${JOB_DIR}/${ds}.sh

done

echo "done ..."
