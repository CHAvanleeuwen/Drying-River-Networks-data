#!/bin/bash

## OAR instructions ##
#OAR -n albarine_dryver_obi3_speed
#OAR --project metabar
#OAR -l nodes=1/core=24,walltime=08:00:00
#OAR -O /bettik/lionnetc/dryver/logs/obi3_speed.%jobid%.stdout
#OAR -E /bettik/lionnetc/dryver/logs/obi3_speed.%jobid%.stderr

## define some bash options
set -e ## exit the script as soon as a function return an error
# set -x ## use to debug show all command

## load ciment environment and required modules
source /applis/site/nix.sh

if [ $# -lt 4 ]
then
  echo "Error, missing argument!"
fi

if [ $# -gt 4 ]
then
  echo "Error, too many argument!"
fi

LIBNAME=$1
PRIMER=$2
LIBDIR=$(realpath $3)
DBINFO=$(realpath $4)

# check to data integrity
echo 'Check the data integrity'

# check if the libraries file exist
if [ -z ${LIBNAME} ]
then
  echo "Error the libraries name is missing!"
  exit 1
fi

if [ -z ${PRIMER} ]
then
  echo "Error the primer name is missing!"
fi

# check if the primer information file exist
if [ ! -d ${LIBDIR} ]
then
  echo "Error the directory "${LIBDIR}" not exist!"
  exit 1
fi

# check if the primer information file exist
if [ ! -f ${DBINFO} ]
then
  echo "Error the primer information file "${DBINFO}" not exist!"
  exit 1
fi

# get the fasta file name of librairy
NBFILE=$(find ${LIBDIR} -name "*${LIBNAME}.filtered.uniq.fasta" | wc -l)
# check if the fasta file is uniq
if [ ${NBFILE} -gt 1 ]
then
  echo "Error, too many fasta file for the library "${LIBNAME}
  exit 1
fi

# process analysis
echo 'Process analysis'
echo `date`
echo "library: "${LIBNAME}
echo 'Primer: '$PRIMER

# get the minimal length of primer sequence
MINLENGTH=$(grep ${PRIMER,,} ${DBINFO} | cut -d',' -f2)
echo 'minimal length of primer: '${MINLENGTH}

if [ -z "${MINLENGTH}" ]
then
  echo "Error, the primer min length is missing!"
  exit 1
fi

# get the fasta file name
FASTAFILE=$(find ${LIBDIR} -name "*${LIBNAME}.filtered.uniq.fasta")
echo 'Fasta file used: '${FASTAFILE}

# generate the file name
FILENAME=${LIBNAME}"-"${PRIMER}
echo 'Prefix of output file name: '${FILENAME}

DMSPATH=${LIBDIR}/${LIBNAME}
# process analysis
obi import --preread --fasta-input ${FASTAFILE} ${DMSPATH}/${FILENAME}
obi annotate -k COUNT -k MERGED_sample --length ${DMSPATH}/${FILENAME} ${DMSPATH}/${FILENAME}_annot
obi grep -p "sequence['COUNT']>1 and sequence['seq_length']>${MINLENGTH}" ${DMSPATH}/${FILENAME}_annot ${DMSPATH}/${FILENAME}_denoise
obi clean --thread-count 24 -s MERGED_sample -r 0.5 ${DMSPATH}/${FILENAME}_denoise ${DMSPATH}/${FILENAME}_clean
obi grep -p "(sequence['obiclean_headcount']>0 or sequence['obiclean_singletoncount'])>0" ${DMSPATH}/${FILENAME}_clean ${DMSPATH}/${FILENAME}_clean_g
obi export --fasta-output -o ${LIBDIR}/${FILENAME}_grep_clean_g.fasta ${DMSPATH}/${FILENAME}_clean_g

set +e
source /applis/ciment/v2/env.bash
source /nfs_scratch/LECA_ENVIRONEMENT/env.bash
module load obitools
set -e

obiannotate -R "COUNT:count" ${LIBDIR}/${FILENAME}_grep_clean_g.fasta > $$; mv $$ ${LIBDIR}/${FILENAME}_grep_clean_g.fasta
obiannotate -R "MERGED_sample:merged_sample" ${LIBDIR}/${FILENAME}_grep_clean_g.fasta > $$; mv $$ ${LIBDIR}/${FILENAME}_grep_clean_g.fasta

echo `date`
## quit the script
exit $?
