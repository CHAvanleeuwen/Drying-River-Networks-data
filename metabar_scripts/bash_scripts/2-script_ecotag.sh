#!/bin/bash

## OAR instructions ##
#OAR -n albarine_dryver_ecotag
#OAR --project metabar
#OAR -l nodes=1/core=12,walltime=12:00:00
#OAR -O /bettik/lionnetc/dryver/logs/ecotag.%jobid%.stdout
#OAR -E /bettik/lionnetc/dryver/logs/ecotag.%jobid%.stderr

set +e 

## load ciment environment and required modules
source /applis/site/nix.sh
source /applis/ciment/v2/env.bash
source /nfs_scratch/LECA_ENVIRONEMENT/env.bash
module load obitools
set -e

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

echo 'Check the data integrity'
if [ -z ${LIBNAME} ]
then
  echo "Error the libraries name is missing!"
  exit 1
fi

if [ -z ${PRIMER} ]
then
  echo "Error the primer name is missing!"
fi

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

FASTAFILE=$(find ${LIBDIR} -name "${LIBNAME}-${PRIMER}_grep_clean_g.fasta")
if [ ! -f "${FASTAFILE}" ]
then
  echo "Error, the fasta file ${LIBNAME}-${PRIMER}_grep_clean_g.fasta not exist!"
fi

DBREF=$(grep "${PRIMER}" ${DBINFO} | cut -d',' -f4)
if [ ! -e ${DBREF} ]
then
  echo 'Error: the reference database file '${DBREF}' not exist!'
  exit 1
fi

# process analysis
echo 'Process analysis'
echo `date`
echo "fasta file: "${FASTAFILE}
echo 'Primer: '$PRIMER

DBREF=$(grep "${PRIMER}" ${DBINFO} | cut -d',' -f4)
echo 'database reference: '${DBREF}

TAXODB=$(grep "${PRIMER}" ${DBINFO} | cut -d',' -f5)
echo 'taxo reference: '${TAXODB}

echo 'Fasta file used: '${FASTAFILE}

FILENAME=$(basename ${FASTAFILE%.*})

mkdir ${LIBDIR}/TMPDIR_ecotag_${LIBNAME}
cat ${FASTAFILE} | obiannotate -S hash:"hash(str(sequence))%100" | obisplit -t hash -p "${LIBDIR}/TMPDIR_ecotag_${LIBNAME}/${FILENAME}_split-"
for FMINI in $(find ${LIBDIR}/TMPDIR_ecotag_${LIBNAME} -name "${FILENAME}_split-*.fasta" | grep -v "ecotag.fasta")
do
	echo "ecotag -d ${TAXODB} -R ${DBREF} ${FMINI} > ${LIBDIR}/TMPDIR_ecotag_${LIBNAME}/$(basename ${FMINI%.*})_ecotag.fasta && rm ${FMINI}"
done | parallel -j 12

cat ${LIBDIR}/TMPDIR_ecotag_${LIBNAME}/*_ecotag.fasta > ${FASTAFILE%.*}_ecotag.fasta
obiannotate --delete-tag=father --delete-tag=fathers --delete-tag=clean ${FASTAFILE%.*}_ecotag.fasta | obiannotate -d ${TAXODB} -S path:'":".join([str(x[3])+"@"+taxonomy.getRank(x[0]) for x in taxonomy.parentalTreeIterator(taxid)][::-1])' > ${FASTAFILE%.*}_ecotag_annot.fasta
obitab -o ${FASTAFILE%.*}_ecotag_annot.fasta > ${FASTAFILE%.*}_ecotag_annot.tab

rm -r ${LIBDIR}/TMPDIR_ecotag_${LIBNAME}

echo `date`

exit $?
