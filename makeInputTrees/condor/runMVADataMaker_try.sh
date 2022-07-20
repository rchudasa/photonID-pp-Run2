#!/bin/bash
set -x
source /cvmfs/cms.cern.ch/cmsset_default.sh
export HOME=/home/rchudasa
export X509_USER_PROXY=/grid_mnt/t3home/rchudasa/grid_proxy/x509up_u56617

dir=/home/rchudasa/work/bsMM/photonID-pp-Run2/makeInputTrees
pwd 
uname -a

suffix="_18April22"

TMPDIR=`mktemp -d`
cd $TMPDIR

sampleName="flatPtPi0" # 10400 files
configPath="$dir/configs/Pi0Sample.cfg" # 10400 files
inputPath=`sed "2q;d" ${dir}/inputFiles/flatPi0.txt`
outputPath="${dir}/workarea/mc_flat_pt_pi0/mvaTrees${suffix}"

cp $configPath .
cp $inputPath .

mkdir -p $outputPath
output="${outputPath}/mvaTrees_${1}.root"

cp $dir/mvaDataMaker.exe .
#cp $dir/configs/Pi0Sample.cfg . 
if [ -s ${output} ]
then
  echo "File already exists, skipping"
else
  echo "Running"
  echo $configPath
  #./mvaDataMaker.exe $configPath $inputPath $output $3
  ./mvaDataMaker.exe $configPath $inputPath $output 1
fi 
