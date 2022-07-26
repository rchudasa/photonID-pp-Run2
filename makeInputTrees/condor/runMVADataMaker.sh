#!/bin/bash
set -x
source /cvmfs/cms.cern.ch/cmsset_default.sh
export HOME=/home/rchudasa
export X509_USER_PROXY=/grid_mnt/t3home/rchudasa/grid_proxy/x509up_u56678

dir=/home/rchudasa/work/bsMM/CMSSW_10_6_29/src/photonID-pp-Run2/makeInputTrees
pwd 
uname -a

inputPath=""
outputPath=""
sampleName=""
configPath=""

suffix="_20July22_v1"

TMPDIR=`mktemp -d`
cd $TMPDIR

if [ $2 -eq 0 ]
then
  sampleName="flatPtPi0" # 10400 files
  configPath="$dir/configs/Pi0Sample.cfg" # 10400 files
  inputPath=`sed "${1}q;d" ${dir}/inputFiles/flatPi0.txt`
  outputPath="${dir}/workarea/mc_flat_pt_pi0/mvaTrees${suffix}"
elif [ $2 -eq 1 ]
then
  sampleName="BsMMG" # 10400 files
  configPath="$dir/configs/BsJpsiGamma.cfg" # 10400 files
  inputPath=`sed "${1}q;d" ${dir}/inputFiles/bsToMuMuGamma.txt`
  outputPath="${dir}/workarea/mc_bsMMG/mvaTrees${suffix}"
elif [ $2 -eq 2 ]
then
  sampleName="BsJpsiGamma" # 10400 files 
  configPath="$dir/configs/BsJpsiGamma.cfg" 
  inputPath=`sed "${1}q;d" ${dir}/inputFiles/bsToJpsiGamma.txt`
  outputPath="${dir}/workarea/mc_bsToJpsiGamma/mvaTrees${suffix}"
elif [ $2 -eq 3 ]
then
  sampleName="Data" # 10400 files
  configPath="$dir/configs/Pi0Sample.cfg" 
  inputPath=`sed "${1}q;d" ${dir}/inputFiles/data2018D.txt`
  outputPath="${dir}/workarea/data/mvaTrees${suffix}"

fi

cp $configPath .

mkdir -p $outputPath
output="mvaTrees_${1}.root"

cp $dir/mvaDataMaker.exe .
if [ -s ${output} ]
then
  echo "File already exists, skipping"
else
  echo "Running"
  echo $inputPath
  ./mvaDataMaker.exe $configPath $inputPath $output $3
  mv $output $outputPath/ 
fi
