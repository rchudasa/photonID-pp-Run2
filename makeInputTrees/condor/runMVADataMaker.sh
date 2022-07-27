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

suffix="_barrel"

TMPDIR=`mktemp -d`
cd $TMPDIR

if [ $2 -eq 0 ]
then
  sampleName="flatPtPi0" # 4 files
  configPath="$dir/configs/Pi0Sample.cfg" 
  inputPath=`sed "${1}q;d" ${dir}/inputFiles/flatPi0.txt`
  outputPath="${dir}/workarea/mc_flat_pt_pi0/mvaTrees${suffix}"
elif [ $2 -eq 1 ]
then
  sampleName="QCD-15To7000" #60 files 
  configPath="$dir/configs/Pi0Sample.cfg" 
  inputPath=`sed "${1}q;d" ${dir}/inputFiles/qcd15To7000.txt`
  outputPath="${dir}/workarea/mc_qcd15To7000/mvaTrees${suffix}"
elif [ $2 -eq 2 ]
then
  sampleName="QCD-30To50" 
  configPath="$dir/configs/Pi0Sample.cfg"
  inputPath=`sed "${1}q;d" ${dir}/inputFiles/qcd30To50.txt`
  outputPath="${dir}/workarea/qcd30To50/mvaTrees${suffix}"
elif [ $2 -eq 3 ]
then
  sampleName="QCD-50To80"
  configPath="$dir/configs/Pi0Sample.cfg"
  inputPath=`sed "${1}q;d" ${dir}/inputFiles/qcd50To80.txt`
  outputPath="${dir}/workarea/qcd50To80/mvaTrees${suffix}"
elif [ $2 -eq 4 ]
then
  sampleName="QCD-20To30EMEnriched" 
  configPath="$dir/configs/Pi0Sample.cfg" 
  inputPath=`sed "${1}q;d" ${dir}/inputFiles/qcd20To30EmEnriched.txt`
  outputPath="${dir}/workarea/qcd20To30EmEnriched/mvaTrees${suffix}"
elif [ $2 -eq 5 ]
then
  sampleName="QCD-30To50EMEnriched" 
  configPath="$dir/configs/Pi0Sample.cfg" 
  inputPath=`sed "${1}q;d" ${dir}/inputFiles/qcd30To50EmEnriched.txt`
  outputPath="${dir}/workarea/qcd30To50EmEnriched/mvaTrees${suffix}"
elif [ $2 -eq 6 ]
then
  sampleName="BsMMG" 
  configPath="$dir/configs/BsJpsiGamma.cfg" 
  inputPath=`sed "${1}q;d" ${dir}/inputFiles/bsToMuMuGamma.txt`
  outputPath="${dir}/workarea/mc_bsMMG/mvaTrees${suffix}"
elif [ $2 -eq 7 ]
then
  sampleName="BsJpsiGamma" 
  configPath="$dir/configs/BsJpsiGamma.cfg" 
  inputPath=`sed "${1}q;d" ${dir}/inputFiles/bsToJpsiGamma.txt`
  outputPath="${dir}/workarea/mc_bsToJpsiGamma/mvaTrees${suffix}"
elif [ $2 -eq 8 ]
then
  sampleName="Data" 
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
