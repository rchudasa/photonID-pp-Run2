#!/bin/bash
set -x
source /cvmfs/cms.cern.ch/cmsset_default.sh
export HOME=/home/rchudasa
export X509_USER_PROXY=/grid_mnt/t3home/rchudasa/grid_proxy/x509up_u56678

dir=/home/rchudasa/work/bsMM/photonID-pp-Run2/makeInputTrees
pwd 
uname -a

#root setup
#source /cvmfs/sft.cern.ch/lcg/contrib/gcc/4.8/x86_64-centos7-gcc48-opt/setup.sh
#source /cvmfs/sft.cern.ch/lcg/app/releases/ROOT/6.24.06/x86_64-centos7-gcc48-opt/bin/thisroot.sh

inputPath=""
outputPath=""
sampleName=""
configPath=""

suffix="_20July22_v7"

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
  sampleName="flatPtPhoton" # 10400 files
  configPath="$dir/configs/SinglePhoton.cfg" # 10400 files
  inputPath="${basePath}/mc_flat_pt_photon/flatPtPhoton_5p02TeV_PbPb/flatPtPhoton_HiForest_v2/220411_102132/0000/flatPtMC_HiForestAOD_${1}.root"
  outputPath="${basePath}/mc_flat_pt_photon/mvaTrees${suffix}"
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
