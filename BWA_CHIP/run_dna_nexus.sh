#!/bin/bash
#dx run cloud_workstation -imax_session_length=10h --ssh --brief -y --name "test"
cd ../..
dx download -r project-G5B07V8JPkg740v9GjfF9PzV:/References/
cd References
dx download -r project-G5B07V8JPkg740v9GjfF9PzV:/References/cloud_references
cd ..

mkdir -p Inputs
cd Inputs
dx download project-G5B07V8JPkg740v9GjfF9PzV:/Bulk/Exome\ sequences/Exome\ OQFE\ CRAM\ files/10/1000184_23143_0_0.cram.crai

dx download project-G5B07V8JPkg740v9GjfF9PzV:/Bulk/Exome\ sequences/Exome\ OQFE\ CRAM\ files/10/1000184_23143_0_0.cram

dx download project-G5B07V8JPkg740v9GjfF9PzV:/Bulk/Exome\ sequences/Exome\ OQFE\ CRAM\ files/10/1000973_23143_0_0.cram.crai

dx download project-G5B07V8JPkg740v9GjfF9PzV:/Bulk/Exome\ sequences/Exome\ OQFE\ CRAM\ files/10/1000973_23143_0_0.cram

dx download project-G5B07V8JPkg740v9GjfF9PzV:/Bulk/Exome\ sequences/Exome\ OQFE\ CRAM\ files/10/1000819_23143_0_0.cram.crai

dx download project-G5B07V8JPkg740v9GjfF9PzV:/Bulk/Exome\ sequences/Exome\ OQFE\ CRAM\ files/10/1000819_23143_0_0.cram

cd ..

#git clone https://github.com/Lab-Jaiswal/workflows
cd References/cloud_references
mv GRCh* ../
mv 201* ../

cd ../..

sudo apt-get update
sudo apt-get install -y parallel

echo "complete"

rm ~/workflows/BWA_CHIP/config.sh
cp ~/workflows/BWA_CHIP/cloud_config.sh workflows/BWA_CHIP/config.sh
rm ~/workflows/BWA_CHIP/cloud_config.sh




