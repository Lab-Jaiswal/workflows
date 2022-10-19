#!/bin/bash
#dx run cloud_workstation -imax_session_length=10h --ssh --brief -y --name "test"
cd ../..
dx download -r project-G5B07V8JPkg740v9GjfF9PzV:/References/
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

sudo apt-get update
sudo apt-get install -y parallel

echo "complete"

rm config.sh
cp cloud_config.sh config.sh
rm cloud_config.sh




