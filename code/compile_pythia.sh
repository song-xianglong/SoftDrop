source /cvmfs/sft.cern.ch/lcg/views/LCG_101_ATLAS_1/x86_64-centos7-gcc8-opt/setup.sh

export FJCONTRIB=/cvmfs/sft.cern.ch/lcg/releases/fjcontrib/1.046-8f6bb/x86_64-centos7-gcc8-opt/

./compile_pythia.sh lundplane_all

./lundplane_all.exe tune random_seeds numberEvents flavor(0-incl, 1-D, 2-B) radius
./lundplane_all.exe 5 0 10000 1 0.4 0
