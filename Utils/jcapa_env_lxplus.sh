#load bashrc
source ~/.bashrc

#load config files
source /afs/cern.ch/sw/lcg/contrib/gcc/4.9.3/x86_64-slc6/setup.sh
source /afs/cern.ch/sw/lcg/app/releases/ROOT/5.34.36/x86_64-slc6-gcc49-opt/root/bin/thisroot.sh
alias cmake='/afs/cern.ch/sw/lcg/external/CMake/3.4.3/Linux-x86_64/bin/cmake'

#change prompt
PS1="(jcapa) [\u@\h \W]\$ "
