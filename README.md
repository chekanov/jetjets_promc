This example shows how to create ROOT histograms with dijet cross sections 
using ProMC files from  the [HepSim Monte Carlo repository](http://atlaswww.hep.anl.gov/hepsim/). The example calculates invariant masses of two jets in events
associated with a lepton. 


 1. Install ProMC (http://atlaswww.hep.anl.gov/asc/promc/), ROOT and FastJet 
 2. Check the installation. The variables: 

```
   echo $FASTJET 
   echo $PROMC
   echo $ROOTSYS
```
  should return the installation paths. 

 3. Compile as "make". You may need to rebuild the library inside "lib", i.e. "cd lib; make; make lib". 

 4. Before running "ana" executable, make sure you use the setup script:

```
export LD_LIBRARY_PATH=$FASTJET/lib:$LD_LIBRARY_PATH
export LD_LIBRARY_PATH=./lib/src/:$LD_LIBRARY_PATH
echo "STORAGE=$STORAGE"
```

where "STORAGE" points to a directory with ProMC files. 
The executable uses "data.in" file that points to the absolute path of ProMC files
(one file per line). To help creating this file, Make_input file is used.

As an example, you can  run scripts that start  with "A_RUN", 
assuming the STORAGE variable is set. 
Each script runs over directories with ProMC files created with different masses of
particles that decay to dijets, and fills ROOT histograms.

S.Chekanov (ANL) 
