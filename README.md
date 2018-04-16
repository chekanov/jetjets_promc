This example shows how to create ROOT histograms with dijet cross sections 
using ProMC files from  the [HepSim Monte Carlo repository](http://atlaswww.hep.anl.gov/hepsim/). The example builds anti-KT jets, and fill ROOT histograms.


 1. Install ProMC (http://atlaswww.hep.anl.gov/asc/promc/), ROOT and FastJet 
 2. Check the installation. The variables: 

```
   echo $FASTJET 
   echo $PROMC
   echo $ROOTSYS
```
  they should return the installation paths. 

 3. Compile as "make"

 4. Run one of the scripts "A_RUN*" assuming $STORAGE variable is set. It points to the location of ProMC files. The script runs over directories with ProMC files and fills ROOT histograms.
