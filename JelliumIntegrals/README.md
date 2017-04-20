Contains a program which will (hopefully) compute all 1- and 2-electron integrals + the self energy and store them in files.

2-electron repulsion integrals ~ half of all integrals are computed -> ERI.dat
1-electron nuclear attraction integrals - symmetry fully exploited  -> NucAttraction.dat
1-electron kinetic energy integrals     - symmetry fully exploited  -> Kinetic.dat
Nuclear Self-Energy - treat same as nuclear repulsion in H2O case   -> SelfEnergy.dat

Tips for compling and running 
(1) Attempt to use compiler optimizers, so compile using
g++ -O2 -o JelliumIntegrals.exe JelliumIntegrals.c

(2) Give the number of grid points and the start and end of the domain as arguments when you execute the program and run in the background
./JelliumIntegrals.exe 50 0 1 >& OUT.txt &

Go enjoy a cup of coffee, its going to take > 20 minutes!
