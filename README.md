# 2D-Ising-MonteCarlo

//************************Monte Carlo Assignment First Part*************************//

In order to use the codes, there is a Folder named Ising_Program, wich contains esentially:
	
	2DIsing.cpp : Program that runs the MC Ising Algorithm with Glauber Dynamics or Metropolis updating
	randomc.h   : Library that contains the mersenne random number generator
        mersenne.cpp: Random Number Generator
        userintf.cpp: User defined properties for the Random Number Generator
	mersenne.o  : Defined mersenne.cpp in order to avoid conflict between C/C++ Libraries (!use this)
	userintf.o  : Defined userintf.cpp in order to avoid conflict between C/C++ libraries (!use this)
	Jack.cpp    : Program that computes different Jacknife estimators
        BIN.cpp     : Program that computes the Binning analysis of the error in the measurements.
	
        make_Ising.sh : Shell script that contains the compiling instructions for 2DIsing.cpp

*************************************************************************************************************************************
In order to run 2D.Ising.cpp you only need to edit the file make_Ising.sh and change the values according to:

	-d (dimension 2) -L (Linear Size) -T (temperature) -nmcs (MC steps) -nmeas (measure every nmeas MC steps) -seed (seed of RNG)

Be aware that the program only works for d=2 up to now. The seed used for the assignment was: -seed 51712


This Will generate The temporary series for E and M in Energy.txt and Magnet.txt respectively


*************************************************************************************************************************************
In order to run the Jacknife estimator, just put the filename in the line 177, run g++ Jack.cpp -o Jack, then ./Jack and in the terminal you will see the options.

In order to run the Binning data, just put the filename in theline 46, run g++ Bin.cpp -o Bin, then ./Bin and it will show the result in the terminal.


Any Further Question:

Andres Henao Aristizabal
andreshenao31@gmail.com
ahenaoa@unal.edu.co
+34 698595499
