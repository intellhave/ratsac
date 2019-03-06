# RATSAC - Random Tree Sampling for Consensus Maximization

Description
---
This package contains the C++ implementation of the following paper

Huu Le, Tat-Jun Chin, and David Suter. *"RATSAC-Random Tree Sampling for Maximum Consensus Estimation."*, In Digital Image Computing: Techniques and Applications (DICTA), 2017


This implementation utilizes the framework of USAC:

Raguram, Rahul, Ondrej Chum, Marc Pollefeys, Jiri Matas, and Jan-Michael Frahm. *"USAC: a universal framework for random sample consensus."* IEEE transactions on pattern analysis and machine intelligence 35, no. 8 (2013): 2022-2038.

This code was tested on an Ubuntu Machine. The performances on Windows and MacOS have not been verified.

Dependencies
---
Since our code requires solving linear programs, the following libraries are required:
* Coin-Or Clp: https://projects.coin-or.org/Clp
* Some of the files also use the Armadillo package: http://arma.sourceforge.net/

To install the dependencies, please follow these steps:
* Install Coin-Or Clp: Follow the instructions in this link: https://projects.coin-or.org/Clp

* Install Armadillo package: 
 ` sudo apt-get install libarmadillo-dev`

Compiling
---
(The instruction is for Ubuntu; Compilation steps for Windows and MacOS should be changed accordingly).
Please install CMake then execute the following commands from the current directory:

~~~
mkdir build

cd build

cmake ..

make
~~~
Running 
---
After the code has been successfully compiled, cd to the `demo_script` folder and execute: `./runhomo.sh` to run a demo of homography estimation.