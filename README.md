# Spaital Heatmap

This is a C++ project for summited paper:

***Quantifying the Competitiveness of a Market in Relation to User Preferences*** 

We supply using cmake to compile our project in a linux machine 
(if you are using windows, you could download cygwin with gcc/g++ and cmake), 
and the default compile setting is using O3 optimization.
(If you are interesting in the compile details, please refer to CMakeLists.txt)

Our project only use standard C++ libraries, which means you don't need to download any additional libs.

Download make and cmake
-----
(usually a linux machine should contain make, you could try with command "make -version")
if your machine doesn't contain make, you can follow the steps shown as below:

```
sudo apt-get install build-essential
```

try "cmake -version" to see if your machine contains cmake, if doesn't, install cmake as: 
``` 
sudo snap install cmake --classic 
```

To compile our project 
-------

Go to where download this project, you could also go to any where else you want.
``` 
cd ~ 
```
download this project:
```
git clone https://github.com/efbeb99/hmp.git 
```
go to the directory of this project:
```
cd klevel
```
make a build directory:
```
mkdir build
``` 
go to build directory:
```
cd build
``` 
List build generators, the system will tell you which generator you should use <br />
(e.g., a linux system would tell you "Unix Makefiles"):
```
cmake --help
``` 
set cmake configurations:
```
cmake -G "Unix Makefiles" .. # sometimes not "Unix Makefiles", depend on your OS
```
build this project:
```
cmake --build . # it may should a lot of warnings because I set "-Wall" in CMakeLists.txt
```
After run command as listed above, you would generate a file called "heatmap" in linux machines or "heatmap.exe" in windows machines. 

Run our project:
-----
a running example:
```
./heatmap -k 10 -d 4 -m MDA+ -h 5 -f U500K4.dat
```
parameter explanations:
- k, topk
- d, dimension of the problem
- f, option data file, an option (a row of this file) can be represented as "id L1 L2 ... Ld U1 U2 ... Ud"
- m, method name, "UTK", "CSA", "CSA+", "MDA", "MDA+"
- h, the height of quad-tree, \lambda=2^h

the full option dataset link: <br />
**TODO** 