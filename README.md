# HMP

This is a C++ project for summited paper:

***Quantifying the Competitiveness of a Dataset in Relation to User Preferences*** 

We supply using cmake to compile our project in a linux machine 
(if you are using windows, you could download cygwin with gcc/g++ and cmake), 
and the default compile setting is using O3 optimization.
(If you are interesting in the compile details, please refer to CMakeLists.txt)

**Our project only uses standard C++ libraries, which means you don't need to download any other additional libs.**

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
cd hmp
```

List build generators, the system will tell you which generator you should use <br />
(e.g., a linux system would tell you "Unix Makefiles"):
```
cmake --help
``` 
set cmake configurations:
```
cmake -G "Unix Makefiles" . # sometimes not "Unix Makefiles", depend on your OS
```
build this project:
```
cmake --build . # it may should a lot of warnings because I set "-Wall" in CMakeLists.txt
```
After run command as listed above, you would generate a file called "heatmap" in linux machines or "heatmap.exe" in windows machines. 

Dataset
------
the full option dataset link: <br />
https://drive.google.com/drive/folders/1lOCzducFXr9herFnH761Re4Oa3EmvNRv?usp=sharing

Run our project:
-----
Before running our project, please download the dataset with the link mentioned above
and place the dataset as something like:<br />
./rtree/  <br />
./data/inde/ <br />
./data/cor/ <br />
./data/anti/ <br />
./data/real/  <br />
./log/ <br />
./doc/ <br />
./CMakeLists.txt <br />
balabalabala <br />
./main.cpp <br />
./heatmap <br />


A running example:
```
./heatmap -k 10 -d 4 -m MDA+ -h 5 -f ./data/inde/U500K4.dat
```
parameter explanations:
- k, topk
- d, dimension of the problem
- f, option data file, an option (a row of this file) can be represented as "id L1 L2 ... Ld U1 U2 ... Ud"
- m, method name, "UTK", "CSA", "CSA+", "MDA", "MDA+"
- h, the height of quad-tree, \lambda=2^h


Experiment reproduce:<br />
----
Figure 8: Effect of |P| 
```
./heatmap -k 10 -d 4 -m BL -h 5 -f ./data/inde/U100K4.dat
./heatmap -k 10 -d 4 -m BL -h 5 -f ./data/inde/U500K4.dat
./heatmap -k 10 -d 4 -m BL -h 5 -f ./data/inde/U1000K4.dat
./heatmap -k 10 -d 4 -m BL -h 5 -f ./data/inde/U5000K4.dat
./heatmap -k 10 -d 4 -m BL -h 5 -f ./data/inde/U10000K4.dat

./heatmap -k 10 -d 4 -m CSA -h 5 -f ./data/inde/U100K4.dat
./heatmap -k 10 -d 4 -m CSA -h 5 -f ./data/inde/U500K4.dat
./heatmap -k 10 -d 4 -m CSA -h 5 -f ./data/inde/U1000K4.dat
./heatmap -k 10 -d 4 -m CSA -h 5 -f ./data/inde/U5000K4.dat
./heatmap -k 10 -d 4 -m CSA -h 5 -f ./data/inde/U10000K4.dat

./heatmap -k 10 -d 4 -m CSA+ -h 5 -f ./data/inde/U100K4.dat
./heatmap -k 10 -d 4 -m CSA+ -h 5 -f ./data/inde/U500K4.dat
./heatmap -k 10 -d 4 -m CSA+ -h 5 -f ./data/inde/U1000K4.dat
./heatmap -k 10 -d 4 -m CSA+ -h 5 -f ./data/inde/U5000K4.dat
./heatmap -k 10 -d 4 -m CSA+ -h 5 -f ./data/inde/U10000K4.dat

./heatmap -k 10 -d 4 -m MDA -h 5 -f ./data/inde/U100K4.dat
./heatmap -k 10 -d 4 -m MDA -h 5 -f ./data/inde/U500K4.dat
./heatmap -k 10 -d 4 -m MDA -h 5 -f ./data/inde/U1000K4.dat
./heatmap -k 10 -d 4 -m MDA -h 5 -f ./data/inde/U5000K4.dat
./heatmap -k 10 -d 4 -m MDA -h 5 -f ./data/inde/U10000K4.dat

./heatmap -k 10 -d 4 -m MDA+ -h 5 -f ./data/inde/U100K4.dat
./heatmap -k 10 -d 4 -m MDA+ -h 5 -f ./data/inde/U500K4.dat
./heatmap -k 10 -d 4 -m MDA+ -h 5 -f ./data/inde/U1000K4.dat
./heatmap -k 10 -d 4 -m MDA+ -h 5 -f ./data/inde/U5000K4.dat
./heatmap -k 10 -d 4 -m MDA+ -h 5 -f ./data/inde/U10000K4.dat
```

Figure 9: Effect of d
```
./heatmap -k 10 -d 2 -m BL -h 5 -f ./data/inde/U500K2.dat
./heatmap -k 10 -d 3 -m BL -h 5 -f ./data/inde/U500K3.dat
./heatmap -k 10 -d 4 -m BL -h 5 -f ./data/inde/U500K4.dat
./heatmap -k 10 -d 5 -m BL -h 5 -f ./data/inde/U500K5.dat
./heatmap -k 10 -d 6 -m BL -h 5 -f ./data/inde/U500K6.dat
./heatmap -k 10 -d 7 -m BL -h 5 -f ./data/inde/U500K7.dat

./heatmap -k 10 -d 2 -m CSA -h 5 -f ./data/inde/U500K2.dat
./heatmap -k 10 -d 3 -m CSA -h 5 -f ./data/inde/U500K3.dat
./heatmap -k 10 -d 4 -m CSA -h 5 -f ./data/inde/U500K4.dat
./heatmap -k 10 -d 5 -m CSA -h 5 -f ./data/inde/U500K5.dat
./heatmap -k 10 -d 6 -m CSA -h 5 -f ./data/inde/U500K6.dat
./heatmap -k 10 -d 7 -m CSA -h 5 -f ./data/inde/U500K7.dat

./heatmap -k 10 -d 2 -m CSA+ -h 5 -f ./data/inde/U500K2.dat
./heatmap -k 10 -d 3 -m CSA+ -h 5 -f ./data/inde/U500K3.dat
./heatmap -k 10 -d 4 -m CSA+ -h 5 -f ./data/inde/U500K4.dat
./heatmap -k 10 -d 5 -m CSA+ -h 5 -f ./data/inde/U500K5.dat
./heatmap -k 10 -d 6 -m CSA+ -h 5 -f ./data/inde/U500K6.dat
./heatmap -k 10 -d 7 -m CSA+ -h 5 -f ./data/inde/U500K7.dat

./heatmap -k 10 -d 2 -m MDA -h 5 -f ./data/inde/U500K2.dat
./heatmap -k 10 -d 3 -m MDA -h 5 -f ./data/inde/U500K3.dat
./heatmap -k 10 -d 4 -m MDA -h 5 -f ./data/inde/U500K4.dat
./heatmap -k 10 -d 5 -m MDA -h 5 -f ./data/inde/U500K5.dat
./heatmap -k 10 -d 6 -m MDA -h 5 -f ./data/inde/U500K6.dat
./heatmap -k 10 -d 7 -m MDA -h 5 -f ./data/inde/U500K7.dat

./heatmap -k 10 -d 2 -m MDA+ -h 5 -f ./data/inde/U500K2.dat
./heatmap -k 10 -d 3 -m MDA+ -h 5 -f ./data/inde/U500K3.dat
./heatmap -k 10 -d 4 -m MDA+ -h 5 -f ./data/inde/U500K4.dat
./heatmap -k 10 -d 5 -m MDA+ -h 5 -f ./data/inde/U500K5.dat
./heatmap -k 10 -d 6 -m MDA+ -h 5 -f ./data/inde/U500K6.dat
./heatmap -k 10 -d 7 -m MAD+ -h 5 -f ./data/inde/U500K7.dat
```

Figure 10: Effect of k
```
./heatmap -k 1 -d 4 -m BL -h 5 -f ./data/inde/U500K4.dat
./heatmap -k 5 -d 4 -m BL -h 5 -f ./data/inde/U500K4.dat
./heatmap -k 10 -d 4 -m BL -h 5 -f ./data/inde/U500K4.dat
./heatmap -k 20 -d 4 -m BL -h 5 -f ./data/inde/U500K4.dat
./heatmap -k 40 -d 4 -m BL -h 5 -f ./data/inde/U500K4.dat

./heatmap -k 1 -d 4 -m CSA -h 5 -f ./data/inde/U500K4.dat
./heatmap -k 5 -d 4 -m CSA -h 5 -f ./data/inde/U500K4.dat
./heatmap -k 10 -d 4 -m CSA -h 5 -f ./data/inde/U500K4.dat
./heatmap -k 20 -d 4 -m CSA -h 5 -f ./data/inde/U500K4.dat
./heatmap -k 40 -d 4 -m CSA -h 5 -f ./data/inde/U500K4.dat

./heatmap -k 1 -d 4 -m CSA+ -h 5 -f ./data/inde/U500K4.dat
./heatmap -k 5 -d 4 -m CSA+ -h 5 -f ./data/inde/U500K4.dat
./heatmap -k 10 -d 4 -m CSA+ -h 5 -f ./data/inde/U500K4.dat
./heatmap -k 20 -d 4 -m CSA+ -h 5 -f ./data/inde/U500K4.dat
./heatmap -k 40 -d 4 -m CSA+ -h 5 -f ./data/inde/U500K4.dat

./heatmap -k 1 -d 4 -m MDA -h 5 -f ./data/inde/U500K4.dat
./heatmap -k 5 -d 4 -m MDA -h 5 -f ./data/inde/U500K4.dat
./heatmap -k 10 -d 4 -m MDA -h 5 -f ./data/inde/U500K4.dat
./heatmap -k 20 -d 4 -m MDA -h 5 -f ./data/inde/U500K4.dat
./heatmap -k 40 -d 4 -m MDA -h 5 -f ./data/inde/U500K4.dat

./heatmap -k 1 -d 4 -m MDA+ -h 5 -f ./data/inde/U500K4.dat
./heatmap -k 5 -d 4 -m MDA+ -h 5 -f ./data/inde/U500K4.dat
./heatmap -k 10 -d 4 -m MDA+ -h 5 -f ./data/inde/U500K4.dat
./heatmap -k 20 -d 4 -m MDA+ -h 5 -f ./data/inde/U500K4.dat
./heatmap -k 40 -d 4 -m MDA+ -h 5 -f ./data/inde/U500K4.dat
```

Figure 11: Effect of \lambda
```
./heatmap -k 10 -d 4 -m BL -h 3 -f ./data/inde/U500K4.dat
./heatmap -k 10 -d 4 -m BL -h 4 -f ./data/inde/U500K4.dat
./heatmap -k 10 -d 4 -m BL -h 5 -f ./data/inde/U500K4.dat
./heatmap -k 10 -d 4 -m BL -h 6 -f ./data/inde/U500K4.dat
./heatmap -k 10 -d 4 -m BL -h 7 -f ./data/inde/U500K4.dat

./heatmap -k 10 -d 4 -m CSA -h 3 -f ./data/inde/U500K4.dat
./heatmap -k 10 -d 4 -m CSA -h 4 -f ./data/inde/U500K4.dat
./heatmap -k 10 -d 4 -m CSA -h 5 -f ./data/inde/U500K4.dat
./heatmap -k 10 -d 4 -m CSA -h 6 -f ./data/inde/U500K4.dat
./heatmap -k 10 -d 4 -m CSA -h 7 -f ./data/inde/U500K4.dat

./heatmap -k 10 -d 4 -m CSA+ -h 3 -f ./data/inde/U500K4.dat
./heatmap -k 10 -d 4 -m CSA+ -h 4 -f ./data/inde/U500K4.dat
./heatmap -k 10 -d 4 -m CSA+ -h 5 -f ./data/inde/U500K4.dat
./heatmap -k 10 -d 4 -m CSA+ -h 6 -f ./data/inde/U500K4.dat
./heatmap -k 10 -d 4 -m CSA+ -h 7 -f ./data/inde/U500K4.dat

./heatmap -k 10 -d 4 -m MDA -h 3 -f ./data/inde/U500K4.dat
./heatmap -k 10 -d 4 -m MDA -h 4 -f ./data/inde/U500K4.dat
./heatmap -k 10 -d 4 -m MDA -h 5 -f ./data/inde/U500K4.dat
./heatmap -k 10 -d 4 -m MDA -h 6 -f ./data/inde/U500K4.dat
./heatmap -k 10 -d 4 -m MDA -h 7 -f ./data/inde/U500K4.dat

./heatmap -k 10 -d 4 -m MDA+ -h 3 -f ./data/inde/U500K4.dat
./heatmap -k 10 -d 4 -m MDA+ -h 4 -f ./data/inde/U500K4.dat
./heatmap -k 10 -d 4 -m MDA+ -h 5 -f ./data/inde/U500K4.dat
./heatmap -k 10 -d 4 -m MDA+ -h 6 -f ./data/inde/U500K4.dat
./heatmap -k 10 -d 4 -m MDA+ -h 7 -f ./data/inde/U500K4.dat
```

Figure 12: Synthetic distributions and real product sets
```
./heatmap -k 10 -d 4 -m BL -h 5 -f ./data/cor/COR500K4.dat
./heatmap -k 10 -d 4 -m BL -h 5 -f ./data/inde/U500K4.dat
./heatmap -k 10 -d 4 -m BL -h 5 -f ./data/anti/ANTI500K4.dat
./heatmap -k 10 -d 4 -m BL -h 5 -f ./data/real/HOTEL4D.dat
./heatmap -k 10 -d 6 -m BL -h 5 -f ./data/real/HOUSE6D.dat
./heatmap -k 10 -d 8 -m BL -h 5 -f ./data/real/NBA8D.dat

./heatmap -k 10 -d 4 -m CSA -h 5 -f ./data/cor/COR500K4.dat
./heatmap -k 10 -d 4 -m CSA -h 5 -f ./data/inde/U500K4.dat
./heatmap -k 10 -d 4 -m CSA -h 5 -f ./data/anti/ANTI500K4.dat
./heatmap -k 10 -d 4 -m CSA -h 5 -f ./data/real/HOTEL4D.dat
./heatmap -k 10 -d 6 -m CSA -h 5 -f ./data/real/HOUSE6D.dat
./heatmap -k 10 -d 8 -m CSA -h 5 -f ./data/real/NBA8D.dat

./heatmap -k 10 -d 4 -m CSA+ -h 5 -f ./data/cor/COR500K4.dat
./heatmap -k 10 -d 4 -m CSA+ -h 5 -f ./data/inde/U500K4.dat
./heatmap -k 10 -d 4 -m CSA+ -h 5 -f ./data/anti/ANTI500K4.dat
./heatmap -k 10 -d 4 -m CAS+ -h 5 -f ./data/real/HOTEL4D.dat
./heatmap -k 10 -d 6 -m CSA+ -h 5 -f ./data/real/HOUSE6D.dat
./heatmap -k 10 -d 8 -m CSA+ -h 5 -f ./data/real/NBA8D.dat

./heatmap -k 10 -d 4 -m MDA -h 5 -f ./data/cor/COR500K4.dat
./heatmap -k 10 -d 4 -m MDA -h 5 -f ./data/inde/U500K4.dat
./heatmap -k 10 -d 4 -m MDA -h 5 -f ./data/anti/ANTI500K4.dat
./heatmap -k 10 -d 4 -m MDA -h 5 -f ./data/real/HOTEL4D.dat
./heatmap -k 10 -d 6 -m MDA -h 5 -f ./data/real/HOUSE6D.dat
./heatmap -k 10 -d 8 -m MDA -h 5 -f ./data/real/NBA8D.dat

./heatmap -k 10 -d 4 -m MDA+ -h 5 -f ./data/cor/COR500K4.dat
./heatmap -k 10 -d 4 -m MDA+ -h 5 -f ./data/inde/U500K4.dat
./heatmap -k 10 -d 4 -m MDA+ -h 5 -f ./data/anti/ANTI500K4.dat
./heatmap -k 10 -d 4 -m MDA+ -h 5 -f ./data/real/HOTEL4D.dat
./heatmap -k 10 -d 6 -m MDA+ -h 5 -f ./data/real/HOUSE6D.dat
./heatmap -k 10 -d 8 -m MDA+ -h 5 -f ./data/real/NBA8D.dat

