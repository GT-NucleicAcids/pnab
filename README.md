# pNAB

## Install instructions

### Dependencies
[Openbabel](https://github.com/openbabel/openbabel), GCC 4+, CMAKE 2.4+\n
(Optional) Doxygen

For an Ubuntu system, run the following
```
sudo apt-get install libxml2-dev zlib1g-dev libeigen2-dev libeigen3-dev libcairo2-dev doxygen
git clone https://github.com/openbabel/openbabel.git
mkdir build
cd build
cmake -DCMAKE_BUILD_TYPE=Release ../openbabel/
sudo make install -j
```
For users who have already installed Openbabel, all dependencies should already be satisfied. Simply run the following:
```
git clone https://github.gatech.edu/jbarnett8/pNAB.git
mkdir build
cd build
cmake -DCMAKE_BUILD_TYPE=Release ../pNAB/
make all -j
```
This should build the executable as well as associated documentation as a webpage
