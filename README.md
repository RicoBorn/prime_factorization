# prime_factorization
This repository contains logic for prime factorization using elliptic curves.

## Build of the project
For building this project you'll need to install the GMP-Library: https://gmplib.org/
  - Follow the instructions for installation (based on your system)
  - For macOS, install via homebrew (https://brew.sh/): ``brew install gmp``

After you have successfully installed the GMP-Library, this project can be built using either a direct ``g++`` command or CMake.
### Option 1: Direct Compilation with g++
For a quick setup, you can directly compile the project using the following command:
```bash
g++ -std=c++17 main.cpp test_numbers.cpp helper.cpp trial_division.cpp elliptic_curve.cpp -lgmp -lgmpxx
```
### Option 2: Build Using CMake
For a more robust and platform-independent build process, you can use CMake. This is the recommended method for larger projects or when working on different systems.
Create a build directory:
```bash
mkdir build
cd build
cmake ..
```
Build the project:
```bash
make
```