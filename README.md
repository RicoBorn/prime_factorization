# Prime Factorization using Elliptic Curves
This repository implements prime factorization algorithms utilizing elliptic curves, a modern and efficient approach to factorizing large numbers.

## Prerequisites
To build this project, you need the GMP (GNU Multiple Precision Arithmetic Library). Follow these steps to install GMP on your system:
- Visit the GMP official website (https://gmplib.org/) for installation instructions tailored to your platform.
- macOS users: You can install GMP easily via Homebrew (https://brew.sh/): ``brew install gmp``

Ensure the library is installed and accessible on your system before proceeding with the build.

## Building the project
It might be necessary to set the following environment variables before building (you may include them in your `.zshrc` (on macOS) or `.bashrc` (on Linux)):
```bash
# for C++ compiling and linking (gmp library)
export LIBRARY_PATH="<lib_path>:$LIBRARY_PATH"
export DYLD_LIBRARY_PATH="<lib_path>:$DYLD_LIBRARY_PATH"
export CPLUS_INCLUDE_PATH="<include>:$CPLUS_INCLUDE_PATH"
```
Replace `<lib_path>` with the path to the GMP library (e.g., `/opt/homebrew/opt/gmp/lib`) and `<include>` with the include path (e.g., `/opt/homebrew/include`).

You can build the project using one of the two following methods: direct compilation with `g++` or using CMake.
### Option 1: Direct Compilation with g++
For a quick setup, you can directly compile the project using the following command:
```bash
g++ -std=c++17 main.cpp test_numbers.cpp helper.cpp trial_division.cpp elliptic_curve.cpp -lgmp -lgmpxx -o prime_factorization
```
This will generate an executable named prime_factorization.
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

## Running the program
After building, you can run the program with:
```bash
./prime_factorization
```
The program expects input for numbers to factorize, and it will compute their prime factors using elliptic curve techniques.