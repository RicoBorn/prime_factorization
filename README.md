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
You can customize the factorization process using several optional command-line arguments.

### Command-line arguments
| Option                     | Description                                                                                                          | Default Value                               |
|----------------------------|----------------------------------------------------------------------------------------------------------------------|---------------------------------------------|
| `--num_mode`, `-m`          | Specifies the type of number to be factored. Supported values: Fermat, Cunningham, Test, RSA. If not provided, the program defaults to `DirectNumber` mode, where the user inputs a custom number. | `DirectNumber`                              |
| `--stage1_bound`, `-b`      | Specifies the stage-1 bound (B) for elliptic curve factorization. This is the bound for prime bases.                | `DEFAULT_B` (per default a function of C)   |
| `--stage2_bound`, `-c`      | Specifies the stage-2 bound (C) for elliptic curve factorization. This is the bound for the smallest prime factor of N. | `DEFAULT_C` (per default a function of N)   |
| `--num_curves`, `-n`        | Specifies the number of elliptic curves used for factorization.                                                      | `DEFAULT_NUM_CURVES`                        |
| `--no_trial_division`, `-nt`| Disables trial division as a preliminary factorization step.                                                         | `false` (trial division enabled by default) |
| `--run_parallel`, `-p`      | Enables parallel execution of the factorization algorithm.                                                           | `false` (disabled by default)               |
| `--num_threads`, `-t`       | Specifies the number of threads to use for parallel factorization. Must be a positive integer and not exceed 100.    | `DEFAULT_NUM_THREADS` (10)                  |
