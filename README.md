
# How it works

## Description

MLXANES - machine learn x-ray absorption near edge structure: A fortran program to predict XANES spectra from an xyz file using machine learning. A multi-dimensional gradient descent is used to minimise the cost function (root-mean-squared distance between hypothesis and true spectra data), giving optimised coefficients to predict XANES from an xyz file.

## Important files

- compile\_mlxanes.sh , shell script to compile the program and all subroutines with gfortran
- mlxanes.f95 , the main code
- mlxanes , the main executable to run the program 
- x2\_gradient\_descent.f95 , the gradient descent algorithm
- functions.f95 , contains functions defining the cost function and other important functions

## Important directories

- xyz/ , directory containing all the xyz files

- train00/ , train01/ , ... , directories for training sets (spectra)

- predicted/ , the program writes all predicted spectra in this directory

- latex/ , directory containing latex and pdf files that describe the algorithm


## Algorithm description

1. Reads the setup.txt file and act according to the defined settings
2. If *learn* is 1 (on), read the training set defined by *itrain*, e.g. if *itrain* is set to 01 the files within the directory 'train01/' will be read. If *learn* is 0, the program skips to step 9. It reads file names 00001\_tddft.txt, 00002\_tddft.txt, ... (or 00001\_tddft\_conv.txt for convoluted spectra), up to a maximum of 10000\_tddft.txt (x = 10000 is hard coded currently).

- The expected format of the spectra files is:

```text
    Energy    <xanes>    
   7097.00 -3.0728018E-25
   7097.10  1.0060670E-23
   7097.20  1.7812559E-23
   7097.30 -1.9822585E-23
   7097.40  1.3765650E-23
   7097.50 -1.7296550E-23
   7097.60  1.6878165E-23
   7097.70 -2.3145454E-23
   7097.80  6.7490110E-24
```

3. Calculates the Coulomb matrix (CM) of size *m^2* for all corresponding xyz files that have spectra in *itrain*. Then either half the non-diagonals are chosen, or the reduced CM (RCM) of size _m_ if *rcm* is set to 1. Any xyz files with corresponding CM or RCM smaller than *m^2* or *m* respectively gets padded with zeros.

The full CM:

```text
| c11 c12 ... c1m |
| c21 c22 ... c2m |
| ...             |
| cm1 cm2 ... cmm |
```

The non-diagonals (rcm\_on = 0):

```text
|                     |
| c21                 |
| c31 c22             |
| ...                 |
| cm1 cm2 ... cm(m-1) |
```

The reduced CM (rcm\_on = 1):

```text
| c11 |
| c21 |
| ... |
| cm1 |
```

4. Either reads the convoluted or unconvoluted spectra in *itrain* (if *conv* is 1 or 0 respectively)
5. Reduces the energy grid to an index range *n1* to *n2* giving a total grid size of *n2*-*n1*+1
6. Checks if the file "coeffs\_train*itrain*.dat" exists, if it does it replaces it, if not it creates it
7. For each energy grid point a gradient descent is performed with the subroutine "x2\_gradient\_descent", with step-size *step*, maximum number of iterations *maxiter*, and convergence criteria is when the cost function (as defined in functions.f95) is less than *cutoff* between the *i*th and (*i*-1)th steps.

The cost function is,

J(θ) = Σⱼ ( hⱼ(θ) - yⱼ )²

where the sum is over the total *p* spectra in the *itrain* training set, yⱼ is the *j*th spectra, and hⱼ(θ) is *j*th hypothesis function, which is a multidimensional equation of a line,

h(θ) = θ₀ + θ₁x₁ + ... + θₘxₘ

there are *m*+1 coefficients (θ₀, θ₁, ..., θₘ) and "observations" (x₀, x₁, ..., xₘ) where formally we define x₀ := 1 (to make the vector θ and x have the same length), and θ₀ is the intercept. In this case the observations are the _m_ elements of the Coulomb matrix. If *percent* is enabled, the cost function is,

J(θ) = Σⱼ ( hⱼ(θ) - yⱼ )² / yⱼ

which does not work in regions where yⱼ is equal or close to 0. At each step the gradient is calculated for each θᵢ,

∂J(θ)/∂θᵢ = Σⱼ ( hⱼ(θ) - yⱼ ) xᵢ

or 

∂J(θ)/∂θᵢ = Σⱼ ( hⱼ(θ) - yⱼ ) xᵢ / yⱼ

if *percent* is enabled. Then each coefficient is iterated in the negative direction of the partial derivative,

θᵢ = θᵢ - α * ∂J(θ)/∂θᵢ

for step-size α. The starting value of α is defined by *step*. If *varystep* is enabled, every 5 steps α increases by a factor of 2, but if J(θ) increases compared to the previous step α decreases by a factor of 2. Every 5 steps and the factor 2 have been briefly tested and are fairly optimal. If *varystep* is disabled a constant step-size is used.

This ends one "step". The next step calculates the hypothesis and cost function and its partial derivatives with the new coefficients.

8. Once convergence is achieved or the maximum iterations are completed, the coefficients are written to the "coeffs\_train*itrain*.dat" file
9. Finally, the coefficients are read from "coeffs\_train*itrain*.dat" (if the file exists) and used to predict the spectra for all xyz files in the xyz directory. 

* The predicted spectrum at energy point *k* is the hypothesis function h(θₖ) where θₖ is the vector of *m*+1 optimised coefficients for the *k*th energy grid point.

The predicted spectra are also convoluted by the "lorenzian\_broadening" function with energy dependent widths calculated by the "fdmnes\_gamma" function (in functions.f95). The raw predictions and convoluted predicitons are written to the predicted/ directory.

## Setup file

Here is an example of the setup.txt file:
```
learn    1
cutoff   1.0e-12
jtarget  5.0e-06
maxiter  2000000
step     1.0e-08
mmax     1000
m        30
kmax     751
n1       1
n2       751
itrain   03
conv     0
jprint   1
rcm      1
varystep 1
percent  1
```
### Explanation

- _learn_: Can be 1 or 0, if 1 (on) the program will read files and perform the optimisation, otherwise it just predicts from the file tik.dat (which is generated coefficients from a previous optimisation). 

- _cutoff_: The cut-off value at which the gradient descent stops, specifically when the cost function differs by < _cutoff_ compared to the previous iteration it will stop (usually set to 1.0e-10). 

- _jtarget_: The absolute value of the cost function at which the gradient descent stops. 
 
Note: It stops when the _cutoff_ and _jtarget_ values have been achieved.

- _maxiter_: The maximum number of iterations until the gradient descent stops regardless of convergence (and it still prints the coefficient to file). 

- _step_: The step-size that the gradient descent takes at each iteration, larger values are faster but can miss the minimum and not converge, and too small values will take a long time or never reach the minimum so this value has to be tweaked occasionally (often is 0.001, or 1e-4, 1e-5 depending on the system). 

- _mmax_: The number of atoms to read from each xyz file (usually set this equal or greater to the largest total number of atoms). 

- _m_: The number of elements of the Coulomb matrix to use, its maximum possible value is *mmax*, but a lower value can be used to speed up the optimisation or help convergence (at the cost of ignoring hydrogens for example). 

- _kmax_: The total number of energy grid points to read from each spectrum file. 

- _n1_ and _n2_: Define the energy grid range the optimisation occurs in, in the above exmaple 1 and 751 define the entire range, as the maxium range is 1 to *kmax*. 

- _itrain_: Defines the training set to use in the optimisation, e.g. _itrain_ set to 03 means the program will read files from the train03/ directory. 

- _conv_: Read from convoluted or unconvoluted spectrum files; 1 (on) or 0 (off).

- _jprint_: Print convergence information at every iteration; 1 (on) or 0 (off).

- _rcm_: Use the reduced Coulomb matrix if on or full CM if off; 1 (on) or 0 (off).

- _varystep_: Use a variable step-size with starting step-size of _step_.

- _percent_: Use a percentage based cost function.

# Usage

### Compile
```
./compile_mlxanes.sh
```
to compile in serial, or edit the line in compile\_mlxanes\_openmp.sh with the number of threads,

```
export OMP_NUM_THREADS=64
```
and
```
./compile_mlxanes_openmp.sh
```
to compile with OpenMP parallisation.

The loop over the energy grid is parallelised. Therefore the maximum number of threads the program can utilise is *n2*-*n1*+1.

Note: It is not necessary to re-compile if setup.txt is changed.

### Run the program(s)

```
./mlxanes

```

# To install BLAS and LAPACK (may not be needed anymore)
```
sudo apt-get install libblas3gf libblas-doc libblas-dev liblapack3gf liblapack-doc liblapack-dev
```

