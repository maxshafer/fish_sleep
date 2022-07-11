## Installing treePL on cluster, ugh so much hard

# First login and make a new directory to install stuff into, because I don't have permission in /usr

mkdir -p ~/.local/{bin,lib,src,include}

# cd and clone the github (change to https)

cd ~/.local/src

git clone https://github.com/blackrim/treePL.git

# Install dependencies first

cd treePL/dep/

tar -xvzf nlopt-2.4.2.tar.gz

cd nlopt-2.4.2

./configure --prefix ~/.local/
make
make install


## Install ADOL-C

cd ~/.local/src/ADOL-C

./configure --prefix ~/.local/ --with-openmp-flag=-fopenmp

*****************************************************************************

  To successfully compile and run programs using the ADOL-C shared library do
  the following things:
     compiling:
        * add "-I/scicore/home/schiera/gizevo30/.local/include" to your compiler call
     linking:
        * add "-L${exec_prefix}/lib64 -ladolc" to your linker call
        * extend your linker call by "-Wl,--rpath -Wl,${exec_prefix}/lib64"
          (if you wish to skip the point "executing")
     executing (do one of the following things):
        * add ${exec_prefix}/lib64 to your LD_LIBRARY_PATH variable
        * ask your system administrator for adding ${exec_prefix}/lib64 to
          the global file containing library search paths (/etc/ld.so.conf)

     (or use the static library by replacing
        "-L${exec_prefix}/lib64 -ladolc" with
        "${exec_prefix}/lib64/libadolc.a" when linking)

  See README for instructions on how to change to a nonlocal installation!

*****************************************************************************

make
make install



# Once the software is properly installed, make sure to include the .local directories in the corresponding paths:

export PATH=~/.local/bin:~/.local:$PATH
export LD_LIBRARY_PATH=~/.local/lib:~/.local/lib64:${exec_prefix}/lib64:$LD_LIBRARY_PATH
export LIBRARY_PATH=~/.local/lib:~/.local/lib64:$LIBRARY_PATH
export CPATH=~/.local/include:~/.local/include/adolc:$CPATH


# The above seems to work, but when installing treePL, it can't find it, BUT
# I can load the correct module, and treePL installation maybe works
# Will likely need to load this each time

# module load NLopt/2.4.2-GCCcore-7.3.0

cd ~/.local/src

git clone https://github.com/blackrim/treePL.git

cd ~/.local/src/treePL/src/

module load foss

# This works, if I preload NLopt

./configure --prefix ~/.local/

# I then have to modify the Makefile, to change the paths (any that look weird)
# This seems to work eventually, not sure of the order
# Then just make sure all the proper paths are exported in the SLURM script

# But then this doesn't work, because it doesn't know to look for adolc in the right place
make














