# Inference-using-FHE
To run the code do the following steps --
1. Install OpenFHE library: https://github.com/openfheorg/openfhe-development following the steps given here: https://github.com/openfheorg/openfhe-development#installation
2. Check that you do not get error messages at the time of running "make install" for OpenFHE. You may need the admin rights to install.
3. Go to a new directory where you will keep “NND.cpp” and "N_N.cpp" file attached here.
4. Copy openfhe-development/CMakeLists.User.txt to the new directory and rename it to CMakeLists.txt.
5. Open CMakeLists.txt for editing and add a line to its end as suggested in the comments in CMakeLists.txt. Something like this:
```
    add_executable( N_N N_N.cpp )
    add_executable( NND NND.cpp )
```
6. ... and after that, execute commands that are very similar to the commands to build and run examples in OpenFHE:
```
    mkdir build
    cd build
    cmake ..
    make
    and, finally, run ./N_N or ./NND
```

