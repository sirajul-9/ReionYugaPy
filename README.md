# ReionYugaPy

This repository contains a Python wrapper for the C code, ReionYuga. The original code is located at https://github.com/rajeshmondal18/ReionYuga. 

For running this code, the outputs of N-body code (https://github.com/rajeshmondal18/N-body) and Friends-of-Friend Halo Finder(https://github.com/rajeshmondal18/FoF-Halo-finder) are needed.



# Required Libraries

CFFI: can be installed with <pre>  pip install cffi </pre>
for local installation (without root access) you may use <pre> pip install --user cffi </pre> or <pre> pip3 install --user cffi </pre>
FFTW-3 library needs to be installed in the system. It can be downloaded with 
<pre> sudo apt-get install libfftw3-dev </pre>

If that does not work download FFTW-3 from http://www.fftw.org/fftw3_doc/Installation-on-Unix.html#Installation-on-Unix Then install it with the flags --enable-float, --enable-shared and --enable-openmp as written below

After extracting the downloaded file go to the directory and do the following.
<pre> 
./configure --enable-float --enable-shared --enable-openmp
 make
 sudo make install
</pre>

For installing locally without root access:
<pre> 
./configure --enable-float --enable-shared --enable-openmp --prefix=/path/to/your/local/install/directory
 make
 make install
</pre>
 

# Instructions for Running the Code

Download all the files in this repository as ZIP and extract. Or use the following command:
<pre>
 git clone https://github.com/sirajul-9/ReionYugaPy
</pre>

If you have installed FFTW locally, edit the "build_funcs.py" file. Add the following arguments to the "ffi.set_source()" function:
<pre>
include_dirs=['/path/to/local/include/directory'],

library_dirs=['/path/to/local/lib/directory'],
</pre>

For example
<pre>
 ffi.set_source(
    "funcs_cffi",  
    """
    #include "source.h"
    """,
    sources=["source.c"],
    libraries=["m", "fftw3f", "fftw3f_omp"],
    extra_compile_args=['-w','-fopenmp', '-std=c99'],
    library_dirs=["/home/sirajul/local_install/lib"],
    include_dirs=["/home/sirajul/local_install/include/"]
)
</pre>

Write the values of all the parameters and N-body and halo_catalog outputs' path in the file input_parameters.py.

Then run the file build_funcs.py with 
<pre>
python3 build_funcs.py
</pre>
If you are running with locally installed libraries run the following command (or write it in your job submission script):
<pre> 
LD_LIBRARY_PATH=/path/to/local/lib/directory
export LD_LIBRARY_PATH
</pre>
Then the main code can be run using:
<pre>
python3 ionz_main.py
</pre>
