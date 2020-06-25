A C++ toolbox for computing Discrete and Fast Fourier Transforms (DFT,FFT), Power Spectral Density (PSD) estimates, and the sound pressure level (SPL) in (dB).

Example for usage:
----------------------

If you have an input file with column formats:\
       time   signal_1     signal_2     signal_3 
      
line1  0.0_____ 0.0          1.0         -0.3     

line2  0.1     2.5          2.0          0.7   

....   ...     ...          ...          ...  


Linux
Inside the cloned directory $ mkdir obj , $ mkdir bin, before compile
1-compile the code by running make inside the code directory after cloning. 
2-run the code as follows:
$ ./bin/fftpsd -i inputfile.dat -n 2048 
for fft only.
or
$ ./bin/fftpsd -i inputfile.dat -n -spl 
for computing the sound pressure level
The above commands by default uses the first two columns in the input file only. 
If you want to compute fft/spl for other signals use the following:
$ ./bin/fftpsd -i inputfile.dat -n 2048 -spl -c 3 -r 10 -q 100
This means use signal_2 and start from row 10 to row 100.

**Windows:
Inside the cloned directory make the two obj and bin directories.
1- open a cmd window and navigate to the cloned directory.
2- compile by typing make in the cmd
3- run the code as the case for Linux above.

** An update will be pushed to make the obj and bin directories automatically. 

For more options please use the following table:

help with implemented options
-------------------------------
       minimum arguments for fft computations
       --------------------------------------
           [-i datafile] 
           [-n N] 

       minimum arguments for psd, power & spl computations
       -----------------------------------------------------
           [-i datafile]
           [-n N] 
           [-psd] or [-pow] or [-spl]

       If DFT is preferred just add to your arguments the flag [-dft]

Detailed option list and discription
---------------------------------------
    [-i datafile]  name of the data file with its address
    [-o directory] output directory, default is the same as input
    [-spec] for output files with headers and spec suffix
    [-r first_row]  of the signal data in the data file, default 1
    [-q last_row ]  of the signal data in the data file, default is the end of the file
    [-c column] of the signal data in the data file, default 2, assuming time is at the first column

    [-n N] number of data in one window subset, or in the whole sample if no shifting/averaging, must be 2^{k}, k is an integer
    [-l Lt_sub] window length in sec or the whole sample length if no shifting/averaging
    [-s shift] a value in [0.,1.] to indicate the ratio of data to be shifted
             shift defaults: 0.5 for psd/pow/spl, 0. for fft
    [-dt time_step] of the window subset or the whole sample in sec
    [-w window] name of the window functions, for psd default is hann, for fft default is rectangular
              windows: hann,hamm,bartlett,welch,blackman,triangular,rectangular
    [-m number] either 1 or 0 for mean substract, default is 1 to substract the mean

    [-psd filename]      compute the power spectral density (psd), filename is optional
    [-pow filename]      compute the power spectrum, filename is optional
    [-spl filename]      compute the sound pressure level (SPL), filename is optional
    [-peak]     peak preserving averaging fft mode, useful for ensemble averaging with shift=1.
    [-variance] variance preserving averaging fft mode for conserving the signal energy, by default it is forced for psd, power and spl computations
    [-pref pressure] reference pressure for SPL computation, default is 2e-5 Pa
    [-bin] bin-averaging for psd and spl
    [-oct number] the number nth octave averaging that is required, default is 1
