!Hard limits on Multiscat parameters included everywhere at compile time

      parameter (NZFIXED_MAX = 550)        !max no of z points in fixed pot.
<<<<<<< Updated upstream
      parameter (NVFCFIXED_MAX = 32768)      !max no of fourier cmpts (fixed, from file)
=======
      parameter (NVFCFIXED_MAX = 100)      !max no of fourier cmpts (fixed, fro/Users/lazyb/Documents/Part-III-Multiscat-and-misc/multiscat.incm file)
>>>>>>> Stashed changes

      parameter (nmax = 1024)    !NOTE: This is set higher in sams repo, this vm eats too much ram!!      diffraction channels
      parameter (mmax = 550)    !         z grid points
      parameter (nfcx = 32768)    !         max number of potential fourier components per z slice
