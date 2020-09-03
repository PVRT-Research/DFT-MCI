 ```
      [.....    [........[... [......   [..       [..    [..   [.. 
      [..   [.. [..           [..       [. [..   [... [..   [..[.. 
      [..    [..[..           [..       [.. [.. [ [..[..       [.. 
      [..    [..[......       [..       [..  [..  [..[..       [.. 
      [..    [..[..           [..    [..[..   [.  [..[..       [.. 
      [..   [.. [..           [..       [..       [.. [..   [..[.. 
      [.....    [..           [..    [..[..       [..   [....  [.. 

        +======================================================+   
        |  Density Functional Theory: Monte Carlo Integration  |   
        |  Version 0.2 (2019.alpha)                            |   
        |  RPV Research                                        |   
        |  Peverati Lab at Florida Tech                        |   
        +======================================================+   
 ```
# Overview
This is the evolution of our initial DFT:IoG (DFT:Integration on a Grid) program for integrating the xc functional in 3D.
It is now capable of performing the integration of the xc functional two ways:
- Monte Carlo Integration (with four different algorithms)
- Lebedev grid with Becke weights
It also has the following capabilities:
- Calculate the Coulomb energy via an inefficiently implemented pseudo-spectral method (it takes time!)
- Calculate -D3, gCP and BSE corrections (via a rearranged versio of the ppot program of Grimme)
- Calculate DFT-C corrections (J. Witte, M. Head-Gordon)

It reads several standard electronic density files from different sources (.chk, .fchk, .molden, .wfn, .wfx, etc.), thanks in part to the gMultiWFN engine https://github.com/stecue/gMultiwfn

# Dependencies:
1. math libraries (currently tested with lapack and lblas standard ones)
2. libxc (in principle this is not a stringent requirement, but for full xc functionality, it is. Available here: https://gitlab.com/libxc/libxc)
3. Cuba library for MC integration (provided, and also available here: http://www.feynarts.de/cuba/)

# Compilation
Hopefully it is just as straightforward as modify the Makefile and compile with make.
(Notice: this is still very untested on several systems, it might require some substantial tinkering around if compilation happens on computers that are different from those in our group!)

# Running a calculation
The simplest way of running a calculation is to provide a .fchk file as input (from Gaussian, q-chem, etc.):
```bash
./dftmci.exe h2.fchk
```
This will read the density from the h2.fchk file and integrate the xc functional using the default settings.
If more control over the setting is needed, a specific .iog input file can be provided, containing the relevant keywords. For an example of the allowed keywords and their explanation, see the in.iog sample file provided.

# Wishlist
Some straightforward (or maybe not so much) updates that might be considered:
- VV nonlocal dispersion (via noloco)
- XDM (via postg)