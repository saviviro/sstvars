# sstvars 1.0.0

* Initial CRAN submission.

# sstvars 1.0.1

* Updated configure script to fix an issue with the installation on Mac OS X.

# sstvars 1.0.2

* Updated readme.
* Updated documentation.

# sstvars 1.1.0

* Implemented independent skewed t distribution as a new conditional distribution.
* Implemented a multiple phase estimation for TVAR models to enhance computational efficiency.
* Fixed a bug in the simulation algorithm for models incorporating independent Student's t conditional distributions
  (the variance of each structural shock was not scaled to one). 
* The argument standard_error_print can now be used directly in the summary function to obtain printout of standard errors. 
* Updated the documentation. 
