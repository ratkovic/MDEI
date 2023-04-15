## R CMD check results

0 errors | 0 warnings | 1 note

* This is a new release.
* I have submitted several times.  As per our last correspondence, I have specified authors/maintainer/contributor via Authors@R
* I have add an example to the main function, MDEI
* I have removed references from the DESCRIPTION file and put them in the manual, using @references with roxygen2
* The random forest function ranger() was using multiple threads at its default that was throwing issues when checking.  The number of threads being used 
	by ranger is now an option in the function.
