## Test environments
Checks were performed on the following platforms using rhub::check_for_cran():
- R-hub windows-x86_64-devel (r-devel)
- R-hub ubuntu-gcc-release (r-release)
- R-hub fedora-clang-devel (r-devel)

## R CMD check results
* On windows-x86_64-devel (r-devel), ubuntu-gcc-release (r-release), fedora-clang-devel (r-devel)
  checking CRAN incoming feasibility ... NOTE
  Maintainer: 'Nils Sturma <nils.sturma@tum.de>'
  New maintainer:
    Nils Sturma <nils.sturma@tum.de>
  Old maintainer(s):
    Luca Weihs <lucaw@uw.edu>
  
  Found the following (possibly) invalid DOIs:
    DOI: 10.1214/10-AOS859
      From: DESCRIPTION
      Status: Internal Server Error
      Message: 500
    DOI: 10.1214/12-AOS1012
      Status: Internal Server Error
      From: DESCRIPTION
  
      Message: 500

0 errors √ | 0 warnings √ | 1 note x

- The maintainer is being changed on purpose.
- None of these DOIs are invalid. They are both correct.


## Downstream dependencies
There are currently no downstream dependencies for this package.