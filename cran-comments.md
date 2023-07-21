## Test environments

-   R-hub windows-x86_64-devel (r-devel)
-   R-hub ubuntu-gcc-release (r-release)
-   R-hub fedora-clang-devel (r-devel)

## R CMD check results

There were no ERRORs or WARNINGs.

There are 4 NOTEs.

#### 1. On windows-x86_64-devel (r-devel)

```         
checking for non-standard things in the check directory ... NOTE 
Found the following files/directories: 
  ''NULL''
```

As noted in [R-hub issue #560](https://github.com/r-hub/rhub/issues/560), this seems to be an Rhub issue and so can likely be ignored.

####2. On windows-x86_64-devel (r-devel)

```         
checking for detritus in the temp directory ... NOTE 
Found the following files/directories: 
  'lastMiKTeXException'
```

As noted in [R-hub issue #503](https://github.com/r-hub/rhub/issues/503), this could be due to a bug/crash in MiKTeX and can likely be ignored.

#### 3. On ubuntu-gcc-release (r-release), fedora-clang-devel (r-devel)

```         
checking examples ... NOTE 
Examples with CPU (user + system) or elapsed time \> 5s 
              user system elapsed 
SEMID-package 3.077 0.107 5.237
```

This is only slightly over 5s.

#### 4. On ubuntu-gcc-release (r-release), fedora-clang-devel (r-devel)

```         
checking HTML version of manual ... NOTE 
Skipping checking HTML validation: no command 'tidy' found
```

As noted in [R-hub issue #548](https://github.com/r-hub/rhub/issues/548), this seems to be an issue related to tidy and can likely be ignored.

## Downstream dependencies

There are currently no downstream dependencies for this package.
