UKF TRACTOGRAPHY
================


About   
-----

We present a framework which uses an unscented Kalman filter for performing
tractography. At each point on the fiber the most consistent direction is found
as a mixture of previous estimates and of the local model.

It is very easy to expand the framework and to implement new fiber representations 
for it. Currently it is possible to tract fibers using two different 1-, 2-, or 3-tensor 
methods. Both methods use a mixture of Gaussian tensors. One limits the diffusion 
ellipsoids to a cylindrical shape (the second and third eigenvalue are assumed to be 
identical) and the other one uses a full tensor representation.

__Authors__:
Yogesh Rathi (yogesh@bwh.harvard.edu), Stefan Lienhard, Yinpeng Li, Martin
Styner, Ipek Oguz, Yundi Shi, Christian Baumgartner (c.f.baumgartner@gmail.com)
Ryan Eckbo



Installation
------------


### 1. From Source

Checkout from github:

    git clone git://github.com/pnlbwh/ukftractography.git

There are 3 ways to build this project from source, as a stand alone
superbuild, against a Slicer 4 build, and as a Slicer 4 extension build (which
is more of a test than anything).


#### a) Standalone Superbuild

    cd <build-dir>
    cmake <path-to-source>
    make
    make test

#### b) Build with Slicer4

    cd <build-dir>
    cmake -DUKFTractography_SUPERBUILD:BOOL=OFF \
      -DSlicer_DIR=<path-to-Slicer4-Superbuild>/Slicer-build <path-to-source>
    make
    make test

#### c) Build as Slicer4 extension

Manual build (for the developer to test that the extension can be installed
successfully)

    mkdir s4ext_build
    cd s4ext_build
    cmake -DCMAKE_BUILD_TYPE:String=Release \
      -DSlicer_DIR:PATH=/path/to/Slicer-SuperBuild-Debug/Slicer-build ../ukftractography 
    make

Extension build, test, package and upload using ExperimentalUpload target
(https://github.com/Slicer/ExtensionsIndex#extension-build-test-package-and-upload-using-experimentalupload-target):

    cd s4ext_build
    cd ukf_tractography-build
    cmake -DMIDAS_PACKAGE_URL:STRING=http://slicer.kitware.com/midas3 -DMIDAS_PACKAGE_EMAIL:STRING=<EMAIL> -DMIDAS_PACKAGE_API_KEY:STRING=<API KEY> .
    make ExperimentalUpload


### 2. As a Slicer 4 Extension

Navigate to the Slicer Extension Manager and download `UKF Tractography` to
install it as a Slicer 4 module.



Basic Usage
-----------

### 1. As Command Line Module

The executable will be called 'UKFTractography', and can be found in the bin directory
of your build folder. 

In order to see all options run.

    ./UKFTractography --help 

In the source directory of the project you will find a shell script called 'sample_run.sh'
It should give you an idea of what a function call could look like. 


### 2. As Slicer 4 module

Navigate to the Slicer Extension Manager and download `UKF Tractography` to
install it as a Slicer 4 module.  There will be 3 modules under
`Diffusion-->Tractography`: `UKF Tractography`, `vtk2mask`, and `vtkFilter`.


Notes
-----

On a Mac, there are rounding errors that affect the accuracy of 2T FW tracts.
This explains why the 2T_FW ctest fails.
