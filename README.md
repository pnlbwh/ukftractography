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


Authors
-------
 
Yogesh Rathi (yogesh@bwh.harvard.edu), Stefan Lienhard, Yinpeng Li, Martin
Styner, Ipek Oguz, Yundi Shi, Christian Baumgartner (c.f.baumgartner@gmail.com)
Ryan Eckbo


Building the Project 
---------------------

### 1. Standalone Superbuild

    cd <build-dir>
    cmake <path-to-source>/superbuild
    make

### 2. Build with Slicer4

    cd <build-dir>
    cmake -DSlicer_DIR=<path-to-Slicer4-Superbuild>/Slicer-build <path-to-source>
    make

### 3. Download as a Slicer4 extension


TODO
----

* Get tests to pass in slicer/s4ext build
* Get `make test` to work with superbuild


-------------------------------
Below are the old instructions

Building the Project 
---------------------

In order to build make a new directory outside of the source tree called for 	
example ukf-build.

There are 3 ways to build the project

### 1. Standalone 

From your build directory run the following commands

    > cmake <path-to-source>/SuperBuild
    > make

This will download all resources and build the project for you.
Note: cmake >= 2.8.4, git, and svn are required for this

### 2. Build with Slicer 

To build against one of the slicer versions run:

    * Slicer 3

        > cmake -DSlicer3_DIR=<path-to-slicer3>/Slicer3-build <path-to-source>
        > make

    * Slicer 4

        > cmake -DSlicer_DIR=<path-to-slicer4-superbuild>/Slicer-build <path-to-source>
        > make

Note: cmake >= 2.6, and an installed version of boost >1.41.1 is required
Also, you must have a Slicer build. Only, executables are not enough.

### 3. Manual Build

You can take care of all dependancies yourself, the easiest
way to do this is using ccmake. Run ccmake from your build directory

    > ccmake <path-to-source>

Provide the Links to the builds of ITK, Teem, GenerateCLP, and Boost

    > make

Note: This is only recommended if the above two methods failed for some reason


Running the Executable
----------------------

### 1. As Command Line Module

The executable will be called 'UKFTractography', and can be found in the bin directory
of your build folder. 

In order to see all options run.

    ./UKFTractography --help 

In the source directory of the project you will find a shell script called 'sample_run.sh'
It should give you an idea of what a function call could look like. 


### 2. As Slicer 3 or 4 module

Open Slicer and in the settings add the '<path-to-build>/bin' directory. When you restart
Slicer the module will be under Diffusion->Tractography.


Everything Else
---------------

Please refer to the wiki page of this project under 
http://www.nitrc.org/plugins/mwiki/index.php/ukftractography:MainPage. 
