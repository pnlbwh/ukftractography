/**
 * \mainpage
 * \brief This module traces fibers in a DWI Volume using the multiple tensor unscented Kalman Filter methology.
 * \image html UKFloop.png "System overview illustrating relation between the neural fibers, scanner signals,\n and the unscented Kalman Filter as it is used to estimate the local model"
 * \authors Yogesh Rathi, Stefan Lienhard, Yinpeng Li, Martin Styner,
 *          Ipek Oguz, Yundi Shi, Christian Baumgartner
 *
 * \section Introduction
 * \par
 * The C++ module for Unscented Kalman Filter (UKF) Tractography enables the user
 * to easily perform multi-tensor filtered tractography on diffusion weighted MR images.
 * \par
 * The model is highly customizable. one-, two-, and three-tensor models are available
 * with, and without free water estimation. In addition, a 'full tensor model' representation
 * is available, where each of the eigenvalues is modeled sperately.
 * The UKF can be run constrained with Quadratic Programming or unconstrained.
 * \par
 * The code can be easily compiled as described below, and can be run either from command line
 * or as Slicer 3 or 4 module.
 *
 * \section Building
 *
 * In order to build make a new directory outside of the source tree called for
 * example ukf-build.
 *
 * There are 3 ways to build the project
 * <ol>
 * <li>
 *  <b>Standalone:</b> From your build directory run the following commands
 *  \code
 *    cmake <path-to-source>/SuperBuild
 *    make
 *  \endcode
 *  This will download all resources and build the project for you.\n
 *  Note: cmake >= 2.8.4, git, and svn are required for this
 * </li>
 * <li>
 *  <b>Build with Slicer:</b> To build against one of the slicer versions run
 *  <ol type='a'>
 *    <li> Slicer 3
 *       \code
 *         cmake -DSlicer3_DIR=<path-to-slicer3>/Slicer3-build <path-to-source>
 *         make
 *       \endcode
 *     </li>
 *    <li> Slicer 4
 *      \code
 *        cmake -DSlicer_DIR=<path-to-slicer4>/Slicer-build <path-to-source>
 *        make
 *      \endcode
 *    </li>
 *  </ol>
 *  Note: cmake >= 2.6, and an installed version of boost >1.41.1 is required
 *  Also, you must have a Slicer build. Only, executables are not enough.
 * </li>
 * <li> <b>Manual Build:</b> You can take care of all dependancies yourself, the easiest
 *    way to do this is using ccmake. Run ccmake from your build directory
 *    \code ccmake <path-to-source> \endcode
 *    Provide the Links to the builds of ITK, Teem, GenerateCLP, and Boost
 *    \code make \endcode
 *    Note: This is only recommended if the above two methods failed for some reason
 * </li>
 *
 * \section Running
 * \par
 * <ol>
 *   <li> <b>As Command Line Module: </b></li>
 *   The executable will be called 'UKFTractography', and can be found in the bin directory
 *   of your build folder.
 *   <br>
 *   In order to see all options run.
 *   \code
 *     ./UKFTractography --help
 *   \endcode
 *   In the source directory of the project you will find a shell script called 'sample_run.sh'
 *   It should give you an idea of what a function call could look like.
 * <li> <b>As Slicer 3 or 4 module: </b></li>
 *   Open Slicer and in the settings add the '<path-to-build>/bin' directory. When you restart
 *   Slicer the module will be under Diffusion->Tractography.
 * </ol>
 * \section Options
 * <p><b>Input/Output (IO)</b></p>
 * <p><table border=0 width=100%>
 * <tr>
 *   <td width=33%>--dwiFile <std::string></td>
 *   <td width=67%>Input DWI Image</td>
 * </tr>
 * <tr>
 *   <td width=33%>--seedsFile <std::string></td>
 *   <td width=67%>Seeds for diffusion. If not specified, full brain tractography will be
 *   performed, and the algorithm will start from every voxel in the brain
 *   mask where the Generalized Anisotropy is bigger than 0.18</td>
 * </tr>
 * <tr>
 *   <td width=33%>--labels <std::vector<int>></td>
 *   <td width=67%>A vector of the ROI labels to be used (default: 1)</td>
 * </tr>
 * <tr>
 *   <td width=33%>--maskFile <std::string></td>
 *   <td width=67%>Mask for diffusion tractography</td>
 * </tr>
 * <tr>
 *   <td width=33%>--tracts <std::string></td>
 *   <td width=67%>Tracts generated, with first tensor output</td>
 * </tr>
 * <tr>
 *   <td width=33%>--tractsWithSecondTensor <std::string></td>
 *   <td width=67%>Tracts generated with second tensor output (if there is one)</td>
 * </tr>
 * </table></p>
 * <p><b>Seeding Options</b></p>
 * <p><table border=0 width=100%>
 * <tr>
 *   <td width=33%>--seedsPerVoxel <int></td>
 *   <td width=67%>Number of seeds per voxel (default: 1)</td>
 * </tr>
 * </table></p>
 * <p><b>Model Options</b></p>
 * <p><table border=0 width=100%>
 * <tr>
 *   <td width=33%>--numTensor <1|2|3></td>
 *   <td width=67%>Number of tensors used (default: 2)</td>
 * </tr>
 * <tr>
 *   <td width=33%>--simpleTensorModel</td>
 *   <td width=67%>Whether to use the simple tensor model. If unchecked, use the full tensor model</td>
 * </tr>
 * <tr>
 *   <td width=33%>--freeWater</td>
 *   <td width=67%>Adds a term for free water difusion to the model. If checked, the 1T
 *   simple model is forced. </td>
 * </tr>
 * </table></p>
 * <p><b>Stopping Criteria</b></p>
 * <p><table border=0 width=100%>
 * <tr>
 *   <td width=33%>--minFA <ukfPrecisionType></td>
 *   <td width=67%>Abort the tractography when the Fractional Anisotropy is less than
 *   this value (default: 0.15)</td>
 * </tr>
 * <tr>
 *   <td width=33%>--minGA <ukfPrecisionType></td>
 *   <td width=67%>Abort the tractography when the Generalized Anisotropy is less than
 *   this value (default: 0.1)</td>
 * </tr>
 * </table></p>
 * <p><b>UKFFiber Scalar Fields</b></p>
 * <p><table border=0>
 * <tr>
 *   <td width=33%>--recordFA</td>
 *   <td width=67%>Whether to store FA. Attaches field 'FA', and 'FA2' for 2-tensor case
 *   to fiber. </td>
 * </tr>
 * <tr>
 *   <td width=33%>--recordNMSE</td>
 *   <td width=67%>Whether to store NMSE. Attaches field 'NMSE' to fiber. </td>
 * </tr>
 * <tr>
 *   <td width=33%>--recordTrace</td>
 *   <td width=67%> Whether to store Trace. Attaches field 'Trace', and 'Trace2' for
 *    2-tensor case to fiber. </td>
 * </tr>
 * <tr>
 *   <td width=33%>--recordFreeWater</td>
 *   <td width=67%>Whether to store the fraction of free water. Attaches field
     'FreeWater' to fiber. </td>
 * </tr>
 * <tr>
 *   <td width=33%>--recordState</td>
 *   <td width=67%>Whether to attach the states to the fiber. Will generate field
 *   'state'. </td>
 * </tr>
 * <tr>
 *   <td width=33%>--recordCovariance</td>
 *   <td width=67%>Whether to store the covariance. Will generate field 'covariance' in
     fiber. </td>
 * </tr>
 * <tr>
 *   <td width=33%>--recordTensors</td>
 *   <td width=67%>Recording the tensors enables Slicer to color the fiber bundles by FA,
 *    orientation, and so on. The fields will be called 'TensorN', where N
 *    is the tensor number. </td>
 * </tr>
 * </table></p>
 * <p><b>Advanced Options</b></p>
 * <p><table border=0>
 * <tr>
 *   <td width=33%>--numThreads <int></td>
 *   <td width=67%>Number of threads used during compuation. Set to the number of cores on your workstation for
 *       optimal speed. If left undefined boost will figure out the number of cores, and hence threads, during runtime.</td>
 * </tr>
 * <tr>
 *   <td width=33%>--normalizedDWIData</td>
 *   <td width=67%>Whether the DWI input data is already normalized</td>
 * </tr>
 * <tr>
 *   <td width=33%>--stepLength <ukfPrecisionType></td>
 *   <td width=67%>Step length of tractography, in millimeters. (If not set, defined during runtime)</td>
 * </tr>
 * <tr>
 *   <td width=33%>--recordLength <ukfPrecisionType></td>
 *   <td width=67%>Record length of tractography, in millimeters. (If not set, defined during runtime)</td>
 * </tr>
 * <tr>
 *   <td width=33%>--sigmaSignal <ukfPrecisionType></td>
 *   <td width=67%>Sigma for gaussian interpolation of signal. (If not set, defined during runtime)</td>
 * </tr>
 * <tr>
 *   <td width=33%>--weightsOnTensors <std::vector<ukfPrecisionType>></td>
 *   <td width=67%>Weights on different tensors when using multiple tensors. There must be one weight for each tensor, and the weights must sum up to 1. Defaults to equally weighted.</td>
 * </tr>
 * <tr>
 *   <td width=33%>--maxHalfFiberLength <ukfPrecisionType></td>
 *   <td width=67%>The max length limit of the half fibers generated during tractography. Here the fiber is "half" because the tractography goes in only one direction from one seed point at a time (default: 10000)</td>
 * </tr>
 * <tr>
 *   <td width=33%>--seedFALimit <ukfPrecisionType></td>
 *   <td width=67%>Seed points whose FA are below this value are excluded. (If not set, defined during runtime)</td>
 * </tr>
 * <tr>
 *   <td width=33%>--Qm <ukfPrecisionType></td>
 *   <td width=67%>Process noise for angles/direction. (If not set, defined during runtime)</td>
 * </tr>
 * <tr>
 *   <td width=33%>--Ql <ukfPrecisionType></td>
 *   <td width=67%>Process noise for eigenvalues. (If not set, defined during runtime)</td>
 * </tr>
 * <tr>
 *   <td width=33%>--Qw <ukfPrecisionType></td>
 *   <td width=67%>Process noise for free water weights, ignored if no free water estimation. (If not set, defined during runtime)</td>
 * </tr>
 * <tr>
 *   <td width=33%>--Rs <ukfPrecisionType></td>
 *   <td width=67%>Measurement noise. (If not set, defined during runtime)</td>
 * </tr>
 * <tr>
 *   <td width=33%>--maxBranchingAngle <ukfPrecisionType></td>
 *   <td width=67%>Maximum branching angle, in degrees. When using multiple tensors, a new branch will be created when the tensors' major directions
 *   form an angle between (minBranchingAngle, maxBranchingAngle). Branching is supressed     when this maxBranchingAngle is set to 0.0</td>
 * </tr>
 * <tr>
 *   <td width=33%>--minBranchingAngle <ukfPrecisionType></td>
 *   <td width=67%>Minimum branching angle, in degrees. When using multiple tensors, a new branch will be created when the tensors' major directions form an angle between (minBranchingAngle, maxBranchingAngle)</td>
 * </tr>
 * </table></p>
 * <p><b>Additional Output Options</b></p>
 * <p><table border=0>
 * <tr>
 *   <td width=33%>--noTransformPosition</td>
 *   <td width=67%>Don't transform Points back from ijk->RAS when writing the output fiber</td>
 * </tr>
 * <tr>
 *   <td width=33%>--branchesOnly</td>
 *   <td width=67%>Only output branches, ignore the primary tracts</td>
 * </tr>
 * <tr>
 *   <td width=33%>--storeGlyphs</td>
 *   <td width=67%>Store tensors' main directions as two-point lines in a separate file named glyphs_{tracts}. When using multiple tensors, only the major tensors' main directions are stored</td>
 * </tr>
 * <tr>
 *   <td width=33%>--outputNormalizedDWIData</td>
 *   <td width=67%>Whether to output the DWI after normalization (i.e. preprocessing)</td>
 * </tr>
 * </table></p>
 * <p><b>Tractography unrelated options</b></p>
 * <p><table border=0 width=100%>
 * <tr>
 *   <td width=33%>--returnparameterfile <std::string></td>
 *   <td width=67%>Filename in which to write simple return parameters (int, float,
 *    int-vector, etc.) as opposed to bulk return parameters (image,
 *    geometry, transform, measurement, table).</td>
 * </tr>
 * <tr>
 *   <td width=33%>--processinformationaddress <std::string></td>
 *   <td width=67%>Address of a structure to store process information (progress, abort,
 *   etc.). (default: 0)</td>
 * </tr>
 * <tr>
 *   <td width=33%>--xml</td>
 *   <td width=67%>Produce xml description of command line arguments.</td>
 * </tr>
 * <tr>
 *   <td width=33%>--echo</td>
 *   <td width=67%>Echo the command line arguments.</td>
 * </tr>
 * <tr>
 *   <td width=33%>--, --ignore_rest</td>
 *   <td width=67%>Ignores the rest of the labeled arguments following this flag.</td>
 * </tr>
 * <tr>
 *   <td width=33%>--version</td>
 *   <td width=67%>Displays version information and exits.</td>
 * </tr>
 * <tr>
 *   <td width=33%>-h, --help</td>
 *   <td width=67%>Displays usage information and exits.</td>
 * </tr>
 * </table></p>
 * \section Dependencies
 * \par Boost
 * The Boost C++ Libraries are a collection of free libraries that extend the functionality of C++. In this project boost is used for multithreading, and the progress bar.
 * \par teem
 * Teem is a coordinated group of libraries for representing, processing, and visualising scientific raster data. We use it to read, and process
 * The input data of type NRRD.
 * \par ITK
 * Insight Segmentation and Registration Toolkit (ITK).  ITK is an open-source, cross-platform system that provides developers
 * with an extensive suite of software tools for image analysis. In this application we only use the underlying Eigen libraries for most of the linear algebra processing.
 * It replaces functionality that was done with LAPACK before.
 * \par "Slicer Execution Model"
 * The Execution Model provides a simple mechanism for incorporating command line programs as Slicer modules. These command line modules
 * are self-describing, emitting an XML description of its command line arguments. Slicer uses this XML description to construct a GUI for the module. We use it to make our application
 * into a Command Line Module that can be run from Slicer, or from a Consolse without being dependand on Slicer.
 * \section Screenshots
 * \image html 2TFW_0.png "Tracing fibers through the anterior limb of the internal capsule"
 * \todo Doxygen can elegantly be generated by cmake as describe e.g. here http://www.cmake.org/pipermail/cmake/2006-August/010794.html.
*/

// Defines content of the main documentation page
