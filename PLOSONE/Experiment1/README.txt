The code for running Experiment 1 has the following software dependencies:
 - mtimesx (downloaded from Matlab Central and included with this code)
    > mtimesx must first be "built" before running the scripts for Experiment 1 (see documentation in the "mtimesx" folder)
    > NOTE: the file mtimesx_build.m has been modified to enable building in Matlab2014a

------------------------
Experiment 1: Comparison of Total Least Squares Point-to-Point Registration Methods

Note: test results are saved to a text file in the local directory, as well as displayed in the MATLAB command window.

Experiment 1A:
  - open App_RunRandomTrials_RigidBodyXfm.m and set:
      isoInit = 0;
      anisBothSets = 1;
  - run App_RunRandomTrials_RigidBodyXfm.m
      
Experiment 1B:
  - open App_RunRandomTrials_RigidBodyXfm.m and set:
      isoInit = 0;
      anisBothSets = 0;
  - run App_RunRandomTrials_RigidBodyXfm.m
  - open App_RunRandomTrials_RigidBodyXfm.m and set:
      isoInit = 1;
      anisBothSets = 0;
  - run App_RunRandomTrials_RigidBodyXfm.m
  
Experiment 1C:
  - run App_RunRandomTrials_Rotation.m