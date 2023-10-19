Philip Mocz, Xinyi Guo
Harvard University
CS 207, Spring 2014

== Prompt ==

  Each team should submit a working demonstration of their final project. With the
  code, submit a README_FINAL.txt that details the features of your project. This should
  include how we should run the code, any run-time commands or interaction, any in-
  teresting things we should try, etc.

== Description ==

Our code extends the Constrained Transport (CT) algorithm for magnetohydrodynamics (MHD) from a static Cartesian grid to an unstructured triangular mesh. The code solves simple MHD test problems in a numerically robust way that maintains the divergence of the magnetic fields to 0 to machine precision. 

The code visualizes the density field using colored triangles (using our extended SDLViewer).

The code also allows the user to interact with the simulation by setting initial conditions and injecting energy into the simulation in real time (using listeners addeded to our extended SDLViewer). 

The code simulates 3 classical test problems:
  [1] Implosion
  [2] BlastWave
  [3] OrszagTang (decaying supersonic turbulence)
and also has an interactive mode to set initial conditions for an Eulerian fluid
  [4] HydroInteractive

All modes support the injection of energy by right-clicking a point in the domain.

== Main Files to Review ==

  mhd.cpp
    main code for solving the MHD equations and setting up initial conditions and listeners.

  CS207/SDLViewer.hpp
    our extended SDLViewer class that supports plotting of shaded triangles and handling user-defined listeners

  Mesh.hpp
   our Mesh class designed in HW4

  final_report/writeup.pdf
    report describing our project

  COLLAB_EVAL.txt
    description of collaboration


== How to Run the Code ==

  $ ./mhd NODES_FILE TRIANGLE_FILE IC_STRING

  right-click to inject energy in the domain in real time

  [1] Implosion

      $ ./mhd data/tub3.nodes data/tub3.tris Implosion


  [2] BlastWave

      $ ./mhd data/tub3.nodes data/tub3.tris BlastWave 


  [3] OrszagTang (decaying supersonic turbulence)

      $ ./mhd data/pond3.nodes data/pond3.tris OrszagTang 


  [4] HydroInteractive

      $ ./mhd data/tub3.nodes data/tub3.tris HydroInteractive 

      Interactive mode will prompt the user to left-click 2 points in the screen to set a fluid-phase dividing line
      The user is then prompted in the terminal to set the density contrast between the two regions
      The user is also asked to set the pressure contrast and also add optional velocity shear
