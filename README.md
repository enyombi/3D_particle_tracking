# 3D_particle_tracking

## Context: 

* PhD reserach on studying random packing of colloids for developing statistical physics models that predict random structure 

* references for hardcore scientists:

   * Dissertation: https://drive.google.com/file/d/0B3T0tfDPwGqXYnczMzV0dUdtZGM/view?usp=sharing
   
   * paper: https://arxiv.org/abs/1708.04702 


------

## Motivation for a general, NON-scientific audience:

(*__You CAN/SHOULD skip this. BUT if you stumbled upon this page and you're wondering what you're looking at. Please don't run away. This is for you__*).  
   
   In elementary school we learn about the phases of matter: solid, liquid, and gas.  
   
   Later on we may learn that there are theories for understanding solids, specifically crystalline solids (i.e., solid state physics/chemistry, crystallography, etc.).  A crystalline solids has an ordered, predictable atomic structure. 
   
   There's also the ideal gas law for explaining ideal gases. 
   
   Recall, an ideal gas has massless and perfectly elastic molecules w/o intermolecular forces. These molecules will fill the entire volume of the space they occupy, blah, blah, blah... 
   
   But ideal gases are just that: _ideal_. 
   
   In the real world, things attract and repel each other (including molecules).  This is most evident in phase-transitions...i.e., if you cool down a gas, then it condenses into a liquid.  And if you cool a liquid it solidifies into a solid.  
   
   So we have 2 extremes: 1) an ordered crystalline solid (real) and 2) an ideal gas (fake/fictional).  
   
   Solid state physics explains and predicts behavior of crystalline solids and the ideal gas law explains/predicts behavior of a _fictional_ ideal gas.  But there's no universal law/theory to explain everything in-between (i.e., liquids, real/non-ideal gases, supercooled liquid like glass, plastics, etc.)
   
   Therefore, let's conduct experiments to visualize and predict random arrangements/structure of molecules like those found in liquids, glasses, etc. 
   
   Then let's analyze these experiments and formulate theories to explain/predict behavior in essentially all types of matter _sans_ crystalline solids. 
   
   Molecules are too small and unwieldy to study. Instead, let's model molecular random structure after the random structure/arrangements found in grains of sand.  (I know, I know it's a stretch but we gotta start somewhere so please bare with me).
   
   _My PhD research studies random arrangements of microscopic particles of 'sand' (aka silica) and then uses statistical physics (i.e., lots of math and probabilities) to try to predict randomness.  I worked closely with physicists to build and to test numerous theories on randomness found in granular matter.  I also designed the experimental system (i.e., 3D images of fluorescent microscopic 'sand') and programmed the machine-learning algorithms used to analyze this system_. 
   
   __*This repo. contains those machine-learning algorithms!*__
   
-----

## Contents of this repository: 

programs and some of the results (i.e., .png, .mov, .m4v) of those programs used...

1) to locate >1000 particles in a 3D image (particleTracker_position.m), 

   * 'Pardon my dust'. This is both a *development*- and final-version of a machine-learning and parallel-computing regression algorithm for locating/mapping 3D images of particles  
   
2) to approximate the size of particles (particleTracker_size.m), and 

3) to identify contacts between particles (contactAnalysis.m)

* __note1__: /particleTrackMatlab contains all the functions needed to run the programs above

* __note2__: the main program w/the machine-learning and parallel-computing algorithm (i.e., particleTracker_position.m) chops up and analyzes 3D images. These files are waaaaay too big (>1GB) and not included in this repo.

## Chronological order of analysis:

1) 3D image of packing (in <font color = 'green'>__Java__</font> using ImageJ software)

    * see /movies_particleTracking/movie_packing_REAL.mov

2) 3D image processing and object-recognition/particle-tracking (in <font color = 'green'>__MATLAB__</font> + _parallel-computing_)

    * see particleTracker_position.m, particleTracker_size.m, contactAnalysis.m, and /particleTrackMatlab

3) 3D visualization (via <font color ='green'>__C/C++__</font>, <font color='green'>__BASH__</font>, POV-ray) 

   * For programs, see /code_movies_and_images
   
      * generateResults.cc --> pngTojpg.sh  --> filenameTonumber.sh --> povray for 3D rendering/visualization (i.e., packing.png, movie_packing_RECONSTRUCTED.mp4, packing_contactNetwork.pdf)
      
   * For results (i.e., reconstructed 3D images and movies), see /movies_particleTracking and /images_particleTracking
