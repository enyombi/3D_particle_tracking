# 3D_particle_tracking

__Context__: PhD research (<font color = 'blue'> https://arxiv.org/abs/1708.04702 </font>) 
studying random packing of collids for developing statistical physics models that predict random structure 

__Content of this repository__: programs and some of the results (i.e., .png, .mov, .m4v) of those programs used...

1) to locate >1000 particles in a 3D image, 
2) to approximate the size of particles, and 
3) to identify contacts between particles

__Order__:

1) 3D image of packing (in <font color = 'green'>__Java__</font> using ImageJ software)

  * see movie_packing_REAL.mov

2) 3D image processing and object-recognition/particle-tracking (in <font color = 'green'>__MATLAB__</font> + _parallel-computing_)

  * see particleTracker_position.m, particleTracker_size.m, contactAnalysis.m, and /particleTrackMatlab

3) 3D visualization (via <font color ='green'>__C__</font>, <font color='green'>__BASH__</font>, POV-ray) 

* generateResults.cc --> pngTojpg.sh  --> filenameTonumber.sh --> povray for 3D rendering/visualization (i.e., packing.png, movie_packing_RECONSTRUCTED.mp4, packing_contactNetwork.pdf)
