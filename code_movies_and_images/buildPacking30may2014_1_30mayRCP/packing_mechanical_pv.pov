//Author: Eru K.
//Objective: To draw a packing of particles with its voronoi network

#version 3.6;

// x-axis == left(-) to right(+)
// y-axis == bottom(-) to top(+)
// z-axis == foreground(-) to background(+)
//image is 1024x1024x795(width x length x depth)
camera {
	location <-1024,2000,-512>
	look_at<512,512,512>
	sky <0,1,0> //orientation of the camera, <0,1,0> means the camera is facing up and isn't slanted
	right x*(image_width/image_height)
}

// White background
background{rgb 1}

// Two lights with slightly different colors
light_source{<800,2000,-300> color rgb <0.77,0.75,0.75>}
light_source{<-2500,-1200,-1200> color rgb <0.38,0.40,0.40>}


//This is .pov file is made by matlab 
// Particles
union{
#include "pks30may2014_1_30mayRCP_4x_mechanical_p.pov"
//	 finish{reflection 0.1 specular 0.3 ambient 0.42}
}


// Voronoi Network
#declare r = 8; //radius of cylinders used to draw voronoi network
union{
#include "pks30may2014_1_30mayRCP_4x_v.pov"
	 pigment{rgb <0.90,0.95,1>} finish{specular 1 ambient 0.42}
}
