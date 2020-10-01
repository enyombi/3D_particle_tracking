//Author: Eru K.
//Objective: To draw a packing of particles with its voronoi network

#version 3.6;

// Right-handed coordinate system in which the z-axis points upwards
camera {
	location <-1024,2000,-512>
	look_at<512,512,512>
}

// White background
background{rgb 1}

// Two lights with slightly different colors
light_source{<800,2000,-300> color rgb <0.77,0.75,0.75>}
light_source{<-2500,-1200,-1200> color rgb <0.38,0.40,0.40>}

/*
//This is .pov file is made by matlab 
// Particles
union{
#include "pks30may2014_1_30mayRCP_4x_network_p.pov"
//	pigment{rgb 0.95} finish{reflection 0.1 specular 0.3 ambient 0.42}
}
*/

// Contact Network according to fluorophore exclusion
#declare r = 8; //radius of cylinders used to draw voronoi network
union{
#include "pks30may2014_1_30mayRCP_4x_mechanical_network.pov"
        pigment{rgb <0.90,0.95,1>} finish{specular 1 ambient 0.42}
}
