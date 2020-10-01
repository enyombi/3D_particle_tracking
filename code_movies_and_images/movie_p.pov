//Author: Eru K.
//Objective: To draw a packing of particles with its voronoi network

#version 3.6;

// x-axis == left(-) to right(+)
// y-axis == bottom(-) to top(+)
// z-axis == foreground(-) to background(+)
//image is 1024x1024x795(width x length x depth)

/*
   +y	
   |
   | /+z
   |/
   -----+x

	   -------	
   C	 /       /|
	/_______/ |
	|   A	| | 
	|	| /
	o_______|/
(0,0,0)

o = origin = (0,0,0)
camera(C) is located at (-1024,2000,-512) wrt origin
camera is looking at pt A i.e., (512, 512,512) wrt origin

d1 = sqrt((C(1)-A(1))^2+(C(2)-A(2))^2+(C(3)-A(3))^2) //distance between C and A

      d3
   C------b
    .	  |
   d1 .	  |d2
        . |
	  A

b = (512,2000,512)
d2 = height of image - A(3) = 1024 - 512 = 512
d1^2 = d3^2 + d2^2
d3 = sqrt(d1^2-d2^2)

so to rotate in a circle centered around pt b...
loacation<b(1)+d3*cos(2*pi*clock),b(2),b(3)+d3*sin(2*pi*clock)>
*/
camera {
       location <512+(2315*cos(2*pi*clock)),2000,512+(2315*sin(2*pi*clock))>
       look_at<512,512,512>
       sky <0,1,0> //orientation of the camera, <0,1,0> means the camera is facing up and isn't slanted
}

// White background
background{rgb 1}

// Two lights with slightly different colors
light_source{<800,2000,-300> color rgb <0.77,0.75,0.75>}
light_source{<-2500,-1200,-1200> color rgb <0.38,0.40,0.40>}


//This is .pov file is made by matlab 
// Particles
union{
#include "pks30may2014_1_30mayRCP_4x_p.pov"
//	pigment{rgb 0.95} finish{reflection 0.1 specular 0.3 ambient 0.42}
}

/*
// Voronoi Network
#declare r = 8; //radius of cylinders used to draw voronoi network
union{
#include "pks30may2014_1_30mayRCP_4x_v.pov"
        pigment{rgb <1,0.4,0.45>} finish{specular 0.5 ambient 0.42}
}
*/