#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "ml6.h"
#include "display.h"
#include "draw.h"
#include "matrix.h"
#include "gmath.h"

int order_points(int d01, int d02, int d12){
  int ans;
  if(d01>0){
    if(d02<0){
      //+--:201
      ans = 201;
    }else{
      if(d12>0){
	//+++:012
	ans = 12;
      }else{
	//++-:021
	//++0:021
	//+0-:021
	ans = 21;
      }
    }
  }else{
    if(d02>0){
      //-++:102
      //0++:102
      ans = 102;
    }else{
      if(d12<0){
	//---:210
	//0--:210
        ans = 210;
      }else{
	//--+:120
	//-0+:120
	//--0:120
	ans = 120;
      }
    }
  }
  return ans;
}


/*======== void add_polygon() ==========
Inputs:   struct matrix *surfaces
         double x0
         double y0
         double z0
         double x1
         double y1
         double z1
         double x2
         double y2
         double z2  
Returns: 
Adds the vertices (x0, y0, z0), (x1, y1, z1)
and (x2, y2, z2) to the polygon matrix. They
define a single triangle surface.

04/16/13 13:05:59
jdyrlandweaver
====================*/
void add_polygon( struct matrix *polygons, 
		  double x0, double y0, double z0, 
		  double x1, double y1, double z1, 
		  double x2, double y2, double z2 ) {
  add_point(polygons, x0, y0, z0);
  add_point(polygons, x1, y1, z1);
  add_point(polygons, x2, y2, z2);
}

/*======== void draw_polygons() ==========
Inputs:   struct matrix *polygons
          screen s
          color c  
Returns: 
Goes through polygons 3 points at a time, drawing 
lines connecting each points to create bounding
triangles

04/16/13 13:13:27
jdyrlandweaver
====================*/
void draw_polygons( struct matrix *polygons, screen s, color c ) {
  
  int i;
  int total;
  int count;
  total = 0;
  count = 0;
  for( i=0; i < polygons->lastcol-2; i+=3 ) {

    if ( calculate_dot( polygons, i ) < 0 ) {
      draw_line( polygons->m[0][i],
		 polygons->m[1][i],
		 polygons->m[0][i+1],
		 polygons->m[1][i+1],
		 s, c);
      draw_line( polygons->m[0][i+1],
		 polygons->m[1][i+1],
		 polygons->m[0][i+2],
		 polygons->m[1][i+2],
		 s, c);
      draw_line( polygons->m[0][i+2],
		 polygons->m[1][i+2],
		 polygons->m[0][i],
		 polygons->m[1][i],
		 s, c);

      /////////////////SCANLINE CONVERSION//////////////////

      ////Designate points TOP, BOTTOM, and MIDDLE////
      int T, M, B, d01, d02, d12, ans;

      d01 = polygons->m[1][i] - polygons->m[1][i+1];
      d02 = polygons->m[1][i] -  polygons->m[1][i+2];
      d12 = polygons->m[1][i+1] -  polygons->m[1][i+2];
     
      ans = order_points(d01, d02, d12);
      T = ans/100;
      M = (ans%100)/10;
      B = (ans%10);

      while(!(polygons->m[1][i+T] >= polygons->m[1][i+M] && polygons->m[1][i+M] >= polygons->m[1][i+B])){
	if(polygons->m[1][i+T] < polygons->m[1][i+M]){
	  int temp = T;
	  T = M;
	  M = temp;
	}
	if(polygons->m[1][i+M] < polygons->m[1][i+B]){
	  int temp = M;
	  M = B;
	  B = temp;
	}
      }
      
      if(polygons->m[1][i+T] == polygons->m[1][i+B]){
	printf("things are going wrong\n");
      }
      if(polygons->m[1][i+T] == polygons->m[1][i+M]){
	if(polygons->m[0][i+T] >  polygons->m[0][i+M]){
	  int temp = T;
	  T = M;
	  M = temp;
	}
      }
      if(polygons->m[1][i+M] == polygons->m[1][i+B]){
	if(polygons->m[0][i+B] >  polygons->m[0][i+M]){
	  int temp = B;
	  B = M;
	  M = temp;
	}
      }

      //printf("T:%d\nM:%d\nB:%d\n",T,M,B);
      total++;
      if(polygons->m[1][i+T] >= polygons->m[1][i+M] && polygons->m[1][i+M] >= polygons->m[1][i+B]){
	//printf("Correct!\n");
	count++;
      }else if(polygons->m[1][i+T]==polygons->m[1][i+B]){
	printf("horizontal line\n");
      }
      //printf("%d/%d\n",count,total);

      ///////////Set x0, y0, x1, y1, and d0////////////
      double x0, y0, x1, y1, d0, d1;
      
      x0 = polygons->m[0][i+B];
      y0 = polygons->m[1][i+B];
      x1 = polygons->m[0][i+B];
      y1 = polygons->m[1][i+B];

      //printf("Top:(%f,%f)\nMiddle:(%f,%f)\nBottom:(%f,%f)\n",polygons->m[0][i+T],polygons->m[1][i+T],polygons->m[0][i+M],polygons->m[1][i+M],x0,y0);
 
      d0 = (polygons->m[0][i+T] - polygons->m[0][i+B])/(polygons->m[1][i+T] - polygons->m[1][i+B]);
      //////////Start moving the points///////////////
      if( (polygons->m[1][i+B] < polygons->m[1][i+M]) && (polygons->m[1][i+M] < polygons->m[1][i+T]) ){
	d1 = (polygons->m[0][i+M] - polygons->m[0][i+B])/(polygons->m[1][i+M] - polygons->m[1][i+B]);
	while(y0 <= polygons->m[1][i+M]){
	  int colors[5] = {0,64,128,192,255};
	  c.red = colors[i%5];
	  c.blue = colors[(i+1)%5];
	  c.green = colors[(i-1)%5];
	  x0+=d0;
	  x1+=d1;
	  y0++;
	  y1++;
	 
	  draw_line(x0,y0,x1,y1,s,c);
	}
	d1 = (polygons->m[0][i+T] - polygons->m[0][i+M])/(polygons->m[1][i+T] - polygons->m[1][i+M]);
	
	int colors[5] = {0,64,128,192,255};
	while(y0 < polygons->m[1][i+T]){
	  c.red = colors[i%5];
	  c.blue = colors[(i+1)%5];
	  c.green = colors[(i-1)%5];
	  x0+=d0;
	  x1+=d1;
	  y0++;
	  y1++;
	 
	  draw_line(x0,y0,x1,y1,s,c);
	}
      }else{
	if(polygons->m[1][i+B] == polygons->m[1][i+M]){
	  //printf("T>M==B\n");
	  //x1 is on MT
	  x1 = polygons->m[0][i+M];
	  y1 = polygons->m[1][i+M];
	  d1 = (polygons->m[0][i+T] - polygons->m[0][i+M])/(polygons->m[1][i+T] - polygons->m[1][i+M]);
	  //printf("d1 is %f\n",d1);
	}else if(polygons->m[1][i+T] == polygons->m[1][i+M]){
	  //printf("T==M>B\n");
	  //x1 is on BM
	  d1 = (polygons->m[0][i+M] - polygons->m[0][i+B])/(polygons->m[1][i+M] - polygons->m[1][i+B]);
	  //printf("d1 is %f\n",d1);
	}
	while(y0 < polygons->m[1][i+T]){
	  int colors[5] = {0,64,128,192,255};
	  c.red = colors[i%5];
	  c.blue = colors[(i+1)%5];
	  c.green = colors[(i-1)%5];
	  x0+=d0;
	  x1+=d1;
	  y0++;
	  y1++;
	  int dif01,dif02,dif12,ans,Rx,Mx,Lx;
	  dif01 = polygons->m[0][i] - polygons->m[0][i+1];
	  dif02 = polygons->m[0][i] - polygons->m[0][i+2];
	  dif12 = polygons->m[0][i+1] - polygons->m[0][i+2];
	  ans = order_points(dif01,dif02,dif12);
	  Rx = ans/100;
	  Mx = (ans%100)/10;
	  Lx = ans%10;
	  
	  if(x0<polygons->m[0][i+Lx] || x0>polygons->m[0][i+Rx]){
	    /*
	    printf("x0 is out of range. d0 is %f. Lx is %f. Rx is %f. x0 is %f\n",d0,polygons->m[0][i+Lx],polygons->m[0][i+Rx],x0);
	    	printf("   Top: (%f,%f)\n",polygons->m[0][i+T],polygons->m[1][i+T]);
		printf("Middle: (%f,%f)\n",polygons->m[0][i+M],polygons->m[1][i+M]);
		printf("Bottom: (%f,%f)\n\n",polygons->m[0][i+B],polygons->m[1][i+B]);
	    */
	    //printf("x0 is out of range\n");
	    if(polygons->m[1][i+T] < polygons->m[1][i+M]){
	      printf("T < M\n");
	    }
	    if(polygons->m[1][i+T] < polygons->m[1][i+B]){
	      printf("T < B\n");
	    }
	    if(polygons->m[1][i+M] < polygons->m[1][i+B]){
	      printf("M < B\n");
	    }
	    if(polygons->m[1][i+T] == polygons->m[1][i+B]){
	      printf("T == B\n");
	    }
	  }
	  if(x1<polygons->m[0][i+Lx] || x1>polygons->m[0][i+Rx]){

	    /*
	    printf("x1 is out of range. d1 is %f. Lx is %f. Rx is %f. x1 is %f\n",d1,polygons->m[0][i+Lx],polygons->m[0][i+Rx],x1);
	    printf("   Top: (%f,%f)\n",polygons->m[0][i+T],polygons->m[1][i+T]);
	    printf("Middle: (%f,%f)\n",polygons->m[0][i+M],polygons->m[1][i+M]);
	    printf("Bottom: (%f,%f)\n\n",polygons->m[0][i+B],polygons->m[1][i+B]);
	    */
	    printf("x1 is out of range\n");
	    if(polygons->m[1][i+T] < polygons->m[1][i+M]){
	      printf("T < M\n");
	    }
	    if(polygons->m[1][i+T] < polygons->m[1][i+B]){
	      printf("T < B\n");
	    }
	    if(polygons->m[1][i+M] < polygons->m[1][i+B]){
	      printf("M < B\n");
	    }
	    if(polygons->m[1][i+T] == polygons->m[1][i+B]){
	      printf("T == B\n");
	    }
	  }
	  
	  draw_line(x0,y0,x1,y1,s,c);
	}
      }
      ////////////////////////////////////////////////
    }
  }
}


/*======== void add_sphere() ==========
  Inputs:   struct matrix * points
            double cx
	    double cy
	    double r
	    double step  
  Returns: 

  adds all the points for a sphere with center 
  (cx, cy) and radius r.

  should call generate_sphere to create the
  necessary points

  jdyrlandweaver
  ====================*/
void add_sphere( struct matrix * points, 
		 double cx, double cy, double cz, double r, 
		 int step ) {

  struct matrix * temp;
  int lat, longt;
  int index;
  int num_steps, num_points;
  double px0, px1, px2, px3;
  double py0, py1, py2, py3;
  double pz0, pz1, pz2, pz3;

  num_steps = MAX_STEPS / step;
  num_points = num_steps * (num_steps + 1);
  
  temp = new_matrix( 4, num_points);
  //generate the points on the sphere
  generate_sphere( temp, cx, cy, cz, r, step );

  int latStop, longStop, latStart, longStart;
  latStart = 0;
  latStop = num_steps;
  longStart = 0;
  longStop = num_steps;

  num_steps++;

  for ( lat = latStart; lat < latStop; lat++ ) {
    for ( longt = longStart; longt < longStop; longt++ ) {
      
      index = lat * num_steps + longt;

      px0 = temp->m[0][ index ];
      py0 = temp->m[1][ index ];
      pz0 = temp->m[2][ index ];
      
      px1 = temp->m[0][ (index + num_steps) % num_points ];
      py1 = temp->m[1][ (index + num_steps) % num_points ];
      pz1 = temp->m[2][ (index + num_steps) % num_points ];

      px3 = temp->m[0][ index + 1 ];
      py3 = temp->m[1][ index + 1 ];
      pz3 = temp->m[2][ index + 1 ];

      if (longt != longStop - 1) {
	px2 = temp->m[0][ (index + num_steps + 1) % num_points ];
	py2 = temp->m[1][ (index + num_steps + 1) % num_points ];
	pz2 = temp->m[2][ (index + num_steps + 1) % num_points ];
      }
      else {
	px2 = temp->m[0][ (index + 1) % num_points ];
	py2 = temp->m[1][ (index + 1) % num_points ];
	pz2 = temp->m[2][ (index + 1) % num_points ];
      }

      if (longt != 0)
	add_polygon( points, px0, py0, pz0, px1, py1, pz1, px2, py2, pz2 );
      if (longt != longStop - 1)
	add_polygon( points, px2, py2, pz2, px3, py3, pz3, px0, py0, pz0 );
    }
  }
}

/*======== void generate_sphere() ==========
  Inputs:   struct matrix * points
            double cx
	    double cy
	    double r
	    double step  
  Returns: 

  Generates all the points along the surface of a 
  sphere with center (cx, cy) and radius r

  Adds these points to the matrix parameter

  03/22/12 11:30:26
  jdyrlandweaver
  ====================*/
void generate_sphere( struct matrix * points, 
		      double cx, double cy, double cz, double r, 
		      int step ) {


  int circle, rotation;
  double x, y, z, circ, rot;

  int rotStart = step * 0;
  int rotStop = MAX_STEPS;
  int circStart = step * 0;
  int circStop = MAX_STEPS;
  
  for ( rotation = rotStart; rotation < rotStop; rotation += step ) {
    rot = (double)rotation / MAX_STEPS;
    for ( circle = circStart; circle <= circStop; circle+= step ) {

      circ = (double)circle / MAX_STEPS;
      x = r * cos( M_PI * circ ) + cx;
      y = r * sin( M_PI * circ ) *
	cos( 2 * M_PI * rot ) + cy;
      z = r * sin( M_PI * circ ) *
	sin( 2 * M_PI * rot ) + cz;

      add_point( points, x, y, z);
    }
  }
}    


/*======== void add_torus() ==========
  Inputs:   struct matrix * points
            double cx
	    double cy
	    double r1
	    double r2
	    double step  
  Returns: 

  adds all the points required to make a torus
  with center (cx, cy) and radii r1 and r2.

  should call generate_torus to create the
  necessary points

  03/22/12 13:34:03
  jdyrlandweaver
  ====================*/
void add_torus( struct matrix * points, 
		double cx, double cy, double cz, double r1, double r2, 
		int step ) {

  struct matrix * temp;
  int lat, longt;
  int index;
  int num_steps;
  
  num_steps = MAX_STEPS / step;

  temp = new_matrix( 4, num_steps * num_steps );
  //generate the points on the torus
  generate_torus( temp, cx, cy, cz, r1, r2, step );
  int num_points = temp->lastcol;

  int latStop, longtStop, latStart, longStart;
  latStart = 0;
  longStart = 0;
  latStop = num_steps;
  longtStop = num_steps;
  for ( lat = latStart; lat < latStop; lat++ )
    for ( longt = longStart; longt < longtStop; longt++ ) {

      index = lat * num_steps + longt;

      if ( longt != num_steps-1) {
	add_polygon( points, temp->m[0][index],
		     temp->m[1][index],
		     temp->m[2][index],
		     temp->m[0][(index+num_steps+1) % num_points],
		     temp->m[1][(index+num_steps+1) % num_points],
		     temp->m[2][(index+num_steps+1) % num_points],
		     temp->m[0][index+1],
		     temp->m[1][index+1],
		     temp->m[2][index+1] );
	add_polygon( points, temp->m[0][index],
		     temp->m[1][index],
		     temp->m[2][index],
		     temp->m[0][(index+num_steps) % num_points],
		     temp->m[1][(index+num_steps) % num_points],
		     temp->m[2][(index+num_steps) % num_points],
		     temp->m[0][(index+num_steps) % num_points + 1],
		     temp->m[1][(index+num_steps) % num_points + 1],
		     temp->m[2][(index+num_steps) % num_points + 1]);
      }
      else {
	add_polygon( points, temp->m[0][index],
		     temp->m[1][index],
		     temp->m[2][index],
		     temp->m[0][(index+1) % num_points],
		     temp->m[1][(index+1) % num_points],
		     temp->m[2][(index+1) % num_points],
		     temp->m[0][index+1-num_steps],
		     temp->m[1][index+1-num_steps],
		     temp->m[2][index+1-num_steps] );
	add_polygon( points, temp->m[0][index],
		     temp->m[1][index],
		     temp->m[2][index],
		     temp->m[0][(index+num_steps) % num_points],
		     temp->m[1][(index+num_steps) % num_points],
		     temp->m[2][(index+num_steps) % num_points],
		     temp->m[0][(index+1) % num_points],
		     temp->m[1][(index+1) % num_points],
		     temp->m[2][(index+1) % num_points]);
      }

    }
}

/*======== void generate_torus() ==========
  Inputs:   struct matrix * points
            double cx
	    double cy
	    double r
	    double step  
  Returns: 

  Generates all the points along the surface of a 
  tarus with center (cx, cy) and radii r1 and r2

  Adds these points to the matrix parameter

  03/22/12 11:30:26
  jdyrlandweaver
  ====================*/
void generate_torus( struct matrix * points, 
		     double cx, double cy, double cz, double r1, double r2, 
		     int step ) {

  double x, y, z, circ, rot;
  int circle, rotation;

  double rotStart = step * 0;
  double rotStop = MAX_STEPS;
  double circStart = step * 0;
  double circStop = MAX_STEPS;

  for ( rotation = rotStart; rotation < rotStop; rotation += step ) {

    rot = (double)rotation / MAX_STEPS;
    for ( circle = circStart; circle < circStop; circle+= step ) {

      circ = (double)circle / MAX_STEPS;
      x = cos( 2 * M_PI * rot ) *
	( r1 * cos( 2 * M_PI * circ ) + r2 ) + cx;
      y = r1 * sin( 2 * M_PI * circ ) + cy;
      z = sin( 2 * M_PI * rot ) *
	( r1 * cos( 2 * M_PI * circ ) + r2 ) + cz;

      add_point( points, x, y, z );
    }
  }
}

/*======== void add_box() ==========
  Inputs:   struct matrix * points
            double x
	    double y
	    double z
	    double width
	    double height
	    double depth
  Returns: 

  add the points for a rectagular prism whose 
  upper-left corner is (x, y, z) with width, 
  height and depth dimensions.

  jdyrlandweaver
  ====================*/
void add_box( struct matrix * polygons,
	      double x, double y, double z,
	      double width, double height, double depth ) {

  double x2, y2, z2;
  x2 = x + width;
  y2 = y - height;
  z2 = z - depth;
  //front
  add_polygon( polygons, 
	       x, y, z, 
	       x, y2, z,
	       x2, y2, z);
  add_polygon( polygons, 
	       x2, y2, z, 
	       x2, y, z,
	       x, y, z);
  //back
  add_polygon( polygons, 
	       x2, y, z2, 
	       x2, y2, z2,
	       x, y2, z2);
  add_polygon( polygons, 
	       x, y2, z2, 
	       x, y, z2,
	       x2, y, z2);
  //top
  add_polygon( polygons, 
	       x, y, z2, 
	       x, y, z,
	       x2, y, z);
  add_polygon( polygons, 
	       x2, y, z, 
	       x2, y, z2,
	       x, y, z2);
  //bottom
  add_polygon( polygons, 
	       x2, y2, z2, 
	       x2, y2, z,
	       x, y2, z);
  add_polygon( polygons, 
	       x, y2, z, 
	       x, y2, z2,
	       x2, y2, z2);
  //right side
  add_polygon( polygons, 
	       x2, y, z, 
	       x2, y2, z,
	       x2, y2, z2);
  add_polygon( polygons, 
	       x2, y2, z2, 
	       x2, y, z2,
	       x2, y, z);
  //left side
  add_polygon( polygons, 
	       x, y, z2, 
	       x, y2, z2,
	       x, y2, z);
  add_polygon( polygons, 
	       x, y2, z, 
	       x, y, z,
	       x, y, z2); 
}
  
/*======== void add_circle() ==========
  Inputs:   struct matrix * points
            double cx
	    double cy
	    double y
	    double step  
  Returns: 


  03/16/12 19:53:52
  jdyrlandweaver
  ====================*/
void add_circle( struct matrix * points, 
		 double cx, double cy, 
		 double r, double step ) {
  
  double x0, y0, x, y, t;
  
  x0 = cx + r;
  y0 = cy;

  for ( t = step; t <= 1; t+= step ) {
    
    x = r * cos( 2 * M_PI * t ) + cx;
    y = r * sin( 2 * M_PI * t ) + cy;
    
    add_edge( points, x0, y0, 0, x, y, 0 );
    x0 = x;
    y0 = y;
  }

  add_edge( points, x0, y0, 0, cx + r, cy, 0 );
}

/*======== void add_curve() ==========
Inputs:   struct matrix *points
         double x0
         double y0
         double x1
         double y1
         double x2
         double y2
         double x3
         double y3
         double step
         int type  
Returns: 

Adds the curve bounded by the 4 points passsed as parameters
of type specified in type (see matrix.h for curve type constants)
to the matrix points

03/16/12 15:24:25
jdyrlandweaver
====================*/
void add_curve( struct matrix *points, 
		double x0, double y0, 
		double x1, double y1, 
		double x2, double y2, 
		double x3, double y3, 
		double step, int type ) {

  double x, y, t;
  struct matrix * xcoefs;
  struct matrix * ycoefs;
  
  //generate the coeficients
  if ( type == BEZIER_MODE ) {
    ycoefs = generate_curve_coefs(y0, y1, y2, y3, BEZIER_MODE);
    xcoefs = generate_curve_coefs(x0, x1, x2, x3, BEZIER_MODE);
  }

  else {
    xcoefs = generate_curve_coefs(x0, x1, x2, x3, HERMITE_MODE);
    ycoefs = generate_curve_coefs(y0, y1, y2, y3, HERMITE_MODE);
  }

  /*
  printf("a = %lf b = %lf c = %lf d = %lf\n", xcoefs->m[0][0],
         xcoefs->m[1][0], xcoefs->m[2][0], xcoefs->m[3][0]);
  */

  for (t=step; t <= 1; t+= step) {
    
    x = xcoefs->m[0][0] * t * t * t + xcoefs->m[1][0] * t * t
      + xcoefs->m[2][0] * t + xcoefs->m[3][0];

    y = ycoefs->m[0][0] * t * t * t + ycoefs->m[1][0] * t * t
      + ycoefs->m[2][0] * t + ycoefs->m[3][0];

    add_edge(points, x0, y0, 0, x, y, 0);
    x0 = x;
    y0 = y;
  }

  free_matrix(xcoefs);
  free_matrix(ycoefs);
}

/*======== void add_point() ==========
Inputs:   struct matrix * points
         int x
         int y
         int z 
Returns: 
adds point (x, y, z) to points and increment points.lastcol
if points is full, should call grow on points
====================*/
void add_point( struct matrix * points, double x, double y, double z) {
  
  if ( points->lastcol == points->cols )
    grow_matrix( points, points->lastcol + 100 );

  points->m[0][points->lastcol] = x;
  points->m[1][points->lastcol] = y;
  points->m[2][points->lastcol] = z;
  points->m[3][points->lastcol] = 1;

  points->lastcol++;
}

/*======== void add_edge() ==========
Inputs:   struct matrix * points
          int x0, int y0, int z0, int x1, int y1, int z1
Returns: 
add the line connecting (x0, y0, z0) to (x1, y1, z1) to points
should use add_point
====================*/
void add_edge( struct matrix * points, 
	       double x0, double y0, double z0, 
	       double x1, double y1, double z1) {
  add_point( points, x0, y0, z0 );
  add_point( points, x1, y1, z1 );
}

/*======== void draw_lines() ==========
Inputs:   struct matrix * points
         screen s
         color c 
Returns: 
Go through points 2 at a time and call draw_line to add that line
to the screen
====================*/
void draw_lines( struct matrix * points, screen s, color c) {

  int i;
 
  if ( points->lastcol < 2 ) {
    
    printf("Need at least 2 points to draw a line!\n");
    return;
  }

  for ( i = 0; i < points->lastcol - 1; i+=2 ) {

    draw_line( points->m[0][i], points->m[1][i], 
	       points->m[0][i+1], points->m[1][i+1], s, c);
    //FOR DEMONSTRATION PURPOSES ONLY
    //draw extra pixels so points can actually be seen    
    /*
    draw_line( points->m[0][i]+1, points->m[1][i], 
	       points->m[0][i+1]+1, points->m[1][i+1], s, c);
    draw_line( points->m[0][i], points->m[1][i]+1, 
	       points->m[0][i+1], points->m[1][i+1]+1, s, c);
    draw_line( points->m[0][i]-1, points->m[1][i], 
	       points->m[0][i+1]-1, points->m[1][i+1], s, c);
    draw_line( points->m[0][i], points->m[1][i]-1, 
	       points->m[0][i+1], points->m[1][i+1]-1, s, c);
    draw_line( points->m[0][i]+1, points->m[1][i]+1, 
	       points->m[0][i+1]+1, points->m[1][i+1]+1, s, c);
    draw_line( points->m[0][i]-1, points->m[1][i]+1, 
	       points->m[0][i+1]-1, points->m[1][i+1]+1, s, c);
    draw_line( points->m[0][i]-1, points->m[1][i]-1, 
	       points->m[0][i+1]-1, points->m[1][i+1]-1, s, c);
    draw_line( points->m[0][i]+1, points->m[1][i]-1, 
	       points->m[0][i+1]+1, points->m[1][i+1]-1, s, c);
    */
  } 	       
}


void draw_line(int x0, int y0, int x1, int y1, screen s, color c) {
 
  int x, y, d, dx, dy;

  x = x0;
  y = y0;
  
  //swap points so we're always draing left to right
  if ( x0 > x1 ) {
    x = x1;
    y = y1;
    x1 = x0;
    y1 = y0;
  }

  //need to know dx and dy for this version
  dx = (x1 - x) * 2;
  dy = (y1 - y) * 2;

  //positive slope: Octants 1, 2 (5 and 6)
  if ( dy > 0 ) {

    //slope < 1: Octant 1 (5)
    if ( dx > dy ) {
      d = dy - ( dx / 2 );
  
      while ( x <= x1 ) {
	plot(s, c, x, y);

	if ( d < 0 ) {
	  x = x + 1;
	  d = d + dy;
	}
	else {
	  x = x + 1;
	  y = y + 1;
	  d = d + dy - dx;
	}
      }
    }

    //slope > 1: Octant 2 (6)
    else {
      d = ( dy / 2 ) - dx;
      while ( y <= y1 ) {

	plot(s, c, x, y );
	if ( d > 0 ) {
	  y = y + 1;
	  d = d - dx;
	}
	else {
	  y = y + 1;
	  x = x + 1;
	  d = d + dy - dx;
	}
      }
    }
  }

  //negative slope: Octants 7, 8 (3 and 4)
  else { 

    //slope > -1: Octant 8 (4)
    if ( dx > abs(dy) ) {

      d = dy + ( dx / 2 );
  
      while ( x <= x1 ) {

	plot(s, c, x, y);

	if ( d > 0 ) {
	  x = x + 1;
	  d = d + dy;
	}
	else {
	  x = x + 1;
	  y = y - 1;
	  d = d + dy + dx;
	}
      }
    }

    //slope < -1: Octant 7 (3)
    else {

      d =  (dy / 2) + dx;

      while ( y >= y1 ) {
	
	plot(s, c, x, y );
	if ( d < 0 ) {
	  y = y - 1;
	  d = d + dx;
	}
	else {
	  y = y - 1;
	  x = x + 1;
	  d = d + dy + dx;
	}
      }
    }
  }
}

