#ifndef DRAW_H
#define DRAW_H

#include "matrix.h"

#define MAX_STEPS 100

double z_values[XRES][YRES];

void set_ambient(color c_ambient);
void set_constants(float tkar,float tkdr,float tksr,
		   float tkag,float tkdg,float tksg,
		   float tkab,float tkdb,float tksb,
		   float tri,float tgi,float tbi);

void add_light(int l0, int l1, int l2, int c0, int c1, int c2);

color calculate_Ia();
color calculate_Id(struct matrix *polygons, screen s, color c, int i);
color calculate_Is(struct matrix *polygons, screen s, color c, int i);

void calculate_surface_normal(struct matrix *polygons, int i);
void normalize_light(struct matrix *polygons, int i, int li);

color calculate_I(struct matrix *polygons, screen s, color c, int i);

void initialize_zs();

void move_points(struct matrix *polygons, screen s, color c,
		 int i, int T, int M, int B);

void move_points_z(struct matrix *polygons, screen s, color c,
		   int i, int T, int M, int B);

int order(struct matrix *polygons, int i, int xy);

void scanline_convert(struct matrix *polygons, screen s, color c, int i);
void scanline_convert_z(struct matrix *polygons, screen s, color c, int i);

void draw_line(int x0, int y0, 
	       int x1, int y1, 
	       screen s, color c);
void draw_line_z(int x0, int y0, int z0, 
		 int x1, int y1, int z1,
	       screen s, color c);
void add_point( struct matrix * points, 
		 double x, double y, double z);
void add_edge( struct matrix * points, 
	       double x0, double y0, double z0, 
	       double x1, double y1, double z1);
void add_polygons( struct matrix * points, 
		   double x0, double y0, double z0, 
		   double x1, double y1, double z1,
		   double x2, double y2, double z2);
void draw_lines( struct matrix * points, screen s, color c);
void draw_polygons( struct matrix * points, screen s, color c);
void draw_polygons_z( struct matrix * points, screen s, color c);

//advanced shapes
void add_circle( struct matrix * points, 
		 double cx, double cy, 
		 double r, double step );
void add_curve( struct matrix *points, 
		double x0, double y0,
		double x1, double y1,
		double x2, double y2,
		double x3, double y3,
		double step, int type );
void add_box( struct matrix *points,
	      double x, double y, double z,
	      double w, double h, double d);
void add_sphere( struct matrix * points, 
		 double cx, double cy, double cz, double r, 
		 int step );
void generate_sphere( struct matrix * points, 
		      double cx, double cy, double cz, double r, 
			   int step );
void add_torus( struct matrix * points, 
		double cx, double cy, double cz, double r1, double r2, 
		     int step );
void generate_torus( struct matrix * points, 
		     double cx, double cy, double cz, double r1, double r2, 
			   int step );
#endif
