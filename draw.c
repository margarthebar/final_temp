#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>

#include "ml6.h"
#include "display.h"
#include "draw.h"
#include "matrix.h"
#include "gmath.h"

float kar,kdr,ksr,kag,kdg,ksg,kab,kdb,ksb,ri,gi,bi;
color ca;
int lights[100][6];
float lights_normalized[100][3];
float L_dot_N[100];
int light_index = 0;
double surface_normal[3];

void set_ambient(color c_ambient){
  ca.red = c_ambient.red;
  ca.green = c_ambient.green;
  ca.blue = c_ambient.blue;
  //printf("Ambient light: r:%d, g:%d, b:%d\n",ca.red,ca.green,ca.blue);
}

void set_constants(float tkar,float tkdr,float tksr,float tkag,float tkdg,
		   float tksg,float tkab,float tkdb,float tksb,float tri,float tgi,float tbi){

  kar = tkar;
  kdr = tkdr;
  ksr = tksr;
  kag = tkag;
  kdg = tkdg;
  ksg = tksg;
  kab = tkab;
  kdb = tkdb;
  ksb = tksb;
  ri = tri;
  gi = tgi;
  bi = tbi;
  /*
  printf("Constants:\n");
  printf("kar:%f, kdr:%f, ksr:%f\n",kar,kdr,ksr);
  printf("kag:%f, kdg:%f, ksg:%f\n",kag,kdg,ksg);
  printf("kab:%f, kdb:%f, ksb:%f\n",kab,kdb,ksb);
  printf("ri:%f, gi:%f, bi:%f\n",ri,gi,bi);
  */
}
void add_light(int l0, int l1, int l2, int c0, int c1, int c2){
  lights[light_index][0] = l0;
  lights[light_index][1] = l1;
  lights[light_index][2] = l2;
  lights[light_index][3] = c0;
  lights[light_index][4] = c1;
  lights[light_index][5] = c2;
  light_index++;
  //printf("Light %d: l:%d,%d,%d c:%d,%d,%d\n",light_index-1,l0,l1,l2,c0,c1,c2);
}

color calculate_Ia(){
  color Ia;
  //Iambient = Ca * Ka
  //print_symtab();
  Ia.red = ca.red * kar;
  Ia.green = ca.green * kag;
  Ia.blue = ca.blue * kab;
  return Ia;
}

color calculate_Id(struct matrix *polygons, screen s, color c, int i){
  color Id, Idt;
  Id.red = 0;
  Id.green = 0;
  Id.blue = 0;
  //Idiffuse = Cp * Kd * (L^ (dot product) N^)
  int li;
  color cp;
  double xdot,ydot,zdot,cos_theta;
  for(li=0; li<light_index; li++){
    cp.red = lights[li][3];
    cp.green = lights[li][4];
    cp.blue = lights[li][5];

    xdot = lights_normalized[li][0]*surface_normal[0];
    ydot = lights_normalized[li][1]*surface_normal[1];
    zdot = lights_normalized[li][2]*surface_normal[2];
    cos_theta = xdot + ydot + zdot;
    //L_dot_N[li] = cos_theta;

    Idt.red = truncf(cp.red * kdr * cos_theta);
    Idt.green = truncf(cp.green * kdg * cos_theta);
    Idt.blue = truncf(cp.blue * kdb * cos_theta);

    if(Idt.red<0){
      Idt.red = 0;
    }else if(Idt.red>255){
      Idt.red = 255;
    }
    if(Idt.green<0){
      Idt.green = 0;
    }else if(Idt.green>255){
      Idt.green = 255;
    }
    if(Idt.blue<0){
      Idt.blue = 0;
    }else if(Idt.blue>255){
      Idt.blue = 255;
    }
    Id.red += Idt.red;
    Id.green += Idt.green;
    Id.blue += Idt.blue;
  }
  if(Id.red>255)
    Id.red=255;
  if(Id.green>255)
    Id.green=255;
  if(Id.blue>255)
    Id.blue=255;
  return Id;
}

color calculate_Is(struct matrix *polygons, screen s, color c, int i){
  color Is, Ist, cp;
  //Ispecular = Cp * Ks * cos(alpha)
  //cos(alpha) = R-> (dot product) V->
  //R-> = P-> + S->
  //R-> = 2[N^ (L^ (dot product) N^)] - L^
  int li,cor;
  double R[3];
  double V[3];
  double cos_alpha_d;
  float cos_alpha;
  V[0] = 0;
  V[1] = 1;
  V[2] = 2;
  double dot;
  for(li=0; li<light_index; li++){
    Is.red = 0;
    Is.green = 0;
    Is.blue = 0;
    cp.red = lights[li][3];
    cp.green = lights[li][4];
    cp.blue = lights[li][5];
    for(cor=0;cor<3;cor++){
      dot = lights_normalized[li][cor]*surface_normal[cor];
      R[cor] = 2*(surface_normal[cor]*dot) - lights_normalized[li][cor];
    }
    cos_alpha_d = (R[0]*V[0]) + (R[1]*V[1]) + (R[2]*V[2]);
    cos_alpha = (float)cos_alpha_d;
    //printf("cos_alpha:%f\n",cos_alpha);
    /*
    Ist.red = truncf(cp.red * ksr * cos_alpha);
    Ist.green = truncf(cp.green * ksg * cos_alpha);
    Ist.blue = truncf(cp.blue * ksb * cos_alpha);
    */
    //printf("ksr:%f, ksg:%f, ksb:%f\n",ksr,ksg,ksb);
    //printf("cp.red * ksr: %f\n",cp.red*ksr);
    //printf("cp.red * ksr * cos_alpha: %f\n",abs(cp.red*ksr*cos_alpha));
    Ist.red = cp.red * ksr * cos_alpha;
    Ist.green = cp.green * ksg * cos_alpha;
    Ist.blue = cp.blue * ksb * cos_alpha;
    /*if(ksr==0){
      Ist.red=0;
    }
    if(ksg==0){
      Ist.green=0;
    }
    if(ksb==0){
      Ist.blue=0;
    }*/

    if(Ist.red<0){
      Ist.red = 0;
    }else if(Ist.red>255){
      Ist.red = 255;
    }
    if(Ist.green<0){
      Ist.green = 0;
    }else if(Ist.green>255){
      Ist.green = 255;
    }
    if(Ist.blue<0){
      Ist.blue = 0;
    }else if(Ist.blue>255){
      Ist.blue = 255;
    }

    Is.red += Ist.red;
    Is.green += Ist.green;
    Is.blue += Ist.blue;

    //printf("Ist: r:%d, g:%d, b:%d\n",Ist.red,Ist.green,Ist.blue);
  }
  if(Is.red<0){
    Is.red = 0;//
  }else if(Is.red>255){
    Is.red=255;
  }
  if(Is.green<0){
    Is.green = 0;
  }else if(Is.green>255){
    Is.green=255;
  }
  if(Is.blue<0){
    Is.blue=0;
  }else if(Is.blue>255){
    Is.blue=255;
  }
  //printf("Is : r:%d, g:%d, b:%d\n\n",Is.red,Is.green,Is.blue);
  return Is;
}

void calculate_surface_normal(struct matrix *polygons, int i){
  double ax,ay,az,bx,by,bz,mag;
  double *normal;
  ax = polygons->m[0][i+1] - polygons->m[0][i];
  bx = polygons->m[0][i+2] - polygons->m[0][i];
  ay = polygons->m[1][i+1] - polygons->m[1][i];
  by = polygons->m[1][i+2] - polygons->m[1][i];
  az = polygons->m[2][i+1] - polygons->m[2][i];
  bz = polygons->m[2][i+2] - polygons->m[2][i];
  
  normal = calculate_normal(ax,ay,az,bx,by,bz);

  mag = sqrt((normal[0]*normal[0])+(normal[1]*normal[1])+(normal[2]*normal[2]));

  surface_normal[0] = normal[0]/mag;
  surface_normal[1] = normal[1]/mag;
  surface_normal[2] = normal[2]/mag;
}

void normalize_light(struct matrix *polygons, int i, int li){
  float Cx,Cy,Cz,Lx,Ly,Lz,L_mag;
  Cx = (polygons->m[0][i] + polygons->m[0][i+1] + polygons->m[0][i+2])/3;
  Cy = (polygons->m[1][i] + polygons->m[1][i+1] + polygons->m[1][i+2])/3;
  Cz = (polygons->m[2][i] + polygons->m[2][i+1] + polygons->m[2][i+2])/3;

  Lx = Cx - lights[li][0];
  Ly = Cy - lights[li][1];
  Lz = Cz - lights[li][2];

  L_mag = sqrt((Lx * Lx)+(Ly * Ly)+(Lz * Lz));

  lights_normalized[li][0] = Lx/L_mag;
  lights_normalized[li][1] = Ly/L_mag;
  lights_normalized[li][2] = Lz/L_mag;
}

color calculate_I(struct matrix *polygons, screen s, color c, int i){
  color Ia, Id, Is;
  int li;
  calculate_surface_normal(polygons, i);

  for(li=0; li<=light_index; li++){
    normalize_light(polygons, i, li);
  }

  Ia = calculate_Ia();
  Id = calculate_Id(polygons, s, c, i);
  Is = calculate_Is(polygons, s, c, i);
  //printf("Ia: r:%d, g:%d, b:%d\n",Ia.red,Ia.green,Ia.blue);
  //printf("Id: r:%d, g:%d, b:%d\n",Id.red,Id.green,Id.blue);
  //printf("Is: r:%d, g:%d, b:%d\n\n",Is.red,Is.green,Is.blue);

  c.red = Ia.red + Id.red + Is.red;
  c.green = Ia.green + Id.green + Is.green;
  c.blue = Ia.blue + Id.blue + Is.blue;

  return c;
}


void initialize_zs(){
  int x,y;
  for(y=0; y<YRES; y++){
    for(x=0; x<XRES; x++){
      //set each pixel to minimum value
      z_values[x][y] = DBL_MAX*-1;
      //printf("z_values[x][y]: %f\n",z_values[x][y]);
    }
  }
}

void move_points(struct matrix *polygons, screen s, color c, int i, int T, int M, int B){
  float x0, y0, x1, y1, d0, d1, Tx, Ty, Mx, My, Bx, By;

  x0 = polygons->m[0][i+B];
  y0 = polygons->m[1][i+B];
  x1 = polygons->m[0][i+B];
  y1 = polygons->m[1][i+B];

  Tx = polygons->m[0][i+T];
  Ty = polygons->m[1][i+T];
  Mx = polygons->m[0][i+M];
  My = polygons->m[1][i+M];
  Bx = polygons->m[0][i+B];
  By = polygons->m[1][i+B];

  int ans;
  float Rx, Cx, Lx;
  ans = order(polygons, i, 0);
  Rx = polygons->m[0][i + ans/100];
  Cx = polygons->m[0][i + (ans%100)/10];
  Lx = polygons->m[0][i + ans%10];

  //printf("(%f,%f),(%f,%f),(%f,%f)\n",Tx,Ty,Mx,My,Bx,By);

  d0 = (Tx - Bx)/(Ty - By);
  //////////Start moving the points///////////////
  if( (By < My) && (My < Ty) ){

    d1 = (Mx - Bx)/(My - By);
    while(y0 < My){
      /*
      int colors[5] = {0,64,128,192,255};
      c.red = colors[i%5];
      c.blue = colors[(i+1)%5];
      c.green = colors[(i-1)%5];
      */
      c = calculate_I(polygons, s, c, i);

      draw_line(x0,y0,x1,y1,s,c);

      d0 = (Tx - x0)/(Ty - y0);
      d1 = (Mx - x1)/(My - y1);
      x0+=d0;
      x1+=d1;
      y0++;
      y1++;

      if(x0 < Lx){
	x0 = Lx;
      }else if(x0 > Rx){
	x0 = Rx;
      }
      if(x1 < Lx){
	x1 = Lx;
      }else if(x1 > Rx){
	x1 = Rx;
      }


    }

    d1 = (Tx - Mx)/(Ty - My);

    int colors[5] = {0,64,128,192,255};
    while(y0 <= Ty){
      /*
      c.red = colors[i%5];
      c.blue = colors[(i+1)%5];
      c.green = colors[(i-1)%5];
      */
      c = calculate_I(polygons, s, c, i);

      draw_line(x0,y0,x1,y1,s,c);

      d0 = (Tx - x0)/(Ty - y0);
      d1 = (Tx - x1)/(Ty - y1);
      x0+=d0;
      x1+=d1;
      y0++;
      y1++;

      if(x0 < Lx){
	x0 = Lx;
      }else if(x0 > Rx){
	x0 = Rx;
      }
      if(x1 < Lx){
	x1 = Lx;
      }else if(x1 > Rx){
	x1 = Rx;
      }

    }
  }else{
    if(By == My){
      //x1 is on MT
      x1 = Mx;
      y1 = My;
      d1 = (Tx - Mx)/(Ty - My);
    }else if(Ty == My){
      //x1 is on BM
      d1 = (Mx - Bx)/(My - By);
    }else{
      printf("STRAIGHT LINE\n");
    }
    while(y0 <= Ty){
      /*
      int colors[5] = {0,64,128,192,255};
      c.red = colors[i%5];
      c.blue = colors[(i+1)%5];
      c.green = colors[(i-1)%5];
      */
      c = calculate_I(polygons, s, c, i);
      
      d0 = (Tx - x0)/(Ty - y0);
      if(By==My){
	d1 = (Tx - x1)/(Ty - y1);
      }else if(Ty==My){
	d1 = (Mx - x1)/(My - y1);
      }
      draw_line(x0,y0,x1,y1,s,c);
      x0+=d0;
      x1+=d1;
      y0++;
      y1++;

      if(x0 < Lx){
	x0 = Lx;
      }else if(x0 > Rx){
	x0 = Rx;
      }
      if(x1 < Lx){
	x1 = Lx;
      }else if(x1 > Rx){
	x1 = Rx;
      }
    }
  }
}

void move_points_z(struct matrix *polygons, screen s, color c, int i, int T, int M, int B){
  float x0, y0, z0, x1, y1, z1, d0, d1, d0z, d1z, Tx, Ty, Tz, Mx, My, Mz, Bx, By, Bz;

  //printf("here!");

  x0 = polygons->m[0][i+B];
  y0 = polygons->m[1][i+B];
  z0 = polygons->m[2][i+B];
  x1 = polygons->m[0][i+B];
  y1 = polygons->m[1][i+B];
  z1 = polygons->m[2][i+B];

  Tx = polygons->m[0][i+T];
  Ty = polygons->m[1][i+T];
  Tz = polygons->m[2][i+T];
  Mx = polygons->m[0][i+M];
  My = polygons->m[1][i+M];
  Mz = polygons->m[2][i+M];
  Bx = polygons->m[0][i+B];
  By = polygons->m[1][i+B];
  Bz = polygons->m[2][i+B];

  int ans;
  float Rx, Cx, Lx, Rz, Cz, Lz;
  ans = order(polygons, i, 0);
  Rx = polygons->m[0][i + ans/100];
  Cx = polygons->m[0][i + (ans%100)/10];
  Lx = polygons->m[0][i + ans%10];
  Rz = polygons->m[2][i + ans/100];
  Cz = polygons->m[2][i + (ans%100)/10];
  Lz = polygons->m[2][i + ans%10];

  //printf("(%f,%f),(%f,%f),(%f,%f)\n",Tx,Ty,Mx,My,Bx,By);

  d0 = (Tx - Bx)/(Ty - By);
  d0z = (Tz - Bz)/(Ty - By);
  //////////Start moving the points///////////////
  if( (By < My) && (My < Ty) ){

    d1 = (Mx - Bx)/(My - By);
    d1z = (Mz - Bz)/(My - By);
    while(y0 < My){
      /*
      int colors[5] = {0,64,128,192,255};
      c.red = colors[i%5];
      c.blue = colors[(i+1)%5];
      c.green = colors[(i-1)%5];
      */
      c = calculate_I(polygons, s, c, i);

      draw_line_z(x0,y0,z0,x1,y1,z1,s,c);
      //draw_line(x0,y0,x1,y1,s,c);

      d0 = (Tx - x0)/(Ty - y0);
      d1 = (Mx - x1)/(My - y1);
      d0z = (Tz - z0)/(Ty - y0);
      d1z = (Mz - z1)/(My - y1);
      x0+=d0;
      x1+=d1;
      z0+=d0z;
      z1+=d1z;
      y0++;
      y1++;

      if(x0 < Lx){
	x0 = Lx;
	z0 = Lz;
      }else if(x0 > Rx){
	x0 = Rx;
	z0 = Rz;
      }
      if(x1 < Lx){
	x1 = Lx;
	z1 = Lz;
      }else if(x1 > Rx){
	x1 = Rx;
	z1 = Rz;
      }


    }

    d1 = (Tx - Mx)/(Ty - My);
    d1z = (Tz - Mz)/(Ty - My);

    int colors[5] = {0,64,128,192,255};
    while(y0 <= Ty){
      /*
      c.red = colors[i%5];
      c.blue = colors[(i+1)%5];
      c.green = colors[(i-1)%5];
      */
      c = calculate_I(polygons, s, c, i);

      draw_line_z(x0,y0,z0,x1,y1,z1,s,c);

      d0 = (Tx - x0)/(Ty - y0);
      d1 = (Tx - x1)/(Ty - y1);
      d0z = (Tz - z0)/(Ty - y0);
      d1z = (Tz - z1)/(Ty - y1);
      x0+=d0;
      x1+=d1;
      z0+=d0z;
      z1+=d1z;
      y0++;
      y1++;

      if(x0 < Lx){
	x0 = Lx;
	z0 = Lz;
      }else if(x0 > Rx){
	x0 = Rx;
	z0 = Rz;
      }
      if(x1 < Lx){
	x1 = Lx;
	z1 = Lz;
      }else if(x1 > Rx){
	x1 = Rx;
	z1 = Rz;
      }

    }
  }else{
    if(By == My){
      //x1 is on MT
      x1 = Mx;
      y1 = My;
      z1 = Mz;
      d1 = (Tx - Mx)/(Ty - My);
      d1z = (Tz - Mz)/(Ty - My);
    }else if(Ty == My){
      //x1 is on BM
      d1 = (Mx - Bx)/(My - By);
      d1z = (Mz - Bz)/(My - By);
    }else{
      printf("STRAIGHT LINE\n");
    }
    while(y0 <= Ty){
      /*
      int colors[5] = {0,64,128,192,255};
      c.red = colors[i%5];
      c.blue = colors[(i+1)%5];
      c.green = colors[(i-1)%5];
      */
      c = calculate_I(polygons, s, c, i);

      d0 = (Tx - x0)/(Ty - y0);
      d0z = (Tz - z0)/(Ty - y0);
      if(By==My){
	d1 = (Tx - x1)/(Ty - y1);
	d1z = (Tz - z1)/(Ty - y1);
      }else if(Ty==My){
	d1 = (Mx - x1)/(My - y1);
	d1z = (Mz - z1)/(My - y1);
      }
      draw_line_z(x0,y0,z0,x1,y1,z1,s,c);
      x0+=d0;
      x1+=d1;
      z0+=d0z;
      z1+=d1z;
      y0++;
      y1++;

      if(x0 < Lx){
	x0 = Lx;
	z0 = Lz;
      }else if(x0 > Rx){
	x0 = Rx;
	z0 = Rz;
      }
      if(x1 < Lx){
	x1 = Lx;
	z1 = Lz;
      }else if(x1 > Rx){
	x1 = Rx;
	z1 = Rz;
      }
    }
  }
}

int order(struct matrix *polygons, int i, int xy){
 int T, M, B;

  T = 0;
  M = 1;
  B = 2;

  int temp;
  if(polygons->m[xy][i+M] < polygons->m[xy][i+B]){
    temp = M;
    M = B;
    B = temp;
  }
  if(polygons->m[xy][i+T] < polygons->m[xy][i+M]){
    temp = T;
    T = M;
    M = temp;
  }
  if(polygons->m[xy][i+M] < polygons->m[xy][i+B]){
    temp = M;
    M = B;
    B = temp;
  }
  if(polygons->m[xy][i+T]==polygons->m[xy][i+M]){
    if(xy==1 && (polygons->m[0][i+T]<polygons->m[0][i+M])){
      temp = T;
      T = M;
      M = temp;
    }
  }
  if(polygons->m[1][i+B]==polygons->m[1][i+M]){
    if(xy==1 && (polygons->m[0][i+M]<polygons->m[0][i+B])){
      temp = B;
      B = M;
      M = temp;
    }
  }

  int ans = (100*T) + (10*M) + B;
  return ans;
}

void scanline_convert(struct matrix *polygons, screen s, color c, int i){
  /////////////////SCANLINE CONVERSION//////////////////
  ////Designate points TOP, BOTTOM, and MIDDLE////
  int order_indices = order(polygons,i,1);
  int T, M, B;

  T = order_indices/100;
  M = (order_indices%100)/10;
  B = order_indices%10;
  ////////////////////////////////////////////////
  //Check to see if everything is in order
  if(!(polygons->m[1][i+T] >= polygons->m[1][i+M] && polygons->m[1][i+M] >= polygons->m[1][i+B])){
    printf("   Top: (%f,%f)\n",polygons->m[0][i+T],polygons->m[1][i+T]);
    printf("Middle: (%f,%f)\n",polygons->m[0][i+M],polygons->m[1][i+M]);
    printf("Bottom: (%f,%f)\n\n",polygons->m[0][i+B],polygons->m[1][i+B]);
  }

  ///////////Set x0, y0, x1, y1, and d0////////////
  move_points(polygons,s,c,i,T,M,B);
}

void scanline_convert_z(struct matrix *polygons, screen s, color c, int i){
   /////////////////SCANLINE CONVERSION//////////////////
  ////Designate points TOP, BOTTOM, and MIDDLE////
  int order_indices = order(polygons,i,1);
  int T, M, B;

  T = order_indices/100;
  M = (order_indices%100)/10;
  B = order_indices%10;
  ////////////////////////////////////////////////
  //Check to see if everything is in order
  if(!(polygons->m[1][i+T] >= polygons->m[1][i+M] && polygons->m[1][i+M] >= polygons->m[1][i+B])){
    printf("   Top: (%f,%f)\n",polygons->m[0][i+T],polygons->m[1][i+T]);
    printf("Middle: (%f,%f)\n",polygons->m[0][i+M],polygons->m[1][i+M]);
    printf("Bottom: (%f,%f)\n\n",polygons->m[0][i+B],polygons->m[1][i+B]);
  }

  ///////////Set x0, y0, x1, y1, and d0////////////
  move_points_z(polygons,s,c,i,T,M,B);
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
  for( i=0; i < polygons->lastcol-2; i+=3 ) {

    if ( calculate_dot( polygons, i ) < 0 ) {
      /*
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
      */
      //c = calculate_I(polygons, s, c, i);
      scanline_convert(polygons, s, c, i);
      //printf("Finished i:%d, total:%d\n",i,polygons->lastcol-2);
    }
  }
  //printf("finished drawing\n");
}

void draw_polygons_z( struct matrix *polygons, screen s, color c ) {
  int i;
  for( i=0; i < polygons->lastcol-2; i+=3 ) {

    if ( calculate_dot( polygons, i ) < 0 ) {
      /*
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
      */
      //c = calculate_I(polygons, s, c, i);
      scanline_convert_z(polygons, s, c, i);
      //printf("Finished i:%d, total:%d\n",i,polygons->lastcol-2);
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

void draw_line_z(int x0, int y0, int z0, int x1, int y1, int z1, screen s, color c) {
 
  double zees[XRES][YRES];
  int z_x,z_y;
  int x, y, z, d, dx, dy, dz, dd;

  x = x0;
  y = y0;
  z = z0;
  
  //swap points so we're always draing left to right
  if ( x0 > x1 ) {
    x = x1;
    y = y1;
    z = z1;
    x1 = x0;
    y1 = y0;
    z1 = z0;
  }

  //need to know dx and dy for this version
  dx = (x1 - x) * 2;
  dy = (y1 - y) * 2;
  dz = (z1 - z) * 2;

  //positive slope: Octants 1, 2 (5 and 6)
  if ( dy > 0 ) {

    //slope < 1: Octant 1 (5)
    if ( dx > dy ) {
      d = dy - ( dx / 2 );
      if( dz > 0){
	dd = dy - ( dz / 2 );
      }else{
	dd = ( dy / 2 ) - dz;
      }
  
      while ( x <= x1 ) {

	for(z_x=0;z_x<500;z_x++){
	  for(z_y=0;z_y<500;z_y++){
	    zees[z_x][z_y] = z_values[z_x][z_y];
	  }
	}
	z_values[x][YRES-1-y] = plot_z(s, c, x, y, z, zees);

	if ( d < 0 ) {
	  x = x + 1;
	  z = z + 1;
	  d = d + dy;
	}
	else {
	  x = x + 1;
	  y = y + 1;
	  z = z + 1;
	  d = d + dy - dx;
	}
	if( dd < 0 ){
	  dd = dd + dy;
	}else{
	  dd = dd + dy - dz;
	}
      }
    }

    //slope > 1: Octant 2 (6)
    else {
      d = ( dy / 2 ) - dx;

      if( dz > 0){
	dd = dy - ( dz / 2 );
      }else{
	dd = ( dz / 2 ) - dz;
      }
      while ( y <= y1 ) {
	for(z_x=0;z_x<500;z_x++){
	  for(z_y=0;z_y<500;z_y++){
	    zees[z_x][z_y] = z_values[z_x][z_y];
	  }
	}
	z_values[x][YRES-1-y] = plot_z(s, c, x, y, z, zees);
	if ( d > 0 ) {
	  y = y + 1;
	  d = d - dx;
	}
	else {
	  y = y + 1;
	  x = x + 1;
	  z = z + 1;
	  d = d + dy - dx;
	}
	if ( dd > 0 ){
	  dd = dd - dz;
	}else{
	  dd = dd + dy - dz;
	}
      }
    }
  }

  //negative slope: Octants 7, 8 (3 and 4)
  else { 

    //slope > -1: Octant 8 (4)
    if ( dx > abs(dy) ) {

      d = dy + ( dx / 2 );
      
      if( dz > abs(dy) ){
	dd = dy + ( dz / 2 );
      }else{
	dd =  (dy / 2) + dz;
      }
  
      while ( x <= x1 ) {
	
	for(z_x=0;z_x<500;z_x++){
	  for(z_y=0;z_y<500;z_y++){
	    zees[z_x][z_y] = z_values[z_x][z_y];
	  }
	}
	z_values[x][YRES-1-y] = plot_z(s, c, x, y, z, zees);

	if ( d > 0 ) {
	  x = x + 1;
	  z = z + 1;
	  d = d + dy;
	}
	else {
	  x = x + 1;
	  y = y - 1;
	  z = z + 1;
	  d = d + dy + dx;
	}
	if ( dd > 0 ){
	  dd = dd + dy;
	}else{
	  dd = dd + dy + dz;
	}
      }
    }

    //slope < -1: Octant 7 (3)
    else {

      d =  (dy / 2) + dx;
      
      if( dz > abs(dy) ){
	dd = dy + ( dz / 2 );
      }else{
	dd =  (dy / 2) + dz;
      }
      
      while ( y >= y1 ) {

	for(z_x=0;z_x<500;z_x++){
	  for(z_y=0;z_y<500;z_y++){
	    zees[z_x][z_y] = z_values[z_x][z_y];
	  }
	}
	z_values[x][YRES-1-y] = plot_z(s, c, x, y, z, zees);
	if ( d < 0 ) {
	  y = y - 1;
	  d = d + dx;
	}
	else {
	  y = y - 1;
	  x = x + 1;
	  z = z + 1;
	  d = d + dy + dx;
	}
	if ( dd < 0 ){
	  dd = dd + dz;
	}else{
	  dd = dd + dy + dz;
	}
      }
    }
  }
}
