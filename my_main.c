/*========== my_main.c ==========

  This is the only file you need to modify in order
  to get a working mdl project (for now).

  my_main.c will serve as the interpreter for mdl.
  When an mdl script goes through a lexer and parser, 
  the resulting operations will be in the array op[].

  Your job is to go through each entry in op and perform
  the required action from the list below:

  frames: set num_frames (in misc_headers.h) for animation

  basename: set name (in misc_headers.h) for animation

  vary: manipluate knob values between two given frames
        over a specified interval

  set: set a knob to a given value
  
  setknobs: set all knobs to a given value

  push: push a new origin matrix onto the origin stack
  
  pop: remove the top matrix on the origin stack

  move/scale/rotate: create a transformation matrix 
                     based on the provided values, then 
		     multiply the current top of the
		     origins stack by it.

  box/sphere/torus: create a solid object based on the
                    provided values. Store that in a 
		    temporary matrix, multiply it by the
		    current top of the origins stack, then
		    call draw_polygons.

  line: create a line based on the provided values. Store 
        that in a temporary matrix, multiply it by the
	current top of the origins stack, then call draw_lines.

  save: call save_extension with the provided filename

  display: view the image live
  
  jdyrlandweaver
  =========================*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <float.h>
#include "parser.h"
#include "symtab.h"
#include "y.tab.h"

#include "misc_headers.h"
#include "matrix.h"
#include "ml6.h"
#include "display.h"
#include "draw.h"
#include "stack.h"

/*======== void first_pass()) ==========
  Inputs:   
  Returns: 

  Checks the op array for any animation commands
  (frames, basename, vary)
  
  Should set num_frames and basename if the frames 
  or basename commands are present

  If vary is found, but frames is not, the entire
  program should exit.

  If frames is found, but basename is not, set name
  to some default value, and print out a message
  with the name being used.

  jdyrlandweaver
  ====================*/
void first_pass() {

  int i;
  short frame_check, vary_check, name_check;

  frame_check = vary_check = name_check = 0;

  for( i=0; i<lastop; i++ ) {

    switch ( op[i].opcode ) {
      
    case FRAMES:
      num_frames = op[i].op.frames.num_frames;
      frame_check = 1;
      break;
      
    case BASENAME:
      strncpy( name, op[i].op.basename.p->name, sizeof( name ) );
      name_check = 1;
      break;

    case VARY:
      vary_check = 1;
      break;
    }
  }

  if ( vary_check && !frame_check ) {
    printf( "Vary command found without setting number of frames!\n");
    exit( 0 );
  }
  
  else if ( frame_check && !name_check ){
    printf("Animation code found but basename was not set, using \"frame\" as basename\n" );
    strncpy( name, "frame", sizeof( name ) );
  }
}

/*======== struct vary_node ** second_pass()) ==========
  Inputs:   
  Returns: An array of vary_node linked lists

  In order to set the knobs for animation, we need to keep
  a seaprate value for each knob for each frame. We can do
  this by using an array of linked lists. Each array index
  will correspond to a frame (eg. knobs[0] would be the first
  frame, knobs[2] would be the 3rd frame and so on).

  Each index should contain a linked list of vary_nodes, each
  node contains a knob name, a value, and a pointer to the
  next node.

  Go through the opcode array, and when you find vary, go 
  from knobs[0] to knobs[frames-1] and add (or modify) the
  vary_node corresponding to the given knob with the
  appropirate value. 

  05/17/12 09:29:31
  jdyrlandweaver
  ====================*/
struct vary_node ** second_pass() {
  
  int i, k;
  double start_frame, end_frame, start_value, end_value;
  int count;
  char knob_names[50][128];
  short found;
  struct vary_node ** knobs;
  struct vary_node * new_node;
  struct vary_node * curr;

  knobs = (struct vary_node **)calloc( num_frames, 
				       sizeof( struct vary_node * ) );

  for ( i=0; i<lastop; i++ ) {
    
    if ( op[i].opcode == VARY ) {
      
      start_frame = op[i].op.vary.start_frame;
      end_frame = op[i].op.vary.end_frame;
      start_value = op[i].op.vary.start_val;
      end_value = op[i].op.vary.end_val;
      
      if ( start_frame < 0 || 
	   end_frame >= num_frames || 
	   end_frame < start_frame ) {

	printf( "Invalid vary command for knob: %s\n", 
		op[i].op.vary.p->name );
	exit( 0 );
      }
	
      //set knob values for each frame
      for ( k=0; k < num_frames; k++ ){
	
	found = 0;
	new_node = knobs[k];
	while ( new_node ) {
	  
	  if ( strcmp( new_node->name, op[i].op.vary.p->name) == 0 ) {
	    found = 1;
	    break;
	  }
	  
	  new_node = new_node->next;
	}
	
	if ( !found ) {
	  new_node = (struct vary_node *)calloc(1, sizeof( struct vary_node));
	  strncpy( new_node->name, op[i].op.vary.p->name, 
		   sizeof( new_node-> name ) );
	  new_node->next = knobs[k];
	  knobs[k] = new_node;
	}

	//if frame is before the knob application,
	//set knob to the start value
	if  ( !found && k < start_frame )  {
	  new_node->value = start_value;
	}
	//if frame is after the knob, set knob to end value
	else if ( !found && k > end_frame ) 
	  new_node->value = end_value;

	else if ( k >= start_frame && k <= end_frame ) {
	  if ( start_value < end_value ) {
	    new_node->value = ( k - start_frame ) *
	      (double)(end_value - start_value) / 
	      (double) (end_frame - start_frame) + start_value;
	  }
	  else 
	    new_node->value = start_value - 
	      ( start_frame - k ) *
	      (double)(end_value - start_value) / 
	      (double) (end_frame - start_frame);
	}
      }
    } //end if vary found
  } //end for loop
  return knobs;
}


/*======== void print_knobs() ==========
Inputs:   
Returns: 

Goes through symtab and display all the knobs and their
currnt values

jdyrlandweaver
====================*/
void print_knobs() {
  
  int i;

  printf( "ID\tNAME\t\tTYPE\t\tVALUE\n" );
  for ( i=0; i < lastsym; i++ ) {

    if ( symtab[i].type == SYM_VALUE ) {
      printf( "%d\t%s\t\t", i, symtab[i].name );

      printf( "SYM_VALUE\t");
      printf( "%6.2f\n", symtab[i].s.value);
    }
  }
}

/*======== void process_knobs() ==========
Inputs:   
Returns: 

Displays the current knob values and provides
an interface for the user to set them

jdyrlandweaver
====================*/
void process_knobs() {
  
  int i;
  double v;

  if ( lastsym == 0 )
    return;

  do {
    printf( "Knob List:\n" );
    print_knobs();
    printf( "Enter knob ID to set (-1 to stop): ");
    scanf( "%d", &i );
    
    if ( i >= lastsym || i < -1 )
      printf( "Invalid entry, please try again.\n" );
    
    else if ( i != -1 ) {
      
      printf( "Enter new value for %s: ", symtab[i].name );
      scanf( "%lf", &v );
      symtab[i].s.value = v;
    }
    printf("\n");
    
  } while ( i != -1 );
}


/*======== void my_main() ==========
  Inputs:   int polygons  
  Returns: 

  This is the main engine of the interpreter, it should
  handle most of the commadns in mdl.

  If frames is not present in the source (and therefore 
  num_frames is 1, then process_knobs should be called.

  If frames is present, the enitre op array must be
  applied frames time. At the end of each frame iteration
  save the current screen to a file named the
  provided basename plus a numeric string such that the
  files will be listed in order, then clear the screen and
  reset any other data structures that need it.

  Important note: you cannot just name your files in 
  regular sequence, like pic0, pic1, pic2, pic3... if that
  is done, then pic1, pic10, pic11... will come before pic2
  and so on. In order to keep things clear, add leading 0s
  to the numeric portion of the name. If you use sprintf, 
  you can use "%0xd" for this purpose. It will add at most
  x 0s in front of a number, if needed, so if used correctly,
  and x = 4, you would get numbers like 0001, 0002, 0011,
  0487

  05/17/12 09:41:35
  jdyrlandweaver
  ====================*/
void my_main( int polygons ) {

  int i, f, j;
  double step;
  double xval, yval, zval, knob_value;
  struct matrix *transform;
  struct matrix *tmp;
  struct stack *s;
  screen t;
  color g;

  color ca;
  float kar,kdr,ksr,kag,kdg,ksg,kab,kdb,ksb,r_intensity,g_intensity,b_intensity;
  int index;

  struct vary_node **knobs;
  struct vary_node *vn;
  char frame_name[128];

  num_frames = 1;
  step = 5;
 
  g.red = 0;
  g.green = 255;
  g.blue = 255;


  first_pass();

  if (num_frames == 1)
    process_knobs();
  else
    knobs = second_pass();
  
  for ( f=0; f < num_frames; f++ ) {

    s = new_stack();
    tmp = new_matrix(4, 1000);
    clear_screen( t );
    initialize_zs();

    //if there are multiple frames, set the knobs
    if ( num_frames > 1 ) {
      
      vn = knobs[f];
      while ( vn ) {
	printf("knob: %s value:%lf\n", vn->name, vn->value);
	set_value( lookup_symbol( vn->name ), vn->value );
	vn = vn-> next;
      }
    }

    /*
    //initialize z values
    //z_values = new_matrix(XRES,YRES);
    int x,y;
    for(y=0; y<YRES; y++){
      for(x=0; x<XRES; x++){
	//set each pixel to minimum value
	z_values->m[x][y] = DBL_MIN;
      }
    }
    */
    
    for (i=0;i<lastop;i++) {
  
      switch (op[i].opcode) {

      case SET:
	set_value( lookup_symbol( op[i].op.set.p->name ), 
		   op[i].op.set.p->s.value );
	break;
	
      case SETKNOBS:
	for ( j=0; j < lastsym; j++ ) 
	  if ( symtab[j].type == SYM_VALUE )
	    symtab[j].s.value = op[i].op.setknobs.value;
	break;

      case SPHERE:
	add_sphere( tmp,op[i].op.sphere.d[0], //cx
		    op[i].op.sphere.d[1],  //cy
		    op[i].op.sphere.d[2],  //cz
		    op[i].op.sphere.r,
		    step);
	//apply the current top origin
	matrix_mult( s->data[ s->top ], tmp );
	//draw_polygons( tmp, t, g );
	draw_polygons_z( tmp, t, g );
	tmp->lastcol = 0;
	break;

      case TORUS:
	add_torus( tmp, op[i].op.torus.d[0], //cx
		   op[i].op.torus.d[1],     //cy
		   op[i].op.torus.d[2],    //cz
		   op[i].op.torus.r0,
		   op[i].op.torus.r1,
		   step);
	matrix_mult( s->data[ s->top ], tmp );
	//draw_polygons( tmp, t, g );
	draw_polygons_z( tmp, t, g );
	tmp->lastcol = 0;
	break;

      case BOX:
	add_box( tmp, op[i].op.box.d0[0],
		 op[i].op.box.d0[1],
		 op[i].op.box.d0[2],
		 op[i].op.box.d1[0],
		 op[i].op.box.d1[1],
		 op[i].op.box.d1[2]);
	matrix_mult( s->data[ s->top ], tmp );
	//draw_polygons( tmp, t, g );
	draw_polygons_z( tmp, t, g );
	tmp->lastcol = 0;
	break;

      case LINE:
	add_edge( tmp, op[i].op.line.p0[0],
		  op[i].op.line.p0[1],
		  op[i].op.line.p0[1],
		  op[i].op.line.p1[0],
		  op[i].op.line.p1[1],
		  op[i].op.line.p1[1]);
	draw_lines( tmp, t, g );
	tmp->lastcol = 0;
	break;

      case MOVE:
	//get the factors
	xval = op[i].op.move.d[0];
	yval =  op[i].op.move.d[1];
	zval = op[i].op.move.d[2];
      
	//get knob if it exists
	if ( op[i].op.move.p != NULL ) {
	  knob_value = lookup_symbol( op[i].op.move.p->name )->s.value;
	  xval = xval * knob_value;
	  yval = yval * knob_value;
	  zval = zval * knob_value;
	}

	transform = make_translate( xval, yval, zval );
	//multiply by the existing origin
	matrix_mult( s->data[ s->top ], transform );
	//put the new matrix on the top
	copy_matrix( transform, s->data[ s->top ] );
	free_matrix( transform );
	break;

      case SCALE:
	xval = op[i].op.scale.d[0];
	yval = op[i].op.scale.d[1];
	zval = op[i].op.scale.d[2];
      
	//get knob if it exists
	if ( op[i].op.scale.p != NULL ) {
	  knob_value = lookup_symbol( op[i].op.scale.p->name )->s.value;
	  xval = xval * knob_value;
	  yval = yval * knob_value;
	  zval = zval * knob_value;
	}

	transform = make_scale( xval, yval, zval );
	matrix_mult( s->data[ s->top ], transform );
	//put the new matrix on the top
	copy_matrix( transform, s->data[ s->top ] );
	free_matrix( transform );
	break;

      case ROTATE:
	xval = op[i].op.rotate.degrees * ( M_PI / 180 );

	//get knob if it exists
	if ( op[i].op.rotate.p != NULL ) {
	  knob_value = lookup_symbol( op[i].op.rotate.p->name )->s.value;
	  xval = xval * knob_value;
	}

	//get the axis
	if ( op[i].op.rotate.axis == 0 ) 
	  transform = make_rotX( xval );
	else if ( op[i].op.rotate.axis == 1 ) 
	  transform = make_rotY( xval );
	else if ( op[i].op.rotate.axis == 2 ) 
	  transform = make_rotZ( xval );

	matrix_mult( s->data[ s->top ], transform );
	//put the new matrix on the top
	copy_matrix( transform, s->data[ s->top ] );
	free_matrix( transform );
	break;

      case PUSH:
	push( s );
	break;
      case POP:
	pop( s );
	break;
      case SAVE:
	save_extension( t, op[i].op.save.p->name );
	break;
      case DISPLAY:
	display( t );
	break;
      case LIGHT:
	//light r g b x y z
	//creates a "light" datastructure with rgb values r,g,b at location x,y,z.
	//This is inserted into the symbol table
        for(j=0; j<lastsym; j++){
	  if ( symtab[j].type == SYM_LIGHT ){
	    index = j;
	  }
	}
        add_light(symtab[index].s.l->l[0],symtab[index].s.l->l[1],symtab[index].s.l->l[2],
		  symtab[index].s.l->c[0],symtab[index].s.l->c[1],symtab[index].s.l->c[2]);
	break;
      case AMBIENT:
	//ambient r g b
	//specifies how much ambient light is in the scene
	//set_value( lookup_symbol( op[i].op.ambient.c[0] ), op[i].op.ambient.c[0]);
	ca.red = op[i].op.ambient.c[0];
	ca.green = op[i].op.ambient.c[1];
	ca.blue = op[i].op.ambient.c[2];
	set_ambient(ca);
	break;
      case CONSTANTS:
	for(j=0; j<lastsym; j++){
	  if ( symtab[j].type == SYM_CONSTANTS ){
	    index = j;
	  }
	}

	kar = symtab[index].s.c->r[0];
	kdr = symtab[index].s.c->r[1];
	ksr = symtab[index].s.c->r[2];

	kag = symtab[index].s.c->g[0];
	kdg = symtab[index].s.c->g[1];
	ksg = symtab[index].s.c->g[2];

	kab = symtab[index].s.c->b[0];
	kdb = symtab[index].s.c->b[1];
	ksb = symtab[index].s.c->b[2];

        r_intensity = symtab[index].s.c->red;
	g_intensity = symtab[index].s.c->green;
	b_intensity = symtab[index].s.c->blue;
	//printf("\nkar:%f\n\n",kar);

	set_constants(kar,kdr,ksr,kag,kdg,ksg,kab,kdb,ksb,r_intensity,g_intensity,b_intensity);
	//printf("\n again kar:%f\n\n",kar);
	break;
      }
    }
  
    free_stack( s );
    free_matrix( tmp );
    //free_matrix( transform );
 
    //save the image with the correct filename
    if ( num_frames > 1 ) {
      printf("Drawing frome: %d\n", f );
      sprintf( frame_name, "anim/%s%03d.png", name, f );
      save_extension( t, frame_name );
    }
  } //end frame loop

}
