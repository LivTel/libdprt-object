/*   
     Copyright 2006, Astrophysics Research Institute, Liverpool John Moores University.

     This file is part of libobject.

     libobject is free software; you can redistribute it and/or modify
     it under the terms of the GNU General Public License as published by
     the Free Software Foundation; either version 2 of the License, or
     (at your option) any later version.

     libobject is distributed in the hope that it will be useful,
     but WITHOUT ANY WARRANTY; without even the implied warranty of
     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
     GNU General Public License for more details.

     You should have received a copy of the GNU General Public License
     along with libobject; if not, write to the Free Software
     Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
*/
/* object.c
** Entry point for Object detection algorithm.
** $Header: /space/home/eng/cjm/cvs/libdprt-object/c/object_jmm.c,v 1.11 2008-05-29 15:35:48 eng Exp $
*/
/**
 * object.c is the main object detection source file.
 * N.B. Derived from getObjectList in dprt_libfits.c. 
 * Note that routine was wrong the following ways:
 * <ul>
 * <li>The old routine (should have) crashed if all the objects were smaller than npix.
 * <li>The old routine seemed to subtract the image_median off the pixel value when calculating 
 *     intensity in calc_object_fwhms, when it had already been subtracted in getObjectList_connect_pixels.
 * </ul>
 * @author Chris Mottram, LJMU
 * @version $Revision: 1.11 $
 */



/*
  $Log: not supported by cvs2svn $
  Revision 1.10  2008/05/06 10:57:15  eng
  Put all new debug logging behind ifdefs (log level 0).

  Revision 1.9  2008/04/30 15:18:19  eng
  Now selects the N brightest objects (default N = 11) sorting by numpix (roughly corresponding
  to area), then resorting by FWHM and taking the median.

  Briefly changed centroid position (xpos,ypos) to (xpos+1,ypos+1) in mistaken belief that
  centroiding code was off by one. Turns out it was difference in reporting position of centroid
  by C code with respect to IRAF code that was causing the confusion (C counts from zero but IRAF
  counts from 1). However, legacy of adjustment code left in (commented out!) as a reminder.

  Revision 1.8  2008/03/06 12:23:06  eng
  Checked in temporarily for laughs.

  Revision 1.7  2008/02/05 18:30:18  eng
  This version (1.6) uses a Gradient Descent (GD) routine in optimise() to curve fit. Should be faster and more accurate - testing will show. Current implementation using default settings. Also added a few lines in Object_Calculate_FWHM to write best fit parameters to new variables in w_object, for diagnostic purposes.

  Revision 1.6  2008/02/05 18:11:10  eng
  Slight tweak to write moffat curve fitting parameters into new w_object variables
  in Object_Calculate_FWHM.

  Revision 1.5  2008/01/09 14:15:32  eng
  Implementation in Object_Calculate_FWHM() of a brute force method of fitting a Moffat curve to the radial
  profile of each object. Manual testing on a few test objects shows it fits well to the data and produces
  a FWHM within 0.2 pixels (0.05") of that calculated by IRAF.

  __Radial profile__ entails converting x,y,z -> r,z where r = sqrt((x-xpos)^2 + (y-ypos)^2), i.e. based
  around the barycentre (1st moment). Using a gsl routine for this as it executes a little faster than
  standard sqrt math lib function and we have a lot of pixels to work on.

  __Moffat curve__ is a better fit to a stellar radial profile than a gaussian, and is given by the equation
  y = k * [ 1 + (x/a)^2 ]^(-b), where k = peak at x=0 and a,b are parameters to find. 'k' can be set to the
  nearest (brightest) pixel to the barycentre as it's close enough.

  This version (coded in December 2007) uses a crude scan through a range of a,b values to find the best fit,
  purely to enable testing on a large number of objects in many frames.

  Revision 1.4  2007/11/23 19:44:49  eng
  Thought some tweaks here were necessary to enable testing of Chris Simpson's idea
  of a FWHM workaround, but turns out it could all be done in object_test_jmm.c
  instead. Initially added code was then deleted again, so apart from some different
  whitespace, the code in this version therefore should be no different from that
  in the previous version.

  Revision 1.3  2007/11/14 13:38:31  eng
  Added extra debugging to print out pixel lists for objects. (CJM)

  Revision 1.2  2007/11/08 17:02:03  jmm
  Lots of changes...!


*/




/**
 * This hash define is needed before including source files give us POSIX.4/IEEE1003.1b-1993 prototypes.
 */
#define _POSIX_SOURCE 1
/**
 * This hash define is needed before including source files give us POSIX.4/IEEE1003.1b-1993 prototypes.
 */
#define _POSIX_C_SOURCE 199309L


#define MIN(a, b)  (((a) < (b)) ? (a) : (b))


#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include "object_jmm.h"
#include <float.h>

/* for new fwhm */
/*#include <gsl/gsl_sf.h> */


/* ------------------------------------------------------- */
/* hash definitions */
/* ------------------------------------------------------- */
/**
 * Value to use for seeing when something goes wrong (in pixels).
 * Should be larger than 10 pixels, so RCS thinks seeing is "bad".
 * RJS wants really large value so archive searches can differentiate between real bad seeing and failed reductions.
 */
#define MAX_N_FWHM            (11)               /* deliberately chosen to be odd */
#define MAX_N_FWHM_MID        (5)                /* median for 11 items           */

#define DEFAULT_BAD_SEEING         (999.0)       /* generic "bad seeing" flag                  */
#define DEFAULT_SEEING_NONSTELLAR  (988.0)       /* too elliptical                             */
#define DEFAULT_SEEING_TOOSMALL    (977.0)       /* less than 0.01 pixels                      */
#define DEFAULT_SEEING_BADFIT      (966.0)       /* fitted moffat curve fwhm > object diameter */

#define STELLAR_ELLIP_LIMIT          (0.2)         /* upper limit of ellipticity for an object to be */
                                                 /* classed as 'stellar' */

/* ------------------------------------------------------- */
/* structure declarations */
/* ------------------------------------------------------- */
/**
 * Structure to hold a pixel position. Lists of these points are used to get round recursion problems.
 * <ul>
 * <li><b>x</b>
 * <li><b>y</b>
 * <li><b>next_point</b>
 * </ul>
 */
struct Point_Struct
{
  int x;
  int y;
  struct Point_Struct *next_point;
};


/**
 * Structure to sort object fwhms by object "size" (numpix).
 * <ul>
 * <li><b>int numpix</b>
 * <li><b>float fwhm</b>
 * </ul>
 */
struct sizefwhm
{
  int numpix;
  float fwhm;
  int objnum;
  float xpos;
  float ypos;
};




/**
 * Structure containing logging data.
 * <ul>
 * <li><b>Log_Handler</b> Function pointer to the routine that will log messages passed to it.
 * <li><b>Log_Filter</b> Function pointer to the routine that will filter log messages passed to it.
 *                   The funtion will return TRUE if the message should be logged, and FALSE if it shouldn't.
 * <li><b>Log_Filter_Level</b> A globally maintained log filter level. 
 *                             This is set using Object_Set_Log_Filter_Level.
 * </ul>
 */
struct Log_Struct
{
  void (*Log_Handler)(int level,char *string);
  int (*Log_Filter)(int level);
  int Log_Filter_Level;
};



/* ------------------------------------------------------- */
/* internal variables */
/* ------------------------------------------------------- */
/**
 * Revision Control System identifier.
 */
/*static char rcsid[] = "$Id: object_jmm.c,v 1.11 2008-05-29 15:35:48 eng Exp $";*/
/**
 * Internal Error Number - set this to a unique value for each location an error occurs.
 */
static int Object_Error_Number = 0;
/**
 * Internal Error String - set this to a descriptive string each place an error occurs.
 * Ensure the string is not longer than OBJECT_ERROR_STRING_LENGTH long.
 * @see #Object_Error_Number
 * @see #OBJECT_ERROR_STRING_LENGTH
 */
static char Object_Error_String[OBJECT_ERROR_STRING_LENGTH] = "";
/**
 * Structure containing details of logging information.
 * @see #Log_Struct
 */
static struct Log_Struct Log_Data;
/**
 * General buffer used for string formatting during logging.
 * @see #OBJECT_ERROR_STRING_LENGTH
 */
static char Object_Buff[OBJECT_ERROR_STRING_LENGTH];

/* ------------------------------------------------------- */
/* internal function declarations */
/* ------------------------------------------------------- */
static int Object_List_Get_Connected_Pixels(int naxis1,int naxis2,float image_median,int x,int y,float thresh,
					    float *image,Object *w_object);
static void Object_Calculate_FWHM(Object *w_object,float BGmedian,int *is_stellar,float *fwhm);
static void Object_Free(Object **w_object);
static int Point_List_Remove_Head(struct Point_Struct **point_list,int *point_count);
static int Point_List_Add(struct Point_Struct **point_list,int *point_count,struct Point_Struct **last_point,
			  int x,int y);


/* ------------------------------------------------------- */
/* external functions */
/* ------------------------------------------------------- */
/*static int Sort_Float(const void *data1,const void *data2);*/
double moffat(double x, double k, double a, double b);
double delta(const double *x, const double *y, const int items, const double parameters[]);
int sign(double x);
double findMax(const double *a, const int items);
double optimize(const double *x, const double *y, int items, double params[]);
int intcmp(const void *v1, const void *v2);
int sizefwhm_cmp_by_numpix(const void *v1, const void *v2);
int sizefwhm_cmp_by_fwhm(const void *v1, const void *v2);

#define MAX_ITERS 2000
#define EARLY_STOP 4.5
#define EPS 10e-10





/*
  ---------------------------------------------------------------------
  ___   _      _           _       _     _      _       ___       _   
 / _ \ | |__  (_) ___  __ | |_    | |   (_) ___| |_    / __| ___ | |_ 
| (_) || '_ \ | |/ -_)/ _||  _|   | |__ | |(_-<|  _|  | (_ |/ -_)|  _|
 \___/ |_.__/_/ |\___|\__| \__|___|____||_|/__/ \__|___\___|\___| \__|
            |__/              |___|                |___|              
*/

/**
 * Routine to get a list of objects on the image.
 * @param image A float array containing the image data.
 *     <b>Note, this function is destructive to the contents of this array.</b>
 * @param naxis1 The length of the first axis.
 * @param naxis2 The length of the second axis.
 * @param thresh The minimum value in the array that is considered 'not background'.
 * @param npix The minimum number of pixels in something that IS an object.
 * @param first_object The address of a pointer to an object, the first in a linked list. This list is filled
 *       with allocated Object's which will need freeing. This list can be NULL, if no objects are found.
 * @param sflag The address of an integer to store a boolean flag. If set to 1, the seeing 
 *   is faked, otherwise it is real.
 * @param seeing The address of a float to return the object's seeing, in <b>pixels</b>.
 * @return Return TRUE on success, FALSE on failure.
 * @see #DEFAULT_BAD_SEEING
 * @see #Sort_Float
 * @see #Object_List_Get_Connected_Pixels
 * @see #Object_Free
 * @see #Object_Calculate_FWHM
 */
int Object_List_Get(float *image,float image_median,int naxis1,int naxis2,float thresh,int npix,
		    Object **first_object,int *sflag,float *seeing,
		    int *initial_count, int *size_count, int *stellar_count, int *fwhm_lt_dia_count)
{
  Object *w_object = NULL;
  Object *last_object = NULL;
  Object *next_object = NULL;
  HighPixel *curpix = NULL;
  float fwhm;
  float *fwhm_list = NULL;
  int y,x,done,is_stellar;
  int fwhmarray_size = 0;
  struct sizefwhm *fwhmarray = NULL;        /* array for objects whose fwhm is smaller than its diameter */
  int obj_area;                             /* number of pixels in object */
  float obj_fwhm;                           /* object fwhm in pixels */
  float obj_dia;                            /* object pseudo-diameter (pixels) */
  int mid_posn;                             /* middle position of fwhmarray, to find median */
  int lower_mid_posn,upper_mid_posn;        /* array positions either side of median, for even-sized fwhmarray */
  float median_fwhm;                        /* median fwhm obtained from fwhmarray */
  int i;                                    /* generic counter */


  /* Initialise object counters */
  *initial_count = 0;               /* initial count of all objects */
  *size_count = 0;                  /* objects bigger than size limit (currently 8 pixels) */
  *stellar_count = 0;               /* objects with ellipticity below limit */
  *fwhm_lt_dia_count = 0;           /* objects where fwhm < diameter (calculated from size) */

  /* Note therefore that: initial_count > size_count > stellar_count > fwhm_lt_dia_count */



  Object_Error_Number = 0;
#ifdef MEMORYCHECK
  if(first_object == NULL)
    {
      Object_Error_Number = 2;
      sprintf(Object_Error_String,"Object_List_Get:first_object was NULL.");
      return FALSE;
    }
  if(sflag == NULL)
    {
      Object_Error_Number = 4;
      sprintf(Object_Error_String,"Object_List_Get:sflag was NULL.");
      return FALSE;
    }
  if(seeing == NULL)
    {
      Object_Error_Number = 5;
      sprintf(Object_Error_String,"Object_List_Get:seeing was NULL.");
      return FALSE;
    }
#endif
  /* ensure top of list is set to NULL - assume is has not been allocated. */
  (*first_object) = NULL;
#if LOGGING > 0
  Object_Log(OBJECT_LOG_BIT_GENERAL,"Object_List_Get:Searching for objects.");
#endif

  for(y=0;y<naxis2;y++)
    {
      for(x=0;x<naxis1;x++)
	{
#if LOGGING > 9
	  Object_Log_Format(OBJECT_LOG_BIT_PIXEL,"Object_List_Get:searching pixel %d,%d.",x,y);
#endif
	  if(image[(y*naxis1)+x] > thresh)
	    {
	      (*initial_count)++;
	      w_object = (Object *) malloc(sizeof(Object));
#ifdef MEMORYCHECK
	      if(w_object == NULL)
		{
		  Object_Error_Number = 1;
		  sprintf(Object_Error_String,"Object_List_Get:Failed to allocate w_object.");
		  return FALSE;
		}
#endif
#if LOGGING > 10
	      Object_Log_Format(OBJECT_LOG_BIT_OBJECT,"Object_List_Get:allocated w_object (%p).",
				w_object);
#endif
	      w_object->nextobject=NULL;
	      w_object->highpixel = NULL;
	      w_object->last_hp = NULL;
	      if((*first_object)==NULL)
		{
		  (*first_object) = w_object;
		  last_object = w_object;
#if LOGGING > 10
		  Object_Log_Format(OBJECT_LOG_BIT_OBJECT,"Object_List_Get:"
				    "set first_object to (%p).",(*first_object));
#endif
		}
	      else
		{
		  last_object->nextobject = w_object;
		  last_object = w_object;
		}
	      w_object->objnum = *initial_count;
#if LOGGING > 3
	      Object_Log_Format(OBJECT_LOG_BIT_OBJECT,"Object_List_Get:"
				"found start of object at %d,%d.",x,y);
#endif
	      /* initialise stats */
	      w_object->total=0;
	      w_object->xpos=0;
	      w_object->ypos=0;
	      w_object->peak=0;
	      w_object->numpix=0;
	      if(!Object_List_Get_Connected_Pixels(naxis1,naxis2,image_median,x,y,thresh,image,
						   w_object))
		{
		  /* diddly free previous objects */
		  return FALSE;
		}
	    }/* end if threshold exceeded for image[x,y] */
	}/* end for on x */
    }/* end for on y */
#if LOGGING > 0
  Object_Log_Format(OBJECT_LOG_BIT_GENERAL,"Object_List_Get:Found %d objects.",*initial_count);
#endif
  /* have we got any objects in the frame? */
  if(*initial_count == 0)
    {
      (*seeing) = DEFAULT_BAD_SEEING;
      (*sflag) = 1; /* the seeing was fudged. */
      (*first_object) = NULL;
      Object_Error_Number = 6;
      sprintf(Object_Error_String,"Object_List_Get:No objects found.");
      Object_Warning();
      /* We used to return FALSE (error) here.
      ** But there are cases where it is OK to have no objects - e.g. Moon images.
      ** We want a fake seeing to be written to the FITS headers,
      ** So we generate a warning message and return TRUE.
      ** Note this means any program using Object_List_Get must be able to cope with
      ** a NULL object list.
      */
      return TRUE;
    }



  /* ---------------------------------- */
  /* REMOVE TOO-SMALL OBJECTS           */
  /*                                    */
  /* go through list of objects,        */
  /* get rid of any with less than npix */
  /* ---------------------------------- */


#if LOGGING > 0 
  Object_Log(OBJECT_LOG_BIT_GENERAL,"Object_List_Get:Finding useful objects.");
#endif



  /* if first object too small, delete it */
  /* ------------------------------------ */

  w_object = (*first_object);
  done = FALSE;
  while(done == FALSE){
    if(w_object->numpix >= npix)             /* if object bigger than limit */
      done = TRUE;                           /* we're done */
    else {                                   /* otherwise */


#if LOGGING > 5
      Object_Log_Format(OBJECT_LOG_BIT_OBJECT,"Object_List_Get:deleting object(1) at %.2f,%.2f(%d).",
			w_object->xpos,w_object->ypos,w_object->numpix);
#endif 


      next_object = w_object->nextobject;    /* take copy of next object pointer */
      Object_Free(&w_object);                /* delete w_object */
      w_object = next_object;                /* set w_object to next object */
      if(w_object == NULL)                   /* if we've reached the end of the list, bail out */
	done = TRUE;
    }
  }
  
  if(w_object == NULL)
    {
      (*seeing) = DEFAULT_BAD_SEEING;
      (*sflag) = 1;                       /* the seeing was fudged. */
      (*first_object) = NULL;
      Object_Error_Number = 7;
      sprintf(Object_Error_String,"Object_List_Get:All objects were too small.");
      Object_Warning();
                                           /* We used to return FALSE (error) here.
					   ** But it is OK to have all objects too small.
					   ** We want  a fake seeing to be written to the FITS headers,
					   ** So we generate a warning message and return TRUE.
					   ** Note this means any program using Object_List_Get must be able to cope with
					   ** a NULL object list.
					   */
      return TRUE;
    }


  /* set new first object */
  /* -------------------- */

#if LOGGING > 10
  Object_Log_Format(OBJECT_LOG_BIT_OBJECT,"Object_List_Get:first_object (%p) set from w_object (%p).",
		    (*first_object),w_object);
#endif


  (*first_object) = w_object;
  w_object->objnum=1;
  last_object = (*first_object);


#if LOGGING > 5
  Object_Log_Format(OBJECT_LOG_BIT_OBJECT,"Object_List_Get:object %d at %.2f,%.2f(%d) is ok(1).",
		    w_object->objnum,w_object->xpos,w_object->ypos,w_object->numpix);
#endif


  w_object = (*first_object)->nextobject;
  *size_count = 1;
  done = FALSE;


  /* go through rest of object list, deleting objects with numpix < npix */
  /* ------------------------------------------------------------------- */

  while(w_object != NULL)
    {
      next_object=w_object->nextobject;         /* take copy of next object to go to */
      if(w_object->numpix < npix)
	{


#if LOGGING > 5
	  Object_Log_Format(OBJECT_LOG_BIT_OBJECT,"Object_List_Get:deleting object(2) at %.2f,%.2f(%d).",
			    w_object->xpos,w_object->ypos,w_object->numpix);
#endif


	  Object_Free(&w_object);
	}
      else
	{
	  (*size_count)++;
	  last_object->nextobject = w_object;   /* tell last object this is its next object */
	  w_object->objnum = *size_count;       /* set objects number */
	  last_object = w_object;               /* set the last object in the list to be this object */


#if LOGGING > 5
	  Object_Log_Format(OBJECT_LOG_BIT_OBJECT,"Object_List_Get:object %d at %.2f,%.2f(%d) is ok(2).",
			    w_object->objnum,w_object->xpos,w_object->ypos,w_object->numpix);
#endif


	}
      w_object = next_object;                   /* change to next object */
    }
  last_object->nextobject=NULL;








  /* -------------------------------- */
  /* Add 1 to each object's xpos,ypos */
  /* -------------------------------- */
/*   w_object = (*first_object); */
/*   while(w_object != NULL){ */
/*     w_object->xpos++; */
/*     w_object->ypos++; */
/*     w_object = w_object->nextobject;		 */
/*   } */










  /* extra debug - list connected pixels in all objects */
  /* -------------------------------------------------- */
#if LOGGING > 5
  w_object = (*first_object);
  while(w_object != NULL)
  {
    Object_Log_Format(OBJECT_LOG_BIT_OBJECT,"Object_List_Get:Printing pixels for object %d at %.2f,%.2f(%d).",
		      w_object->objnum,w_object->xpos,w_object->ypos,w_object->numpix);
    curpix = w_object->highpixel;
    while(curpix != NULL)
      {
	Object_Log_Format(OBJECT_LOG_BIT_OBJECT,"Object_List_Get:Printing pixels:object:%d pixel %d,%d value %.2f.",
		      w_object->objnum,curpix->x,curpix->y,curpix->value);
	 curpix = curpix->next_pixel;           /* goto next pixel */
      }
      w_object = w_object->nextobject;	        /* goto next object */	
  }
#endif




  /* --------------------------- */
  /* CREATE LIST OF OBJECT FWHMS */
  /* --------------------------- */
  
#if LOGGING > 0
  Object_Log(OBJECT_LOG_BIT_GENERAL,"Object_List_Get:Finding FWHM of objects.");
#endif


  /* set first object */
  /* ---------------- */  
  w_object = (*first_object);

#if LOGGING > 10
  Object_Log_Format(OBJECT_LOG_BIT_OBJECT,"Object_List_Get:w_object (%p) set from first_object (%p).",
		    w_object,(*first_object));
#endif


  /* run through list of objects */
  /* --------------------------- */
  while(w_object != NULL)
    {

#if LOGGING > 5
      Object_Log_Format(OBJECT_LOG_BIT_FWHM,"Object_List_Get:Calculating FWHM for object at %.2f,%.2f.",
			w_object->xpos,w_object->ypos);
#endif


      /* calculate FWHM of object */
      /* ------------------------ */
      Object_Calculate_FWHM(w_object,image_median,&is_stellar,&fwhm);

#if LOGGING > 5
      Object_Log_Format(OBJECT_LOG_BIT_FWHM,"Object_List_Get:"
			"object at %.2f,%.2f has FWHM %.2f pixels and is_stellar = %d.",
			w_object->xpos,w_object->ypos,fwhm,is_stellar);
#endif


      /* if object stellar, add its FWHM to list of FWHMs */
      /* ------------------------------------------------ */
      if(is_stellar)
	{
	  if(fwhm_list == NULL)
	    fwhm_list = (float*)malloc(sizeof(float));
	  else
	    fwhm_list = (float*)realloc(fwhm_list,((*stellar_count)+1)*sizeof(float));
#ifdef MEMORYCHECK
	  if(fwhm_list == NULL)
	    {
	      Object_Error_Number = 8;
	      sprintf(Object_Error_String,"Object_List_Get:Failed to allocate FWHM list(%d).",
		      *stellar_count);
	      return FALSE;
	    }
#endif
	  fwhm_list[(*stellar_count)] = fwhm;
	  (*stellar_count)++;
	}
      w_object = w_object->nextobject;		
    }







#if LOGGING > 0
  Object_Log(OBJECT_LOG_BIT_GENERAL,"Object_List_Get:Calculating final seeing.");
#endif




  /* ----------------------------- */     /* CJM: diddly In the original code theres */
  /* CALCULATE FINAL FWHM (MEDIAN) */     /* a bit here concerned with numobjs that */
  /* ----------------------------- */     /* makes no sense at all to me */


  /* If any fwhms at all */
  /* ------------------- */
#if LOGGING > 0
  Object_Log_Format(OBJECT_LOG_BIT_GENERAL,"Number of objects potentially available: %d\n", *stellar_count);
#endif

  if(*stellar_count > 0) {

    /* Create array of numpix & fwhm structs              */
    /* NB: fwhm_count >= fwhm_lt_dia_count defined below. */
    /* Will realloc array after it's been filled.         */
    /* -------------------------------------------------- */

    fwhmarray = (struct sizefwhm *) malloc((*stellar_count) * sizeof(struct sizefwhm));



    /* populate array of structs from w_object,            */
    /* but ONLY IF fwhm < diameter of object;              */
    /* fwhm_lt_dia_count = "fwhm-less-than-diameter-count" */
    /* --------------------------------------------------- */
#if LOGGING > 0
    Object_Log_Format(OBJECT_LOG_BIT_GENERAL,"         \tfwhm\tdia\n");
#endif
    w_object = (*first_object);                           /* set w_object to first obj in list */
    while(w_object != NULL){                              /* start running through all objects */
      obj_area = w_object->numpix;
      obj_fwhm = w_object->fwhmx;	  
      obj_dia = sqrt( 1.2732 * obj_area);                 /* pseudo-diameter. 1.2732 = 4/pi    */
      if (obj_fwhm < obj_dia){                            /* add object to array only if fwhm < dia... */
	fwhmarray[*fwhm_lt_dia_count].numpix = obj_area;
	fwhmarray[*fwhm_lt_dia_count].fwhm = obj_fwhm;
	fwhmarray[*fwhm_lt_dia_count].objnum = w_object->objnum;
	fwhmarray[*fwhm_lt_dia_count].xpos = w_object->xpos;
	fwhmarray[*fwhm_lt_dia_count].ypos = w_object->ypos;
	(*fwhm_lt_dia_count)++;                              /* ... and increment counter */
      }
#if LOGGING > 0
      Object_Log_Format(OBJECT_LOG_BIT_GENERAL,"Object %d\t%f\t%f\t(%d)\n",w_object->objnum,obj_fwhm,obj_dia,*fwhm_lt_dia_count);
#endif
      w_object = w_object->nextobject;                    /* go to next object */		
    }
      

#if LOGGING > 0
    Object_Log_Format(OBJECT_LOG_BIT_GENERAL,"Number of objects with FWHM < diameter: %d\n\n", *fwhm_lt_dia_count);
#endif


    /* If some fwhms were < diameter, then */
    /* find median fwhm of top N brightest */
    /* objects in this array               */
    /* ----------------------------------- */
    if ( *fwhm_lt_dia_count > 0 ){

      /* trim fwhmarray array to exact number of objects (fwhm_lt_dia_count) */
      fwhmarray = (struct sizefwhm *) realloc(fwhmarray,(*fwhm_lt_dia_count)*sizeof(struct sizefwhm));
      
      /* set standard array size descriptor (necessary for later on) */
      fwhmarray_size = (int) (*fwhm_lt_dia_count);


#if LOGGING > 0
      Object_Log_Format(OBJECT_LOG_BIT_GENERAL,"original object list\n--------------------\n");
      for (i=0;i<fwhmarray_size;i++)
	Object_Log_Format(OBJECT_LOG_BIT_GENERAL,"[%d] (%d)\t%d\t%f\n",i,fwhmarray[i].objnum,fwhmarray[i].numpix,fwhmarray[i].fwhm);
      Object_Log_Format(OBJECT_LOG_BIT_GENERAL,"\n\n");
#endif

      /* sort array (LARGEST FIRST) by 1st struct member (numpix) */
      qsort (fwhmarray, fwhmarray_size, sizeof(struct sizefwhm), sizefwhm_cmp_by_numpix);

#if LOGGING > 0
      Object_Log_Format(OBJECT_LOG_BIT_GENERAL,"sorted by numpix\n--------------------\n");
      for (i=0;i<fwhmarray_size;i++)
	Object_Log_Format(OBJECT_LOG_BIT_GENERAL,"[%d] (%d)\t%d\t%f\n",i,fwhmarray[i].objnum,fwhmarray[i].numpix,fwhmarray[i].fwhm);
      Object_Log_Format(OBJECT_LOG_BIT_GENERAL,"\n\n");
#endif     
      
      /* now they're sorted by size, if the number of objects is greater than the maximum
	 we're going to use to find the median (i.e. the "Top N") then we need to truncate the
	 array even further, i.e. reallocate again, this time to N objects (i.e. MAX_N_FWHM).
	 NB: MAX_N_FWHM is deliberately chosen to be odd so that the median position MAX_N_FWHM_MID
	 can be stated straightaway. */
      if (fwhmarray_size > MAX_N_FWHM){
	fwhmarray = (struct sizefwhm *) realloc(fwhmarray,MAX_N_FWHM*sizeof(struct sizefwhm));
	fwhmarray_size = MAX_N_FWHM;
	
#if LOGGING > 0
	Object_Log_Format(OBJECT_LOG_BIT_GENERAL,"N > MAX_N_FWHM:\n");
#endif

	/* sort array by 2nd struct member (fwhm) SMALLEST FIRST */
	qsort (fwhmarray, fwhmarray_size, sizeof(struct sizefwhm), sizefwhm_cmp_by_fwhm);
#if LOGGING > 0
	Object_Log_Format(OBJECT_LOG_BIT_GENERAL,"truncated & sorted by fwhm\n--------------------\n");
	for (i=0;i<fwhmarray_size;i++)
	  Object_Log_Format(OBJECT_LOG_BIT_GENERAL,"toplist\t%d\t%d\t%d\t%f\t%f\t%f\n",
		 i,fwhmarray[i].objnum,
		 fwhmarray[i].numpix,fwhmarray[i].fwhm,
		 fwhmarray[i].xpos,fwhmarray[i].ypos);
	Object_Log_Format(OBJECT_LOG_BIT_GENERAL,"\n\n");
#endif
	
	/* find median */
	mid_posn = MAX_N_FWHM_MID;
	median_fwhm = fwhmarray[mid_posn].fwhm;
      }

      /* otherwise, if fwhmarray_size < MAX_N_FWHM */
      else {

#if LOGGING > 0	
	Object_Log_Format(OBJECT_LOG_BIT_GENERAL,"N < MAX_N_FWHM:\n");
#endif
	/* sort array by 2nd struct member (fwhm) SMALLEST FIRST */
	qsort (fwhmarray, fwhmarray_size, sizeof(struct sizefwhm), sizefwhm_cmp_by_fwhm);
#if LOGGING > 0	
	Object_Log_Format(OBJECT_LOG_BIT_GENERAL,"sorted by FWHM\n---------\n");
	for (i=0;i<fwhmarray_size;i++)
	  Object_Log_Format(OBJECT_LOG_BIT_GENERAL,"[%d] (%d)\t%d\t%f\t%f\t%f\n",
		 i,fwhmarray[i].objnum,
		 fwhmarray[i].numpix,fwhmarray[i].fwhm,
		 fwhmarray[i].xpos,fwhmarray[i].ypos);
	Object_Log_Format(OBJECT_LOG_BIT_GENERAL,"\n\n");
#endif
	
	/* find median */
	/* fwhmarray_size EVEN */
	if (fwhmarray_size % 2 == 0){
	  lower_mid_posn = (int) ((fwhmarray_size - 1)/2);
	  upper_mid_posn = (int) (fwhmarray_size/2);
	  median_fwhm = (fwhmarray[lower_mid_posn].fwhm + fwhmarray[upper_mid_posn].fwhm)/2.0;
	}
	
	/* fwhmarray_size ODD */
	else {
	  mid_posn = (int) (fwhmarray_size/2);
	  median_fwhm = fwhmarray[mid_posn].fwhm;
	}
      }

#if LOGGING > 0	
      if ( fwhmarray_size % 2 == 0 ) /* if EVEN */
	Object_Log_Format(OBJECT_LOG_BIT_GENERAL,"median_fwhm = [%d,%d] (%d,%d) %f\n",
	       lower_mid_posn,upper_mid_posn,
	       fwhmarray[lower_mid_posn].objnum,fwhmarray[upper_mid_posn].objnum,
	       median_fwhm);
      else /* if ODD */
	Object_Log_Format(OBJECT_LOG_BIT_GENERAL,"median_fwhm = [%d] (%d) %f\n",mid_posn,fwhmarray[mid_posn].objnum,median_fwhm);
#endif

      /* Set seeing to median_fwhm */
      /* ------------------------- */
      (*seeing) = median_fwhm;
      (*sflag) = 0;                                  /* set sflag to zero (obviously) */
      
      
      /* If the seeing is less than 0.01.... */
      if ((*seeing)<=0.01){
	(*seeing) = DEFAULT_SEEING_TOOSMALL;        /* set the seeing to DEFAULT_SEEING_TOOSMALL */
	(*sflag) = 1;                               /* and set sflag to show the seeing was fudged */
      }

    } /* end fwhm_lt_dia_count > 0 */



    /* If NO object fwhms < diameter....    */
    /* i.e. if fwhm_lt_dia_count == 0       */
    /* ------------------------------------ */
    else { 
      (*seeing) = DEFAULT_SEEING_BADFIT;             /* set the seeing to show no object's fwhm was less than its dia */
      (*sflag) = 1;                                  /* show fudged */
    }



  } /* end "if any fwhms at all" */



  /* If NO fwhms at all */
  /* ------------------ */
  else {
    (*seeing) = DEFAULT_BAD_SEEING;                  /* set the seeing to DEFAULT_BAD_SEEING (pixels) */
    (*sflag) = 1;                                    /* and set sflag to show the seeing was fudged */
  }




  /* ----------- */
  /* FREE MEMORY */
  /* ----------- */
  if(fwhm_list != NULL)
    free(fwhm_list);

  if(fwhmarray != NULL)
    free(fwhmarray);





#if LOGGING > 0
  Object_Log_Format(OBJECT_LOG_BIT_GENERAL,"Object_List_Get:number of objects > %d pixels = %d",npix,*size_count);
  Object_Log_Format(OBJECT_LOG_BIT_GENERAL,"Object_List_Get:number of objects identified as stellar = %d",*stellar_count);
  if ((*sflag)==0)
    {
      Object_Log_Format(OBJECT_LOG_BIT_GENERAL,"Object_List_Get:seeing derived from stellar sources "
			"= %.2f pixels.",(*seeing));
    }
  else
    {
      Object_Log_Format(OBJECT_LOG_BIT_GENERAL,"Object_List_Get:Unable to derive seeing, "
			"faking result = %.2f pixels.",(*seeing));
    }
#endif


  return TRUE;
}









/*
---------------------------------------------------------------------
  ___   _      _           _     _     _      _     ___               
 / _ \ | |__  (_) ___  __ | |_  | |   (_) ___| |_  | __|_ _  ___  ___ 
| (_) || '_ \ | |/ -_)/ _||  _| | |__ | |(_-<|  _| | _|| '_|/ -_)/ -_)
 \___/ |_.__/_/ |\___|\__| \__| |____||_|/__/ \__| |_| |_|  \___|\___|
            |__/                                                      
*/
/**
 * Routine to free the list allocated in Object_List_Get.
 * @param list The address of a pointer to the first element in the list.
 * @return Return TRUE on success, FALSE on failure.
 * @see #Object_Free
 */
int Object_List_Free(Object **list)
{
  Object *this_object;
  Object *next_object;

  this_object = (*list);
  if(this_object == NULL)
    return TRUE;
  while(this_object != NULL)
    {
      next_object = this_object->nextobject;
      Object_Free(&this_object);
      this_object = next_object;
    }
  return TRUE;
}




/*
---------------------------------------------------------------------
  ___   _      _           _     ___                     
 / _ \ | |__  (_) ___  __ | |_  | __| _ _  _ _  ___  _ _ 
| (_) || '_ \ | |/ -_)/ _||  _| | _| | '_|| '_|/ _ \| '_|
 \___/ |_.__/_/ |\___|\__| \__| |___||_|  |_|  \___/|_|  
            |__/                                         
*/
/**
 * The error routine that reports any errors occuring in object in a standard way.
 * @see object.html#Object_Get_Current_Time_String
 */
void Object_Error(void)
{
  char time_string[32];

  Object_Get_Current_Time_String(time_string,32);
  /* if the error number is zero an error message has not been set up
  ** This is in itself an error as we should not be calling this routine
  ** without there being an error to display */
  if(Object_Error_Number == 0)
    sprintf(Object_Error_String,"Logic Error:No Error defined");
  fprintf(stderr,"%s Object:Error(%d) : %s\n",time_string,Object_Error_Number,Object_Error_String);
}





/*
---------------------------------------------------------------------
  ___   _      _           _     ___                       _____     
 / _ \ | |__  (_) ___  __ | |_  | __| _ _  _ _  ___  _ _  |_   _|___ 
| (_) || '_ \ | |/ -_)/ _||  _| | _| | '_|| '_|/ _ \| '_|   | | / _ \
 \___/ |_.__/_/ |\___|\__| \__| |___||_|  |_|  \___/|_|     |_| \___/
            |__/                                                     
 ___  _         _             
/ __|| |_  _ _ (_) _ _   __ _ 
\__ \|  _|| '_|| || ' \ / _` |
|___/ \__||_|  |_||_||_|\__, |
                        |___/ 
*/
/**
 * The error routine that reports any errors occuring in object in a standard way. This routine places the
 * generated error string at the end of a passed in string argument.
 * @param error_string A string to put the generated error in. This string should be initialised before
 * being passed to this routine. The routine will try to concatenate it's error string onto the end
 * of any string already in existance.
 * @see object.html#Object_Get_Current_Time_String
 */
void Object_Error_To_String(char *error_string)
{
  char time_string[32];

  Object_Get_Current_Time_String(time_string,32);
  /* if the error number is zero an error message has not been set up
  ** This is in itself an error as we should not be calling this routine
  ** without there being an error to display */
  if(Object_Error_Number == 0)
    sprintf(Object_Error_String,"Logic Error:No Error defined");
  sprintf(error_string+strlen(error_string),"%s Object:Error(%d) : %s\n",time_string,
	  Object_Error_Number,Object_Error_String);
}






/*
---------------------------------------------------------------------
  ___   _      _           _      ___       _     ___                     
 / _ \ | |__  (_) ___  __ | |_   / __| ___ | |_  | __| _ _  _ _  ___  _ _ 
| (_) || '_ \ | |/ -_)/ _||  _| | (_ |/ -_)|  _| | _| | '_|| '_|/ _ \| '_|
 \___/ |_.__/_/ |\___|\__| \__|  \___|\___| \__| |___||_|  |_|  \___/|_|  
            |__/                                                          
 _  _               _              
| \| | _  _  _ __  | |__  ___  _ _ 
| .` || || || '  \ | '_ \/ -_)| '_|
|_|\_| \_,_||_|_|_||_.__/\___||_|  
                                   
*/
/**
 * Routine to return the object error number.
 * @return The object error number.
 * @see #Object_Error_Number
 */
int Object_Get_Error_Number(void)
{
  return Object_Error_Number;
}





/*
---------------------------------------------------------------------
  ___   _      _           _    __      __                 _             
 / _ \ | |__  (_) ___  __ | |_  \ \    / /__ _  _ _  _ _  (_) _ _   __ _ 
| (_) || '_ \ | |/ -_)/ _||  _|  \ \/\/ // _` || '_|| ' \ | || ' \ / _` |
 \___/ |_.__/_/ |\___|\__| \__|   \_/\_/ \__,_||_|  |_||_||_||_||_|\__, |
            |__/                                                   |___/ 
*/
/**
 * The warning routine that reports any warnings occuring in object in a standard way.
 * @see object.html#Object_Get_Current_Time_String
 */
void Object_Warning(void)
{
  char time_string[32];

  Object_Get_Current_Time_String(time_string,32);
  /* if the error number is zero an warning message has not been set up
  ** This is in itself an error as we should not be calling this routine
  ** without there being an warning to display */
  if(Object_Error_Number == 0)
    sprintf(Object_Error_String,"Logic Error:No Warning defined");
  fprintf(stderr,"%s Object:Warning(%d) : %s\n",time_string,Object_Error_Number,Object_Error_String);
}




/*
---------------------------------------------------------------------
  ___   _      _           _      ___       _   
 / _ \ | |__  (_) ___  __ | |_   / __| ___ | |_ 
| (_) || '_ \ | |/ -_)/ _||  _| | (_ |/ -_)|  _|
 \___/ |_.__/_/ |\___|\__| \__|  \___|\___| \__|
            |__/                                
  ___                           _     _____  _             
 / __|_  _  _ _  _ _  ___  _ _ | |_  |_   _|(_) _ __   ___ 
| (__| || || '_|| '_|/ -_)| ' \|  _|   | |  | || '  \ / -_)
 \___|\_,_||_|  |_|  \___||_||_|\__|   |_|  |_||_|_|_|\___|
                                                           
 ___  _         _             
/ __|| |_  _ _ (_) _ _   __ _ 
\__ \|  _|| '_|| || ' \ / _` |
|___/ \__||_|  |_||_||_|\__, |
                        |___/ 
*/
/**
 * Routine to get the current time in a string. The string is returned in the format
 * '01/01/2000 13:59:59', or the string "Unknown time" if the routine failed.
 * The time is in UTC.
 * @param time_string The string to fill with the current time.
 * @param string_length The length of the buffer passed in. It is recommended the length is at least 20 characters.
 */
void Object_Get_Current_Time_String(char *time_string,int string_length)
{
  time_t current_time;
  struct tm *utc_time = NULL;

  if(time(&current_time) > -1)
    {
      utc_time = gmtime(&current_time);
      strftime(time_string,string_length,"%d/%m/%Y %H:%M:%S",utc_time);
    }
  else
    strncpy(time_string,"Unknown time",string_length);
}



/*
---------------------------------------------------------------------
  ___   _      _           _     _               
 / _ \ | |__  (_) ___  __ | |_  | |    ___  __ _ 
| (_) || '_ \ | |/ -_)/ _||  _| | |__ / _ \/ _` |
 \___/ |_.__/_/ |\___|\__| \__| |____|\___/\__, |
            |__/                           |___/ 
 ___                        _   
| __|___  _ _  _ __   __ _ | |_ 
| _|/ _ \| '_|| '  \ / _` ||  _|
|_| \___/|_|  |_|_|_|\__,_| \__|
                                
*/
/**
 * Routine to log a message to a defined logging mechanism. This routine has an arbitary number of arguments,
 * and uses vsprintf to format them i.e. like fprintf. The Global_Buff is used to hold the created string,
 * therefore the total length of the generated string should not be longer than OBJECT_ERROR_STRING_LENGTH.
 * Object_Log is then called to handle the log message.
 * @param level An integer, used to decide whether this particular message has been selected for
 * 	logging or not.
 * @param format A string, with formatting statements the same as fprintf would use to determine the type
 * 	of the following arguments.
 * @see #Object_Log
 * @see #Object_Buff
 * @see #OBJECT_ERROR_STRING_LENGTH
 */
void Object_Log_Format(int level,char *format,...)
{
  va_list ap;

  /* Note the first two tests below were copied from Object_Log.
  ** This means for logs that do occur, they are tested twice.
  ** BUT, for logs that are to be filtered, the var args sprintf is not done.
  ** This should improve performance, if not delete Log_Handler and Log_Filter test HERE. */ 
  /* If there is no log handler, return */
  if(Log_Data.Log_Handler == NULL)
    return;
  /* If there's a log filter, check it returns TRUE for this message */
  if(Log_Data.Log_Filter != NULL)
    {
      if(Log_Data.Log_Filter(level) == FALSE)
	return;
    }
  /* format the arguments */
  va_start(ap,format);
  vsprintf(Object_Buff,format,ap);
  va_end(ap);
  /* call the log routine to log the results */
  Object_Log(level,Object_Buff);
}




/*
---------------------------------------------------------------------
  ___   _      _           _     _               
 / _ \ | |__  (_) ___  __ | |_  | |    ___  __ _ 
| (_) || '_ \ | |/ -_)/ _||  _| | |__ / _ \/ _` |
 \___/ |_.__/_/ |\___|\__| \__| |____|\___/\__, |
            |__/                           |___/ 
*/
/**
 * Routine to log a message to a defined logging mechanism. If the string or Log_Data.Log_Handler are NULL
 * the routine does not log the message. If the Log_Data.Log_Filter function pointer is non-NULL, the
 * message is passed to it to determoine whether to log the message.
 * @param level An integer, used to decide whether this particular message has been selected for
 * 	logging or not.
 * @param string The message to log.
 * @see #Log_Data
 */
void Object_Log(int level,char *string)
{
  /* If the string is NULL, don't log. */
  if(string == NULL)
    return;
  /* If there is no log handler, return */
  if(Log_Data.Log_Handler == NULL)
    return;
  /* If there's a log filter, check it returns TRUE for this message */
  if(Log_Data.Log_Filter != NULL)
    {
      if(Log_Data.Log_Filter(level) == FALSE)
	return;
    }
  /* We can log the message */
  (*Log_Data.Log_Handler)(level,string);
}



/*
---------------------------------------------------------------------
  ___   _      _           _     ___       _     _               
 / _ \ | |__  (_) ___  __ | |_  / __| ___ | |_  | |    ___  __ _ 
| (_) || '_ \ | |/ -_)/ _||  _| \__ \/ -_)|  _| | |__ / _ \/ _` |
 \___/ |_.__/_/ |\___|\__| \__| |___/\___| \__| |____|\___/\__, |
            |__/                                           |___/ 
 _  _                 _  _             ___                 _    _            
| || | __ _  _ _   __| || | ___  _ _  | __|_  _  _ _   __ | |_ (_) ___  _ _  
| __ |/ _` || ' \ / _` || |/ -_)| '_| | _|| || || ' \ / _||  _|| |/ _ \| ' \ 
|_||_|\__,_||_||_|\__,_||_|\___||_|   |_|  \_,_||_||_|\__| \__||_|\___/|_||_|
                                                                             
*/
/**
 * Routine to set the Log_Data.Log_Handler used by Object_Log.
 * @param log_fn A function pointer to a suitable handler.
 * @see #Log_Data
 * @see #Object_Log
 */
void Object_Set_Log_Handler_Function(void (*log_fn)(int level,char *string))
{
  Log_Data.Log_Handler = log_fn;
}



/*
---------------------------------------------------------------------
  ___   _      _           _     ___       _     _               
 / _ \ | |__  (_) ___  __ | |_  / __| ___ | |_  | |    ___  __ _ 
| (_) || '_ \ | |/ -_)/ _||  _| \__ \/ -_)|  _| | |__ / _ \/ _` |
 \___/ |_.__/_/ |\___|\__| \__| |___/\___| \__| |____|\___/\__, |
            |__/                                           |___/ 
 ___  _  _  _               ___                 _    _            
| __|(_)| || |_  ___  _ _  | __|_  _  _ _   __ | |_ (_) ___  _ _  
| _| | || ||  _|/ -_)| '_| | _|| || || ' \ / _||  _|| |/ _ \| ' \ 
|_|  |_||_| \__|\___||_|   |_|  \_,_||_||_|\__| \__||_|\___/|_||_|
                                                                  
*/
/**
 * Routine to set the Log_Data.Log_Filter used by Object_Log.
 * @param log_fn A function pointer to a suitable filter function.
 * @see #Log_Data
 * @see #Object_Log
 */
void Object_Set_Log_Filter_Function(int (*filter_fn)(int level))
{
  Log_Data.Log_Filter = filter_fn;
}




/*
---------------------------------------------------------------------
  ___   _      _           _     _               
 / _ \ | |__  (_) ___  __ | |_  | |    ___  __ _ 
| (_) || '_ \ | |/ -_)/ _||  _| | |__ / _ \/ _` |
 \___/ |_.__/_/ |\___|\__| \__| |____|\___/\__, |
            |__/                           |___/ 
 _  _                 _  _             ___  _       _             _   
| || | __ _  _ _   __| || | ___  _ _  / __|| |_  __| | ___  _  _ | |_ 
| __ |/ _` || ' \ / _` || |/ -_)| '_| \__ \|  _|/ _` |/ _ \| || ||  _|
|_||_|\__,_||_||_|\__,_||_|\___||_|   |___/ \__|\__,_|\___/ \_,_| \__|
                                                                      
*/

/**
 * A log handler to be used for the Log_Handler function.
 * Just prints the message to stdout, terminated by a newline.
 * @param level The log level for this message.
 * @param string The log message to be logged. 
 * @see #Log_Handler
 */
void Object_Log_Handler_Stdout(int level,char *string)
{
  char time_string[32];

  if(string == NULL)
    return;
  Object_Get_Current_Time_String(time_string,32);
  fprintf(stdout,"%s %s\n",time_string,string);
}



/*
---------------------------------------------------------------------
  ___   _      _           _     ___       _     _               
 / _ \ | |__  (_) ___  __ | |_  / __| ___ | |_  | |    ___  __ _ 
| (_) || '_ \ | |/ -_)/ _||  _| \__ \/ -_)|  _| | |__ / _ \/ _` |
 \___/ |_.__/_/ |\___|\__| \__| |___/\___| \__| |____|\___/\__, |
            |__/                                           |___/ 
 ___  _  _  _               _                    _ 
| __|(_)| || |_  ___  _ _  | |    ___ __ __ ___ | |
| _| | || ||  _|/ -_)| '_| | |__ / -_)\ V // -_)| |
|_|  |_||_| \__|\___||_|   |____|\___| \_/ \___||_|
                                                   
*/

/**
 * Routine to set the Log_Data.Log_Filter_Level.
 * @see #Log_Data
 */
void Object_Set_Log_Filter_Level(int level)
{
  Log_Data.Log_Filter_Level = level;
}





/*
---------------------------------------------------------------------
  ___   _      _           _     _                 ___  _  _  _             
 / _ \ | |__  (_) ___  __ | |_  | |    ___  __ _  | __|(_)| || |_  ___  _ _ 
| (_) || '_ \ | |/ -_)/ _||  _| | |__ / _ \/ _` | | _| | || ||  _|/ -_)| '_|
 \___/ |_.__/_/ |\___|\__| \__| |____|\___/\__, | |_|  |_||_| \__|\___||_|  
            |__/                           |___/                            
 _                    _     _    _              _        _        
| |    ___ __ __ ___ | |   /_\  | |__  ___ ___ | | _  _ | |_  ___ 
| |__ / -_)\ V // -_)| |  / _ \ | '_ \(_-</ _ \| || || ||  _|/ -_)
|____|\___| \_/ \___||_| /_/ \_\|_.__//__/\___/|_| \_,_| \__|\___|

*/
/**
 * A log message filter routine, to be used for Log_Data.Log_Filter function pointer.
 * @param level The log level of the message to be tested.
 * @return The routine returns TRUE if the level is less than or equal to the Log_Data.Log_Filter_Level,
 * 	otherwise it returns FALSE.
 * @see #Log_Data
 */
int Object_Log_Filter_Level_Absolute(int level)
{
  return (level <= Log_Data.Log_Filter_Level);
}




/*
---------------------------------------------------------------------
  ___   _      _           _     _                 ___  _  _  _             
 / _ \ | |__  (_) ___  __ | |_  | |    ___  __ _  | __|(_)| || |_  ___  _ _ 
| (_) || '_ \ | |/ -_)/ _||  _| | |__ / _ \/ _` | | _| | || ||  _|/ -_)| '_|
 \___/ |_.__/_/ |\___|\__| \__| |____|\___/\__, | |_|  |_||_| \__|\___||_|  
            |__/                           |___/                            
 _                    _   ___  _  _            _          
| |    ___ __ __ ___ | | | _ )(_)| |_ __ __ __(_) ___ ___ 
| |__ / -_)\ V // -_)| | | _ \| ||  _|\ V  V /| |(_-</ -_)
|____|\___| \_/ \___||_| |___/|_| \__| \_/\_/ |_|/__/\___|

*/
/**
 * A log message filter routine, to be used for the Log_Data.Log_Filter function pointer.
 * @param level The log level of the message to be tested.
 * @return The routine returns TRUE if the level has bits set that are also set in the 
 * 	Log_Data.Log_Filter_Level, otherwise it returns FALSE.
 * @see #Log_Data
 */
int Object_Log_Filter_Level_Bitwise(int level)
{
  return ((level & Log_Data.Log_Filter_Level) > 0);
}







/*
---------------------------------------------------------------------
  ___   _      _           _     _     _      _      ___       _   
 / _ \ | |__  (_) ___  __ | |_  | |   (_) ___| |_   / __| ___ | |_ 
| (_) || '_ \ | |/ -_)/ _||  _| | |__ | |(_-<|  _| | (_ |/ -_)|  _|
 \___/ |_.__/_/ |\___|\__| \__| |____||_|/__/ \__|  \___|\___| \__|
            |__/                                                   
  ___                            _            _   ___  _            _     
 / __| ___  _ _   _ _   ___  __ | |_  ___  __| | | _ \(_)__ __ ___ | | ___
| (__ / _ \| ' \ | ' \ / -_)/ _||  _|/ -_)/ _` | |  _/| |\ \ // -_)| |(_-<
 \___|\___/|_||_||_||_|\___|\__| \__|\___|\__,_| |_|  |_|/_\_\\___||_|/__/
                                                                          
*/

/**
 * Routine to get connected pixels starting at the specified location.
 * @param naxis1 The number of columns in the image.
 * @param naxis2 The number of rows in the image.
 * @param image_median The median pixel value in the image.
 * @param x The position in x of a pixel above the threshold.
 * @param y The position in y of a pixel above the threshold.
 * @param thresh The threshold pixel value, above which a pixel is deemed to be part of an object.
 * @param image The image data array.
 * @param w_object A pointer to a previously allocated Object, holding all data about it.
 * @return The routine returns TRUE on success and FALSE on failire.
 */

/*
  Comments added and/or tweaked by JMM 4/3/08
*/

static int Object_List_Get_Connected_Pixels(int naxis1,int naxis2,float image_median,int x,int y,float thresh,
					    float *image,Object *w_object)
{
  
  /* ------------- */
  /* SET VARIABLES */
  /* ------------- */
  int x1,y1,cx,cy;                            /* Don't know what these do */
  HighPixel *temp_hp = NULL;
  HighPixel *curpix = NULL;
  struct Point_Struct *point_list = NULL;
  struct Point_Struct *last_point = NULL;
  int point_count=0;


  /* ------------------------------- */
  /* ADD FIRST POINT TO BE PROCESSED */
  /* ------------------------------- */

#if LOGGING > 9
  Object_Log_Format(OBJECT_LOG_BIT_POINT,"Object_List_Get_Connected_Pixels:"
		    "adding point %d,%d to list.",x,y);
#endif

  if(!Point_List_Add(&point_list,&point_count,&last_point,x,y))
    return FALSE;
  

  

  /* -------------------------------- */
  /* RUN THROUGH POINTS ON POINT LIST */
  /* -------------------------------- */
  while(point_count > 0){
    

    /* start of per pixel stuff */
    /* ------------------------ */
    cx = point_list->x;
    cy = point_list->y;
    

    /* check this point hasn't been added to the object since it was added */
    /* ------------------------------------------------------------------- */
    if (image[(cy*naxis1)+cx] > thresh){   /* if pixel above threshold */      
#if LOGGING > 9
      Object_Log_Format(OBJECT_LOG_BIT_PIXEL,"Object_List_Get_Connected_Pixels:"
			"adding pixel %d,%d to object.",cx,cy);
#endif

      /* allocate new object pixel */
      temp_hp=(HighPixel*)malloc(sizeof(HighPixel));
#ifdef MEMORYCHECK
      if(temp_hp == NULL){
	Object_Error_Number = 3;
	sprintf(Object_Error_String,"Object_List_Get_Connected_Pixels:"
		"Failed to allocate temp_hp.");
	return FALSE;
      }
#endif


      /* don't know what this bit does */
      temp_hp->next_pixel=NULL;
      if(w_object->highpixel == NULL){
	w_object->highpixel = temp_hp;
	w_object->last_hp = temp_hp;
	}
      else{
	w_object->last_hp->next_pixel = temp_hp;
	w_object->last_hp = temp_hp;
      }

      
      /* set currentx, current y from point list element */
      w_object->last_hp->x=cx;
      w_object->last_hp->y=cy;
      w_object->last_hp->value=image[(cy*naxis1)+cx] - image_median;
       

      /* important for recursion - stops infinite loops */
      image[(cy*naxis1)+cx]=0.0;
	

      /* end of per-pixel stuff: do some overall w_object stats */
      curpix = w_object->last_hp;
      w_object->total=(w_object->total)+(temp_hp->value);
      w_object->xpos=(w_object->xpos)+((temp_hp->x)*(temp_hp->value)); 
      w_object->ypos=(w_object->ypos)+((temp_hp->y)*(temp_hp->value));
      if ((w_object->peak)<(temp_hp->value))
	w_object->peak=(temp_hp->value);  
      w_object->numpix ++;
      
      /* start of recursive stuff */
      for (x1 = cx-1; x1<=cx+1; x1++){
	for (y1 = cy-1; y1<=cy+1; y1++){
	  if (x1 >= naxis1 || y1 >= naxis2 || x1<0 || y1<0)  
	    continue;                                           /* set a flag here to say crap object? */
	  if (image[(y1*naxis1)+x1] > thresh){
	    /* add this point to be processed */
#if LOGGING > 9
	    Object_Log_Format(OBJECT_LOG_BIT_POINT,
			      "Object_List_Get_Connected_Pixels:"
			      "adding point %d,%d to list.",x1,y1);
#endif
	    if(!Point_List_Add(&point_list,&point_count,&last_point,x1,y1))
	      return FALSE;
	  }
	}/* end for on y1 */
      }/* end for on x1 */
    }

    /* if already added to object */
    /* -------------------------- */
    else {
#if LOGGING > 9
      Object_Log_Format(OBJECT_LOG_BIT_PIXEL,"Object_List_Get_Connected_Pixels:"
			"pixel %d,%d already added to object, ignoring.",cx,cy);
#endif
    }
    
    /* delete processed point */
    /* ---------------------- */
#if LOGGING > 9
    Object_Log_Format(OBJECT_LOG_BIT_POINT,"Object_List_Get_Connected_Pixels:"
		      "deleting point %d,%d from list.",cx,cy);
#endif
    if(!Point_List_Remove_Head(&point_list,&point_count))
      return FALSE;

  }/* end while on point list */
  



  /* -------------------------------- */
  /* CALCULATE MEAN X AND Y POSITIONS */
  /* -------------------------------- */
  w_object->xpos=(w_object->xpos)/(w_object->total);
  w_object->ypos=(w_object->ypos)/(w_object->total);








  return TRUE;
}




/*
---------------------------------------------------------------------
  ___   _      _           _     ___               
 / _ \ | |__  (_) ___  __ | |_  | __|_ _  ___  ___ 
| (_) || '_ \ | |/ -_)/ _||  _| | _|| '_|/ -_)/ -_)
 \___/ |_.__/_/ |\___|\__| \__| |_| |_|  \___|\___|
            |__/                                   
*/
/**
 * Routine to free an Object. The pointer contents themselves are freed, after
 * freeing the highpixel list inside Object. The last_hp element is NOT freed, as this should
 * point to the last element in the highpixel list that IS freed.
 * @param w_object The address of a pointer to an Object. The pointer is set to NULL. 
 */
static void Object_Free(Object **w_object)
{
  HighPixel *high_pixel,*next_pixel;

#if LOGGING > 10
  Object_Log_Format(OBJECT_LOG_BIT_OBJECT,"Object_Free:w_object (%p).",(*w_object));
#endif
  /* free highpixel list */
  high_pixel = (*w_object)->highpixel;
  while(high_pixel != NULL)
    {
      next_pixel = high_pixel->next_pixel;
      free(high_pixel);
      high_pixel = next_pixel;
    }/* end while */
  /* free w_object, and set pointer to NULL */
  free((*w_object));
  (*w_object) = NULL;
}





/*
---------------------------------------------------------------------
 ___       _       _     _     _      _     ___                            
| _ \ ___ (_) _ _ | |_  | |   (_) ___| |_  | _ \ ___  _ __   ___ __ __ ___ 
|  _// _ \| || ' \|  _| | |__ | |(_-<|  _| |   // -_)| '  \ / _ \\ V // -_)
|_|  \___/|_||_||_|\__| |____||_|/__/ \__| |_|_\\___||_|_|_|\___/ \_/ \___|
                                                                           
 _  _                _ 
| || | ___  __ _  __| |
| __ |/ -_)/ _` |/ _` |
|_||_|\___|\__,_|\__,_|
                       
*/
/**
 * Routine to remove the first point in the list
 * @param point_list The address of the pointer pointing to the first item in the list.
 * @param point_count The address of an integer holding the number of elements in the list.
 * @return The routine returns TRUE on success and FALSE on failire.
 */
static int Point_List_Remove_Head(struct Point_Struct **point_list,int *point_count)
{
  struct Point_Struct *old_head = NULL;

  /* check parameters */
#ifdef MEMORYCHECK
  if(point_list == NULL)
    {
      Object_Error_Number = 9;
      sprintf(Object_Error_String,"Point_List_Remove_Head:point_list was NULL.");
      return FALSE;
    }
  if(point_count == NULL)
    {
      Object_Error_Number = 10;
      sprintf(Object_Error_String,"Point_List_Remove_Head:point_count was NULL.");
      return FALSE;
    }
  /* check we have something to remove */
  if((*point_list) == NULL)
    {
      Object_Error_Number = 11;
      sprintf(Object_Error_String,"Point_List_Remove_Head:point_list was NULL (%d).",(*point_count));
      return FALSE;
    }
#endif
  /* keep old_head */
  old_head = (*point_list);
  /* delete off front of list */
  (*point_list) = (*point_list)->next_point;
  /* free deleted point */
  free(old_head);
  /* now less items in list */
  (*point_count)--;
  return TRUE;
}




/*
---------------------------------------------------------------------
 ___       _       _     _     _      _       _       _     _ 
| _ \ ___ (_) _ _ | |_  | |   (_) ___| |_    /_\   __| | __| |
|  _// _ \| || ' \|  _| | |__ | |(_-<|  _|  / _ \ / _` |/ _` |
|_|  \___/|_||_||_|\__| |____||_|/__/ \__| /_/ \_\\__,_|\__,_|
                                                              
*/
/**
 * This routine is used to add the specified point to the point_list. 
 * @param point_list A linked list of points, this pointer points to the head of the list.
 * @param point_count The number of points in the list.
 * @param last_point This can be NULL, if so the routine will find the last point in point_list and add
 *       the new point to that. If it is non-null, it should be set to point to the LAST element in point_list.
 *       Then the new point is added to that (speeds up adding point). last_point is then set to point to
 *       the new point (as it is last in the list).
 * @param x The x of the new point.
 * @param y The y of the new point.
 * @return The routine returns TRUE on success and FALSE on failire.
 */
static int Point_List_Add(struct Point_Struct **point_list,int *point_count,struct Point_Struct **last_point,
			  int x,int y)
{
  struct Point_Struct *new_point = NULL;
  struct Point_Struct *a_point = NULL;

  /* check parameters */
#ifdef MEMORYCHECK
  if(point_list == NULL)
    {
      Object_Error_Number = 12;
      sprintf(Object_Error_String,"Point_List_Add:point_list was NULL.");
      return FALSE;
    }
  if(point_count == NULL)
    {
      Object_Error_Number = 13;
      sprintf(Object_Error_String,"Point_List_Add:point_count was NULL.");
      return FALSE;
    }
#endif
  /* last_point can be null - see below */
  /* allocate new point */
  new_point = (struct Point_Struct *)malloc(sizeof(struct Point_Struct));
#ifdef MEMORYCHECK
  if(new_point == NULL)
    {
      Object_Error_Number = 14;
      sprintf(Object_Error_String,"Point_List_Add:Failed to allocate new_point.");
      return FALSE;
    }
#endif
  /* new point is at end of list */
  new_point->next_point = NULL;
  /* add new_point to end of list */
  if((last_point != NULL) && ((*last_point) != NULL))
    {
      /* if last_point was correctly set,
      ** but is not at the end of the list throw an error */
      if((*last_point)->next_point != NULL)
	{
	  Object_Error_Number = 15;
	  sprintf(Object_Error_String,"Point_List_Add:last_point specified that was not "
		  "last in the list.");
	  return FALSE;
	}
      /* add new point to end of list */
      (*last_point)->next_point = new_point;
    }
  else
    {
      /* if last_point not set, find end of point list */
      if((*point_list) == NULL)
	{
	  /* if point list is NULL, new point is start of point list */
	  (*point_list) = new_point;
	}
      else
	{
	  /* find end of point_list */
	  a_point = (*point_list);
	  while(a_point->next_point != NULL)
	    {
	      a_point = a_point->next_point;
	    }
	  /* add new point to end of list */
	  a_point->next_point = new_point;
	}
    }
  (*point_count)++;
  /* fill in new point's data */
  new_point->x = x;
  new_point->y = y;
  /* reset last_point to new_point */
  (*last_point) = new_point;
  return TRUE;
}




/*
-----------------------------------------------------------------------------
  ___   _      _           _       ___        _            _        _        
 / _ \ | |__  (_) ___  __ | |_    / __| __ _ | | __  _  _ | | __ _ | |_  ___ 
| (_) || '_ \ | |/ -_)/ _||  _|  | (__ / _` || |/ _|| || || |/ _` ||  _|/ -_)
 \___/ |_.__/_/ |\___|\__| \__|___\___|\__,_||_|\__| \_,_||_|\__,_| \__|\___|
            |__/              |___|                                          
     ___ __      __ _  _  __  __ 
    | __|\ \    / /| || ||  \/  |
    | _|  \ \/\/ / | __ || |\/| |
 ___|_|    \_/\_/  |_||_||_|  |_|
|___|                            
*/

/**
 * Routine to calculate the FWHM of the specified object.
 * @param w_object The object to calculate the FWHM from.
 * @param is_stellar The address of an integer to store a boolean. On exit of the routine,
 *        will be TRUE if stellar, FALSE if non-stellar.
 * @param fwhm An address to store the calculated full width half maximum, in pixels.
 */
static void Object_Calculate_FWHM(Object *w_object,float BGmedian,int *is_stellar,float *fwhm)
{

  /* ---------------- */
  /* VARIABLES CONFIG */
  /* ---------------- */

                                  /* object ellipticity */
                                  /* ------------------ */
  float x2nd=0,y2nd=0,xy2nd=0;    /* first and second order moments of the ellipse data */
  float xoff=0,yoff=0;            /* offsets from the centre of the object determind by */
                                  /*   the centroid curpix->x = spatial pixel values */
  float intensity=0;              /* intensity in each pixel corrected for background */
  float aux=0,aux2=0;             /* auxiliary variable */
  float minor=0,major=0;          /* semi-minor and semi-major axes of the object ellipse */
  float theta;                    /* angle of the ellipse, measured from Naxis1 to semi-major */
                                  /*   axis (ie anti clockwise) */
  float diff;                     /* difference between major and minor axes */
  float ellip;                    /* ellipticity = (major-minor)/major */
  HighPixel *curpix;              /* pixel pointer */
  char stellarflag[32];           /* stellar flag string for diagnostics */

                                  /* Moffat curve fitting optimisation */
                                  /* --------------------------------- */
  double *pixr, *pixz;            /* optimisation r,z arrays */
  int ipix;                       /* pixel counter */
  double dx,dy;                   /* pixel offset from 1st moment */
  double params[3];               /* Moffat curve parameters (k,a,b) */
  int m;                          /* optimisation infinite counter */  
  double k,a,b;                   /* Moffat curve parameters (k,a,b) */
  double FWHM;                    /* FWHM */




  

  /* --------------------------------------------- */
  /* CALCULATE OBJECT ELLIPTICITY (VIA 2ND MOMENT) */
  /*                                               */
  /* uses equations from SEXtractor and PISA       */
  /* manual - see SUN109.1                         */
  /* --------------------------------------------- */

  curpix = w_object->highpixel;
  while(curpix !=NULL)
    {
      intensity=(curpix->value);
      xoff=(w_object->xpos)-(curpix->x);
      yoff=(w_object->ypos)-(curpix->y);
      x2nd  += xoff*xoff*intensity;
      y2nd  += yoff*yoff*intensity;
      xy2nd += xoff*yoff*intensity;
      aux += intensity;
      curpix=curpix->next_pixel;
    }
  x2nd  /= aux;
  y2nd  /= aux;
  xy2nd /= aux;

  aux=0;
  aux=2*(x2nd+y2nd);
  aux2=2*(sqrt( ((x2nd-y2nd)*(x2nd-y2nd)) + 4*xy2nd*xy2nd ));
  major=sqrt(aux+aux2);
  minor=sqrt(aux-aux2);
  if (major != minor)
    theta=atan2((2*xy2nd),(x2nd-y2nd));
  else
    theta=0;
  
  diff=major-minor;
  ellip = (major-minor)/major;

  if (ellip <= STELLAR_ELLIP_LIMIT)
    {
      (*is_stellar) = TRUE;            /* object is STELLAR */
      w_object->is_stellar = TRUE;
      sprintf(stellarflag,"stellar");
    }
  else
    {
      (*is_stellar) = FALSE;           /* object is NONSTELLAR */
      w_object->is_stellar = FALSE;
      sprintf(stellarflag,"NON-stellar");
    }

  /* extra debug - show ellipticity of object */
  /* ---------------------------------------- */
#if LOGGING > 5
  Object_Log_Format(OBJECT_LOG_BIT_OBJECT,
		    "Object_Calculate_FWHM: (%d) a = %.2f, b = %.2f\tellip = %.2f\t%s",
		    w_object->objnum,
		    major,minor,
		    ellip,
		    stellarflag);
#endif
  






  /* -------------- */
  /* CALCULATE FWHM */
  /* -------------- */

  /* define default fwhms if meet rejection criteria */
  /* ----------------------------------------------- */


  if (w_object->is_stellar == FALSE){
    w_object->fwhmx = DEFAULT_SEEING_NONSTELLAR;
    w_object->fwhmy = DEFAULT_SEEING_NONSTELLAR;
    (*fwhm) = DEFAULT_SEEING_NONSTELLAR;
  }


  /*
    ---------------
    FWHM IF STELLAR
    ---------------
  */
  else {


    /*
      ------------------------------
      create optimisation r,z arrays
      ------------------------------
    */
    pixr = (double *)malloc(sizeof(double) * w_object->numpix);
    pixz = (double *)malloc(sizeof(double) * w_object->numpix);


    /*
      ---------------------------------------------
      create radial profile of object in r,z arrays
      ---------------------------------------------
    */
    ipix = 0;
    curpix = w_object->highpixel;
    while(curpix !=NULL){
      dx = curpix->x - w_object->xpos;
      dy = curpix->y - w_object->ypos;
      pixr[ipix] = sqrt( dx*dx + dy*dy );
      pixz[ipix] = curpix->value;

/*       fprintf(stdout,"%d:%d\t%f,%f,%f\t%f,%f\n", */
/* 	     w_object->objnum,ipix, */
/* 	     dx,dy,curpix->value, */
/* 	     pixr[ipix],pixz[ipix]); */

      ipix++;
      curpix=curpix->next_pixel;
    }


    /*
      --------------------------------
      first guess at Moffat parameters
      --------------------------------
    */
    params[0] = findMax(pixz, w_object->numpix);
    params[1] = 4.9;
    params[2] = 3.0;


    /*
      --------------------------
      optimise Moffat parameters
      --------------------------
    */
    for(m = 0; m < 1; m++) {
      optimize(pixr, pixz, w_object->numpix, params);
    }


    /*
      --------------
      calculate FWHM
      --------------
    */
    k = params[0]; /* k is z(r=0) */
    a = params[1];
    b = params[2];
    FWHM = 2.0 * a * sqrt( pow(2.0,(1.0/b)) -1.0 );

    w_object->fwhmx = FWHM;
    w_object->fwhmy = FWHM;
    (*fwhm) = FWHM;

    w_object->moffat_k = k;    /*     added     */
    w_object->moffat_a = a;    /*      for      */
    w_object->moffat_b = b;    /*  diagnostics  */
    
    if (pixr != NULL)
      free(pixr);
    if (pixz != NULL)
      free(pixz);

  } /* end of FWHM IF STELLAR */

}



/*
---------------------------------------------------------------------
 ___            _     ___  _             _   
/ __| ___  _ _ | |_  | __|| | ___  __ _ | |_ 
\__ \/ _ \| '_||  _| | _| | |/ _ \/ _` ||  _|
|___/\___/|_|   \__| |_|  |_|\___/\__,_| \__|
                                             
*/
/*
static int Sort_Float(const void *data1,const void *data2)
{
  return *((float *) data1) < * ((float *) data2);
}
*/
/* ---------------------------------------------------------------------
            __  __      _   
 _ __  ___ / _|/ _|__ _| |_ 
| '  \/ _ \  _|  _/ _` |  _|
|_|_|_\___/_| |_| \__,_|\__|
                            
y = k(1+(x/a)^2)^(-b)

*/

double moffat(double x, double k, double a, double b) {
  return k * pow(1.0 + pow(x / a, 2.0), -b);
}


/* ---------------------------------------------------------------------
    _     _ _        
 __| |___| | |_ __ _ 
/ _` / -_) |  _/ _` |
\__,_\___|_|\__\__,_|
                     
*/
double delta(const double *x, const double *y, const int items, const double parameters[]) {
  double sum = 0.0;
  int i;
  for(i = 0; i < items; i++) {
    sum += pow(y[i] - moffat(x[i], parameters[0], parameters[1], parameters[2]), 2.0);
  }
  return sqrt(sum);
}


/* ---------------------------------------------------------------------
  __ _         _ __  __          
 / _(_)_ _  __| |  \/  |__ ___ __
|  _| | ' \/ _` | |\/| / _` \ \ /
|_| |_|_||_\__,_|_|  |_\__,_/_\_\
                                 
*/
double findMax(const double *a, const int items) {
  double m = a[0];
  int i;
  for(i = 1; i < items; i++) {
    if(a[i] > m) {
      m = a[i];
    }
  }
  return m;
}


/* ---------------------------------------------------------------------
          _   _       _         
 ___ _ __| |_(_)_ __ (_)___ ___ 
/ _ \ '_ \  _| | '  \| (_-</ -_)
\___/ .__/\__|_|_|_|_|_/__/\___|
    |_|                         
*/

/*                          pixr             pixz     numpix        params */
double optimize(const double *x, const double *y, int items, double params[]) {
  double p[3];
  double p2[3];
  double momentum[3] = {0.01, 0.1, 0.1};
  double minDelta = 10e10;
  double best[3];
  int bestIter = 0;
  int i,j;

  
  p[0] = findMax(y, items);
  p[1] = 5.0;
  p[2] = 3.0;


  for(i = 0; i < MAX_ITERS; i++) {                 /* MAX_ITERS here */
    double d1, d2, m;
    double sm = 0.0;                          /* sum of slopes */
    p2[0] = p[0];
    p2[1] = p[1];
    p2[2] = p[2];
    d1 = delta(x, y, items, p);

    for(j = 0; j < 3; j++) {                  /* optimising all 3 parameters */
      double tmp = p[j];
      p2[j] = p[j] + EPS;
      d2 = delta(x, y, items, p2);
      m = (d2 - d1) / EPS;
      sm += (1 + m) * (1 + m);
                                              /* This is stolen from iRprop */
      if(momentum[j] * m > 0) {               /* momentum & slope same signs? */
	momentum[j] *= 1.25;                  /* accelerate forwards */
      } else {
	momentum[j] *= -0.8;                  /* slow down & backwards */
      }
      
      /* Threshold the momentum before applying */
      momentum[j] = sign(momentum[j]) * MIN(fabs(momentum[j]),0.2);

      p[j] = p2[j] - momentum[j];
      p2[j] = tmp;

      /* STOPPING CONDITIONS */

      /* by delta */
      if(d2 < minDelta) {
	minDelta = d2;
	best[0] = p2[0];
	best[1] = p2[1];
	best[2] = p2[2];
	bestIter = i;
      } else {
	if(i - bestIter > 20) {
	  break;
	}
      }
    }

    /* by slope */
    if(sm < EARLY_STOP) {
      break;
    }
  }
  
  for(i = 0; i < 3; i++) {
    params[i] = best[i];
  }
  

  return delta(x, y, items, params);
}



/* ---- */
/* SIGN */
/* ---- */
int sign(double x){
  if ( x == 0.0 )
    return 0;
  else
    return (int) (x/fabs(x)); 
}



/* ------ */
/* INTCMP */
/* ------ */
int intcmp(const void *v1, const void *v2)
{
return(*(int *)v1- *(int *)v2);
}



/* ---------------------- */
/* SIZEFWHM_CMP_BY_NUMPIX */
/* sorts by LARGEST first */
/* ---------------------- */
int sizefwhm_cmp_by_numpix(const void *v1, const void *v2)
{
  struct sizefwhm *sf1, *sf2;
  sf1 = (struct sizefwhm *) v1;
  sf2 = (struct sizefwhm *) v2;
  
  return(sf2->numpix - sf1->numpix);  /* numpix is an int */
}


/* -------------------- */
/* SIZEFWHM_CMP_BY_FWHM */
/* -------------------- */
int sizefwhm_cmp_by_fwhm(const void *v1, const void *v2)
{
  struct sizefwhm *sf1, *sf2;
  sf1 = (struct sizefwhm *) v1;
  sf2 = (struct sizefwhm *) v2;
  
  if (sf1->fwhm > sf2->fwhm)      /* this needed because fwhm is a float */
    return 1;
  else if (sf1->fwhm < sf2->fwhm)
    return -1;
  else
    return 0;
}



/*
** $Log: not supported by cvs2svn $
** Revision 1.10  2008/05/06 10:57:15  eng
** Put all new debug logging behind ifdefs (log level 0).
**
** Revision 1.9  2008/04/30 15:18:19  eng
** Now selects the N brightest objects (default N = 11) sorting by numpix (roughly corresponding
** to area), then resorting by FWHM and taking the median.
**
** Briefly changed centroid position (xpos,ypos) to (xpos+1,ypos+1) in mistaken belief that
** centroiding code was off by one. Turns out it was difference in reporting position of centroid
** by C code with respect to IRAF code that was causing the confusion (C counts from zero but IRAF
** counts from 1). However, legacy of adjustment code left in (commented out!) as a reminder.
**
** Revision 1.8  2008/03/06 12:23:06  eng
** Checked in temporarily for laughs.
**
** Revision 1.7  2008/02/05 18:30:18  eng
** This version (1.6) uses a Gradient Descent (GD) routine in optimise() to curve fit. Should be faster and more accurate - testing will show. Current implementation using default settings. Also added a few lines in Object_Calculate_FWHM to write best fit parameters to new variables in w_object, for diagnostic purposes.
**
** Revision 1.6  2008/02/05 18:11:10  eng
** Slight tweak to write moffat curve fitting parameters into new w_object variables
** in Object_Calculate_FWHM.
**
** Revision 1.5  2008/01/09 14:15:32  eng
** Implementation in Object_Calculate_FWHM() of a brute force method of fitting a Moffat curve to the radial
** profile of each object. Manual testing on a few test objects shows it fits well to the data and produces
** a FWHM within 0.2 pixels (0.05") of that calculated by IRAF.
**
** __Radial profile__ entails converting x,y,z -> r,z where r = sqrt((x-xpos)^2 + (y-ypos)^2), i.e. based
** around the barycentre (1st moment). Using a gsl routine for this as it executes a little faster than
** standard sqrt math lib function and we have a lot of pixels to work on.
**
** __Moffat curve__ is a better fit to a stellar radial profile than a gaussian, and is given by the equation
** y = k * [ 1 + (x/a)^2 ]^(-b), where k = peak at x=0 and a,b are parameters to find. 'k' can be set to the
** nearest (brightest) pixel to the barycentre as it's close enough.
**
** This version (coded in December 2007) uses a crude scan through a range of a,b values to find the best fit,
** purely to enable testing on a large number of objects in many frames.
**
** Revision 1.4  2007/11/23 19:44:49  eng
** Thought some tweaks here were necessary to enable testing of Chris Simpson's idea
** of a FWHM workaround, but turns out it could all be done in object_test_jmm.c
** instead. Initially added code was then deleted again, so apart from some different
** whitespace, the code in this version therefore should be no different from that
** in the previous version.
**
** Revision 1.3  2007/11/14 13:38:31  eng
** Added extra debugging to print out pixel lists for objects. (CJM)
**
** Revision 1.2  2007/11/08 17:02:03  jmm
** Lots of changes...!
**
** Revision 1.1  2007/09/18 17:24:41  jmm
** Initial revision
**
** Revision 1.13  2006/09/28 09:58:04  cjm
** Changed "All objects were too small" to be a warning - TRUE returned from Object_List_Get
** with a NULL object list, and fake seeing etc.
** This gives it the same response as "No objects found".
**
** Revision 1.12  2006/06/29 20:25:14  cjm
** Object_Calculate_FWHM now calculates FWHM even if is not stellar (AG change).
**
** Revision 1.11  2006/05/25 13:36:51  cjm
** Added is_stellar,fwhmx, and fwhmy to the Object_Struct for better seeing determination in client code.
**
** Revision 1.10  2006/05/16 18:47:54  cjm
** gnuify: Added GNU General Public License.
**
** Revision 1.9  2005/03/04 14:38:52  cjm
** Added DEFAULT_BAD_SEEING and changed bad seeing from 10 pixels to 999 pixels,
** for archive searches.
**
** Revision 1.8  2005/01/27 17:54:55  cjm
** Changed faked seeing results to 10 pixels from 2 - gives bad seeing measure even
** in binning 1.
** Added documentation and prints to say calculated seeing in pixels - conversion done
** in rjs_dprt.c:dprt_process.
**
** Revision 1.7  2004/11/22 12:47:13  cjm
** If Object_Get_List finds no objects in the frame, it now generates a warning and returns TRUE.
** The returned object list is NULL.
**
** Revision 1.6  2004/02/06 17:04:26  cjm
** Added reset of Object_Error_Number.
**
** Revision 1.5  2004/02/06 11:55:26  cjm
** Added protection for no stellar sources.
**
** Revision 1.4  2004/01/30 16:19:51  cjm
** Ensure top of list is set to NULL.
**
** Revision 1.3  2004/01/29 13:28:38  cjm
** Added Object_Get_Error_Number.
**
** Revision 1.2  2004/01/29 12:27:19  cjm
** Performance optimisation, of logging routines.
**
** Revision 1.1  2004/01/26 15:16:35  cjm
** Initial revision
**
*/
