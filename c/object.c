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
** $Header: /space/home/eng/cjm/cvs/libdprt-object/c/object.c,v 1.17 2023-11-10 16:47:04 cjm Exp $
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
 * @version $Revision: 1.17 $
 */




/*
  $Log: not supported by cvs2svn $
  Revision 1.16  2014/07/30 21:38:06  eng
  *** empty log message ***

  Revision 1.15  2014/07/30 18:26:21  eng
  Added calculation of ellip_theta, the long axis orientaion of ellipticity in the
  frame. The value is not yet passed back to the calling application. That comes
  in the next update.

  Revision 1.14  2014/07/30 17:42:09  eng
  Continuation of changes from v1.13. Finished adding teh function that allows you to set
  the saturation level and created a DEFAULT #def to set the satuation to 63000 ADU if not
  set explicitly by the calling application.

  When determining the list of candidate objects which will be used to derive the
  median seeing, added a check to reject any sources with peak > saturation limit. Only
  unsaturated stars to be used to measure seeing.

  Revision 1.13  2014/07/30 17:26:41  eng
  Added Object_Saturation_Limit_Set() to allow the calling application to define
  a saturation limit for teh detector that overrides the built in default.

  Revision 1.12  2009/08/19 17:54:00  eng
  Tweaked to apply the MARGIN criteria instead to the object centroid position
  (xpos,ypos) rather than pixels used in object creation. This is because
  there's nothing to stop the finding-of-connected-pixels process going into
  the margin area in some cases.

  Revision 1.11  2009/08/11 14:24:32  cjm
  Turned down per-pixel logging.
  Reenabled rcsid.

  Revision 1.10  2009/08/06 13:34:03  eng
  Defining a "margin" around the image frame to prevent objects being made from
  pixels within a distance N pixels from the frame edge. This is primarily to
  avoid objects made from the "hot rows" at the bottom of RATCam images. NB: I
  would have thought however that such an object would have been rejected as
  very much non-stellar anyway later in the code, but so be it. (JMM)

  Revision 1.9  2009/06/12 10:32:06  cjm
  Added Object_Stellar_Ellipticity_Limit_Set/Stellar_Ellipticity_Limit to make
  the stellar ellipticity limit configurable.

  Revision 1.8  2009/01/30 15:20:31  cjm
  Replaced log bit levels with log_udp.h verbositys.

  Revision 1.7  2009/01/28 14:18:53  cjm
  Added extra parameters to logging routines for GLS support.

  Revision 1.6  2008/10/07 13:55:32  cjm
  Fixed object total counts measure by uncommenting increment.

  Revision 1.5  2008/10/07 13:17:59  cjm
  Changed some pixel logging from OBJECT to PIXEL bit.

  Revision 1.4  2008/10/07 12:57:32  eng
  Added HighPixel *curpix; (jmm)

  Revision 1.3  2008/10/01 15:53:36  eng
  Some invisible non-printing character was preventing compilation with error
  syntax error: line 1810, token = " ", char = '270'
  I fixed it by deleteing a couple of apparently blank lines.

  Revision 1.2  2008/09/30 19:42:44  eng
  changed object_jmm.h to object.h

  Revision 1.1  2008/09/30 19:26:27  eng
  Initial revision

  Revision 1.12.2.10  2008/09/26 12:51:05  eng
  Minor tweaks:
  * altered some default seeing flag names
  * in Object_Calculate_FWHM:
    - removed moffat & 2M methods
    - removed theta & diff in ellipticity calculation
    - moved indirections to object-specific constants

  Revision 1.12.2.9  2008/09/25 16:18:58  eng
  Object_List_Get_New has been renamed back to Object_List_Get. Also, as we now only
  use one blanket threshold, "thresh1" has been renamed back to "thresh". The extra args
  to Object_List_Get, namely thresh2 and the object counters, have been removed,
  putting it back in its original state. The object counters themselves are now local
  variables inside Object_List_Get and are not passed out of the function. Also tidied
  up some of the comments.

  Revision 1.12.2.8  2008/09/25 15:30:36  eng
  Added RJS' Object_Find_Peak to find the local peak of each individual object and
  then set a local-to-object 1/5th peak threshold (replacing the blanket "thresh2"
  of earlier versions).

  Revision 1.12.2.7  2008/09/25 13:27:30  eng
  Fix of memory leak whereby the call to Object_Calculate_FWHM only produced sane
  results if a mystery 'printf' was added just before or after the function call.
  The variable 'xy2I' in Object_Calculate_FWHM was uninitialised and the best
  explanation is that the extra printf happened to push the pointer into the right
  place. It is initialised now and everything seems to run now.

  Revision 1.12.2.6  2008/09/15 13:12:09  eng
  Changed fwhm_lt_dia_count --> usable_count
  Deleted Object_List_Get
  Added new default FWHM settings in case of bad data
  Added ellipticity to sizefwhm struct, for diagnostic purposes.
  In Object_List_Get_New:
    - Removed fwhmlist
    - Added stellarity & +ve fwhm to selection for fwhmarray:
        if ((obj_fwhm < obj_dia) && (w_object->is_stellar == 1) && (obj_fwhm > 0.0))
  In Object_Calculate_FWHM (sExtractor-adapted-code section):
    - Added default fwhm settings in case fwhm -ve or sex_d -ve.

  Revision 1.12.2.5  2008/09/09 11:08:10  eng
  This is basically the same two-threshold code as before but with a few extra tweaks
  based on/adapted from SExtractor code (not the manual). This uses only the top 80%
  of pixels, i.e. ignores the "wings", and then uses a factor to compensate.

  Revision 1.12.2.4  2008/08/30 15:37:21  eng
    Closing 1.12.2.3 - FWHM code was fixed, part typo error. New code from scratch worked but
    produced FWHMs that were now much larger than the original.

  Revision 1.12.2.3  2008/08/22 15:58:17  eng
  Working version of two-threshold ("TT") algorithm. Emphasis on applying the TT to the old
  2nd-moment ("2M")method of calculating FWHM again. Why? Because the Moffatt-fitting algorithm
  barfs out during a focus run, which admittedly may be a coding problem rather than
  fundamental method problem, AND it is very slow, which again could just be down to
  breadboard (not efficinet) coding practices. Looked again at the 2m version because
  it looked more robust and a lot faster. However this version still produces wrong
  FWHMs - a fault which might be down to the maths used for the 1st & 2nd moment calc.
  Therefore, going to try coding the 1st/2nd moments again from scratch rather than try to
  second-guess the original existing code. 1st moment will still be calculated in
  Object_List_Get_New & 2nd in Object_Calculate_FWHM, but the maths involved will be sourced
  from fresh locations.

  Revision 1.12.2.2  2008/06/05 13:56:51  eng
  *** empty log message ***

  Revision 1.12.2.1  2008/06/03 14:25:48  eng
  Another branch-off from 1.12, this time to persue the concept of just
  using two thresholds; one high-sigma to use when detecting objects,
  the other low-sigma to use when getting connected pixels. This way we
  get well-sampled objects out to the "wings", but only for bright objects, not
  wasting time doing it for every object a low-sigma value above the median.

  Revision 1.12  2008/05/29 22:26:09  eng
  Checked in by RJS to enable JMM to roll back on Monday any changes he disapproves of.
  THis is just after Jon has reverted teh calling parameters to the old 9 parameter format.

  Revision 1.11  2008/05/29 15:35:48  eng
  Having to check this version in unfinished because I need to look at an older
  version for RJS' benefit. If I just had a different file with a different name
  then this would be easy, but because the earlier version is buried in the depths
  of RCS with no way of even looking at it, I have to actually check out the older
  version just to look at it. This will apparently destroy the current version unless
  it's checked in, so to keep it safe, I have to temporarily check it in even though
  it's unfinished. RCS: useful in some ways, bloody annoying in others.

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


#include <float.h>
#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include "object.h"
#include "log_udp.h"

/* for new fwhm */
/*#include <gsl/gsl_sf.h> */


/* ------------------------------------------------------- */
/* hash definitions */
/* ------------------------------------------------------- */
/**
 * Values to use for seeing when something goes wrong (in pixels).
 * Should be larger than 10 pixels, so RCS thinks seeing is "bad".
 * RJS wants really large value so archive searches can differentiate
 * between real bad seeing and failed reductions.
 */
#define DEFAULT_BAD_SEEING         (999.0)       /* generic "bad seeing" flag                  */
#define DEFAULT_SEEING_NONSTELLAR  (988.0)       /* too elliptical                             */
#define DEFAULT_SEEING_TOOSMALL    (977.0)       /* less than 0.01 pixels                      */
#define DEFAULT_SEEING_TOOBIG      (966.0)       /* fitted moffat curve fwhm > object diameter */
#define DEFAULT_SEEING_SEXD_ZERO   (950.0)       /* if sex_d not > 0                           */
#define DEFAULT_SEEING_ZERO        (951.0)       /* if fwhm not > 0                            */

/**
 * Default upper limit of ellipticity for an object to be classed as 'stellar'.
 */
#define DEFAULT_STELLAR_ELLIP_LIMIT   (0.3)

/**
 * Default saturation limit (ADU) for the detector 
 */
#define DEFAULT_SATURATION_LIMIT   (63000)


#define MAX_N_FWHM            (17)          /* deliberately chosen to be odd */
#define MAX_N_FWHM_MID        (8)           /* median for 11 items with c-like counting from zero */


/**
 * Define a "margin" around the image frame to prevent objects being made from
 * pixels within a distance N pixels from the frame edge.
 */
#define MARGIN (5) /* pixels */




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
  float ellipticity;
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
  void (*Log_Handler)(char *sub_system,char *source_filename,char *function,int level,char *category,char *string);
  int (*Log_Filter)(char *sub_system,char *source_filename,char *function,int level,char *category);
  int Log_Filter_Level;
};



/* ------------------------------------------------------- */
/* internal variables */
/* ------------------------------------------------------- */
/**
 * Revision Control System identifier.
 */
static char rcsid[] = "$Id: object.c,v 1.17 2023-11-10 16:47:04 cjm Exp $";
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
/**
 * Upper limit of ellipticity for an object to be classed as 'stellar'.
 * @see #DEFAULT_STELLAR_ELLIP_LIMIT
 */
static float Stellar_Ellipticity_Limit = DEFAULT_STELLAR_ELLIP_LIMIT;
/**
 * Max good value of detector array in ADU before a source is considered to 
 * be saturated. 
 * @see #DEFAULT_SATURATION_LIMIT
 */
static float Saturation_Limit = DEFAULT_SATURATION_LIMIT;

/* ------------------------------------------------------- */
/* internal function declarations */
/* ------------------------------------------------------- */
static int Object_Find_Peak(int naxis1,int naxis2,int x,int y,float *image,Object *w_object);
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
  ___  _     _           _       _     _     _       ____      _   
 / _ \| |__ (_) ___  ___| |_    | |   (_)___| |_    / ___| ___| |_ 
| | | | '_ \| |/ _ \/ __| __|   | |   | / __| __|  | |  _ / _ \ __|
| |_| | |_) | |  __/ (__| |_    | |___| \__ \ |_   | |_| |  __/ |_ 
 \___/|_.__// |\___|\___|\__|   |_____|_|___/\__|   \____|\___|\__|
          |__/             
                     
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
int Object_List_Get(float *image,float image_median,int naxis1,int naxis2,float thresh,
			int npix,Object **first_object,int *sflag,float *seeing)
{
  Object *w_object = NULL;
  Object *last_object = NULL;
  Object *next_object = NULL;
  HighPixel *curpix;
  float fwhm = 0.0;
  float thresh2 = 0.0;                      /* individual object 2nd threshold (1/5th peak) to build object */
  int y,x,done,is_stellar;
  int local_peak_x,local_peak_y;	    /* Location of the peak returned by Object_Find_Peak() */
  int fwhmarray_size = 0;
  struct sizefwhm *fwhmarray = NULL;        /* array for objects whose fwhm is smaller than its diameter */
  int obj_area;                             /* number of pixels in object */
  float obj_fwhm;                           /* object fwhm in pixels */
  float obj_dia;                            /* object pseudo-diameter (pixels) */
  float obj_peak;			    /* ADU of brightest pixel in the image */
  int mid_posn;                             /* middle position of fwhmarray, to find median */
  int lower_mid_posn,upper_mid_posn;        /* array positions either side of median, for even-sized fwhmarray */
  float median_fwhm;                        /* median fwhm obtained from fwhmarray */


  /* Initialise object counters - now internal to this function only as of 1.12.2.9. 
     Note therefore that: initial_count > size_count > stellar_count > usable_count */
  int initial_count = 0;               /* initial count of all objects */
  int size_count = 0;                  /* objects bigger than size limit (currently 8 pixels) */
  int stellar_count = 0;               /* objects with ellipticity below limit (i.e. "stellar") */
  int usable_count = 0;                /* stellar objects where fwhm < diameter (calculated from size) */

  int i = 0; /* needed in logging */


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


  /* ensure top of list is set to NULL - assume it has not been allocated. */
  (*first_object) = NULL;



/*
                     _       __               _     _        _      
 ___ ___ __ _ _ _ __| |_    / _|___ _ _   ___| |__ (_)___ __| |_ ___
(_-</ -_) _` | '_/ _| ' \  |  _/ _ \ '_| / _ \ '_ \| / -_) _|  _(_-<
/__/\___\__,_|_| \__|_||_| |_| \___/_|   \___/_.__// \___\__|\__/__/
                                                 |__/               
*/

#if LOGGING > 0
  Object_Log("object","object.c","Object_List_Get",LOG_VERBOSITY_TERSE,NULL,"Searching for objects.");
#endif



  /* ---------------------- */
  /* RUN THROUGH ALL PIXELS */
  /* ---------------------- */



#if LOGGING > 7
  Object_Log_Format("object","object.c","Object_List_Get",LOG_VERBOSITY_INTERMEDIATE,NULL,"(AGD) All pixels above threshold %.2f",thresh);
#endif



  for(y=0;y<naxis2;y++)
    {
      for(x=0;x<naxis1;x++)
	{

#if LOGGING > 9
	  Object_Log_Format("object","object.c","Object_List_Get",LOG_VERBOSITY_VERY_VERBOSE,NULL,"searching pixel %d,%d.",x,y);
#endif

	  /* ---------------------------- */
	  /* IF PIXEL ABOVE THRESHOLD -1- */
	  /* ---------------------------- */
    
	  if(image[(y*naxis1)+x] > thresh)  
	    {
	      initial_count++;
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
	      Object_Log_Format("object","object.c","Object_List_Get",LOG_VERBOSITY_VERY_VERBOSE,NULL,
				"allocated w_object (%p).",w_object);
#endif

	      w_object->nextobject=NULL;
	      w_object->highpixel = NULL;
	      w_object->last_hp = NULL;
	      if((*first_object)==NULL)
		{
		  (*first_object) = w_object;
		  last_object = w_object;

#if LOGGING > 10
		  Object_Log_Format("object","object.c","Object_List_Get",LOG_VERBOSITY_VERY_VERBOSE,NULL,
				    "set first_object to (%p).",(*first_object));
#endif

		}
	      else
		{
		  last_object->nextobject = w_object;
		  last_object = w_object;
		}
	      w_object->objnum = initial_count;


#if LOGGING > 3
	      Object_Log_Format("object","object.c","Object_List_Get",LOG_VERBOSITY_INTERMEDIATE,NULL,
				"found start of object at %d,%d,%.2f",x,y,image[(y*naxis1)+x]);
#endif


	      /* --------------------------------------------------------- */
	      /* GET ALL CONNECTED PIXELS ABOVE LOCAL 1/5th PEAK THRESHOLD */
	      /* --------------------------------------------------------- */

	      /*
		initialise stats
		----------------
	      */
	      w_object->total=0;
	      w_object->xpos=0;
	      w_object->ypos=0;
	      w_object->peak=0;
	      w_object->numpix=0;

	      /*
		find object peak value
		----------------------
	      */
#if LOGGING > 7
	      Object_Log_Format("object","object.c","Object_List_Get",LOG_VERBOSITY_INTERMEDIATE,NULL,
				"(AGD) calling Object_Find_Peak to find local 1/5th peak value",x,y);
#endif

	      Object_Find_Peak(naxis1,naxis2,x,y,image,w_object);

	      /* 
		set local peak coordinates
		--------------------------
		Do not reassign x,y to the peak because once we have extracted
		this source we want to go back to searching from where we left
		off, so we keep x,y to be where we first found this object. It
		is however much more efficient to start the extraction from the
		peak, so we create these local_peak_x,local_peak_y coords and
		start from there.
	      */
	      local_peak_x = w_object->xpos;
	      local_peak_y = w_object->ypos;


	      /*
		set 1/5th peak level
		--------------------
		Object_Find_Peak does not background subtract
	      */
	      thresh2 = image_median + ( (w_object->peak-image_median) / 5); 

#if LOGGING > 7
	      Object_Log_Format("object","object.c","Object_List_Get",LOG_VERBOSITY_INTERMEDIATE,NULL,
				"(AGD) Found object peak at %d,%d,%.2f so setting thresh2 = median + (peak/5) = %.2f",
				local_peak_x,local_peak_y,image[(local_peak_y*naxis1)+local_peak_x],thresh2);
#endif
	      
	      /* 
		check if thresh2 > thresh
		-------------------------
		You must extract at least down to thresh. Never let thresh2 be
		above thresh otherwise you will rediscover this object a second
		time and extract its halo as a second object after you have
		extracted the core above thresh2 as a first obejct.
	      */
	      if (thresh2 > thresh){
		thresh2 = thresh;    

#if LOGGING > 7
		Object_Log_Format("object","object.c","Object_List_Get",LOG_VERBOSITY_INTERMEDIATE,NULL,
				  "(AGD) thresh2 must be below thresh, but here thresh2 > thresh, so setting thresh2 = thresh");
#endif
	      }


	      /* 
		 reset stats
		 -----------
		 as if we have not run the Object_List_Get_Connected_Pixels() above
	      */
	      w_object->xpos=0;
	      w_object->ypos=0;
	      w_object->peak=0;
	      w_object->numpix=0;



#if LOGGING > 7
		Object_Log_Format("object","object.c","Object_List_Get",LOG_VERBOSITY_INTERMEDIATE,NULL,
				  "(AGD) calling Object_List_Get_Connected_Pixels to build object now, using all appropriate pixels");
#endif

	      if(!Object_List_Get_Connected_Pixels(naxis1,naxis2,image_median,x,y,thresh2,image,
						   w_object))
		{
		  return FALSE;
		}

	    }/* end if threshold exceeded for image[x,y] */
       }/* end for on x */
    }/* end for on y */




#if LOGGING > 0
  Object_Log_Format("object","object.c","Object_List_Get",LOG_VERBOSITY_TERSE,NULL,"Found %d objects.",
		    initial_count);
#endif




  /*
     _  __                  _     _        _      
    (_)/ _|  _ _  ___   ___| |__ (_)___ __| |_ ___
    | |  _| | ' \/ _ \ / _ \ '_ \| / -_) _|  _(_-<
    |_|_|   |_||_\___/ \___/_.__// \___\__|\__/__/
                               |__/               
  */
  if(initial_count == 0)
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




  /*
                                      _     _        _      
     _ _ ___ _ __  _____ _____    ___| |__ (_)___ __| |_ ___
    | '_/ -_) '  \/ _ \ V / -_)  / _ \ '_ \| / -_) _|  _(_-<
    |_| \___|_|_|_\___/\_/\___|  \___/_.__// \___\__|\__/__/
                                         |__/               

    Go through list of objects, getting rid of:
    - objects with less than npix
    - objects where xpos,ypos are within MARGIN pixels of the frame edge
  */



#if LOGGING > 0
  Object_Log("object","object.c","Object_List_Get",LOG_VERBOSITY_VERY_VERBOSE,NULL,"Finding useful objects.");
#endif




  /* ------------------------------------------------- */
  /* IF FIRST OBJECT TOO SMALL OR IN MARGIN, DELETE IT */
  /* ------------------------------------------------- */
  
  w_object = (*first_object);
  done = FALSE;
  while(done == FALSE){
    if ((w_object->numpix >= npix) 
	&& (w_object->xpos > MARGIN) && (w_object->xpos <(naxis1-MARGIN)) 
	&& (w_object->ypos > MARGIN) && (w_object->ypos <(naxis2-MARGIN)))
      done = TRUE;                           /* we're done */
    else {                                   /* otherwise */




#if LOGGING > 5
      Object_Log_Format("object","object.c","Object_List_Get",LOG_VERBOSITY_VERY_VERBOSE,NULL,
			"deleting object(1) at %.2f,%.2f(%d).",
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
      sprintf(Object_Error_String,"Object_List_Get: All objects were too small.");
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


  /* -------------------- */
  /* SET NEW FIRST OBJECT */
  /* -------------------- */



#if LOGGING > 10
  Object_Log_Format("object","object.c","Object_List_Get",LOG_VERBOSITY_VERY_VERBOSE,NULL,
		    "first_object (%p) set from w_object (%p).",(*first_object),w_object);
#endif



  (*first_object) = w_object;
  w_object->objnum=1;
  last_object = (*first_object);


#if LOGGING > 5
  Object_Log_Format("object","object.c","Object_List_Get",LOG_VERBOSITY_VERY_VERBOSE,NULL,
		    "object %d at %.2f,%.2f(%d) is ok(1).",
		    w_object->objnum,w_object->xpos,w_object->ypos,w_object->numpix);
#endif




  w_object = (*first_object)->nextobject;
  size_count = 1;
  done = FALSE;


  /* ------------------------------------------------------------------- */
  /* GO THROUGH REST OF OBJECT LIST, DELETING OBJECTS WITH NUMPIX < NPIX */
  /* OR IF OBJECT IN MARGIN                                              */
  /* ------------------------------------------------------------------- */

  while(w_object != NULL)
    {
      next_object=w_object->nextobject;         /* take copy of next object to go to */
      if((w_object->numpix < npix)
	 || (w_object->xpos < MARGIN) || (w_object->xpos >(naxis1-MARGIN)) 
	 || (w_object->ypos < MARGIN) || (w_object->ypos >(naxis2-MARGIN)))
	{




#if LOGGING > 5
	  Object_Log_Format("object","object.c","Object_List_Get",LOG_VERBOSITY_VERY_VERBOSE,NULL,
			    "deleting object(2) at %.2f,%.2f(%d).",
			    w_object->xpos,w_object->ypos,w_object->numpix);
#endif



	  Object_Free(&w_object);
	}
      else
	{
	  size_count++;
	  last_object->nextobject = w_object;   /* tell last object this is its next object */
	  w_object->objnum = size_count;        /* set objects number */
	  last_object = w_object;               /* set the last object in the list to be this object */




#if LOGGING > 5
	  Object_Log_Format("object","object.c","Object_List_Get",LOG_VERBOSITY_VERY_VERBOSE,NULL,
			    "object %d at %.2f,%.2f(%d) is ok(2).",
			    w_object->objnum,w_object->xpos,w_object->ypos,w_object->numpix);
#endif




	}
      w_object = next_object;                   /* change to next object */
    }
  last_object->nextobject=NULL;




  /* -------------------------------------------------- */
  /* EXTRA DEBUG - LIST CONNECTED PIXELS IN ALL OBJECTS */
  /* -------------------------------------------------- */

#if LOGGING > 12
  w_object = (*first_object);
  while(w_object != NULL)
  {
    Object_Log_Format("object","object.c","Object_List_Get",LOG_VERBOSITY_VERY_VERBOSE,NULL,
		      "Printing pixels for object %d at %.2f,%.2f(%d).",
		      w_object->objnum,w_object->xpos,w_object->ypos,w_object->numpix);
    curpix = w_object->highpixel;
    while(curpix != NULL)
      {
	Object_Log_Format("object","object.c","Object_List_Get",LOG_VERBOSITY_VERY_VERBOSE,NULL,
			  "Printing pixels:object:%d pixel %d,%d value %.2f.",
			  w_object->objnum,curpix->x,curpix->y,curpix->value);
	 curpix = curpix->next_pixel;           /* goto next pixel */
      }
      w_object = w_object->nextobject;	        /* goto next object */
  }
#endif



  /*
                     _         _ _    _          __       _     _        _   
     __ _ _ ___ __ _| |_ ___  | (_)__| |_   ___ / _|  ___| |__ (_)___ __| |_ 
    / _| '_/ -_) _` |  _/ -_) | | (_-<  _| / _ \  _| / _ \ '_ \| / -_) _|  _|
    \__|_| \___\__,_|\__\___| |_|_/__/\__| \___/_|   \___/_.__// \___\__|\__|
                                                             |__/            
     _____      ___  _ __  __    
    | __\ \    / / || |  \/  |___
    | _| \ \/\/ /| __ | |\/| (_-<
    |_|   \_/\_/ |_||_|_|  |_/__/
                             
  */
  
#if LOGGING > 0
  Object_Log("object","object.c","Object_List_Get",LOG_VERBOSITY_INTERMEDIATE,NULL,"Finding FWHM of objects.");
#endif


  /* ---------------- */
  /* SET FIRST OBJECT */
  /* ---------------- */
  w_object = (*first_object);
#if LOGGING > 10
  Object_Log_Format("object","object.c","Object_List_Get",LOG_VERBOSITY_VERY_VERBOSE,NULL,
		    "w_object (%p) set from first_object (%p).",w_object,(*first_object));
#endif


  /* --------------------------- */
  /* RUN THROUGH LIST OF OBJECTS */
  /* CALCULATING FWHM            */
  /* --------------------------- */
  while(w_object != NULL){
#if LOGGING > 5
    Object_Log_Format("object","object.c","Object_List_Get",LOG_VERBOSITY_VERBOSE,NULL,
		      "Calculating FWHM for object (%d) at %.2f,%.2f.",
		      w_object->objnum,w_object->xpos,w_object->ypos);
#endif
    

    /* calculate FWHM of object */
    /* ------------------------ */
    Object_Calculate_FWHM(w_object,image_median,&is_stellar,&fwhm);



#if LOGGING > 5
    Object_Log_Format("object","object.c","Object_List_Get",LOG_VERBOSITY_VERBOSE,NULL,
		      "object (%d) at %.2f,%.2f has FWHM %.2f pixels and is_stellar = %d.",
		      w_object->objnum,w_object->xpos,w_object->ypos,fwhm,is_stellar);
#endif


    /* if object stellar, increment stellar_count */
    /* ------------------------------------------ */
    if(is_stellar)
      stellar_count++;

    w_object = w_object->nextobject;		
  }









  /*
             _       __ _           _   _____      ___  _ __  __ 
     __ __ _| |__   / _(_)_ _  __ _| | | __\ \    / / || |  \/  |
    / _/ _` | / _| |  _| | ' \/ _` | | | _| \ \/\/ /| __ | |\/| |
    \__\__,_|_\__| |_| |_|_||_\__,_|_| |_|   \_/\_/ |_||_|_|  |_|
                                                             
  */


  /* If any stellar objects at all */
  /* ----------------------------- */
#if LOGGING > 0
  Object_Log("object","object.c","Object_List_Get",LOG_VERBOSITY_INTERMEDIATE,NULL,"Calculating final seeing.");
  Object_Log_Format("object","object.c","Object_List_Get",LOG_VERBOSITY_INTERMEDIATE,NULL,
		    "Number of stellar objects: %d", stellar_count);
#endif


  if(stellar_count > 0) {
#if LOGGING > 0
    Object_Log("object","object.c","Object_List_Get",LOG_VERBOSITY_VERY_VERBOSE,NULL,"Creating fwhmarray");
#endif

    fwhmarray = (struct sizefwhm *) malloc((stellar_count) * sizeof(struct sizefwhm));

    /* populate array of structs from w_object, */
    /* but ONLY IF:                             */
    /*   (a) object is stellar,                 */
    /*   (b) fwhm < diameter of object,         */ /* <-- this will remove any 900+ codes */
    /*   (c) fwhm > 0                           */
    /*   (d) no pixel exceeds Saturation_Limit  */
    /* ---------------------------------------- */

    w_object = (*first_object);                                      /* set w_object to first obj in list */
    while(w_object != NULL){                                         /* start running through all objects */
      obj_area = w_object->numpix;                                   /* area == numpix                    */
      obj_fwhm = (w_object->fwhmx + w_object->fwhmy)/2.0;            /* calc mean fwhm                    */  
      obj_dia = sqrt( 1.2732 * obj_area);                            /* pseudo-diameter. 1.2732 = 4/pi    */
      obj_peak = w_object->peak;				     /* Brightest pixel in object         */


      /* if fwhm < dia & is stellar & fwhm is +ve & not saturated */
      if ((obj_fwhm < obj_dia) && (w_object->is_stellar == 1) && (obj_fwhm > 0.0) && (obj_peak < Saturation_Limit) ) { 
	fwhmarray[usable_count].numpix = obj_area;
	fwhmarray[usable_count].fwhm = obj_fwhm;
	fwhmarray[usable_count].objnum = w_object->objnum;
	fwhmarray[usable_count].xpos = w_object->xpos;
	fwhmarray[usable_count].ypos = w_object->ypos;
	fwhmarray[usable_count].ellipticity = w_object->ellipticity;
	usable_count++;                                           /* increment counter of usable objects */
      }


#if LOGGING > 0
      Object_Log_Format("object","object.c","Object_List_Get",LOG_VERBOSITY_VERY_VERBOSE,NULL,
			"Object %d\t%f\t%f\t(%d)",w_object->objnum,obj_fwhm,obj_dia,usable_count);
#endif
      w_object = w_object->nextobject;                               /* go to next object */		
    }
      

#if LOGGING > 0
    Object_Log_Format("object","object.c","Object_List_Get",LOG_VERBOSITY_VERBOSE,NULL,
		      "Number of usable objects: %d", usable_count);
#endif


    /* If some stellar fwhms were < diameter, then */
    /* find median fwhm of top N brightest objects */
    /* in this array                               */
    /* ------------------------------------------- */
    if ( usable_count > 0 ){

      /* trim fwhmarray array to exact number of objects (fwhm_lt_dia_count) */
      fwhmarray = (struct sizefwhm *) realloc(fwhmarray,(usable_count)*sizeof(struct sizefwhm));
      
      /* set standard array size descriptor (necessary for later on) */
      fwhmarray_size = (int) (usable_count);


#if LOGGING > 0
      Object_Log_Format("object","object.c","Object_List_Get",LOG_VERBOSITY_VERY_VERBOSE,NULL,
			"original object list\n[n] (objnum)\tnumpix\tfwhm\tellip\n--------------------");
      for (i=0;i<fwhmarray_size;i++)
	Object_Log_Format("object","object.c","Object_List_Get",LOG_VERBOSITY_VERY_VERBOSE,NULL,
			  "[%d] (%d)\t%d\t%f\t%f",
			  i,fwhmarray[i].objnum,fwhmarray[i].numpix,fwhmarray[i].fwhm,fwhmarray[i].ellipticity);
#endif

      /* sort array (LARGEST FIRST) by 1st struct member (numpix) */
      qsort (fwhmarray, fwhmarray_size, sizeof(struct sizefwhm), sizefwhm_cmp_by_numpix);

#if LOGGING > 0
      Object_Log_Format("object","object.c","Object_List_Get",LOG_VERBOSITY_VERY_VERBOSE,NULL,
			"sorted by numpix\n[n] (objnum)\tnumpix\tfwhm\tellip\n--------------------");
      for (i=0;i<fwhmarray_size;i++)
	Object_Log_Format("object","object.c","Object_List_Get",LOG_VERBOSITY_VERY_VERBOSE,NULL,
			  "[%d] (%d)\t%d\t%f\t%f",
			  i,fwhmarray[i].objnum,fwhmarray[i].numpix,fwhmarray[i].fwhm,fwhmarray[i].ellipticity);
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
	Object_Log_Format("object","object.c","Object_List_Get",LOG_VERBOSITY_VERY_VERBOSE,NULL,"N > MAX_N_FWHM:");
#endif

	/* sort array by 2nd struct member (fwhm) SMALLEST FIRST */
	qsort (fwhmarray, fwhmarray_size, sizeof(struct sizefwhm), sizefwhm_cmp_by_fwhm);
#if LOGGING > 0
	Object_Log_Format("object","object.c","Object_List_Get",LOG_VERBOSITY_VERY_VERBOSE,NULL,
		 "truncated & sorted by fwhm\n[n] (objnum)\tnumpix\tfwhm\txpos\typos\tellip\n--------------------");
	for (i=0;i<fwhmarray_size;i++)
	  Object_Log_Format("object","object.c","Object_List_Get",LOG_VERBOSITY_VERY_VERBOSE,NULL,
			    "fwhmsort: [%d] (%d)\t%d\t%f\t%f\t%f\t%f",
			    i,fwhmarray[i].objnum,
			    fwhmarray[i].numpix,fwhmarray[i].fwhm,
			    fwhmarray[i].xpos,fwhmarray[i].ypos,
			    fwhmarray[i].ellipticity);
#endif
	
	


	/* find median */
	mid_posn = MAX_N_FWHM_MID;
	median_fwhm = fwhmarray[mid_posn].fwhm;
      }

      /* otherwise, if fwhmarray_size < MAX_N_FWHM */
      else {

#if LOGGING > 0
	Object_Log_Format("object","object.c","Object_List_Get",LOG_VERBOSITY_VERY_VERBOSE,NULL,"N < MAX_N_FWHM:");
#endif
	/* sort array by 2nd struct member (fwhm) SMALLEST FIRST */
	qsort (fwhmarray, fwhmarray_size, sizeof(struct sizefwhm), sizefwhm_cmp_by_fwhm);
#if LOGGING > 0
	Object_Log_Format("object","object.c","Object_List_Get",LOG_VERBOSITY_VERY_VERBOSE,NULL,
			  "sorted by FWHM\n[n] (objnum)\tnumpix\tfwhm\txpos\typos\tellip\n-----------------");
	for (i=0;i<fwhmarray_size;i++)
	  Object_Log_Format("object","object.c","Object_List_Get",LOG_VERBOSITY_VERY_VERBOSE,NULL,
			    "fwhmsort: [%d] (%d)\t%d\t%f\t%f\t%f\t%f",
			    i,fwhmarray[i].objnum,
			    fwhmarray[i].numpix,fwhmarray[i].fwhm,
			    fwhmarray[i].xpos,fwhmarray[i].ypos,
			    fwhmarray[i].ellipticity);
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
	Object_Log_Format("object","object.c","Object_List_Get",LOG_VERBOSITY_VERY_VERBOSE,NULL,
			  "median_fwhm = [%d,%d] (%d,%d) %f",
			  lower_mid_posn,upper_mid_posn,
			  fwhmarray[lower_mid_posn].objnum,fwhmarray[upper_mid_posn].objnum,
			  median_fwhm);
      else /* if ODD */
	Object_Log_Format("object","object.c","Object_List_Get",LOG_VERBOSITY_VERY_VERBOSE,NULL,
			  "median_fwhm = [%d] (%d) %f",
			  mid_posn,fwhmarray[mid_posn].objnum,median_fwhm);
#endif

      /* Set seeing to median_fwhm */
      /* ------------------------- */
      (*seeing) = median_fwhm;
      (*sflag) = 0;                                  /* set sflag to zero (obviously) */
      
      
      /* If the seeing is less than 0.01 */
      /* ------------------------------- */
      if ((*seeing)<=0.01){
	(*seeing) = DEFAULT_SEEING_TOOSMALL;        /* set the seeing to DEFAULT_SEEING_TOOSMALL */
	(*sflag) = 1;                               /* and set sflag to show the seeing was fudged */
      }

    } /* end of "if usable_count > 0" */


    /* if usable_count = 0 (i.e. if there are no objects whose fwhm < dia) */
    /* ------------------------------------------------------------------- */
    else {
      (*seeing) = DEFAULT_SEEING_TOOBIG;  
      (*sflag) = 1;                       
    }



  } /* end 'if any fwhms at all' */



  /* If NO fwhms at all */
  /* ------------------ */
  else {
#if LOGGING > 0
    Object_Log_Format("object","object.c","Object_List_Get",LOG_VERBOSITY_VERY_VERBOSE,NULL,
			"No objects found - defaulting seeing to BAD_SEEING");
#endif
    (*seeing) = DEFAULT_BAD_SEEING;                  /* set the seeing to DEFAULT_BAD_SEEING (pixels) */
    (*sflag) = 1;                                    /* and set sflag to show the seeing was fudged */
  }


  /* ----------- */
  /* FREE MEMORY */
  /* ----------- */
  if(fwhmarray != NULL)
    free(fwhmarray);



#if LOGGING > 0
  Object_Log_Format("object","object.c","Object_List_Get",LOG_VERBOSITY_VERBOSE,NULL,
		    "number of objects > %d pixels = %d",npix,size_count);
  Object_Log_Format("object","object.c","Object_List_Get",LOG_VERBOSITY_VERBOSE,NULL,
		    "number of objects identified as stellar = %d",stellar_count);
  Object_Log_Format("object","object.c","Object_List_Get",LOG_VERBOSITY_VERBOSE,NULL,
		    "number of stellar objects with fwhm < dia (\"usable\") = %d",usable_count);
  if ((*sflag)==0)
    {
      Object_Log_Format("object","object.c","Object_List_Get",LOG_VERBOSITY_INTERMEDIATE,NULL,
			"seeing derived from stellar sources = %.2f pixels.",(*seeing));
    }
  else
    {
      Object_Log_Format("object","object.c","Object_List_Get",LOG_VERBOSITY_INTERMEDIATE,NULL,
			"Unable to derive seeing, faking result = %.2f pixels.",(*seeing));
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


/**
 * Set the stellar ellipticity limit. If the computed ellipticity of an object is above the limit,
 * the object is flagged non-stellar.
 * @param limit The ellipticity limit. This must be a positive number.
 * @return The routine returns TRUE on success and false on failure.
 */
int Object_Stellar_Ellipticity_Limit_Set(float limit)
{
	if(limit <= 0.0)
	{
		Object_Error_Number = 8;
		sprintf(Object_Error_String,"Object_Stellar_Ellipticity_Limit_Set:ellipticity %.2f out of range.",
			limit);
		return FALSE;
	}
	Stellar_Ellipticity_Limit = limit;
	return TRUE;
}

/**
 * Set the detector saturation limit. All detected sources are analyzed and ertuned in the
 * list. This saturation threshold is only used to prevent saturated stars from being
 * included in the overall median seeing estimate for the frame.
 * @param saturation The detector saturation in ADU. This must be a positive number.
 * @return The routine returns TRUE on success and false on failure.
 */
int Object_Saturation_Limit_Set(float saturation)
{
	if(saturation <= 0.0)
	{
		Object_Error_Number = 16;
		sprintf(Object_Error_String,"Object_Saturation_Limit_Set:saturation %.2f out of range.", saturation);
		return FALSE;
	}
	Saturation_Limit = saturation;
	return TRUE;
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
 * @param sub_system The sub system. Can be NULL.
 * @param source_file The source filename. Can be NULL.
 * @param function The function calling the log. Can be NULL.
 * @param level At what level is the log message (TERSE/high level or VERBOSE/low level), 
 *         a valid member of LOG_VERBOSITY.
 * @param category What sort of information is the message. Designed to be used as a filter. Can be NULL.
 * @param format A string, with formatting statements the same as fprintf would use to determine the type
 * 	of the following arguments.
 * @see #Object_Log
 * @see #Object_Buff
 * @see #OBJECT_ERROR_STRING_LENGTH
 */
void Object_Log_Format(char *sub_system,char *source_filename,char *function,int level,char *category,char *format,...)
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
      if(Log_Data.Log_Filter(sub_system,source_filename,function,level,category) == FALSE)
	return;
    }
  /* format the arguments */
  va_start(ap,format);
  vsprintf(Object_Buff,format,ap);
  va_end(ap);
  /* call the log routine to log the results */
  Object_Log(sub_system,source_filename,function,level,category,Object_Buff);
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
 * @param sub_system The sub system. Can be NULL.
 * @param source_file The source filename. Can be NULL.
 * @param function The function calling the log. Can be NULL.
 * @param level At what level is the log message (TERSE/high level or VERBOSE/low level), 
 *         a valid member of LOG_VERBOSITY.
 * @param category What sort of information is the message. Designed to be used as a filter. Can be NULL.
 * @param string The message to log.
 * @see #Log_Data
 */
void Object_Log(char *sub_system,char *source_filename,char *function,int level,char *category,char *string)
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
      if(Log_Data.Log_Filter(sub_system,source_filename,function,level,category) == FALSE)
	return;
    }
  /* We can log the message */
  (*Log_Data.Log_Handler)(sub_system,source_filename,function,level,category,string);
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
void Object_Set_Log_Handler_Function(void (*log_fn)(char *sub_system,char *source_filename,char *function,
						    int level,char *category,char *string))
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
void Object_Set_Log_Filter_Function(int (*filter_fn)(char *sub_system,char *source_filename,char *function,
						     int level,char *category))
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
 * @param sub_system The sub system. Can be NULL.
 * @param source_file The source filename. Can be NULL.
 * @param function The function calling the log. Can be NULL.
 * @param level At what level is the log message (TERSE/high level or VERBOSE/low level), 
 *         a valid member of LOG_VERBOSITY.
 * @param category What sort of information is the message. Designed to be used as a filter. Can be NULL.
 * @param string The log message to be logged. 
 * @see #Log_Handler
 */
void Object_Log_Handler_Stdout(char *sub_system,char *source_filename,char *function,
			       int level,char *category,char *string)
{
  char time_string[32];

  if(string == NULL)
    return;
  Object_Get_Current_Time_String(time_string,32);
  fprintf(stdout,"%s %s:%s\n",time_string,function,string);
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
 * @param sub_system The sub system. Can be NULL.
 * @param source_file The source filename. Can be NULL.
 * @param function The function calling the log. Can be NULL.
 * @param level At what level is the log message (TERSE/high level or VERBOSE/low level), 
 *         a valid member of LOG_VERBOSITY.
 * @param category What sort of information is the message. Designed to be used as a filter. Can be NULL.
 * @return The routine returns TRUE if the level is less than or equal to the Log_Data.Log_Filter_Level,
 * 	otherwise it returns FALSE.
 * @see #Log_Data
 */
int Object_Log_Filter_Level_Absolute(char *sub_system,char *source_filename,char *function,int level,char *category)
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
 * @param sub_system The sub system. Can be NULL.
 * @param source_file The source filename. Can be NULL.
 * @param function The function calling the log. Can be NULL.
 * @param level At what level is the log message (TERSE/high level or VERBOSE/low level), 
 *         a valid member of LOG_VERBOSITY.
 * @param category What sort of information is the message. Designed to be used as a filter. Can be NULL.
 * @return The routine returns TRUE if the level has bits set that are also set in the 
 * 	Log_Data.Log_Filter_Level, otherwise it returns FALSE.
 * @see #Log_Data
 */
int Object_Log_Filter_Level_Bitwise(char *sub_system,char *source_filename,char *function,int level,char *category)
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
  

  /* Note the threshold value passed to this function is actually what's referred to as "thresh2" outside of this
   function. It's still called "thresh" in here though. */


  /* ------------- */
  /* SET VARIABLES */
  /* ------------- */
  int x1,y1,cx,cy;                            /* Don't know what these do */
  HighPixel *temp_hp = NULL;
  HighPixel *curpix = NULL;
  struct Point_Struct *point_list = NULL;
  struct Point_Struct *last_point = NULL;
  int point_count=0;

  float SumXI = 0.0;                   /* Running totals for moment calculation */
  float SumYI = 0.0;                   /* Running totals for moment calculation */
  float SumI = 0.0;                    /* Running totals for moment calculation */

  /* float orig_pixelvalue = 0.0; */

  /* ------------------------------- */
  /* ADD FIRST POINT TO BE PROCESSED */
  /* ------------------------------- */

#if LOGGING > 9
  Object_Log_Format("object","object.c","Object_List_Get_Connected_Pixels",LOG_VERBOSITY_VERY_VERBOSE,NULL,
		    "adding point %d,%d to list.",x,y);
#endif

#if LOGGING > 7
  Object_Log_Format("object","object.c","Object_List_Get_Connected_Pixels",LOG_VERBOSITY_INTERMEDIATE,NULL,
		    "(AGD) adding first point (%d,%d,%.2f) to new object (count %d)",
		    x,y,image[(y*naxis1)+x],w_object->numpix);
#endif
  




  if(!Point_List_Add(&point_list,&point_count,&last_point,x,y))
    return FALSE;
  




  /* -------------------------------- */
  /* RUN THROUGH POINTS ON POINT LIST */
  /* -------------------------------- */
#if LOGGING > 7
  Object_Log_Format("object","object.c","Object_List_Get_Connected_Pixels",LOG_VERBOSITY_INTERMEDIATE,NULL,
		    "(AGD) running through points on point list");
#endif

  while(point_count > 0){
    

    /* start of per pixel stuff */
    /* ------------------------ */
    cx = point_list->x;
    cy = point_list->y;
    

    /* if pixel value above threshold */
    /* ------------------------------ */
    if (image[(cy*naxis1)+cx] > thresh){     
      
 
#if LOGGING > 9
      Object_Log_Format("object","object.c","Object_List_Get_Connected_Pixels",LOG_VERBOSITY_VERY_VERBOSE,NULL,
			"(AGD) pixel %d,%d,%f above thresh2 (%f) (count %d)",
			cx,cy,image[(cy*naxis1)+cx],thresh,point_count);
#endif


#if LOGGING > 9
      Object_Log_Format("object","object.c","Object_List_Get_Connected_Pixels",LOG_VERBOSITY_VERY_VERBOSE,NULL,
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
      w_object->last_hp->value=image[(cy*naxis1)+cx] - image_median;  /* leave as is for now but return to this
									 - (later) decided it's OK */




      /* important for recursion - stops infinite loops */
#if LOGGING > 9
      orig_pixelvalue = image[(cy*naxis1)+cx];
#endif	

      /*image[(cy*naxis1)+cx]=0.0; */  /* (AGD) Change to stop looping with negative thresh2 (11/4/12) */

      image[(cy*naxis1)+cx]=-1e9;


/* #if LOGGING > 7 */
/*       Object_Log_Format("object","object.c","Object_List_Get_Connected_Pixels",LOG_VERBOSITY_INTERMEDIATE,NULL, */
/* 			"(AGD) anti-recursion: zeroing point %d,%d,%f -> %d,%d,%f (count %d)", */
/* 			cx,cy,orig_pixelvalue, */
/* 			cx,cy,image[(cy*naxis1)+cx],point_count); */
/* #endif	 */
	



      /* end of per-pixel stuff: do some overall w_object stats */
      curpix = w_object->last_hp;

       w_object->total=(w_object->total)+(temp_hp->value);
/*       w_object->xpos=(w_object->xpos)+((temp_hp->x)*(temp_hp->value));  */
/*       w_object->ypos=(w_object->ypos)+((temp_hp->y)*(temp_hp->value)); */

      SumXI += (temp_hp->x)*(temp_hp->value);
      SumYI += (temp_hp->y)*(temp_hp->value);
      SumI += temp_hp->value;

      if ((w_object->peak)<(temp_hp->value))
	w_object->peak=(temp_hp->value);  

      w_object->numpix ++;


#if LOGGING > 7
      /* print 1-5 then every 10000th point */
      if ((w_object->numpix) <= 5)
	Object_Log_Format("object","object.c","Object_List_Get_Connected_Pixels",LOG_VERBOSITY_INTERMEDIATE,NULL,
			  "(AGD) object first 5: point (%d,%d,%.2f) > thresh2 (%.2f), adding to object (size %d)",
			  cx,cy,w_object->last_hp->value,thresh,(w_object->numpix));
	    

      if (((w_object->numpix) % 10000 ) == 0)
	Object_Log_Format("object","object.c","Object_List_Get_Connected_Pixels",LOG_VERBOSITY_INTERMEDIATE,NULL,
			  "(AGD) object runaway: point (%d,%d,%.2f) > thresh2 (%.2f), adding to object (size %d)",
			  cx,cy,w_object->last_hp->value,thresh,(w_object->numpix));
#endif	









      /* start of recursive stuff */
      for (x1 = cx-1; x1<=cx+1; x1++){
	for (y1 = cy-1; y1<=cy+1; y1++){
	  if (x1 >= naxis1 || y1 >= naxis2 || x1<0 || y1<0)  
	    continue;                                           /* set a flag here to say crap object? */
	  if (image[(y1*naxis1)+x1] > thresh){
	    /* add this point to be processed */
#if LOGGING > 9
	    Object_Log_Format("object","object.c","Object_List_Get_Connected_Pixels",LOG_VERBOSITY_VERY_VERBOSE,NULL,
		  	      "adding point %d,%d to list.",x1,y1);
#endif



	    if(!Point_List_Add(&point_list,&point_count,&last_point,x1,y1))
	      return FALSE;
	  }
	}/* end for on y1 */
      }/* end for on x1 */
    }/* end of if pixel value above threshold */



    /* log if already added to object i.e./or if pixel value below threshold */
    /* --------------------------------------------------------------------- */
    else {




#if LOGGING > 9
      Object_Log_Format("object","object.c","Object_List_Get_Connected_Pixels",LOG_VERBOSITY_VERY_VERBOSE,NULL,
			"pixel %d,%d already added to object, ignoring.",cx,cy);
#endif



    } /* end of log if pixel value below threshold */
    


    /* delete processed point */
    /* ---------------------- */
#if LOGGING > 9
    Object_Log_Format("object","object.c","Object_List_Get_Connected_Pixels",LOG_VERBOSITY_VERY_VERBOSE,NULL,
		      "deleting point %d,%d from list.",cx,cy);
#endif



    if(!Point_List_Remove_Head(&point_list,&point_count))
      return FALSE;

  }/* end while on point list */
  



  /* -------------------------------- */
  /* CALCULATE MEAN X AND Y POSITIONS */
  /* -------------------------------- */
  w_object->xpos = SumXI/SumI;
  w_object->ypos = SumYI/SumI;


  return TRUE;
}




/*
  ---------------------------------------------------------------------
  ___  _     _        _     ___ _         _   ___          _   
 / _ \| |__ (_)___ __| |_  | __(_)_ _  __| | | _ \___ __ _| |__
| (_) | '_ \| / -_) _|  _| | _|| | ' \/ _` | |  _/ -_) _` | / /
 \___/|_.__// \___\__|\__| |_| |_|_||_\__,_| |_| \___\__,_|_\_\
          |__/                                                 
*/

/**
 * This is strongly based on the connect pixel finder, maininly for symmetry and ease of maintenance reasons.
 * The function being performed here does not absolutely need to use the recursion, but it is a relatively
 * neat way to do it.
 * 
 * This function runs basically the same was as the connected pixel finder except that it looks for
 * pixels greater then the local peak rather than just any greater than a threshold. It them climbs to
 * the peak of each source. Essentially it is same as the connected pixel finder but with an adaptive
 * threshold which is reset to the local max on ever step.
 * 
 * Other differences are 
 *	we do not subtract off the sky background.
 * 	we do not create the linked list of pixels because we would only have to free it again
 *	we do not mask pixels to 0.0 once they have been found.
 * 
 * This function returns its reults back in w_object, but the meansings are slightly different.
 * w_object->xpos	Integer X coord of brightest pixel rather than a true centroid.
 * w_object->ypos	Integer Y coord of brightest pixel rather than a true centroid.
 * w_object->peak	Counts in peak pixel. Not sky subtracted.
 * w_object->numpix	Number of steps taken in ascendng to the peak. Not the total number in the object.
 */

static int Object_Find_Peak(int naxis1,int naxis2,int x,int y, float *image,Object *w_object)
{
  
  /* ------------- */
  /* SET VARIABLES */
  /* ------------- */
  int x1,y1,cx,cy;     
  struct Point_Struct *point_list = NULL;
  struct Point_Struct *last_point = NULL;
  int point_count=0;
  float curr_peak,new_peak;


  /* ------------------------------- */
  /* ADD FIRST POINT TO BE PROCESSED */
  /* ------------------------------- */

/* #if LOGGING > 9 */
/*   Object_Log_Format("object","object.c","Object_Find_Peak",LOG_VERBOSITY_VERY_VERBOSE,NULL, */
/* 		    "adding point %d,%d to list.",x,y); */
/* #endif */



#if LOGGING > 7
  Object_Log_Format("object","object.c","Object_Find_Peak",LOG_VERBOSITY_INTERMEDIATE,NULL,
		    "(AGD) adding first point %d,%d,%.2f to peak-finding list (count was %d)",
		    x,y,image[(y*naxis1)+x],(w_object->numpix));
#endif


  if(!Point_List_Add(&point_list,&point_count,&last_point,x,y))
    return FALSE;
  

  /* -------------------------------- */
  /* RUN THROUGH POINTS ON POINT LIST */
  /* -------------------------------- */
#if LOGGING > 7
  Object_Log_Format("object","object.c","Object_Find_Peak",LOG_VERBOSITY_INTERMEDIATE,NULL,
		    "(AGD) running through points on point list");
#endif


  while(point_count > 0){
    

    /* start of per pixel stuff */
    /* ------------------------ */
    cx = point_list->x;
    cy = point_list->y;
    

    /* if pixel value above current peak */
    /* --------------------------------- */
    if (image[(cy*naxis1)+cx] > w_object->peak){   
  

#if LOGGING > 9
      Object_Log_Format("object","object.c","Object_Find_Peak",LOG_VERBOSITY_VERY_VERBOSE,NULL,
			"adding point %d,%d to object.",cx,cy);
#endif



      curr_peak = w_object->peak;

      w_object->peak = image[(cy*naxis1)+cx] ;  
      w_object->xpos = cx ;
      w_object->ypos = cy ;
      w_object->numpix ++;

      new_peak = w_object->peak;




#if LOGGING > 7
      Object_Log_Format("object","object.c","Object_Find_Peak",LOG_VERBOSITY_INTERMEDIATE,NULL,
			"(AGD) point %d,%d,%.2f > current peak %.2f, new peak %.2f, starting recursive loop (count %d)",
			cx,cy,image[(cy*naxis1)+cx],curr_peak,new_peak,(w_object->numpix));
#endif






      /* start of recursive stuff */
      for (x1 = cx-1; x1<=cx+1; x1++){
	for (y1 = cy-1; y1<=cy+1; y1++){
	  if (x1 >= naxis1 || y1 >= naxis2 || x1<0 || y1<0)  
	    continue;                                           /* set a flag here to say crap object? */

	  if (image[(y1*naxis1)+x1] > w_object->peak){


	    /* add this point to be processed */
#if LOGGING > 9
	    Object_Log_Format("object","object.c","Object_Find_Peak",LOG_VERBOSITY_VERY_VERBOSE,NULL,
			      "adding point %d,%d to list.",x1,y1);
#endif

	    
	    

#if LOGGING > 7
	    /* print first 5 then every 10000th point */
	    if ((w_object->numpix) <= 5)
	      Object_Log_Format("object","object.c","Object_Find_Peak",LOG_VERBOSITY_INTERMEDIATE,NULL,
				"(AGD) peak first 5: adding point (%d,%d,%.2f) to list (count %d)",
				x1,y1,image[(y1*naxis1)+x1],(w_object->numpix));


	    if (((w_object->numpix) % 10000 ) == 0) 
	      Object_Log_Format("object","object.c","Object_Find_Peak",LOG_VERBOSITY_INTERMEDIATE,NULL,
				"(AGD) peak runaway: adding point (%d,%d,%.2f) to list (count %d)",
				x1,y1,image[(y1*naxis1)+x1],(w_object->numpix));	      
#endif

	    if(!Point_List_Add(&point_list,&point_count,&last_point,x1,y1))
	      return FALSE;



	  }
	}/* end for on y1 */
      }/* end for on x1 */
    } /* end of if pixel value above current peak */


    /* if already added to object */
    /* -------------------------- */
    else {


#if LOGGING > 9
      Object_Log_Format("object","object.c","Object_Find_Peak",LOG_VERBOSITY_VERY_VERBOSE,NULL,
			"pixel %d,%d already added to object, ignoring.",cx,cy);
#endif
    }
    
    /* delete processed point */
    /* ---------------------- */

#if LOGGING > 9
    Object_Log_Format("object","object.c","Object_Find_Peak",LOG_VERBOSITY_VERY_VERBOSE,NULL,
		      "deleting point %d,%d from list.",cx,cy);
#endif
    if(!Point_List_Remove_Head(&point_list,&point_count))
      return FALSE;

  }/* end while on point list */
 
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
  Object_Log_Format("object","object.c","Object_Free",LOG_VERBOSITY_VERY_VERBOSE,NULL,"w_object (%p).",(*w_object));
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
      (*last_point)->next_point = new_point;          /* <---- new pixel added here */

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






  /* Log size of object as it's being created here, by noting every nth pixel */
  /* Don't want log level as high as that which lists every pixel but still higher than 
     normal level which I think is currently set to 5 */

/* #if LOGGING > 7 */
/*   if ( ((*point_count) % 10000 ) == 0 ) */
/*     Object_Log_Format("object","object.c","Point_List_Add",LOG_VERBOSITY_VERY_VERBOSE,NULL, */
/* 		      "Adding point %d %d to latest object, size now %ld\n",x,y,*point_count); */
/* #endif */







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
 * @see #Stellar_Ellipticity_Limit
 */
static void Object_Calculate_FWHM(Object *w_object,float BGmedian,int *is_stellar,float *fwhm)
{

  /* ---------------- */
  /* VARIABLES CONFIG */
  /* ---------------- */

  /* object-specific constants */
  /* ------------------------- */
  float object_xpos,object_ypos;  /* w_object->xpos,ypos */
  float object_peak;              /* peak pixel in object */
  float sex_onefifthpeak;         /* inner threshold to only deal with pixels brighter than
				     one-fifth of the peak pixel */

  /* object ellipticity */
  /* ------------------ */
  float x2I=0,y2I=0,xy2I=0;       /* running total of offset^2 x intensity in x, y & xy */
  float x2nd=0,y2nd=0,xy2nd=0;    /* first and second order moments of the ellipse data */
  float xoff=0,yoff=0;            /* offsets from the centre of the object determined by */
                                  /*   the centroid curpix->x = spatial pixel values */
  float intensity=0;              /* intensity in each pixel corrected for background */
  float SumI=0;                   /* scalar sum of pixel intensities */
  float aux=0,aux2=0;             /* auxiliary variable */
  float minor=0,major=0;          /* semi-minor and semi-major axes of the object ellipse */
  float ellip;                    /* ellipticity = (major-minor)/major */
  float ellip_theta;		  /* Orientation in the frame of the ellipticity major axis */

  HighPixel *curpix;              /* pixel pointer */
  char stellarflag[32];           /* stellar flag string for diagnostics */


  /* SExtractor code & 1/5th peak for fwhm calc */
  /* ------------------------------------------ */
  float sex_s,sex_sx,sex_sy,sex_sxx,sex_sxy;
  float sex_pix;                  /* pixel value */
  float sex_dx,sex_dy;            /* pixel x,y offsets from centroid */
  float sex_lpix;                 /* natural log of pixel value */
  float sex_inverr2;              /* (pixel value)^2 - why's it called inverr2? */
  float sex_d2;                   /* pixel radial distance from centroid (squared) */
  float sex_d;
  float sex_b;
  float sex_fwhm;




  /* ------------------------- */
  /* OBJECT-SPECIFIC CONSTANTS */
  /* ------------------------- */
  object_xpos = w_object->xpos;
  object_ypos = w_object->ypos;
  object_peak = w_object->peak;
  sex_onefifthpeak = object_peak / 5.0; 



  /* --------------------------------------------- */
  /* CALCULATE OBJECT ELLIPTICITY (VIA 2ND MOMENT) */
  /*                                               */
  /* uses equations from SEXtractor and PISA       */
  /* manual - see SUN109.1                         */
  /* --------------------------------------------- */

  /*
    2nd moment
    ----------
  */
  curpix = w_object->highpixel;
  while(curpix !=NULL)
    {
      intensity=(curpix->value);
      xoff = object_xpos - (curpix->x);
      yoff = object_ypos - (curpix->y);
      x2I  += xoff*xoff*intensity;
      y2I  += yoff*yoff*intensity;
      xy2I += xoff*yoff*intensity;
      SumI += intensity;
      curpix=curpix->next_pixel;
    }
  x2nd = x2I/SumI;
  y2nd = y2I/SumI;
  xy2nd = xy2I/SumI;


  /*
    ellipse axes
    ------------
  */
  aux = 0;
  aux = 2*(x2nd+y2nd);
  aux2 = 2*(sqrt( ((x2nd-y2nd)*(x2nd-y2nd)) + 4*xy2nd*xy2nd ));
  major = sqrt(aux+aux2); /* FULL ellipse major axis - not semimajor */
  minor = sqrt(aux-aux2); /* FULL ellipse minor axis - not semiminor */


  /*
    ellipticity
    -----------
  */
  ellip = (major-minor)/major;
  w_object->ellipticity = ellip;

  /* Measure orientation of the ellipticity axis.
   * I am not sure how we will want to define orientaion ultimately but for now
   * I use the standard cartesian convention of measuring anticlockwise from the X
   * axis so that a horizontal elongation is 0 (or 180) and a vertical elongation is 90 up (or -90 down) */
  if ( x2I == y2I ) 
    /* (M_PI / 2 == 90deg)  sets orientation 'vertical' if undefined for a _perfectly_ circular source */
    ellip_theta = 1.5707963268 ; 
  else {
    ellip_theta = 0.5 * atan2( (2.0*xy2I) , (x2I-y2I) ); 
    /* Keep to 0 < ellip_theta < 180deg */
    if (ellip_theta < 0) ellip_theta += 3.14159265359;
  }
  w_object->ellip_theta = ellip_theta;

#if LOGGING > 5
  Object_Log_Format("object","object.c","Object_Calculate_FWHM",LOG_VERBOSITY_VERY_VERBOSE,NULL,
		    "(%d) a = %.2f, b = %.2f\tellip = %.2f\ttheta = %.2f",w_object->objnum,major,minor,ellip,ellip_theta);
#endif
 

  /*
    set stellar flag
    ----------------
  */
  if (ellip <= Stellar_Ellipticity_Limit){
    (*is_stellar) = TRUE;            /* object is STELLAR */
    w_object->is_stellar = TRUE;
    sprintf(stellarflag,"stellar");
  }
  else {
    (*is_stellar) = FALSE;           /* object is NONSTELLAR */
    w_object->is_stellar = FALSE;
    sprintf(stellarflag,"NON-stellar");
  }

#if LOGGING > 5
  Object_Log_Format("object","object.c","Object_Calculate_FWHM",LOG_VERBOSITY_VERBOSE,NULL,
		    "object (%d) stellarflag %s, is_stellar [%d]",
		    w_object->objnum,stellarflag,w_object->is_stellar);
#endif





  /* -------------------------------- */
  /* CALCULATE FWHM IF OBJECT STELLAR */
  /* -------------------------------- */

  if (w_object->is_stellar == TRUE){ 

    /* 
       initialise
       ----------
    */
    sex_s = sex_sx = sex_sy = sex_sxx = sex_sxy = 0.0;


    /* 
       for each pixel IN OBJECT i.e. AFTER thresholding
       ------------------------------------------------
    */      
    curpix = w_object->highpixel;            /* set current pixel to object's first pixel */
    while(curpix != NULL){                   /* start looping through pixels in object */

      /* pixel value */
      sex_pix = curpix->value;                   /* This should have had BGmedian subtracted already */


      /* reject if sex_pix < 1/5 of peak*/
      if (sex_pix < sex_onefifthpeak){
	curpix=curpix->next_pixel;               /* next pixel in object */
	continue;                                /* go to next iteration of while loop */
      }

      /* X,Y offsets from centroid */
      sex_dx = curpix->x - object_xpos;
      sex_dy = curpix->y - object_ypos;

      /* natural log */
      sex_lpix = log(sex_pix);
      sex_inverr2 = sex_pix*sex_pix;
      sex_s += sex_inverr2;

      /* find radial distance of pixel from centroid */
      sex_d2 = sex_dx*sex_dx+sex_dy*sex_dy;             /* d^2 = dx^2 + dy^2 */
      sex_sx += sex_d2*sex_inverr2;           
      sex_sxx += sex_d2*sex_d2*sex_inverr2;
      sex_sy += sex_lpix*sex_inverr2;
      sex_sxy += sex_lpix*sex_d2*sex_inverr2;

      curpix=curpix->next_pixel;                        /* next pixel in object */
    }

/* RJS making first attempt at ellipse orientation */
/*       if ( sex_sxx == sex_syy ) 
        theta = PI/4.0;
      else
        theta = 0.5 * atan2( (2.0*sex_sxy) , (sex_sxx-sex_syy) ); */
/* End RJS */

    sex_d = sex_s*sex_sxx-sex_sx*sex_sx;
    if (fabs(sex_d) > 0.0) {
      sex_b = -(sex_s*sex_sxy-sex_sx*sex_sy)/sex_d;
      sex_fwhm = (float)(1.6651/sqrt(sex_b));
      if (sex_fwhm > 0.0) {
	sex_fwhm -= 1.0/(4.0*sex_fwhm);
	w_object->fwhmx = sex_fwhm;
	w_object->fwhmy = sex_fwhm;
	(*fwhm) = sex_fwhm;
      }
      else {
	w_object->fwhmx = DEFAULT_SEEING_ZERO;
	w_object->fwhmy = DEFAULT_SEEING_ZERO;
	(*fwhm) = DEFAULT_SEEING_ZERO;	  
      }
    }
    else {
      w_object->fwhmx = DEFAULT_SEEING_SEXD_ZERO;
      w_object->fwhmy = DEFAULT_SEEING_SEXD_ZERO;
      (*fwhm) = DEFAULT_SEEING_SEXD_ZERO;
    }



   
  } /* end of CALCULATE FWHM IF OBJECT IS STELLAR */



  
  /* ------------------------------------- */
  /* SET DEFAULT FWHM IF OBJECT NONSTELLAR */
  /* ------------------------------------- */

  else {
    w_object->fwhmx = DEFAULT_SEEING_NONSTELLAR;
    w_object->fwhmy = DEFAULT_SEEING_NONSTELLAR;
    (*fwhm) = DEFAULT_SEEING_NONSTELLAR;

#if LOGGING > 5
    Object_Log_Format("object","object.c","Object_Calculate_FWHM",LOG_VERBOSITY_VERBOSE,NULL,
		      "object (%d) is %s, setting FWHM to %f",
		      w_object->objnum,stellarflag,DEFAULT_SEEING_NONSTELLAR);
#endif
    
  }
  
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
** Revision 1.16  2014/07/30 21:38:06  eng
** *** empty log message ***
**
** Revision 1.15  2014/07/30 18:26:21  eng
** Added calculation of ellip_theta, the long axis orientaion of ellipticity in the
** frame. The value is not yet passed back to the calling application. That comes
** in the next update.
**
** Revision 1.14  2014/07/30 17:42:09  eng
** Continuation of changes from v1.13. Finished adding teh function that allows you to set
** the saturation level and created a DEFAULT #def to set the satuation to 63000 ADU if not
** set explicitly by the calling application.
**
** When determining the list of candidate objects which will be used to derive the
** median seeing, added a check to reject any sources with peak > saturation limit. Only
** unsaturated stars to be used to measure seeing.
**
** Revision 1.13  2014/07/30 17:26:41  eng
** Added Object_Saturation_Limit_Set() to allow the calling application to define
** a saturation limit for teh detector that overrides the built in default.
**
** Revision 1.12  2009/08/19 17:54:00  eng
** Tweaked to apply the MARGIN criteria instead to the object centroid position
** (xpos,ypos) rather than pixels used in object creation. This is because
** there's nothing to stop the finding-of-connected-pixels process going into
** the margin area in some cases.
**
** Revision 1.11  2009/08/11 14:24:32  cjm
** Turned down per-pixel logging.
** Reenabled rcsid.
**
** Revision 1.10  2009/08/06 13:34:03  eng
** Defining a "margin" around the image frame to prevent objects being made from
** pixels within a distance N pixels from the frame edge. This is primarily to
** avoid objects made from the "hot rows" at the bottom of RATCam images. NB: I
** would have thought however that such an object would have been rejected as
** very much non-stellar anyway later in the code, but so be it. (JMM)
**
** Revision 1.9  2009/06/12 10:32:06  cjm
** Added Object_Stellar_Ellipticity_Limit_Set/Stellar_Ellipticity_Limit to make
** the stellar ellipticity limit configurable.
**
** Revision 1.8  2009/01/30 15:20:31  cjm
** Replaced log bit levels with log_udp.h verbositys.
**
** Revision 1.7  2009/01/28 14:18:53  cjm
** Added extra parameters to logging routines for GLS support.
**
** Revision 1.6  2008/10/07 13:55:32  cjm
** Fixed object total counts measure by uncommenting increment.
**
** Revision 1.5  2008/10/07 13:17:59  cjm
** Changed some pixel logging from OBJECT to PIXEL bit.
**
** Revision 1.4  2008/10/07 12:57:32  eng
** Added HighPixel *curpix; (jmm)
**
** Revision 1.3  2008/10/01 15:53:36  eng
** Some invisible non-printing character was preventing compilation with error
** syntax error: line 1810, token = " ", char = '270'
** I fixed it by deleteing a couple of apparently blank lines.
**
** Revision 1.2  2008/09/30 19:42:44  eng
** changed object_jmm.h to object.h
**
** Revision 1.1  2008/09/30 19:26:27  eng
** Initial revision
**
** Revision 1.12.2.10  2008/09/26 12:51:05  eng
** Minor tweaks:
** * altered some default seeing flag names
** * in Object_Calculate_FWHM:
**   - removed moffat & 2M methods
**   - removed theta & diff in ellipticity calculation
**   - moved indirections to object-specific constants
**
** Revision 1.12.2.9  2008/09/25 16:18:58  eng
** Object_List_Get_New has been renamed back to Object_List_Get. Also, as we now only
** use one blanket threshold, "thresh1" has been renamed back to "thresh". The extra args
** to Object_List_Get, namely thresh2 and the object counters, have been removed,
** putting it back in its original state. The object counters themselves are now local
** variables inside Object_List_Get and are not passed out of the function. Also tidied
** up some of the comments.
**
** Revision 1.12.2.8  2008/09/25 15:30:36  eng
** Added RJS' Object_Find_Peak to find the local peak of each individual object and
** then set a local-to-object 1/5th peak threshold (replacing the blanket "thresh2"
** of earlier versions).
**
** Revision 1.12.2.7  2008/09/25 13:27:30  eng
** Fix of memory leak whereby the call to Object_Calculate_FWHM only produced sane
** results if a mystery 'printf' was added just before or after the function call.
** The variable 'xy2I' in Object_Calculate_FWHM was uninitialised and the best
** explanation is that the extra printf happened to push the pointer into the right
** place. It is initialised now and everything seems to run now.
**
** Revision 1.12.2.6  2008/09/15 13:12:09  eng
** Changed fwhm_lt_dia_count --> usable_count
** Deleted Object_List_Get
** Added new default FWHM settings in case of bad data
** Added ellipticity to sizefwhm struct, for diagnostic purposes.
** In Object_List_Get_New:
**   - Removed fwhmlist
**   - Added stellarity & +ve fwhm to selection for fwhmarray:
**       if ((obj_fwhm < obj_dia) && (w_object->is_stellar == 1) && (obj_fwhm > 0.0))
** In Object_Calculate_FWHM (sExtractor-adapted-code section):
**   - Added default fwhm settings in case fwhm -ve or sex_d -ve.
**
** Revision 1.12.2.5  2008/09/09 11:08:10  eng
** This is basically the same two-threshold code as before but with a few extra tweaks
** based on/adapted from SExtractor code (not the manual). This uses only the top 80%
** of pixels, i.e. ignores the "wings", and then uses a factor to compensate.
**
** Revision 1.12.2.4  2008/08/30 15:37:21  eng
**   Closing 1.12.2.3 - FWHM code was fixed, part typo error. New code from scratch worked but
**   produced FWHMs that were now much larger than the original.
**
** Revision 1.12.2.3  2008/08/22 15:58:17  eng
** Working version of two-threshold ("TT") algorithm. Emphasis on applying the TT to the old
** 2nd-moment ("2M")method of calculating FWHM again. Why? Because the Moffatt-fitting algorithm
** barfs out during a focus run, which admittedly may be a coding problem rather than
** fundamental method problem, AND it is very slow, which again could just be down to
** breadboard (not efficinet) coding practices. Looked again at the 2m version because
** it looked more robust and a lot faster. However this version still produces wrong
** FWHMs - a fault which might be down to the maths used for the 1st & 2nd moment calc.
** Therefore, going to try coding the 1st/2nd moments again from scratch rather than try to
** second-guess the original existing code. 1st moment will still be calculated in
** Object_List_Get_New & 2nd in Object_Calculate_FWHM, but the maths involved will be sourced
** from fresh locations.
**
** Revision 1.12.2.2  2008/06/05 13:56:51  eng
** *** empty log message ***
**
** Revision 1.12.2.1  2008/06/03 14:25:48  eng
** Another branch-off from 1.12, this time to persue the concept of just
** using two thresholds; one high-sigma to use when detecting objects,
** the other low-sigma to use when getting connected pixels. This way we
** get well-sampled objects out to the "wings", but only for bright objects, not
** wasting time doing it for every object a low-sigma value above the median.
**
** Revision 1.12  2008/05/29 22:26:09  eng
** Checked in by RJS to enable JMM to roll back on Monday any changes he disapproves of.
** THis is just after Jon has reverted teh calling parameters to the old 9 parameter format.
**
** Revision 1.11  2008/05/29 15:35:48  eng
** Having to check this version in unfinished because I need to look at an older
** version for RJS' benefit. If I just had a different file with a different name
** then this would be easy, but because the earlier version is buried in the depths
** of RCS with no way of even looking at it, I have to actually check out the older
** version just to look at it. This will apparently destroy the current version unless
** it's checked in, so to keep it safe, I have to temporarily check it in even though
** it's unfinished. RCS: useful in some ways, bloody annoying in others.
**
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
