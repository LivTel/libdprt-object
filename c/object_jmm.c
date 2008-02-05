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
** $Header: /space/home/eng/cjm/cvs/libdprt-object/c/object_jmm.c,v 1.7 2008-02-05 18:30:18 eng Exp $
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
 * @version $Revision: 1.7 $
 */



/*
  $Log: not supported by cvs2svn $
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
#include <gsl/gsl_sf.h>


/* ------------------------------------------------------- */
/* hash definitions */
/* ------------------------------------------------------- */
/**
 * Value to use for seeing when something goes wrong (in pixels).
 * Should be larger than 10 pixels, so RCS thinks seeing is "bad".
 * RJS wants really large value so archive searches can differentiate between real bad seeing and failed reductions.
 */
#define DEFAULT_BAD_SEEING    (999.0)

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
/*static char rcsid[] = "$Id: object_jmm.c,v 1.7 2008-02-05 18:30:18 eng Exp $";*/
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
static int Sort_Float(const void *data1,const void *data2);
double moffat(double x, double k, double a, double b);
double delta(const double *x, const double *y, const int items, const double parameters[]);
int sign(double x);
double findMax(const double *a, const int items);
double optimize(const double *x, const double *y, int items, double params[]);
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
		    Object **first_object,int *sflag,float *seeing)
{
  Object *w_object = NULL;
  Object *last_object = NULL;
  Object *next_object = NULL;
  HighPixel *curpix = NULL;
  float fwhm;
  float *fwhm_list = NULL;
  int y,x,count,done,is_stellar,mid;
  int fwhm_count = 0;

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
  count = 0;
  for(y=0;y<naxis2;y++)
    {
      for(x=0;x<naxis1;x++)
	{
#if LOGGING > 9
	  Object_Log_Format(OBJECT_LOG_BIT_PIXEL,"Object_List_Get:searching pixel %d,%d.",x,y);
#endif
	  if(image[(y*naxis1)+x] > thresh)
	    {
	      count++;
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
	      w_object->objnum = count;
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
  Object_Log_Format(OBJECT_LOG_BIT_GENERAL,"Object_List_Get:Found %d objects.",count);
#endif
  /* have we got any objects in the frame? */
  if(count == 0)
    {
      (*seeing) = DEFAULT_BAD_SEEING;
      (*sflag) = 1; /* the seeing was fudged. */
      (*first_object) = NULL;
      Object_Error_Number = 6;
      sprintf(Object_Error_String,"Object_List_Get:No objects found.");
      Object_Warning();
      /* We used to return FALSE (error) here.
      ** But there are cases where it is OK to have no objects - e.g. Moon images.
      ** We want  a fake seeing to be written to the FITS headers,
      ** So we generate a warning message and return TRUE.
      ** Note this means any program using Object_List_Get must be able to cope with
      ** a NULL object list.
      */
      return TRUE;
    }
  /* go through list of objects, get rid of any with less than npix */
  /* diddly why got get rid of any with lots of pixels i.e. non-stellar stuff? */
  /* find new first_object with numpix > npix */
#if LOGGING > 0
  Object_Log(OBJECT_LOG_BIT_GENERAL,"Object_List_Get:Finding useful objects.");
#endif
  w_object = (*first_object);
  count = 0;
  done = FALSE;
  while(done == FALSE)
    {
      if(w_object->numpix >= npix)
	done = TRUE;
      else
	{
#if LOGGING > 5
	  Object_Log_Format(OBJECT_LOG_BIT_OBJECT,"Object_List_Get:deleting object(1) at %.2f,%.2f(%d).",
			    w_object->xpos,w_object->ypos,w_object->numpix);
#endif
	  /* take copy of next object pointer */
	  next_object = w_object->nextobject;
	  /* delete w_object */
	  Object_Free(&w_object);
	  /* set w_object to next object */
	  w_object = next_object;
	  /* if we've reached the end of the list, bail out */
	  if(w_object == NULL)
	    done = TRUE;
	}
    }
  if(w_object == NULL)
    {
      (*seeing) = DEFAULT_BAD_SEEING;
      (*sflag) = 1; /* the seeing was fudged. */
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
  count = 1;
  done = FALSE;
  /* go through rest of object list, deleting objects with numpix < npix */
  while(w_object != NULL)
    {
      /* take copy of next object to go to */
      next_object=w_object->nextobject;
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
	  count++;
	  /* tell last object this is it's next object */
	  last_object->nextobject=w_object;
	  /* set objects number */
	  w_object->objnum=count;
	  /* set the last object in the list to be this object */
	  last_object = w_object;
#if LOGGING > 5
	  Object_Log_Format(OBJECT_LOG_BIT_OBJECT,"Object_List_Get:object %d at %.2f,%.2f(%d) is ok(2).",
			    w_object->objnum,w_object->xpos,w_object->ypos,w_object->numpix);
#endif
	}
      /* change to next object */
      w_object = next_object;
    }
  last_object->nextobject=NULL;
  /* extra debug - list connected pixels in all objects */
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
	/* goto next pixel */
	 curpix = curpix->next_pixel;
      }
    /* goto next object */
      w_object = w_object->nextobject;		
  }
#endif

  /* find fwhm */
#if LOGGING > 0
  Object_Log(OBJECT_LOG_BIT_GENERAL,"Object_List_Get:Finding FWHM of objects.");
#endif
  w_object = (*first_object);
#if LOGGING > 10
  Object_Log_Format(OBJECT_LOG_BIT_OBJECT,"Object_List_Get:w_object (%p) set from first_object (%p).",
		    w_object,(*first_object));
#endif
  while(w_object != NULL)
    {
#if LOGGING > 5
      Object_Log_Format(OBJECT_LOG_BIT_FWHM,"Object_List_Get:Calculating FWHM for object at %.2f,%.2f.",
			w_object->xpos,w_object->ypos);
#endif
      Object_Calculate_FWHM(w_object,image_median,&is_stellar,&fwhm);
#if LOGGING > 5
      Object_Log_Format(OBJECT_LOG_BIT_FWHM,"Object_List_Get:"
			"object at %.2f,%.2f has FWHM %.2f pixels and is_stellar = %d.",
			w_object->xpos,w_object->ypos,fwhm,is_stellar);
#endif
      if(is_stellar)
	{
	  if(fwhm_list == NULL)
	    fwhm_list = (float*)malloc(sizeof(float));
	  else
	    fwhm_list = (float*)realloc(fwhm_list,(fwhm_count+1)*sizeof(float));
#ifdef MEMORYCHECK
	  if(fwhm_list == NULL)
	    {
	      Object_Error_Number = 8;
	      sprintf(Object_Error_String,"Object_List_Get:Failed to alllocate FWHM list(%d).",
		      fwhm_count);
	      return FALSE;
	    }
#endif
	  fwhm_list[fwhm_count] = fwhm;
	  fwhm_count++;
	}
      w_object = w_object->nextobject;		
    }
#if LOGGING > 0
  Object_Log(OBJECT_LOG_BIT_GENERAL,"Object_List_Get:Calculating final seeing.");
#endif
  /* diddly In the original code theres a bit here concerned with numobjs
  ** that makes no sense at all to me */
  if(fwhm_count > 0)
    {
      qsort(fwhm_list,fwhm_count,sizeof(float),Sort_Float);
      mid = (fwhm_count-1)/2;
      (*seeing) = fwhm_list[mid];
      (*sflag) = 0;
      /* If the seeing is less than 0.01 set the seeing to DEFAULT_BAD_SEEING (pixels) 
      ** and sflag to show the seeing was fudged */
      if ((*seeing)<=0.01)
	{
	  (*seeing) = DEFAULT_BAD_SEEING;
	  (*sflag) = 1;
	}
    }
  else
    {
      (*seeing) = DEFAULT_BAD_SEEING;
      (*sflag) = 1; /* the seeing was fudged. */
    }
  if(fwhm_list != NULL)
    free(fwhm_list);
#if LOGGING > 0
  Object_Log_Format(OBJECT_LOG_BIT_GENERAL,"Object_List_Get:number of objects > %d pixels = %d",npix,count);
  Object_Log_Format(OBJECT_LOG_BIT_GENERAL,"Object_List_Get:number of objects identified as stellar = %d",
		    fwhm_count);
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




/**
 * Routine to return the object error number.
 * @return The object error number.
 * @see #Object_Error_Number
 */
int Object_Get_Error_Number(void)
{
  return Object_Error_Number;
}


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

/**
 * Routine to set the Log_Data.Log_Filter_Level.
 * @see #Log_Data
 */
void Object_Set_Log_Filter_Level(int level)
{
  Log_Data.Log_Filter_Level = level;
}

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

/* ------------------------------------------------------- */
/* internal functions */
/* ------------------------------------------------------- */
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
static int Object_List_Get_Connected_Pixels(int naxis1,int naxis2,float image_median,int x,int y,float thresh,
					    float *image,Object *w_object)
{
  int x1,y1,cx,cy;
  HighPixel *temp_hp=NULL;
  HighPixel *curpix = NULL;
  struct Point_Struct *point_list = NULL;
  struct Point_Struct *last_point = NULL;
  int point_count=0;

  /* add first point to be processed */
#if LOGGING > 9
  Object_Log_Format(OBJECT_LOG_BIT_POINT,"Object_List_Get_Connected_Pixels:"
		    "adding point %d,%d to list.",x,y);
#endif
  if(!Point_List_Add(&point_list,&point_count,&last_point,x,y))
    return FALSE;
  /* while there are points to process, process them */
  while(point_count > 0)
    {
      /* start of per pixel stuff */
      cx = point_list->x;
      cy = point_list->y;
      /* check this point hasn't been added to the object since it was added. */
      if (image[(cy*naxis1)+cx] > thresh)
	{
#if LOGGING > 9
	  Object_Log_Format(OBJECT_LOG_BIT_PIXEL,"Object_List_Get_Connected_Pixels:"
			    "adding pixel %d,%d to object.",cx,cy);
#endif
	  /* allocate new object pixel */
	  temp_hp=(HighPixel*)malloc(sizeof(HighPixel));
#ifdef MEMORYCHECK
	  if(temp_hp == NULL)
	    {
	      Object_Error_Number = 3;
	      sprintf(Object_Error_String,"Object_List_Get_Connected_Pixels:"
		      "Failed to allocate temp_hp.");
	      return FALSE;
	    }
#endif
	  temp_hp->next_pixel=NULL;
	  if(w_object->highpixel == NULL)
	    {
	      w_object->highpixel = temp_hp;
	      w_object->last_hp = temp_hp;
	    }
	  else
	    {
	      w_object->last_hp->next_pixel = temp_hp;
	      w_object->last_hp = temp_hp;
	    }
	  /* set currentx, current y from point list element. */
	  w_object->last_hp->x=cx;
	  w_object->last_hp->y=cy;
	  w_object->last_hp->value=image[(cy*naxis1)+cx] - image_median;
	  /* important for recursion - stops infinite loops. */
	  image[(cy*naxis1)+cx]=0.0;
	  /* end of per-pixel stuff*/
	  /* do some overall w_object stats.
	  ** Copied from getObjectList, again per pixel stuff. */
	  curpix = w_object->last_hp;
	  w_object->total=(w_object->total)+(temp_hp->value);
	  w_object->xpos=(w_object->xpos)+((temp_hp->x)*(temp_hp->value)); 
	  w_object->ypos=(w_object->ypos)+((temp_hp->y)*(temp_hp->value));
	  if ((w_object->peak)<(temp_hp->value))
	    w_object->peak=(temp_hp->value);  
	  w_object->numpix ++;
	  /* start of recursive stuff */
	  for (x1 = cx-1; x1<=cx+1; x1++)
	    {
	      for (y1 = cy-1; y1<=cy+1; y1++)
		{
		  if (x1 >= naxis1 || y1 >= naxis2 || x1<0 || y1<0)  
		    continue; /*set a flag here to say crap object?*/
		  if (image[(y1*naxis1)+x1] > thresh)
		    {
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
      else
	{
#if LOGGING > 9
	  Object_Log_Format(OBJECT_LOG_BIT_PIXEL,"Object_List_Get_Connected_Pixels:"
			    "pixel %d,%d already added to object, ignoring.",cx,cy);
#endif
	}
      /* delete processed point */
#if LOGGING > 9
      Object_Log_Format(OBJECT_LOG_BIT_POINT,"Object_List_Get_Connected_Pixels:"
			"deleting point %d,%d from list.",cx,cy);
#endif
      if(!Point_List_Remove_Head(&point_list,&point_count))
	return FALSE;
    }/* end while on point list */
  /* calculate mean x and y positions.
  ** Copied from getObjectList */
  w_object->xpos=(w_object->xpos)/(w_object->total);
  w_object->ypos=(w_object->ypos)/(w_object->total);
  return TRUE;
}

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
---------------------------------------------------------------------
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
  /* first and second order moments of the ellipse data */
  float x2nd=0,y2nd=0,xy2nd=0;
  /* off sets from the centre of the object determind by the centroid */
  /* curpix->x  = spatial pixel values */
  float xoff=0,yoff=0;
  /* intensity = intensity in each pixel corrected for background*/
  float intensity=0;
  /* auxillary variable */
  float aux=0,aux2=0;
  /* semi-minor and semi-major axes of the ellipse */
  float minor=0,major=0;
  /* inclination angle of the ellipse, measured from Naxis1 to semi-major */
  /* axis (ie anti clockwise) */
  float theta;
  /* difference between major and minor axes */
  float diff;
  /* pixel pointer */
  HighPixel *curpix;


  /* Moffat curve fitting optimisation */
  /* --------------------------------- */
  double *pixr, *pixz;                 /* optimisation r,z arrays */
  int ipix;                            /* pixel counter */
  double dx,dy;                        /* pixel offset from 1st moment */
  double params[3];                    /* Moffat curve parameters (k,a,b) */
  int m;                               /* optimisation infinite counter */  
  double k,a,b;                        /* Moffat curve parameters (k,a,b) */
  double FWHM;                         /* FWHM */



  /*
    ----------
    2ND MOMENT
    ----------
  */

  /* calculate the 2nd moments of the ellipse in x y and xy  */
  /* uses equations from SEXtractor manual and PISA manual -
     SUN109.1 */

  /* should really check w_object, fwhm and is_stellar here for NULL. */
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
  if (major!= minor)
    theta=atan2((2*xy2nd),(x2nd-y2nd));
  else
    theta=0;
  
  diff=major-minor;
  if (diff<=major*0.3)
    {
      (*is_stellar) = TRUE;            /* object is STELLAR */
      w_object->is_stellar = TRUE;
    }
  else
    {
      (*is_stellar) = FALSE;           /* object is NONSTELLAR */
      w_object->is_stellar = FALSE;
    }


  /*
    ---------------
    FWHM CONDITIONS
    ---------------
  */

  /*
    -----------
    non-stellar
    -----------
  */
  
  if (w_object->is_stellar == FALSE){
    w_object->fwhmx = 999.0;
    w_object->fwhmy = 999.0;
    (*fwhm) = 999.0;
  }

  /*
    ---------
    too faint
    ---------
  */
/*   else if ((w_object->peak/BGmedian) < 10.0){ */
/*     w_object->fwhmx = 999.0; */
/*     w_object->fwhmy = 999.0; */
/*     (*fwhm) = 999.0;     */
/*   } */


  /*
    ------------------------
    FWHM IF BRIGHT & STELLAR
    ------------------------
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





















/*    for(ii = 1; ii < w_object->numpix; ii++) { */
/*       if(p[ii].z > peak) */
/* 	peak = p[ii].z; */
/*     } */

/*     params[0] = peak; */
/*     params[1] = mf_a_start; */
/*     params[2] = mf_b_start; */
    


/*     /\* OPTIMIZE *\/ */
/*     for (i=0; i<=2000; i++){ */
/*       sm = 0.0;                         /\* sum of slopes *\/ */
/*       params2[0] = params[0];           /\* initialise 2nd parameter array *\/ */
/*       params2[1] = params[1]; */
/*       params2[2] = params[2]; */
/*       d1 = delta(p, w_object->numpix, params);    /\* calc goodness of fit (delta) for initial params *\/ */
 
/*       for(j = 1; j <= 2; j++) {                     /\* for parameters a and b (not Io) *\/ */
/* 	tmp = params[j];                            /\* current parameter value *\/ */
/* 	params2[j] = params[j] + EPS;               /\* increment new parameter by EPS *\/ */
/* 	d2 = delta(p, w_object->numpix, params2);   /\* calc new goodness of fit (delta) *\/ */
/* 	m = (d2 - d1) / EPS;                        /\* calc slope in delta space *\/ */
/* 	sm += (1 + m) * (1 + m);                    /\* add to sum of slopes *\/ */
	
/* 	/\* This is stolen from iRprop. *\/ */
/* 	if(momentum[j] * m > 0) {         /\* momentum & slope same signs? *\/ */
/* 	  momentum[j] *= 1.1;             /\* accelerate forwards *\/ */
/* 	} else { */
/* 	  momentum[j] *= -0.9091;         /\* slow down & backwards *\/ */
/* 	} */
/* 	/\* Threshold the momentum before applying *\/ */
/* 	momentum[j] = sign(momentum[j]) * MIN(fabs(momentum[j]),0.2); */

/* 	params[j] = params2[j] - momentum[j];  /\* change original parameter *\/ */
/* 	params2[j] = tmp;                      /\* new parameter becomes old one *\/ */


/* 	/\* STOPPING CONDITIONS *\/ */
	
/* 	/\* -- by delta *\/ */

/* 	if (d2 < minDelta) {              /\* if difference small enough *\/ */
/* 	  printf("minDelta stop at iteration %d\n", i); */
/* 	  minDelta = d2;                  /\* note values & carry on a bit.... *\/ */
/* 	  best[0] = params2[0]; */
/* 	  best[1] = params2[1]; */
/* 	  best[2] = params2[2]; */
/* 	  bestIter = i; */
/* 	} else { */
/* 	  if(i - bestIter > 20)           /\* ... but eventually stop after   *\/ */
/* 	    break;                        /\* 20 more goes if no new bestIter *\/ */
/* 	} */
	
	
/* 	/\* -- by slope *\/ */
/* 	if(sm < earlystop) { */
/* 	  printf("Early Stop at iteration %d\n", i); */
/* 	  break; */
/* 	} */

/*       } /\* end of parameter loop *\/ */
      
/*       /\* Write final parameters to params array *\/ */
/*       params[0] = best[0]; */
/*       params[1] = best[1]; */
/*       params[2] = best[2]; */

/*     } */

/*     /\* printf("-- Optimized moffat parameters for object %d: Io = %f\ta = %f\tb = %f\n", */
/*        w_object->objnum,params[0], params[1], params[2]); *\/ */

/*     moffat_a = params[1]; */
/*     moffat_b = params[2]; */


/*     /\*  */
/*        -------------- */
/*        CALCULATE FWHM   */
/*        -------------- */
/*     *\/ */
  
/*     /\* printf("-- calc FWHM\n"); *\/ */
/*     /\* Calculate FWHM from ma, mb *\/ */
/*     fw = 2.0 * moffat_a * sqrt( pow(2.0,(1.0/moffat_b)) - 1.0 ); */

/*     /\* printf("-- write to object struct\n"); *\/ */
/*     /\* Put in object struct requirements *\/ */
/*     w_object->fwhmx = fw; */
/*     w_object->fwhmy = fw; */
/*     (*fwhm) = fw; */


/*     /\*printf("%d\t%f\t%f\t%f\n",w_object->objnum,min_a,min_b,fw);*\/ */

/*     /\* printf("-- clear struct\n"); *\/ */
/*     /\* Clear struct 'p' memory *\/ */
/*     if (p != NULL) */
/*       free(p); */
/*   } /\* end of if stellar flag TRUE *\/ */

/* } */


/* Sort_Float */
static int Sort_Float(const void *data1,const void *data2)
{
  return *((float *) data1) < * ((float *) data2);
}

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


/* SIGN */
int sign(double x){
  if ( x == 0.0 )
    return 0;
  else
    return (int) (x/fabs(x)); 
}





/*
** $Log: not supported by cvs2svn $
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
