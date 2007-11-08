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
** $Header: /space/home/eng/cjm/cvs/libdprt-object/c/object_jmm.c,v 1.2 2007-11-08 17:02:03 jmm Exp $
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
 * @version $Revision: 1.2 $
 */



/*
  $Log: not supported by cvs2svn $

*/




/**
 * This hash define is needed before including source files give us POSIX.4/IEEE1003.1b-1993 prototypes.
 */
#define _POSIX_SOURCE 1
/**
 * This hash define is needed before including source files give us POSIX.4/IEEE1003.1b-1993 prototypes.
 */
#define _POSIX_C_SOURCE 199309L

#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include "object_jmm.h"

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
static char rcsid[] = "$Id: object_jmm.c,v 1.2 2007-11-08 17:02:03 jmm Exp $";
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
static void Object_Calculate_FWHM(Object *w_object,int *is_stellar,float *fwhm);
static void Object_Free(Object **w_object);
static int Point_List_Remove_Head(struct Point_Struct **point_list,int *point_count);
static int Point_List_Add(struct Point_Struct **point_list,int *point_count,struct Point_Struct **last_point,
			  int x,int y);
static int Sort_Float(const void *data1,const void *data2);

/* ------------------------------------------------------- */
/* external functions */
/* ------------------------------------------------------- */
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
      Object_Calculate_FWHM(w_object,&is_stellar,&fwhm);
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

/**
 * Routine to calculate the FWHM of the specified object.
 * @param w_object The object to calculate the FWHM from.
 * @param is_stellar The address of an integer to store a boolean. On exit of the routine,
 *        will be TRUE if stellar, FALSE if non-stellar.
 * @param fwhm An address to store the calculated full width half maximum, in pixels.
 */
static void Object_Calculate_FWHM(Object *w_object,int *is_stellar,float *fwhm)
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


  float xpos,ypos;                            /* barycentre coordinates */
  int int_xpos,int_ypos;                      /* nearest pixel to barycentre */
  float delta_x,delta_y;                      /* dist. btwn pixel and real (float) barycentre */
  float dummy_x,dummy_y,dummy_v;              /* dummy vars for sigma calc */
  float Sum_I_horiz_cut,Sum_I_vert_cut;       /* Sum of pixel values along transects */
  float Sigma_x,Sigma_y;
  float pix_I;                                /* pixel value */
  int pix_x, pix_y;                           /* pixel x position, y position */
  float pix_v;                                /* pixel value */


  /* TEST ONLY   */
  FILE *test_of;
  test_of = fopen("blob.txt", "w");
	



  /* calculate the 2nd moments of the ellipse in x y and xy  */
  /* uses equations from SEXtractor manual and PISA manual -
     SUN109.1 */

  /* should really check w_object, fwhm and is_stellar here for NULL. */
  curpix = w_object->highpixel;
  while(curpix !=NULL)
    {
      /* diddly this was:
      ** intensity=(curpix->value)-(fitsimage->median);
      ** However, Object_List_Get_Connected_Pixels already subtracts the median. 
      ** See below, however, for futher potential screw-up. */
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
  
  /* ********************** */
 
  diff=major-minor;
  if (diff<=major*0.3)
    {
      /* object is stellar; */
      (*is_stellar) = TRUE;
      w_object->is_stellar = TRUE;
    }/* endif */
  else
    {
      (*is_stellar) = FALSE;
      w_object->is_stellar = FALSE;
    }











	
  /*
    __        _          
    / _|_ __ _| |_  _ __  
    |  _\ V  V / ' \| '  \ 
    |_|  \_/\_/|_||_|_|_|_|
	  
  */


  /*
    ----------
    INITIALISE
    ----------
  */
  xpos = w_object->xpos;
  ypos = w_object->ypos;
  curpix = w_object->highpixel;
  delta_x = 0.0;                                  /* dist. btwn x and real barycentre x position (float) */
  delta_y = 0.0;                                  /* dist. btwn y and real barycentre y position (float) */
  dummy_x = 0;                                    
  dummy_y = 0;
  dummy_v = 0;
  Sum_I_horiz_cut = 0;                            /* Sum of pixel values through horizontal transect */
  Sum_I_vert_cut = 0;                             /* Sum of pixel values through vertical transect */
  Sigma_x = 0;                                    
  Sigma_y = 0;

  /*
    ---------------------------
    NEAREST PIXEL TO BARYCENTRE
    ---------------------------
  */
  int_ypos = floor((w_object->ypos)+0.5); 	/* ypos nearest int */
  int_xpos = floor((w_object->xpos)+0.5); 	/* xpos nearest int */


  /*
    --------------------------------
    RUN THROUGH ALL PIXELS IN OBJECT
    --------------------------------
  */
  while(curpix !=NULL)
    {


      /* if on horizontal transect
	 ------------------------- */
      /* 	  if ((curpix->y) == int_ypos){                 /\* if pixel's y position is same as horizontal transect's *\/ */
      /* 	    pix_I = curpix->value;                      /\* pixel value *\/    */
      /* 	    delta_x = curpix->x - xpos;                 /\* dist. btwn x and real barycentre x position (float) *\/ */
      /* 	    dummy_x += (pix_I * (delta_x*delta_x)); */
      /* 	    Sum_I_horiz_cut += pix_I;                   /\* Sum of pixel values through horizontal transect *\/ */
      /* 	  } */

      /* 	  /\* if on vertical transect */
      /* 	     ----------------------- *\/ */
      /* 	  if((curpix->x) == int_xpos){ */
      /* 	    pix_I = curpix->value;                      /\* pixel value *\/  */
      /* 	    delta_y = curpix->y - ypos;                 /\* dist. btwn y and real barycentre y position (float) *\/ */
      /* 	    dummy_y += (pix_I * (delta_y*delta_y)); */
      /* 	    Sum_I_vert_cut += pix_I;                    /\* Sum of pixel values through vertical transect *\/ */

      /* 	  } */
	   

      pix_x = curpix->x;
      pix_y = curpix->y;
      pix_v = curpix->value;
	  
      dummy_x += (pix_v * pow((pix_x-xpos),2));
      dummy_y += (pix_v * pow((pix_y-ypos),2));
      dummy_v += pix_v;
	  
      curpix=curpix->next_pixel;
    }


  /*
    ---------------
    CALCULATE SIGMA
    ---------------
  */
  Sigma_x = sqrt( dummy_x / dummy_v );
  Sigma_y = sqrt( dummy_y / dummy_v );


  /*
    --------------
    CALCULATE FWHM
    --------------
  */
  w_object->fwhmx = 2.3548 * Sigma_x;
  w_object->fwhmy = 2.3548 * Sigma_y;
  (*fwhm) = (w_object->fwhmx + w_object->fwhmy)/2;


	
  fprintf(test_of,"X: %f\t%f\t%f\n",dummy_x,dummy_v,Sigma_x);
  fprintf(test_of,"Y: %f\t%f\t%f\n",dummy_y,dummy_v,Sigma_y);
	
	
  fclose (test_of);

} 












/**
 * Routine, used in conjunction with qsort, to sort floats (FWHM)s.
 */
static int Sort_Float(const void *data1,const void *data2)
{
  return *((float *) data1) < * ((float *) data2);
}

/*
** $Log: not supported by cvs2svn $
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
