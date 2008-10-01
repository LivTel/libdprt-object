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
/* object.h
** $Header: /space/home/eng/cjm/cvs/libdprt-object/include/object.h,v 1.7 2008-10-01 15:54:59 eng Exp $
*/
#ifndef OBJECT_H
#define OBJECT_H

/* hash definitions */
/**
 * TRUE is the value usually returned from routines to indicate success.
 */
#ifndef TRUE
#define TRUE 1
#endif
/**
 * FALSE is the value usually returned from routines to indicate failure.
 */
#ifndef FALSE
#define FALSE 0
#endif

/**
 * This is the length of error string of modules in the library.
 */
#define OBJECT_ERROR_STRING_LENGTH	256

/**
 * The number of nanoseconds in one second. A struct timespec has fields in nanoseconds.
 */
#define ONE_SECOND_NS       (1000000000)

/**
 * The number of nanoseconds in one millisecond. A struct timespec has fields in nanoseconds.
 */
#define ONE_MILLISECOND_NS   (1000000)                             

/**
 * The number of milliseconds in one second.
 */
#define ONE_SECOND_MS        (1000)                                

/* Logging bits.
** Note the DpRt java layer reserves bits 0-7 inclusive.
** Note libdprt_rjs reserves bits 8-15 inclusive (currently unused).
** Note libdprt_object reserves bits 16-23 inclusive.
*/
/**
 * Value to pass into logging calls, used to log general execution.
 */
#define OBJECT_LOG_BIT_GENERAL       (1<<16)
/**
 * Value to pass into logging calls, used to log object changing code.
 */
#define OBJECT_LOG_BIT_OBJECT        (1<<17)
/**
 * Value to pass into logging calls, used to log FWHM code.
 */
#define OBJECT_LOG_BIT_FWHM          (1<<18)
/**
 * Value to pass into logging calls, used to log pixel (of object) code.
 */
#define OBJECT_LOG_BIT_PIXEL         (1<<19)
/**
 * Value to pass into logging calls, used to log point list (internal list recursion) code.
 */
#define OBJECT_LOG_BIT_POINT         (1<<20)

/* structures */
/**
 * A structure containing high pixels. These are pixels thats make up an object.
 * The structure contains:
 * <ul>
 * <li><b>x</b> X location of pixel.
 * <li><b>y</b> Y location of pixel.
 * <li><b>value</b> Value above the image median of this pixel.
 * <li><b>next_pixel</b> Pointer to the next pixel in this (linked) list.
 * </ul>
 */
struct HighPixel_Struct
{
	int x;
	int y;
	float value;
	struct HighPixel_Struct *next_pixel;
};

typedef struct HighPixel_Struct HighPixel;

/**
 * A structure containing an object. This object is a collection of contiguous pixels above the threshold.
 * <ul>
 * <li><b>objnum</b> Which number object this is, i.e. each object has a number.
 * <li><b>xpos</b> The x position of the centre of the object, weighted by pixel * value in pixel.
 * <li><b>ypos</b> The y position of the centre of the object, weighted by pixel * value in pixel.
 * <li><b>total</b> The total number of counts above the median for all pixels in the object.
 * <li><b>numpix</b> The number of pixels this object covers.
 * <li><b>peak</b> The number of counts above the median for the brightest pixel in the object.
 * <li><b>is_stellar</b> Boolean determining whether the object is stellar or not.
 * <li><b>fwhmx</b> Full width Half Maximum in X in pixels. Only valid if is_stellar is TRUE. (See below)
 * <li><b>fwhmy</b> Full width Half Maximum in Y in pixels. Only valid if is_stellar is TRUE.  (See below)
 * <li><b>ellipticity</b> A simple ellipticity measure. (A-B)/A.  (See below)
 * <li><b>nextobject</b> A pointer to the next object in the (linked) list.
 * <li><b>highpixel</b> A pointer to a linked list of pixels in the object.
 * <li><b>last_hp</b> A pointer to the end element in the highpixel list.
 * </ul>
 * When using the SExtractor-derived half-flux-radius method of FWHM measures, then both
 * fhhmx and fwhmy will be the same and simly be the object fwhm. Separate values are not
 * returned. The ellipticity value is still valid because A,B are derived internally in the 
 * code.
 */
struct Object_Struct
{
	int objnum;
	float xpos;
	float ypos;
	float total;
	int numpix;
	float peak;
	int is_stellar;
	float fwhmx;
	float fwhmy;
	float ellipticity;
	struct Object_Struct *nextobject;
	HighPixel *highpixel;
	HighPixel *last_hp;
};

/**
 * Object typedef.
 */
typedef struct Object_Struct Object;

/* function declarations */
extern int Object_List_Get(float *image,float image_median,int naxis1,int naxis2,float thresh,int npix,
			   Object **first_object,int *sflag,float *seeing);
extern int Object_List_Free(Object **list);
extern void Object_Error(void);
extern void Object_Error_To_String(char *error_string);
extern int Object_Get_Error_Number(void);
extern void Object_Warning(void);
extern void Object_Get_Current_Time_String(char *time_string,int string_length);
extern void Object_Log_Format(int level,char *format,...);
extern void Object_Log(int level,char *string);
extern void Object_Set_Log_Handler_Function(void (*log_fn)(int level,char *string));
extern void Object_Set_Log_Filter_Function(int (*filter_fn)(int level));
extern void Object_Log_Handler_Stdout(int level,char *string);
extern void Object_Set_Log_Filter_Level(int level);
extern int Object_Log_Filter_Level_Absolute(int level);
extern int Object_Log_Filter_Level_Bitwise(int level);

#endif
