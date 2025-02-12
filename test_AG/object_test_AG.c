/*
         _     _           _      _            _          _    ____ 
    ___ | |__ (_) ___  ___| |_   | |_ ___  ___| |_       / \  / ___|
   / _ \| '_ \| |/ _ \/ __| __|  | __/ _ \/ __| __|     / _ \| |  _ 
  | (_) | |_) | |  __/ (__| |_   | ||  __/\__ \ |_     / ___ \ |_| |
   \___/|_.__// |\___|\___|\__|___\__\___||___/\__|___/_/   \_\____|
            |__/             |_____|             |_____|            

     Copyright 2009, Astrophysics Research Institute, Liverpool John Moores University.

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
/* object_test_AG.c
** $Header: /space/home/eng/cjm/cvs/libdprt-object/test_AG/object_test_AG.c,v 1.2 2025-02-12 11:17:24 eng Exp $
*/


/**
 * object_test_AG.c Tests libdprt_object.so.
 *
 * Autoguider (AG) files do not contain any median and standard deviation headers. Therefore the
 * object_test.c code cannot run on them as it expects these headers to exist. Therefore a separate
 * test code, based on object_test, needs to be made.
 *
 * The only headers in an AG fits are:
 * SIMPLE, BITPIX, NAXIS, NAXIS1, NAXIS2, EXTEND, COMMENT, COMMENT, END
 *
 * If the output filename is specified a FITS image is written with each object mask having it's object number
 * written in.
 * <pre>
 * object_test <filename> [-output <FITS filename>] [-h|--help] [-v|--verbose]
 * </pre>
 */


/*
  $Log: not supported by cvs2svn $
  Revision 1.1  2009/03/04 11:21:55  eng
  Initial revision

  *
*/



#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <math.h>
#include "fitsio.h"
#include "object.h"



/* ------------------------------------------------------- */
/* internal hash definitions */
/* ------------------------------------------------------- */
/**
 * Number of axes in a valid FITS file.
 */
#define FITS_GET_DATA_NAXIS (2)

/* ------------------------------------------------------- */
/* internal functions declarations */
/* ------------------------------------------------------- */
static void Help(void);
static int Parse_Args(int argc,char *argv[]);
static int Load(void);
static int Save(void);
static int Object_Mask_Create(Object *object_list);
static int difftimems(struct timespec start_time,struct timespec stop_time);
static int fltcmp(const void *v1, const void *v2);

/* ------------------------------------------------------- */
/* internal variables */
/* ------------------------------------------------------- */

/* Revision Control System identifier */
static char rcsid[] = "$Id: object_test_AG.c,v 1.2 2025-02-12 11:17:24 eng Exp $";
static char Input_Filename[256] = "";                      /* Filename of file to be processed. */
static char Output_Filename[256] = "";                     /* Filename of file to be output. */
static float *Image_Data = NULL;                           /* Data in image array. */
static float *Image_Data_Sort = NULL;                     /* Data in image array - destructively sorted. */
static unsigned short *Object_Mask_Data = NULL;            /* Data created from object pixel list showing extent of each object. */
static int Naxis1;                                         /* Dimensions of data array. */
static int Naxis2;                                         /* Dimensions of data array. */
static long int ImageSize;                                 /* Data array size (pixels) */
static float Median;                                       /* Background median counts */
static float BGSigma = 7.0;                                /* Threshold sigma level (was 10.0 on 20120403) */
static float Background_SD;                                /* Standard deviation of background */
static float PixelScale = 0.21920;                         /* Pixel scale of binned image (arcsec per binned pixel) */
                                                           /* (was 0.27837 on 20120403) */
static int Log_Level = 0;                                  /* Log level */
static int verbose = FALSE;                                /* Verbose flag (off by default) */

/* ------------------------------------------------------- */
/* external functions */
/* ------------------------------------------------------- */
/**
 * The main program.
 */
int main(int argc, char *argv[])
{
  Object *object_list = NULL;
  Object *tmp_object = NULL;
  Object *object = NULL;
  struct timespec start_time,stop_time;
  int seeing_flag;
  float seeing,thresh;
  int obj_count_init;                   /* initial count of all objects */
  int obj_count_size;                   /* objects bigger than size limit (currently 8 pixels) */
  int obj_count_stellar;                /* objects with ellipticity below limit */
  int obj_count_dia;                    /* objects where fwhm < diameter (calculated from size) */
  int retval;
  float BGSD_factor;
  float peak_abs;
  float fwhmx2,fwhmy2;
  float bx,by,bc;
  float brightest_x,brightest_y,brightest_count;
  long int median_posn;
  long int median_posn_upper,median_posn_lower;
  float total,total2,total3;
  long int i;
  float mean;


  /* TEST ONLY   */
  FILE *fPtr;
  fPtr = fopen("objectfile.txt", "w");


  /*
    ----
    ARGS
    ----
  */

  if(argc < 2){
    Help();
    return 0;
  }

  if(!Parse_Args(argc,argv))
    return 1;
  if(strcmp(Input_Filename,"")==0){
    fprintf(stderr,"object_test: No filename specified.\n");
    return 2;
  }


  Object_Set_Log_Handler_Function(Object_Log_Handler_Stdout);
  Object_Set_Log_Filter_Function(Object_Log_Filter_Level_Absolute);
  Object_Set_Log_Filter_Level(Log_Level);



  /*
    ----------
    LOAD IMAGE
    ----------
    Load the FITS image into a float array.
  */
  if (verbose)
    fprintf(stdout,"object_test: loading image from %s\n",Input_Filename);
  if(!Load())
    return 3;


  /*
    -------------------
    CALCULATE THRESHOLD
    -------------------
    Median & Background_SD are not found automatically
    in AG fitsfiles and have to be calculated. BGSigma
    is set to default to 10.0 sigma.
  */
  if (verbose)
    fprintf(stdout,"object_test: calculating threshold\n");
  
  /*
    Median
    ------
  */  
  ImageSize = (long) Naxis1 * (long) Naxis2;
  qsort (Image_Data_Sort, ImageSize, sizeof(float), fltcmp);
  if (ImageSize % 2){ /* if ODD */
    median_posn = (int) (ImageSize/2);
    Median = (float) Image_Data_Sort[median_posn];
  }
  else { /* if EVEN */
    median_posn_upper = (int) (ImageSize/2);
    median_posn_lower = median_posn_upper - 1;
    Median = (float) ((Image_Data_Sort[median_posn_upper] + Image_Data_Sort[median_posn_lower])/2.0);
  }

  /*
    Background_SD - Standard Deviation of background
    ----
    (might as well use Image_Data_Sort)
  */
  total = 0;
  for (i=0;i<ImageSize;i++)
    total += Image_Data_Sort[i];
  mean = total / (float) ImageSize;

  total2 = 0;
  for (i=0;i<ImageSize;i++)
    total2 += (Image_Data_Sort[i] * Image_Data_Sort[i]);

  total3 = 0;
  for (i=0;i<ImageSize;i++)
    total3 += ((Image_Data_Sort[i] - mean) * (Image_Data_Sort[i] - mean));
  Background_SD = sqrt((total3/(float) ImageSize));

  thresh = Median + (BGSigma * Background_SD);

  if (verbose){
    fprintf(stdout,"object_test: measured background: Median = %.2f and StDev (1 sigma) = %.2f\n",Median,Background_SD);
    fprintf(stdout,"object_test: threshold hard-coded to %.2f sigma in AG config, so thresh passed to object.c = %.2f\n",BGSigma,thresh);
  }


  /*
    -------------
    OBJECT DETECT
    -------------
  */
  if (verbose)
    fprintf(stdout,"object_test: Running object detection....\n");
  clock_gettime(CLOCK_REALTIME,&start_time);
  retval = Object_List_Get(Image_Data,Median,Naxis1,Naxis2,thresh,8,&object_list,&seeing_flag,&seeing);
  clock_gettime(CLOCK_REALTIME,&stop_time);
  if(retval == FALSE){
    Object_Error();
    return 4;
  }


  /* --------------------------------- */
  /* FIND BRIGHTEST OBJECT COORDINATES */
  /* --------------------------------- */
  brightest_x = 0.0;
  brightest_y = 0.0;
  brightest_count = 0.0;
  tmp_object = object_list; /* loop through this copy & leave object_list itself alone */
  while (tmp_object != NULL){
    bc = tmp_object->total;
    if (bc > brightest_count){
      brightest_count = bc;
      brightest_x = tmp_object->xpos;
      brightest_y = tmp_object->ypos;
    }
    tmp_object = tmp_object->nextobject;
  }


  /* ----------------- */
  /* PRINT OUT RESULTS */
  /* ----------------- */
  if (verbose){
    fprintf(stdout,"object_test: The object detection procedure took %d ms.\n",difftimems(start_time,stop_time));
    fprintf(stdout,"object_test: The seeing was %.2f pixels (%.2f arcsec) with seeing_flag = %d (0 is good).\n",
	    seeing,seeing*PixelScale,seeing_flag);
    fprintf(stdout,"object_test: The brightest object was at %.2f,%.2f with %.2f counts.\n",
	    brightest_x,brightest_y,brightest_count);
  }

  /*
    ----------
    FREE IMAGE
    ----------
  */
  if(Image_Data != NULL)
    free(Image_Data);

  /* do object mask output?
     ---------------------- */
  if(strcmp(Output_Filename,"") != 0){
    if(!Object_Mask_Create(object_list))
      return 5;
    if(!Save())
      return 6;
    if(Object_Mask_Data != NULL)
      free(Object_Mask_Data);
  }

  if (verbose)
    fprintf(stdout,"object_test: Freeing object data.\n");

  if(!Object_List_Free(&object_list))
    {
      Object_Error();
      return 7;
    }

  fclose(fPtr);

  if (verbose)
    fprintf(stdout,"object_test: All done.\n");

  return 0;
}






/* ------------------------------------------------------- */
/* internal functionss */
/* ------------------------------------------------------- */


/* ---------------------------------------------------------------------------------------------- */

/**
 * Routine to parse arguments.
 * @param argc The argument count.
 * @param argv The argument list.
 * @return Returns TRUE if the program can proceed, FALSE if it should stop (the user requested help).
 * @see #Input_Filename
 * @see #Output_Filename
 * @see #Median
 * @see #Background_SD
 */
static int Parse_Args(int argc,char *argv[])
{
  int i,retval;
  int call_help = FALSE;

  strcpy(Input_Filename,"");
  strcpy(Output_Filename,"");

  for(i=1;i<argc;i++){

    /* ---- */
    /* HELP */
    /* ---- */
   
    if(strcmp(argv[i],"--help")==0||(strcmp(argv[i],"-h")==0)){
      call_help = TRUE;
    }

    /* ------------ */
    /* VERBOSE FLAG */
    /* ------------ */
    else if ((strcmp(argv[i],"--verbose")==0)||(strcmp(argv[i],"-v")==0)){
      verbose = TRUE;
      fprintf(stdout,"object_test: Parse_Args: verbose ON\n");
    }


    /* ------------- */
    /* LOGGING LEVEL */
    /* ------------- */
    
    else if((strcmp(argv[i],"-log_level")==0)||(strcmp(argv[i],"-l")==0)){
      if((i+1) < argc){
	retval = sscanf(argv[i+1],"%d",&Log_Level);
	if(retval != 1){
	  fprintf(stderr,"object_test: Parse_Args: log level parameter %s not an integer.\n",argv[i+1]);
	  return FALSE;
	}
	i++;
      }
      else {
	fprintf(stderr,"object_test: Parse_Args: log level parameter missing.\n");
	return FALSE;
      }
    }
    

    /* ----------------------- */
    /* OBJECT FITS OUTPUT FILE */
    /* ----------------------- */

    else if((strcmp(argv[i],"-output")==0)||(strcmp(argv[i],"-o")==0)){
      if((i+1) < argc){
	strcpy(Output_Filename,argv[i+1]);
	i++;
      }
      else {
	fprintf(stderr,"object_test: Parse_Args: output filename missing.\n");
	return FALSE;
      }
    }



    /* --------------- */
    /* INPUT FITS FILE */
    /* --------------- */

    else
      strcpy(Input_Filename,argv[i]);
  }
  

  /* ------- */
  /* DO HELP */
  /* ------- */
  if(call_help){
    Help();
    return FALSE;
  }
  


  return TRUE;
}


/* ---------------------------------------------------------------------------------------------- */


/**
 * Routine to produce some help.
 */
static void Help(void)
{
  fprintf(stdout,"object_test: Tests the object finding routine in libdprt_object.\n");
  fprintf(stdout,"object_test [-h[elp]] [-v[erbose]] [-l[og_level] <level>]\n");
  fprintf(stdout,"\t<FITS filename>  [-o[utput] <FITS filename>]\n");
  fprintf(stdout,"-help prints this help message and exits.\n");
  fprintf(stdout,"-verbose prints progress to stdout (default off).\n");
  fprintf(stdout,"-log_level sets the amount of logging produced.\n");
  fprintf(stdout,"-threshold sets the threshold level in counts\n");
  fprintf(stdout,"-sigma sets the threshold level in sigma (default 10.0)\n");
  fprintf(stdout,"-output writes an object mask to the specified FITS filename.\n");
  fprintf(stdout,"You must always specify a filename to reduce.\n");
  fprintf(stdout,"Ideally, pass in a flat-fielded de-biased image.\n");
}



/* ---------------------------------------------------------------------------------------------- */

/**
 * Load the FITS image into a float array.
 * @return TRUE on success, FALSE on failure.
 * @see #Input_Filename
 * @see #Naxis1
 * @see #Naxis2
 * @see #Image_Data
 */
static int Load(void)
{
  fitsfile *fits_fp = NULL;
  int retval=0,status=0,integer_value;

  /* 
     open file 
     ---------
  */
  retval = fits_open_file(&fits_fp,Input_Filename,READONLY,&status);
  if(retval)
    {
      fits_report_error(stderr,status);
      return FALSE;
    }

  /*
    check axes
    ----------
  */
  retval = fits_read_key(fits_fp,TINT,"NAXIS",&integer_value,NULL,&status);
  if(retval)
    {
      fits_report_error(stderr,status);
      fits_close_file(fits_fp,&status);
      return FALSE;
    }
  if(integer_value != FITS_GET_DATA_NAXIS)
    {
      fprintf(stderr,"object_test: Wrong NAXIS value(%d).\n",integer_value);
      fits_close_file(fits_fp,&status);
      return FALSE;
    }


  /*
    get naxis1,naxis2
    -----------------
  */
  retval = fits_read_key(fits_fp,TINT,"NAXIS1",&Naxis1,NULL,&status);
  if(retval)
    {
      fits_report_error(stderr,status);
      fits_close_file(fits_fp,&status);
      return FALSE;
    }
  retval = fits_read_key(fits_fp,TINT,"NAXIS2",&Naxis2,NULL,&status);
  if(retval)
    {
      fits_report_error(stderr,status);
      fits_close_file(fits_fp,&status);
      return FALSE;
    }


  /*
    allocate image memory
    ---------------------
  */

  Image_Data = (float *)malloc(Naxis1*Naxis2*sizeof(float));
  Image_Data_Sort = (float *)malloc(Naxis1*Naxis2*sizeof(float));
  if(Image_Data == NULL)
    {
      fprintf(stderr,"object_test: failed to allocate memory (%d,%d).\n",
	      Naxis1,Naxis2);
      fits_close_file(fits_fp,&status);
      return FALSE;
    }
  if(Image_Data_Sort == NULL)
    {
      fprintf(stderr,"object_test: failed to allocate memory (%d,%d).\n",
	      Naxis1,Naxis2);
      fits_close_file(fits_fp,&status);
      return FALSE;
    }


  /*
    read data
    ---------
  */
  retval = fits_read_img(fits_fp,TFLOAT,1,Naxis1*Naxis2,NULL,Image_Data,NULL,&status);
  if(retval)
    {
      fits_report_error(stderr,status);
      fits_close_file(fits_fp,&status);
      fprintf(stderr,"object_test:fits_read_img:Failed to read FITS (%d,%d).\n",Naxis1,Naxis2);
      return FALSE;
    }
  retval = fits_read_img(fits_fp,TFLOAT,1,Naxis1*Naxis2,NULL,Image_Data_Sort,NULL,&status);
  if(retval)
    {
      fits_report_error(stderr,status);
      fits_close_file(fits_fp,&status);
      fprintf(stderr,"object_test:fits_read_img:Failed to read FITS (%d,%d).\n",Naxis1,Naxis2);
      return FALSE;
    }


  /* 
     close file 
     ----------
  */
  retval = fits_close_file(fits_fp,&status);
  if(retval)
    {
      fits_report_error(stderr,status);
      return FALSE;
    }
  return TRUE;
}



/* ---------------------------------------------------------------------------------------------- */


/**
 * Save the object mask data into a FITS image.
 * @return TRUE on success, FALSE on failure.
 * @see #Output_Filename
 * @see #Naxis1
 * @see #Naxis2
 * @see #Object_Mask_Data
 */
static int Save(void)
{
  fitsfile *fits_fp = NULL;
  int retval=0,status=0,ivalue;
  double dvalue;

  if(Object_Mask_Data == NULL)
    {
      fprintf(stderr,"Object mask is NULL.\n");
      return FALSE;
    }
  /* open file */
  retval = fits_create_file(&fits_fp,Output_Filename,&status);
  if(retval)
    {
      fits_report_error(stderr,status);
      return FALSE;
    }
  /* ensure the simple keywords are correct */
  /* SIMPLE */
  ivalue = TRUE;
  retval = fits_update_key(fits_fp,TLOGICAL,(char*)"SIMPLE",&ivalue,NULL,&status);
  if(retval != 0)
    {
      fits_report_error(stderr, status);
      fprintf(stderr,"SIMPLE keyword.\n");
      return FALSE;
    }
  /* BITPIX */
  ivalue = 16;
  retval = fits_update_key(fits_fp,TINT,(char*)"BITPIX",&ivalue,NULL,&status);
  if(retval != 0)
    {
      fits_report_error(stderr, status);
      fprintf(stderr,"BITPIX keyword.\n");
      return FALSE;
    }
  /* NAXIS */
  ivalue = 2;
  retval = fits_update_key(fits_fp,TINT,(char*)"NAXIS",&ivalue,NULL,&status);
  if(retval != 0)
    {
      fits_report_error(stderr, status);
      fprintf(stderr,"NAXIS keyword.\n");
      return FALSE;
    }
  /* NAXIS1 */
  ivalue = Naxis1;
  retval = fits_update_key(fits_fp,TINT,(char*)"NAXIS1",&ivalue,NULL,&status);
  if(retval != 0)
    {
      fits_report_error(stderr, status);
      fprintf(stderr,"NAXIS1 keyword = %d.\n",ivalue);
      return FALSE;
    }
  /* NAXIS2 */
  ivalue = Naxis2;
  retval = fits_update_key(fits_fp,TINT,(char*)"NAXIS2",&ivalue,NULL,&status);
  if(retval != 0)
    {
      fits_report_error(stderr, status);
      fprintf(stderr,"NAXIS2 keyword = %d.\n",ivalue);
      return FALSE;
    }
  /* BZERO */
  dvalue = 32768.0;
  retval = fits_update_key_fixdbl(fits_fp,"BZERO",dvalue,6,NULL,&status);
  if(retval != 0)
    {
      fits_report_error(stderr, status);
      fprintf(stderr,"BZERO keyword = %.2f.\n",dvalue);
      return FALSE;
    }
  /* BSCALE */
  dvalue = 1.0;
  retval = fits_update_key_fixdbl(fits_fp,"BSCALE",dvalue,6,NULL,&status);
  if(retval != 0)
    {
      fits_report_error(stderr, status);
      fprintf(stderr,"BSCALE keyword = %.2f.\n",dvalue);
      return FALSE;
    }
  /* write image */
  retval = fits_write_img(fits_fp,TUSHORT,1,Naxis1*Naxis2,Object_Mask_Data,&status);
  if(retval)
    {
      fits_report_error(stderr,status);
      return FALSE;
    }
  /* close file */
  retval = fits_close_file(fits_fp,&status);
  if(retval)
    {
      fits_report_error(stderr,status);
      return FALSE;
    }
  return TRUE;
}



/* ---------------------------------------------------------------------------------------------- */


/**
 * Create an output object mask from pixel data in object list.
 * @param object_list The objects to create.
 * @return TRUE on success, FALSE on failure.
 * @see #Naxis1
 * @see #Naxis2
 * @see #Object_Mask_Data
 */
static int Object_Mask_Create(Object *object_list)
{
  Object *object = NULL;
  HighPixel *high_pixel = NULL;
  int x,y,objnum;

  Object_Mask_Data = (unsigned short*)malloc(Naxis1*Naxis2*sizeof(unsigned short));
  if(Object_Mask_Data == NULL)
    {
      fprintf(stderr,"Failed to allocate Object Mask(%d,%d).\n",Naxis1,Naxis2);
      return FALSE;
    }
  object = object_list;
  while(object != NULL)
    {
      objnum = object->objnum;
      objnum = (objnum % 65536);
      high_pixel = object->highpixel;
      while(high_pixel != NULL)
	{
	  x = high_pixel->x;
	  y = high_pixel->y;
	  if((x > -1) &&(x < Naxis1)&&(y > -1) &&(y < Naxis2))
	    Object_Mask_Data[(y*Naxis1)+x] = objnum;
	  /* and to next pixel */
	  high_pixel = high_pixel->next_pixel;
	}

      /* highlight centre-point */
      x = (int)(object->xpos);
      y = (int)(object->ypos);
      if((x > -1) &&(x < Naxis1)&&(y > -1) &&(y < Naxis2))
	Object_Mask_Data[(y*Naxis1)+x] = 65535;
      /* and to next object */
      object = object->nextobject;
    }
  return TRUE;
}


/* ---------------------------------------------------------------------------------------------- */

/*
     _ _  __  __ _   _                              
  __| (_)/ _|/ _| |_(_)_ __ ___   ___ _ __ ___  ___ 
 / _` | | |_| |_| __| | '_ ` _ \ / _ \ '_ ` _ \/ __|
| (_| | |  _|  _| |_| | | | | | |  __/ | | | | \__ \
 \__,_|_|_| |_|  \__|_|_| |_| |_|\___|_| |_| |_|___/
                                                    
*/
/**
 * Routine to calculate the difference between start_time and stop_time, and to return 
 * the number of milliseconds difference.
 * @param start_time The start time.
 * @param stop_time The end time.
 * @return The number of milliseconds between start_time and stop_time.
 * @see #ONE_SECOND_NS
 * @see #ONE_SECOND_MS
 * @see #ONE_MILLISECOND_NS
 */
static int difftimems(struct timespec start_time,struct timespec stop_time)
{
  int sec,ns,ms;

  sec = stop_time.tv_sec-start_time.tv_sec;
  ns = stop_time.tv_nsec-start_time.tv_nsec;
  if(ns < 0)
    {
      ns += ONE_SECOND_NS;
      sec -= 1;
    }
  ms = (sec*ONE_SECOND_MS)+(ns/ONE_MILLISECOND_NS);
  return ms;
}


/*
  __ _ _                       
 / _| | |_ ___ _ __ ___  _ __  
| |_| | __/ __| '_ ` _ \| '_ \ 
|  _| | || (__| | | | | | |_) |
|_| |_|\__\___|_| |_| |_| .__/ 
                        |_|    
*/
int fltcmp(const void *v1, const void *v2){
  return (*(int *)v1 - *(int *)v2);
}
