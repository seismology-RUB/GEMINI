/* these are mimicking functions to implement sunfortran routines */

#include <stddef.h>
#include <time.h>

void idate_(int *iarray);
void itime_(int *iarray);
/*
 *--------------------------------------
 *  function mimicking fortran idate
 *--------------------------------------
*/
void idate_(iarray)
	int *iarray;
{
	struct tm *datime;
	time_t zeit;
	zeit = time((time_t *) NULL);
	datime = localtime(&zeit);
	iarray[0] = datime->tm_mday;
	iarray[1] = datime->tm_mon + 1;
	iarray[2] = datime->tm_year;
}
/*
 *--------------------------------------
 *  function mimicking fortran itime
 *--------------------------------------
*/
void itime_(iarray)
	int *iarray;
{
	struct tm *datime;
	time_t zeit;
	zeit = time((time_t *) NULL);
	datime = localtime(&zeit);
	iarray[0] = datime->tm_hour;
	iarray[1] = datime->tm_min;
	iarray[2] = datime->tm_sec;
}
