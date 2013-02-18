/* --------------------------------------------------------------------------
c	Copyright 2013 Wolfgang Friederich
c
c	This file is part of Gemini II.
c
c	Gemini II is free software: you can redistribute it and/or modify
c	it under the terms of the GNU General Public License as published by
c	the Free Software Foundation, either version 2 of the License, or
c	any later version.
c
c	Gemini II is distributed in the hope that it will be useful,
c	but WITHOUT ANY WARRANTY; without even the implied warranty of
c	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
c	GNU General Public License for more details.
c
c	You should have received a copy of the GNU General Public License
c	along with Gemini II.  If not, see <http://www.gnu.org/licenses/>.
c---------------------------------------------------------------------------- */
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
