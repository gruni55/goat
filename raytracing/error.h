/***************************************************************************
                          error.h  -  description
                             -------------------
    begin                : Fri Jul 27 2001
    copyright            : (C) 2001 by Thomas Weigel
    email                : weigel@lat.ruhr-uni-bochum.de
 ***************************************************************************/

/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/

namespace GOAT
 {
   namespace raytracing
   {
 #define N_ERRORS 8
 #define NO_ERRORS  0
 #define MALLOC_ERR 1
 #define CALLOC_ERR 2
 #define REALLOC_ERR 3
 #define SUPERGITTER 4
 #define IMHOST 5
 #define KUGELGITTER_K 6
 #define NOT_FOUND 7
 
 #ifndef ERROR_H
 #define ERROR_H
 void error (int nerr, char *msg);
 #endif
   }
 }
