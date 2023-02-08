/***************************************************************************
                          error.cpp  -  description
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
#include "error.h"
#include <iostream>
namespace GOAT
{
   namespace raytracing
   {
      const char *error_type[N_ERRORS]={
                                    "Speicherzuordnungsfehler (MALLOC)",
                                    "Speicherzuordnungsfehler (CALLOC)",
                                    "Speicherzuordnungsfehler (REALLOC)",
                                    "allg. Fehler in SuperArray"
                                  };
      void error (int nerr, std::string msg)
      {
         std::cerr << "Programmabbruch !!!!!n" << std::endl;
         std::cerr << "Fehler: " << error_type[nerr] << "  Info: " << msg << std::endl;
         exit (1);
      }
   }
}
