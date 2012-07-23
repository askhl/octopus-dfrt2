/*
** Copyright (C) 2010-2012 X. Andrade <xavier@tddft.org>
** 
** FortranCL is free software: you can redistribute it and/or modify
** it under the terms of the GNU Lesser General Public License as published by
** the Free Software Foundation, either version 3 of the License, or
** (at your option) any later version.
**
** FortranCL is distributed in the hope that it will be useful,
** but WITHOUT ANY WARRANTY; without even the implied warranty of
** MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
** GNU Lesser General Public License for more details.
**
** You should have received a copy of the GNU Lesser General Public License
** along with this program.  If not, see <http://www.gnu.org/licenses/>.
**
** $Id$
*/

#include <config.h>
#include <stdlib.h>
#include <stdio.h>
#include "localcl.h"

#include <string_f.h>

typedef void* ptrtype;

/* -----------------------------------------------------------------------*/

void FC_FUNC_(fortrancl_get_component, FORTRANCL_GET_COMPONENT)
     (const ptrtype * array, const int * index, ptrtype * component){
  *component = array[*index];
}


/* -----------------------------------------------------------------------*/

void FC_FUNC_(fortrancl_set_component, FORTRANCL_SET_COMPONENT)
     (ptrtype * array, const int * index, const ptrtype * component){
  array[*index] = *component;
}

/* -----------------------------------------------------------------------*/

void FC_FUNC_(fortrancl_set_null, FORTRANCL_SET_NULL)
     (ptrtype * ptr){
  *ptr = NULL;
}
