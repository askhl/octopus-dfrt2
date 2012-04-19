/* bravais_art.c */
/* Copyright (C) 2008 Atsushi Togo */

/* This program is free software; you can redistribute it and/or */
/* modify it under the terms of the GNU General Public License */
/* as published by the Free Software Foundation; either version 2 */
/* of the License, or (at your option) any later version. */

/* This program is distributed in the hope that it will be useful, */
/* but WITHOUT ANY WARRANTY; without even the implied warranty of */
/* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the */
/* GNU General Public License for more details. */

/* You should have received a copy of the GNU General Public License */
/* along with this program; if not, write to the Free Software */
/* Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA. */

#include <stdio.h>
#include <stdlib.h>
#include "bravais.h"
#include "cell.h"
#include "debug.h"
#include "mathfunc.h"
#include "spacegroup.h"
#include "symmetry.h"

static double bcc_axes[13][3] = {
  { 1.0, 0.0, 0.0},
  { 0.0, 1.0, 0.0},
  { 0.0, 0.0, 1.0},
  { 0.0, 1.0, 1.0},
  { 0.0, 1.0,-1.0},
  { 1.0, 0.0, 1.0},
  {-1.0, 0.0, 1.0},
  { 1.0, 1.0, 0.0},
  { 1.0,-1.0, 0.0},
  {-0.5, 0.5, 0.5},
  { 0.5,-0.5, 0.5},
  { 0.5, 0.5,-0.5},
  { 0.5, 0.5, 0.5},
};

static double fcc_axes[13][3] = {
  { 1.0, 0.0, 0.0},
  { 0.0, 1.0, 0.0},
  { 0.0, 0.0, 1.0},
  { 0.0, 0.5, 0.5},
  { 0.0, 0.5,-0.5},
  { 0.5, 0.0, 0.5},
  {-0.5, 0.0, 0.5},
  { 0.5, 0.5, 0.0},
  { 0.5,-0.5, 0.0},
  { 1.0, 1.0, 1.0},
  {-1.0, 1.0, 1.0},
  { 1.0,-1.0, 1.0},
  { 1.0, 1.0,-1.0},
};

static double primitive_axes[13][3] = {
  { 1.0, 0.0, 0.0},
  { 0.0, 1.0, 0.0},
  { 0.0, 0.0, 1.0},
  { 0.0, 1.0, 1.0},
  { 0.0, 1.0,-1.0},
  { 1.0, 0.0, 1.0},
  {-1.0, 0.0, 1.0},
  { 1.0, 1.0, 0.0},
  { 1.0,-1.0, 0.0},
  { 1.0, 1.0, 1.0},
  {-1.0, 1.0, 1.0},
  { 1.0,-1.0, 1.0},
  { 1.0, 1.0,-1.0},
};

static int is_holohedry(Bravais *bravais, const Cell *cell, const Holohedry holohedry, const double symprec);
static int get_rotation_axis(const int rot[3][3], const int axis_num);
static int is_monocli(Bravais *bravais, const Symmetry *conv_sym, const double symprec);
static int is_monocli_from_I(Bravais *bravais, const Symmetry *conv_sym, const double symprec);
static int is_monocli_from_F(Bravais *bravais, const Symmetry *conv_sym, const double symprec);
static int is_monocli_from_P(Bravais *bravais, const Symmetry *conv_sym, const double symprec);
static int get_monocli_bravais(double lattice[3][3], const Symmetry *conv_sym,
			       const double relative_vol,  const
			       double monocli_axes[13][3],
			       const double symprec);
static void get_monocli_bcc_to_c_center(double lattice[3][3]);
static int is_monocli_orthogonal(const int b_axis, const int naxis, const double monocli_axes[13][3], const double symprec);
static void get_monocli_relative_axes(double relative_axis[3][3], const double new_lattice[3][3],
				      const double old_lattice[3][3], const double symprec);
static int is_tetra(Bravais *bravais, const Symmetry *conv_sym, const double symprec);
static int is_ortho(Bravais *bravais, const Symmetry *conv_sym, const double symprec);
static int is_ortho_from_H(Bravais *bravais, const Cell *cell, const Symmetry *symmetry, const double symprec);
static int is_ortho_from_H_axis(const Bravais *bravais, const Cell *cell, const Symmetry *symmetry, const double symprec);
static int is_ortho_from_I(Bravais *bravais, const Symmetry *conv_sym, const double symprec);
static int is_ortho_from_F(Bravais *bravais, const Symmetry *conv_sym, const double symprec);
static int is_ortho_from_P(Bravais *bravais, const Symmetry *conv_sym, const double symprec);
static int get_ortho_axis(int naxis[3], const Symmetry *conv_sym);

/* bravais is going to be changed. */
int art_get_artificial_bravais(Bravais *bravais, const Symmetry *symmetry,
			       const Cell *cell, const Holohedry holohedry,
			       const double symprec)
{
  int i, j;
  double lattice[3][3];
  Symmetry conv_sym, tmp_sym;
  
  /* Triogonal */
  if (holohedry == TRIGO && bravais->holohedry == HEXA) {
    /* Just change the holohedry. */
    goto found;
  }

  /* Rhombohedral */
  if (holohedry == TRIGO && bravais->holohedry == RHOMB) {
    /* Do nothing */
    goto end;
  }

  /* Rhombohedral from Cubic */
  if (holohedry == TRIGO && bravais->holohedry == CUBIC ) {
    if (is_holohedry(bravais, cell, RHOMB, symprec)) {
      goto end;
    }
    goto not_found;
  }

  /* Orthorhombic from Hexagonal */
  if (holohedry == ORTHO && bravais->holohedry == HEXA ) {

    if (is_ortho_from_H(bravais, cell, symmetry, symprec)) {
      goto end;
    }
    goto not_found;
  }


  conv_sym = tbl_get_conventional_symmetry(bravais, cell, symmetry, symprec);

#ifdef DEBUG
  debug_print("Rotation of conventional cell of the following Bravais lattice\n");
  debug_print_matrix_d3(bravais->lattice);
  for ( i = 0; i < conv_sym.size; i++ ) {
    debug_print("---- %d ----\n", i+1 );
    debug_print_matrix_i3(conv_sym.rot[i]);
  }
#endif

  /* Monoclinic */
  if ( holohedry == MONOCLI ) {
    if (is_monocli(bravais, &conv_sym, symprec)) {
      goto found_and_deallocate;
    }
    goto not_found_and_deallocate;
  }

  /* Tetragonal */
  if (holohedry == TETRA) {
    if (is_tetra(bravais, &conv_sym, symprec)) {
      goto found_and_deallocate;
    }
    goto not_found_and_deallocate;
  }

  /* Orthorhombic */
  if (holohedry == ORTHO) {
    if (is_ortho(bravais, &conv_sym, symprec)) {
      goto found_and_deallocate;
    }

    goto not_found_and_deallocate;
  }

  sym_delete_symmetry(&conv_sym);

  /* Triclinic */
  if (is_holohedry(bravais, cell, holohedry, symprec)) {
    goto end;
  }

  goto not_found;



  /*************/
  /*** Found ***/
  /*************/

 found_and_deallocate:
  sym_delete_symmetry(&conv_sym);
 found:
  bravais->holohedry = holohedry;
 end:
  /* Check if right hand system  */
  if (mat_get_determinant_d3(bravais->lattice) < 0) {
    for (i = 0; i < 3; i++)
      for (j = 0; j < 3; j++)
	bravais->lattice[i][j] = -bravais->lattice[i][j];
  }
  return 1;


  /*****************/
  /*** Not found ***/
  /*****************/

 not_found_and_deallocate:
  sym_delete_symmetry(&conv_sym);
 not_found:
  return 0;
}

static int is_ortho(Bravais *bravais, const Symmetry *conv_sym, const double symprec)
{
  /* I-Ortho and F-Ortho from I-Cbuic or I-Tetragonal */
  if (bravais->centering == BODY) {
    if (is_ortho_from_I(bravais, conv_sym, symprec))
      goto found;
  }  

  if (bravais->centering == FACE) {
    if (is_ortho_from_F(bravais, conv_sym, symprec))
      goto found;
  }

  if (bravais->centering == NO_CENTER) {
    if (is_ortho_from_P(bravais, conv_sym, symprec))
      goto found;
  }

  return 0;

 found:
  return 1;
}

static int is_ortho_from_H(Bravais *bravais, const Cell *cell,
			   const Symmetry *symmetry, const double symprec)
{
  int i;
  double bravais_lattice[3][3];

  mat_copy_matrix_d3(bravais_lattice, bravais->lattice);
  bravais->holohedry = ORTHO;
  bravais->centering = C_FACE;

  /* Try three kinds of C-base orthorhombic lattice */
  /* Type 1 (no rotation) */
  for ( i = 0; i < 3; i++ ) {
    bravais->lattice[i][0] = bravais_lattice[i][0] - bravais_lattice[i][1];
    bravais->lattice[i][1] = bravais_lattice[i][0] + bravais_lattice[i][1];
    bravais->lattice[i][2] = bravais_lattice[i][2];
  }
  if ( is_ortho_from_H_axis(bravais, cell, symmetry, symprec) )
    goto found;

  /* Type 2 (60 degs) */
  for ( i = 0; i < 3; i++ ) {
    bravais->lattice[i][0] = 2* bravais_lattice[i][0] + bravais_lattice[i][1];
    bravais->lattice[i][1] = bravais_lattice[i][1];
    bravais->lattice[i][2] = bravais_lattice[i][2];
  }
  if ( is_ortho_from_H_axis(bravais, cell, symmetry, symprec) )
    goto found;

  /* Type 3 (120 degs) */
  for ( i = 0; i < 3; i++ ) {
    bravais->lattice[i][0] = bravais_lattice[i][0] + 2 * bravais_lattice[i][1];
    bravais->lattice[i][1] = - bravais_lattice[i][0];
    bravais->lattice[i][2] = bravais_lattice[i][2];
  }
  if ( is_ortho_from_H_axis(bravais, cell, symmetry, symprec) )
    goto found;


  /* Not found */
  mat_copy_matrix_d3(bravais->lattice, bravais_lattice);
  return 0;

  /* Found */
 found:
  return 1;
}

static int is_ortho_from_H_axis(const Bravais *bravais, const Cell *cell,
			   const Symmetry *symmetry, const double symprec)
{
  int naxis[3];
  Symmetry conv_sym;

  conv_sym = tbl_get_conventional_symmetry(bravais, cell, symmetry, symprec);

  if ( get_ortho_axis(naxis, &conv_sym) ) {
    /* Found */
    sym_delete_symmetry(&conv_sym);
    return 1;
  }

  /* Not found */
  sym_delete_symmetry(&conv_sym);
  return 0;
}

static int get_ortho_axis(int naxis[3], const Symmetry *conv_sym)
{
  int i, tmp_naxis, num_ortho_axis = 0;

  for (i = 0; i < conv_sym->size; i++) {
    tmp_naxis = get_rotation_axis(conv_sym->rot[i], 2);

    if (num_ortho_axis == 0) {
      if (tmp_naxis > -1) {
	naxis[0] = tmp_naxis;
	num_ortho_axis++;
      }
    }

    if (num_ortho_axis == 1) {
      if (tmp_naxis != naxis[0]) {
	naxis[1] = tmp_naxis;
	num_ortho_axis++;
      }
    }

    if (num_ortho_axis == 2) {
      if (tmp_naxis != naxis[0] && tmp_naxis != naxis[1]) {
	naxis[2] = tmp_naxis;
	num_ortho_axis++;
      }
    }

    if (num_ortho_axis > 2) {
      debug_print("Ortho axes(%d): %d %d %d\n", num_ortho_axis, naxis[0], naxis[1], naxis[2]);
      return 1;
    }
  }

  return 0;
}

static int is_ortho_from_I(Bravais *bravais, const Symmetry *conv_sym,
			   const double symprec)
{
  int i, j, naxis[3];
  double relative_axis[3][3];

  if (! get_ortho_axis(naxis, conv_sym))
    goto not_found;

  /* Each axis has at least 2-fold rotation. */
  for (i = 0; i < 3; i++)
    for (j = 0; j < 3; j++)
      relative_axis[j][i] = bcc_axes[naxis[i]][j];

  debug_print("is_ortho_from_I\n");
  debug_print_matrix_d3(relative_axis);

  mat_multiply_matrix_d3(bravais->lattice, bravais->lattice, relative_axis);

  return 1;

 not_found:
  return 0;
}

static int is_ortho_from_F(Bravais *bravais, const Symmetry *conv_sym,
			   const double symprec)
{
  int i, j, naxis[3];
  double relative_axis[3][3];

  if (! get_ortho_axis(naxis, conv_sym))
    goto not_found;
      
  /* Each axis has at least 2-fold rotation. */
  for (i = 0; i < 3; i++)
    for (j = 0; j < 3; j++)
      relative_axis[j][i] = fcc_axes[naxis[i]][j];

  debug_print("is_ortho_from_F\n");
  debug_print_matrix_d3(relative_axis);
  mat_multiply_matrix_d3(bravais->lattice, bravais->lattice, relative_axis);

  /* I case */
  if (mat_Dabs(mat_Dabs(mat_get_determinant_d3(relative_axis)) - 0.5) < symprec) {
    bravais->centering = BODY;
  }

  return 1;

 not_found:
  return 0;
}

static int is_ortho_from_P(Bravais *bravais, const Symmetry *conv_sym,
			   const double symprec)
{
  int i, j, naxis[3];
  double relative_axis[3][3];

  if (! get_ortho_axis(naxis, conv_sym))
    goto not_found;
      
  /* Each axis has at least 2-fold rotation. */
  for (i = 0; i < 3; i++)
    for (j = 0; j < 3; j++)
      relative_axis[j][i] = primitive_axes[naxis[i]][j];
  mat_multiply_matrix_d3(bravais->lattice, bravais->lattice, relative_axis);
  if (naxis[0] < 3 &&  naxis[1] > 2 && naxis[2] > 2)
    bravais->centering = A_FACE;
  if (naxis[0] > 2 &&  naxis[1] < 3 && naxis[2] > 2)
    bravais->centering = B_FACE;
  if (naxis[0] > 2 &&  naxis[1] > 2 && naxis[2] < 3)
    bravais->centering = C_FACE;

  debug_print("is_ortho_from_P\n");
  debug_print_matrix_d3(relative_axis);

  return 1;

 not_found:
  return 0;
}

static int is_tetra(Bravais *bravais, const Symmetry *conv_sym,
		    const double symprec)
{
  int i, tmp_naxis;
  double permutate_axis_a[3][3] = {
    { 0.0, 0.0, 1.0},
    { 1.0, 0.0, 0.0},
    { 0.0, 1.0, 0.0},
  };
  double permutate_axis_b[3][3] = {
    { 0.0, 1.0, 0.0},
    { 0.0, 0.0, 1.0},
    { 1.0, 0.0, 0.0},
  };
  double permutate_axis_F_a[3][3] = {
    { 0.0, 0.0, 1.0},
    { 0.5, 0.5, 0.0},
    {-0.5, 0.5, 0.0},
  };
  double permutate_axis_F_b[3][3] = {
    {-0.5, 0.5, 0.0},
    { 0.0, 0.0, 1.0},
    { 0.5, 0.5, 0.0},
  };
  double permutate_axis_F_c[3][3] = {
    { 0.5, 0.5, 0.0},
    {-0.5, 0.5, 0.0},
    { 0.0, 0.0, 1.0},
  };

  /* Get 4 fold axis */
  for (i = 0; i < conv_sym->size; i++) {
    tmp_naxis = get_rotation_axis(conv_sym->rot[i], 4);
    if (tmp_naxis > -1)
      break;
  }

  if (tmp_naxis < 0)
    goto not_found;

  /* I from F-Cubic */
  if (bravais->centering == FACE) {
    if (tmp_naxis == 0 || tmp_naxis == 1 || tmp_naxis == 2) {
      if (tmp_naxis == 0) { /* b c a*/
	mat_multiply_matrix_d3(bravais->lattice, bravais->lattice, permutate_axis_F_a);
      }
      if (tmp_naxis == 1) { /* c a b*/
	mat_multiply_matrix_d3(bravais->lattice, bravais->lattice, permutate_axis_F_b);
      }
      if (tmp_naxis == 2) { /* a b c*/
	mat_multiply_matrix_d3(bravais->lattice, bravais->lattice, permutate_axis_F_c);
      }
      
      bravais->centering = BODY;
      goto found;
    }
  }

  /* P or I: just exchange the 4-fold axis to c-axis */
  if (bravais->centering != FACE) {
    if (tmp_naxis == 0 || tmp_naxis == 1 || tmp_naxis == 2) {
      if (tmp_naxis == 0) { /* c a b*/
	mat_multiply_matrix_d3(bravais->lattice, bravais->lattice, permutate_axis_a);
      }
      if (tmp_naxis == 1) { /* b c a*/
	mat_multiply_matrix_d3(bravais->lattice, bravais->lattice, permutate_axis_b);
      }

      debug_print("is_tetra\n");
      debug_print_matrix_d3(bravais->lattice);
      goto found;
    }
  }

 not_found: /* This may not happen. */
  return 0;

 found:
  return 1;
}

static int is_monocli(Bravais *bravais, const Symmetry *conv_sym,
		      const double symprec)
{
  /* Monoclinic from I-cubic, I-tetra, and I-ortho */
  if (bravais->centering == BODY) {
    if (is_monocli_from_I(bravais, conv_sym, symprec))
      goto found;
  }

  /* Monoclinic from F-cubic */
  if (bravais->centering == FACE) {
    if (is_monocli_from_F(bravais, conv_sym, symprec))
      goto found;
  }

  /* Monoclinic from P */
  if (bravais->centering == NO_CENTER || A_FACE || B_FACE || C_FACE) {
    if (is_monocli_from_P(bravais, conv_sym, symprec))
      goto found;
  }

  return 0;

 found:
  return 1;
}

static int is_monocli_from_I(Bravais *bravais, const Symmetry *conv_sym,
			     const double symprec)
{
  double lattice[3][3], relative_axis[3][3];

  mat_copy_matrix_d3(lattice, bravais->lattice);
  
  /* base center */
  if (get_monocli_bravais(bravais->lattice, conv_sym, 1.0, bcc_axes, symprec)) {
    get_monocli_relative_axes(relative_axis, bravais->lattice, lattice, symprec);
    if (mat_Dabs(mat_Dabs(relative_axis[0][0]+relative_axis[0][1]) - 1.0) < symprec &&
	mat_Dabs(mat_Dabs(relative_axis[1][0]+relative_axis[1][1]) - 1.0) < symprec &&
	mat_Dabs(mat_Dabs(relative_axis[2][0]+relative_axis[2][1]) - 1.0) < symprec) {
      bravais->centering = C_FACE;
      goto found;
    }
    if (mat_Dabs(mat_Dabs(relative_axis[0][1]+relative_axis[0][2]) - 1.0) < symprec &&
	mat_Dabs(mat_Dabs(relative_axis[1][1]+relative_axis[1][2]) - 1.0) < symprec &&
	mat_Dabs(mat_Dabs(relative_axis[2][1]+relative_axis[2][2]) - 1.0) < symprec) {
      bravais->centering = A_FACE;
      goto found;
    }	

    /* otherwise bcc */
    get_monocli_bcc_to_c_center(bravais->lattice);
    bravais->centering = C_FACE;
    goto found;
  }

  /* primitive */
  if (get_monocli_bravais(bravais->lattice, conv_sym, 0.5, bcc_axes, symprec)) {
    bravais->centering = NO_CENTER;
    goto found;
  }

  return 0;

 found:
  return 1;
}

static int is_monocli_from_F(Bravais *bravais, const Symmetry *conv_sym,
			     const double symprec)
{
  double lattice[3][3], relative_axis[3][3];

  mat_copy_matrix_d3(lattice, bravais->lattice);
  
  /* base center */
  if (get_monocli_bravais(bravais->lattice, conv_sym, 0.5, fcc_axes, symprec)) {
    get_monocli_relative_axes(relative_axis, bravais->lattice, lattice, symprec);

    if ((mat_Dabs(mat_Dabs(relative_axis[0][0]+relative_axis[0][1]) - 0.0) < symprec &&
	 mat_Dabs(mat_Dabs(relative_axis[1][0]+relative_axis[1][1]) - 1.0) < symprec &&
	 mat_Dabs(mat_Dabs(relative_axis[2][0]+relative_axis[2][1]) - 1.0) < symprec) ||
	(mat_Dabs(mat_Dabs(relative_axis[0][0]+relative_axis[0][1]) - 1.0) < symprec &&
	 mat_Dabs(mat_Dabs(relative_axis[1][0]+relative_axis[1][1]) - 0.0) < symprec &&
	 mat_Dabs(mat_Dabs(relative_axis[2][0]+relative_axis[2][1]) - 1.0) < symprec) ||
	(mat_Dabs(mat_Dabs(relative_axis[0][0]+relative_axis[0][1]) - 1.0) < symprec &&
	 mat_Dabs(mat_Dabs(relative_axis[1][0]+relative_axis[1][1]) - 1.0) < symprec &&
	 mat_Dabs(mat_Dabs(relative_axis[2][0]+relative_axis[2][1]) - 0.0) < symprec)) {
      bravais->centering = C_FACE;
      goto found;
    }
    if ((mat_Dabs(mat_Dabs(relative_axis[0][1]+relative_axis[0][2]) - 0.0) < symprec &&
	 mat_Dabs(mat_Dabs(relative_axis[1][1]+relative_axis[1][2]) - 1.0) < symprec &&
	 mat_Dabs(mat_Dabs(relative_axis[2][1]+relative_axis[2][2]) - 1.0) < symprec) ||
	(mat_Dabs(mat_Dabs(relative_axis[0][1]+relative_axis[0][2]) - 1.0) < symprec &&
	 mat_Dabs(mat_Dabs(relative_axis[1][1]+relative_axis[1][2]) - 0.0) < symprec &&
	 mat_Dabs(mat_Dabs(relative_axis[2][1]+relative_axis[2][2]) - 1.0) < symprec) ||
	(mat_Dabs(mat_Dabs(relative_axis[0][1]+relative_axis[0][2]) - 1.0) < symprec &&
	 mat_Dabs(mat_Dabs(relative_axis[1][1]+relative_axis[1][2]) - 1.0) < symprec &&
	 mat_Dabs(mat_Dabs(relative_axis[2][1]+relative_axis[2][2]) - 0.0) < symprec)) {
      bravais->centering = A_FACE;
      goto found;
    }
    /* otherwise bcc */
    get_monocli_bcc_to_c_center(bravais->lattice);
    bravais->centering = C_FACE;
    goto found;
  }

  /* primitive */
  if (get_monocli_bravais(bravais->lattice, conv_sym, 0.25, fcc_axes, symprec)) {
    bravais->centering = NO_CENTER;
    goto found;
  }

  return 0;

 found:
  return 1;
}

static int is_monocli_from_P(Bravais *bravais, const Symmetry *conv_sym, 
			     const double symprec)
{
  double lattice[3][3], relative_axis[3][3];

  /* primitive */
  if (get_monocli_bravais(bravais->lattice, conv_sym, 1.0, primitive_axes, symprec)) {
    bravais->centering = NO_CENTER;
    goto found;
  }

  /* base center */
  mat_copy_matrix_d3(lattice, bravais->lattice);
  if (get_monocli_bravais(bravais->lattice, conv_sym, 2.0, primitive_axes, symprec)) {

    get_monocli_relative_axes(relative_axis, bravais->lattice, lattice, symprec);

    debug_print("is_monocli_from_P\n");
    debug_print_matrix_d3(relative_axis);

    if ((mat_Dabs(mat_Dabs(relative_axis[0][0]+relative_axis[0][1]) - 2.0) < symprec &&
	 mat_Dabs(mat_Dabs(relative_axis[1][0]+relative_axis[1][1]) - 0.0) < symprec &&
	 mat_Dabs(mat_Dabs(relative_axis[2][0]+relative_axis[2][1]) - 0.0) < symprec) ||
	(mat_Dabs(mat_Dabs(relative_axis[0][0]+relative_axis[0][1]) - 0.0) < symprec &&
	 mat_Dabs(mat_Dabs(relative_axis[1][0]+relative_axis[1][1]) - 2.0) < symprec &&
	 mat_Dabs(mat_Dabs(relative_axis[2][0]+relative_axis[2][1]) - 0.0) < symprec) ||
	(mat_Dabs(mat_Dabs(relative_axis[0][0]+relative_axis[0][1]) - 0.0) < symprec &&
	 mat_Dabs(mat_Dabs(relative_axis[1][0]+relative_axis[1][1]) - 0.0) < symprec &&
	 mat_Dabs(mat_Dabs(relative_axis[2][0]+relative_axis[2][1]) - 2.0) < symprec)) {
      bravais->centering = C_FACE;
      goto found;
    }

    if ((mat_Dabs(mat_Dabs(relative_axis[0][1]+relative_axis[0][2]) - 2.0) < symprec &&
	 mat_Dabs(mat_Dabs(relative_axis[1][1]+relative_axis[1][2]) - 0.0) < symprec &&
	 mat_Dabs(mat_Dabs(relative_axis[2][1]+relative_axis[2][2]) - 0.0) < symprec) ||
	(mat_Dabs(mat_Dabs(relative_axis[0][1]+relative_axis[0][2]) - 0.0) < symprec &&
	 mat_Dabs(mat_Dabs(relative_axis[1][1]+relative_axis[1][2]) - 2.0) < symprec &&
	 mat_Dabs(mat_Dabs(relative_axis[2][1]+relative_axis[2][2]) - 0.0) < symprec) ||
	(mat_Dabs(mat_Dabs(relative_axis[0][1]+relative_axis[0][2]) - 0.0) < symprec &&
	 mat_Dabs(mat_Dabs(relative_axis[1][1]+relative_axis[1][2]) - 0.0) < symprec &&
	 mat_Dabs(mat_Dabs(relative_axis[2][1]+relative_axis[2][2]) - 2.0) < symprec)) {
      bravais->centering = A_FACE;
      goto found;
    }

    /* otherwise bcc */
    get_monocli_bcc_to_c_center(bravais->lattice);
    bravais->centering = C_FACE;
    goto found;
  }

  return 0;

 found:
  return 1;
}

static void get_monocli_relative_axes(double relative_axis[3][3],
				      const double new_lattice[3][3],
				      const double old_lattice[3][3],
				      const double symprec)
{
  double inv_lattice[3][3];

  if (!(mat_inverse_matrix_d3(inv_lattice, old_lattice, symprec)))
    fprintf(stderr, "spglib: BUG in spglib in __LINE__, __FILE__.");

  mat_multiply_matrix_d3(relative_axis, inv_lattice, new_lattice);
}

static void get_monocli_bcc_to_c_center(double lattice[3][3])
{
  int i;
  double tmp_lattice[3][3];

  mat_copy_matrix_d3(tmp_lattice, lattice);
  for (i = 0; i < 3; i++)
    lattice[i][0] = tmp_lattice[i][0] + tmp_lattice[i][2];
}

static int get_monocli_bravais(double lattice[3][3], const Symmetry *conv_sym,
			       const double relative_vol,
			       const double monocli_axes[13][3],
			       const double symprec)
{
  int i, j, k, l, naxis[3];
  double relative_axis[3][3], volume;

  /* get 2-fold axis */
  for (i = 0; i < conv_sym->size; i++) {
    naxis[1] = get_rotation_axis(conv_sym->rot[i], 2);
    if (naxis[1] > -1) {
      break;
    }
  }

  /* get the other axes orthonormal to the 2-fold axis */
  for (i = 0; i < 13; i++) {
    if (is_monocli_orthogonal(naxis[1], i, monocli_axes, symprec)) {
      naxis[0] = i;
      for (j = 0; j < 13; j++) {
	if (is_monocli_orthogonal(naxis[1], j, monocli_axes, symprec)) {
	  naxis[2] = j;
	  for (k = 0; k < 3; k++)
	    for (l = 0; l < 3; l++)
	      relative_axis[l][k] = monocli_axes[naxis[k]][l];

	  volume = mat_Dabs(mat_get_determinant_d3(relative_axis));
	  debug_print("axes: %d %d %d, volume: %f \n", naxis[0], naxis[1], naxis[2], volume);
	  if (mat_Dabs(volume - relative_vol) < symprec) {
	    goto found;
	  }
	}
      }
    }
  }

  /* not found */
  return 0;

 found:
  for (i = 0; i < 3; i++)
    for (j = 0; j < 3; j++)
      relative_axis[j][i] = monocli_axes[naxis[i]][j];
  mat_multiply_matrix_d3(lattice, lattice, relative_axis);
  return 1;
}

static int is_monocli_orthogonal(const int b_axis, const int naxis,
				 const double monocli_axes[13][3],
				 const double symprec)
{
  if (mat_Dabs(monocli_axes[b_axis][0] * monocli_axes[naxis][0] +
	       monocli_axes[b_axis][1] * monocli_axes[naxis][1] +
	       monocli_axes[b_axis][2] * monocli_axes[naxis][2]) < symprec) {
    return 1;
  }

  return 0;
}

/* bravais is going to be modified */
static int is_holohedry(Bravais *bravais, const Cell *cell,
			const Holohedry holohedry,
			const double symprec)
{
  Bravais temp_bravais;
  double min_lattice[3][3];

  temp_bravais = *bravais;
  temp_bravais.holohedry = holohedry;
  brv_smallest_lattice_vector(min_lattice, cell->lattice, symprec);
  
  if (brv_get_brv_lattice_in_loop(&temp_bravais, min_lattice, symprec)) {
    *bravais = temp_bravais;
    return 1;
  }
  else {
    return 0;
  }
}

static int get_rotation_axis(const int rot[3][3], const int axis_num)
{
  int i, axis = -1, tmp_rot[3][3], test_rot[3][3], vec[3];
  int rot_axes[13][3] = {
    { 1, 0, 0},
    { 0, 1, 0},
    { 0, 0, 1},
    { 0, 1, 1},
    { 0, 1,-1},
    { 1, 0, 1},
    {-1, 0, 1},
    { 1, 1, 0},
    { 1,-1, 0},
    {-1, 1, 1},
    { 1,-1, 1},
    { 1, 1,-1},
    { 1, 1, 1},
  };

  int identity[3][3] = {
    { 1, 0, 0},
    { 0, 1, 0},
    { 0, 0, 1},
  };

  int inversion[3][3] = {
    {-1, 0, 0},
    { 0,-1, 0},
    { 0, 0,-1},
  };

  mat_copy_matrix_i3(test_rot, identity);
  mat_copy_matrix_i3(tmp_rot, rot);

  /* If improper, multiply inversion and get proper rotation */
  if (mat_get_determinant_i3(rot) != 1)
    mat_multiply_matrix_i3(tmp_rot, inversion, rot);

  /* Look for rotation axis */
  if (!(mat_check_identity_matrix_i3(tmp_rot, identity))) {

    for (i = 0; i < axis_num; i++) {
      mat_multiply_matrix_i3(test_rot, tmp_rot, test_rot);
      if (mat_check_identity_matrix_i3(test_rot, identity) && i < axis_num - 1)
	goto end;
    }

    if (mat_check_identity_matrix_i3(test_rot, identity)) {

      /* Look for eigenvector = rotation axis */
      for (i = 0; i < 13; i++) {
	mat_multiply_matrix_vector_i3(vec, tmp_rot, rot_axes[i]);
	
	if (vec[0] == rot_axes[i][0] &&
	    vec[1] == rot_axes[i][1] &&
	    vec[2] == rot_axes[i][2]) {

	  axis = i;

	  debug_print("--------\n");
	  debug_print("invariant axis: %d %d %d\n", rot_axes[i][0], rot_axes[i][1], rot_axes[i][2]);
	  debug_print_matrix_i3(tmp_rot);
	  break;
	}
      }
    }
  }
  
 end:
  return axis;
}
