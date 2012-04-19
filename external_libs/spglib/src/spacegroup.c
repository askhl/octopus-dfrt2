/* spacegroup.c */
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
#include "spacegroup.h"
#include "bravais.h"
#include "bravais_art.h"
#include "cell.h"
#include "debug.h"
#include "mathfunc.h"
#include "pointgroup.h"
#include "primitive.h"
#include "spacegroup_data.h"
#include "spacegroup_database.h"
#include "symmetry.h"

typedef enum {
  NO_AXIS = 0,
  PRIMARY_AXIS = 1,
  SECONDARY_AXIS = 2,
  TERTIARY_AXIS = 3
} Axis_type;



static int get_spacegroup_number(const Bravais *bravais, const Cell *primitive,
				 const Symmetry *primitive_sym, const double symprec);
static int get_rotation_class(const Bravais * bravais, const int order, const int rot[3][3],
			      const double trans[3], const double symprec);
static int get_improper_rotation_class_2axis(const Bravais * bravais, const int rot[3][3],
					     const double raw_trans[3], const double symprec);
static int get_proper_rotation_class(const Bravais * bravais, const int order, const int rot[3][3],
				     const double raw_trans[3], const double symprec);
static int get_proper_rotation_class_6axis(const int rot[3][3], const int trans_order,
					   const int nonsymmorphic);
static int get_proper_rotation_class_4axis(const Centering centering, const int rot[3][3],
					   const int trans_order, const int nonsymmorphic[3]);
static int get_proper_rotation_class_3axis(const Holohedry holohedry, const int rot[3][3],
					   const int trans_order, const int nonsymmorphic);
static int get_proper_rotation_class_2axis(const Bravais * bravais, const int rot[3][3],
					   const int trans_order);
static int get_translation_order(const Centering centering, const int order, const double raw_trans[3],
				 const int rot[3][3], const double symprec);
static int get_class_order(int rot[3][3]);
static double get_abs_diff_two_vectors(double a[3], double b[3]);
static void get_order_times_translation(double trans[3], const int order,
					const int rot[3][3],
					const double raw_trans[3]);
static int get_class_order(int rot[3][3]);
#ifdef DEBUG
static void print_rotation_class(int rot_class);
#endif

Spacegroup tbl_get_spacegroup(const Cell * cell, const double symprec)
{
  int spacegroup_number;
  Bravais bravais;
  Cell primitive;
  Holohedry holohedry;
  Symmetry symmetry;
  Spacegroup spacegroup;
  Pointgroup pointgroup;

  /* get primitive cell */
  if (sym_get_multiplicity(cell, symprec) > 1) {

    debug_print("Cell is not primitive.\n");
    primitive = prm_get_primitive(cell, symprec);
    spacegroup = tbl_get_spacegroup(&primitive, symprec);

    cel_delete_cell(&primitive);
    return spacegroup;
  }

  /* symmetry check after getting correct Bravais lattice */
  bravais = brv_get_brv_lattice(cell->lattice, symprec);
  symmetry = sym_get_operation(&bravais, cell, symprec);

  /* Get correct Bravais lattice including internal symmetry */
  holohedry = ptg_get_holohedry(bravais.holohedry, &symmetry);
  /* When holohedry from lattice doesn't correspond to holohedry   */
  /* determined by pointgroup, careful treatment is required. */
  if (holohedry < bravais.holohedry) {
    if (!(art_get_artificial_bravais(&bravais, &symmetry, cell, holohedry, symprec))) {
      spacegroup.number = 0;
      goto end;
    }
  }

  /* get space group */
  spacegroup_number = get_spacegroup_number(&bravais, cell, &symmetry, symprec);
  if (spacegroup_number) {
    spacegroup = tbl_get_spacegroup_database(spacegroup_number, 1, 1);
    pointgroup = tbl_get_pointgroup_database(spacegroup_number);
    spacegroup.pointgroup = pointgroup;
  }
  else {
    spacegroup.number = 0;
  }

 end:
  if (spacegroup.number == 0) {
    fprintf(stderr, "spglib BUG: Space group could not be found.\n");
  }
  sym_delete_symmetry(&symmetry);
  return spacegroup;
}

Symmetry tbl_get_conventional_symmetry(const Bravais *bravais,
				       const Cell *primitive,
				       const Symmetry *primitive_sym,
				       const double symprec)
{
  int i, j, k, multi, size;
  Centering centering;
  Holohedry holohedry;
  double tmp_trans;
  double coordinate_inv[3][3], tmp_matrix_d3[3][3], shift[4][3];
  double symmetry_rot_d3[3][3], primitive_sym_rot_d3[3][3];
  Symmetry symmetry;

  centering = bravais->centering;
  holohedry = bravais->holohedry;
  size = primitive_sym->size;

  if (centering == FACE) {
    symmetry = sym_new_symmetry(size * 4);
  }
  else {
    if (centering) {
      symmetry = sym_new_symmetry(size * 2);
    } else {
      symmetry = sym_new_symmetry(size);
    }
  }

  /* C^-1 = P^-1*B */
  mat_inverse_matrix_d3(tmp_matrix_d3, primitive->lattice, symprec);
  mat_multiply_matrix_d3(coordinate_inv, tmp_matrix_d3, bravais->lattice);

  for (i = 0; i < size; i++) {
    mat_cast_matrix_3i_to_3d(primitive_sym_rot_d3, primitive_sym->rot[i]);

    /* C*S*C^-1: recover conventional cell symmetry operation */
    mat_get_similar_matrix_d3(symmetry_rot_d3, primitive_sym_rot_d3,
                              coordinate_inv, symprec);

    mat_cast_matrix_3d_to_3i(symmetry.rot[i], symmetry_rot_d3);

    /* translation in conventional cell: C = B^-1*P */
    mat_inverse_matrix_d3(tmp_matrix_d3, coordinate_inv, symprec);
    mat_multiply_matrix_vector_d3(symmetry.trans[i], tmp_matrix_d3,
                                  primitive_sym->trans[i]);

    for (j = 0; j < 3; j++) {
      symmetry.trans[i][j] = symmetry.trans[i][j]
	- mat_Nint(symmetry.trans[i][j] - symprec * 2);
    }
  }


  if (centering) {
    if (centering != FACE) {
      for (i = 0; i < 3; i++)
	shift[0][i] = 0.5;

      if (centering == A_FACE)
	shift[0][0] = 0;
      if (centering == B_FACE)
	shift[0][1] = 0;
      if (centering == C_FACE)
	shift[0][2] = 0;

      multi = 2;
    }

    if (centering == FACE) {
      shift[0][0] = 0;
      shift[0][1] = 0.5;
      shift[0][2] = 0.5;
      shift[1][0] = 0.5;
      shift[1][1] = 0;
      shift[1][2] = 0.5;
      shift[2][0] = 0.5;
      shift[2][1] = 0.5;
      shift[2][2] = 0;

      multi = 4;
    }

    for (i = 0; i < multi - 1; i++) {
      for (j = 0; j < size; j++) {
	mat_copy_matrix_i3(symmetry.rot[(i+1) * size + j], symmetry.rot[j]);
	/* Translations are confined to -0 < trans < 1.0 .*/
	/* This helps reducing options (or giving definite search) */
	/* in spacegroup search. */
	for (k = 0; k < 3; k++) {
	  tmp_trans = symmetry.trans[j][k] + shift[i][k];
	  tmp_trans -= mat_Nint(tmp_trans);
	  if (tmp_trans < -symprec)
	    tmp_trans += 1.0;
	  symmetry.trans[(i+1) * size + j][k] = tmp_trans;
	}
      }
    }
  }

#ifdef DEBUG
  debug_print("*** get_conventional_symmetry ***\n");
  debug_print("--- Bravais lattice ---\n");
  debug_print_matrix_d3(bravais->lattice);
  debug_print("--- Prim lattice ---\n");
  debug_print_matrix_d3(primitive->lattice);
  debug_print("--- Coordinate inv ---\n");
  debug_print_matrix_d3(coordinate_inv);
  
  for (i = 0; i < symmetry.size; i++) {
    debug_print("--- %d ---\n", i + 1);
    debug_print_matrix_i3(symmetry.rot[i]);
    debug_print("%f %f %f\n", symmetry.trans[i][0],
		symmetry.trans[i][1], symmetry.trans[i][2]);
  }
#endif

  return symmetry;
}

static int get_spacegroup_number(const Bravais *bravais, const Cell *cell,
				 const Symmetry *prim_sym, const double symprec)
{
  int i, order, spacegroup_data;
  Symmetry symmetry;

  symmetry = tbl_get_conventional_symmetry(bravais, cell, prim_sym, symprec);

  int rot_class[symmetry.size];

  debug_print("*** get_spacegroup_number ***\n");
  for (i = 0; i < symmetry.size; i++) {

    order = get_class_order(symmetry.rot[i]);

    debug_print("--- %d ---\n", i + 1);
    debug_print("Order of symmetry: %d\n", order);

    rot_class[i] = get_rotation_class(bravais, order, symmetry.rot[i],
				      symmetry.trans[i], symprec);

    debug_print("Symmetry operation No.%d is ", i + 1);
#ifdef DEBUG
    print_rotation_class(rot_class[i]);
#endif

    /* Error to find class */
    if (rot_class[i] == 0)
      return 0;
  }

  debug_print("holohedry: %d\n", bravais->holohedry);
  debug_print("centering: %d\n", bravais->centering);

  spacegroup_data = tbl_get_spacegroup_data(&symmetry, bravais, rot_class, symprec);

  sym_delete_symmetry(&symmetry);
  
  return spacegroup_data;
}

static int get_class_order(int rot[3][3])
{
  int i, order = 0, test_rot[3][3];
  int identity[3][3] = {
    {1, 0, 0},
    {0, 1, 0},
    {0, 0, 1}
  };
  int inversion[3][3] = {
    {-1, 0, 0},
    {0, -1, 0},
    {0, 0, -1}
  };


  mat_copy_matrix_i3(test_rot, identity);

  for (i = 0; i < 6; i++) {

    mat_multiply_matrix_i3(test_rot, rot, test_rot);

    if (mat_check_identity_matrix_i3(test_rot, identity) ||
	(mat_check_identity_matrix_i3(test_rot, inversion))) {

      order = i + 1;
      break;
    }
  }

  if (!order || order > 6) {
    fprintf(stderr, "spglib BUG: Invalid rotation matrix\n");
    return 0;
  }

  return order;
}

static int get_rotation_class(const Bravais * bravais, const int order,
			      const int rot[3][3], const double trans[3],
			      const double symprec)
{
  int i;
  double sum;

  if (mat_get_determinant_i3(rot) == 1) {

    if (order == 1) {

      sum = 0;
      for (i = 0; i < 3; i++)
	sum = sum + mat_Dabs(trans[i] - mat_Nint(trans[i]));

      if (sum < symprec * 3) {
	return IDENTITY;
      }
      else {
	return PURE_TRANSLATION;
      }

    }
    else {
      return get_proper_rotation_class(bravais, order, rot, trans, symprec);
    }
  }

  if (mat_get_determinant_i3(rot) == -1) {

    if (order == 1)
      return INVERSION;

    if (order == 2)
      return get_improper_rotation_class_2axis(bravais, rot, trans,
					       symprec);

    if (order == 3)
      return PRIMARY_M3_AXIS;

    if (order == 4)
      return PRIMARY_M4_AXIS;

    if (order == 6)
      return PRIMARY_M6_AXIS;
  }

  return 0;
}

static int get_improper_rotation_class_2axis(const Bravais * bravais,
					     const int rot[3][3],
					     const double raw_trans[3],
					     const double symprec)
{
  int i, sum=0;
  Holohedry holohedry;
  Centering centering;
  double trans[3], trans2[3];

  int mirror_x[3][3] = {
    {-1, 0, 0},
    {0, 1, 0},
    {0, 0, 1}
  };
  int mirror_y[3][3] = {
    {1, 0, 0},
    {0, -1, 0},
    {0, 0, 1}
  };
  int mirror_z[3][3] = {
    {1, 0, 0},
    {0, 1, 0},
    {0, 0, -1}
  };

  double halves[7][3] = {
    {0.5, 0, 0},            /* 0 */
    {0, 0.5, 0},            /* 1 */
    {0, 0, 0.5},            /* 2 */
    {0, 0.5, 0.5},          /* 3 */
    {0.5, 0, 0.5},          /* 4 */
    {0.5, 0.5, 0},          /* 5 */
    {0.5, 0.5, 0.5},        /* 6 */
  };

  Axis_type axis_type = NO_AXIS;

  sum = rot[0][0] + rot[0][1] + rot[0][2] +
    rot[1][0] + rot[1][1] + rot[1][2] +
    rot[2][0] + rot[2][1] + rot[2][2];

  holohedry = bravais->holohedry;
  centering = bravais->centering;

  /* determine axis_type */
  if (sum == 1){

    if (mat_check_identity_matrix_i3(rot, mirror_z)) {
      axis_type = PRIMARY_AXIS;
    }
    else {
      if (holohedry == MONOCLI || holohedry == ORTHO || holohedry == CUBIC)
	axis_type = PRIMARY_AXIS;

      if (holohedry == TETRA || holohedry == HEXA || holohedry == TRIGO)
	axis_type = SECONDARY_AXIS;
    }

  }

  if (sum == 2)
    axis_type = SECONDARY_AXIS;

  if (sum == 3 || sum == 0) {
    axis_type = TERTIARY_AXIS;
  }

  if (sum == -1) {

    if (holohedry == TETRA || holohedry == CUBIC)
      axis_type = TERTIARY_AXIS;

    if (holohedry == HEXA || holohedry == TRIGO)
      axis_type = SECONDARY_AXIS;
  }


  /* operate twice from origin and extract one translation */
  mat_multiply_matrix_vector_id3(trans, rot, raw_trans);

  for (i = 0; i < 3; i++) {
    trans[i] = (trans[i] + raw_trans[i]) * 0.5;
    trans[i] = trans[i] - mat_Nint(trans[i] - symprec);
  }


  debug_print("Improper 2-axis\n");
  debug_print_matrix_i3(rot);
  debug_print("raw: %f %f %f\n", raw_trans[0], raw_trans[1], raw_trans[2]);
  debug_print("new: %f %f %f\n", trans[0], trans[1], trans[2]);
  debug_print("holohedry: %d\n", holohedry);
  debug_print("centering: %d\n", centering);
  debug_print("axis: %d\n", axis_type);


  if (mat_Dabs(trans[0]) + mat_Dabs(trans[1]) + mat_Dabs(trans[2]) < symprec*3 &&
      (!(holohedry == HEXA || holohedry == TRIGO)))
    return M_PLANE;

  if (holohedry == TETRA && centering == NO_CENTER) {

    if (axis_type == PRIMARY_AXIS)
      return PRIMARY_N;

    if (axis_type == SECONDARY_AXIS) {

      if (get_abs_diff_two_vectors(trans, halves[0]) < symprec*3 ||
	  get_abs_diff_two_vectors(trans, halves[1]) < symprec*3)
	return SECONDARY_A_B;

      if (get_abs_diff_two_vectors(trans, halves[2]) < symprec*3)
	return SECONDARY_C;

      return SECONDARY_N;
    }

    if (axis_type == TERTIARY_AXIS) {

      if (mat_Dabs(trans[2]) < symprec)
	return TERTIARY_M;

      if (mat_Dabs(trans[2] - 0.5) < symprec)
	return TERTIARY_C;
    }
  }

  if (holohedry == TETRA && centering == BODY) {

    if (axis_type == PRIMARY_AXIS || axis_type == SECONDARY_AXIS) {

      if (get_abs_diff_two_vectors(trans, halves[0]) < symprec*3 ||
	  get_abs_diff_two_vectors(trans, halves[1]) < symprec*3 ||
	  get_abs_diff_two_vectors(trans, halves[2]) < symprec*3)
	return A_B_C_PLANE;

      if (get_abs_diff_two_vectors(trans, halves[3]) < symprec*3 ||
	  get_abs_diff_two_vectors(trans, halves[4]) < symprec*3 ||
	  get_abs_diff_two_vectors(trans, halves[5]) < symprec*3)
	return M_PLANE;
    }

    if (axis_type == TERTIARY_AXIS) {

      if (mat_Dabs(trans[2]) < symprec || mat_Dabs(trans[2] - 0.5) < symprec)
	return TERTIARY_M;

      return TERTIARY_D;

    }
  }

  if (holohedry == RHOMB) {

    if (mat_Dabs(mat_Dabs(trans[0]) + mat_Dabs(trans[1]) + mat_Dabs(trans[2]) - 1)
	< symprec*3)
      return SECONDARY_M;

    if (mat_Dabs(trans[0]) + mat_Dabs(trans[1]) + mat_Dabs(trans[2]) - 0.5 <
	symprec*3
	|| mat_Dabs(trans[0]) + mat_Dabs(trans[1]) + mat_Dabs(trans[2]) - 1.5 <
	symprec*3)
      return SECONDARY_C_RHOMB;
  }

  if (holohedry == HEXA || holohedry == TRIGO) {

    if (axis_type == PRIMARY_AXIS)

      if (mat_Dabs(trans[2]) < symprec)
	return PRIMARY_M;

    if (axis_type == SECONDARY_AXIS) {

      if (mat_Dabs(trans[2]) < symprec)
	return SECONDARY_M;

      if (mat_Dabs(trans[2] - 0.5) < symprec)
	return SECONDARY_C_HEXA;
    }

    if (axis_type == TERTIARY_AXIS) {

      if (mat_Dabs(trans[2]) < symprec) 
	return TERTIARY_M_HEXA;

      if (mat_Dabs(trans[2]-0.5) < symprec)
	return TERTIARY_C_HEXA;

    }                
  }

  if (holohedry == CUBIC && centering == NO_CENTER) {

    if (axis_type == PRIMARY_AXIS)
      return PRIMARY_N;

    if (mat_Dabs(mat_Dabs(trans[0]) + mat_Dabs(trans[1]) + mat_Dabs(trans[2]) - 0.5)
	< symprec*3
	|| mat_Dabs(mat_Dabs(trans[0]) + mat_Dabs(trans[1]) + mat_Dabs(trans[2]) - 1.5)
	< symprec*3)
      return TERTIARY_N;

    if (mat_Dabs(mat_Dabs(trans[0]) + mat_Dabs(trans[1]) + mat_Dabs(trans[2]) - 1)
	< symprec*3)
      return TERTIARY_M;
  }

  /* other cases */
  if (get_abs_diff_two_vectors(trans, halves[0]) < symprec*3 ||
      get_abs_diff_two_vectors(trans, halves[1]) < symprec*3 ||
      get_abs_diff_two_vectors(trans, halves[2]) < symprec*3) {
    return A_B_C_PLANE;
  }

  if ((axis_type == PRIMARY_AXIS || axis_type == SECONDARY_AXIS) &&
      (get_abs_diff_two_vectors(trans, halves[3]) < symprec*3 ||
       get_abs_diff_two_vectors(trans, halves[4]) < symprec*3 ||
       get_abs_diff_two_vectors(trans, halves[5]) < symprec*3)) {
    return N_PLANE;
  }

  if (axis_type == TERTIARY_AXIS &&
      get_abs_diff_two_vectors(trans, halves[6]) < symprec*3) {
    return N_PLANE;
  }

  for (i = 0; i < 3; i++)
    trans2[i] = trans[i] * 2;

  if ((mat_check_identity_matrix_i3(rot, mirror_x) &&
       get_abs_diff_two_vectors(trans2, halves[3]) < symprec*3) ||
      (mat_check_identity_matrix_i3(rot, mirror_y) &&
       get_abs_diff_two_vectors(trans2, halves[4]) < symprec*3) ||
      (mat_check_identity_matrix_i3(rot, mirror_z) &&
       get_abs_diff_two_vectors(trans2, halves[5]) < symprec*3)) {
    return D_PLANE;
  }

  if (axis_type == TERTIARY_AXIS &&
      get_abs_diff_two_vectors(trans2, halves[6]) < symprec*3) {
    return D_PLANE;
  }

  return G_PLANE;
}

static double get_abs_diff_two_vectors(double a[3], double b[3])
{
  int i;
  double sum = 0;

  for (i = 0; i < 3; i++)
    sum = sum + mat_Dabs(mat_Dabs(a[i]) - mat_Dabs(b[i]));

  return sum;
}

static int get_proper_rotation_class(const Bravais * bravais, const int order,
				     const int rot[3][3], const double raw_trans[3],
				     const double symprec)
{
  int i, rot_class, trans_order, nonsymmorphic[3];
  Holohedry holohedry;
  Centering centering;
  double trans[3];

  centering = bravais->centering;
  holohedry = bravais->holohedry;

  trans_order =
    get_translation_order(centering, order, raw_trans, rot, symprec);

  get_order_times_translation(trans, order, rot, raw_trans);
  /* trans[i] should be integer unless centering */
  for (i = 0; i < 3; i++)
    nonsymmorphic[i] = mat_Nint(trans[i]);

  if (order == 2)
    rot_class =
      get_proper_rotation_class_2axis(bravais, rot, trans_order);

  if (order == 3)
    rot_class =
      get_proper_rotation_class_3axis(bravais->holohedry, rot,
				      trans_order, nonsymmorphic[2]);

  if (order == 4)
    rot_class =
      get_proper_rotation_class_4axis(bravais->centering, rot,
				      trans_order, nonsymmorphic);

  if (order == 6)
    rot_class = get_proper_rotation_class_6axis(rot, trans_order,
						nonsymmorphic[2]);


  debug_print("Order of translation: %d\n", trans_order);
  debug_print("Nonsymmorphic translation*order of symmetry: %d %d %d\n",
	      nonsymmorphic[0], nonsymmorphic[1], nonsymmorphic[2]);
  debug_print("Nonsymmorphic translation*order of symmetry: %f %f %f\n",
	      trans[0], trans[1], trans[2]);

  return rot_class;
}

int get_proper_rotation_class_6axis(const int rot[3][3], const int trans_order,
                                    const int nonsymmorphic)
{
  if (trans_order == 1)
    return PRIMARY_6_AXIS;

  if (trans_order == 2)
    return PRIMARY_6_3_AXIS;

  if (trans_order == 3) {

    /* +pi/3 */
    if (rot[0][0] == 1) {

      if (nonsymmorphic == 2)
	return PRIMARY_6_2_AXIS;

      if (nonsymmorphic == -2)
	return PRIMARY_6_4_AXIS;
    }

    /* -pi/3 */
    if (rot[0][0] == 0) {

      if (nonsymmorphic == 2)
	return PRIMARY_6_4_AXIS;

      if (nonsymmorphic == -2)
	return PRIMARY_6_2_AXIS;
    }
  }

  /* +pi/3 */
  if (rot[0][0] == 1) {

    if (nonsymmorphic == 1)
      return PRIMARY_6_1_AXIS;

    if (nonsymmorphic == -1)
      return PRIMARY_6_5_AXIS;
  }

  /* -pi/3 */
  if (rot[0][0] == 0) {

    if (nonsymmorphic == 1)
      return PRIMARY_6_5_AXIS;

    if (nonsymmorphic == -1)
      return PRIMARY_6_1_AXIS;
  }

  debug_print("Invalid 6-axis\n");
  debug_print("[0][0]: %d, trans_order: %d\n", rot[0][0], trans_order);
  debug_print("nonsymmorphic: %d\n", nonsymmorphic);
  debug_print_matrix_i3(rot);
  fprintf(stderr, "spglib BUG: Invalid 6-axis.\n");
  return 0;
}

int get_proper_rotation_class_4axis(const Centering centering, const  int rot[3][3],
                                    const int trans_order,
				    const int nonsymmorphic[3])
{
  int i;

  if (trans_order == 1)
    return PRIMARY_4_AXIS;

  if (trans_order == 2)
    return PRIMARY_4_2_AXIS;

  if (centering)
    return PRIMARY_4_1_4_3_AXIS;

  for (i = 0; i < 3; i++)

    if (rot[i][i] == 1) {

      /* +pi/2 rotation */
      if ((i == 0 && rot[1][2] == -1) ||
	  (i == 1 && rot[2][0] == -1) ||
	  (i == 2 && rot[0][1] == -1)) {

	if (nonsymmorphic[i] == 1)
	  return PRIMARY_4_1_AXIS;

	if (nonsymmorphic[i] == -1)
	  return PRIMARY_4_3_AXIS;
      }

      /* -pi/2 rotation */
      if ((i == 0 && rot[1][2] == 1) ||
	  (i == 1 && rot[2][0] == 1) ||
	  (i == 2 && rot[0][1] == 1)) {

	if (nonsymmorphic[i] == 1)
	  return PRIMARY_4_3_AXIS;

	if (nonsymmorphic[i] == -1)
	  return PRIMARY_4_1_AXIS;
      }
    }

  fprintf(stderr, "spglib BUG: Invalid 4-axis.\n");
  return 0;
}

static int get_proper_rotation_class_3axis(const Holohedry holohedry,
					   const int rot[3][3],
					   const int trans_order,
					   const int nonsymmorphic)
{
  if (trans_order == 1)
    return PRIMARY_3_AXIS;

  if (holohedry == RHOMB || holohedry == CUBIC)
    return PRIMARY_3_3_1_3_2_AXIS;

  /* +2pi/3 rotation */
  if (!rot[0][0]) {

    if (nonsymmorphic == 1)
      return PRIMARY_3_1_AXIS;

    if (nonsymmorphic == -1)
      return PRIMARY_3_2_AXIS;
  }

  /* -2pi/3 rotation */
  if (rot[0][0] == -1) {

    if (nonsymmorphic == 1)
      return PRIMARY_3_2_AXIS;

    if (nonsymmorphic == -1)
      return PRIMARY_3_1_AXIS;
  }

  fprintf(stderr, "spglib BUG: Invalid 3-axis.\n");
  return 0;

}

static int get_proper_rotation_class_2axis(const Bravais * bravais,
					   const int rot[3][3],
					   const int trans_order)
{
  Holohedry holohedry;
  Centering centering;
  int sum;
  Axis_type axis_type = NO_AXIS;

  centering = bravais->centering;
  holohedry = bravais->holohedry;
  axis_type = PRIMARY_AXIS;

  /* determine axis type */
  if (holohedry == TETRA || holohedry == CUBIC)

    if (abs(rot[0][0]) + abs(rot[1][1]) + abs(rot[2][2]) == 1)
      axis_type = TERTIARY_AXIS;

  if (holohedry == HEXA || holohedry == TRIGO) {

    sum = rot[0][0] + rot[0][1] + rot[0][2] +
      rot[1][0] + rot[1][1] + rot[1][2] +
      rot[2][0] + rot[2][1] + rot[2][2];

    if (sum != -1)
      axis_type = SECONDARY_AXIS;

    if (sum == 0 || sum == -3)
      axis_type = TERTIARY_AXIS;
  }

  debug_print("axis type: %d\n", axis_type);

  /* get rotation class */
  if (axis_type == SECONDARY_AXIS)
    return SECONDARY_2_AXIS;

  if (axis_type == TERTIARY_AXIS && holohedry == TETRA)
    return TERTIARY_2_AXIS;

  if (axis_type == TERTIARY_AXIS && centering == NO_CENTER &&
      (holohedry == HEXA || holohedry == TRIGO || holohedry == CUBIC))
    return TERTIARY_2_AXIS;

  if (trans_order == 1 ||
      (holohedry == TETRA && centering == BODY) || holohedry == RHOMB)
    return PRIMARY_2_AXIS;

  /* other */
  return PRIMARY_2_1_AXIS;
}


static int get_translation_order(const Centering centering, const int order,
				 const double raw_trans[3], const int rot[3][3],
				 const double symprec)
{
  int i, j, k;
  double sum, reduced_trans[3], trans[3];
  double face_center[3][3] = {
    {0, 0.5, 0.5},
    {0.5, 0, 0.5},
    {0.5, 0.5, 0}
  };

  get_order_times_translation(trans, order, rot, raw_trans);

  for (i = 0; i < 3; i++)
    trans[i] = trans[i] / order - mat_Nint(trans[i] / order - symprec);

  for (i = 1; i < order + 1; i++) {

    for (j = 0; j < 3; j++) {
      reduced_trans[j] =
	trans[j] * i - mat_Nint(trans[j] * i - symprec);
    }

    debug_print_matrix_i3(rot);
    debug_print("order: %d\n", order);
    debug_print("raw_trans: %f %f %f\n",
		raw_trans[0], raw_trans[1], raw_trans[2]);
    debug_print("order times trans: %f %f %f\n",
		trans[0], trans[1], trans[2]);
    debug_print("reduced order times trans: %f %f %f\n",
		reduced_trans[0], reduced_trans[1], reduced_trans[2]);

    /* Origine */
    sum = 0;
    for (j = 0; j < 3; j++)
      sum = sum + mat_Dabs(reduced_trans[j]);

    if (sum < symprec * 3)
      return i;

    /* A-, B-, C-face or face-centered */
    for (j = 0; j < 3; j++) {

      sum = 0;
      for (k = 0; k < 3; k++)
	sum = sum + mat_Dabs(reduced_trans[k] - face_center[j][k]);

      if ((centering == j + 1 || centering == FACE) && sum < symprec * 3)
	return i;
    }

    /* Body-centered */
    sum = 0;
    for (j = 0; j < 3; j++)
      sum = sum + mat_Dabs(reduced_trans[j] - 0.5);

    if (centering == BODY && sum < symprec * 3)
      return i;

  }

  fprintf(stderr, "spglib BUG: Could not find translation order.");
  return 0;
}

/* measure translation after 'symmetry order' times symmetry operations */
static void get_order_times_translation(double trans[3], const int order,
					const int rot[3][3],
					const double raw_trans[3])
{
  int i, j;

  for (i = 0; i < 3; i++)
    trans[i] = 0;
    
  for (i = 0; i < order; i++) {

    mat_multiply_matrix_vector_id3(trans, rot, trans);

    for (j = 0; j < 3; j++)
      trans[j] = trans[j] + raw_trans[j];
  }

}

#ifdef DEBUG
static void print_rotation_class(int rot_class)
{

  switch (rot_class) {

  case IDENTITY:
    printf("identity\n");
    break;
  case PURE_TRANSLATION:
    printf("pure translation\n");
    break;
  case SECONDARY_2_AXIS:
    printf("secondary 2 axis\n");
    break;
  case TERTIARY_2_AXIS:
    printf("tertiary 2 axis\n");
    break;
  case PRIMARY_2_AXIS:
    printf("2 axis\n");
    break;
  case PRIMARY_2_1_AXIS:
    printf("2_1 axis\n");
    break;
  case PRIMARY_3_AXIS:
    printf("3 axis\n");
    break;
  case PRIMARY_3_3_1_3_2_AXIS:
    printf("3, 3_1 or 3_2 axis\n");
    break;
  case PRIMARY_3_1_AXIS:
    printf("3_1 axis\n");
    break;
  case PRIMARY_3_2_AXIS:
    printf("3_2 axis\n");
    break;
  case PRIMARY_4_AXIS:
    printf("4 axis\n");
    break;
  case PRIMARY_4_2_AXIS:
    printf("4_2 axis\n");
    break;
  case PRIMARY_4_1_4_3_AXIS:
    printf("4_1 or 4_3 axis\n");
    break;
  case PRIMARY_4_1_AXIS:
    printf("4_1 axis\n");
    break;
  case PRIMARY_4_3_AXIS:
    printf("4_3 axis\n");
    break;
  case PRIMARY_6_AXIS:
    printf("6 axis\n");
    break;
  case PRIMARY_6_3_AXIS:
    printf("6_3 axis\n");
    break;
  case PRIMARY_6_2_AXIS:
    printf("6_2 axis\n");
    break;
  case PRIMARY_6_4_AXIS:
    printf("6_4 axis\n");
    break;
  case PRIMARY_6_1_AXIS:
    printf("6_1 axis\n");
    break;
  case PRIMARY_6_5_AXIS:
    printf("6_5 axis\n");
    break;
  case INVERSION:
    printf("inversion\n");
    break;
  case PRIMARY_M3_AXIS:
    printf("-3 axis\n");
    break;
  case PRIMARY_M4_AXIS:
    printf("-4 axis\n");
    break;
  case PRIMARY_M6_AXIS:
    printf("-6 axis\n");
    break;
  case M_PLANE:
    printf("mirror plane\n");
    break;
  case PRIMARY_N:
    printf("primary n plane\n");
    break;
  case SECONDARY_A_B:
    printf("secondary a or b plane\n");
    break;
  case SECONDARY_C:
    printf("secondary c plane\n");
    break;
  case SECONDARY_N:
    printf("secondary n plane\n");
    break;
  case TERTIARY_M:
    printf("tertiary m plane\n");
    break;
  case TERTIARY_C:
    printf("tertiary c plane\n");
    break;
  case A_B_C_PLANE:
    printf("a, b or c plane\n");
    break;
  case TERTIARY_D:
    printf("tertiary m plane\n");
    break;
  case SECONDARY_M:
    printf("secondary m plane\n");
    break;
  case SECONDARY_C_RHOMB:
    printf("secondary c plane\n");
    break;
  case PRIMARY_M:
    printf("primary m plane\n");
    break;
  case TERTIARY_M_HEXA:
    printf("tertiary m plane\n");
    break;
  case TERTIARY_C_HEXA:
    printf("tertiary c plane\n");
    break;
  case SECONDARY_C_HEXA:
    printf("secondary c plane\n");
    break;
  case TERTIARY_N:
    printf("tertiary n plane\n");
    break;
  case N_PLANE:
    printf("n plane\n");
    break;
  case D_PLANE:
    printf("d plane\n");
    break;
  case G_PLANE:
    printf("g plane\n");
    break;
  default:
    printf("bug \n");
    break;
  }
}
#endif
