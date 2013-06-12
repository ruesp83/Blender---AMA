/*
 * ***** BEGIN GPL LICENSE BLOCK *****
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software  Foundation,
 * Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 *
 * The Original Code is Copyright (C) 2005 by the Blender Foundation.
 * All rights reserved.
 *
 * Contributor(s): Daniel Dunbar
 *                 Ton Roosendaal,
 *                 Ben Batt,
 *                 Brecht Van Lommel,
 *                 Campbell Barton
 *
 * ***** END GPL LICENSE BLOCK *****
 *
 */

/** \file blender/modifiers/intern/MOD_array.c
 *  \ingroup modifiers
 */


/* Array modifier: duplicates the object multiple times along an axis */

#include "MEM_guardedalloc.h"

#include "BLI_math.h"
#include "BLI_utildefines.h"
#include "BLI_string.h"
#include "BLI_ghash.h"
#include "BLI_edgehash.h"

#include "DNA_curve_types.h"
#include "DNA_group_types.h"
#include "DNA_meshdata_types.h"
#include "DNA_object_types.h"
#include "DNA_scene_types.h"

#include "BKE_ama.h"
#include "BKE_anim.h"
#include "BKE_cdderivedmesh.h"
#include "BKE_displist.h"
#include "BKE_mesh.h"
#include "BKE_modifier.h"
#include "BKE_object.h"

#include "bmesh.h"

#include "depsgraph_private.h"

#include <ctype.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

static void initData(ModifierData *md)
{
	ArrayModifierData *amd = (ArrayModifierData *) md;

	/* default to 2 duplicates distributed along the x-axis by an
	 * offset of 1 object-width
	 */
	amd->start_cap = amd->mid_cap = amd->end_cap = amd->curve_cap = amd->curve_ob = amd->offset_ob = NULL;
	amd->count = 2;
	zero_v3(amd->offset);
	amd->scale[0] = 1;
	amd->scale[1] = amd->scale[2] = 0;
	amd->length = 0;
	amd->merge_dist = 0.01;
	amd->fit_type = MOD_ARR_FIXEDCOUNT;
	amd->offset_type = MOD_ARR_OFF_RELATIVE;
	amd->flags = 0;
	
	amd->seed[0] = amd->seed[1] = amd->seed[2] = 0;
	amd->type = MOD_ARR_MOD_NRM;
	amd->mode = !MOD_ARR_MOD_ADV;
	amd->displays = !MOD_ARR_DIS_ADV;
	zero_v3(amd->loc_offset);
	zero_v3(amd->rot_offset);
	amd->scale_offset[0] = amd->scale_offset[1] = amd->scale_offset[2] = 1;
	amd->sign = MOD_ARR_SIGN_P + MOD_ARR_SIGN_L;
	amd->lock = MOD_ARR_LOCK;
	amd->flag_offset = MOD_ARR_PROP + MOD_ARR_LOCAL;
	amd->rays = 1;
	amd->rand_mat = MOD_ARR_MAT;
	amd->mat_ob = MOD_ARR_AR_MAT_RND;
	amd->cont_mat = 1;
	amd->count_mc = 1;
	amd->dist_mc = MOD_ARR_DIST_SEQ;
	amd->rays_dir = MOD_ARR_RAYS_X;
	amd->arr_group = NULL;
	amd->rand_group = !MOD_ARR_RAND_GROUP;
	amd->dist_cu = MOD_ARR_DIST_EVENLY;
	amd->outer_cp = !MOD_ARR_CP_FIRST;
	amd->Mem_Mat_Ob.start_cap = 0;
	amd->Mem_Mat_Ob.end_cap = 0;

	zero_v3(amd->tmp.loc);
	zero_v3(amd->tmp.rot);
	zero_v3(amd->tmp.scale);
	amd->tmp.seed[0] = amd->tmp.seed[1] = amd->tmp.seed[2] = 0;
}

static void copyData(ModifierData *md, ModifierData *target)
{
	ArrayModifierData *amd = (ArrayModifierData *) md;
	ArrayModifierData *tamd = (ArrayModifierData *) target;

	tamd->start_cap = amd->start_cap;
	tamd->mid_cap = amd->mid_cap;
	tamd->end_cap = amd->end_cap;
	tamd->curve_cap = amd->curve_cap;
	tamd->curve_ob = amd->curve_ob;
	tamd->offset_ob = amd->offset_ob;
	tamd->count = amd->count;
	copy_v3_v3(tamd->offset, amd->offset);
	copy_v3_v3(tamd->scale, amd->scale);
	tamd->length = amd->length;
	tamd->merge_dist = amd->merge_dist;
	tamd->fit_type = amd->fit_type;
	tamd->offset_type = amd->offset_type;
	tamd->flags = amd->flags;
	
	copy_v3_v3_int(tamd->seed, amd->seed);
	tamd->type = amd->type;
	tamd->mode = amd->mode;
	tamd->displays = amd->displays;
	copy_v3_v3(tamd->loc_offset, amd->loc_offset);
	copy_v3_v3(tamd->rot_offset, amd->rot_offset);
	copy_v3_v3(tamd->scale_offset, amd->scale_offset);
	tamd->sign = amd->sign;	
	tamd->lock = amd->lock;
	tamd->flag_offset = amd->flag_offset;
	tamd->rays = amd->rays;
	tamd->Mem_Ob = MEM_dupallocN(amd->Mem_Ob);
	tamd->Mem_Mat_Ob.start_cap = amd->Mem_Mat_Ob.start_cap;
	tamd->Mem_Mat_Ob.end_cap = amd->Mem_Mat_Ob.end_cap;
	tamd->Mem_Mat_Ob.mid_cap = MEM_dupallocN(amd->Mem_Mat_Ob.mid_cap);
	tamd->rand_mat = amd->rand_mat;
	tamd->mat_ob = amd->mat_ob;
	tamd->cont_mat = amd->cont_mat;
	tamd->count_mc = amd->count_mc;
	tamd->dist_mc = amd->dist_mc;
	tamd->rays_dir = amd->rays_dir;
	tamd->arr_group = amd->arr_group;
	tamd->rand_group = amd->rand_group;
	tamd->dist_cu = amd->dist_cu;
	tamd->outer_cp = amd->outer_cp;
}

static void foreachObjectLink(
        ModifierData *md, Object *ob,
        void (*walk)(void *userData, Object *ob, Object **obpoin),
        void *userData)
{
	ArrayModifierData *amd = (ArrayModifierData *) md;

	walk(userData, ob, &amd->start_cap);
	walk(userData, ob, &amd->mid_cap);
	walk(userData, ob, &amd->end_cap);
	walk(userData, ob, &amd->curve_cap);
	walk(userData, ob, &amd->curve_ob);
	walk(userData, ob, &amd->offset_ob);
}

static void updateDepgraph(ModifierData *md, DagForest *forest,
                           struct Scene *UNUSED(scene), Object *UNUSED(ob), DagNode *obNode)
{
	ArrayModifierData *amd = (ArrayModifierData *) md;

	if (amd->start_cap) {
		DagNode *curNode = dag_get_node(forest, amd->start_cap);

		dag_add_relation(forest, curNode, obNode,
		                 DAG_RL_DATA_DATA | DAG_RL_OB_DATA, "Array Modifier");
	}
	if (amd->mid_cap) {
		DagNode *curNode = dag_get_node(forest, amd->mid_cap);

		dag_add_relation(forest, curNode, obNode,
						 DAG_RL_DATA_DATA | DAG_RL_OB_DATA, "Array Modifier");
	}
	if (amd->end_cap) {
		DagNode *curNode = dag_get_node(forest, amd->end_cap);

		dag_add_relation(forest, curNode, obNode,
		                 DAG_RL_DATA_DATA | DAG_RL_OB_DATA, "Array Modifier");
	}
	if (amd->curve_cap) {
		DagNode *curNode = dag_get_node(forest, amd->curve_cap);

		dag_add_relation(forest, curNode, obNode,
		                 DAG_RL_DATA_DATA | DAG_RL_OB_DATA, "Array Modifier");
	}
	if (amd->curve_ob) {
		DagNode *curNode = dag_get_node(forest, amd->curve_ob);

		dag_add_relation(forest, curNode, obNode,
		                 DAG_RL_DATA_DATA | DAG_RL_OB_DATA, "Array Modifier");
	}
	if (amd->offset_ob) {
		DagNode *curNode = dag_get_node(forest, amd->offset_ob);

		dag_add_relation(forest, curNode, obNode,
		                 DAG_RL_DATA_DATA | DAG_RL_OB_DATA, "Array Modifier");
	}
}

static int *find_doubles_index_map(BMesh *bm, BMOperator *dupe_op,
                                   const ArrayModifierData *amd,
                                   int *index_map_length)
{
	BMOperator find_op;
	BMOIter oiter;
	BMVert *v, *v2;
	BMElem *ele;
	int *index_map, i;

	BMO_op_initf(bm, &find_op, (BMO_FLAG_DEFAULTS & ~BMO_FLAG_RESPECT_HIDE),
	             "find_doubles verts=%av dist=%f keep_verts=%s",
	             amd->merge_dist, dupe_op, "geom");

	BMO_op_exec(bm, &find_op);

	i = 0;
	BMO_ITER (ele, &oiter, dupe_op->slots_in, "geom", BM_ALL) {
		BM_elem_index_set(ele, i); /* set_dirty */
		i++;
	}

	BMO_ITER (ele, &oiter, dupe_op->slots_out, "geom.out", BM_ALL) {
		BM_elem_index_set(ele, i); /* set_dirty */
		i++;
	}
	/* above loops over all, so set all to dirty, if this is somehow
	 * setting valid values, this line can be removed - campbell */
	bm->elem_index_dirty |= BM_VERT | BM_EDGE | BM_FACE;

	(*index_map_length) = i;
	index_map = MEM_callocN(sizeof(int) * (*index_map_length), "index_map");

	/*element type argument doesn't do anything here*/
	BMO_ITER (v, &oiter, find_op.slots_out, "targetmap.out", 0) {
		v2 = BMO_iter_map_value_p(&oiter);

		index_map[BM_elem_index_get(v)] = BM_elem_index_get(v2) + 1;
	}

	BMO_op_finish(bm, &find_op);

	return index_map;
}

/* Used for start/mid/end cap.
 *
 * this function expects all existing vertices to be tagged,
 * so we can know new verts are not tagged.
 *
 * All verts will be tagged on exit.
 */
static void bm_merge_dm_transform(BMesh *bm, DerivedMesh *dm, float mat[4][4],
                                  const ArrayModifierData *amd,
                                  BMOperator *dupe_op,
                                  BMOpSlot dupe_op_slot_args[BMO_OP_MAX_SLOTS], const char *dupe_slot_name,
                                  BMOperator *weld_op)
{
	const int is_input = (dupe_op->slots_in == dupe_op_slot_args);
	BMVert *v, *v2, *v3;
	BMIter iter;

	/* Add the DerivedMesh's elements to the BMesh. The pre-existing
	 * elements were already tagged, so the new elements can be
	 * identified by not having the BM_ELEM_TAG flag set. */
	DM_to_bmesh_ex(dm, bm, false);

	if (amd->flags & MOD_ARR_MERGE) {
		/* if merging is enabled, find doubles */
		
		BMOIter oiter;
		BMOperator find_op;
		BMOpSlot *slot_targetmap;

		BMO_op_initf(bm, &find_op, (BMO_FLAG_DEFAULTS & ~BMO_FLAG_RESPECT_HIDE),
		             is_input ?  /* ugh */
		             "find_doubles verts=%Hv dist=%f keep_verts=%s" :
		             "find_doubles verts=%Hv dist=%f keep_verts=%S",
		             BM_ELEM_TAG, amd->merge_dist,
		             dupe_op, dupe_slot_name);

		/* append the dupe's geom to the findop input verts */
		if (is_input) {
			BMO_slot_buffer_append(&find_op, slots_in, "verts",
			                       dupe_op,  slots_in, dupe_slot_name);
		}
		else if (dupe_op->slots_out == dupe_op_slot_args) {
			BMO_slot_buffer_append(&find_op, slots_in,  "verts",
			                       dupe_op,  slots_out, dupe_slot_name);
		}
		else {
			BLI_assert(0);
		}

		/* transform and tag verts */
		BM_ITER_MESH (v, &iter, bm, BM_VERTS_OF_MESH) {
			if (!BM_elem_flag_test(v, BM_ELEM_TAG)) {
				mul_m4_v3(mat, v->co);
				BM_elem_flag_enable(v, BM_ELEM_TAG);
			}
		}

		BMO_op_exec(bm, &find_op);

		slot_targetmap = BMO_slot_get(weld_op->slots_in, "targetmap");

		/* add new merge targets to weld operator */
		BMO_ITER (v, &oiter, find_op.slots_out, "targetmap.out", 0) {
			v2 = BMO_iter_map_value_p(&oiter);
			/* check in case the target vertex (v2) is already marked
			 * for merging */
			while ((v3 = BMO_slot_map_elem_get(slot_targetmap, v2))) {
				v2 = v3;
			}
			BMO_slot_map_elem_insert(weld_op, slot_targetmap, v, v2);
		}

		BMO_op_finish(bm, &find_op);
	}
	else {
		/* transform and tag verts */
		BM_ITER_MESH (v, &iter, bm, BM_VERTS_OF_MESH) {
			if (!BM_elem_flag_test(v, BM_ELEM_TAG)) {
				mul_m4_v3(mat, v->co);
				BM_elem_flag_enable(v, BM_ELEM_TAG);
			}
		}
	}
}

static void assign_mat(BMesh *bm, int id_mat)
{
	BMFace *f;
	BMIter iter;

	BM_ITER_MESH (f, &iter, bm, BM_FACES_OF_MESH) {
		if (!BM_elem_flag_test(f, BM_ELEM_TAG)) {
			f->mat_nr = id_mat;
			BM_elem_flag_enable(f, BM_ELEM_TAG);
		}
	}
}

static void merge_first_last(BMesh *bm,
                             const ArrayModifierData *amd,
                             BMOperator *dupe_first,
                             BMOperator *dupe_last,
                             BMOperator *weld_op)
{
	BMOperator find_op;
	BMOIter oiter;
	BMVert *v, *v2;
	BMOpSlot *slot_targetmap;

	BMO_op_initf(bm, &find_op, (BMO_FLAG_DEFAULTS & ~BMO_FLAG_RESPECT_HIDE),
	             "find_doubles verts=%s dist=%f keep_verts=%s",
	             dupe_first, "geom", amd->merge_dist,
	             dupe_first, "geom");

	/* append the last dupe's geom to the findop input verts */
	BMO_slot_buffer_append(&find_op,  slots_in,  "verts",
	                       dupe_last, slots_out, "geom.out");

	BMO_op_exec(bm, &find_op);

	/* add new merge targets to weld operator */
	slot_targetmap = BMO_slot_get(weld_op->slots_in, "targetmap");
	BMO_ITER (v, &oiter, find_op.slots_out, "targetmap.out", 0) {
		v2 = BMO_iter_map_value_p(&oiter);
		BMO_slot_map_elem_insert(weld_op, slot_targetmap, v, v2);
	}

	BMO_op_finish(bm, &find_op);
}

static DerivedMesh *arrayModifier_doArray(ArrayModifierData *amd,
                                          Scene *scene, Object *ob, DerivedMesh *dm,
                                          int UNUSED(initFlags))
{
	DerivedMesh *result;
	BMesh *bm = DM_to_bmesh(dm, false);
	BMOperator first_dupe_op, dupe_op, old_dupe_op, weld_op;
	BMVert **first_geom = NULL;
	int i, j;
	int index_len = -1;  /* initialize to an invalid value */
	int start, start_mc;
	float alpha = 0, d_alp = 0, circle;
	float f_o;
	/* offset matrix */
	float offset[4][4], rot[4][4];
	float final_offset[4][4];
	float mid_offset[4][4], half_offset[4][4];
	float tmp_mat[4][4], prec_mid[4][4];
	float length = amd->length;
	int count = amd->count, maxVerts;
	int *indexMap = NULL;
	DerivedMesh *start_cap = NULL, *mid_cap = NULL, *end_cap = NULL;
	MVert *src_mvert;
	BMOpSlot *slot_targetmap = NULL;  /* for weldop */
	Nurb *nu = NULL;

	/* need to avoid infinite recursion here */
	if (amd->start_cap && amd->start_cap != ob && amd->start_cap->type == OB_MESH)
		start_cap = mesh_get_derived_final(scene, amd->start_cap, CD_MASK_MESH);
	if (amd->mid_cap && amd->mid_cap != ob && amd->mid_cap->type == OB_MESH)
		mid_cap = mesh_get_derived_final(scene, amd->mid_cap, CD_MASK_MESH);
	if (amd->end_cap && amd->end_cap != ob && amd->end_cap->type == OB_MESH)
		end_cap = mesh_get_derived_final(scene, amd->end_cap, CD_MASK_MESH);

	if (amd->count_mc < 1)
		amd->count_mc = 1;
	
	unit_m4(offset);

	src_mvert = dm->getVertArray(dm);
	maxVerts = dm->getNumVerts(dm);

	if (amd->offset_type & MOD_ARR_OFF_CONST)
		add_v3_v3v3(offset[3], offset[3], amd->offset);

	if (amd->offset_type & MOD_ARR_OFF_RELATIVE) {
		for (j = 0; j < 3; j++)
			offset[3][j] += amd->scale[j] * BKE_vertarray_size(src_mvert, maxVerts, j);
	}

	if ((amd->offset_type & MOD_ARR_OFF_OBJ) && (amd->offset_ob)) {
		float obinv[4][4];
		float result_mat[4][4];

		if (ob)
			invert_m4_m4(obinv, ob->obmat);
		else
			unit_m4(obinv);

		mul_serie_m4(result_mat, offset,
		             obinv, amd->offset_ob->obmat,
		             NULL, NULL, NULL, NULL, NULL);
		copy_m4_m4(offset, result_mat);
	}

	if ((amd->type & MOD_ARR_MOD_CURVE) && (amd->curve_ob))
		length = BKE_length_fitcurve(amd, scene);

	/* calculate the maximum number of copies which will fit within the
	 * prescribed length */
	if (amd->fit_type & MOD_ARR_FITLENGTH) { // || (amd->type & MOD_ARR_MOD_CURVE)) {
		if (amd->type & MOD_ARR_MOD_NRM) {
			count = BKE_length_to_count(length, offset[3]);
			amd->count = count;
		}
		else {
			count = BKE_length_to_count(amd->length, offset[3]);
			amd->count = count;
		}
	}
	else {
		if ((amd->type & MOD_ARR_MOD_CURVE) && (amd->curve_ob))
			amd->length = length;
		else {
			amd->length = BKE_count_to_length(count, offset[3]);
			length = amd->length;
		}
	}

	if (amd->fit_type & MOD_ARR_FITBETWEEN) {
		count = count + 2;

		if (amd->type & MOD_ARR_MOD_NRM) {
			if ((amd->offset_type & MOD_ARR_OFF_OBJ) && (amd->offset_ob)) {
				//float dist = sqrt(dot_v3v3(amd->offset_ob->obmat[3], amd->offset_ob->obmat[3]));
				offset[3][0] = amd->offset_ob->obmat[3][0] / (count - 1);
				offset[3][1] = amd->offset_ob->obmat[3][1] / (count - 1);
				offset[3][2] = amd->offset_ob->obmat[3][2] / (count - 1);
			}
			else {
				offset[3][0] = offset[3][0] / (count - 1);
				offset[3][1] = offset[3][1] / (count - 1);
				offset[3][2] = offset[3][2] / (count - 1);
			}
		} 
		else { /* amd->type & MOD_ARR_MOD_CURVE*/
			Curve *cu = (Curve *)amd->curve_ob->data;
			nu = (Nurb *)cu->nurb.first;

			if (nu->bezt)
				copy_v3_v3(offset[3], nu->bezt[nu->pntsu - 1].vec[1]);
			else
				copy_v3_v3(offset[3], nu->bp[nu->pntsu * nu->pntsv - 1].vec);

			offset[3][0] = offset[3][0] / (count - 1);
			offset[3][1] = offset[3][1] / (count - 1);
			offset[3][2] = offset[3][2] / (count - 1);
		}
	}

	if (count < 1)
		count = 1;

	if ((amd->dist_mc & MOD_ARR_DIST_CURVE) && (amd->curve_cap)) {
		int dec = 0;

		Curve *cu = (Curve *)amd->curve_cap->data;
		nu = (Nurb *)cu->nurb.first;

		unit_m4(mid_offset);

		if ((amd->outer_cp & MOD_ARR_CP_FIRST) && amd->start_cap)
			dec = dec + 1;

		if ((amd->outer_cp & MOD_ARR_CP_LAST) && amd->end_cap)
			dec = dec + 1;

		if (nu->bezt) {
			if ((amd->outer_cp & MOD_ARR_CP_FIRST) && amd->start_cap)
				copy_v3_v3(mid_offset[3], nu->bezt[1].vec[1]);
			else
				copy_v3_v3(mid_offset[3], nu->bezt[0].vec[1]);

			if (amd->count_mc > (nu->pntsu - dec))
				amd->count_mc = nu->pntsu - dec;
		}
		else {
			if ((amd->outer_cp & MOD_ARR_CP_FIRST) && amd->start_cap)
				copy_v3_v3(mid_offset[3], nu->bp[1].vec);
			else
				copy_v3_v3(mid_offset[3], nu->bp[0].vec);

			if (amd->count_mc > (nu->pntsu * nu->pntsv - dec))
				amd->count_mc = nu->pntsu * nu->pntsv - dec;
		}
	}
	else 
		if (amd->count_mc >= count)
			amd->count_mc = count - 1;

	if ((amd->mode & MOD_ARR_MOD_ADV) || (amd->mode & MOD_ARR_MOD_ADV_MAT)) {
		start = 0;
		start_mc = 0;

		if (!amd->Mem_Ob)
			amd->Mem_Ob = MEM_callocN(sizeof(*amd->Mem_Ob) * count, "Mem_Ob");
		else {
			int dim = 0;

			dim = MEM_allocN_len(amd->Mem_Ob) / sizeof(*amd->Mem_Ob);
			if (dim != count) {
				amd->Mem_Ob = MEM_reallocN(amd->Mem_Ob, sizeof(*amd->Mem_Ob) * count);
				if (dim < count)
					start = dim;
				else if (dim > count)
					start = amd->count;
			}
		}

		if (amd->mid_cap) {
			if (!amd->Mem_Mat_Ob.mid_cap)
				amd->Mem_Mat_Ob.mid_cap = MEM_callocN(sizeof(int) * (amd->count_mc), "Mem_Mat_Ob");
			else {
				int dim = 0;

				dim = MEM_allocN_len(amd->Mem_Mat_Ob.mid_cap) / sizeof(*amd->Mem_Mat_Ob.mid_cap);
				if (dim != amd->count_mc) {
					amd->Mem_Mat_Ob.mid_cap = MEM_reallocN(amd->Mem_Mat_Ob.mid_cap, sizeof(*amd->Mem_Mat_Ob.mid_cap) * amd->count_mc);
					if (dim < amd->count_mc)
						start_mc = dim;
					else if (dim > amd->count_mc)
						start_mc = amd->count_mc;
				}
			}
			/* Inizializzare i nuovi cloni creati */
			if (amd->Mem_Mat_Ob.mid_cap) {
				if ((start_mc != 0) && (start_mc != amd->count_mc))
					BKE_init_mat_object_cap(start_mc, amd->count_mc, amd->Mem_Mat_Ob.mid_cap);
			}
		}

		if ((start != 0) && (start != count))
			BKE_init_offset(start, count, amd);
		BKE_create_offset(count, ob->totcol, amd, ob);
	}

	f_o = count-1;

	if ((amd->mode & MOD_ARR_MOD_ADV_CLONE) && (amd->rays > 1)) {
		alpha = (float)6.2831 / amd->rays;
		f_o = ceil((float)(count - 1) / amd->rays);
	}

	/* calculate the offset matrix of the final copy (for merging) */
	unit_m4(final_offset);

	for (j = 0; j < f_o; j++) {
		float tmp_mat[4][4];
		mul_m4_m4m4(tmp_mat, offset, final_offset);
		copy_m4_m4(final_offset, tmp_mat);
	}

	if ((amd->mid_cap) && (amd->mode & MOD_ARR_MOD_ADV_MID)) {
		if (amd->dist_mc & MOD_ARR_DIST_SEQ) {
			copy_m4_m4(mid_offset, offset);
			mid_offset[3][0] = mid_offset[3][0] / 2;
			mid_offset[3][1] = mid_offset[3][1] / 2;
			mid_offset[3][2] = mid_offset[3][2] / 2;

			if (amd->mode & MOD_ARR_MOD_ADV) {
				if (amd->Mem_Ob[0].transform) {
					float app[4][4];

					unit_m4(app);
					loc_eul_size_to_mat4(app, amd->Mem_Ob[0].loc, amd->Mem_Ob[0].rot, amd->Mem_Ob[0].scale);
					copy_m4_m4(prec_mid, mid_offset);
					copy_m4_m4(tmp_mat, mid_offset);
					mul_m4_m4m4(mid_offset, tmp_mat, app);
				}
			}
		}
		else if (amd->dist_mc & MOD_ARR_DIST_HALF) {
			copy_m4_m4(half_offset, final_offset);
			half_offset[3][0] = half_offset[3][0] / (amd->count_mc + 1);
			half_offset[3][1] = half_offset[3][1] / (amd->count_mc + 1);
			half_offset[3][2] = half_offset[3][2] / (amd->count_mc + 1);
			copy_m4_m4(mid_offset, half_offset);
		}
	}
	
	copy_m4_m4(amd->delta, offset);

	/* BMESH_TODO: bumping up the stack level avoids computing the normals
	 * after every top-level operator execution (and this modifier has the
	 * potential to execute a *lot* of top-level BMOps. There should be a
	 * cleaner way to do this. One possibility: a "mirror" BMOp would
	 * certainly help by compressing it all into one top-level BMOp that
	 * executes a lot of second-level BMOps. */
	BM_mesh_elem_toolflags_ensure(bm);
	BMO_push(bm, NULL);
	bmesh_edit_begin(bm, 0);

	if (amd->flags & MOD_ARR_MERGE) {
		BMO_op_init(bm, &weld_op, (BMO_FLAG_DEFAULTS & ~BMO_FLAG_RESPECT_HIDE),
		            "weld_verts");

		slot_targetmap = BMO_slot_get(weld_op.slots_in, "targetmap");
	}

	BMO_op_initf(bm, &dupe_op, (BMO_FLAG_DEFAULTS & ~BMO_FLAG_RESPECT_HIDE),
	             "duplicate geom=%avef");
	first_dupe_op = dupe_op;

	if ((amd->mode & MOD_ARR_MOD_ADV_CLONE) && (amd->rays > 1)) {
		d_alp = 0;
		circle = 0;

		unit_m4(rot);
		if (amd->rays_dir == MOD_ARR_RAYS_X)
			rotate_m4(rot, 'X', alpha);
		else if (amd->rays_dir == MOD_ARR_RAYS_Y)
			rotate_m4(rot, 'Y', alpha);
		else
			rotate_m4(rot, 'Z', alpha);
	}

	for (j = 0; j < count - 1; j++) {
		BMVert *v, *v2, *v3;
		BMOpSlot *geom_slot;
		BMOpSlot *geom_out_slot;
		BMOIter oiter;
		float loc[4][4];
		float ray[4][4];
		float orig[4][4];
		
		if ((amd->mode & MOD_ARR_MOD_ADV) && (amd->Mem_Ob[j].transform)) {
			if (j > 0) {
				copy_m4_m4(orig, offset);
				BKE_reset_offset(orig, j, amd->Mem_Ob[j - 1].loc, amd->Mem_Ob[j - 1].rot, amd->Mem_Ob[j - 1].scale);
			}
			
			if (amd->flag_offset & MOD_ARR_LOCAL) {
				copy_m4_m4(loc, offset);
				loc[3][0] *= j + 1;
				loc[3][1] *= j + 1;
				loc[3][2] *= j + 1;
			}
		}

		if (j != 0) {
			BMO_op_initf(bm, &dupe_op,
			             (BMO_FLAG_DEFAULTS & ~BMO_FLAG_RESPECT_HIDE),
			             "duplicate geom=%S", &old_dupe_op, "geom.out");
		}
		BMO_op_exec(bm, &dupe_op);

		geom_slot   = BMO_slot_get(dupe_op.slots_in,  "geom");
		geom_out_slot = BMO_slot_get(dupe_op.slots_out, "geom.out");

		if ((amd->flags & MOD_ARR_MERGEFINAL) && (j == 0)) {
			int first_geom_bytes = sizeof(BMVert *) * geom_slot->len;
				
			/* make a copy of the initial geometry ordering so the
			 * last duplicate can be merged into it */
			first_geom = MEM_mallocN(first_geom_bytes, "first_geom");
			memcpy(first_geom, geom_slot->data.buf, first_geom_bytes);
		}

		if ((amd->mode & MOD_ARR_MOD_ADV_CLONE) && (amd->rays > 1))
			copy_m4_m4(ray, offset);

		/* apply transformation matrix */
		BMO_ITER (v, &oiter, dupe_op.slots_out, "geom.out", BM_VERT) {

			if ((amd->mode & MOD_ARR_MOD_ADV_CLONE) && (amd->rays > 1)) {
				if (d_alp == 0) {
					if (circle > 0)
						mul_m4_v3(rot, v->co);

					mul_m4_v3(ray, v->co);
				}
				else
					mul_m4_v3(rot, v->co);
			}
			else
				mul_m4_v3(offset, v->co);

			if (amd->mode & MOD_ARR_MOD_ADV) {	
				if (amd->Mem_Ob[j].transform) {
					float app[4][4];

					unit_m4(app);
					loc_eul_size_to_mat4(app, amd->Mem_Ob[j].loc, amd->Mem_Ob[j].rot, amd->Mem_Ob[j].scale);
					if (amd->flag_offset & MOD_ARR_LOCAL) {
						/* recalculates the offset of the clone */
						loc[3][0] *= -1;
						loc[3][1] *= -1;
						loc[3][2] *= -1;
						mul_m4_v3(loc, v->co);
						/* calculates the new coordinates of the clone */
						mul_m4_v3(app, v->co);
						loc[3][0] *= -1;
						loc[3][1] *= -1;
						loc[3][2] *= -1;
						mul_m4_v3(loc, v->co);
					}
					else {
						if (j > 0)
							mul_m4_v3(orig, v->co);
						//print_m4("app", app);
						if (j == 0)
							mul_m4_v3(app, v->co);
					}
				}
			}
				
			if ((amd->type & MOD_ARR_MOD_CURVE) && (amd->curve_ob)) {
				if (amd->dist_cu & MOD_ARR_DIST_EVENLY) {

				}
				else { /* MOD_ARR_DIST_SEGMENT */
					
				}
			}
		}

		if (amd->mode & MOD_ARR_MOD_ADV_MAT) {
			if ((amd->rand_mat & MOD_ARR_MAT) && (amd->mat_ob & MOD_ARR_AR_MAT_RND)) {
				BMFace *f;
				BMO_ITER (f, &oiter, dupe_op.slots_out, "geom.out", BM_FACE) {
					f->mat_nr = amd->Mem_Ob[j].id_mat;
				}
			}
		}

		if ((amd->mode & MOD_ARR_MOD_ADV_CLONE) && (amd->rays > 1)) {
			d_alp = d_alp + alpha;
			if (d_alp > 6.2831) {
				d_alp = 0;
				circle++;
				ray[3][0] *= circle;
				ray[3][1] *= circle;
				ray[3][2] *= circle;
			}
		}

		if (amd->flags & MOD_ARR_MERGE) {
			/*calculate merge mapping*/
			if (j == 0) {
				indexMap = find_doubles_index_map(bm, &dupe_op,
				                                  amd, &index_len);
			}

			#define _E(s, i) ((BMVert **)(s)->data.buf)[i]

			/* ensure this is set */
			BLI_assert(index_len != -1);

			for (i = 0; i < index_len; i++) {
				if (!indexMap[i]) continue;

				/* merge v (from 'geom.out') into v2 (from old 'geom') */
				v = _E(geom_out_slot, i - geom_slot->len);
				v2 = _E(geom_slot, indexMap[i] - 1);

				/* check in case the target vertex (v2) is already marked
				 * for merging */
				while ((v3 = BMO_slot_map_elem_get(slot_targetmap, v2))) {
					v2 = v3;
				}

				BMO_slot_map_elem_insert(&weld_op, slot_targetmap, v, v2);
			}

			#undef _E
		}

		/* already copied earlier, but after executation more slot
		 * memory may be allocated */
		if (j == 0)
			first_dupe_op = dupe_op;
		
		if (j >= 2)
			BMO_op_finish(bm, &old_dupe_op);
		old_dupe_op = dupe_op;
	}

	if ((amd->flags & MOD_ARR_MERGE) &&
	    (amd->flags & MOD_ARR_MERGEFINAL) &&
	    (count > 1))
	{
		/* Merge first and last copies. Note that we can't use the
		 * indexMap for this because (unless the array is forming a
		 * loop) the offset between first and last is different from
		 * dupe X to dupe X+1. */

		merge_first_last(bm, amd, &first_dupe_op, &dupe_op, &weld_op);
	}

	/* start capping */
	if (start_cap || mid_cap || end_cap) {
		BM_mesh_elem_hflag_enable_all(bm, BM_VERT, BM_ELEM_TAG, FALSE);

		if (start_cap) {
			float startoffset[4][4];

			invert_m4_m4(startoffset, offset);
			if (amd->dist_mc & MOD_ARR_DIST_CURVE) {
				if (amd->outer_cp & MOD_ARR_CP_FIRST && amd->start_cap && amd->curve_cap) {
					if (nu->bezt)
						copy_v3_v3(startoffset[3], nu->bezt[0].vec[1]);
					else
						copy_v3_v3(startoffset[3], nu->bp[0].vec);
				}
			}
			bm_merge_dm_transform(bm, start_cap, startoffset, amd, 
			                      &first_dupe_op, first_dupe_op.slots_in, "geom", &weld_op);

			if (amd->mode & MOD_ARR_MOD_ADV_MAT)
				if ((amd->rand_mat & MOD_ARR_MAT) && (amd->mat_ob & MOD_ARR_SC_MAT_RND))
					assign_mat(bm, amd->Mem_Mat_Ob.start_cap);
		}

		if ((mid_cap) && ((amd->dist_mc & MOD_ARR_DIST_SEQ) || (amd->dist_mc & MOD_ARR_DIST_HALF) ||
			((amd->dist_mc & MOD_ARR_DIST_CURVE) && (amd->curve_cap))))
		{
			int count_midcap, inc, j;
			/* d_alp_md = 0;
			 * if (amd->rays > 1)
			 *	  alpha_md = (float)6.2831 / amd->rays; */
			if (amd->mode & MOD_ARR_MOD_ADV_MID)
				count_midcap = amd->count_mc;
			else {
				count_midcap = 1;
				copy_m4_m4(mid_offset, offset);
				mid_offset[3][0] = mid_offset[3][0] / 2;
				mid_offset[3][1] = mid_offset[3][1] / 2;
				mid_offset[3][2] = mid_offset[3][2] / 2;
			}

			for (j = 0; j < count_midcap; j++) {
				bm_merge_dm_transform(bm, mid_cap, mid_offset, amd,
					                  &first_dupe_op, first_dupe_op.slots_in, "geom", &weld_op);

				if (amd->mode & MOD_ARR_MOD_ADV_MAT)
					if ((amd->rand_mat & MOD_ARR_MAT) && (amd->mat_ob & MOD_ARR_MC_MAT_RND))
						assign_mat(bm, amd->Mem_Mat_Ob.mid_cap[j]);

				if (amd->dist_mc & MOD_ARR_DIST_SEQ) {
					/********/
					/*if (amd->rays > 1) {
						float rot[4][4];
						unit_m4(rot);
						if (amd->rays_dir == MOD_ARR_RAYS_X)
							rotate_m4(rot,'X',d_alp_md);
						else if (amd->rays_dir == MOD_ARR_RAYS_Y)
							rotate_m4(rot,'Y',d_alp_md);
						else
							rotate_m4(rot,'Z',d_alp_md);
						if (d_alp_md == 0){
							mult_m4_m4m4(mid_offset, tmat, offset);

							copy_m4_m4(tmat, mid_offset);
							mult_m4_m4m4(mid_offset, tmat, rot);
						}
						else{
							mult_m4_m4m4(mid_offset, tmat, offset);

							copy_m4_m4(tmat, mid_offset);
							mult_m4_m4m4(mid_offset, tmat, rot);
						}
					}
					else
						mult_m4_m4m4(mid_offset, tmat, offset);
					*/
					/********/
					if (amd->mode & MOD_ARR_MOD_ADV) {
						if (amd->Mem_Ob[j + 1].transform)
							copy_m4_m4(mid_offset, prec_mid);
					}
					mul_m4_m4m4(tmp_mat, offset, mid_offset);
					copy_m4_m4(mid_offset, tmp_mat);
					if (amd->mode & MOD_ARR_MOD_ADV) {
						if (amd->Mem_Ob[j + 1].transform) {
							float app[4][4];
							
							unit_m4(app);
							loc_eul_size_to_mat4(app, amd->Mem_Ob[j + 1].loc, amd->Mem_Ob[j + 1].rot, 
												 amd->Mem_Ob[j + 1].scale);
							copy_m4_m4(prec_mid, mid_offset);
							copy_m4_m4(tmp_mat, mid_offset);
							mul_m4_m4m4(mid_offset, tmp_mat, app);
						}
					}
				}
				else if (amd->dist_mc & MOD_ARR_DIST_HALF) {
					mul_m4_m4m4(tmp_mat, half_offset, mid_offset);
					copy_m4_m4(mid_offset, tmp_mat);
				}
				else if (amd->dist_mc & MOD_ARR_DIST_CURVE) {
					if (amd->outer_cp & MOD_ARR_CP_FIRST && amd->start_cap)
						inc = 1;
					else 
						inc = 0;
					if ((j + 1) < amd->count_mc){
						if (nu->bezt)
							copy_v3_v3(mid_offset[3], nu->bezt[j + 1 + inc].vec[1]);
						else
							copy_v3_v3(mid_offset[3], nu->bp[j + 1 + inc].vec);
					}
				}
			}
		}

		if (end_cap) {
			float endoffset[4][4];

			mul_m4_m4m4(endoffset, offset, final_offset);
			if (amd->dist_mc & MOD_ARR_DIST_CURVE) {
				if (amd->outer_cp & MOD_ARR_CP_LAST && amd->end_cap && amd->curve_cap) {
					if (nu->bezt)
						copy_v3_v3(endoffset[3], nu->bezt[nu->pntsu - 1].vec[1]);
					else
						copy_v3_v3(endoffset[3], nu->bp[nu->pntsu * nu->pntsv - 1].vec);
				}
			}
			bm_merge_dm_transform(bm, end_cap, endoffset, amd,
			                      &dupe_op, (count == 1) ? dupe_op.slots_in : dupe_op.slots_out,
			                      (count == 1) ? "geom" : "geom.out", &weld_op);
			if (amd->mode & MOD_ARR_MOD_ADV_MAT)
				if ((amd->rand_mat & MOD_ARR_MAT) && (amd->mat_ob & MOD_ARR_EC_MAT_RND))
					assign_mat(bm, amd->Mem_Mat_Ob.end_cap);
		}
	}
	/* done capping */

	/* free remaining dupe operators */
	BMO_op_finish(bm, &first_dupe_op);
	if (count > 2)
		BMO_op_finish(bm, &dupe_op);

	/* run merge operator */
	if (amd->flags & MOD_ARR_MERGE) {
		BMO_op_exec(bm, &weld_op);
		BMO_op_finish(bm, &weld_op);
	}

	/* Bump the stack level back down to match the adjustment up above */
	BMO_pop(bm);

	result = CDDM_from_bmesh(bm, FALSE);

	if ((dm->dirty & DM_DIRTY_NORMALS) ||
	    ((amd->offset_type & MOD_ARR_OFF_OBJ) && (amd->offset_ob)))
	{
		/* Update normals in case offset object has rotation. */
		result->dirty |= DM_DIRTY_NORMALS;
	}

	BM_mesh_free(bm);

	if (indexMap)
		MEM_freeN(indexMap);
	if (first_geom)
		MEM_freeN(first_geom);

	return result;
}

static DerivedMesh *applyModifier(ModifierData *md, Object *ob,
                                  DerivedMesh *dm,
                                  ModifierApplyFlag UNUSED(flag))
{
	DerivedMesh *result;
	ArrayModifierData *amd = (ArrayModifierData *) md;

	result = arrayModifier_doArray(amd, md->scene, ob, dm, 0);

	if(amd->arr_group!=NULL) {
		if (amd->mode & MOD_ARR_MOD_ADV_CLONE) {
			ob->transflag = OB_DUPLIARRAY;
			ob->dup_group = amd->arr_group;
		}
		else {
			ob->transflag = 0;
			ob->dup_group = NULL;
		}
	}
	else {
		ob->transflag = 0;
		ob->dup_group = NULL;
	}
	/* } */
	return result;
}


static void freeData(ModifierData *md)
{
	ArrayModifierData *amd = (ArrayModifierData*) md;
	
	if (amd) {
		if (amd->Mem_Ob)
			MEM_freeN(amd->Mem_Ob);

		if (amd->Mem_Mat_Ob.mid_cap)
			MEM_freeN(amd->Mem_Mat_Ob.mid_cap);
	}
}

ModifierTypeInfo modifierType_Array = {
	/* name */              "Array",
	/* structName */        "ArrayModifierData",
	/* structSize */        sizeof(ArrayModifierData),
	/* type */              eModifierTypeType_Constructive,
	/* flags */             eModifierTypeFlag_AcceptsMesh |
	                        eModifierTypeFlag_SupportsMapping |
	                        eModifierTypeFlag_SupportsEditmode |
	                        eModifierTypeFlag_EnableInEditmode |
	                        eModifierTypeFlag_AcceptsCVs,

	/* copyData */          copyData,
	/* deformVerts */       NULL,
	/* deformMatrices */    NULL,
	/* deformVertsEM */     NULL,
	/* deformMatricesEM */  NULL,
	/* applyModifier */     applyModifier,
	/* applyModifierEM */   NULL,
	/* initData */          initData,
	/* requiredDataMask */  NULL,
	/* freeData */          freeData,
	/* isDisabled */        NULL,
	/* updateDepgraph */    updateDepgraph,
	/* dependsOnTime */     NULL,
	/* dependsOnNormals */	NULL,
	/* foreachObjectLink */ foreachObjectLink,
	/* foreachIDLink */     NULL,
	/* foreachTexLink */    NULL,
};
