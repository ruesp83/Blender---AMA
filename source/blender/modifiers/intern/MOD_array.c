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

#include <time.h>
#include <math.h>
#include "MEM_guardedalloc.h"

#include "BLI_math.h"
#include "BLI_utildefines.h"
#include "BLI_ghash.h"
#include "BLI_edgehash.h"

#include "DNA_curve_types.h"
#include "DNA_group_types.h"
#include "DNA_meshdata_types.h"
#include "DNA_object_types.h"

#include "BKE_cdderivedmesh.h"
#include "BKE_displist.h"
#include "BKE_mesh.h"
#include "BKE_modifier.h"
#include "BKE_object.h"
#include "BKE_anim.h"

//#include "BLI_rand.h"

#include "depsgraph_private.h"
#include "BKE_ama.h"
#include "MOD_util.h"

static void initData(ModifierData *md)
{
	ArrayModifierData *amd = (ArrayModifierData*) md;

	/* default to 2 duplicates distributed along the x-axis by an
	offset of 1 object-width
	*/
	amd->start_cap = amd->mid_cap = amd->end_cap = amd->curve_cap = amd->curve_ob = amd->offset_ob = NULL;
	amd->count = 2;
	amd->offset[0] = amd->offset[1] = amd->offset[2] = 0;
	amd->scale[0] = 1;
	amd->scale[1] = amd->scale[2] = 0;
	amd->length = 0;
	amd->merge_dist = 0.01;
	amd->fit_type = MOD_ARR_FIXEDCOUNT;
	amd->offset_type = MOD_ARR_OFF_RELATIVE;
	amd->flags = 0;

	amd->type = MOD_ARR_MOD_NRM;
	amd->mode = !MOD_ARR_MOD_ADV;
	amd->loc_offset[0] = amd->loc_offset[1] = amd->loc_offset[2] = 0;
	amd->rot_offset[0] = amd->rot_offset[1] = amd->rot_offset[2] = 0;
	amd->scale_offset[0] = amd->scale_offset[1] = amd->scale_offset[2] = 1;
	amd->sign = MOD_ARR_SIGN_P;
	amd->lock = !MOD_ARR_LOCK;
	amd->flag_offset = MOD_ARR_PROP;
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
}

static void copyData(ModifierData *md, ModifierData *target)
{
	ArrayModifierData *amd = (ArrayModifierData*) md;
	ArrayModifierData *tamd = (ArrayModifierData*) target;

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

	tamd->type = amd->type;
	tamd->mode = amd->mode;
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
	ArrayModifierData *amd = (ArrayModifierData*) md;

	walk(userData, ob, &amd->start_cap);
	walk(userData, ob, &amd->mid_cap);
	walk(userData, ob, &amd->end_cap);
	walk(userData, ob, &amd->curve_ob);
	walk(userData, ob, &amd->curve_cap);
	walk(userData, ob, &amd->offset_ob);
}

static void updateDepgraph(ModifierData *md, DagForest *forest,
	struct Scene *UNUSED(scene), Object *UNUSED(ob), DagNode *obNode)
{
	ArrayModifierData *amd = (ArrayModifierData*) md;

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
	if (amd->curve_ob) {
		DagNode *curNode = dag_get_node(forest, amd->curve_ob);

		dag_add_relation(forest, curNode, obNode,
		                 DAG_RL_DATA_DATA | DAG_RL_OB_DATA, "Array Modifier");
	}
	if (amd->curve_cap) {
		DagNode *curNode = dag_get_node(forest, amd->curve_cap);

		dag_add_relation(forest, curNode, obNode,
		                 DAG_RL_DATA_DATA | DAG_RL_OB_DATA, "Array Modifier");
	}
	if (amd->offset_ob) {
		DagNode *curNode = dag_get_node(forest, amd->offset_ob);

		dag_add_relation(forest, curNode, obNode,
		                 DAG_RL_DATA_DATA | DAG_RL_OB_DATA, "Array Modifier");
	}
}

static DerivedMesh *arrayModifier_doArray(ArrayModifierData *amd,
					  struct Scene *scene, Object *ob, DerivedMesh *dm,
	   int initFlags)
{
	/* offset matrix */
	float offset[4][4];
	float final_offset[4][4];
	float mid_offset[4][4], half_offset[4][4];
	float tmp_mat[4][4], prec_mid[4][4];
	float length = amd->length;
	float alpha = 0, d_alp = 0, circle;
	float f_o;
	int i, j, flag;
	int start, start_mc;
	int count = amd->count;
	int numVerts, numEdges, numFaces;
	int maxVerts, maxEdges, maxFaces;
	int finalVerts, finalEdges, finalFaces;
	DerivedMesh *result, *start_cap = NULL, *mid_cap = NULL, *end_cap = NULL;
	MVert *mvert, *src_mvert;
	MEdge *medge;
	MFace *mface;
	Nurb *nu = NULL;

	IndexMapEntry *indexMap;

	EdgeHash *edges;
	
	/* need to avoid infinite recursion here */
	if(amd->start_cap && amd->start_cap != ob)
		start_cap = amd->start_cap->derivedFinal;
	if(amd->mid_cap && amd->mid_cap != ob)
		mid_cap = amd->mid_cap->derivedFinal;
	if(amd->end_cap && amd->end_cap != ob)
		end_cap = amd->end_cap->derivedFinal;
	
	if (amd->count_mc < 1)
		amd->count_mc = 1;

	unit_m4(offset);

	indexMap = MEM_callocN(sizeof(*indexMap) * dm->getNumVerts(dm),
				   "indexmap");

	src_mvert = dm->getVertArray(dm);

	maxVerts = dm->getNumVerts(dm);

	if(amd->offset_type & MOD_ARR_OFF_CONST)
		add_v3_v3(offset[3], amd->offset);
	if(amd->offset_type & MOD_ARR_OFF_RELATIVE) {
		for(j = 0; j < 3; j++)
			offset[3][j] += amd->scale[j] * vertarray_size(src_mvert,
					maxVerts, j);
	}

	if((amd->offset_type & MOD_ARR_OFF_OBJ) && (amd->offset_ob)) {
		float obinv[4][4];
		float result_mat[4][4];

		if(ob)
			invert_m4_m4(obinv, ob->obmat);
		else
			unit_m4(obinv);

		mul_serie_m4(result_mat, offset,
		             obinv, amd->offset_ob->obmat,
		             NULL, NULL, NULL, NULL, NULL);
		copy_m4_m4(offset, result_mat);
	}
	
	if ((amd->type & MOD_ARR_MOD_CURVE) && (amd->curve_ob))
		length = length_fitcurve(amd, scene);

	/* calculate the maximum number of copies which will fit within the
	prescribed length */
	if (amd->fit_type & MOD_ARR_FITLENGTH) { // || (amd->type & MOD_ARR_MOD_CURVE)) {
		if (amd->type & MOD_ARR_MOD_NRM) {
			count = length_to_count(length, offset[3]);
			amd->count = count;
		}
		else {
			count = length_to_count(amd->length, offset[3]);
			amd->count = count;
		}
	}
	else {
		if ((amd->type & MOD_ARR_MOD_CURVE) && (amd->curve_ob)) {
			amd->length = length;
		}
		else {
			amd->length = count_to_length(count, offset[3]);
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
			offset[3][0] = offset[3][0] * (count - 1);
			offset[3][1] = offset[3][1] * (count - 1);
			offset[3][2] = offset[3][2] * (count - 1);
		}
	}

	if(count < 1)
		count = 1;
	
	/**/
	if ((amd->dist_mc & MOD_ARR_DIST_CURVE) && (amd->curve_cap)){
		int dec = 0;
		Curve *cu = amd->curve_cap->data;
		nu = cu->nurb.first;
		unit_m4(mid_offset);
		if (amd->outer_cp & MOD_ARR_CP_FIRST && amd->start_cap)
			dec = dec + 1;
		if (amd->outer_cp & MOD_ARR_CP_LAST && amd->end_cap)
			dec = dec + 1;
		if (nu->bezt) {
			if (amd->outer_cp & MOD_ARR_CP_FIRST && amd->start_cap)
				copy_v3_v3(mid_offset[3], nu->bezt[1].vec[1]);
			else
				copy_v3_v3(mid_offset[3], nu->bezt[0].vec[1]);
			if (amd->count_mc > (nu->pntsu - dec))
				amd->count_mc = nu->pntsu - dec;
		}
		else {
			if (amd->outer_cp & MOD_ARR_CP_FIRST && amd->start_cap)
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

	/* allocate memory for count duplicates (including original) plus
		  * start, mid and end caps
	*/
	finalVerts = dm->getNumVerts(dm) * count;
	finalEdges = dm->getNumEdges(dm) * count;
	finalFaces = dm->getNumFaces(dm) * count;
	
	if(start_cap) {
		finalVerts += start_cap->getNumVerts(start_cap);
		finalEdges += start_cap->getNumEdges(start_cap);
		finalFaces += start_cap->getNumFaces(start_cap);
	}
	if(mid_cap) {
		finalVerts += mid_cap->getNumVerts(mid_cap) * amd->count_mc;
		finalEdges += mid_cap->getNumEdges(mid_cap) * amd->count_mc;
		finalFaces += mid_cap->getNumFaces(mid_cap) * amd->count_mc;
	}
	if(end_cap) {
		finalVerts += end_cap->getNumVerts(end_cap);
		finalEdges += end_cap->getNumEdges(end_cap);
		finalFaces += end_cap->getNumFaces(end_cap);
	}
	
	result = CDDM_from_template(dm, finalVerts, finalEdges, finalFaces);
	
	flag = 0;
	if ((amd->mode & MOD_ARR_MOD_ADV) || (amd->mode & MOD_ARR_MOD_ADV_MAT)){
		start = 0;
		start_mc = 0;
		if (!amd->Mem_Ob)
			amd->Mem_Ob = MEM_callocN(sizeof(*amd->Mem_Ob) * count, "Mem_Ob");
		else {
			int dim = 0;
			dim = MEM_allocN_len(amd->Mem_Ob) / sizeof(*amd->Mem_Ob);
			if (dim < count) {
				amd->Mem_Ob = MEM_reallocN(amd->Mem_Ob, sizeof(*amd->Mem_Ob) * count);
				start = dim;
			}
			else if (dim > count) {
				amd->Mem_Ob = MEM_reallocN(amd->Mem_Ob, sizeof(*amd->Mem_Ob) * count);
				start = amd->count_mc;
			}
		}

		if (!amd->Mem_Mat_Ob.mid_cap)
			amd->Mem_Mat_Ob.mid_cap = MEM_callocN(sizeof(int) * (amd->count_mc), "Mem_Mat_Ob");
		else {
			int dim = 0;
			dim = MEM_allocN_len(amd->Mem_Mat_Ob.mid_cap) / sizeof(*amd->Mem_Mat_Ob.mid_cap);
			if (dim < amd->count_mc) {
				amd->Mem_Mat_Ob.mid_cap = MEM_reallocN(amd->Mem_Mat_Ob.mid_cap, sizeof(*amd->Mem_Mat_Ob.mid_cap) * amd->count_mc);
				start_mc = dim;
			}
			else if (dim > amd->count_mc) {
				amd->Mem_Mat_Ob.mid_cap = MEM_reallocN(amd->Mem_Mat_Ob.mid_cap, sizeof(*amd->Mem_Mat_Ob.mid_cap) * amd->count_mc);
				start_mc = amd->count_mc;
			}
		}
		//Inizializzare i nuovi cloni creati
		if (!amd->lock) {
			if (amd->Mem_Mat_Ob.mid_cap) {
				if (start_mc != amd->count_mc)
					init_mat_oc(start_mc, amd->count_mc, amd->Mem_Mat_Ob.mid_cap);
			}
			if (start != count)
				init_offset(start, count, amd);
			create_offset(count, ob->totcol, amd, ob);
			amd->lock = 1;
			flag = 1;
		}
		else {
			if ((start != 0) && (start != count)) // || (amd->mode & MOD_ARR_MOD_ADV_MAT))
				init_offset(start, count, amd);
			if ((start_mc != 0) && (start_mc != amd->count_mc))
 				init_mat_oc(start_mc, amd->count_mc, amd->Mem_Mat_Ob.mid_cap);
		}
	}
	
	f_o = count-1;

	if (amd->rays>1) {
		alpha = (float)6.2831 / amd->rays;
		circle = (float)(count - 1) / amd->rays;
		f_o = ceil(circle);
	}

	/* calculate the offset matrix of the final copy (for merging) */
	unit_m4(final_offset);
	for(j=0; j < f_o; j++) {
		mult_m4_m4m4(tmp_mat, offset, final_offset);
		copy_m4_m4(final_offset, tmp_mat);
	}
	
	if (amd->dist_mc & MOD_ARR_DIST_SEQ){
		copy_m4_m4(mid_offset, offset);
		mid_offset[3][0] = mid_offset[3][0] / 2;
		mid_offset[3][1] = mid_offset[3][1] / 2;
		mid_offset[3][2] = mid_offset[3][2] / 2;
		if (amd->mode & MOD_ARR_MOD_ADV) {
			if (amd->Mem_Ob[0].transform == 1) {
				float app[4][4];
				unit_m4(app);
				loc_eul_size_to_mat4(app, amd->Mem_Ob[0].loc, amd->Mem_Ob[0].rot, amd->Mem_Ob[0].scale);
				copy_m4_m4(prec_mid, mid_offset);
				copy_m4_m4(tmp_mat, mid_offset);
				mult_m4_m4m4(mid_offset, tmp_mat, app);
			}
		}
	}
	else if (amd->dist_mc & MOD_ARR_DIST_HALF){
		copy_m4_m4(half_offset, final_offset);
		half_offset[3][0] = half_offset[3][0] / (amd->count_mc + 1);
		half_offset[3][1] = half_offset[3][1] / (amd->count_mc + 1);
		half_offset[3][2] = half_offset[3][2] / (amd->count_mc + 1);
		copy_m4_m4(mid_offset, half_offset);
	}
	
	copy_m4_m4(amd->delta, offset);
	numVerts = numEdges = numFaces = 0;
	mvert = CDDM_get_verts(result);

	for (i = 0; i < maxVerts; i++) {
		indexMap[i].merge = -1; /* default to no merge */
		indexMap[i].merge_final = 0; /* default to no merge */
	}
	
	for (i = 0; i < maxVerts; i++) {
		MVert *inMV;
		MVert *mv = &mvert[numVerts];
		MVert *mv2;
		float co[3];
		
		inMV = &src_mvert[i];

		DM_copy_vert_data(dm, result, i, numVerts, 1);
		*mv = *inMV;
		numVerts++;

		indexMap[i].new = numVerts - 1;

		copy_v3_v3(co, mv->co);
		
		/* Attempts to merge verts from one duplicate with verts from the
			  * next duplicate which are closer than amd->merge_dist.
			  * Only the first such vert pair is merged.
			  * If verts are merged in the first duplicate pair, they are merged
			  * in all pairs.
		*/
		if((count > 1) && (amd->flags & MOD_ARR_MERGE)) {
			float tmp_co[3];
			mul_v3_m4v3(tmp_co, offset, mv->co);

			for(j = 0; j < maxVerts; j++) {
				/* if vertex already merged, don't use it */
				if( indexMap[j].merge != -1 ) continue;

				inMV = &src_mvert[j];
				/* if this vert is within merge limit, merge */
				if(compare_len_v3v3(tmp_co, inMV->co, amd->merge_dist)) {
					indexMap[i].merge = j;

					/* test for merging with final copy of merge target */
					if(amd->flags & MOD_ARR_MERGEFINAL) {
						copy_v3_v3(tmp_co, inMV->co);
						inMV = &src_mvert[i];
						mul_m4_v3(final_offset, tmp_co);
						if(compare_len_v3v3(tmp_co, inMV->co, amd->merge_dist))
							indexMap[i].merge_final = 1;
					}
					break;
				}
			}
		}
		/* if no merging, generate copies of this vert */
		if(indexMap[i].merge < 0) {
			if (amd->rays>1)
				d_alp=0;

			for(j=0; j < count - 1; j++) {
				float rot[4][4];
				mv2 = &mvert[numVerts];

				DM_copy_vert_data(result, result, numVerts - 1, numVerts, 1);
				*mv2 = *mv;
				numVerts++;
				/*//Aggiungendo questa parte ed eliminando 1) e 2) 
				//si ottiene una spirale
				mul_m4_v3(offset, co);
				copy_v3_v3(mv2->co, co);*/
				if (amd->rays>1) {
					float ro[3];
					unit_m4(rot);
					if (amd->rays_dir == MOD_ARR_RAYS_X)
						rotate_m4(rot,'X',d_alp);
					else if (amd->rays_dir == MOD_ARR_RAYS_Y)
						rotate_m4(rot,'Y',d_alp);
					else
						rotate_m4(rot,'Z',d_alp);
					if (d_alp == 0){
						/*1)*/
						mul_m4_v3(offset, co);
						copy_v3_v3(mv2->co, co);
						/******/
						copy_v3_v3(ro, mv2->co);
						mul_m4_v3(rot, ro);
						copy_v3_v3(mv2->co, ro);
					}
					else {
						copy_v3_v3(ro,co);
						mul_m4_v3(rot, ro);
						copy_v3_v3(mv2->co, ro);
					}
					d_alp = d_alp + alpha;
					if (d_alp>6.2831)
						d_alp=0;
				}
				else {
					/*2)*/
					mul_m4_v3(offset, co);
					copy_v3_v3(mv2->co, co);
					/******/
				}

				if (amd->mode & MOD_ARR_MOD_ADV) {	
					if (amd->Mem_Ob[j].transform) {
						float fo[3];
						float app[4][4];
						unit_m4(app);

						if (amd->flag_offset & MOD_ARR_LOCAL) {
							float loc[4][4];

							copy_v3_v3(fo, mv->co);
							copy_v3_v3(app[3], fo);
							loc_eul_size_to_mat4(app, amd->Mem_Ob[j].loc, amd->Mem_Ob[j].rot, amd->Mem_Ob[j].scale);
							mul_m4_v3(app, fo);
							copy_m4_m4(loc, offset);
							mul_v3_fl(loc[3], j+1);
							mul_m4_v3(loc, fo);
							copy_v3_v3(mv2->co, fo);
						}
						else {
							loc_eul_size_to_mat4(app, amd->Mem_Ob[j].loc, amd->Mem_Ob[j].rot, amd->Mem_Ob[j].scale);
							copy_v3_v3(fo, mv2->co);
							mul_m4_v3(app, fo);
							copy_v3_v3(mv2->co, fo);
						}
					}
				}
				
				if ((amd->type & MOD_ARR_MOD_CURVE) && (amd->curve_ob)) {
					if (amd->dist_cu & MOD_ARR_DIST_EVENLY){
						if (i==0){
							float fo[3];
							float cent[3];
							float loc[4];
							float app[4][4];

							unit_m4(app);
							copy_v3_v3(fo, mv2->co);
							array_to_curve(scene, amd->curve_ob, ob, fo, loc, cent);
							copy_v3_v3(amd->Mem_Ob[j].cu_cent, cent);
							copy_v4_v4(amd->Mem_Ob[j].cu_loc, loc);
							copy_v3_v3(app[3], amd->Mem_Ob[j].cu_loc);
							mul_m4_v3(app, fo);
							//add_v3_v3v3(fo, cent, loc);
							copy_v3_v3(mv2->co, fo);
						}
						else {
							float fo[3];
							float app[4][4];
							unit_m4(app);
							copy_v3_v3(fo, mv2->co);
							/*print_v3("cent", amd->Mem_Ob[j].cu_cent);
							print_v4("loc", amd->Mem_Ob[j].cu_loc);*/
							//add_v3_v3v3(fo, amd->Mem_Ob[j].cu_cent, amd->Mem_Ob[j].cu_loc);
							copy_v3_v3(app[3], amd->Mem_Ob[j].cu_loc);
							mul_m4_v3(app, fo);
							print_v3("fo", fo);
							copy_v3_v3(mv2->co, fo);
						}
					}
				}
			}
		} else if(indexMap[i].merge != i && indexMap[i].merge_final) {
			/* if this vert is not merging with itself, and it is merging
				  * with the final copy of its merge target, remove the first copy
			*/
			numVerts--;
			DM_free_vert_data(result, numVerts, 1);
		}
	}

	/* make a hashtable so we can avoid duplicate edges from merging */
	edges = BLI_edgehash_new();

	maxEdges = dm->getNumEdges(dm);
	medge = CDDM_get_edges(result);
	for(i = 0; i < maxEdges; i++) {
		MEdge inMED;
		MEdge med;
		MEdge *med2;
		int vert1, vert2;

		dm->getEdge(dm, i, &inMED);

		med = inMED;
		med.v1 = indexMap[inMED.v1].new;
		med.v2 = indexMap[inMED.v2].new;

		/* if vertices are to be merged with the final copies of their
			  * merge targets, calculate that final copy
		*/
		if(indexMap[inMED.v1].merge_final) {
			med.v1 = calc_mapping(indexMap, indexMap[inMED.v1].merge,
			count - 1);
		}
		if(indexMap[inMED.v2].merge_final) {
			med.v2 = calc_mapping(indexMap, indexMap[inMED.v2].merge,
			count - 1);
		}

		if(med.v1 == med.v2) continue;

		/* XXX Unfortunately the calc_mapping returns sometimes numVerts... leads to bad crashes */
		if(med.v1 >= numVerts)
			med.v1= numVerts-1;
		if(med.v2 >= numVerts)
			med.v2= numVerts-1;

		if (initFlags) {
			med.flag |= ME_EDGEDRAW | ME_EDGERENDER;
		}

		if(!BLI_edgehash_haskey(edges, med.v1, med.v2)) {
			DM_copy_edge_data(dm, result, i, numEdges, 1);
			medge[numEdges] = med;
			numEdges++;

			BLI_edgehash_insert(edges, med.v1, med.v2, NULL);
		}

		for(j = 1; j < count; j++)
		{
			vert1 = calc_mapping(indexMap, inMED.v1, j);
			vert2 = calc_mapping(indexMap, inMED.v2, j);

			/* edge could collapse to single point after mapping */
			if(vert1 == vert2) continue;

			/* XXX Unfortunately the calc_mapping returns sometimes numVerts... leads to bad crashes */
			if(vert1 >= numVerts)
				vert1= numVerts-1;
			if(vert2 >= numVerts)
				vert2= numVerts-1;

			/* avoid duplicate edges */
			if(!BLI_edgehash_haskey(edges, vert1, vert2)) {
				med2 = &medge[numEdges];

				DM_copy_edge_data(dm, result, i, numEdges, 1);
				*med2 = med;
				numEdges++;

				med2->v1 = vert1;
				med2->v2 = vert2;

				BLI_edgehash_insert(edges, med2->v1, med2->v2, NULL);
			}
		}
	}

	maxFaces = dm->getNumFaces(dm);
	mface = CDDM_get_faces(result);
	for (i=0; i < maxFaces; i++) {
		MFace inMF;
		MFace *mf = &mface[numFaces];

		dm->getFace(dm, i, &inMF);

		DM_copy_face_data(dm, result, i, numFaces, 1);
		*mf = inMF;

		mf->v1 = indexMap[inMF.v1].new;
		mf->v2 = indexMap[inMF.v2].new;
		mf->v3 = indexMap[inMF.v3].new;
		if(inMF.v4)
			mf->v4 = indexMap[inMF.v4].new;

		/* if vertices are to be merged with the final copies of their
			  * merge targets, calculate that final copy
		*/
		if(indexMap[inMF.v1].merge_final)
			mf->v1 = calc_mapping(indexMap, indexMap[inMF.v1].merge, count-1);
		if(indexMap[inMF.v2].merge_final)
			mf->v2 = calc_mapping(indexMap, indexMap[inMF.v2].merge, count-1);
		if(indexMap[inMF.v3].merge_final)
			mf->v3 = calc_mapping(indexMap, indexMap[inMF.v3].merge, count-1);
		if(inMF.v4 && indexMap[inMF.v4].merge_final)
			mf->v4 = calc_mapping(indexMap, indexMap[inMF.v4].merge, count-1);

		if(test_index_face_maxvert(mf, &result->faceData, numFaces, inMF.v4?4:3, numVerts) < 3)
			continue;

		numFaces++;

		/* if the face has fewer than 3 vertices, don't create it */
		if(mf->v3 == 0 || (mf->v1 && (mf->v1 == mf->v3 || mf->v1 == mf->v4))) {
			numFaces--;
			DM_free_face_data(result, numFaces, 1);
		}

		for(j = 1; j < count; j++)
		{
			MFace *mf2 = &mface[numFaces];

			DM_copy_face_data(dm, result, i, numFaces, 1);
			*mf2 = *mf;

			mf2->v1 = calc_mapping(indexMap, inMF.v1, j);
			mf2->v2 = calc_mapping(indexMap, inMF.v2, j);
			mf2->v3 = calc_mapping(indexMap, inMF.v3, j);
			if (inMF.v4)
				mf2->v4 = calc_mapping(indexMap, inMF.v4, j);

			numFaces++;

			/*Rand Material*/
			if ((flag = 1) && (amd->mode & MOD_ARR_MOD_ADV_MAT) && (amd->mat_ob & MOD_ARR_AR_MAT_RND) && (ob->totcol>1))
			{
					mf2->mat_nr = amd->Mem_Ob[j-1].id_mat;
					flag = 0;
			}
			/* if the face has fewer than 3 vertices, don't create it */
			if(test_index_face_maxvert(mf2, &result->faceData, numFaces-1, inMF.v4?4:3, numVerts) < 3) {
				numFaces--;
				DM_free_face_data(result, numFaces, 1);
			}
		}
	}

	/* add start, mid and end caps */
	if(start_cap) {
		float startoffset[4][4];
		MVert *cap_mvert;
		MEdge *cap_medge;
		MFace *cap_mface;
		int *origindex;
		int *vert_map;
		int capVerts, capEdges, capFaces;
		
		capVerts = start_cap->getNumVerts(start_cap);
		capEdges = start_cap->getNumEdges(start_cap);
		capFaces = start_cap->getNumFaces(start_cap);
		cap_mvert = start_cap->getVertArray(start_cap);
		cap_medge = start_cap->getEdgeArray(start_cap);
		cap_mface = start_cap->getFaceArray(start_cap);

		invert_m4_m4(startoffset, offset);

		if (amd->dist_mc & MOD_ARR_DIST_CURVE) {
			if (amd->outer_cp & MOD_ARR_CP_FIRST && amd->start_cap && amd->curve_cap) {
				if (nu->bezt)
					copy_v3_v3(startoffset[3], nu->bezt[0].vec[1]);
				else
					copy_v3_v3(startoffset[3], nu->bp[0].vec);
			}
		}
		vert_map = MEM_callocN(sizeof(*vert_map) * capVerts,
		"arrayModifier_doArray vert_map");

		origindex = result->getVertDataArray(result, CD_ORIGINDEX);
		for(i = 0; i < capVerts; i++) {
			MVert *mv = &cap_mvert[i];
			short merged = 0;

			if(amd->flags & MOD_ARR_MERGE) {
				float tmp_co[3];
				MVert *in_mv;
				int j;

				copy_v3_v3(tmp_co, mv->co);
				mul_m4_v3(startoffset, tmp_co);

				for(j = 0; j < maxVerts; j++) {
					in_mv = &src_mvert[j];
					/* if this vert is within merge limit, merge */
					if(compare_len_v3v3(tmp_co, in_mv->co, amd->merge_dist)) {
						vert_map[i] = calc_mapping(indexMap, j, 0);
						merged = 1;
						break;
					}
				}
			}

			if(!merged) {
				DM_copy_vert_data(start_cap, result, i, numVerts, 1);
				mvert[numVerts] = *mv;
				mul_m4_v3(startoffset, mvert[numVerts].co);
				origindex[numVerts] = ORIGINDEX_NONE;

				vert_map[i] = numVerts;

				numVerts++;
			}
		}

		origindex = result->getEdgeDataArray(result, CD_ORIGINDEX);
		for(i = 0; i < capEdges; i++) {
			int v1, v2;

			v1 = vert_map[cap_medge[i].v1];
			v2 = vert_map[cap_medge[i].v2];

			if(!BLI_edgehash_haskey(edges, v1, v2)) {
				DM_copy_edge_data(start_cap, result, i, numEdges, 1);
				medge[numEdges] = cap_medge[i];
				medge[numEdges].v1 = v1;
				medge[numEdges].v2 = v2;
				origindex[numEdges] = ORIGINDEX_NONE;

				numEdges++;
			}
		}
		
		origindex = result->getFaceDataArray(result, CD_ORIGINDEX);
		for(i = 0; i < capFaces; i++) {
			DM_copy_face_data(start_cap, result, i, numFaces, 1);
			mface[numFaces] = cap_mface[i];
			mface[numFaces].v1 = vert_map[mface[numFaces].v1];
			mface[numFaces].v2 = vert_map[mface[numFaces].v2];
			mface[numFaces].v3 = vert_map[mface[numFaces].v3];
			if(mface[numFaces].v4) {
				mface[numFaces].v4 = vert_map[mface[numFaces].v4];

				test_index_face_maxvert(&mface[numFaces], &result->faceData,
				numFaces, 4, numVerts);
			}
			else
			{
				test_index_face(&mface[numFaces], &result->faceData,
				numFaces, 3);
			}
			if ((amd->mode & MOD_ARR_MOD_ADV_MAT) && (amd->mat_ob & MOD_ARR_SC_MAT_RND) && (ob->totcol>1))
				mface[numFaces].mat_nr = amd->Mem_Mat_Ob.start_cap;
			origindex[numFaces] = ORIGINDEX_NONE;

			numFaces++;
		}

		MEM_freeN(vert_map);
		start_cap->release(start_cap);
	}
	
	if ((mid_cap) && ((amd->dist_mc & MOD_ARR_DIST_SEQ) || (amd->dist_mc & MOD_ARR_DIST_HALF) ||
		((amd->dist_mc & MOD_ARR_DIST_CURVE) && (amd->curve_cap)))) {
		MVert *cap_mvert;
		MEdge *cap_medge;
		MFace *cap_mface;
		int *origindex;
		int *vert_map;
		int capVerts, capEdges, capFaces;
		int inc;

		capVerts = mid_cap->getNumVerts(mid_cap);
		capEdges = mid_cap->getNumEdges(mid_cap);
		capFaces = mid_cap->getNumFaces(mid_cap);
		cap_mvert = mid_cap->getVertArray(mid_cap);
		cap_medge = mid_cap->getEdgeArray(mid_cap);
		cap_mface = mid_cap->getFaceArray(mid_cap);

		vert_map = MEM_callocN(sizeof(*vert_map) * capVerts,
				"arrayModifier_doArray vert_map");
		for(j=0; j < amd->count_mc; j++) {
			
			origindex = result->getVertDataArray(result, CD_ORIGINDEX);
			for(i = 0; i < capVerts; i++) {
				MVert *mv = &cap_mvert[i];
				short merged = 0;

				if(amd->flags & MOD_ARR_MERGE) {
					float tmp_co[3];
					MVert *in_mv;
					int k;

					copy_v3_v3(tmp_co, mv->co);
					mul_m4_v3(mid_offset, tmp_co);
					for(k = 0; k < maxVerts; k++) {
						in_mv = &src_mvert[k];
						/* if this vert is within merge limit, merge */
						if(compare_len_v3v3(tmp_co, in_mv->co, amd->merge_dist)) {
							vert_map[i] = calc_mapping(indexMap, k, 0);
							merged = 1;
							break;
						}
					}
				}

				if(!merged) {
					DM_copy_vert_data(mid_cap, result, i, numVerts, 1);
					mvert[numVerts] = *mv;
					mul_m4_v3(mid_offset, mvert[numVerts].co);
					origindex[numVerts] = ORIGINDEX_NONE;
					vert_map[i] = numVerts;
					numVerts++;
				}
			}

			origindex = result->getEdgeDataArray(result, CD_ORIGINDEX);
			for(i = 0; i < capEdges; i++) {
				int v1, v2;

				v1 = vert_map[cap_medge[i].v1];
				v2 = vert_map[cap_medge[i].v2];

				if(!BLI_edgehash_haskey(edges, v1, v2)) {
					DM_copy_edge_data(mid_cap, result, i, numEdges, 1);
					medge[numEdges] = cap_medge[i];
					medge[numEdges].v1 = v1;
					medge[numEdges].v2 = v2;
					origindex[numEdges] = ORIGINDEX_NONE;

					numEdges++;
				}
			}
			
			origindex = result->getFaceDataArray(result, CD_ORIGINDEX);
			for(i = 0; i < capFaces; i++) {
				DM_copy_face_data(mid_cap, result, i, numFaces, 1);
				mface[numFaces] = cap_mface[i];
				mface[numFaces].v1 = vert_map[mface[numFaces].v1];
				mface[numFaces].v2 = vert_map[mface[numFaces].v2];
				mface[numFaces].v3 = vert_map[mface[numFaces].v3];
				if(mface[numFaces].v4) {
					mface[numFaces].v4 = vert_map[mface[numFaces].v4];

					test_index_face_maxvert(&mface[numFaces], &result->faceData,
					numFaces, 4, numVerts);
				}
				else
				{
					test_index_face(&mface[numFaces], &result->faceData,
					numFaces, 3);
				}
				if ((amd->mode & MOD_ARR_MOD_ADV_MAT) && (amd->mat_ob & MOD_ARR_MC_MAT_RND) && (ob->totcol>1))
					mface[numFaces].mat_nr = amd->Mem_Mat_Ob.mid_cap[j];
				origindex[numFaces] = ORIGINDEX_NONE;
				numFaces++;
			}

			if (amd->dist_mc & MOD_ARR_DIST_SEQ){
				if (amd->mode & MOD_ARR_MOD_ADV) {
					if (amd->Mem_Ob[j+1].transform == 1)
						copy_m4_m4(mid_offset, prec_mid);
				}
				mult_m4_m4m4(tmp_mat, offset, mid_offset);
				copy_m4_m4(mid_offset, tmp_mat);
				if (amd->mode & MOD_ARR_MOD_ADV) {
					if (amd->Mem_Ob[j+1].transform == 1) {
						float app[4][4];
						unit_m4(app);
						loc_eul_size_to_mat4(app, amd->Mem_Ob[j+1].loc, amd->Mem_Ob[j+1].rot, amd->Mem_Ob[j+1].scale);
						copy_m4_m4(prec_mid, mid_offset);
						copy_m4_m4(tmp_mat, mid_offset);
						mult_m4_m4m4(mid_offset, tmp_mat, app);
					}
				}
			}
			else if (amd->dist_mc & MOD_ARR_DIST_HALF){
				mult_m4_m4m4(tmp_mat, half_offset, mid_offset);
				copy_m4_m4(mid_offset, tmp_mat);
			}
			else if (amd->dist_mc & MOD_ARR_DIST_CURVE){
				if (amd->outer_cp & MOD_ARR_CP_FIRST && amd->start_cap)
					inc = 1;
				else 
					inc = 0;
				if (j+1 < amd->count_mc){
					if (nu->bezt)
						copy_v3_v3(mid_offset[3], nu->bezt[j+1+inc].vec[1]);
					else
						copy_v3_v3(mid_offset[3], nu->bp[j+1+inc].vec);
				}
			}
		}
		MEM_freeN(vert_map);
		mid_cap->release(mid_cap);
	}

	if(end_cap) {
		float endoffset[4][4];
		MVert *cap_mvert;
		MEdge *cap_medge;
		MFace *cap_mface;
		int *origindex;
		int *vert_map;
		int capVerts, capEdges, capFaces;

		capVerts = end_cap->getNumVerts(end_cap);
		capEdges = end_cap->getNumEdges(end_cap);
		capFaces = end_cap->getNumFaces(end_cap);
		cap_mvert = end_cap->getVertArray(end_cap);
		cap_medge = end_cap->getEdgeArray(end_cap);
		cap_mface = end_cap->getFaceArray(end_cap);

		mult_m4_m4m4(endoffset, offset, final_offset);
		if (amd->dist_mc & MOD_ARR_DIST_CURVE) {
			if (amd->outer_cp & MOD_ARR_CP_LAST && amd->end_cap && amd->curve_cap) {
				if (nu->bezt)
					copy_v3_v3(endoffset[3], nu->bezt[nu->pntsu - 1].vec[1]);
				else
					copy_v3_v3(endoffset[3], nu->bp[nu->pntsu * nu->pntsv - 1].vec);
			}
		}
		vert_map = MEM_callocN(sizeof(*vert_map) * capVerts,
		"arrayModifier_doArray vert_map");

		origindex = result->getVertDataArray(result, CD_ORIGINDEX);
		for(i = 0; i < capVerts; i++) {
			MVert *mv = &cap_mvert[i];
			short merged = 0;

			if(amd->flags & MOD_ARR_MERGE) {
				float tmp_co[3];
				MVert *in_mv;
				int j;

				copy_v3_v3(tmp_co, mv->co);
				mul_m4_v3(offset, tmp_co);

				for(j = 0; j < maxVerts; j++) {
					in_mv = &src_mvert[j];
					/* if this vert is within merge limit, merge */
					if(compare_len_v3v3(tmp_co, in_mv->co, amd->merge_dist)) {
						vert_map[i] = calc_mapping(indexMap, j, count - 1);
						merged = 1;
						break;
					}
				}
			}

			if(!merged) {
				DM_copy_vert_data(end_cap, result, i, numVerts, 1);
				mvert[numVerts] = *mv;
				mul_m4_v3(endoffset, mvert[numVerts].co);
				origindex[numVerts] = ORIGINDEX_NONE;

				vert_map[i] = numVerts;

				numVerts++;
			}
		}

		origindex = result->getEdgeDataArray(result, CD_ORIGINDEX);
		for(i = 0; i < capEdges; i++) {
			int v1, v2;

			v1 = vert_map[cap_medge[i].v1];
			v2 = vert_map[cap_medge[i].v2];

			if(!BLI_edgehash_haskey(edges, v1, v2)) {
				DM_copy_edge_data(end_cap, result, i, numEdges, 1);
				medge[numEdges] = cap_medge[i];
				medge[numEdges].v1 = v1;
				medge[numEdges].v2 = v2;
				origindex[numEdges] = ORIGINDEX_NONE;

				numEdges++;
			}
		}

		origindex = result->getFaceDataArray(result, CD_ORIGINDEX);
		for(i = 0; i < capFaces; i++) {
			DM_copy_face_data(end_cap, result, i, numFaces, 1);
			mface[numFaces] = cap_mface[i];
			mface[numFaces].v1 = vert_map[mface[numFaces].v1];
			mface[numFaces].v2 = vert_map[mface[numFaces].v2];
			mface[numFaces].v3 = vert_map[mface[numFaces].v3];
			if(mface[numFaces].v4) {
				mface[numFaces].v4 = vert_map[mface[numFaces].v4];

				test_index_face(&mface[numFaces], &result->faceData,
				numFaces, 4);
			}
			else
			{
				test_index_face(&mface[numFaces], &result->faceData,
				numFaces, 3);
			}
			if ((amd->mode & MOD_ARR_MOD_ADV_MAT) && (amd->mat_ob & MOD_ARR_EC_MAT_RND) && (ob->totcol>1))
				mface[numFaces].mat_nr = amd->Mem_Mat_Ob.end_cap;
			origindex[numFaces] = ORIGINDEX_NONE;

			numFaces++;
		}

		MEM_freeN(vert_map);
		end_cap->release(end_cap);
	}

	/*if ((amd->type & MOD_ARR_MOD_CURVE) && (amd->curve_ob)) {
		if (amd->dist_cu & MOD_ARR_DIST_EVENLY){
			array_to_curve(scene, amd->curve_ob, ob, mvert, numVerts);
		}
	}*/

	BLI_edgehash_free(edges, NULL);
	MEM_freeN(indexMap);

	CDDM_lower_num_verts(result, numVerts);
	CDDM_lower_num_edges(result, numEdges);
	CDDM_lower_num_faces(result, numFaces);

	return result;
}

static DerivedMesh *applyModifier(ModifierData *md, Object *ob,
						DerivedMesh *dm,
						int UNUSED(useRenderParams),
						int UNUSED(isFinalCalc))
{
	DerivedMesh *result;
	ArrayModifierData *amd = (ArrayModifierData*) md;

	result = arrayModifier_doArray(amd, md->scene, ob, dm, 0);

	if(result != dm) {
		CDDM_calc_normals(result);
		
		if(amd->arr_group!=NULL) {	
			ob->transflag = OB_DUPLIARRAY;
			ob->dup_group = amd->arr_group;
		}
		else {
			ob->transflag = 0;
			ob->dup_group = NULL;
		}
	}
	return result;
}

static DerivedMesh *applyModifierEM(ModifierData *md, Object *ob,
						struct EditMesh *UNUSED(editData),
						DerivedMesh *dm)
{
	return applyModifier(md, ob, dm, 0, 1);
}

static void freeData(ModifierData *md)
{
	ArrayModifierData *amd = (ArrayModifierData*) md;
	
	if (amd) 
	{
		if (amd->Mem_Ob)
			MEM_freeN(amd->Mem_Ob);
	}
}

ModifierTypeInfo modifierType_Array = {
	/* name */              "Array",
	/* structName */        "ArrayModifierData",
	/* structSize */        sizeof(ArrayModifierData),
	/* type */              eModifierTypeType_Constructive,
	/* flags */             eModifierTypeFlag_AcceptsMesh
							| eModifierTypeFlag_SupportsMapping
							| eModifierTypeFlag_SupportsEditmode
							| eModifierTypeFlag_EnableInEditmode
							| eModifierTypeFlag_AcceptsCVs,

	/* copyData */          copyData,
	/* deformVerts */       NULL,
	/* deformMatrices */    NULL,
	/* deformVertsEM */     NULL,
	/* deformMatricesEM */  NULL,
	/* applyModifier */     applyModifier,
	/* applyModifierEM */   applyModifierEM,
	/* initData */          initData,
	/* requiredDataMask */  NULL,
	/* freeData */          freeData,
	/* isDisabled */        NULL,
	/* updateDepgraph */    updateDepgraph,
	/* dependsOnTime */     NULL,
	/* dependsOnNormals */	NULL,
	/* foreachObjectLink */ foreachObjectLink,
	/* foreachIDLink */     NULL,
};
