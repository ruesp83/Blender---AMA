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
 * along with this program; if not, write to the Free Software Foundation,
 * Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 *
 * The Original Code is Copyright (C) Blender Foundation.
 * All rights reserved.
 *
 * The Original Code is: all of this file.
 *
 * Contributor(s): none yet.
 *
 * ***** END GPL LICENSE BLOCK *****
 */
#ifndef BKE_AMA_H
#define BKE_AMA_H

/** \file BKE_ama.h
 *  \ingroup bke
 */

struct ArrayModifierData;
struct Object;
struct Scene;
struct MVert;
struct MFace;
struct CustomData;

/*typedef struct IndexMapEntry {
 *	// the new vert index that this old vert index maps to 
 *	int new;
 *	// -1 if this vert isn't merged, otherwise the old vert index it
 *	//  should be replaced with
 *	int merge;
 *	// 1 if this vert's first copy is merged with the last copy of its
 *	// merge target, otherwise 0
 *	short merge_final;
 *} IndexMapEntry;*/


/* calculations is in local space of deformed object
 * so we store in latmat transform from path coord inside object 
 */
typedef struct {
	float dmin[3], dmax[3], dsize, dloc[3];
	float curvespace[4][4], objectspace[4][4], objectspace3[3][3];
	int no_rot_axis;
} CurveDeform;

/*typedef struct Temp {
 *	float loc[3];
 *	float rot[3];
 *	float scale[3];
 *	int seed;
 *} Temp;*/

float BKE_vertarray_size(struct MVert *mvert, int numVerts, int axis);

float BKE_length_fitcurve(struct ArrayModifierData *amd, struct Scene *scene);
int BKE_length_to_count(float length, const float offset[3]);
float BKE_count_to_length(int count, const float offset[3]);
//float rand_max_fl(float max);

//void array_scale_offset(const float max_off[3], float rit[3], int prop, int sign, int seed);
//void array_offset(const float max_off[3], float rit[3], int sign, int seed);
void BKE_init_mat_object_cap(const int start, const int end, int *vet_mc);
void BKE_init_offset(const int start, const int end, struct ArrayModifierData *ar);
void BKE_create_offset(const int n, const int totmat, struct ArrayModifierData *ar, struct Object *ob);
void BKE_reset_offset(float offset[4][4], const int count, float loc[3], float rot[3], float scale[3]);
//void array_to_curve(struct Scene *scene, struct Object *cuOb, float (*vertexCos)[3], int numVerts);
void array_to_curve(struct Scene *scene, struct Object *cuOb, struct Object *target, float *vertexCos, float *vec, 
					float *cent);

#endif
