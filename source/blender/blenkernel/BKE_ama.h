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


float length_fitcurve(struct ArrayModifierData *amd, struct Scene *scene);
int length_to_count(float length, const float offset[3]);
float count_to_length(int count, const float offset[3]);
float f_rand_max(float max);
void array_scale_offset(const float max_off[3], float rit[3],int prop);
void array_offset(const float max_off[3], float rit[3],int sign);
void init_offset(const int start, const int end, struct ArrayModifierData *ar);
void create_offset(const int n, const int totmat, struct ArrayModifierData *ar, struct Object *ob);
void array_to_curve(struct Scene *scene, struct Object *cuOb, struct Object *target, float (*vertexCos)[3], int numVerts);

#endif

