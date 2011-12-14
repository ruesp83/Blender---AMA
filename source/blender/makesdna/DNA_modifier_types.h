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
 * ***** END GPL LICENSE BLOCK *****
 */

#ifndef DNA_MODIFIER_TYPES_H
#define DNA_MODIFIER_TYPES_H

/** \file DNA_modifier_types.h
 *  \ingroup DNA
 */

#include "DNA_listBase.h"


#define MODSTACK_DEBUG 1

/* WARNING ALERT! TYPEDEF VALUES ARE WRITTEN IN FILES! SO DO NOT CHANGE!
 * (ONLY ADD NEW ITEMS AT THE END) */

typedef enum ModifierType {
	eModifierType_None = 0,
	eModifierType_Subsurf,
	eModifierType_Lattice,
	eModifierType_Curve,
	eModifierType_Build,
	eModifierType_Mirror,
	eModifierType_Decimate,
	eModifierType_Wave,
	eModifierType_Armature,
	eModifierType_Hook,
	eModifierType_Softbody,
	eModifierType_Boolean,
	eModifierType_Array,
	eModifierType_EdgeSplit,
	eModifierType_Displace,
	eModifierType_UVProject,
	eModifierType_Smooth,
	eModifierType_Cast,
	eModifierType_MeshDeform,
	eModifierType_ParticleSystem,
	eModifierType_ParticleInstance,
	eModifierType_Explode,
	eModifierType_Cloth,
	eModifierType_Collision,
	eModifierType_Bevel,
	eModifierType_Shrinkwrap,
	eModifierType_Fluidsim,
	eModifierType_Mask,
	eModifierType_SimpleDeform,
	eModifierType_Multires,
	eModifierType_Surface,
	eModifierType_Smoke,
	eModifierType_ShapeKey,
	eModifierType_Solidify,
	eModifierType_Screw,
	eModifierType_Warp,
	eModifierType_WeightVGEdit,
	eModifierType_WeightVGMix,
	eModifierType_WeightVGProximity,
	eModifierType_Ocean,
	eModifierType_DynamicPaint,
	NUM_MODIFIER_TYPES
} ModifierType;

typedef enum ModifierMode {
	eModifierMode_Realtime = (1<<0),
	eModifierMode_Render = (1<<1),
	eModifierMode_Editmode = (1<<2),
	eModifierMode_OnCage = (1<<3),
	eModifierMode_Expanded = (1<<4),
	eModifierMode_Virtual = (1<<5),
	eModifierMode_ApplyOnSpline = (1<<6),
	eModifierMode_DisableTemporary = (1 << 31)
} ModifierMode;

typedef struct ModifierData {
	struct ModifierData *next, *prev;

	int type, mode;
	int stackindex, pad;
	char name[32];
	
	/* XXX for timing info set by caller... solve later? (ton) */
	struct Scene *scene;
	
	char *error;
} ModifierData;

typedef enum {
	eSubsurfModifierFlag_Incremental = (1<<0),
	eSubsurfModifierFlag_DebugIncr = (1<<1),
	eSubsurfModifierFlag_ControlEdges = (1<<2),
	eSubsurfModifierFlag_SubsurfUv = (1<<3)
} SubsurfModifierFlag;

/* not a real modifier */
typedef struct MappingInfoModifierData {
	ModifierData modifier;

	struct Tex *texture;
	struct Object *map_object;
	char uvlayer_name[32];
	int uvlayer_tmp;
	int texmapping;
} MappingInfoModifierData;

typedef struct SubsurfModifierData {
	ModifierData modifier;

	short subdivType, levels, renderLevels, flags;

	void *emCache, *mCache;
} SubsurfModifierData;

typedef struct LatticeModifierData {
	ModifierData modifier;

	struct Object *object;
	char name[32];			/* optional vertexgroup name */
} LatticeModifierData;

typedef struct CurveModifierData {
	ModifierData modifier;

	struct Object *object;
	char name[32];			/* optional vertexgroup name */
	short defaxis;			/* axis along which curve deforms */
	char pad[6];
} CurveModifierData;

/* CurveModifierData->defaxis */
#define MOD_CURVE_POSX	1
#define MOD_CURVE_POSY	2
#define MOD_CURVE_POSZ	3
#define MOD_CURVE_NEGX	4
#define MOD_CURVE_NEGY	5
#define MOD_CURVE_NEGZ	6

typedef struct BuildModifierData {
	ModifierData modifier;

	float start, length;
	int randomize, seed;
} BuildModifierData;

/* Mask Modifier */
typedef struct MaskModifierData {
	ModifierData modifier;
	
	struct Object *ob_arm;	/* armature to use to in place of hardcoded vgroup */
	char vgroup[32];		/* name of vertex group to use to mask */
	
	int mode;				/* using armature or hardcoded vgroup */
	int flag;				/* flags for various things */
} MaskModifierData;

/* Mask Modifier -> mode */
#define MOD_MASK_MODE_VGROUP		0
#define MOD_MASK_MODE_ARM			1

/* Mask Modifier -> flag */
#define MOD_MASK_INV			(1<<0)


typedef struct ArrayChangeObject {
	//float location[4][4];
	float rot[3];
	float scale[3];
	float loc[3];
	/*0 not modified
	  1 modified
	*/
	float cu_cent[3];
	float cu_loc[4];
	int transform;
	int id_mat;
	int rand_group_obj;
} ArrayChangeObject;

/*typedef struct ObjectCap {
	struct Object *ob_cap;
	struct Object *curve;
	int count;
	int distribution;
} ObjectCap;*/

typedef struct Mat_ObCap {
	int *mid_cap;
	int start_cap;
	int end_cap;
} Mat_ObCap;

typedef struct ArrayModifierData {
	ModifierData modifier;

	/* the object with which to cap the start of the array  */
	struct Object *start_cap;
	/* the object with which to cap the mid of the array  */
	struct Object *mid_cap;
	/* the object with which to cap the end of the array  */
	struct Object *end_cap;
	/* the curve object to use for MOD_ARR_MOD_CURVE */
	struct Object *curve_ob;
	/* the curve object to use for Object Cap */
	struct Object *curve_cap;
	/* the object to use for object offset */
	struct Object *offset_ob;
	/* stores information about objects subject to random */
	struct ArrayChangeObject *Mem_Ob;
	/* stores the material of the "Object Cap" */
	struct Mat_ObCap Mem_Mat_Ob;
	/* group linked to the modifier */
	struct Group *arr_group;
	/* a constant duplicate offset;
	   1 means the duplicates are 1 unit apart
	*/
	float offset[3];
	/* a scaled factor for duplicate offsets;
	   1 means the duplicates are 1 object-width apart
	*/
	float scale[3];
	/* Offset for clone group */
	float delta[4][4];
	/* the length over which to distribute the duplicates */
	float length;
	/* the limit below which to merge vertices in adjacent duplicates */
	float merge_dist;
	/* randomize Offset */
	float loc_offset[3];
	float rot_offset[3];
	float scale_offset[3];
	/* determines how duplicate count is calculated; one of:
		  MOD_ARR_FIXEDCOUNT -> fixed
		  MOD_ARR_FITLENGTH  -> calculated to fit a set length
		  MOD_ARR_FITBETWEEN   -> Number of duplicates between two objects
	*/
	int fit_type;
	/* flags specifying how total offset is calculated; binary OR of:
		 MOD_ARR_OFF_CONST    -> total offset += offset
		 MOD_ARR_OFF_RELATIVE -> total offset += relative * object width
		 MOD_ARR_OFF_OBJ      -> total offset += offset_ob's matrix
	   total offset is the sum of the individual enabled offsets
	*/
	int offset_type;
	/* general flags:
		  MOD_ARR_MERGE -> merge vertices in adjacent duplicates
	*/
	int flags;
	/* the number of duplicates to generate for MOD_ARR_FIXEDCOUNT */
	int count;
	/* the number of duplicates of Mid Cap and type of distribution */
	int count_mc;
	int dist_mc;
	/* Indicates whether you use a random material or the original*/
	int mat_ob;
	/* Indicates whether to use the Start Cap and End Cap on the curve of the Mid Cap */
	int outer_cp;
	/* number of rays and the direction of the Clones */
	int rays;
	int rays_dir;
	/* Number of copies of the same material */
	int cont_mat;
	/* randomizes the materials */
	int rand_mat;
	/* lock the noise offset */
	int lock;
	/* normal Mode-Advanced Mode */
	int mode;
	/* indicates how the modifier should be used */
	int type;
	/* direction Offset */
	int sign;
	/* keeps the ratio on the scale */
	int flag_offset;
	/* how to distribute the clones on a curve */
	int dist_cu;
	/* ability to randomization of objects belonging to the group linked */
	int rand_group;
	//int pad1;
} ArrayModifierData;

/* ArrayModifierData->fit_type */
#define MOD_ARR_FIXEDCOUNT		0
#define MOD_ARR_FITLENGTH		1
#define MOD_ARR_FITBETWEEN		2

/* ArrayModifierData->offset_type */
#define MOD_ARR_OFF_CONST		(1<<0)
#define MOD_ARR_OFF_RELATIVE	(1<<1)
#define MOD_ARR_OFF_OBJ			(1<<2)

/* ArrayModifierData->flags */
#define MOD_ARR_MERGE      		(1<<0)
#define MOD_ARR_MERGEFINAL 		(1<<1)

/* ArrayModifierData->mat_ob */
/* Array */
#define MOD_ARR_AR_MAT_RND     	(1<<0)
/* Start Cap */
#define MOD_ARR_SC_MAT_RND     	(1<<1)
#define MOD_ARR_SC_MAT_OR 		(1<<2)
/* Mid Cap */
#define MOD_ARR_MC_MAT_RND     	(1<<3)
#define MOD_ARR_MC_MAT_OR 		(1<<4)
/* End Cap */
#define MOD_ARR_EC_MAT_RND     	(1<<5)
#define MOD_ARR_EC_MAT_OR 		(1<<6)

/* ArrayModifierData->type */
#define MOD_ARR_MOD_NRM			(1<<0)
#define MOD_ARR_MOD_CURVE		(1<<1)

/* ArrayModifierData->outer_cp */
#define MOD_ARR_CP_FIRST		(1<<0)
#define MOD_ARR_CP_LAST			(1<<1)

/* ArrayModifierData->mode */
#define MOD_ARR_MOD_ADV			(1<<0)
#define MOD_ARR_MOD_ADV_MAT		(1<<1)
#define MOD_ARR_MOD_ADV_MID		(1<<3)
#define MOD_ARR_MOD_ADV_CLONE	(1<<4)

/* ArrayModifierData->sign */
#define MOD_ARR_SIGN_P			(1<<0)
#define MOD_ARR_SIGN_L			(1<<1)

/* ArrayModifierData->lock */
#define MOD_ARR_LOCK			(1<<0)

/* ArrayModifierData->flag_offset */
#define MOD_ARR_PROP			(1<<0)
#define MOD_ARR_LOCAL			(1<<1)

/* ArrayModifierData->rand_mat */
#define MOD_ARR_MAT				(1<<0)
#define MOD_ARR_SEQ				(1<<1)

/* ArrayModifierData->rays_dir */
#define MOD_ARR_RAYS_X			0
#define MOD_ARR_RAYS_Y			1
#define MOD_ARR_RAYS_Z			2

/* ArrayModifierData->rand_group */
#define MOD_ARR_RAND_GROUP		(1<<0)
#define MOD_ARR_RAND_MAT_GROUP  (1<<1)

/* ArrayModifierData->dist_cu */
#define MOD_ARR_DIST_EVENLY		(1<<0)
#define MOD_ARR_DIST_SEGMENT	(1<<1)

/* ArrayModifierData->dist_mc */
#define MOD_ARR_DIST_SEQ		(1<<0)
#define MOD_ARR_DIST_HALF		(1<<1)
#define MOD_ARR_DIST_CURVE		(1<<2)

typedef struct MirrorModifierData {
	ModifierData modifier;

	short axis; /* deprecated, use flag instead */
	short flag;
	float tolerance;
	struct Object *mirror_ob;
} MirrorModifierData;

/* MirrorModifierData->flag */
#define MOD_MIR_CLIPPING	(1<<0)
#define MOD_MIR_MIRROR_U	(1<<1)
#define MOD_MIR_MIRROR_V	(1<<2)
#define MOD_MIR_AXIS_X		(1<<3)
#define MOD_MIR_AXIS_Y		(1<<4)
#define MOD_MIR_AXIS_Z		(1<<5)
#define MOD_MIR_VGROUP		(1<<6)
#define MOD_MIR_NO_MERGE	(1<<7)

typedef struct EdgeSplitModifierData {
	ModifierData modifier;

	float split_angle;    /* angle above which edges should be split */
	int flags;
} EdgeSplitModifierData;

/* EdgeSplitModifierData->flags */
#define MOD_EDGESPLIT_FROMANGLE   (1<<1)
#define MOD_EDGESPLIT_FROMFLAG    (1<<2)

typedef struct BevelModifierData {
	ModifierData modifier;

	float value;          /* the "raw" bevel value (distance/amount to bevel) */
	int res;              /* the resolution (as originally coded, it is the number of recursive bevels) */
	int pad;
	short flags;          /* general option flags */
	short val_flags;      /* flags used to interpret the bevel value */
	short lim_flags;      /* flags to tell the tool how to limit the bevel */
	short e_flags;        /* flags to direct how edge weights are applied to verts */
	float bevel_angle;    /* if the BME_BEVEL_ANGLE is set, this will be how "sharp" an edge must be before it gets beveled */
	char defgrp_name[32]; /* if the BME_BEVEL_VWEIGHT option is set, this will be the name of the vert group */
} BevelModifierData;

typedef struct BMeshModifierData {
	ModifierData modifier;

	float pad;
	int type;
} BMeshModifierData;


/* Smoke modifier flags */
#define MOD_SMOKE_TYPE_DOMAIN (1 << 0)
#define MOD_SMOKE_TYPE_FLOW (1 << 1)
#define MOD_SMOKE_TYPE_COLL (1 << 2)

typedef struct SmokeModifierData {
	ModifierData modifier;

	struct SmokeDomainSettings *domain;
	struct SmokeFlowSettings *flow; /* inflow, outflow, smoke objects */
	struct SmokeCollSettings *coll; /* collision objects */
	float time;
	int type;  /* domain, inflow, outflow, ... */
} SmokeModifierData;

typedef struct DisplaceModifierData {
	ModifierData modifier;

	/* keep in sync with MappingInfoModifierData */
	struct Tex *texture;
	struct Object *map_object;
	char uvlayer_name[32];
	int uvlayer_tmp;
	int texmapping;
	int pad10;
	/* end MappingInfoModifierData */

	float strength;
	int direction;
	char defgrp_name[32];
	float midlevel;
} DisplaceModifierData;

/* DisplaceModifierData->direction */
enum {
	MOD_DISP_DIR_X,
	MOD_DISP_DIR_Y,
	MOD_DISP_DIR_Z,
	MOD_DISP_DIR_NOR,
	MOD_DISP_DIR_RGB_XYZ,
};

/* DisplaceModifierData->texmapping */
enum {
	MOD_DISP_MAP_LOCAL,
	MOD_DISP_MAP_GLOBAL,
	MOD_DISP_MAP_OBJECT,
	MOD_DISP_MAP_UV,
};

typedef struct UVProjectModifierData {
	ModifierData modifier;

	/* the objects which do the projecting */
	struct Object *projectors[10]; /* MOD_UVPROJECT_MAX */
	struct Image *image;      /* the image to project */
	int flags;
	int num_projectors;
	float aspectx, aspecty;
	float scalex, scaley;												
	char uvlayer_name[32];
	int uvlayer_tmp, pad;
} UVProjectModifierData;

#define MOD_UVPROJECT_MAXPROJECTORS 10

/* UVProjectModifierData->flags */
#define MOD_UVPROJECT_OVERRIDEIMAGE (1<<0)

typedef struct DecimateModifierData {
	ModifierData modifier;

	float percent;
	int faceCount;
} DecimateModifierData;

/* Smooth modifier flags */
#define MOD_SMOOTH_X (1<<1)
#define MOD_SMOOTH_Y (1<<2)
#define MOD_SMOOTH_Z (1<<3)

typedef struct SmoothModifierData {
	ModifierData modifier;
	float fac;
	char defgrp_name[32];
	short flag, repeat;

} SmoothModifierData;

/* Cast modifier flags */
#define MOD_CAST_X (1<<1)
#define MOD_CAST_Y (1<<2)
#define MOD_CAST_Z (1<<3)
#define MOD_CAST_USE_OB_TRANSFORM (1<<4)
#define MOD_CAST_SIZE_FROM_RADIUS (1<<5)

/* Cast modifier projection types */
#define MOD_CAST_TYPE_SPHERE 0
#define MOD_CAST_TYPE_CYLINDER 1
#define MOD_CAST_TYPE_CUBOID 2

typedef struct CastModifierData {
	ModifierData modifier;

	struct Object *object;
	float fac;
	float radius;
	float size;
	char defgrp_name[32];
	short flag, type;
} CastModifierData;

enum {
	MOD_WAV_MAP_LOCAL,
	MOD_WAV_MAP_GLOBAL,
	MOD_WAV_MAP_OBJECT,
	MOD_WAV_MAP_UV,
};

/* WaveModifierData.flag */
#define MOD_WAVE_X      (1<<1)
#define MOD_WAVE_Y      (1<<2)
#define MOD_WAVE_CYCL   (1<<3)
#define MOD_WAVE_NORM   (1<<4)
#define MOD_WAVE_NORM_X (1<<5)
#define MOD_WAVE_NORM_Y (1<<6)
#define MOD_WAVE_NORM_Z (1<<7)

typedef struct WaveModifierData {
	ModifierData modifier;

	struct Object *objectcenter;
	char defgrp_name[32];
	struct Tex *texture;
	struct Object *map_object;

	short flag, pad;

	float startx, starty, height, width;
	float narrow, speed, damp, falloff;

	int texmapping, uvlayer_tmp;

	char uvlayer_name[32];

	float timeoffs, lifetime;
	float pad1;
} WaveModifierData;

typedef struct ArmatureModifierData {
	ModifierData modifier;

	short deformflag, multi;		/* deformflag replaces armature->deformflag */
	int pad2;
	struct Object *object;
	float *prevCos;		/* stored input of previous modifier, for vertexgroup blending */
	char defgrp_name[32];
} ArmatureModifierData;

typedef struct HookModifierData {
	ModifierData modifier;

	struct Object *object;
	char subtarget[32];		/* optional name of bone target */
	
	float parentinv[4][4];	/* matrix making current transform unmodified */
	float cent[3];			/* visualization of hook */
	float falloff;			/* if not zero, falloff is distance where influence zero */
	
	int *indexar;			/* if NULL, it's using vertexgroup */
	int totindex;
	float force;
	char name[32];			/* optional vertexgroup name */
} HookModifierData;

typedef struct SoftbodyModifierData {
	ModifierData modifier;
} SoftbodyModifierData;

typedef struct ClothModifierData {
	ModifierData		modifier;

	struct Scene *scene;			/* the context, time etc is here */
	struct Cloth *clothObject;		/* The internal data structure for cloth. */
	struct ClothSimSettings *sim_parms; /* definition is in DNA_cloth_types.h */
	struct ClothCollSettings *coll_parms; /* definition is in DNA_cloth_types.h */
	struct PointCache *point_cache;	/* definition is in DNA_object_force.h */
	struct ListBase ptcaches;
} ClothModifierData;

typedef struct CollisionModifierData {
	ModifierData	modifier;
	
	struct MVert *x; /* position at the beginning of the frame */
	struct MVert *xnew; /* position at the end of the frame */
	struct MVert *xold; /* unsued atm, but was discussed during sprint */
	struct MVert *current_xnew; /* new position at the actual inter-frame step */
	struct MVert *current_x; /* position at the actual inter-frame step */
	struct MVert *current_v; /* (xnew - x) at the actual inter-frame step */
	
	struct MFace *mfaces; /* object face data */
	
	unsigned int numverts;
	unsigned int numfaces;
	float time_x, time_xnew;		/* cfra time of modifier */
	struct BVHTree *bvhtree; /* bounding volume hierarchy for this cloth object */
} CollisionModifierData;

typedef struct SurfaceModifierData {
	ModifierData	modifier;

	struct MVert *x; /* old position */
	struct MVert *v; /* velocity */

	struct DerivedMesh *dm;

	struct BVHTreeFromMesh *bvhtree; /* bounding volume hierarchy of the mesh faces */

	int cfra, numverts;
} SurfaceModifierData;

typedef enum {
	eBooleanModifierOp_Intersect,
	eBooleanModifierOp_Union,
	eBooleanModifierOp_Difference,
} BooleanModifierOp;
typedef struct BooleanModifierData {
	ModifierData modifier;

	struct Object *object;
	int operation, pad;
} BooleanModifierData;

#define MOD_MDEF_INVERT_VGROUP	(1<<0)
#define MOD_MDEF_DYNAMIC_BIND	(1<<1)

#define MOD_MDEF_VOLUME			0
#define MOD_MDEF_SURFACE		1

typedef struct MDefInfluence {
	int vertex;
	float weight;
} MDefInfluence;

typedef struct MDefCell {
	int offset;
	int totinfluence;
} MDefCell;

typedef struct MeshDeformModifierData {
	ModifierData modifier;

	struct Object *object;			/* mesh object */
	char defgrp_name[32];			/* optional vertexgroup name */

	short gridsize, flag, mode, pad;

	/* result of static binding */
	MDefInfluence *bindinfluences;	/* influences */
	int *bindoffsets;				/* offsets into influences array */
	float *bindcagecos;				/* coordinates that cage was bound with */
	int totvert, totcagevert;		/* total vertices in mesh and cage */

	/* result of dynamic binding */
	MDefCell *dyngrid;				/* grid with dynamic binding cell points */
	MDefInfluence *dyninfluences;	/* dynamic binding vertex influences */
	int *dynverts, *pad2;			/* is this vertex bound or not? */
	int dyngridsize;				/* size of the dynamic bind grid */
	int totinfluence;				/* total number of vertex influences */
	float dyncellmin[3];			/* offset of the dynamic bind grid */
	float dyncellwidth;				/* width of dynamic bind cell */
	float bindmat[4][4];			/* matrix of cage at binding time */

	/* deprecated storage */
	float *bindweights;				/* deprecated inefficient storage */
	float *bindcos;					/* deprecated storage of cage coords */

	/* runtime */
	void (*bindfunc)(struct Scene *scene,
		struct MeshDeformModifierData *mmd,
		float *vertexcos, int totvert, float cagemat[][4]);
} MeshDeformModifierData;

typedef enum {
	eParticleSystemFlag_Pars =			(1<<0),
	eParticleSystemFlag_psys_updated =	(1<<1),
	eParticleSystemFlag_file_loaded =	(1<<2),
} ParticleSystemModifierFlag;

typedef struct ParticleSystemModifierData {
	ModifierData modifier;
	struct ParticleSystem *psys;
	struct DerivedMesh *dm;
	int totdmvert, totdmedge, totdmface;
	short flag, rt;
} ParticleSystemModifierData;

typedef enum {
	eParticleInstanceFlag_Parents =		(1<<0),
	eParticleInstanceFlag_Children =	(1<<1),
	eParticleInstanceFlag_Path =		(1<<2),
	eParticleInstanceFlag_Unborn =		(1<<3),
	eParticleInstanceFlag_Alive =		(1<<4),
	eParticleInstanceFlag_Dead =		(1<<5),
	eParticleInstanceFlag_KeepShape =	(1<<6),
	eParticleInstanceFlag_UseSize =		(1<<7),
} ParticleInstanceModifierFlag;

typedef struct ParticleInstanceModifierData {
	ModifierData modifier;
	struct Object *ob;
	short psys, flag, axis, rt;
	float position, random_position;
} ParticleInstanceModifierData;

typedef enum {
	eExplodeFlag_CalcFaces =	(1<<0),
	eExplodeFlag_PaSize =		(1<<1),
	eExplodeFlag_EdgeCut =		(1<<2),
	eExplodeFlag_Unborn =		(1<<3),
	eExplodeFlag_Alive =		(1<<4),
	eExplodeFlag_Dead =			(1<<5),
} ExplodeModifierFlag;

typedef struct ExplodeModifierData {
	ModifierData modifier;
	int *facepa;
	short flag, vgroup;
	float protect;
	char uvname[32];
} ExplodeModifierData;

typedef struct MultiresModifierData {
	ModifierData modifier;

	char lvl, sculptlvl, renderlvl, totlvl;
	char simple, flags, pad[2];
} MultiresModifierData;

typedef enum {
	eMultiresModifierFlag_ControlEdges = (1<<0),
	eMultiresModifierFlag_PlainUv = (1<<1),
} MultiresModifierFlag;

typedef struct FluidsimModifierData {
	ModifierData modifier;
	
	struct FluidsimSettings *fss; /* definition is in DNA_object_fluidsim.h */
	struct PointCache *point_cache;	/* definition is in DNA_object_force.h */
} FluidsimModifierData;

typedef struct ShrinkwrapModifierData {
	ModifierData modifier;

	struct Object *target;	/* shrink target */
	struct Object *auxTarget; /* additional shrink target */
	char vgroup_name[32];	/* optional vertexgroup name */
	float keepDist;			/* distance offset to keep from mesh/projection point */
	short shrinkType;		/* shrink type projection */
	short shrinkOpts;		/* shrink options */
	char projAxis;			/* axis to project over */

	/*
	 * if using projection over vertex normal this controls the
	 * the level of subsurface that must be done before getting the
	 * vertex coordinates and normal
	 */
	char subsurfLevels;

	char pad[6];

} ShrinkwrapModifierData;

/* Shrinkwrap->shrinkType */
#define MOD_SHRINKWRAP_NEAREST_SURFACE	0
#define MOD_SHRINKWRAP_PROJECT			1
#define MOD_SHRINKWRAP_NEAREST_VERTEX	2

/* Shrinkwrap->shrinkOpts */
#define MOD_SHRINKWRAP_PROJECT_ALLOW_POS_DIR	(1<<0)	/* allow shrinkwrap to move the vertex in the positive direction of axis */
#define MOD_SHRINKWRAP_PROJECT_ALLOW_NEG_DIR	(1<<1)	/* allow shrinkwrap to move the vertex in the negative direction of axis */

#define MOD_SHRINKWRAP_CULL_TARGET_FRONTFACE	(1<<3)	/* ignore vertex moves if a vertex ends projected on a front face of the target */
#define MOD_SHRINKWRAP_CULL_TARGET_BACKFACE		(1<<4)	/* ignore vertex moves if a vertex ends projected on a back face of the target */

#define MOD_SHRINKWRAP_KEEP_ABOVE_SURFACE		(1<<5)	/* distance is measure to the front face of the target */

#define MOD_SHRINKWRAP_PROJECT_OVER_X_AXIS		(1<<0)
#define MOD_SHRINKWRAP_PROJECT_OVER_Y_AXIS		(1<<1)
#define MOD_SHRINKWRAP_PROJECT_OVER_Z_AXIS		(1<<2)
#define MOD_SHRINKWRAP_PROJECT_OVER_NORMAL			0	/* projection over normal is used if no axis is selected */

typedef struct SimpleDeformModifierData {
	ModifierData modifier;

	struct Object *origin;	/* object to control the origin of modifier space coordinates */
	char vgroup_name[32];	/* optional vertexgroup name */
	float factor;			/* factors to control simple deforms */
	float limit[2];			/* lower and upper limit */		

	char mode;				/* deform function */
	char axis;				/* lock axis (for taper and strech) */
	char originOpts;		/* originOptions */
	char pad;

} SimpleDeformModifierData;

#define MOD_SIMPLEDEFORM_MODE_TWIST		1
#define MOD_SIMPLEDEFORM_MODE_BEND		2
#define MOD_SIMPLEDEFORM_MODE_TAPER		3
#define MOD_SIMPLEDEFORM_MODE_STRETCH	4

#define MOD_SIMPLEDEFORM_LOCK_AXIS_X			(1<<0)
#define MOD_SIMPLEDEFORM_LOCK_AXIS_Y			(1<<1)

/* indicates whether simple deform should use the local
   coordinates or global coordinates of origin */
#define MOD_SIMPLEDEFORM_ORIGIN_LOCAL			(1<<0)

#define MOD_UVPROJECT_MAX				10

typedef struct ShapeKeyModifierData {
	ModifierData modifier;
} ShapeKeyModifierData;

typedef struct SolidifyModifierData {
	ModifierData modifier;

	char defgrp_name[32];	/* name of vertex group to use */
	float offset;			/* new surface offset level*/
	float offset_fac;		/* midpoint of the offset  */
	float offset_fac_vg;	/* factor for the minimum weight to use when vgroups are used, avoids 0.0 weights giving duplicate geometry */
	float crease_inner;
	float crease_outer;
	float crease_rim;
	int flag;
	short mat_ofs;
	short mat_ofs_rim;
} SolidifyModifierData;

#define MOD_SOLIDIFY_RIM			(1<<0)
#define MOD_SOLIDIFY_EVEN			(1<<1)
#define MOD_SOLIDIFY_NORMAL_CALC	(1<<2)
#define MOD_SOLIDIFY_VGROUP_INV		(1<<3)
#define MOD_SOLIDIFY_RIM_MATERIAL	(1<<4) /* deprecated, used in do_versions */

typedef struct ScrewModifierData {
	ModifierData modifier;
	struct Object *ob_axis;
	int		steps;
	int		render_steps;
	int		iter;
	float	screw_ofs;
	float	angle;
	short	axis;
	short	flag;
} ScrewModifierData;

#define MOD_SCREW_NORMAL_FLIP	(1<<0)
#define MOD_SCREW_NORMAL_CALC	(1<<1)
#define MOD_SCREW_OBJECT_OFFSET	(1<<2)
// #define MOD_SCREW_OBJECT_ANGLE	(1<<4)

typedef struct OceanModifierData {
	ModifierData modifier;
	
	struct Ocean *ocean;
	struct OceanCache *oceancache;
	
	int		resolution;
	int		spatial_size;
	
	float	wind_velocity;
	
	float	damp;
	float	smallest_wave;
	float	depth;
	
	float	wave_alignment;
	float	wave_direction;
	float	wave_scale;
	
	float	chop_amount;
	float	foam_coverage;
	float	time;
	
	int		bakestart;
	int		bakeend;
	
	char	cachepath[240];	// FILE_MAX
	char	foamlayername[32];
	char	cached;
	char	geometry_mode;

	char	flag;
	char	refresh;

	short	repeat_x;
	short	repeat_y;

	int		seed;

	float	size;
	
	float	foam_fade;

	int pad;

} OceanModifierData;

#define MOD_OCEAN_GEOM_GENERATE	0
#define MOD_OCEAN_GEOM_DISPLACE	1
#define MOD_OCEAN_GEOM_SIM_ONLY	2

#define MOD_OCEAN_REFRESH_RESET			1
#define MOD_OCEAN_REFRESH_SIM			2
#define MOD_OCEAN_REFRESH_ADD			4
#define MOD_OCEAN_REFRESH_CLEAR_CACHE	8
#define MOD_OCEAN_REFRESH_TOPOLOGY		16

#define MOD_OCEAN_GENERATE_FOAM	1
#define MOD_OCEAN_GENERATE_NORMALS	2


typedef struct WarpModifierData {
	ModifierData modifier;

	/* keep in sync with MappingInfoModifierData */
	struct Tex *texture;
	struct Object *map_object;
	char uvlayer_name[32];
	int uvlayer_tmp;
	int texmapping;
	int pad10;
	/* end MappingInfoModifierData */

	float strength;

	struct Object *object_from;
	struct Object *object_to;
	struct CurveMapping *curfalloff;
	char defgrp_name[32];			/* optional vertexgroup name */
	float falloff_radius;
	char flag; /* not used yet */
	char falloff_type;
	char pad[2];
} WarpModifierData;

#define MOD_WARP_VOLUME_PRESERVE 1

typedef enum {
	eWarp_Falloff_None =		0,
	eWarp_Falloff_Curve =		1,
	eWarp_Falloff_Sharp =		2, /* PROP_SHARP */
	eWarp_Falloff_Smooth =		3, /* PROP_SMOOTH */
	eWarp_Falloff_Root =		4, /* PROP_ROOT */
	eWarp_Falloff_Linear =		5, /* PROP_LIN */
	eWarp_Falloff_Const =		6, /* PROP_CONST */
	eWarp_Falloff_Sphere =		7, /* PROP_SPHERE */
	/* PROP_RANDOM not used */
} WarpModifierFalloff;

typedef struct WeightVGEditModifierData {
	ModifierData modifier;

	/* Note: I tried to keep everything logically ordered - provided the
	 * alignment constraints... */

	char	defgrp_name[32];      /* Name of vertex group to edit. */

	short	edit_flags;     /* Using MOD_WVG_EDIT_* flags. */
	short	falloff_type;   /* Using MOD_WVG_MAPPING_* defines. */
	float	default_weight; /* Weight for vertices not in vgroup. */

	/* Mapping stuff. */
	struct CurveMapping *cmap_curve;  /* The custom mapping curve! */

	/* The add/remove vertices weight thresholds. */
	float	add_threshold, rem_threshold;

	/* Masking options. */
	float	mask_constant; /* The global "influence", if no vgroup nor tex is used as mask. */
	/* Name of mask vertex group from which to get weight factors. */
	char	mask_defgrp_name[32];

	/* Texture masking. */
	int		mask_tex_use_channel;      /* Which channel to use as weightf. */
	struct Tex *mask_texture;          /* The texture. */
	struct Object *mask_tex_map_obj;   /* Name of the map object. */
	/* How to map the texture (using MOD_DISP_MAP_* constants). */
	int		mask_tex_mapping;
	char	mask_tex_uvlayer_name[32]; /* Name of the UV map. */

	/* Padding... */
	int pad_i1;
} WeightVGEditModifierData;

/* WeightVGEdit flags. */
/* Use parametric mapping. */
//#define MOD_WVG_EDIT_MAP					(1 << 0)
/* Use curve mapping. */
//#define MOD_WVG_EDIT_CMAP					(1 << 1)
/* Reverse weights (in the [0.0, 1.0] standard range). */
//#define MOD_WVG_EDIT_REVERSE_WEIGHTS		(1 << 2)
/* Add vertices with higher weight than threshold to vgroup. */
#define MOD_WVG_EDIT_ADD2VG					(1 << 3)
/* Remove vertices with lower weight than threshold from vgroup. */
#define MOD_WVG_EDIT_REMFVG					(1 << 4)
/* Clamp weights. */
//#define MOD_WVG_EDIT_CLAMP					(1 << 5)

typedef struct WeightVGMixModifierData {
	ModifierData modifier;

	/* XXX Note: I tried to keep everything logically ordered – provided the
	 *           alignment constraints... */

	char	defgrp_name_a[32];      /* Name of vertex group to modify/weight. */
	char	defgrp_name_b[32];     /* Name of other vertex group to mix in. */
	float	default_weight_a;       /* Default weight value for first vgroup. */
	float	default_weight_b;      /* Default weight value to mix in. */
	char	mix_mode;             /* How second vgroups weights affect first ones */
	char	mix_set;              /* What vertices to affect. */

	char	pad_c1[6];

	/* Masking options. */
	float	mask_constant; /* The global "influence", if no vgroup nor tex is used as mask. */
	/* Name of mask vertex group from which to get weight factors. */
	char	mask_defgrp_name[32];

	/* Texture masking. */
	int		mask_tex_use_channel;      /* Which channel to use as weightf. */
	struct Tex *mask_texture;          /* The texture. */
	struct Object *mask_tex_map_obj;   /* Name of the map object. */
	int		mask_tex_mapping;          /* How to map the texture! */
	char	mask_tex_uvlayer_name[32]; /* Name of the UV map. */

	/* Padding... */
	int pad_i1;
} WeightVGMixModifierData;

/* How second vgroup's weights affect first ones. */
#define MOD_WVG_MIX_SET			1 /* Second weights replace weights. */
#define MOD_WVG_MIX_ADD			2 /* Second weights are added to weights. */
#define MOD_WVG_MIX_SUB			3 /* Second weights are subtracted from weights. */
#define MOD_WVG_MIX_MUL			4 /* Second weights are multiplied with weights. */
#define MOD_WVG_MIX_DIV			5 /* Second weights divide weights. */
#define MOD_WVG_MIX_DIF			6 /* Difference between second weights and weights. */
#define MOD_WVG_MIX_AVG			7 /* Average of both weights. */

/* What vertices to affect. */
#define MOD_WVG_SET_ALL			1 /* Affect all vertices. */
#define MOD_WVG_SET_A			2 /* Affect only vertices in first vgroup. */
#define MOD_WVG_SET_B			3 /* Affect only vertices in second vgroup. */
#define MOD_WVG_SET_OR			4 /* Affect only vertices in one vgroup or the other. */
#define MOD_WVG_SET_AND			5 /* Affect only vertices in both vgroups. */

typedef struct WeightVGProximityModifierData {
	ModifierData modifier;

	/* Note: I tried to keep everything logically ordered - provided the
	 * alignment constraints... */

	char	defgrp_name[32];      /* Name of vertex group to modify/weight. */

	/* Proximity modes. */
	int		proximity_mode;
	int		proximity_flags;

	/* Target object from which to calculate vertices distances. */
	struct Object *proximity_ob_target;

	/* Masking options. */
	float	mask_constant; /* The global "influence", if no vgroup nor tex is used as mask. */
	/* Name of mask vertex group from which to get weight factors. */
	char	mask_defgrp_name[32];

	/* Texture masking. */
	int		mask_tex_use_channel;      /* Which channel to use as weightf. */
	struct Tex *mask_texture;          /* The texture. */
	struct Object *mask_tex_map_obj;   /* Name of the map object. */
	int		mask_tex_mapping;          /* How to map the texture! */
	char	mask_tex_uvlayer_name[32]; /* Name of the UV Map. */

	float	min_dist, max_dist;        /* Distances mapping to 0.0/1.0 weights. */

	/* Put here to avoid breaking existing struct... */
	short	falloff_type;              /* Using MOD_WVG_MAPPING_* defines. */

	/* Padding... */
	short pad_s1;
} WeightVGProximityModifierData;

/* Modes of proximity weighting. */
/* Dist from target object to affected object. */
#define MOD_WVG_PROXIMITY_OBJECT			1 /* source vertex to other location */
/* Dist from target object to vertex. */
#define MOD_WVG_PROXIMITY_GEOMETRY			2 /* source vertex to other geometry */

/* Flags options for proximity weighting. */
/* Use nearest vertices of target obj, in MOD_WVG_PROXIMITY_GEOMETRY mode. */
#define MOD_WVG_PROXIMITY_GEOM_VERTS		(1 << 0)
/* Use nearest edges of target obj, in MOD_WVG_PROXIMITY_GEOMETRY mode. */
#define MOD_WVG_PROXIMITY_GEOM_EDGES		(1 << 1)
/* Use nearest faces of target obj, in MOD_WVG_PROXIMITY_GEOMETRY mode. */
#define MOD_WVG_PROXIMITY_GEOM_FACES		(1 << 2)

/* Defines common to all WeightVG modifiers. */
/* Mapping modes. */
#define MOD_WVG_MAPPING_NONE				0
#define MOD_WVG_MAPPING_CURVE				1
#define MOD_WVG_MAPPING_SHARP				2 /* PROP_SHARP */
#define MOD_WVG_MAPPING_SMOOTH				3 /* PROP_SMOOTH */
#define MOD_WVG_MAPPING_ROOT				4 /* PROP_ROOT */
/* PROP_LIN not used (same as NONE, here...). */
/* PROP_CONST not used. */
#define MOD_WVG_MAPPING_SPHERE				7 /* PROP_SPHERE */
#define MOD_WVG_MAPPING_RANDOM				8 /* PROP_RANDOM */
#define MOD_WVG_MAPPING_STEP				9 /* Median Step. */

/* Tex channel to be used as mask. */
#define MOD_WVG_MASK_TEX_USE_INT			1
#define MOD_WVG_MASK_TEX_USE_RED			2
#define MOD_WVG_MASK_TEX_USE_GREEN			3
#define MOD_WVG_MASK_TEX_USE_BLUE			4
#define MOD_WVG_MASK_TEX_USE_HUE			5
#define MOD_WVG_MASK_TEX_USE_SAT			6
#define MOD_WVG_MASK_TEX_USE_VAL			7
#define MOD_WVG_MASK_TEX_USE_ALPHA			8

/* Dynamic paint modifier flags */
#define MOD_DYNAMICPAINT_TYPE_CANVAS (1 << 0)
#define MOD_DYNAMICPAINT_TYPE_BRUSH (1 << 1)

typedef struct DynamicPaintModifierData {
	ModifierData modifier;

	struct DynamicPaintCanvasSettings *canvas;
	struct DynamicPaintBrushSettings *brush;
	int type;  /* ui display: canvas / brush */
	int pad;
} DynamicPaintModifierData;

#endif
