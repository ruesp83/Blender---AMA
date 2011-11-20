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
 * The Original Code is Copyright (C) 2011 Blender Foundation.
 * All rights reserved.
 *
 *
 * Contributor(s): Blender Foundation,
 *                 Sergey Sharybin
 *
 * ***** END GPL LICENSE BLOCK *****
 */

/** \file blender/editors/space_clip/clip_header.c
 *  \ingroup spclip
 */

#include <string.h>

#include "DNA_windowmanager_types.h"

#include "MEM_guardedalloc.h"

#include "BLI_blenlib.h"
#include "BLI_utildefines.h"

#include "BKE_context.h"
#include "BKE_screen.h"

#include "ED_screen.h"
#include "ED_util.h"

#include "WM_types.h"
#include "WM_api.h"

#include "UI_interface.h"
#include "UI_resources.h"

#include "clip_intern.h" /* own include */

/* ************************ header area region *********************** */

/************************** properties ******************************/

static ARegion *clip_has_properties_region(ScrArea *sa)
{
	ARegion *ar, *arnew;

	ar= BKE_area_find_region_type(sa, RGN_TYPE_UI);
	if(ar)
		return ar;

	/* add subdiv level; after header */
	ar= BKE_area_find_region_type(sa, RGN_TYPE_HEADER);

	/* is error! */
	if(ar==NULL)
		return NULL;

	arnew= MEM_callocN(sizeof(ARegion), "clip properties region");

	BLI_insertlinkafter(&sa->regionbase, ar, arnew);
	arnew->regiontype= RGN_TYPE_UI;
	arnew->alignment= RGN_ALIGN_RIGHT;

	arnew->flag= RGN_FLAG_HIDDEN;

	return arnew;
}

static int properties_poll(bContext *C)
{
	return (CTX_wm_space_clip(C) != NULL);
}

static int properties_exec(bContext *C, wmOperator *UNUSED(op))
{
	ScrArea *sa= CTX_wm_area(C);
	ARegion *ar= clip_has_properties_region(sa);

	if(ar)
		ED_region_toggle_hidden(C, ar);

	return OPERATOR_FINISHED;
}

void CLIP_OT_properties(wmOperatorType *ot)
{
	/* identifiers */
	ot->name= "Properties";
	ot->description= "Toggle clip properties panel";
	ot->idname= "CLIP_OT_properties";

	/* api callbacks */
	ot->exec= properties_exec;
	ot->poll= properties_poll;
}

/************************** tools ******************************/

static ARegion *clip_has_tools_region(ScrArea *sa)
{
	ARegion *ar, *artool=NULL, *arprops=NULL, *arhead;

	for(ar= sa->regionbase.first; ar; ar= ar->next) {
		if(ar->regiontype==RGN_TYPE_TOOLS)
			artool= ar;
		if(ar->regiontype==RGN_TYPE_TOOL_PROPS)
			arprops= ar;
	}

	/* tool region hide/unhide also hides props */
	if(arprops && artool)
		return artool;

	if(artool==NULL) {
		/* add subdiv level; after header */
		arhead= BKE_area_find_region_type(sa, RGN_TYPE_HEADER);

		/* is error! */
		if(arhead==NULL)
			return NULL;

		artool= MEM_callocN(sizeof(ARegion), "clip tools region");

		BLI_insertlinkafter(&sa->regionbase, arhead, artool);
		artool->regiontype= RGN_TYPE_TOOLS;
		artool->alignment= RGN_ALIGN_LEFT;

		artool->flag= RGN_FLAG_HIDDEN;
	}

	if(arprops==NULL) {
		/* add extra subdivided region for tool properties */
		arprops= MEM_callocN(sizeof(ARegion), "tool props for clip");

		BLI_insertlinkafter(&sa->regionbase, artool, arprops);
		arprops->regiontype= RGN_TYPE_TOOL_PROPS;
		arprops->alignment= RGN_ALIGN_BOTTOM|RGN_SPLIT_PREV;
	}

	return artool;
}

static int tools_poll(bContext *C)
{
	return (CTX_wm_space_clip(C) != NULL);
}

static int tools_exec(bContext *C, wmOperator *UNUSED(op))
{
	ScrArea *sa= CTX_wm_area(C);
	ARegion *ar= clip_has_tools_region(sa);

	if(ar)
		ED_region_toggle_hidden(C, ar);

	return OPERATOR_FINISHED;
}

void CLIP_OT_tools(wmOperatorType *ot)
{
	/* identifiers */
	ot->name= "Tools";
	ot->description= "Toggle clip tools panel";
	ot->idname= "CLIP_OT_tools";

	/* api callbacks */
	ot->exec= tools_exec;
	ot->poll= tools_poll;
}

/************************** redo panel ******************************/

static void clip_panel_operator_redo_buts(const bContext *C, Panel *pa, wmOperator *op)
{
	uiLayoutOperatorButs(C, pa->layout, op, NULL, 'V', 0);
}

static void clip_panel_operator_redo_header(const bContext *C, Panel *pa)
{
	wmOperator *op= WM_operator_last_redo(C);

	if(op) BLI_strncpy(pa->drawname, op->type->name, sizeof(pa->drawname));
	else BLI_strncpy(pa->drawname, "Operator", sizeof(pa->drawname));
}

static void clip_panel_operator_redo_operator(const bContext *C, Panel *pa, wmOperator *op)
{
	if(op->type->flag & OPTYPE_MACRO) {
		for(op= op->macro.first; op; op= op->next) {
			uiItemL(pa->layout, op->type->name, ICON_NONE);
			clip_panel_operator_redo_operator(C, pa, op);
		}
	}
	else {
		clip_panel_operator_redo_buts(C, pa, op);
	}
}

static void clip_panel_operator_redo(const bContext *C, Panel *pa)
{
	wmOperator *op= WM_operator_last_redo(C);
	uiBlock *block;

	if(op==NULL)
		return;
	if(WM_operator_poll((bContext*)C, op->type) == 0)
		return;

	block= uiLayoutGetBlock(pa->layout);

	if(ED_undo_valid(C, op->type->name)==0)
		uiLayoutSetEnabled(pa->layout, 0);

	/* note, blockfunc is a default but->func, use Handle func to allow button callbacks too */
	uiBlockSetHandleFunc(block, ED_undo_operator_repeat_cb_evt, op);

	clip_panel_operator_redo_operator(C, pa, op);
}

void ED_clip_tool_props_register(ARegionType *art)
{
	PanelType *pt;

	pt= MEM_callocN(sizeof(PanelType), "spacetype clip panel last operator");
	strcpy(pt->idname, "CLIP_PT_last_operator");
	strcpy(pt->label, "Operator");
	pt->draw_header= clip_panel_operator_redo_header;
	pt->draw= clip_panel_operator_redo;
	BLI_addtail(&art->paneltypes, pt);
}