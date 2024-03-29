/*
    Global definitions.

    Copyright (C) 1999 Mads Bak

    This file is part of the SIMPSON General NMR Simulation Package

    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program; if not, write to the Free Software
    Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
    
*/

#ifndef __DEFS_H
#define __DEFS_H

#ifdef __cplusplus
extern "C" {
#endif

#define VERBOSE_NOT       0
#define VERBOSE_SPINSYS   (1 << 0)
#define VERBOSE_PROGRESS  (1 << 1)
#define VERBOSE_PAR       (1 << 2)
#define VERBOSE_SIMINFO   (1 << 3)
#define VERBOSE_OPER      (1 << 4)
#define VERBOSE_POWDER    (1 << 5)
#define VERBOSE_RFPROF    (1 << 6) 

extern int verbose;

#define VARIOUS_NOT         0
#define VARIOUS_NODIAGPROP  (1 << 0)
#define VARIOUS_NOCHECKPROPTIME  (2 << 0)

extern int various;

#define MAGIC_ANGLE  54.735610317245345868286676704883575
#define SQRT3BY8      0.612372435695794470333908066095318
#define SQRT1BY3      0.577350269189625731058868041145615
#define SQRT2BY3      0.816496580927726034460079063137527
#define SQRT2_2BY3    0.942809041582063356301546264148783
#define DEG2RAD       0.017453292519943295474371680597869
#define RAD2DEG      57.295779513082322864647721871733665

#ifdef __cplusplus
}
#endif

#endif /* __DEFS_H */   
