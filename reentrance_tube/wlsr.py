""" WLSR placement with TPB as parent volume and TTX as daughter """
from __future__ import annotations

import numpy as np
import pyg4ometry.geant4 as g4

from .profiles import (
    WLSR_TPB_THICKNESS,
    WLSR_TTX_THICKNESS,
    PROTECTION_GAP,
    make_outer_wlsr_atmospheric_profiles,
    make_inner_wlsr_argon_profiles,
)

def _add_wls_surfaces(materials, reg, tpb_pv, tetratex_pv, lar_pv, prefix=""):
    """Add optical border surfaces for WLS layers."""
    g4.BorderSurface(
        f"bsurface_{prefix}tpb_ttx", tpb_pv, tetratex_pv,
        materials.surfaces.to_tetratex, reg
    )
    g4.BorderSurface(
        f"bsurface_{prefix}lar_tpb", lar_pv, tpb_pv,
        materials.surfaces.lar_to_tpb, reg
    )
    g4.BorderSurface(
        f"bsurface_{prefix}tpb_lar", tpb_pv, lar_pv,
        materials.surfaces.lar_to_tpb, reg
    )

def place_inner_wlsr_in_argon(
    materials,
    registry,
    lar_cavity_lv,
    lar_cavity_pv,
    neckradius,
    tubeheight,
    totalheight,
    curvefraction,
    wls_height,
    inner_z,
    inner_r,
    outer_z,
    outer_r,
):
    """
    Place inner WLS layers in the underground argon cavity.
    NEW: TPB is parent, TTX is daughter inside with 1 μm margin.
    """
    result = make_inner_wlsr_argon_profiles(
        neckradius, tubeheight, totalheight, curvefraction, wls_height, inner_z, inner_r, outer_z, outer_r
    )
    
    # Unpack profiles
    (tpb_outer_z, tpb_outer_r, tpb_inner_z, tpb_inner_r,
     ttx_outer_z, ttx_outer_r, ttx_inner_z, ttx_inner_r) = result

    if len(tpb_outer_z) > 1:
        # Create TPB polycones (PARENT/MOTHER volume)
        tpb_outer_bound = g4.solid.GenericPolycone(
            "tpb_inner_argon_outer_bound", 0, 2 * np.pi, tpb_outer_r, tpb_outer_z, registry, "mm"
        )
        tpb_inner_bound = g4.solid.GenericPolycone(
            "tpb_inner_argon_inner_bound", 0, 2 * np.pi, tpb_inner_r, tpb_inner_z, registry, "mm"
        )
        tpb_solid = g4.solid.Subtraction(
            "tpb_inner_argon_solid", tpb_outer_bound, tpb_inner_bound,
            [[0, 0, 0], [0, 0, 0, "mm"]], registry
        )
        
        wls_tpb_inner_lv = g4.LogicalVolume(
            tpb_solid, materials.tpb_on_tetratex, "wls_tpb_inner_argon_lv", registry
        )
        wls_tpb_inner_lv.pygeom_color_rgba = [0.0, 0.5, 1.0, 0.7]
        tpb_inner_pv = g4.PhysicalVolume(
            [0, 0, 0], [0, 0, 0, "mm"], wls_tpb_inner_lv,
            "wls_tpb_inner_argon", lar_cavity_lv, registry=registry
        )

        if len(ttx_outer_z) > 1:
            # Create TTX polycones (DAUGHTER volume inside TPB)
            tetratex_outer_bound = g4.solid.GenericPolycone(
                "tetratex_inner_argon_outer_bound", 0, 2 * np.pi, ttx_outer_r, ttx_outer_z, registry, "mm"
            )
            tetratex_inner_bound = g4.solid.GenericPolycone(
                "tetratex_inner_argon_inner_bound", 0, 2 * np.pi, ttx_inner_r, ttx_inner_z, registry, "mm"
            )
            tetratex_solid = g4.solid.Subtraction(
                "wls_tetratex_inner_argon_solid", tetratex_outer_bound, tetratex_inner_bound,
                [[0, 0, 0], [0, 0, 0, "mm"]], registry
            )
            
            wls_tetratex_inner_lv = g4.LogicalVolume(
                tetratex_solid, materials.tetratex, "wls_tetratex_inner_argon_lv", registry
            )
            wls_tetratex_inner_lv.pygeom_color_rgba = [1.0, 0.5, 0.0, 0.7]
            
            # Place TTX INSIDE TPB (parent is TPB logical volume)
            tetratex_inner_pv = g4.PhysicalVolume(
                [0, 0, 0], [0, 0, 0, "mm"], wls_tetratex_inner_lv,
                "wls_tetratex_inner_argon", wls_tpb_inner_lv, registry=registry
            )
            
            _add_wls_surfaces(materials, registry, tpb_inner_pv, tetratex_inner_pv, lar_cavity_pv, prefix="inner_")

    print(f"\nInner WLS layer (Steel -> 5nm gap -> TPB 1μm coating (parent) -> TTX core (daughter) -> UAr):")
    print(f"  Height: {wls_height}mm from bottom")
    print(f"  Gap to tube: {PROTECTION_GAP*1e6:.1f} nm")
    print(f"  TPB uniformly coats TTX with {WLSR_TPB_THICKNESS*1e3:.1f} μm on all sides")

def place_outer_wlsr_in_atmospheric(
    materials,
    registry,
    lar_mother_lv,
    lar_mother_pv,
    neckradius,
    tubeheight,
    totalheight,
    curvefraction,
    wls_height,
    outer_z,
    outer_r,
):
    """
    Place outer WLS layers in the atmospheric argon.
    NEW: TPB is parent, TTX is daughter inside with 1 μm margin.
    """
    (tpb_outer_z, tpb_outer_r, tpb_inner_z, tpb_inner_r,
     ttx_outer_z, ttx_outer_r, ttx_inner_z, ttx_inner_r) = make_outer_wlsr_atmospheric_profiles(
        neckradius, tubeheight, totalheight, curvefraction, wls_height, outer_z, outer_r
    )

    if len(tpb_outer_z) > 1:
        # Create TPB polycones (PARENT/MOTHER volume)
        tpb_outer_bound = g4.solid.GenericPolycone(
            "tpb_outer_atmospheric_outer_bound", 0, 2 * np.pi, tpb_outer_r, tpb_outer_z, registry, "mm"
        )
        tpb_inner_bound = g4.solid.GenericPolycone(
            "tpb_outer_atmospheric_inner_bound", 0, 2 * np.pi, tpb_inner_r, tpb_inner_z, registry, "mm"
        )
        tpb_solid = g4.solid.Subtraction(
            "tpb_outer_atmospheric_solid", tpb_outer_bound, tpb_inner_bound,
            [[0, 0, 0], [0, 0, 0, "mm"]], registry
        )
        
        wls_tpb_outer_lv = g4.LogicalVolume(
            tpb_solid, materials.tpb_on_tetratex, "wls_tpb_outer_atmospheric_lv", registry
        )
        wls_tpb_outer_lv.pygeom_color_rgba = [0.0, 1.0, 0.0, 0.7]
        tpb_outer_pv = g4.PhysicalVolume(
            [0, 0, 0], [0, 0, 0, "mm"], wls_tpb_outer_lv,
            "wls_tpb_outer_atmospheric", lar_mother_lv, registry=registry
        )

        if len(ttx_outer_z) > 1:
            # Create TTX polycones (DAUGHTER volume inside TPB)
            tetratex_outer_bound = g4.solid.GenericPolycone(
                "tetratex_outer_atmospheric_outer_bound", 0, 2 * np.pi, ttx_outer_r, ttx_outer_z, registry, "mm"
            )
            tetratex_inner_bound = g4.solid.GenericPolycone(
                "tetratex_outer_atmospheric_inner_bound", 0, 2 * np.pi, ttx_inner_r, ttx_inner_z, registry, "mm"
            )
            tetratex_solid = g4.solid.Subtraction(
                "wls_tetratex_outer_atmospheric_solid", tetratex_outer_bound, tetratex_inner_bound,
                [[0, 0, 0], [0, 0, 0, "mm"]], registry
            )
            
            wls_tetratex_outer_lv = g4.LogicalVolume(
                tetratex_solid, materials.tetratex, "wls_tetratex_outer_atmospheric_lv", registry
            )
            wls_tetratex_outer_lv.pygeom_color_rgba = [1.0, 0.0, 0.0, 0.7]
            
            # Place TTX INSIDE TPB (parent is TPB logical volume)
            tetratex_outer_pv = g4.PhysicalVolume(
                [0, 0, 0], [0, 0, 0, "mm"], wls_tetratex_outer_lv,
                "wls_tetratex_outer_atmospheric", wls_tpb_outer_lv, registry=registry
            )
            
            _add_wls_surfaces(materials, registry, tpb_outer_pv, tetratex_outer_pv, lar_mother_pv, prefix="outer_")

    print(f"\nOuter WLS layer (Tube -> 5nm gap -> TPB 1μm coating (parent) -> TTX core (daughter) -> Atmospheric Ar):")
    print(f"  Height: {wls_height}mm from bottom")
    print(f"  Gap to tube: {PROTECTION_GAP*1e6:.1f} nm")
    print(f"  TPB uniformly coats TTX with {WLSR_TPB_THICKNESS*1e3:.1f} μm on all sides")
