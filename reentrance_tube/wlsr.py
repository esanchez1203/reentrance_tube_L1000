"""WLS reflector construction for reentrance tube."""
from __future__ import annotations

import numpy as np
import pyg4ometry.geant4 as g4

from .profiles import (
    WLSR_TPB_THICKNESS,
    WLSR_TTX_THICKNESS,
    make_tetratex_inner_profiles,
    make_tetratex_outer_profiles,
    make_wls_inner_profiles,
    make_wls_outer_profiles,
)


def _add_wls_surfaces(materials, reg, tpb_pv, tetratex_pv, lar_pv, prefix=""):
    """
    Add optical border surfaces for WLS layers.
    
    Parameters
    ----------
    materials : OpticalMaterialRegistry
        Material registry with surfaces
    reg : g4.Registry
        pyg4ometry registry
    tpb_pv : g4.PhysicalVolume
        TPB physical volume
    tetratex_pv : g4.PhysicalVolume
        Tetratex physical volume
    lar_pv : g4.PhysicalVolume
        Liquid argon physical volume
    prefix : str
        Prefix for surface names (e.g., "outer_" or "inner_")
    """
    # Surface between TPB and Tetratex
    g4.BorderSurface(
        f"bsurface_{prefix}tpb_ttx",
        tpb_pv,
        tetratex_pv,
        materials.surfaces.to_tetratex,
        reg,
    )

    # Surface between LAr and TPB (bidirectional)
    g4.BorderSurface(
        f"bsurface_{prefix}lar_tpb",
        lar_pv,
        tpb_pv,
        materials.surfaces.lar_to_tpb,
        reg,
    )
    g4.BorderSurface(
        f"bsurface_{prefix}tpb_lar",
        tpb_pv,
        lar_pv,
        materials.surfaces.lar_to_tpb,
        reg,
    )


def place_outer_wlsr(
    materials,
    registry,
    steel_tube_lv,
    neckradius,
    tubeheight,
    totalheight,
    curvefraction,
    wls_height,
    outer_z,
    outer_r,
    lar_mother_pv,
):
    """
    Place outer WLS layers on the exterior steel surface.
    
    Structure: Steel -> Tetratex (254um) -> TPB (1um) -> LAr
    
    Parameters
    ----------
    materials : OpticalMaterialRegistry
        Material registry
    registry : g4.Registry
        pyg4ometry registry
    steel_tube_lv : g4.LogicalVolume
        Steel tube logical volume (mother for WLS)
    neckradius : float
        Tube neck radius in mm
    tubeheight : float
        Tube height in mm
    totalheight : float
        Total height in mm
    curvefraction : float
        Fraction of curved section
    wls_height : float
        Height of WLS from bottom in mm
    outer_z : list
        Steel outer z coordinates
    outer_r : list
        Steel outer r coordinates
    lar_mother_pv : g4.PhysicalVolume
        LAr physical volume for surface
    """
    # Get WLS outer profile (for TPB outer surface)
    wls_outer_z, wls_outer_r, _, _ = make_wls_outer_profiles(
        neckradius, tubeheight, totalheight, curvefraction, wls_height, outer_z, outer_r
    )

    # Get Tetratex profiles
    ttx_outer_z, ttx_outer_r, ttx_inner_z, ttx_inner_r = make_tetratex_outer_profiles(
        neckradius, tubeheight, totalheight, curvefraction, wls_height, outer_z, outer_r
    )

    if len(ttx_outer_z) > 1:
        # Create Tetratex solid (inner layer)
        tetratex_outer_bound = g4.solid.GenericPolycone(
            "tetratex_outer_bound", 0, 2 * np.pi, ttx_outer_r, ttx_outer_z, registry, "mm"
        )
        tetratex_inner_bound = g4.solid.GenericPolycone(
            "tetratex_inner_bound", 0, 2 * np.pi, ttx_inner_r, ttx_inner_z, registry, "mm"
        )

        tetratex_solid = g4.solid.Subtraction(
            "wls_tetratex_solid",
            tetratex_outer_bound,
            tetratex_inner_bound,
            [[0, 0, 0], [0, 0, 0, "mm"]],
            registry,
        )

        # Create Tetratex logical volume
        wls_tetratex_lv = g4.LogicalVolume(
            tetratex_solid, materials.tetratex, "wls_tetratex_lv", registry
        )
        wls_tetratex_lv.pygeom_color_rgba = [1.0, 0.0, 0.0, 0.7]  # Red, semi-transparent

        # Place Tetratex in steel tube
        tetratex_pv = g4.PhysicalVolume(
            [0, 0, 0], [0, 0, 0, "mm"], wls_tetratex_lv, "wls_tetratex", steel_tube_lv, registry=registry
        )

    if len(wls_outer_z) > 1:
        # Create TPB solid (outer layer)
        wls_outer_bound = g4.solid.GenericPolycone(
            "wls_outer_bound", 0, 2 * np.pi, wls_outer_r, wls_outer_z, registry, "mm"
        )
        wls_inner_bound = g4.solid.GenericPolycone(
            "wls_inner_bound", 0, 2 * np.pi, ttx_outer_r, ttx_outer_z, registry, "mm"
        )

        wls_solid = g4.solid.Subtraction(
            "wls_solid", wls_outer_bound, wls_inner_bound, [[0, 0, 0], [0, 0, 0, "mm"]], registry
        )

        # Create TPB logical volume
        wls_outer_lv = g4.LogicalVolume(wls_solid, materials.tpb_on_tetratex, "wls_outer_lv", registry)
        wls_outer_lv.pygeom_color_rgba = [0.0, 1.0, 0.0, 1.0]  # Green, opaque

        # Place TPB in steel tube (sibling to Tetratex)
        wls_outer_pv = g4.PhysicalVolume(
            [0, 0, 0], [0, 0, 0, "mm"], wls_outer_lv, "tpb_coating", steel_tube_lv, registry=registry
        )

        # Add optical surfaces
        _add_wls_surfaces(materials, registry, wls_outer_pv, tetratex_pv, lar_mother_pv, prefix="outer_")

    print(f"\nOuter WLS layer (two-layer structure):")
    print(f"  Height: {wls_height}mm from bottom")
    print(f"  Inner layer (Tetratex): {WLSR_TTX_THICKNESS:.3f}mm thick, on steel surface")
    print(f"  Outer layer (TPB): {WLSR_TPB_THICKNESS:.3f}mm thick, in contact with Ar")


def place_inner_wlsr(
    materials,
    registry,
    steel_tube_lv,
    neckradius,
    tubeheight,
    totalheight,
    curvefraction,
    wls_height,
    inner_z,
    inner_r,
    outer_z,
    outer_r,
    lar_cavity_pv,
):
    """
    Place inner WLS layers on the interior steel surface.
    
    Structure: Steel <- TPB (1um) <- Tetratex (254um) <- LAr
    (layers grow INTO the steel)
    
    Parameters
    ----------
    materials : OpticalMaterialRegistry
        Material registry
    registry : g4.Registry
        pyg4ometry registry
    steel_tube_lv : g4.LogicalVolume
        Steel tube logical volume (mother for WLS)
    neckradius : float
        Tube neck radius in mm
    tubeheight : float
        Tube height in mm
    totalheight : float
        Total height in mm
    curvefraction : float
        Fraction of curved section
    wls_height : float
        Height of WLS from bottom in mm
    inner_z : list
        Steel inner z coordinates
    inner_r : list
        Steel inner r coordinates
    outer_z : list
        Steel outer z coordinates
    outer_r : list
        Steel outer r coordinates
    lar_cavity_pv : g4.PhysicalVolume
        LAr cavity physical volume for surface
    """
    # Get WLS profiles (inner layer, on steel inner surface)
    wls_outer_z, wls_outer_r, wls_inner_z, wls_inner_r = make_wls_inner_profiles(
        neckradius, tubeheight, totalheight, curvefraction, wls_height, inner_z, inner_r, outer_z, outer_r
    )

    # Get Tetratex profiles (outer layer, coating on WLS)
    ttx_outer_z, ttx_outer_r, ttx_inner_z, ttx_inner_r = make_tetratex_inner_profiles(
        neckradius, tubeheight, totalheight, curvefraction, wls_height, inner_z, inner_r, outer_z, outer_r
    )

    if len(wls_outer_z) > 1:
        # Create WLS/TPB solid (inner layer, on steel inner surface)
        wls_outer_bound = g4.solid.GenericPolycone(
            "wls_inner_outer_bound", 0, 2 * np.pi, wls_outer_r, wls_outer_z, registry, "mm"
        )
        wls_inner_bound = g4.solid.GenericPolycone(
            "wls_inner_inner_bound", 0, 2 * np.pi, wls_inner_r, wls_inner_z, registry, "mm"
        )

        wls_solid = g4.solid.Subtraction(
            "wls_inner_solid",
            wls_outer_bound,
            wls_inner_bound,
            [[0, 0, 0], [0, 0, 0, "mm"]],
            registry,
        )

        wls_inner_lv = g4.LogicalVolume(wls_solid, materials.tpb_on_tetratex, "wls_inner_lv", registry)
        wls_inner_lv.pygeom_color_rgba = [0.0, 0.5, 1.0, 1.0]  # Light blue

        # Place in steel tube
        wls_inner_pv = g4.PhysicalVolume(
            [0, 0, 0], [0, 0, 0, "mm"], wls_inner_lv, "wls_inner_coating", steel_tube_lv, registry=registry
        )

    if len(ttx_outer_z) > 1:
        # Create Tetratex solid (outer layer, on WLS)
        tetratex_outer_bound = g4.solid.GenericPolycone(
            "tetratex_inner_outer_bound", 0, 2 * np.pi, ttx_outer_r, ttx_outer_z, registry, "mm"
        )
        tetratex_inner_bound = g4.solid.GenericPolycone(
            "tetratex_inner_inner_bound", 0, 2 * np.pi, ttx_inner_r, ttx_inner_z, registry, "mm"
        )

        tetratex_solid = g4.solid.Subtraction(
            "wls_tetratex_inner_solid",
            tetratex_outer_bound,
            tetratex_inner_bound,
            [[0, 0, 0], [0, 0, 0, "mm"]],
            registry,
        )

        wls_tetratex_inner_lv = g4.LogicalVolume(
            tetratex_solid, materials.tetratex, "wls_tetratex_inner_lv", registry
        )
        wls_tetratex_inner_lv.pygeom_color_rgba = [1.0, 0.5, 0.0, 0.7]  # Orange

        # Place in steel tube (sibling to WLS)
        tetratex_inner_pv = g4.PhysicalVolume(
            [0, 0, 0],
            [0, 0, 0, "mm"],
            wls_tetratex_inner_lv,
            "wls_tetratex_inner",
            steel_tube_lv,
            registry=registry,
        )

        # Add optical surfaces
        _add_wls_surfaces(materials, registry, wls_inner_pv, tetratex_inner_pv, lar_cavity_pv, prefix="inner_")

    print(f"\nInner WLS layer (two-layer structure):")
    print(f"  Height: {wls_height}mm from bottom")
    print(f"  Inner layer (WLS/TPB): {WLSR_TPB_THICKNESS:.3f}mm thick, on steel inner surface")
    print(f"  Outer layer (Tetratex): {WLSR_TTX_THICKNESS:.3f}mm thick, coating on WLS")
