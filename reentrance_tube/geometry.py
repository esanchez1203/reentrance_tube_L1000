"""Main reentrance tube geometry construction."""
from __future__ import annotations

import numpy as np
import pyg4ometry.geant4 as g4

from .profiles import (
    make_316l_ss_profiles,
    make_inner_profile,
    make_ofhc_cu_profiles,
    make_outer_profile,
)
from .wlsr import place_inner_wlsr, place_outer_wlsr


def construct_reentrance_tube(
    materials,
    registry,
    mother_lv,
    mother_pv,
    neckradius=1015,
    tubeheight=6247,
    totalheight=8000,
    curvefraction=0.05,
    with_outer_wls=False,
    with_inner_wls=False,
    wls_height=2179,
    with_ofhc_cu=False,
    with_316l_ss=False,
    ofhc_start_height=2179,
    ofhc_end_height=4184,
    ss_start_height=4184,
):
    
    print(f"\nReentrance tube: radius={neckradius}mm, height={tubeheight}mm")

    # Generate steel profiles
    outer_z, outer_r = make_outer_profile(neckradius, tubeheight, totalheight, curvefraction)
    inner_z, inner_r = make_inner_profile(neckradius, tubeheight, totalheight, curvefraction)

    print(f"Generated {len(outer_z)} profile points")

    # Steel tube
    tube_solid = g4.solid.GenericPolycone(
        "reentrance_tube_solid", 0, 2 * np.pi, outer_r, outer_z, registry, "mm"
    )
    tube_lv = g4.LogicalVolume(tube_solid, materials.metal_steel, "reentrance_tube_lv", registry)
    tube_lv.pygeom_color_rgba = [0.5, 0.5, 0.5, 0.8]
    tube_pv = g4.PhysicalVolume(
        [0, 0, 0], [0, 0, 0, "mm"], tube_lv, "reentrance_tube", mother_lv, registry=registry
    )

    # Underground argon cavity
    cavity_solid = g4.solid.GenericPolycone(
        "underground_argon_solid", 0, 2 * np.pi, inner_r, inner_z, registry, "mm"
    )
    cavity_lv = g4.LogicalVolume(cavity_solid, materials.liquidargon, "underground_argon_lv", registry)
    cavity_lv.pygeom_color_rgba = [0.1, 0.8, 0.3, 0.3]
    cavity_pv = g4.PhysicalVolume(
        [0, 0, 0], [0, 0, 0, "mm"], cavity_lv, "underground_argon", tube_lv, registry=registry
    )

    # Add outer WLS if requested
    if with_outer_wls:
        place_outer_wlsr(
            materials,
            registry,
            tube_lv,
            neckradius,
            tubeheight,
            totalheight,
            curvefraction,
            wls_height,
            outer_z,
            outer_r,
            mother_pv,
        )

    # Add inner WLS if requested
    if with_inner_wls:
        place_inner_wlsr(
            materials,
            registry,
            tube_lv,
            neckradius,
            tubeheight,
            totalheight,
            curvefraction,
            wls_height,
            inner_z,
            inner_r,
            outer_z,
            outer_r,
            cavity_pv,
        )

    # Add OFHC copper if requested
    if with_ofhc_cu:
        ofhc_outer_z, ofhc_outer_r, ofhc_inner_z, ofhc_inner_r = make_ofhc_cu_profiles(
            neckradius,
            tubeheight,
            totalheight,
            curvefraction,
            ofhc_start_height,
            ofhc_end_height,
            outer_z,
            outer_r,
            inner_z,
            inner_r,
        )

        if len(ofhc_outer_z) > 1:
            # Create outer boundary polycone
            ofhc_outer_bound = g4.solid.GenericPolycone(
                "ofhc_cu_outer_bound", 0, 2 * np.pi, ofhc_outer_r, ofhc_outer_z, registry, "mm"
            )
            # Create inner boundary polycone
            ofhc_inner_bound = g4.solid.GenericPolycone(
                "ofhc_cu_inner_bound", 0, 2 * np.pi, ofhc_inner_r, ofhc_inner_z, registry, "mm"
            )

            # Subtract to create solid copper between surfaces
            ofhc_cu_solid = g4.solid.Subtraction(
                "ofhc_cu_solid", ofhc_outer_bound, ofhc_inner_bound, [[0, 0, 0], [0, 0, 0, "mm"]], registry
            )

            ofhc_cu_lv = g4.LogicalVolume(ofhc_cu_solid, materials.metal_copper, "ofhc_cu_lv", registry)
            ofhc_cu_lv.pygeom_color_rgba = [1.0, 0.5, 0.0, 1.0]

            g4.PhysicalVolume([0, 0, 0], [0, 0, 0, "mm"], ofhc_cu_lv, "ofhc_cu", tube_lv, registry=registry)

            print(f"\nOFHC Copper solid cylinder:")
            print(f"  Height: {ofhc_start_height}mm to {ofhc_end_height}mm from bottom")
            print(f"  Solid volume between steel inner and outer surfaces")
            print(f"  Sample thickness at mid-height:")
            mid_idx = len(ofhc_outer_r) // 2
            print(f"    Outer r: {ofhc_outer_r[mid_idx]:.1f}mm")
            print(f"    Inner r: {ofhc_inner_r[mid_idx]:.1f}mm")
            print(f"    Copper thickness: {ofhc_outer_r[mid_idx] - ofhc_inner_r[mid_idx]:.1f}mm")

    # Add 316L stainless steel if requested
    if with_316l_ss:
        ss_outer_z, ss_outer_r, ss_inner_z, ss_inner_r = make_316l_ss_profiles(
            neckradius, tubeheight, totalheight, curvefraction, ss_start_height, outer_z, outer_r, inner_z, inner_r
        )

        if len(ss_outer_z) > 1:
            # Create outer and inner boundaries
            ss_outer_bound = g4.solid.GenericPolycone(
                "ss_316l_outer_bound", 0, 2 * np.pi, ss_outer_r, ss_outer_z, registry, "mm"
            )
            ss_inner_bound = g4.solid.GenericPolycone(
                "ss_316l_inner_bound", 0, 2 * np.pi, ss_inner_r, ss_inner_z, registry, "mm"
            )

            # Create shell via subtraction
            ss_316l_solid = g4.solid.Subtraction(
                "ss_316l_solid", ss_outer_bound, ss_inner_bound, [[0, 0, 0], [0, 0, 0, "mm"]], registry
            )

            # Create logical volume
            ss_316l_lv = g4.LogicalVolume(ss_316l_solid, materials.metal_steel, "ss_316l_lv", registry)
            ss_316l_lv.pygeom_color_rgba = [0.7, 0.7, 0.8, 1.0]  # Light blue-gray for 316L

            # Place in steel tube
            g4.PhysicalVolume([0, 0, 0], [0, 0, 0, "mm"], ss_316l_lv, "ss_316l", tube_lv, registry=registry)

            print(f"\n316L Stainless Steel cylinder:")
            print(f"  Height: {ss_start_height}mm from bottom to top of cylindrical section")
            print(f"  Material: {materials.metal_steel.name}")
            print(f"  Sample thickness at mid-height:")
            mid_idx = len(ss_outer_r) // 2
            print(f"    Outer r: {ss_outer_r[mid_idx]:.1f}mm")
            print(f"    Inner r: {ss_inner_r[mid_idx]:.1f}mm")
            print(f"    Thickness: {ss_outer_r[mid_idx] - ss_inner_r[mid_idx]:.1f}mm")

    return tube_lv, cavity_lv, tube_pv, cavity_pv
