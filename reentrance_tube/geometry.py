"""Modified reentrance tube geometry with inner WLSR in underground argon."""
from __future__ import annotations

import numpy as np
import pyg4ometry.geant4 as g4

from .profiles import (
    make_316l_ss_profiles,
    make_inner_profile,
    make_ofhc_cu_profiles,
    make_outer_profile,
)
from .wlsr import place_inner_wlsr_in_argon, place_outer_wlsr_in_atmospheric


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
    """
    Construct reentrance tube with modified inner WLSR placement.
    
    MODIFICATION: Inner WLSR is now placed in underground argon instead of steel tube.
    Structure: LAr -> TPB -> TTX -> Gap -> Tube inner surface
    """
    
    print(f"\nReentrance tube: radius={neckradius}mm, height={tubeheight}mm")

    # Generate steel profiles - PASS wls_height to add critical boundary points
    outer_z, outer_r = make_outer_profile(neckradius, tubeheight, totalheight, curvefraction, wls_height)
    inner_z, inner_r = make_inner_profile(neckradius, tubeheight, totalheight, curvefraction, wls_height)

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
    cavity_lv.pygeom_color_rgba = [1.0, 0.0, 1.0, 0.3]
    cavity_pv = g4.PhysicalVolume(
        [0, 0, 0], [0, 0, 0, "mm"], cavity_lv, "underground_argon", tube_lv, registry=registry
    )

    # Inner WLSR in underground argon
    if with_inner_wls:
        place_inner_wlsr_in_argon(
            materials,
            registry,
            cavity_lv,      
            cavity_pv,
            neckradius,
            tubeheight,
            totalheight,
            curvefraction,
            wls_height,
            inner_z,
            inner_r,
            outer_z,
            outer_r,
        )

    # Outer WLSR 
    if with_outer_wls:
        place_outer_wlsr_in_atmospheric(
            materials,
            registry,
            mother_lv,  
            mother_pv,
            neckradius,
            tubeheight,
            totalheight,
            curvefraction,
            wls_height,
            outer_z,
            outer_r,
        )

    # OFHC copper (unchanged)
    if with_ofhc_cu:
        ofhc_outer_z, ofhc_outer_r, ofhc_inner_z, ofhc_inner_r = make_ofhc_cu_profiles(
            neckradius, tubeheight, totalheight, curvefraction, 
            ofhc_start_height, ofhc_end_height, outer_z, outer_r, inner_z, inner_r
        )

        if len(ofhc_outer_z) > 1:
            ofhc_outer_bound = g4.solid.GenericPolycone(
                "ofhc_cu_outer_bound", 0, 2 * np.pi, ofhc_outer_r, ofhc_outer_z, registry, "mm"
            )
            ofhc_inner_bound = g4.solid.GenericPolycone(
                "ofhc_cu_inner_bound", 0, 2 * np.pi, ofhc_inner_r, ofhc_inner_z, registry, "mm"
            )

            ofhc_solid = g4.solid.Subtraction(
                "ofhc_cu_solid",
                ofhc_outer_bound,
                ofhc_inner_bound,
                [[0, 0, 0], [0, 0, 0, "mm"]],
                registry,
            )

            ofhc_lv = g4.LogicalVolume(ofhc_solid, materials.metal_copper, "ofhc_cu_lv", registry)
            ofhc_lv.pygeom_color_rgba = [1.0, 0.5, 0.0, 1.0]
            ofhc_pv = g4.PhysicalVolume(
                [0, 0, 0], [0, 0, 0, "mm"], ofhc_lv, "ofhc_cu", tube_lv, registry=registry
            )

            print(f"\nOFHC copper layer:")
            print(f"  Height range: {ofhc_start_height}-{ofhc_end_height}mm from bottom")
            print(f"  Replaces steel in this region")

    # 316L stainless steel 
    if with_316l_ss:
        ss_outer_z, ss_outer_r, ss_inner_z, ss_inner_r = make_316l_ss_profiles(
            neckradius, tubeheight, totalheight, curvefraction, 
            ss_start_height, outer_z, outer_r, inner_z, inner_r
        )

        if len(ss_outer_z) > 1:
            ss_outer_bound = g4.solid.GenericPolycone(
                "ss_316l_outer_bound", 0, 2 * np.pi, ss_outer_r, ss_outer_z, registry, "mm"
            )
            ss_inner_bound = g4.solid.GenericPolycone(
                "ss_316l_inner_bound", 0, 2 * np.pi, ss_inner_r, ss_inner_z, registry, "mm"
            )

            ss_solid = g4.solid.Subtraction(
                "ss_316l_solid",
                ss_outer_bound,
                ss_inner_bound,
                [[0, 0, 0], [0, 0, 0, "mm"]],
                registry,
            )

            ss_lv = g4.LogicalVolume(ss_solid, materials.metal_steel, "ss_316l_lv", registry)
            ss_lv.pygeom_color_rgba = [0.7, 0.7, 0.8, 1.0]
            ss_pv = g4.PhysicalVolume(
                [0, 0, 0], [0, 0, 0, "mm"], ss_lv, "ss_316l", tube_lv, registry=registry
            )

            print(f"\n316L stainless steel layer:")
            print(f"  Height range: {ss_start_height}mm from bottom to top")
            print(f"  Replaces steel in this region")

    return tube_lv, cavity_lv, tube_pv, cavity_pv
