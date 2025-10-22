"""Standalone script to create and test reentrance tube geometry."""
from __future__ import annotations

import argparse

import numpy as np
import pyg4ometry as pg4
import pyg4ometry.geant4 as g4
import pint

from materials import OpticalMaterialRegistry
from reentrance_tube import construct_reentrance_tube
from pygeomtools import RemageDetectorInfo
from overlap_checker import check_overlaps, print_volume_hierarchy, summarize_geometry


def create_world_and_argon(reg, materials):
    """Create world volume with nested blackhole shell containing atmospheric argon."""
    
    # World - keeps vacuum as before
    world_solid = g4.solid.Box("world_solid", 10000, 10000, 10000, reg, "mm")
    world_lv = g4.LogicalVolume(world_solid, materials.vacuum, "world_logical", reg)
    reg.setWorld(world_lv.name)
    world_lv.pygeom_color_rgba = [1.0, 1.0, 1.0, 0.0]
    
    # Create blackhole material
    u = pint.get_application_registry()
    
    blackhole_mat = g4.Material(
        name="blackhole",
        density=1e-25,
        number_of_components=1,
        registry=reg,
    )
    hydrogen = g4.ElementSimple("H_blackhole", "H", 1, 1.008, reg)
    blackhole_mat.add_element_natoms(hydrogen, natoms=1)
    
    # Add optical properties to absorb photons
    wavelengths = np.array([100, 800]) * u.nm
    abs_length = np.array([0.001, 0.001]) * u.mm  # 1 micron absorption
    rindex = np.array([1.0, 1.0])
    
    with u.context("sp"):
        blackhole_mat.addVecPropertyPint("ABSLENGTH", wavelengths.to("eV"), abs_length)
        blackhole_mat.addVecPropertyPint("RINDEX", wavelengths.to("eV"), rindex)
    
    # Blackhole containment volume - slightly larger than AAr
    # This will be the mother volume for AAr
    blackhole_solid = g4.solid.Tubs(
        "blackhole_container_solid",
        0,          # Inner radius = 0 (solid cylinder)
        4000,       # Outer radius = 4m 
        9000,       # height = 9m 
        0,
        2 * np.pi,
        reg,
        "mm"
    )
    blackhole_lv = g4.LogicalVolume(
        blackhole_solid,
        blackhole_mat,
        "blackhole_container_lv",
        reg
    )
    blackhole_lv.pygeom_color_rgba = [0.0, 0.0, 0.0, 0.5]  # Semi-transparent black
    
    # Place blackhole container in world
    blackhole_pv = g4.PhysicalVolume(
        [0, 0, 0], [0, 0, 0, "mm"],
        blackhole_lv,
        "blackhole_container",
        world_lv,
        registry=reg
    )
    
    # Register blackhole as optical detector to count escaped photons
    blackhole_pv.set_pygeom_active_detector(RemageDetectorInfo("optical", 999))
    
    # Atmospheric argon - now placed INSIDE blackhole container
    atm_argon_solid = g4.solid.Tubs(
        "atmospheric_argon_solid",
        0,
        3000,       # 3m radius
        8000,       # 8m height
        0,
        2 * np.pi,
        reg,
        "mm"
    )
    atm_argon_lv = g4.LogicalVolume(
        atm_argon_solid,
        materials.liquidargon,
        "atmospheric_argon_lv",
        reg
    )
    atm_argon_lv.pygeom_color_rgba = [1.0, 0.75, 0.8, 0.3]
    
    # Place AAr inside blackhole container (not in world!)
    atm_argon_pv = g4.PhysicalVolume(
        [0, 0, 0], [0, 0, 0, "mm"],
        atm_argon_lv,
        "atmospheric_argon",
        blackhole_lv,  # Mother is now blackhole container!
        registry=reg
    )

    return world_lv, atm_argon_lv, atm_argon_pv


def main():
    parser = argparse.ArgumentParser(description="Generate reentrance tube with WLS, OFHC Cu, and 316L SS")
    parser.add_argument("--output", default="gdml/reentrance_tube.gdml", help="Output GDML filename")
    parser.add_argument("--outer-wls", action="store_true", help="Include outer WLS layer")
    parser.add_argument("--inner-wls", action="store_true", help="Include inner WLS layer")
    parser.add_argument("--ofhc-cu", action="store_true", help="Include OFHC copper cylinder")
    parser.add_argument("--ss-316l", action="store_true", help="Include 316L stainless steel cylinder")
    parser.add_argument(
        "--wls-height", type=float, default=2179, help="WLS height from bottom in mm (default: 2179)"
    )
    parser.add_argument("--check-overlaps", action="store_true", help="Check for geometry overlaps")
    parser.add_argument("--overlap-tolerance", type=float, default=0.01, help="Overlap tolerance in mm (default: 0.01)")
    parser.add_argument("--print-hierarchy", action="store_true", help="Print volume hierarchy")
    parser.add_argument("--summary", action="store_true", help="Print geometry summary")
    args = parser.parse_args()

    print("Creating reentrance tube geometry...")

    # Create registry and materials
    reg = g4.Registry()
    materials = OpticalMaterialRegistry(reg)

    # Create world and argon
    world_lv, atm_argon_lv, atm_argon_pv = create_world_and_argon(reg, materials)

    # Construct reentrance tube
    tube_lv, cavity_lv, tube_pv, cavity_pv = construct_reentrance_tube(
        materials=materials,
        registry=reg,
        mother_lv=atm_argon_lv,
        mother_pv=atm_argon_pv,
        with_outer_wls=args.outer_wls,
        with_inner_wls=args.inner_wls,
        wls_height=args.wls_height,
        with_ofhc_cu=args.ofhc_cu,
        with_316l_ss=args.ss_316l,
    )

    # Geometry analysis
    if args.summary:
        summarize_geometry(reg)
    
    if args.print_hierarchy:
        print_volume_hierarchy(reg, max_depth=5)
    
    if args.check_overlaps:
        print("\nChecking for overlaps...")
        overlaps = check_overlaps(reg, tolerance=args.overlap_tolerance, verbose=True)
        if len(overlaps) > 0:
            print(f"\n⚠ WARNING: {len(overlaps)} overlap(s) detected!")
            print("Review the overlap details above.")
        else:
            print("\n✓ No overlaps detected!")

    print(f"\nExporting to {args.output}...")
    w = pg4.gdml.Writer()
    w.addDetector(reg)
    w.write(args.output)
    print("GDML file created successfully!")
    
    # Additional note about Geant4 overlap checking
    if args.check_overlaps:
        print("\nNote: For comprehensive overlap checking, you can also:")
        print("  1. Load the GDML file in Geant4")
        print("  2. Run: /geometry/test/run")
        print("  This will perform a more thorough overlap check.")


if __name__ == "__main__":
    main()
