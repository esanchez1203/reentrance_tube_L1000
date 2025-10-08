"""Standalone script to create and test reentrance tube geometry."""
from __future__ import annotations

import argparse

import numpy as np
import pyg4ometry as pg4
import pyg4ometry.geant4 as g4

from materials import OpticalMaterialRegistry
from reentrance_tube import construct_reentrance_tube


def create_world_and_argon(reg, materials):
    """Create world volume and atmospheric argon."""
    # World
    world_solid = g4.solid.Box("world_solid", 10000, 10000, 10000, reg, "mm")
    world_lv = g4.LogicalVolume(world_solid, materials.vacuum, "world_logical", reg)
    reg.setWorld(world_lv.name)
    world_lv.pygeom_color_rgba = [1.0, 1.0, 1.0, 0.0]

    # Atmospheric argon
    atm_argon_solid = g4.solid.Tubs("atmospheric_argon_solid", 0, 3000, 8000, 0, 2 * np.pi, reg, "mm")
    atm_argon_lv = g4.LogicalVolume(atm_argon_solid, materials.liquidargon, "atmospheric_argon_lv", reg)
    atm_argon_lv.pygeom_color_rgba = [1.0, 0.75, 0.8, 0.3]
    atm_argon_pv = g4.PhysicalVolume(
        [0, 0, 0], [0, 0, 0, "mm"], atm_argon_lv, "atmospheric_argon", world_lv, registry=reg
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

    print(f"\nExporting to {args.output}...")
    w = pg4.gdml.Writer()
    w.addDetector(reg)
    w.write(args.output)
    print("GDML file created successfully!")


if __name__ == "__main__":
    main()
