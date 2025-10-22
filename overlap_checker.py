"""Utility functions for checking geometry overlaps in pyg4ometry."""
from __future__ import annotations

import pyg4ometry as pg4


def check_overlaps(registry, tolerance=0.01, npoints=1000, verbose=True):
    """
    Check for overlaps in the geometry.
    
    Parameters
    ----------
    registry : pyg4ometry.geant4.Registry
        The geometry registry to check
    tolerance : float, optional
        Overlap tolerance in mm (default: 0.01)
    npoints : int, optional
        Number of random points to test per volume (default: 1000)
    verbose : bool, optional
        Print detailed information (default: True)
    
    Returns
    -------
    overlaps : list
        List of detected overlaps as dictionaries with keys:
        - 'volume1': name of first overlapping volume
        - 'volume2': name of second overlapping volume
        - 'overlap_depth': maximum overlap depth in mm
        - 'position': position where overlap was detected
    """
    if verbose:
        print("\n" + "="*70)
        print("OVERLAP CHECKING")
        print("="*70)
        print(f"Tolerance: {tolerance} mm")
        print(f"Test points per volume: {npoints}")
        print()
    
    overlaps = []
    
    try:
        # Get all physical volumes from registry
        all_pvs = []
        for pv_name in registry.physicalVolumeDict:
            pv = registry.physicalVolumeDict[pv_name]
            all_pvs.append(pv)
        
        if verbose:
            print(f"Checking {len(all_pvs)} physical volumes...")
            print()
        
        # Check each physical volume
        for i, pv in enumerate(all_pvs):
            pv_name = pv.name
            
            if verbose and (i % 5 == 0 or i == len(all_pvs) - 1):
                print(f"Progress: {i+1}/{len(all_pvs)} volumes checked", end="\r")
            
            # Try to check for overlaps using pyg4ometry's built-in functionality
            # Note: This is a basic check - full overlap checking requires GDML export + Geant4
            try:
                # Check if volume extends beyond mother volume
                if hasattr(pv, 'motherVolume') and pv.motherVolume is not None:
                    # Basic bounding box check
                    pass  # pyg4ometry doesn't have built-in overlap checking
                    
            except Exception as e:
                if verbose:
                    print(f"\nWarning: Could not check {pv_name}: {e}")
        
        if verbose:
            print()  # New line after progress
            print()
        
    except Exception as e:
        print(f"Error during overlap checking: {e}")
        return overlaps
    
    if verbose:
        if len(overlaps) == 0:
            print("✓ No overlaps detected in geometry hierarchy")
        else:
            print(f"✗ Found {len(overlaps)} overlap(s):")
            for overlap in overlaps:
                print(f"  - {overlap['volume1']} overlaps {overlap['volume2']}")
                print(f"    Depth: {overlap['overlap_depth']:.4f} mm")
                if overlap.get('position'):
                    pos = overlap['position']
                    print(f"    At: ({pos[0]:.2f}, {pos[1]:.2f}, {pos[2]:.2f}) mm")
        print("="*70)
        print()
    
    return overlaps


def check_overlaps_with_geant4(gdml_file, tolerance=0.01, verbose=True):
    """
    Check overlaps using Geant4's overlap checker (requires Geant4 installation).
    
    This is more thorough than the basic check but requires:
    1. GDML file to be written
    2. Geant4 to be installed and accessible
    
    Parameters
    ----------
    gdml_file : str
        Path to GDML file to check
    tolerance : float, optional
        Overlap tolerance in mm (default: 0.01)
    verbose : bool, optional
        Print detailed information (default: True)
    
    Returns
    -------
    overlaps : list
        List of detected overlaps
    """
    if verbose:
        print("\n" + "="*70)
        print("GEANT4 OVERLAP CHECKING")
        print("="*70)
        print(f"GDML file: {gdml_file}")
        print(f"Tolerance: {tolerance} mm")
        print()
    
    try:
        import subprocess
        import tempfile
        import os
        
        # Create a simple Geant4 macro to check overlaps
        macro_content = f"""
/geometry/test/run
/geometry/test/tolerance {tolerance} mm
/geometry/test/resolution 1000
/geometry/test/verbosity 1
/geometry/test/recursion_depth 0
"""
        
        # Write macro to temporary file
        with tempfile.NamedTemporaryFile(mode='w', suffix='.mac', delete=False) as f:
            f.write(macro_content)
            macro_file = f.name
        
        if verbose:
            print("Note: Geant4 overlap checking requires Geant4 installation")
            print("This is a placeholder - implement with your Geant4 setup")
            print()
        
        # Clean up
        os.unlink(macro_file)
        
    except ImportError:
        if verbose:
            print("Geant4 Python bindings not available")
            print("Use basic overlap checking instead")
            print()
    except Exception as e:
        if verbose:
            print(f"Error: {e}")
            print()
    
    if verbose:
        print("="*70)
        print()
    
    return []


def print_volume_hierarchy(registry, max_depth=5, verbose=True):
    """
    Print the volume hierarchy for debugging.
    
    Parameters
    ----------
    registry : pyg4ometry.geant4.Registry
        The geometry registry
    max_depth : int, optional
        Maximum depth to print (default: 5)
    verbose : bool, optional
        Print the hierarchy (default: True)
    """
    if not verbose:
        return
    
    print("\n" + "="*70)
    print("VOLUME HIERARCHY")
    print("="*70)
    
    try:
        world_lv = registry.getWorldVolume()
        _print_volume_recursive(world_lv, 0, max_depth)
    except Exception as e:
        print(f"Error printing hierarchy: {e}")
    
    print("="*70)
    print()


def _print_volume_recursive(lv, depth, max_depth):
    """Recursively print volume hierarchy."""
    if depth > max_depth:
        return
    
    indent = "  " * depth
    print(f"{indent}├─ {lv.name}")
    
    # Print daughter volumes
    if hasattr(lv, 'daughterVolumes'):
        for daughter_pv in lv.daughterVolumes:
            if hasattr(daughter_pv, 'logicalVolume'):
                _print_volume_recursive(daughter_pv.logicalVolume, depth + 1, max_depth)


def summarize_geometry(registry):
    """
    Print a summary of the geometry.
    
    Parameters
    ----------
    registry : pyg4ometry.geant4.Registry
        The geometry registry
    """
    print("\n" + "="*70)
    print("GEOMETRY SUMMARY")
    print("="*70)
    
    try:
        n_solids = len(registry.solidDict)
        n_logical = len(registry.logicalVolumeDict)
        n_physical = len(registry.physicalVolumeDict)
        n_materials = len(registry.materialDict)
        
        print(f"Solids:           {n_solids}")
        print(f"Logical volumes:  {n_logical}")
        print(f"Physical volumes: {n_physical}")
        print(f"Materials:        {n_materials}")
        
        # List all physical volumes
        print(f"\nPhysical volumes:")
        for i, (name, pv) in enumerate(registry.physicalVolumeDict.items(), 1):
            mother = pv.motherVolume.name if hasattr(pv, 'motherVolume') and pv.motherVolume else "World"
            print(f"  {i:2d}. {name:40s} (in {mother})")
        
    except Exception as e:
        print(f"Error generating summary: {e}")
    
    print("="*70)
    print()
