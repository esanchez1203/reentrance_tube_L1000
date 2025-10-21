"""Profile generation functions for reentrance tube geometry - TPB as parent with 1μm coating"""
from __future__ import annotations

# Constants
WLSR_TPB_THICKNESS = 1 * 1e-3  # 1 um TPB coating (in mm)
WLSR_TTX_THICKNESS = 254 * 1e-3  # 254 um Tetratex foil (in mm)
WLSR_THICKNESS = WLSR_TPB_THICKNESS + WLSR_TTX_THICKNESS
PROTECTION_GAP = 5 * 1e-6  # 5 nm gap in mm
PROTECTION_GAP_LAYER = 100 * 1e-6  # 5 nm gap in mm


def get_thickness_at_distance(
    distance_from_top: float,
    include_wlsr: bool = False,
    wlsr_height: float = 2179,
    total_tube_height: float = 6247,
) -> float:
    """Calculate steel thickness at given distance from top of tube."""
    if distance_from_top <= 4067:
        progress = distance_from_top / 4067
        steel_thickness = 6.0 - progress * 0.0  # Constant thickness
    else:
        remaining = distance_from_top - 4067
        max_remaining = total_tube_height - 4067
        if max_remaining > 0:
            progress = remaining / max_remaining
            steel_thickness = max(1.5 - progress * 0.0, 1.5)
        else:
            steel_thickness = 1.5

    return steel_thickness

def ensure_closed_bottom(z_list, r_list, bottom_z, closure_thickness=None, gap_threshold=0.000001):
    """
    Ensure proper bottom closure by adding intermediate conical points.
    
    Args:
        z_list: list of z coordinates
        r_list: list of r coordinates  
        bottom_z: z coordinate where closure should occur
        closure_thickness: if provided, close at bottom_z + closure_thickness
        gap_threshold: minimum gap size (in mm) to trigger point insertion (default 0.0001 = 100 nm)
    """
    if closure_thickness is not None:
        closure_z = bottom_z + closure_thickness
    else:
        closure_z = bottom_z
    
    # Find first point with r > 0
    first_nonzero_idx = None
    for i, r in enumerate(r_list):
        if r > 0:
            first_nonzero_idx = i
            break
    
    if first_nonzero_idx is None:
        return z_list, r_list
    
    first_z = z_list[first_nonzero_idx]
    first_r = r_list[first_nonzero_idx]
    
    # If there's a gap between closure and first point, add intermediate points
    z_gap = first_z - closure_z
    if z_gap > gap_threshold:
        # Add conical transition points with adaptive spacing
        # For small gaps (< 0.01mm), use finer spacing; for large gaps use 0.1mm
        if z_gap < 0.01:  # Less than 10 μm
            point_spacing = z_gap / 5  # At least 5 points
        else:
            point_spacing = 0.01  # 10 μm spacing for larger gaps
        
        n_points = max(3, int(z_gap / point_spacing))
        new_points_z = []
        new_points_r = []
        
        for i in range(1, n_points):
            frac = i / n_points
            z_intermediate = closure_z + frac * z_gap
            r_intermediate = frac * first_r  # Linear interpolation for cone
            new_points_z.append(z_intermediate)
            new_points_r.append(r_intermediate)
        
        # Insert new points after closure point (index 0)
        for i, (z, r) in enumerate(zip(new_points_z, new_points_r)):
            z_list.insert(i + 1, z)
            r_list.insert(i + 1, r)
    
    return z_list, r_list

def make_outer_profile(
    neckradius: float,
    tubeheight: float,
    totalheight: float,
    curvefraction: float,
    wls_height: float = 2179,
) -> tuple[list, list]:
    """Create complete outer profile with filled cylindrical section."""
    z = []
    r = []

    top_z = totalheight - 1
    curve_start_z = totalheight - tubeheight * (1 - curvefraction)
    bottom_z = totalheight - tubeheight

    z.append(top_z)
    r.append(0)
    z.append(top_z - 0.00001)
    r.append(neckradius)
    z.append(curve_start_z)
    r.append(neckradius)

    curve_z_fractions = [
        16.39, 26.23, 33.61, 40.98, 47.54, 53.28, 55.74, 59.02,
        64.75, 67.21, 70.49, 74.59, 78.69, 80.33, 82.79, 86.07,
        88.52, 90.98, 92.62, 94.26, 95.90, 97.54, 99.18, 99.88,
    ]

    curve_r_fractions = [
        1.89, 4.17, 6.82, 10.23, 14.02, 18.56, 21.21, 23.86,
        28.79, 31.82, 34.47, 39.39, 44.32, 46.97, 49.62, 54.55,
        59.47, 64.39, 68.56, 73.11, 76.52, 82.58, 88.64, 93.94,
    ]

    for i in range(len(curve_z_fractions)):
        z_frac = curve_z_fractions[i] * 0.01
        r_frac = curve_r_fractions[i] * 0.01
        z_pos = curve_start_z - (z_frac * curvefraction * tubeheight)
        r_pos = neckradius * (1 - r_frac)
        z.append(z_pos)
        r.append(r_pos)

    z.append(bottom_z)
    r.append(0)

    # Adjust z coordinates
    z = [zi - 5000 for zi in z]
    z.reverse()
    r.reverse()

    # Sort by z coordinate
    combined = sorted(zip(z, r), key=lambda x: x[0])
    z_sorted, r_sorted = zip(*combined)
    z_out, r_out = list(z_sorted), list(r_sorted)

    # Fill large gaps in cylindrical section
    for i in range(len(z_out) - 1):
        gap = z_out[i + 1] - z_out[i]
        if gap > 100:
            z_filled = list(z_out[: i + 1])
            r_filled = list(r_out[: i + 1])
            current_z = z_out[i] + 50
            while current_z < z_out[i + 1]:
                z_filled.append(current_z)
                r_filled.append(neckradius)
                current_z += 50
            z_filled.extend(z_out[i + 1 :])
            r_filled.extend(r_out[i + 1 :])
            z_out, r_out = z_filled, r_filled
            break

    # DISABLED: Adding critical points causes geometry crossing issues
    # The WLSR profiles should work with existing tube points
    # If needed, enable this section and debug the crossing issue
    
    # Add points at critical WLS boundaries
    wls_top_z = bottom_z + wls_height
    wls_margin = 0.001  # 1 μm in mm
    
    critical_z_values = [
        wls_top_z - wls_margin,
        wls_top_z,
        wls_top_z + wls_margin,
    ]
    
    for z_crit in critical_z_values:
        if bottom_z < z_crit < top_z:
            exists = any(abs(z_existing - z_crit) < 0.0001 for z_existing in z_out)
            if not exists:
                r_crit = None
                for i in range(len(z_out) - 1):
                    if z_out[i] < z_crit < z_out[i + 1]:
                        frac = (z_crit - z_out[i]) / (z_out[i + 1] - z_out[i])
                        r_crit = r_out[i] + frac * (r_out[i + 1] - r_out[i])
                        break
                
                if r_crit is not None and r_crit > 0:
                    z_out.append(z_crit)
                    r_out.append(r_crit)
    
    if len(z_out) > len(list(z_sorted)):
        combined = sorted(zip(z_out, r_out), key=lambda x: x[0])
        z_sorted, r_sorted = zip(*combined)
        z_out, r_out = list(z_sorted), list(r_sorted)
    
    z_out, r_out = ensure_closed_bottom(z_out, r_out, z_out[0])
    
    return z_out, r_out


def make_inner_profile(
    neckradius: float,
    tubeheight: float,
    totalheight: float,
    curvefraction: float,
    wls_height: float = 2179,
) -> tuple[list, list]:
    """
    Create inner profile with constant wall thickness.
    Inner surface closes at higher z than outer to maintain thickness.
    """
    outer_z, outer_r = make_outer_profile(neckradius, tubeheight, totalheight, curvefraction, wls_height)
    top_z_original = totalheight - 1
    bottom_z = totalheight - tubeheight
    
    # Calculate the thickness offset needed at the bottom
    bottom_thickness = get_thickness_at_distance(totalheight - 1 - bottom_z)
    
    inner_z, inner_r = [], []

    for z, r in zip(outer_z, outer_r):
        z_original = z + 5000
        dist_from_top = top_z_original - z_original
        thickness = get_thickness_at_distance(dist_from_top)

        if r == 0:
            # Skip r=0 points during iteration - closures added separately
            continue
        elif r == neckradius:
            # Cylindrical section
            inner_r_value = max(0, r - thickness)
            inner_z.append(z)
            inner_r.append(inner_r_value)
        else:
            # Curved section - scale thickness by radius ratio
            radius_ratio = r / neckradius
            scaled_thickness = thickness * radius_ratio
            inner_r_value = max(0, r - scaled_thickness)
            inner_z.append(z)
            inner_r.append(inner_r_value)

    # Add closures explicitly
    # Bottom closure: inner closes ABOVE outer by thickness offset
    bottom_z_adjusted = (bottom_z - 5000)
    inner_bottom_z = bottom_z_adjusted + bottom_thickness
    
    inner_z.insert(0, inner_bottom_z)
    inner_r.insert(0, 0)
    
    #close tube with 6 mm 
    top_thickness = get_thickness_at_distance(0)  # 6.0 mm at the top
    
    inner_z.append(top_z_original - 5000 - top_thickness)
    inner_r.append(0)
    
    
    # Ensure proper conical closure at bottom
    inner_z, inner_r = ensure_closed_bottom(inner_z, inner_r, inner_bottom_z)

    return inner_z, inner_r

def make_inner_wlsr_argon_profiles(
    neckradius: float,
    tubeheight: float,
    totalheight: float,
    curvefraction: float,
    wls_height: float,
    inner_z: list,
    inner_r: list,
    outer_z: list,
    outer_r: list,
) -> tuple:
    """
    Structure: Steel inner -> Gap (5nm) -> TPB (mother) -> TTX (daughter, 1μm inside TPB) -> UAr
    
    Creates closed volume with TPB coating on all sides (top, bottom, and radial).
    """
    
    bottom_z = (totalheight - tubeheight) - 5000
    top_wls_z = bottom_z + wls_height
    
    # TPB profiles (mother volume - full thickness)
    tpb_outer_z, tpb_outer_r = [], []
    tpb_inner_z, tpb_inner_r = [], []
    
    # TTX profiles (daughter volume - 1μm inside TPB on all sides)
    ttx_outer_z, ttx_outer_r = [], []
    ttx_inner_z, ttx_inner_r = [], []

    # Collect points in WLSR region
    points_collected = 0
    for z, r in zip(inner_z, inner_r):
        # For negative z-coordinates: top_wls_z is HIGHER than bottom_z
        # Collect all points in the WLS region
        if top_wls_z > z >= bottom_z and r > 0:
            # Start from steel inner surface, offset by protection gap
            base = r - PROTECTION_GAP
            
            # TPB (mother) - contains TTX plus 1μm coating on all sides
            # Total TPB thickness = TTX + 2×1μm
            tpb_total_thickness = WLSR_TTX_THICKNESS + 2 * WLSR_TPB_THICKNESS
            tpb_outer = base
            tpb_inner = base - tpb_total_thickness
            
            # TTX (daughter) - 1μm inside TPB on all sides
            ttx_outer = base - WLSR_TPB_THICKNESS
            ttx_inner = base - tpb_total_thickness + WLSR_TPB_THICKNESS
            
            # Add points if valid
            if tpb_inner > 0 and ttx_inner > 0 and ttx_inner < ttx_outer:
                tpb_outer_z.append(z)
                tpb_outer_r.append(tpb_outer)
                tpb_inner_z.append(z)
                tpb_inner_r.append(tpb_inner)
                
                ttx_outer_z.append(z)
                ttx_outer_r.append(ttx_outer)
                ttx_inner_z.append(z)
                ttx_inner_r.append(ttx_inner)
                
                points_collected += 1
    
    print(f"\nInner WLSR: Collected {points_collected} data points")
    print(f"  WLS region: bottom_z={bottom_z:.3f} mm to top_wls_z={top_wls_z:.3f} mm")
    if tpb_outer_z:
        print(f"  Actual data z range: [{min(tpb_outer_z):.3f}, {max(tpb_outer_z):.3f}]")
    else:
        print(f"  WARNING: No data points collected!")

    # Add bottom closures with z-offsets
    # For inner: layers grow INWARD, so close at progressively HIGHER z
    if tpb_outer_z:
        # Find where steel inner closes
        steel_inner_bottom_z = min([z for z, r in zip(inner_z, inner_r) if r == 0])
        
        # TPB outer closes first (lowest z) - right after the gap
        tpb_outer_bottom_z = steel_inner_bottom_z + PROTECTION_GAP
        tpb_outer_z.insert(0, tpb_outer_bottom_z)
        tpb_outer_r.insert(0, 0)
        
        # TPB inner closes higher by total thickness
        tpb_total_thickness = WLSR_TTX_THICKNESS + 2 * WLSR_TPB_THICKNESS
        tpb_inner_bottom_z = tpb_outer_bottom_z + tpb_total_thickness
        tpb_inner_z.insert(0, tpb_inner_bottom_z)
        tpb_inner_r.insert(0, 0)
        
        # TTX outer is 1μm (TPB thickness) above TPB outer
        ttx_outer_bottom_z = tpb_outer_bottom_z + WLSR_TPB_THICKNESS
        ttx_outer_z.insert(0, ttx_outer_bottom_z)
        ttx_outer_r.insert(0, 0)
        
        # TTX inner is 1μm (TPB thickness) below TPB inner
        ttx_inner_bottom_z = tpb_inner_bottom_z - WLSR_TPB_THICKNESS
        ttx_inner_z.insert(0, ttx_inner_bottom_z)
        ttx_inner_r.insert(0, 0)

        tpb_inner_z, tpb_inner_r = ensure_closed_bottom(tpb_inner_z, tpb_inner_r, tpb_inner_bottom_z)
        tpb_outer_z, tpb_outer_r = ensure_closed_bottom(tpb_outer_z, tpb_outer_r, tpb_outer_bottom_z)
        ttx_inner_z, ttx_inner_r = ensure_closed_bottom(ttx_inner_z, ttx_inner_r, ttx_inner_bottom_z)
        ttx_outer_z, ttx_outer_r = ensure_closed_bottom(ttx_outer_z, ttx_outer_r, ttx_outer_bottom_z)


        # Debug output
        print(f"\nBottom closure z-positions:")
        print(f"  Steel inner closes at: {steel_inner_bottom_z:.6f} mm")
        print(f"  TPB outer: {tpb_outer_bottom_z:.6f} mm (offset: {(tpb_outer_bottom_z-steel_inner_bottom_z)*1e3:.3f} μm)")
        print(f"  TPB inner: {tpb_inner_bottom_z:.6f} mm (offset: {(tpb_inner_bottom_z-tpb_outer_bottom_z)*1e3:.3f} μm)")
        print(f"  TTX outer: {ttx_outer_bottom_z:.6f} mm (offset from TPB outer: {(ttx_outer_bottom_z-tpb_outer_bottom_z)*1e3:.3f} μm)")
        print(f"  TTX inner: {ttx_inner_bottom_z:.6f} mm (offset from TPB inner: {(tpb_inner_bottom_z-ttx_inner_bottom_z)*1e3:.3f} μm)")

 
        # All layers close at the same z at the top
    tpb_outer_z.append(top_wls_z)
    tpb_outer_r.append(0)
    tpb_inner_z.append(top_wls_z+10.0)
    tpb_inner_r.append(0)
    ttx_outer_z.append(top_wls_z-WLSR_TPB_THICKNESS)
    ttx_outer_r.append(0)
    ttx_inner_z.append(top_wls_z+10.0)
    ttx_inner_r.append(0)
        
    return (tpb_outer_z, tpb_outer_r, tpb_inner_z, tpb_inner_r,
            ttx_outer_z, ttx_outer_r, ttx_inner_z, ttx_inner_r)

def make_outer_wlsr_atmospheric_profiles(
    neckradius: float,
    tubeheight: float,
    totalheight: float,
    curvefraction: float,
    wls_height: float,
    outer_z: list,
    outer_r: list,
):
    """
    Create WLS profiles for outer (atmospheric) region.
    TPB is 1 μm thick shell coating that surrounds TTX completely.
    
    Structure: Steel outer -> Gap (5nm, argon) -> TPB shell (1 μm, parent) -> TTX core (daughter, 1 μm smaller on all sides) -> Atmospheric Argon
    """
    bottom_z = (totalheight - tubeheight) - 5000
    top_wls_z = bottom_z + wls_height

    # TPB profiles (parent/mother volume)
    tpb_outer_z, tpb_outer_r = [], []
    tpb_inner_z, tpb_inner_r = [], []
    
    # TTX profiles (daughter volume - smaller by 1 μm on all sides)
    ttx_outer_z, ttx_outer_r = [], []
    ttx_inner_z, ttx_inner_r = [], []

    # Collect points in WLSR region
    points_collected = 0
    for z, r in zip(outer_z, outer_r):
        # For negative z-coordinates: top_wls_z is HIGHER than bottom_z
        # Collect all points in the WLS region
        if top_wls_z > z >= bottom_z and r > 0:
            # Base position: steel outer + gap
            base = r + PROTECTION_GAP

            # TPB inner surface at base
            tpb_inner = base
            # TPB total thickness to accommodate TTX + coating on both sides
            total_thickness = WLSR_TTX_THICKNESS + 2 * WLSR_TPB_THICKNESS
            tpb_outer = tpb_inner + total_thickness

            # TTX sits inside with 1 μm clearance on all sides
            # TTX inner surface is 1 μm outward from TPB inner
            ttx_inner = tpb_inner + WLSR_TPB_THICKNESS
            # TTX outer surface is 1 μm inward from TPB outer
            ttx_outer = tpb_outer - WLSR_TPB_THICKNESS

            # Add points
            tpb_inner_z.append(z)
            tpb_inner_r.append(tpb_inner)
            tpb_outer_z.append(z)
            tpb_outer_r.append(tpb_outer)
            
            ttx_inner_z.append(z)
            ttx_inner_r.append(ttx_inner)
            ttx_outer_z.append(z)
            ttx_outer_r.append(ttx_outer)
            
            points_collected += 1
    
    print(f"\nOuter WLSR: Collected {points_collected} data points")
    print(f"  WLS region: bottom_z={bottom_z:.3f} mm to top_wls_z={top_wls_z:.3f} mm")
    if tpb_outer_z:
        print(f"  Actual data z range: [{min(tpb_outer_z):.3f}, {max(tpb_outer_z):.3f}]")
    else:
        print(f"  WARNING: No data points collected!")

    # Add bottom closures with z-offsets
    # For outer: layers grow OUTWARD, so close at progressively HIGHER z
    if tpb_outer_z:
        # Find where steel outer closes
        steel_outer_bottom_z = min([z for z, r in zip(outer_z, outer_r) if r == 0])
        
        # TPB inner closes first (lowest z) - right after the gap
        tpb_inner_bottom_z = steel_outer_bottom_z + PROTECTION_GAP
        tpb_inner_z.insert(0, tpb_inner_bottom_z)
        tpb_inner_r.insert(0, 0)
        
        # TPB outer closes higher by total thickness
        tpb_total_thickness = WLSR_TTX_THICKNESS + 2 * WLSR_TPB_THICKNESS
        tpb_outer_bottom_z = tpb_inner_bottom_z + tpb_total_thickness
        tpb_outer_z.insert(0, tpb_outer_bottom_z)
        tpb_outer_r.insert(0, 0)
        
        # TTX inner is 1μm (TPB thickness) above TPB inner
        ttx_inner_bottom_z = tpb_inner_bottom_z + WLSR_TPB_THICKNESS
        ttx_inner_z.insert(0, ttx_inner_bottom_z)
        ttx_inner_r.insert(0, 0)
        
        # TTX outer is 1μm (TPB thickness) below TPB outer
        ttx_outer_bottom_z = tpb_outer_bottom_z - WLSR_TPB_THICKNESS
        ttx_outer_z.insert(0, ttx_outer_bottom_z)
        ttx_outer_r.insert(0, 0)

        # Add ensure_closed_bottom for all layers (same as inner)
        tpb_inner_z, tpb_inner_r = ensure_closed_bottom(tpb_inner_z, tpb_inner_r, tpb_inner_bottom_z)
        tpb_outer_z, tpb_outer_r = ensure_closed_bottom(tpb_outer_z, tpb_outer_r, tpb_outer_bottom_z)
        ttx_inner_z, ttx_inner_r = ensure_closed_bottom(ttx_inner_z, ttx_inner_r, ttx_inner_bottom_z)
        ttx_outer_z, ttx_outer_r = ensure_closed_bottom(ttx_outer_z, ttx_outer_r, ttx_outer_bottom_z)

        # Debug output
        print(f"\nBottom closure z-positions:")
        print(f"  Steel outer closes at: {steel_outer_bottom_z:.6f} mm")
        print(f"  TPB inner: {tpb_inner_bottom_z:.6f} mm (offset: {(tpb_inner_bottom_z-steel_outer_bottom_z)*1e3:.3f} μm)")
        print(f"  TPB outer: {tpb_outer_bottom_z:.6f} mm (offset: {(tpb_outer_bottom_z-tpb_inner_bottom_z)*1e3:.3f} μm)")
        print(f"  TTX inner: {ttx_inner_bottom_z:.6f} mm (offset from TPB inner: {(ttx_inner_bottom_z-tpb_inner_bottom_z)*1e3:.3f} μm)")
        print(f"  TTX outer: {ttx_outer_bottom_z:.6f} mm (offset from TPB outer: {(tpb_outer_bottom_z-ttx_outer_bottom_z)*1e3:.3f} μm)")

    # Add top closures - same style as inner WLSR
    # TPB inner closes at top_wls_z
    tpb_inner_z.append(top_wls_z + 10.0)
    tpb_inner_r.append(0)
    
    # TPB outer extends beyond by 10mm 
    tpb_outer_z.append(top_wls_z)
    tpb_outer_r.append(0)
    
    # TTX inner closes 1μm below TPB inner
    ttx_inner_z.append(top_wls_z + 10.0)
    ttx_inner_r.append(0)
    
    # TTX outer extends beyond by 10mm
    ttx_outer_z.append(top_wls_z - WLSR_TPB_THICKNESS)
    ttx_outer_r.append(0)
    
    print(f"\nTop closure: Same style as inner WLSR")
    print(f"  TPB inner: {top_wls_z:.6f} mm")
    print(f"  TTX inner: {top_wls_z - WLSR_TPB_THICKNESS:.6f} mm")
    print(f"  TPB outer: {top_wls_z + 10.0:.6f} mm (extended)")
    print(f"  TTX outer: {top_wls_z + 10.0:.6f} mm (extended)")
    print(f"  (Subtraction creates open top automatically)")

    return (
        tpb_outer_z, tpb_outer_r,
        tpb_inner_z, tpb_inner_r,
        ttx_outer_z, ttx_outer_r,
        ttx_inner_z, ttx_inner_r,
    )

def make_ofhc_cu_profiles(
    neckradius: float,
    tubeheight: float,
    totalheight: float,
    curvefraction: float,
    ofhc_start_height: float,
    ofhc_end_height: float,
    outer_z: list,
    outer_r: list,
    inner_z: list,
    inner_r: list,
) -> tuple:
    """Create OFHC copper profiles as SOLID volume with protection gap from steel."""
    bottom_z = (totalheight - tubeheight) - 5000
    ofhc_start_z = bottom_z + ofhc_start_height
    ofhc_end_z = bottom_z + ofhc_end_height

    ofhc_outer_z, ofhc_outer_r = [], []
    ofhc_inner_z, ofhc_inner_r = [], []

    # Build a lookup for inner_r by z
    inner_r_by_z = {z: r for z, r in zip(inner_z, inner_r)}

    for z_out, r_out in zip(outer_z, outer_r):
        if ofhc_start_z <= z_out <= ofhc_end_z:
            # Find corresponding inner r at this z
            if z_out in inner_r_by_z:
                # OFHC outer is inward from steel outer by protection gap
                ofhc_outer_z.append(z_out)
                ofhc_outer_r.append(r_out - PROTECTION_GAP_LAYER)
                
                # OFHC inner is outward from steel inner by protection gap
                ofhc_inner_z.append(z_out)
                ofhc_inner_r.append(inner_r_by_z[z_out] + PROTECTION_GAP_LAYER)

    # Ensure boundaries
    if ofhc_outer_z and ofhc_outer_z[0] > ofhc_start_z + 1:
        ofhc_outer_z.insert(0, ofhc_start_z)
        ofhc_outer_r.insert(0, neckradius - PROTECTION_GAP_LAYER)
        inner_start_r = neckradius - get_thickness_at_distance(tubeheight - ofhc_start_height)
        ofhc_inner_z.insert(0, ofhc_start_z)
        ofhc_inner_r.insert(0, inner_start_r + PROTECTION_GAP_LAYER)

    if ofhc_outer_z and ofhc_outer_z[-1] < ofhc_end_z - 1:
        ofhc_outer_z.append(ofhc_end_z)
        ofhc_outer_r.append(neckradius - PROTECTION_GAP_LAYER)
        inner_end_r = neckradius - get_thickness_at_distance(tubeheight - ofhc_end_height)
        ofhc_inner_z.append(ofhc_end_z)
        ofhc_inner_r.append(inner_end_r + PROTECTION_GAP_LAYER)

    # Add r=0 closure caps at both ends
    ofhc_outer_z.insert(0, ofhc_start_z + PROTECTION_GAP_LAYER )
    ofhc_outer_r.insert(0, 0)
    ofhc_inner_z.insert(0, ofhc_start_z)
    ofhc_inner_r.insert(0, 0)

    ofhc_outer_z.append(ofhc_end_z - PROTECTION_GAP_LAYER )
    ofhc_outer_r.append(0)
    ofhc_inner_z.append(ofhc_end_z)
    ofhc_inner_r.append(0)

    return ofhc_outer_z, ofhc_outer_r, ofhc_inner_z, ofhc_inner_r

def make_316l_ss_profiles(
    neckradius: float,
    tubeheight: float,
    totalheight: float,
    curvefraction: float,
    ss_start_height: float,
    outer_z: list,
    outer_r: list,
    inner_z: list,
    inner_r: list,
) -> tuple:
    """Create 316L stainless steel profiles with protection gap on all sides."""
    bottom_z = (totalheight - tubeheight) - 5000
    ss_start_z = bottom_z + ss_start_height
    ss_end_z = (totalheight - 1) - 5000

    ss_outer_z, ss_outer_r = [], []
    ss_inner_z, ss_inner_r = [], []

    # Build a lookup for inner_r by z
    inner_r_by_z = {z: r for z, r in zip(inner_z, inner_r)}

    for z_out, r_out in zip(outer_z, outer_r):
        if ss_start_z <= z_out <= ss_end_z and r_out > 0:
            # Find corresponding inner r at this z
            if z_out in inner_r_by_z:
                # SS outer is inward from outer profile by protection gap
                ss_outer_z.append(z_out)
                ss_outer_r.append(r_out - PROTECTION_GAP_LAYER)
                
                # SS inner is outward from inner profile by protection gap
                ss_inner_z.append(z_out)
                ss_inner_r.append(inner_r_by_z[z_out] + PROTECTION_GAP_LAYER)

    # Ensure boundaries
    if not ss_outer_z or ss_outer_z[0] > ss_start_z + 1:
        inner_start_r = neckradius - get_thickness_at_distance(tubeheight - ss_start_height)
        ss_outer_z.insert(0, ss_start_z)
        ss_outer_r.insert(0, neckradius - PROTECTION_GAP_LAYER)
        ss_inner_z.insert(0, ss_start_z)
        ss_inner_r.insert(0, inner_start_r + PROTECTION_GAP_LAYER)

    if ss_outer_z and ss_outer_z[-1] < ss_end_z - 1:
        inner_end_r = neckradius - get_thickness_at_distance(0)
        ss_outer_z.append(ss_end_z)
        ss_outer_r.append(neckradius - PROTECTION_GAP_LAYER)
        ss_inner_z.append(ss_end_z)
        ss_inner_r.append(inner_end_r + PROTECTION_GAP_LAYER)

    # Add r=0 closure caps at BOTH ends with protection gap
    ss_outer_z.insert(0, ss_start_z + PROTECTION_GAP_LAYER)
    ss_outer_r.insert(0, 0)
    ss_inner_z.insert(0, ss_start_z)
    ss_inner_r.insert(0, 0)

    ss_outer_z.append(ss_end_z - PROTECTION_GAP_LAYER)
    ss_outer_r.append(0)
    ss_inner_z.append(ss_end_z)
    ss_inner_r.append(0)

    return ss_outer_z, ss_outer_r, ss_inner_z, ss_inner_r
