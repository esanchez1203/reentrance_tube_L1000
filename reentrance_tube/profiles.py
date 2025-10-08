"""Profile generation functions for reentrance tube geometry."""
from __future__ import annotations

# Constants
WLSR_TPB_THICKNESS = 1 * 1e-3  # 1 um TPB coating (in mm)
WLSR_TTX_THICKNESS = 254 * 1e-3  # 254 um Tetratex foil (in mm)
WLSR_THICKNESS = WLSR_TPB_THICKNESS + WLSR_TTX_THICKNESS


def get_thickness_at_distance(
    distance_from_top: float,
    include_wlsr: bool = False,
    wlsr_height: float = 2179,
    total_tube_height: float = 6247,
) -> float:
    """
    Calculate steel thickness at given distance from top of tube.
    
    Parameters
    ----------
    distance_from_top : float
        Distance from top of tube in mm
    include_wlsr : bool
        Whether to include WLSR thickness
    wlsr_height : float
        Height of WLSR region from bottom
    total_tube_height : float
        Total height of tube
        
    Returns
    -------
    float
        Steel thickness at this position in mm
    """
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

    # Add WLSR thickness if requested
    if include_wlsr:
        distance_from_bottom = total_tube_height - distance_from_top
        if distance_from_bottom <= wlsr_height:
            wlsr_total_thickness = WLSR_TTX_THICKNESS + WLSR_TPB_THICKNESS
            return steel_thickness + 2 * wlsr_total_thickness

    return steel_thickness


def make_outer_profile(
    neckradius: float,
    tubeheight: float,
    totalheight: float,
    curvefraction: float,
) -> tuple[list, list]:
    """
    Create complete outer profile with filled cylindrical section.
    
    Parameters
    ----------
    neckradius : float
        Radius of tube neck in mm
    tubeheight : float
        Total height of tube in mm
    totalheight : float
        Total height including top clearance in mm
    curvefraction : float
        Fraction of tube height that is curved at bottom
        
    Returns
    -------
    tuple[list, list]
        z and r coordinates for GenericPolycone
    """
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
            return z_filled, r_filled

    return z_out, r_out


def make_inner_profile(
    neckradius: float,
    tubeheight: float,
    totalheight: float,
    curvefraction: float,
) -> tuple[list, list]:
    """
    Create inner profile with variable thickness.
    
    Parameters
    ----------
    neckradius : float
        Radius of tube neck in mm
    tubeheight : float
        Total height of tube in mm
    totalheight : float
        Total height including top clearance in mm
    curvefraction : float
        Fraction of tube height that is curved at bottom
        
    Returns
    -------
    tuple[list, list]
        z and r coordinates for GenericPolycone
    """
    outer_z, outer_r = make_outer_profile(neckradius, tubeheight, totalheight, curvefraction)
    top_z_original = totalheight - 1

    inner_z, inner_r = [], []

    for z, r in zip(outer_z, outer_r):
        z_original = z + 5000
        dist_from_top = top_z_original - z_original
        thickness = get_thickness_at_distance(dist_from_top)

        if r == 0:
            inner_r_value = 0
        elif r == neckradius:
            inner_r_value = max(0, r - thickness)
        else:
            radius_ratio = r / neckradius
            scaled_thickness = thickness * radius_ratio
            inner_r_value = max(0, r - scaled_thickness)

        inner_z.append(z)
        inner_r.append(inner_r_value)

    return inner_z, inner_r


def make_wls_outer_profiles(
    neckradius: float,
    tubeheight: float,
    totalheight: float,
    curvefraction: float,
    wls_height: float,
    outer_z: list,
    outer_r: list,
) -> tuple:
    """
    Create outer WLS layer profiles (from steel outer surface inward).
    
    Returns properly closed polycone profiles with flat top cap.
    
    Returns
    -------
    tuple
        (wls_outer_z, wls_outer_r, wls_inner_z, wls_inner_r)
    """
    bottom_z = (totalheight - tubeheight) - 5000
    top_wls_z = bottom_z + wls_height

    wls_outer_z, wls_outer_r = [], []
    wls_inner_z, wls_inner_r = [], []

    # Extract points in WLS region
    for z, r in zip(outer_z, outer_r):
        if bottom_z <= z <= top_wls_z:
            wls_outer_z.append(z)
            wls_outer_r.append(r)

            if r == 0:
                inner_r_value = 0
            elif r == neckradius:
                inner_r_value = max(0, r - WLSR_THICKNESS)
            else:
                radius_ratio = r / neckradius
                scaled_thickness = WLSR_THICKNESS * radius_ratio
                inner_r_value = max(0, r - scaled_thickness)

            wls_inner_z.append(z)
            wls_inner_r.append(inner_r_value)

    # Ensure point at exact top height
    if wls_outer_z and wls_outer_z[-1] < top_wls_z - 1:
        wls_outer_z.append(top_wls_z)
        wls_outer_r.append(neckradius)
        wls_inner_z.append(top_wls_z)
        wls_inner_r.append(max(0, neckradius - WLSR_THICKNESS))

    # Add r=0 point at top to create flat end cap
    wls_outer_z.append(top_wls_z)
    wls_outer_r.append(0)
    wls_inner_z.append(top_wls_z)
    wls_inner_r.append(0)

    return wls_outer_z, wls_outer_r, wls_inner_z, wls_inner_r


def make_tetratex_outer_profiles(
    neckradius: float,
    tubeheight: float,
    totalheight: float,
    curvefraction: float,
    wls_height: float,
    outer_z: list,
    outer_r: list,
) -> tuple:
    """
    Create Tetratex layer profiles (inner layer of WLS, on steel surface).
    
    Returns properly closed polycone profiles with flat top cap.
    
    Returns
    -------
    tuple
        (ttx_outer_z, ttx_outer_r, ttx_inner_z, ttx_inner_r)
    """
    # Get WLS profiles (which span the full thickness)
    wls_outer_z, wls_outer_r, wls_inner_z, wls_inner_r = make_wls_outer_profiles(
        neckradius, tubeheight, totalheight, curvefraction, wls_height, outer_z, outer_r
    )

    # Tetratex inner = WLS inner (steel outer surface)
    ttx_inner_z = wls_inner_z.copy()
    ttx_inner_r = wls_inner_r.copy()

    # Tetratex outer = WLS inner + ttx_thickness
    ttx_outer_z = []
    ttx_outer_r = []

    for z, r in zip(wls_inner_z, wls_inner_r):
        ttx_outer_z.append(z)
        if r == 0:
            ttx_outer_r.append(0)
        elif r == neckradius:
            ttx_outer_r.append(r + WLSR_TTX_THICKNESS)
        else:
            radius_ratio = r / neckradius
            scaled_thickness = WLSR_TTX_THICKNESS * radius_ratio
            ttx_outer_r.append(r + scaled_thickness)

    return ttx_outer_z, ttx_outer_r, ttx_inner_z, ttx_inner_r


def make_wls_inner_profiles(
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
    Create inner WLS layer profiles (growing into steel from inner surface).
    
    Returns
    -------
    tuple
        (wls_outer_z, wls_outer_r, wls_inner_z, wls_inner_r)
    """
    bottom_z = (totalheight - tubeheight) - 5000
    top_wls_z = bottom_z + wls_height

    wls_inner_z, wls_inner_r = [], []
    wls_outer_z, wls_outer_r = [], []

    for z, r in zip(inner_z, inner_r):
        if bottom_z <= z <= top_wls_z:
            wls_inner_z.append(z)
            wls_inner_r.append(r)  # Steel inner surface
            wls_outer_z.append(z)

            if r == 0:
                wls_outer_r.append(0)
            else:
                # Find corresponding outer radius
                r_outer_steel = neckradius
                for z_steel, r_steel in zip(outer_z, outer_r):
                    if abs(z_steel - z) < 0.01:
                        r_outer_steel = r_steel
                        break

                # Grow INTO steel, clamped to not exceed outer surface
                r_target = r + WLSR_TPB_THICKNESS
                r_max = r_outer_steel - 0.05  # Small safety margin
                wls_outer_r.append(min(r_target, r_max))

    if wls_inner_z and wls_inner_z[-1] < top_wls_z - 1:
        last_r = wls_inner_r[-1]

        r_outer_steel = neckradius
        for z_steel, r_steel in zip(outer_z, outer_r):
            if abs(z_steel - top_wls_z) < 0.01:
                r_outer_steel = r_steel
                break

        wls_inner_z.append(top_wls_z)
        wls_inner_r.append(last_r)
        wls_outer_z.append(top_wls_z)
        r_target = last_r + WLSR_TPB_THICKNESS
        r_max = r_outer_steel - 0.05
        wls_outer_r.append(min(r_target, r_max))

    wls_inner_z.append(top_wls_z)
    wls_inner_r.append(0)
    wls_outer_z.append(top_wls_z)
    wls_outer_r.append(0)

    return wls_outer_z, wls_outer_r, wls_inner_z, wls_inner_r


def make_tetratex_inner_profiles(
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
    Create inner Tetratex layer profiles (outer layer, coating on TPB).
    
    Returns
    -------
    tuple
        (ttx_outer_z, ttx_outer_r, ttx_inner_z, ttx_inner_r)
    """
    # Get WLS profiles
    wls_outer_z, wls_outer_r, wls_inner_z, wls_inner_r = make_wls_inner_profiles(
        neckradius, tubeheight, totalheight, curvefraction, wls_height, inner_z, inner_r, outer_z, outer_r
    )

    # Tetratex inner = WLS outer
    ttx_inner_z = wls_outer_z.copy()
    ttx_inner_r = wls_outer_r.copy()

    # Tetratex outer = WLS outer + ttx_thickness (further into steel, CLAMPED)
    ttx_outer_z = []
    ttx_outer_r = []

    for z, r_wls in zip(wls_outer_z, wls_outer_r):
        ttx_outer_z.append(z)
        if r_wls == 0:
            ttx_outer_r.append(0)
        else:
            # Find corresponding outer radius
            r_outer_steel = neckradius
            for z_steel, r_steel in zip(outer_z, outer_r):
                if abs(z_steel - z) < 0.01:
                    r_outer_steel = r_steel
                    break

            # Grow further into steel with safety margin and clamping
            r_target = r_wls + WLSR_TTX_THICKNESS
            r_max = r_outer_steel - 0.05  # 0.05mm safety margin
            ttx_outer_r.append(min(r_target, r_max))

    return ttx_outer_z, ttx_outer_r, ttx_inner_z, ttx_inner_r


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
    """
    Create OFHC copper profiles as SOLID volume between steel inner and outer surfaces.
    
    Returns
    -------
    tuple
        (ofhc_outer_z, ofhc_outer_r, ofhc_inner_z, ofhc_inner_r)
    """
    bottom_z = (totalheight - tubeheight) - 5000
    ofhc_start_z = bottom_z + ofhc_start_height
    ofhc_end_z = bottom_z + ofhc_end_height

    ofhc_outer_z, ofhc_outer_r = [], []
    ofhc_inner_z, ofhc_inner_r = [], []

    # Process outer and inner together
    for i, (z_out, r_out) in enumerate(zip(outer_z, outer_r)):
        if ofhc_start_z <= z_out <= ofhc_end_z:
            ofhc_outer_z.append(z_out)
            ofhc_outer_r.append(r_out)
            ofhc_inner_z.append(inner_z[i])
            ofhc_inner_r.append(inner_r[i])

    # Ensure boundaries
    if ofhc_outer_z and ofhc_outer_z[0] > ofhc_start_z + 1:
        ofhc_outer_z.insert(0, ofhc_start_z)
        ofhc_outer_r.insert(0, neckradius)
        inner_start_r = neckradius - get_thickness_at_distance(tubeheight - ofhc_start_height)
        ofhc_inner_z.insert(0, ofhc_start_z)
        ofhc_inner_r.insert(0, inner_start_r)

    if ofhc_outer_z and ofhc_outer_z[-1] < ofhc_end_z - 1:
        ofhc_outer_z.append(ofhc_end_z)
        ofhc_outer_r.append(neckradius)
        inner_end_r = neckradius - get_thickness_at_distance(tubeheight - ofhc_end_height)
        ofhc_inner_z.append(ofhc_end_z)
        ofhc_inner_r.append(inner_end_r)

    # Add r=0 closure caps at both ends
    ofhc_outer_z.insert(0, ofhc_start_z)
    ofhc_outer_r.insert(0, 0)
    ofhc_inner_z.insert(0, ofhc_start_z)
    ofhc_inner_r.insert(0, 0)

    ofhc_outer_z.append(ofhc_end_z)
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
    """
    Create 316L stainless steel from ss_start_height to absolute top of tube.
    
    Returns
    -------
    tuple
        (ss_outer_z, ss_outer_r, ss_inner_z, ss_inner_r)
    """
    bottom_z = (totalheight - tubeheight) - 5000
    ss_start_z = bottom_z + ss_start_height
    ss_end_z = (totalheight - 1) - 5000  # Absolute top

    ss_outer_z, ss_outer_r = [], []
    ss_inner_z, ss_inner_r = [], []

    # Get all points in range, skip r=0 for now
    for i, (z_out, r_out) in enumerate(zip(outer_z, outer_r)):
        if ss_start_z <= z_out <= ss_end_z and r_out > 0:
            ss_outer_z.append(z_out)
            ss_outer_r.append(r_out)
            ss_inner_z.append(inner_z[i])
            ss_inner_r.append(inner_r[i])

    # Ensure boundaries
    if not ss_outer_z or ss_outer_z[0] > ss_start_z + 1:
        inner_start_r = neckradius - get_thickness_at_distance(tubeheight - ss_start_height)
        ss_outer_z.insert(0, ss_start_z)
        ss_outer_r.insert(0, neckradius)
        ss_inner_z.insert(0, ss_start_z)
        ss_inner_r.insert(0, inner_start_r)

    if ss_outer_z and ss_outer_z[-1] < ss_end_z - 1:
        inner_end_r = neckradius - get_thickness_at_distance(0)
        ss_outer_z.append(ss_end_z)
        ss_outer_r.append(neckradius)
        ss_inner_z.append(ss_end_z)
        ss_inner_r.append(inner_end_r)

    # Add r=0 closure caps at BOTH ends (bottom and top)
    ss_outer_z.insert(0, ss_start_z)
    ss_outer_r.insert(0, 0)
    ss_inner_z.insert(0, ss_start_z)
    ss_inner_r.insert(0, 0)

    ss_outer_z.append(ss_end_z)
    ss_outer_r.append(0)
    ss_inner_z.append(ss_end_z)
    ss_inner_r.append(0)

    print(f"DEBUG 316L: {len(ss_outer_z)} points")

    return ss_outer_z, ss_outer_r, ss_inner_z, ss_inner_r
