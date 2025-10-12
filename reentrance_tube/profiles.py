"""Profile generation functions for reentrance tube geometry."""
from __future__ import annotations

# Constants
WLSR_TPB_THICKNESS = 1 * 1e-3  # 1 um TPB coating (in mm)
WLSR_TTX_THICKNESS = 254 * 1e-3  # 254 um Tetratex foil (in mm)
WLSR_THICKNESS = WLSR_TPB_THICKNESS + WLSR_TTX_THICKNESS
PROTECTION_GAP = 5 * 1e-6  # 5 nm gap in mm


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
    #if include_wlsr:
    #    distance_from_bottom = total_tube_height - distance_from_top
    #    if distance_from_bottom <= wlsr_height:
    #        wlsr_total_thickness = WLSR_TTX_THICKNESS + WLSR_TPB_THICKNESS
    #        return steel_thickness + 2 * wlsr_total_thickness

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
    Create inner WLS layer profiles for placement in underground argon.
    
    Structure: Steel -> Gap -> TTX -> TPB -> UAr
    
    TTX outer surface = Tube inner - gap
    TTX inner surface = Tube inner - gap - TTX_thickness - TPB_thickness
    TPB inner surface = Tube inner - gap - TTX_thickness
    TPB outer surface = Tube inner - gap - TTX_thickness - TPB_thickness
    
    Returns
    -------
    tuple
        (tpb_outer_z, tpb_outer_r, tpb_inner_z, tpb_inner_r,
         ttx_outer_z, ttx_outer_r, ttx_inner_z, ttx_inner_r)
    """
    bottom_z = (totalheight - tubeheight) - 5000
    top_wls_z = bottom_z + wls_height

    ttx_outer_z, ttx_outer_r = [], []
    ttx_inner_z, ttx_inner_r = [], []
    tpb_outer_z, tpb_outer_r = [], []
    tpb_inner_z, tpb_inner_r = [], []

    for z, r_inner in zip(inner_z, inner_r):
        if bottom_z <= z <= top_wls_z:
            # Skip r=0 points except at the very top
            if r_inner == 0 and z < top_wls_z - 0.01:
                continue
                
            if r_inner == 0:
                ttx_outer_z.append(z)
                ttx_outer_r.append(0)
                ttx_inner_z.append(z)
                ttx_inner_r.append(0)
                tpb_outer_z.append(z)
                tpb_outer_r.append(0)
                tpb_inner_z.append(z)
                tpb_inner_r.append(0)
            else:
                r_ttx_outer_desired = r_inner - PROTECTION_GAP 
                r_ttx_outer = min(r_ttx_outer_desired, r_inner)
                r_ttx_outer = max(0, r_ttx_outer)
                
                r_ttx_inner = max(0, r_inner - PROTECTION_GAP  - WLSR_TTX_THICKNESS - WLSR_TPB_THICKNESS)
                r_tpb_outer = max(0, r_inner - PROTECTION_GAP  - WLSR_TTX_THICKNESS)
                r_tpb_inner = max(0, r_inner - PROTECTION_GAP  - WLSR_TTX_THICKNESS - WLSR_TPB_THICKNESS)
                
                ttx_outer_z.append(z)
                ttx_outer_r.append(r_ttx_outer)
                ttx_inner_z.append(z)
                ttx_inner_r.append(r_ttx_inner)
                tpb_outer_z.append(z)
                tpb_outer_r.append(r_tpb_outer)
                tpb_inner_z.append(z)
                tpb_inner_r.append(r_tpb_inner)

    # Ensure point at exact top height
    if ttx_outer_z and ttx_outer_z[-1] < top_wls_z - 1:
        last_r_inner = inner_r[-1]
        for i, z in enumerate(inner_z):
            if abs(z - top_wls_z) < 0.1:
                last_r_inner = inner_r[i]
                break
        
        r_ttx_outer_desired = last_r_inner - PROTECTION_GAP 
        r_ttx_outer = min(r_ttx_outer_desired, last_r_inner )
        r_ttx_outer = max(0, r_ttx_outer)
        r_ttx_inner = max(0, last_r_inner - PROTECTION_GAP  - WLSR_TTX_THICKNESS - WLSR_TPB_THICKNESS)
        r_tpb_outer = max(0, last_r_inner - PROTECTION_GAP  - WLSR_TTX_THICKNESS)
        r_tpb_inner = max(0, last_r_inner - PROTECTION_GAP  - WLSR_TTX_THICKNESS - WLSR_TPB_THICKNESS)
        
        ttx_outer_z.append(top_wls_z)
        ttx_outer_r.append(r_ttx_outer)
        ttx_inner_z.append(top_wls_z)
        ttx_inner_r.append(r_ttx_inner)
        tpb_outer_z.append(top_wls_z)
        tpb_outer_r.append(r_tpb_outer)
        tpb_inner_z.append(top_wls_z)
        tpb_inner_r.append(r_tpb_inner)

    # Add r=0 closure at top
    ttx_outer_z.append(top_wls_z)
    ttx_outer_r.append(0)
    ttx_inner_z.append(top_wls_z)
    ttx_inner_r.append(0)
    tpb_outer_z.append(top_wls_z)
    tpb_outer_r.append(0)
    tpb_inner_z.append(top_wls_z)
    tpb_inner_r.append(0)

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
    - TTX mother volume.
    - TPB is a 1 µm daughter coating on its outer surface.
    """

    global WLSR_TTX_THICKNESS, WLSR_TPB_THICKNESS

    bottom_z = (totalheight - tubeheight) - 5000
    top_wls_z = bottom_z + wls_height

    ttx_outer_z, ttx_outer_r = [], []
    ttx_inner_z, ttx_inner_r = [], []
    tpb_outer_z, tpb_outer_r = [], []
    tpb_inner_z, tpb_inner_r = [], []

    for z, r in zip(outer_z, outer_r):
        if bottom_z <= z <= top_wls_z and r > 0:
            base = r + PROTECTION_GAP  # tube outer + gap

            # TTX (mother)
            ttx_inner = base
            ttx_outer = base + WLSR_TTX_THICKNESS + WLSR_TPB_THICKNESS

            # TPB (daughter)
            tpb_inner = base + WLSR_TTX_THICKNESS
            tpb_outer = tpb_inner + WLSR_TPB_THICKNESS

            ttx_inner_z.append(z); ttx_inner_r.append(ttx_inner)
            ttx_outer_z.append(z); ttx_outer_r.append(ttx_outer)
            tpb_inner_z.append(z); tpb_inner_r.append(tpb_inner)
            tpb_outer_z.append(z); tpb_outer_r.append(tpb_outer)

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
