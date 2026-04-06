import re
import numpy as np
import math


def ensure_unique_name(name: str,
                       existing: set[str],
                       suffix_str: str = '_',
                       max_length: int | None = None) -> str:
    suffix_pattern = rf"{suffix_str}$"
    if name not in existing:
        existing.add(name)
        return name
    else:
        base = name
        m = re.search(suffix_pattern, base)
        start = int(m.group(1)) + 1 if m else 2
        i = start
        while True:
            candidate = f"{base}{suffix_str}{i}"
            if max_length:
                if len(candidate) > max_length:
                    overflow = len(candidate) - max_length
                    candidate = f"{base[:max_length]}{suffix_str}{i}"
            if candidate not in existing:
                existing.add(candidate)
                return candidate
            i += 1


def is_valid_regex_pattern(pattern: str) -> str:
    try:
        re.compile(pattern)
        return True
    except re.error:
        return False


def sanitize_name(name: str,
                  invalid_chars: str | None = None,
                  invalid_replace_char: str = '_',
                  max_length: int | None = None) -> str:
    """
    Santize Excel sheet name to satisfy excel constraints.
    """
    # Basic validation
    if not (isinstance(name, str)):
        raise ValueError(
            f"name is not a string-type value"
        )

    if not (isinstance(invalid_replace_char, str)):
        raise ValueError(
            f"invalid_replace_char is not a string-type value"
        )

    if not isinstance(invalid_chars, (str, type(None))):
        raise ValueError(
            f"Type of name is not a string-type value"
        )
    else:
        if not isinstance(invalid_chars, type(None)):
            if not (is_valid_regex_pattern(invalid_chars)):
                return None
        else:
            pass
    if not isinstance(max_length, (int, type(None))):
        raise ValueError(
            f"max_legnth is not a integer value"
        )
    # Removing comas
    name = name.strip().strip("'").strip()

    # Replacing invalid chars
    if invalid_chars:
        name = re.sub(invalid_chars, invalid_replace_char, name)

    # Limiting length
    if max_length:
        if len(name) > max_length:
            name = name[:max_length]

    return name


def transform_stress_global_to_local(stress_global: np.ndarray, theta_deg: float) -> np.ndarray:
    '''
    Performs a rotation transformation of the components of a plane-stress vector
    {sigma_x, sigma_y, sigma_xy} from a global cordinate system into a local.

    Parameters
    ----------
    stress_global : np.ndarray

    theta_deg : float
        Angle of rotation about the 3-axis

    Returns
    -------
    stress_local : np.ndarray
        DESCRIPTION.

    '''

    theta = math.radians(theta_deg)
    c = math.cos(theta)
    s = math.sin(theta)
    c2 = c**2
    s2 = s**2
    cs = c*s
    sx, sy, txy = stress_global
    s1 = c2 * sx + s2 * sy + 2 * cs * txy
    s2 = s2 * sx + c2 * sy - 2 * cs * txy
    t12 = -cs * sx + cs * sy + (c2 - s2) * txy
    stress_local = np.array([s1, s2, t12])
    return stress_local


def transform_strain_global_to_local(strain_global: np.ndarray, theta_deg: float) -> np.ndarray:
    theta = math.radians(theta_deg)
    c = math.cos(theta)
    s = math.sin(theta)
    c2 = c**2
    s2 = s**2
    cs = c*s
    ex, ey, gxy = strain_global
    e1 = c2 * ex + s2 * ey + cs * gxy
    e2 = s2 * ex + c2 * ey - cs * gxy
    g12 = -2.0 * cs * ex + 2.0 * cs * ey + (c2 - s2) * gxy
    strain_local = np.array([e1, e2, g12])
    return strain_local
