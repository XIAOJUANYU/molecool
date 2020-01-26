"""
measure.py

this module is for functions for mesurement
"""
import numpy as np


def calculate_distance(rA, rB):
    """
    Calculate the distance between two points

    Parameters
    ----------
    rA, rB: np.adarray
        The coordinates of each points.

    Returns
    -------
    distance : float        
        The distance between the two points.

    Examples
    --------
    >>> r1 = np.array([0, 0, 1])
    >>> r2 = np.array({0, 0.1, 0})
    >>> calculate_distance(r1, r2)
    0.1
    """

    if not isinstance(rA, np.ndarray) or not isinstance(rB, np.ndarray):
        raise TypeError("Input must be type np.ndarray for calculate_distance")


    # This function calculates the distance between two points given as numpy arrays.
    distance_vector = (rA - rB)
    distance_scalar = np.linalg.norm(distance_vector)
    return distance_scalar

def calculate_angle(rA, rB, rC, degrees=False):
    # Calculate the angle between three points. Answer is given in radians by default, but can be given in degrees
    # by setting degrees=True
    AB = rB - rA
    BC = rB - rC
    theta=np.arccos(np.dot(AB, BC)/(np.linalg.norm(AB)*np.linalg.norm(BC)))

    if degrees:
        return np.degrees(theta)
    else:
        return theta    

