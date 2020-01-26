"""
"""

import numpy as numpy
import molecool

from .atom_data import atomic_weights

def test_molecular_mass():
    symbols = ['C', 'H', 'H', 'H', 'H']
    
    calculated_mass = molecool.calculate_molecular_mass(symbols)

    actual_mass = molecool.atom_data.atomic_weights['C'] + molecool.atom_data.atomic_weights['H'] +\
         molecool.atom_data.atomic_weights['H'] + molecool.atom_data.atomic_weights['H'] + molecool.atom_data.atom_weights['H']
    
    assert actual_mass == calculated_mass


def test_center_of_mass():
    symbols = np.array(['C', 'H', 'H', 'H', 'H'])
    coordinates = np.array([[1,1,1], [2.4,1,1], [-0.4, 1, 1], [1, 1, 2.4], [1, 1, -0.4]])

    center_of_mass = molecool.calculate_center_of_mass(symbols, coordinates)

    expected_center = np.array([1,1,1])

    assert np.allclose(expected_center, center_of_mass)   

def calculate_molecular_mass(symbols):
   """Calculate the mass of a molecule.
   
   Parameters
   ----------
   symbols : list
       A list of elements.
   
   Returns
   -------
   mass : float
       The mass of the molecule
   """
   pass

   mass = 0
   for atom in symbol:
       if atom not in atomic_weights.key():
           raise ValueError(F"Element")




    def calculate_center_of_mass(symbols, coordinates):
   """Calculate the center of mass of a molecule.
   
   The center of mass is weighted by each atom's weight.
   
   Parameters
   ----------
   symbols : list
       A list of elements for the molecule
   coordinates : np.ndarray
       The coordinates of the molecule.
   
   Returns
   -------
   center_of_mass: np.ndarray
       The center of mass of the molecule.
   Notes
   -----
   The center of mass is calculated with the formula
   
   .. math:: \\vec{R}=\\frac{1}{M} \\sum_{i=1}^{n} m_{i}\\vec{r_{}i}
   
   """
    total_mass = calculate_molecular_mass(symbols)
    mass_array = np.zeros([len(symbols), 1])

    for i in range(len(symbols)):
        mass_array[i] = atomic_weights[symbols[i]]

    center_of_mass = sum(coordinates * mass_array) / total_mass

    return center_of_mass




   return np.array([])