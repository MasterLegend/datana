# -*- coding: utf-8 -*-

import pyvista as pv
from pyvistaqt import BackgroundPlotter
from structure import Structure

def plot_structure(s):
    
    plotter = BackgroundPlotter()
    
    for i in range(len(s._coords_cartesian)):
         _ = plotter.add_mesh(pv.Sphere(radius=s.get_atom(i).get_radius(), center = s._coords_cartesian[i]))
    
    _ = plotter.add_mesh(pv.Line(pointa=(0.0, 0.0, 0.0), pointb=s._lattice[0], resolution=1), color='white')
    _ = plotter.add_mesh(pv.Line(pointa=(0.0, 0.0, 0.0), pointb=s._lattice[1], resolution=1), color='white')
    _ = plotter.add_mesh(pv.Line(pointa=s._lattice[0], pointb=s._lattice[0]+s._lattice[1], resolution=1), color='white')
    _ = plotter.add_mesh(pv.Line(pointa=s._lattice[1], pointb=s._lattice[0]+s._lattice[1], resolution=1), color='white')
    _ = plotter.add_mesh(pv.Line(pointa=(0.0, 0.0, 0.0), pointb=s._lattice[2], resolution=1), color='white')
    _ = plotter.add_mesh(pv.Line(pointa=s._lattice[0], pointb=s._lattice[0]+s._lattice[2], resolution=1), color='white')
    _ = plotter.add_mesh(pv.Line(pointa=s._lattice[1], pointb=s._lattice[1]+s._lattice[2], resolution=1), color='white')
    _ = plotter.add_mesh(pv.Line(pointa=s._lattice[0]+s._lattice[1], pointb=s._lattice[0]+s._lattice[1]+s._lattice[2], resolution=1), color='white')
    _ = plotter.add_mesh(pv.Line(pointa=s._lattice[2], pointb=s._lattice[0]+s._lattice[2], resolution=1), color='white')
    _ = plotter.add_mesh(pv.Line(pointa=s._lattice[2], pointb=s._lattice[1]+s._lattice[2], resolution=1), color='white')
    _ = plotter.add_mesh(pv.Line(pointa=s._lattice[0]+s._lattice[2], pointb=s._lattice[0]+s._lattice[1]+s._lattice[2], resolution=1), color='white')
    _ = plotter.add_mesh(pv.Line(pointa=s._lattice[1]+s._lattice[2], pointb=s._lattice[0]+s._lattice[1]+s._lattice[2], resolution=1), color='white')
    
    return