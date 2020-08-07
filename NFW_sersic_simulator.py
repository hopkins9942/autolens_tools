#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug  4 19:12:36 2020

@author: matthew

Simulates image of lens system with an NFW + Sersic lens, with parameters entered as command line arguments.

To add:
    save positions of lensed images
"""

# Configuration
from autoconf import conf
import os
import sys

workspace_path = "/home/matthew/Durham2020/AutoLensWorkspace/autolens_workspace/"

# Params - from command line arguments
sersic_elliptical_comps_0 = float(sys.argv[1])
sersic_elliptical_comps_1 = float(sys.argv[2])
NFW_mass = float(sys.argv[3])
einstein_radius = float(sys.argv[4])
source_centre_0 = float(sys.argv[5])
source_centre_1 = float(sys.argv[6])
source_elliptical_comps_0 = float(sys.argv[7])
source_elliptical_comps_1 = float(sys.argv[8])
source_radius = float(sys.argv[9])

sersic_elliptical_comps = (sersic_elliptical_comps_0, sersic_elliptical_comps_1)
source_centre = (source_centre_0, source_centre_1)
source_elliptical_comps = (source_elliptical_comps_0, source_elliptical_comps_1)


dataset_name = "_".join(
    (
        str(sersic_elliptical_comps),
        str(NFW_mass),
        str(einstein_radius),
        str(source_centre),
        str(source_elliptical_comps),
        str(source_radius),
    )
)

save_path = (
    "/home/matthew/Durham2020/MyScripts/autolens_tools/ExampleData/NFW_sersic/"
    + dataset_name
    + "/"
)

conf.instance = conf.Config(
    config_path=f"{workspace_path}config/", output_path=save_path,
)

if not os.path.exists(save_path):
    os.makedirs(save_path)


# Script
import autolens as al
import autolens.plot as aplt
import datetime
from scipy.optimize import fsolve

print("Generating:")
print(datetime.datetime.now())

redshift_lens = 0.5
redshift_source = 1.0

# Creating profiles and galaxies:
dark_profile = al.mp.SphericalTruncatedNFWMCRLudlow(
    centre=(0.0, 0.0),
    mass_at_200=10 ** NFW_mass,
    redshift_object=redshift_lens,
    redshift_source=redshift_source,
)


def sersic_profile(intensity):
    return al.mp.EllipticalSersic(
        centre=(0.0, 0.0),
        elliptical_comps=sersic_elliptical_comps,
        effective_radius=1.2,
        sersic_index=4,
        intensity=intensity,
    )


def intensity_setter(target_e_radius):
    """
    For given profiles, calculates intensity that gives correct Einstein radius
    """
    print("\nsetting intensity:")
    print(datetime.datetime.now())

    def func(intensity):
        value = (
            al.Galaxy(
                0.5, sersic=sersic_profile(intensity), dark=dark_profile
            ).einstein_radius_in_units()
            - target_e_radius
        )
        print(f"Intensity: {intensity}")
        print(f"Value: {value}")
        return value

    return fsolve(func, 0.1, xtol=1e-3)


lens_galaxy = al.Galaxy(
    redshift=redshift_lens,
    sersic=sersic_profile(intensity_setter(einstein_radius)),
    dark=dark_profile,
)

source_galaxy = al.Galaxy(
    redshift=redshift_source,
    light=al.lp.EllipticalExponential(
        centre=source_centre,
        elliptical_comps=source_elliptical_comps,
        intensity=0.15,
        effective_radius=source_radius,
    ),
)

# Creating Tracer and Simualtor, then writing data
tracer = al.Tracer.from_galaxies(galaxies=[lens_galaxy, source_galaxy])
tracer.save(file_path=save_path, filename="true_tracer")
grid = al.Grid.uniform(shape_2d=(200, 200), pixel_scales=0.05,)


print("\nplotting observed image")
print(datetime.datetime.now())
psf = al.Kernel.from_gaussian(
    shape_2d=(11, 11), sigma=0.1, pixel_scales=grid.pixel_scales
)
simulator = al.SimulatorImaging(
    exposure_time_map=al.Array.full(fill_value=7500.0, shape_2d=grid.shape_2d),
    psf=psf,
    background_sky_map=al.Array.full(fill_value=0.1, shape_2d=grid.shape_2d),
    add_noise=True,
)
imaging = simulator.from_tracer_and_grid(tracer=tracer, grid=grid, name=dataset_name)
imaging.output_to_fits(
    image_path=f"{save_path}image.fits",
    psf_path=f"{save_path}psf.fits",
    noise_map_path=f"{save_path}noise_map.fits",
    overwrite=True,
)

observed_plotter = aplt.Plotter(
    output=aplt.Output(path=save_path, filename="ObservedImage", format="png")
)
aplt.Imaging.image(imaging, plotter=observed_plotter)
