#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug  4 19:13:59 2020

@author: hopkins9942

Simulates image of lens system with an NFW + 2*Sersic lens, with parameters entered as command line arguments.

To add:
    Test set values such as effective radii of sersics
    tweak starting value in intensity setter
    save positions of lensed images
"""

import autolens as al
import autolens.plot as aplt
import os
import datetime
import sys

#Params
sersic1_elliptical_comps_0=float(sys.argv[1])
sersic1_elliptical_comps_1=float(sys.argv[2])
sersic2_elliptical_comps_0=float(sys.argv[3])
sersic2_elliptical_comps_1=float(sys.argv[4])
NFW_mass=float(sys.argv[5])
einstein_radius=float(sys.argv[6])
structure_mass_ratio=float(sys.argv[7])
# equivalent to structure instensity ratio as same mass/light ratio is assumed, equal to ratio of structure2 intensity to structure1 intensity
source_centre_0=float(sys.argv[8])
source_centre_1=float(sys.argv[9])
source_elliptical_comps_0=float(sys.argv[10])
source_elliptical_comps_1=float(sys.argv[11])
source_radius=float(sys.argv[12])

sersic1_elliptical_comps = (sersic1_elliptical_comps_0, sersic1_elliptical_comps_1)
sersic2_elliptical_comps = (sersic2_elliptical_comps_0, sersic2_elliptical_comps_1)
source_centre = (source_centre_0, source_centre_1)
source_elliptical_comps = (source_elliptical_comps_0, source_elliptical_comps_1)


print("Generating:")
print(datetime.datetime.now())


dataset_name = "_".join((
                      str(sersic1_elliptical_comps),
                      str(sersic2_elliptical_comps),
                      str(NFW_mass),
                      str(einstein_radius),
                      str(structure_mass_ratio),
                      str(source_centre),
                      str(source_elliptical_comps),
                      str(source_radius)))

save_path = "/home/matthew/Durham2020/MyScripts/autolens_tools/ExampleData/NFW_double_sersic/" + dataset_name + "/"

if not os.path.exists(save_path):
    os.makedirs(save_path)
    
redshift_lens = 0.5
redshift_source = 1.0

dark_profile=al.mp.SphericalTruncatedNFWMCRLudlow(
        centre = (0.0, 0.0),
        mass_at_200 = 10**NFW_mass,
        redshift_object = redshift_lens,
        redshift_source = redshift_source,
    )

def structure_profile1(intensity1):
    return al.mp.EllipticalSersic(
        centre = (0.0, 0.0),
        elliptical_comps = sersic1_elliptical_comps,
        effective_radius = 1.2,
        sersic_index = 4,
        intensity = intensity1,
        )
def structure_profile2(intensity1): 
    return al.mp.EllipticalSersic(
        centre = (0.0, 0.0),
        elliptical_comps = sersic2_elliptical_comps,
        effective_radius = 0.4,
        sersic_index = 1,
        intensity = structure_mass_ratio*intensity1,
        )

def intensity_setter(target):
    """
    for given profiles, calculates intensity of structure1 that gives correct Einstein radius
    """
    print("\nsetting intensity:")
    print(datetime.datetime.now())
    intensity = 0.1 
    difference = 0.01
    tolerance = 1e-2
    e_radius_list = [
        al.Galaxy(0.5, light1=structure_profile1(intensity - difference), light2=structure_profile2(intensity - difference), dark=dark_profile).einstein_radius_in_units(),
        al.Galaxy(0.5, light1=structure_profile1(intensity             ), light2=structure_profile2(intensity - difference), dark=dark_profile).einstein_radius_in_units(),
        al.Galaxy(0.5, light1=structure_profile1(intensity + difference), light2=structure_profile2(intensity - difference), dark=dark_profile).einstein_radius_in_units(),
        ]
    for i in range(100):
        print(f"\ni = {i}, intensity = {intensity}, difference = {difference},\ne_radius_list=\n{e_radius_list}")
        print(datetime.datetime.now())
        if (e_radius_list[0]<target) and (target<e_radius_list[2]):
            if abs(e_radius_list[1]-e_radius_list[0])>tolerance or abs(e_radius_list[1]-e_radius_list[2])>tolerance:
                difference *= 0.1
                e_radius_list[0] = al.Galaxy(0.5, light1=structure_profile1(intensity - difference), light2=structure_profile2(intensity - difference), dark=dark_profile).einstein_radius_in_units()
                e_radius_list[2] = al.Galaxy(0.5, light1=structure_profile1(intensity + difference), light2=structure_profile2(intensity + difference), dark=dark_profile).einstein_radius_in_units()
            else:
                break
        elif (e_radius_list[1]<target):
            intensity += difference
            e_radius_list[0:2] = e_radius_list[1:3]
            e_radius_list[2] = al.Galaxy(0.5, light1=structure_profile1(intensity + difference), light2=structure_profile2(intensity + difference), dark=dark_profile).einstein_radius_in_units()
        else:
            intensity -= difference
            e_radius_list[1:3] = e_radius_list[0:2]
            e_radius_list[0] = al.Galaxy(0.5, light1=structure_profile1(intensity - difference), light2=structure_profile2(intensity - difference), dark=dark_profile).einstein_radius_in_units()
    print(f"\nfinal intensity: {intensity}, final e_radius_list: \n{e_radius_list}")
    print(datetime.datetime.now())
    return intensity

intensity1 = intensity_setter(einstein_radius)

lens_galaxy = al.Galaxy(
    redshift=redshift_lens,
    structure1=structure_profile1(intensity1), #structure_profile(mlr_setter(einstein_radius)),
    structure2=structure_profile2(intensity1),
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

tracer = al.Tracer.from_galaxies(galaxies=[lens_galaxy, source_galaxy])
tracer.save(file_path=save_path, filename="true_tracer")
grid = al.Grid.uniform(
    shape_2d=(200, 200), 
    pixel_scales=0.05,
)

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
imaging = simulator.from_tracer_and_grid(tracer=tracer, grid=grid)
imaging.output_to_fits(
    image_path=f"{save_path}image.fits",
    psf_path=f"{save_path}psf.fits",
    noise_map_path=f"{save_path}noise_map.fits",
    overwrite=True,
)

observed_plotter = aplt.Plotter(output=aplt.Output(path = save_path, filename="ObservedImage", format = "png"))
aplt.Imaging.image(imaging,plotter=observed_plotter)


