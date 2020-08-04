#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug  4 19:13:59 2020

@author: matthew
"""

import autolens as al
import autolens.plot as aplt
import os
import datetime
import sys

#Params
NFW_mass=12
sersic1_elliptical_comps=(0.0,0.111)
sersic1_intensity=0.1
sersic1_effective_radius = 1.0
sersic2_elliptical_comps=(0.0, - 0.333)
sersic2_intensity=0.3
sersic2_effective_radius = 1.0
source_centre=(0.1,0.2)
source_elliptical_comps=(0.1,0.2)
source_radius=0.1

# sersic_elliptical_comps_0=float(sys.argv[1])
# sersic_elliptical_comps_1=float(sys.argv[2])
# NFW_mass=float(sys.argv[3])
# einstein_radius=float(sys.argv[4])
# source_centre_0=float(sys.argv[5])
# source_centre_1=float(sys.argv[6])
# source_elliptical_comps_0=float(sys.argv[7])
# source_elliptical_comps_1=float(sys.argv[8])
# source_radius=float(sys.argv[9])

# sersic_elliptical_comps = (sersic_elliptical_comps_0, sersic_elliptical_comps_1)
# source_centre = (source_centre_0, source_centre_1)
# source_elliptical_comps = (source_elliptical_comps_0, source_elliptical_comps_1)


print("Generating:")
print(datetime.datetime.now())

save_path= "_".join(("/home/matthew/Durham2020/Outputs/Run3/",
                     str(NFW_mass),
                     str(sersic1_elliptical_comps),
                     str(sersic1_intensity),
                     str(sersic1_effective_radius),
                     str(sersic2_elliptical_comps),
                     str(sersic2_intensity),
                     str(sersic2_effective_radius),
                     str(source_centre),
                     str(source_elliptical_comps),
                     str(source_radius))) + "/"
#save_path="/home/matthew/Durham2020/"
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

def structure_profile1(mass_to_light_ratio):
    return al.mp.EllipticalSersic(
        centre = (0.0, 0.0),
        elliptical_comps = sersic1_elliptical_comps,
        effective_radius = sersic1_effective_radius,
        sersic_index = 4,
        intensity = sersic1_intensity,
        )
def structure_profile2(mass_to_light_ratio): # for now use same mass light ratio, just pick an effective radius
    return al.mp.EllipticalSersic(
        centre = (0.0, 0.0),
        elliptical_comps = sersic2_elliptical_comps,
        effective_radius = sersic2_effective_radius,
        sersic_index = 1,
        intensity = sersic2_intensity,
        )

# def mlr_setter(target):# for given profiles, calculates mass to light ratio that gives correct Einstein radius
#     print("\nsetting mlr:")
#     print(datetime.datetime.now())
#     mlr = 1.6 #best value for target=1
#     difference = 0.1
#     tolerance = 1e-2
#     e_radius_list = [
#         al.Galaxy(0.5, light=structure_profile(mlr - difference), dark=dark_profile).einstein_radius_in_units(),
#         al.Galaxy(0.5, light=structure_profile(mlr             ), dark=dark_profile).einstein_radius_in_units(),
#         al.Galaxy(0.5, light=structure_profile(mlr + difference), dark=dark_profile).einstein_radius_in_units(),
#         ]
#     for i in range(100):
#         print(f"\ni = {i}, mlr = {mlr}, difference = {difference},\ne_radius_list=\n{e_radius_list}")
#         print(datetime.datetime.now())
#         if (e_radius_list[0]<target) and (target<e_radius_list[2]):
#             if abs(e_radius_list[1]-e_radius_list[0])>tolerance or abs(e_radius_list[1]-e_radius_list[2])>tolerance:
#                 difference *= 0.1
#                 e_radius_list[0] = al.Galaxy(0.5, light=structure_profile(mlr - difference), dark=dark_profile).einstein_radius_in_units()
#                 e_radius_list[2] = al.Galaxy(0.5, light=structure_profile(mlr + difference), dark=dark_profile).einstein_radius_in_units()
#             else:
#                 break
#         elif (e_radius_list[1]<target):
#             mlr += difference
#             e_radius_list[0:2] = e_radius_list[1:3]
#             e_radius_list[2] = al.Galaxy(0.5, light=structure_profile(mlr + difference), dark=dark_profile).einstein_radius_in_units()
#         else:
#             mlr -= difference
#             e_radius_list[1:3] = e_radius_list[0:2]
#             e_radius_list[0] = al.Galaxy(0.5, light=structure_profile(mlr - difference), dark=dark_profile).einstein_radius_in_units()
#     print(f"\nfinal mlr: {mlr}, final e_radius_list: \n{e_radius_list}")
#     print(datetime.datetime.now())
#     return mlr

lens_galaxy = al.Galaxy(
    redshift=redshift_lens,
    structure1=structure_profile1(1.0), #structure_profile(mlr_setter(einstein_radius)),
    structure2=structure_profile2(1.0),
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
tracer.save(file_path=save_path, filename="true_tracer") # save tracer
grid = al.Grid.uniform(
    shape_2d=(100, 100), 
    pixel_scales=0.1,
)

# print("\nplotting true image")
# print(datetime.datetime.now())
# true_plotter = aplt.Plotter(output=aplt.Output(path = save_path, filename="TrueImage", format = "png"))
#aplt.Tracer.image(tracer,grid=grid,plotter=true_plotter) takes a while, might as well just plot realistic while producing FITS

print("\nplotting observed image")
print(datetime.datetime.now())
psf = al.Kernel.from_gaussian(
    shape_2d=(11, 11), sigma=0.1, pixel_scales=grid.pixel_scales
)
simulator = al.SimulatorImaging(
    exposure_time_map=al.Array.full(fill_value=300.0, shape_2d=grid.shape_2d),
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
# do this later
# positions_list=tracer.image_plane_multiple_image_positions_of_galaxies(grid=grid)
# # with open(f"{save_path}positions","w") as f:
# #     f.write(str(positions_list))
# # NOTE: for next time, use GridCoordinates.output_to_file()
# al.GridCoordinates.output_to_file()

# observed_plotter = aplt.Plotter(output=aplt.Output(path = save_path, filename="ObservedImage", format = "png"))
# if len(positions_list[0])!=0:
#     positions = al.GridCoordinates(positions_list)
#     aplt.Imaging.image(imaging,plotter=observed_plotter,positions=positions)
# else:
#     aplt.Imaging.image(imaging,plotter=observed_plotter)


