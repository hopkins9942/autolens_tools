#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug  4 19:14:18 2020

@author: matthew

Fits Lens image with power law profile using two phases

To Do:
    test
    utilise positions
"""



from autoconf import conf
import os
import sys

#Configuration
data_path = sys.argv[1]
workspace_path = "/home/matthew/Durham2020/AutoLensWorkspace/autolens_workspace/" 


dataset_name = sys.argv[2]
output_path = os.path.join(data_path, dataset_name) + '/'
mask_radius = sys.argv[3]
    
conf.instance = conf.Config(
    config_path=f"{workspace_path}config/", output_path=output_path,
)


#Script
import autolens as al
import autolens.plot as aplt
import datetime
import autofit as af

imaging = al.Imaging.from_fits(
    image_path=f"{output_path}image.fits",
    noise_map_path=f"{output_path}noise_map.fits",
    psf_path=f"{output_path}psf.fits",
    pixel_scales=0.05,
    name = dataset_name
) 

pickle_files = [f"{output_path}true_tracer.pickle"]

info = {"dataset_name": dataset_name, "pixel_scales": imaging.pixel_scales}



#Fitting - using imaging object from fits
print(f"\nFitting {dataset_name} at {datetime.datetime.now()}")

settings = al.PhaseSettingsImaging(grid_class=al.Grid, sub_size=2)

# First phase
print(f"\n{dataset_name} creating phase 1 at {datetime.datetime.now()}")

lens_model_1 = al.GalaxyModel(
    redshift=0.5, mass=al.mp.EllipticalPowerLaw
)
source_model_1 = al.GalaxyModel(redshift=1.0, light=al.lp.EllipticalExponential)

mask = al.Mask.circular(
    shape_2d=imaging.shape_2d, pixel_scales=imaging.pixel_scales, radius=mask_radius
    )


# simplifying assumptions - check against plots from generator
lens_model_1.mass.centre_0 = 0.0
lens_model_1.mass.centre_1 = 0.0
lens_model_1.mass.slope = 2 #ie isothermal


phase_1 = al.PhaseImaging(
    phase_name="phase_1",
    settings=settings,
    galaxies=dict(lens_model_1=lens_model_1, source_model_1=source_model_1),
    search=af.DynestyStatic(n_live_points=50, evidence_tolerance=5.0),
)

print(f"\n{dataset_name} running phase 1 at {datetime.datetime.now()}")

phase_1_result = phase_1.run(dataset=imaging, mask=mask)

phase_1_sub_plotter = aplt.SubPlotter(output=aplt.Output(path = output_path, filename="Phase1", format = "png"))
aplt.FitImaging.subplot_fit_imaging(fit=phase_1_result.max_log_likelihood_fit, sub_plotter=phase_1_sub_plotter)


# Second phase
print(f"\n{dataset_name} creating phase 2 at {datetime.datetime.now()}")


lens_model_2 = al.GalaxyModel(redshift=0.5,
                              mass=af.PriorModel(al.mp.EllipticalPowerLaw),
                              )
lens_model_2.mass.elliptical_comps = phase_1_result.model.galaxies.lens_model_1.mass.elliptical_comps
lens_model_2.mass.einstein_radius = phase_1_result.model.galaxies.lens_model_1.mass.einstein_radius


phase_2 = al.PhaseImaging(
    phase_name="phase_2",
    settings=settings,
    galaxies=dict(lens_model_2=lens_model_2, 
                  source_model_2=phase_1_result.model.galaxies.source_model_1),
    search=af.DynestyStatic(n_live_points=50, evidence_tolerance=5.0),
)
print(f"\n{dataset_name} running phase 2 at {datetime.datetime.now()}")

phase_2_result=phase_2.run(dataset=imaging, mask=mask)

phase_2_sub_plotter = aplt.SubPlotter(output=aplt.Output(path = output_path, filename="Phase2", format = "png"))
aplt.FitImaging.subplot_fit_imaging(fit=phase_2_result.max_log_likelihood_fit, sub_plotter=phase_2_sub_plotter)

print(f"\n{dataset_name} done at {datetime.datetime.now()}")

