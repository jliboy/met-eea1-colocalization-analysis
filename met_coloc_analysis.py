#This code is for estimating the fraction of protein MET that colocalizes with an endosome marker. 
#It first segments the nuclei, MET, and endosomal marker signal (EEA1). 
#Then it calculates the colocalization between MET and EEA1 signals.
#It also determines properties such as the number of MET foci, the total area of the MET signal, and the intensity of the MET signal per nucleus.
#It also determines the count of EEA1 foci. 

import glob
import tkinter as tk 
from tkinter import filedialog
import os

import sys
sys.path.append('src') #Calling needed functions

root = tk.Tk()
root.withdraw()

directory_path = filedialog.askdirectory()
images_path = glob.glob(os.path.join(directory_path, "HGF*")) 
#assuming images are named with HGF at the beginning, adjust as needed to match your file naming convention.

import nd2reader as nd2
import matplotlib.pyplot as plt
import numpy as np
import re

filename = []
overlap_signal = []
normalized_overlap_signal = []
nuclear_count = []

met_count = []
normalized_met_count = []

eea1_count = []
normalized_eea1_count = []

normalized_met_area = []
met_intensity_per_nuclei = []

treatment = []
time = []

masked_met = []
masked_eea1 = []

for img in images_path:
    print(f"DEBUG: img path = {img}") # Print the path of the current image for debugging purposes.
    raw_img = nd2.ND2Reader(img)
    ##############################################################
    #Segment nuclei
    ##############################################################
    
    nuclei = raw_img.get_frame_2D(c=1)

    from nuclear_segmentation import nuclear_segment
    min_nuclear_size = 500
    min_distance = 15 #minimum distance between nuclei centers to be considered separate nuclei
    gaussian_nuclei = 5
    dilation_radius = 5 #radius for dilation to ensure that the nuclei masks cover the entire nuclei and not just the centers
    label_nuclei = nuclear_segment(nuclei, min_nuclear_size, min_distance, gaussian_nuclei, dilation_radius)

    ##############################################################
    #Segment met signal
    ##############################################################

    met = raw_img.get_frame_2D(c=2)
    
    # Applying a gaussian filter
    from skimage.filters import gaussian, threshold_multiotsu
    transformed_met = gaussian(met, sigma = 2)

    from skimage.morphology import disk
    footprint = disk(5) # Create a disk-shaped footprint for the morphological operation

    from skimage.morphology import white_tophat
    tophat_met = white_tophat(transformed_met, footprint = footprint) # Apply the white top-hat filter to remove background noise

    # Thresholding and masking MET signal
    met_threshold = threshold_multiotsu(tophat_met, classes=3)[1] # Calculate a threshold value
    mask_met = np.zeros(tophat_met.shape)
    mask_met = tophat_met > met_threshold

    # Label the MET mask
    from skimage.measure import label
    labeled_met = label(mask_met)

    from skimage.measure import regionprops
    regions_met = regionprops(labeled_met) # Get properties of each region in the labeled MET regions.
    met_focicount = len(regions_met) # Count the number of MET foci

    regions_nuclei = regionprops(label_nuclei) # Get properties of each region in the labeled nuclei.
    nuclei_count = len(regions_nuclei) # Count the number of nuclei
    nuclear_count.append(nuclei_count) # Append the number of nuclei to the nuclear_count list.
    
    #Normalize MET foci count to nuclei count to account for differences in cell number across images.
    normalized_foci_count = met_focicount / nuclei_count
    met_count.append(met_focicount)
    normalized_met_count.append(normalized_foci_count)

    from skimage.color import label2rgb
    met_overlay = label2rgb(labeled_met, image=met*10, bg_label=0) # Create an overlay of the labeled MET signal and the original MET image.
    masked_met.append(met_overlay) # Append the MET overlay to the masked_met list.

    ##############################################################
    #Segment endosome eea1 marker signal
    ##############################################################
    
    #reading channel corresponding to endosome marker
    eea1 = raw_img.get_frame_2D(c=0)

    # Applying a gaussian filter
    transformed_eea1 = gaussian(eea1, sigma=2)
    
    # Apply the white top-hat filter to remove background noise
    tophat_eea1 = white_tophat(transformed_eea1, footprint = footprint)

    # Thresholding and masking endosome marker
    eea1_threshold = threshold_multiotsu(tophat_eea1, classes=3)[1]  # Calculate a threshold value
    mask_eea1 = np.zeros(tophat_eea1.shape)
    mask_eea1 = tophat_eea1 > eea1_threshold

    labeled_eea1 = label(mask_eea1) # Label the EEA1 mask
    regions_eea1 = regionprops(labeled_eea1) # Get properties of each region in the labeled EEA1 regions.
    eea1_focicount = len(regions_eea1) # Count the number of EEA1 foci
    normalized_eea1foci_count = eea1_focicount / nuclei_count # Normalize the EEA1 foci count to the number of nuclei

    eea1_count.append(eea1_focicount) # Append the number of EEA1 foci to the eea1_count list.
    normalized_eea1_count.append(normalized_eea1foci_count) # Append the normalized EEA1 foci count to the normalized_eea1_count list.
    
    eea1_overlay = label2rgb(labeled_eea1, image=eea1, bg_label=0) # Create an overlay of the labeled EEA1 signal and the original EEA1 image.
    masked_eea1.append(eea1_overlay) # Append the EEA1 overlay to the masked_eea1 list.

    from skimage.measure import intersection_coeff

    #Calculating intersection of fraction of MET signal colocalizing with endosome marker.
    overlap_img_level = intersection_coeff(mask_met, mask_eea1, mask = None)
    norm_overlap = overlap_img_level / nuclei_count # Normalize the overlap to the number of nuclei
    
    overlap_signal.append(overlap_img_level)
    normalized_overlap_signal.append(norm_overlap)
    
    filename.append(os.path.basename(img)) # Append the filename to the filename list.
    

    areas = []
    # Calculate the area of each MET focus
    for regions in regions_met:
        areas.append(regions.area)
    
    sum_area = sum(areas) # Sum the areas of all MET foci
    normalized_area = sum_area / nuclei_count # Normalize the total area to the number of nuclei
    normalized_met_area.append(normalized_area) # Append the normalized MET area to the normalized_met_area list.


    eea1_props = regionprops(labeled_eea1, intensity_image = met) # Get properties of each region in the labeled EEA1 regions, including intensity information.
    
    import numpy as np
    pixel_intensity = []
    
    # Calculate the pixel intensity of each EEA1 focus
    for eea1_prop in eea1_props:
            pixel_intensity.append(np.sum(eea1_prop.image_intensity))
    
    
    total_met_intensity = sum(pixel_intensity) # Sum the pixel intensities of all MET foci
    total_met_intensity_per_nuclei = total_met_intensity / nuclei_count # Normalize the total MET intensity to the number of nuclei
    met_intensity_per_nuclei.append(total_met_intensity_per_nuclei) # Append the normalized MET intensity to the met_intensity_per_nuclei list.


    #Selecting treatment groups based on file name, adjust as needed to match your file naming convention.
    if re.search('sNRP', os.path.basename(img), re.IGNORECASE):
        treatment.append('HGF + SNRP')
    elif re.search('noligand', os.path.basename(img), re.IGNORECASE):
        treatment.append('No Ligand')
    else:
        treatment.append('HGF only')

    # Extract the time point from the filename
    t = re.search(r'(\d+)\s*min', os.path.basename(img))
    if t:
        time.append(t.group(1))
    else:
        time.append(0)
    
import pandas as pd
#Saving data including masks and colocalization results.
coloc_data = pd.DataFrame({'Sample':filename, 'Colocalization_fraction':overlap_signal, 'Norm_colocalization': normalized_overlap_signal, 'Nuclei_count':nuclear_count, 
                           'MET_foci_count':met_count, 'Normalized_MET_foci_count':normalized_met_count, 'EEA1_focicount': eea1_count, 'Norm_EEA1_focicount': normalized_eea1_count, 
                           'Normalized_MET_area':normalized_met_area, 'MET_total_pixel_intensity_per_nuclei': met_intensity_per_nuclei, 'Treatment':treatment, 'Time':time})
coloc_data_path = os.path.join(directory_path, 'met_eea1_coloc.csv')
coloc_data.to_csv(coloc_data_path, index=False)

eea1_array = np.array(masked_eea1)
eea1_path = os.path.join(directory_path, 'masked_eea1.npy')
np.save(eea1_path, eea1_array) 

met_array = np.array(masked_met)
met_path = os.path.join(directory_path, 'masked_met.npy')
np.save(met_path, met_array)
