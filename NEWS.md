# gcxgclab 1.1.0

This is a minor release adding new additional capabilities and providing some improvements and bug fixes.

## Bug fixes

* Fixes erros in 'align()' caused if all peaks found are above the given threshold, THR.

* Inconsequencial warnings are suppressed in 'align()'.

* Modification to the color scale ledgend for data of various scales.
  
## Improvements
  
* 'align()' has two additional optional inputs: use_ref_peak and ref_peak. Previously, 'align()' assumed the presense of toluene (present with themal desorption tubes using Tenax TA sorbent material). For users that do not have toluene present, a different initial reference m/z value can be input for ref_peak to provide an initial shift to the data for alignment. Alternatively, the boolean use_ref_peak can be set to FALSE to skip the initial shift to a reference peak completely. 

* 'plot_chr()' has two additional option inputs: xlab and ylab. This allows modification of the axes labeling on the hoizontal (x) and vertical (x) axis. 

* 'plot_peak()' has two additional option inputs: xlab and ylab. This allows modification of the axes labeling on the hoizontal (x) and vertical (x) axis. 

* 'batch_preprocess()' now has three addtional optional inputs: do_align, use_ref_peak, and ref_peak. do_align allows the user to change a boolean to skip alignment of the given list of files to preprocess if alignment is not necessary. It also includes the two new inputs use_ref_peak and ref_peak for 'align()'.

* 'extract_data()' has been modified to address irregularly shaped data. If the last column of data is incomplete (not a full modulation period), this data is now removed to allow plotting and analysis.

* Addtional commenting added to 'align()' for code organiziation and understandablility.
  
## New capabilities

* New function 'batch_extract()' which extracts data from all CFD files in a given folder path, using the methods of 'extract_data()'. Resulting data frames are returned in a list.

* New function 'TIC_integrate()' allows the user to integrate the area under the entire TIC curve. Previously area could only be calculated for a given peak. The user may also indicate starting and ending times for the region of integration.

# gcxgclab 1.0.1

This is a small release, adding vignette to the initial package submission, 1.0.0.
