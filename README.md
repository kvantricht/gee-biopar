# Sentinel-2 biopar computation in Google Earth Engine

## Background

This repository contains an implementation of Sentinel-2 biophysical parameters in Google Earth Engine.
Computation of these parameters is based on the S2ToolBox methodology. For algorithm details, see the [original ATBD](https://step.esa.int/docs/extra/ATBD_S2ToolBox_L2B_V1.1.pdf).
Porting of the algorithm to Google Earth Engine is mostly based on the implementation from the [Sentinelhub custom scripts repository](https://github.com/sentinel-hub/custom-scripts/tree/master/sentinel-2).

Currently supported parameters are:
- 3-band FAPAR
- 8-band FAPAR
- 3-band LAI
- 8-band LAI

fCOVER, CCC and CWC can be done as well. Input for the biophysical parameter computation should always be Sentinel-2 L2A products.

## Example

A Google Earth Engine example script showing how to use the Sentinel-2 biopar computation method can be consulted [here](https://code.earthengine.google.com/497f61f10af651c585254c14f24de020)
