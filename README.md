# Cadaref

Cadaref is a tool for georeferencing scanned cadastral maps.

The tool was developed in 2024 as part of an Innovation Project
for the City of Zürich. For about a century, any changes to land
parcels were documented on paper, eventually adding up to some 100K
cadastral plans. Instead of browsing through moldy paper, it would be
nice to make these historical plans digitally accessible.

![scan](./doc/sample.png)

The Cadaref tool is part of a
larger [pipeline](https://github.com/brawer/cadaref-zurich)
that uses Computer Vision to recognize map symbols on the scanned images.
The pipeline also extracts the text (such as parcel numbers) written
on the maps, trying to guess the rough geographic area. Ultimately,
the pipeline calls Cadaref with the scanned image, a set of symbols
that were recognized on the image with Computer Vision, and a set of
known geographic points (such as survey markers) that may or may not
be depicted on the map. The pipeline also tries to guess the scale
of the map (typically from Optical Character Recognition, using a
fallback in case the map did not happen to indicate its scale).
Cadaref matches the two point sets against each other, and produces
a Cloud-Optimized GeoTIFF image with embedded transformation parameters.
While this processing pipeline is very specific to the City of Zürich,
the matching tool might be usefus on its own, so it is getting released
separately.


## Build

See [Continuous Build](.github/workflows/ci.yml).


## Usage

The cadaref tool takes the following arguments on the command line:

* `--image` File path to the input image, in TIFF format
    like [this](testdata/HG3099.tif).

* `--page` Page number to process, in case the input TIFF has multiple pages.

* `--symbols` A set of map symbols, in CSV format like
  [this](testdata/symbols.csv). Typically detected by Computer Vision.
  Symbol locations are passed in pixel coordinates,
  relative to the top left of the image being processed.

* `--symbols` A set of points on the globe, in CSV format like
  [this](testdata/points.csv). Typically extracted from a database
  of survey markers, or whatever else the paper maps may depict.
  Locations are passed in geographic coordinates. (The tool currently
  assumes the Swiss CH1903+/LV95 spatial reference system, but this would
  be trivial to change).

* `--output` File path to the output image that will be written.
  The output will be a Cloud-Optimized GeoTIFF with the transformation
  parameters embedded.

If everything goes well, the tool returns with status code 0.
In case of failures, in particular if the matching algorithm could
not find a good enough transformation, the tool returns with a non-zero
status code.
