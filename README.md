# Cadaref

Cadaref is a key component of a [larger
pipeline](https://github.com/brawer/cadaref-zurich), initially built
for the City of Zürich, for automatic
[georeferencing](https://en.wikipedia.org/wiki/Georeferencing)
of historic [cadastral maps](https://en.wikipedia.org/wiki/Cadastre).

![scan](./doc/sample.png)

## Context

For most of the 20th century, each change to a land parcel in Zürich
got documented on a paper map. Likewise, as in the image above,
changes to building footprints were documented in a similar way.
Today, this record keeping is fully digital, but the City of Zürich
still has about 100K paper maps in its archive. This collection
documents a significant part of the city’s construction history.

To preserve this heritage, the archive was scanned to PDF, and then
each scan got processed by a computer system that tries to find the
precise geographic location for each map. As its output, Cadaref
produces [Cloud-Optimized GeoTIFF](http://cogeo.org/), an
industry-standard file format that is understood by Geographic
Information Systems and similar tools.


## Pipeline

After some image pre-processing (such as resolution enhancement and
thresholding), our [driver
pipeline](https://github.com/brawer/cadaref-zurich) uses Computer
Vision to recognize certain [cartographic
symbols](https://github.com/brawer/cadasym).  For example, in the map
above, the small white circles stand for one particular type of
boundary point; on the ground, their location is marked [small metal
disks](https://en.wikipedia.org/wiki/Survey_marker). When the Computer
Vision part of the pipeline looks at a map, it produces a [CSV file
with the recognized map symbols](testdata/symbols.csv), listing the
symbol type and image coordinates. (Since our Computer Vision
component works on an enhanced-resolution version of the image, the
resulting pixel coordinates can be fractional).

Separately, the processing pipeline also tries to find the
approximate geographic area of the map in question. For example,
[Optical Character Recognition](https://en.wikipedia.org/wiki/Optical_character_recognition) is used to extract text like parcel numbers.
Next, the pipeline queries today’s land survey database
for the survey points in the presumed map area, and stores
the result in another [CSV file with geographic points](testdata/points.csv).

Often enough, but not always, the symbols on the historical map
correspond to things (survey markers, fix points) that still
exist today. However, it is not obvious which map symbol corresponds
to what geographic point. We can only guess the rough geographic
area; many of the maps are rotated; neither symbol recognition nor
OCR are perfectly accurate; and the physical points may have been
demolished some time between the map was drawn and today.
To georeference the historical maps despite all these complications,
we need some kind of fuzzy match. This is exactly what Cadaref does.

GeoTIFF image with embedded transformation parameters.  While this
processing pipeline is very specific to the City of Zürich, the
matching tool might be useful for other projects, which is why we
release it independently. (See
[here](
pipeline that calls Cadaref for the Zürich project).


## Build

For performance reasons, the core algorithm of Cadaref is written
in the Rust programming language. To set up development on Linux,
have a look at the [Continuous Build](.github/workflows/ci.yml).


## Usage

Cadaref is a command-line tool that takes the following arguments:

* `--image` File path to the input image, in TIFF format
    like [this](testdata/HG3099.tif).

* `--page` Page number to process, in case the input TIFF has multiple pages.

* `--scales` Comma-separated list of map scales, such as `1:200,1:500`.
  For the Zürich project, we use OCR to extract the map scale that was
  printed on the map, with a fallback in case the map scale is missing.

* `--symbols` A set of map symbols, in CSV format like
  [this](testdata/symbols.csv). Typically detected by Computer Vision.
  Symbol locations are passed in (possibly fractional) pixel coordinates,
  relative to the top left of the image being processed.

* `--points` A set of points on the globe, in CSV format like
  [this](testdata/points.csv). Typically extracted from a database
  of survey markers, or whatever else the paper maps may depict.
  Locations are passed in geographic coordinates. (The tool currently
  assumes the Swiss CH1903+/LV95 spatial reference system, but this would
  be trivial to change).

* `--output` File path to the output image that will be written.
  The output will be a Cloud-Optimized GeoTIFF with embedded transformation
  parameters.

If everything goes well, the tool returns with status code 0.
In case of failures, in particular if the matching algorithm could
not find a good enough transformation, the tool returns with a non-zero
status code.
