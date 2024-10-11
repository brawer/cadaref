# Algorithm

The matching of map symbols against geographic points is inspired by
[Random Sample Consensus
(RANSAC)](https://en.wikipedia.org/wiki/Random_sample_consensus), a
classic algorithm developed in the ealy 1980s.

However, other than with classic RANSAC, we do not randomly select
potential inliers. Instead, we generate and test our candidates from
point pairs that are located maximally far away from each other.  On a
map, it is intuitively better to interpolate between far-away points
than it is to extrapolate from points very near to each other; this
reduces measurement errors.

## Input

* A map image [like this](../testdata/HG3099.tif). The matching
  algorithm only needs the image *resolution* in pixels per inch or
  centimeter. With the TIFF image format, this is part of the image
  metadata.

* A set of cartographic symbols [like this](../testdata/symbols.csv).
  In the project for the City of Zürich, we automatically detect
  cartographic symbols on historical maps by means of Computer Vision,
  but this does not matter for the matching algorithm.  The pixel
  coordinates may also be fractional.

* A set of geographic points [like this](../testdata/points.csv).  In
  the Zürich project, we took these points from a 2007 backup file of
  the land survey register, which was the year when the city switched
  to a completely digital workflow. For each [survey
  marker](https://en.wikipedia.org/wiki/Survey_marker), the database
  lists its construction date, so we only pass those points to the
  algorithm that already existed at the time the historical map was
  made. However, the matching algorithm does not really care about the
  source or meaning of these points; it just works on coordinates. (In
  the above CSV file, the points have an ID, but this is only used for
  tracing and debugging).

* A set of map scales, such as 1:200 or 1:5000. In the archived
  historical maps for Zürich, most maps indicated their scale, and we
  could extract the text by means of Optical Character
  Recognition. Some maps, however, had no scale indication, and
  sometimes OCR failed to read the text. In those cases, we assumed a
  set of typical map scales as a fallback.  However, all this does not
  matter for the matching algorithm, it just needs a set of candidate
  scales as its input.


## Output

If the matching algorithm is successful, it returns an [affine
transformation](https://en.wikipedia.org/wiki/Affine_transformation)
that projects pixel coordinates of the scanned historical map to
geographic coordinates of the real world. For some inputs, however,
the algorithm will find no result and give up.


## Preparation

In a preparation step, we build the following data structures:

1. An array with cartographic symbols, sorted by
Euclidean distance to the top left corner (origin)
of the scanned map image.

2. A sorted table with the distance (in millimeters) between every pair
of geographic points. We will use this table to find all point pairs
whose distance is within a buffer (±1 meter) around a query distance.
To implement this range query, we do a binary search for the lower limit.
Then, we iterate over the table entries until we raech the upper limit.

3. An [R*-tree](https://en.wikipedia.org/wiki/R*-tree) for efficiently
finding geographic points near a query location.


## Matching

After preparation, the algorithm works like this. In the implementation,
many steps are interleaved; we do not actually need to expand all those
various pairs in memory. But the algorithm is easier to understand if the
steps are described separately.

1. We pair the 30 symbols that are closest to the map’s top left
corner with the 30 symbols farthest away. This gives 900 symbols
pairs. (Should the input map contain less than 30 symbols, we use a
smaller limit).

2. For each of the 900 symbol pairs on the historical map, we measure
their distance in meters. We can do this because we know both the map
scale and the image resolution. For example, on a 1:500 map scanned at
300 dpi, the width of a single map pixel corresponds to a real-world
distance of 42.3 millimeters; if two symbols on this map have a
distance of 791 raster pixels, the depicted locations would be
33486 millimeters away from each other.

3. Next, we look for pairs of geographic points with roughly this
distance.  To allow for drawing and scanning inaccuracies, we search
for distances within a range of plus or minus one meter of the value
computed in step 2.  In our example, we might find 7 point pair whose
distance is between 32486 and 34486 millimeters.

4. Given a particular map symbol pair of step 2 (let’s call them A and
B), and a geographic point pair of step 3 (let’s call them P and Q),
we compute an
[affine transformation](https://en.wikipedia.org/wiki/Affine_transformation)
for projecting A/B to P/Q. The historical maps can be arbitrarily
rotated, so we do not know which point matches which. Therefore, we
generate _two_ hypothetical affine transformations in this step: One
for projecting A⟶P and B⟶Q, and another for A⟶Q and B⟶P.

5. In RANSAC terminology, the affine transformations of step 4 are
estimation parameters for “models”. Let’s find out how good each
model is, and pick the best. We apply each affine transformation
(“model hypothesis”) of step 4 to every map symbol, which gives use a
projected location in geodetic coordinates. Then, we query the R*-tree
to find the survey point that is closest to the projected location.
If the distance is less than 1 meter, we have a match or (in RANSAC
terms) an “inlier”.

6. We prefer the model with the most “inliers”. In other words, we take
the affine transformation that projects the most cartographic symbols
to a location within 1 meter of a survey marker that already existed
at the time the historical map had been drawn. If two models produce
the same number of inliers, we prefer the model with the smaller summed
squared distances between projected location and matching geographic point.

7. To prevent garbage results, we require that at least 50% of map symbols
match some survey point. Otherwise, if there’s too many outliers,
we give up.

8. By construction, our transformation will be perfect for the matches
A–P and B–Q (or A–Q and B–P) chosen in step 4, and less than 1 meter
off for all other inliers. We try to redistribute the error residual
to all inliers. Actually, now that we have established which map symbol
corresponds to what survey marker, we have reduced our problem to
georeferencing imagery with Ground Control Points. This is a standard
task in geographic information processing, and we can simply
delegate this to the [GDAL library](https://gdal.org/) which has
an implementation in
[GCPsToGeoTransform()](https://gdal.org/en/release-3.9/doxygen/gdal_8h.html#ae6bc0eeea40d1645fbd44d7431c8db07). We compare the quality of GDAL’s
re-estimation with our own transform and take the better of the two.
