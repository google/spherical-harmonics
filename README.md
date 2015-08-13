3D Graphics-oriented Spherical Harmonics Library
================================================

Spherical harmonics can be a tricky thing to wrap your head around.
Even once the basic theories are understood, there's some surprisingly
finicky implementation work to get the functions coded properly.
This is especially true when it comes to rotations of spherical harmonics
(much of the literature is math-dense and contains errata).
Additionally, different literature sources use slightly different
conventions when defining the basis functions.

This library is a collection of useful functions for working with
spherical harmonics. It is not restricted to a maximum order of basis
function, using recursive definitions for both the SH basis functions
and SH rotation matrices. This library uses the convention of 
including the Condon-Shortely phase function ((-1)^m) in the definition of 
the basis function.

This is not an official Google project.

**Dependencies**

This library depends on [Eigen3](http://eigen.tuxfamily.org) to for its
underlying linear algebra primitives. Colors are represented as 
`Eigen::Array3f`, where the components are ordered red, green, and blue.
[Google Test](https://code.google.com/p/googletest) is used for unit 
testing.

The [Bazel](http://bazel.io) build tool is used to build the library.
This is responsible for downloading and configuring Eigen3 and the
testing framework. You may build the library by executing in the root directory:

    bazel build //sh:spherical_harmonics

**General Functions**

See documentation in `sh/spherical_harmonics.h` for details on specific
functions. `sh/image.h` provides a very generic and simple image interface
that can be used to adapt this library with any actual imaging toolkit
already in use.

*Core SH Functions*

`EvalSH` - Evaluate the SH basis function of the given degree and order
at the provided position on a unit sphere. The position is described as
either a unit vector or as spherical coordinates.

`EvalSHSum` - Evaluate the approximation of a spherical function that
has already been converted to a vector of basis function coefficients.

*Projection Functions*

Used to estimate coefficients applied to basis functions to approximate
complex spherical functions as a weighted sum of the spherical harmonic
basis functions. Once projected, the returned coefficients can be
passed into `EvalSHSum`.

`ProjectFunction` - Project an analytic spherical function into every
basis function up to the specified order. This uses Monte Carlo 
integration to estimate the coefficient for each basis function.

`ProjectEnvironment` - Project an environment map image arranged in
a latitude-longitude projection into the basis functions up to the
specified order. This is a specialization of `ProjectFunction` that
is more efficient when the spherical function is described as an
image containing an environment.

`ProjectSparseSamples` - Project a spherical function that has only 
been sparsely evaluated (i.e. 10-50 times). Unlike the analytic 
function, this uses a least-squares fitting to best estimate the
coefficients for each basis function. This works well when fitting
to photographic data where there can only be so many photos captured.

*Diffuse Irradiance Functions*

Diffuse irradiance can be efficiently represented in low-order
spherical harmonics. It can be computed quickly by estimating
the standard diffuse cosine-lobe as a vector of coefficients,
and the environment as spherical harmonics. Diffuse irradiance
is simply the dot product of the two coefficient vectors.

`RenderDiffuseIrradiance` - Compute diffuse irradiance for a given
unit normal vector and SH coefficients that describe the environment
illumination (i.e. from `ProjectEnvironment`).

`RenderDiffuseIrradianceMap` - Compute diffuse irradiance for every
normal vector described by the texels of the provided latitude-longitude
image. This can be useful for computing a texture map of diffuse
irradiance and then transferring it to the GPU for shader-based rendering.

*Spherical Harmonic Rotations*

If a complex spherical function is rotated, and a set of spherical
harmonic coefficients is needed for this new function, it's possible
to rotate the spherical harmonic coefficients of the original approximation
rather than re-projecting the rotated function. This is often much more
efficient and is used in `RenderDiffuseIrradiance` to transform the cosine
lobe function for the unit z-axis to any other normal vector.

`Rotation` - Object type that computes the transformation matrices that
suitably transform spherical harmonic coefficients given a quaternion
rotation.

*Utility Functions*

`GetCoefficientCount` - Return the total number of coefficients needed to
represent all basis functions up to a given order.

`GetIndex` - Return a 1-dimensional index (suitable for accessing the returned
vectors from all the project functions) given a degree and order.

`ToVector` - Transform spherical coordinates into a unit vector.

`ToSphericalCoords` - Transform a unit vector into spherical coordinates.

`ImageXToPhi` - Transform a pixel's x coordinate in an image of a specific width
to the phi spherical coordinate.

`ImageYToTheta` - Transform a pixel's y coordinate in an image of a specific height
to the theta spherical coordinate.

`ToImageCoords` - Transform spherical coordinates into floating-point image coordinates
given particular image dimensions. The coordinates can be used to bilinearly 
interpolate an environment map, or cast to integers to access direct pixels.

**Literature**
The general spherical harmonic functions and fitting methods are from [1], the
environment map related functions are based on methods in [2] and [3], and 
spherical harmonic rotations are from [4] and [5]:

1. R. Green, "Spherical Harmonic Lighting: The Gritty Details", GDC 2003,
   http://www.research.scea.com/gdc2003/spherical-harmonic-lighting.pdf
2. R. Ramamoorthi and P. Hanrahan, "An Efficient Representation for
   Irradiance Environment Maps",. , P., SIGGRAPH 2001, 497-500
3. R. Ramamoorthi and P. Hanrahan, “On the Relationship between Radiance and
   Irradiance: Determining the Illumination from Images of a Convex
   Lambertian Object,” J. Optical Soc. Am. A, vol. 18, no. 10, pp. 2448-2459,
   2001.
4. J. Ivanic and K. Ruedenberg, "Rotation Matrices for Real Spherical
   Harmonics. Direct Determination by Recursion", J. Phys. Chem., vol. 100,
   no. 15, pp. 6342-6347, 1996. http://pubs.acs.org/doi/pdf/10.1021/jp953350u
5. Corrections to [4]: http://pubs.acs.org/doi/pdf/10.1021/jp9833350

