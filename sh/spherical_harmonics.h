// Copyright 2015 Google Inc. All Rights Reserved.
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//    http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.

// The general spherical harmonic functions and fitting methods are from:
// 1. R. Green, "Spherical Harmonic Lighting: The Gritty Details", GDC 2003,
//    http://www.research.scea.com/gdc2003/spherical-harmonic-lighting.pdf
//
// The environment map related functions are based on the methods in:
// 2. R. Ramamoorthi and P. Hanrahan, "An Efficient Representation for
//    Irradiance Environment Maps",. , P., SIGGRAPH 2001, 497-500
// 3. R. Ramamoorthi and P. Hanrahan, “On the Relationship between Radiance and
//    Irradiance: Determining the Illumination from Images of a Convex
//    Lambertian Object,” J. Optical Soc. Am. A, vol. 18, no. 10, pp. 2448-2459,
//    2001.
//
// Spherical harmonic rotations are implemented using the recurrence relations
// described by:
// 4. J. Ivanic and K. Ruedenberg, "Rotation Matrices for Real Spherical
//    Harmonics. Direct Determination by Recursion", J. Phys. Chem., vol. 100,
//    no. 15, pp. 6342-6347, 1996.
//    http://pubs.acs.org/doi/pdf/10.1021/jp953350u
// 4b. Corrections to initial publication:
//    http://pubs.acs.org/doi/pdf/10.1021/jp9833350

#ifndef SH_SPHERICAL_HARMONICS_H_
#define SH_SPHERICAL_HARMONICS_H_

#include <array>
#include <vector>
#include <functional>
#include <memory>

#include "sh/image.h"

namespace sh {

// A spherical function, the first argument is phi, the second is theta.
// See EvalSH(int, int, double, double) for a description of these terms.
typedef std::function<double(double, double)> SphericalFunction;

const int kDefaultSampleCount = 10000;

// Get the total number of coefficients for a function represented by
// all spherical harmonic basis of degree <= @order (it is a point of
// confusion that the order of an SH refers to its degree and not the order).
constexpr int GetCoefficientCount(int order) {
  return (order + 1) * (order + 1);
}

// Get the one dimensional index associated with a particular degree @l
// and order @m. This is the index that can be used to access the Coeffs
// returned by SHSolver.
constexpr int GetIndex(int l, int m) {
  return l * (l + 1) + m;
}

// Convert from spherical coordinates to a direction vector. @phi represents
// the rotation about the Z axis and is from [0, 2pi]. @theta represents the
// angle down from the Z axis, from [0, pi].
Eigen::Vector3d ToVector(double phi, double theta);

// Convert from a direction vector to its spherical coordinates. The
// coordinates are written out to @phi and @theta. This is the inverse of
// ToVector.
// Check will fail if @dir is not unit.
void ToSphericalCoords(const Eigen::Vector3d& dir, double* phi, double* theta);

// Convert the (x, y) pixel coordinates into spherical coordinates (phi, theta)
// suitable for use with spherical harmonic evaluation. The x image axis maps
// to phi (0 to 2pi) and the y image axis maps to theta (0 to pi). A pixel index
// maps to the center of the pixel, so phi = 2pi (x + 0.5) / width and
// theta = pi (y + 0.5) / height. This is consistent with ProjectEnvironmentMap.
//
// x and y are not bounds checked against the image, but given the repeated
// nature of trigonometry functions, out-of-bounds x/y values produce reasonable
// phi and theta values (e.g. extending linearly beyond 0, pi, or 2pi).
// Results are undefined if the image dimensions are less than or equal to 0.
//
// The x and y functions are separated because they can be computed
// independently, unlike ToImageCoords.
double ImageXToPhi(int x, int width);
double ImageYToTheta(int y, int height);

// Convert the (phi, theta) spherical coordinates (using the convention
// defined spherical_harmonics.h) to pixel coordinates (x, y). The pixel
// coordinates are floating point to allow for later subsampling within the
// image. This is the inverse of ImageCoordsToSphericalCoords. It properly
// supports angles outside of the standard (0, 2pi) or (0, pi) range by mapping
// them back into it.
Eigen::Vector2d ToImageCoords(double phi, double theta, int width, int height);

// Evaluate the spherical harmonic basis function of degree @l and order @m
// for the given spherical coordinates, @phi and @theta.
// For low values of @l this will use a hard-coded function, otherwise it
// will fallback to EvalSHSlow that uses a recurrence relation to support all l.
double EvalSH(int l, int m, double phi, double theta);

// Evaluate the spherical harmonic basis function of degree @l and order @m
// for the given direction vector, @dir.
// Check will fail if @dir is not unit.
// For low values of @l this will use a hard-coded function, otherwise it
// will fallback to EvalSHSlow that uses a recurrence relation to support all l.
double EvalSH(int l, int m, const Eigen::Vector3d& dir);

// As EvalSH, but always uses the recurrence relationship. This is exposed
// primarily for testing purposes to ensure the hard-coded functions equal the
// recurrence relation version.
double EvalSHSlow(int l, int m, double phi, double theta);

// As EvalSH, but always uses the recurrence relationship. This is exposed
// primarily for testing purposes to ensure the hard-coded functions equal the
// recurrence relation version.
// Check will fail if @dir is not unit.
double EvalSHSlow(int l, int m, const Eigen::Vector3d& dir);

// Fit the given analytical spherical function to the SH basis functions
// up to @order. This uses Monte Carlo sampling to estimate the underlying
// integral. @sample_count determines the number of function evaluations
// performed. @sample_count is rounded to the greatest perfect square that
// is less than or equal to it.
//
// The samples are distributed uniformly over the surface of a sphere. The
// number of samples required to get a reasonable sampling of @func depends on
// the frequencies within that function. Lower frequency will not require as
// many samples. The recommended default kDefaultSampleCount should be
// sufficiently high for most functions, but is also likely overly conservative
// for many applications.
std::unique_ptr<std::vector<double>> ProjectFunction(
    int order, const SphericalFunction& func, int sample_count);

// Fit the given environment map to the SH basis functions up to @order.
// It is assumed that the environment map is parameterized by theta along
// the x-axis (ranging from 0 to 2pi after normalizing out the resolution),
// and phi along the y-axis (ranging from 0 to pi after normalization).
//
// This fits three different functions, one for each color channel. The
// coefficients for these functions are stored in the respective indices
// of the Array3f values of the returned vector.
std::unique_ptr<std::vector<Eigen::Array3f>> ProjectEnvironment(
    int order, const Image& env);

// Fit the given samples of a spherical function to the SH basis functions
// up to @order. This variant is used when there are relatively sparse
// evaluations or samples of the spherical function that must be fit and a
// regression is performed.
// @dirs and @values must have the same size. The directions in @dirs are
// assumed to be unit.
std::unique_ptr<std::vector<double>> ProjectSparseSamples(
    int order, const std::vector<Eigen::Vector3d>& dirs,
    const std::vector<double>& values);

// Evaluate the already computed coefficients for the SH basis functions up
// to @order, at the coordinates @phi and @theta. The length of the @coeffs
// vector must be equal to GetCoefficientCount(order).
// There are explicit instantiations for double, float, and Eigen::Array3f.
template <typename T>
T EvalSHSum(int order, const std::vector<T>& coeffs, double phi, double theta);

// As EvalSHSum, but inputting a direction vector instead of spherical coords.
// Check will fail if @dir is not unit.
template <typename T>
T EvalSHSum(int order, const std::vector<T>& coeffs, 
            const Eigen::Vector3d& dir);

// Render into @diffuse_out the diffuse irradiance for every normal vector
// representable in @diffuse_out, given the luminance stored in @env_map.
// Both @env_map and @diffuse_out use the latitude-longitude projection defined
// specified in ImageX/YToPhi/Theta. They may be of different
// resolutions. The resolution of @diffuse_out must be set before invoking this
// function.
void RenderDiffuseIrradianceMap(const Image& env_map, 
                                Image* diffuse_out);

// Render into @diffuse_out diffuse irradiance for every normal vector
// representable in @diffuse_out, for the environment represented as the given
// spherical harmonic coefficients, @sh_coeffs. The resolution of
// @diffuse_out must be set before calling this function. Note that a high
// resolution is not necessary (64 x 32 is often quite sufficient).
// See RenderDiffuseIrradiance for how @sh_coeffs is interpreted.
void RenderDiffuseIrradianceMap(const std::vector<Eigen::Array3f>& sh_coeffs,
                                Image* diffuse_out);

// Compute the diffuse irradiance for @normal given the environment represented
// as the provided spherical harmonic coefficients, @sh_coeffs. Check will
// fail if @normal is not unit length. @sh_coeffs can be of any length. Any
// coefficient provided beyond what's used internally to represent the diffuse
// lobe (9) will be ignored. If @sh_coeffs is less than 9, the remaining
// coefficients are assumed to be 0. This naturally means that providing an
// empty coefficient array (e.g. when the environment is assumed to be black and
// not provided in calibration) will evaluate to 0 irradiance.
Eigen::Array3f RenderDiffuseIrradiance(
    const std::vector<Eigen::Array3f>& sh_coeffs, 
    const Eigen::Vector3d& normal);

class Rotation {
 public:
  // Create a new Rotation that can applies @rotation to sets of coefficients
  // for the given @order. @order must be at least 0.
  static std::unique_ptr<Rotation> Create(int order,
                                          const Eigen::Quaterniond& rotation);

  // Create a new Rotation that applies the same rotation as @rotation. This
  // can be used to efficiently calculate the matrices for the same 3x3
  // transform when a new order is necessary.
  static std::unique_ptr<Rotation> Create(int order, const Rotation& rotation);

  // Transform the SH basis coefficients in @coeff by this rotation and store
  // them into @result. These may be the same vector. The @result vector will
  // be resized if necessary, but @coeffs must have its size equal to
  // GetCoefficientCount(order()).
  //
  // This rotation transformation produces a set of coefficients that are equal
  // to the coefficients found by projecting the original function rotated by
  // the same rotation matrix.
  //
  // There are explicit instantiations for double, float, and Array3f.
  template <typename T>
  void Apply(const std::vector<T>& coeffs, std::vector<T>* result) const;

  // The order (0-based) that the rotation was constructed with. It can only
  // transform coefficient vectors that were fit using the same order.
  int order() const;

  // Return the rotation that is effectively applied to the inputs of the
  // original function.
  Eigen::Quaterniond rotation() const;

  // Return the (2l+1)x(2l+1) matrix for transforming the coefficients within
  // band @l by the rotation. @l must be at least 0 and less than or equal to
  // the order this rotation was initially constructed with.
  const Eigen::MatrixXd& band_rotation(int l) const;

 private:
  explicit Rotation(int order, const Eigen::Quaterniond& rotation);

  const int order_;
  const Eigen::Quaterniond rotation_;

  std::vector<Eigen::MatrixXd> band_rotations_;
};

}  // namespace sh

#endif  // SH_SPHERICAL_HARMONICS_H_
