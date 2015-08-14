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

#include "sh/default_image.h"
#include "sh/spherical_harmonics.h"
#include "gtest/gtest.h"

namespace sh {

namespace {

#define EXPECT_TUPLE3_NEAR(expected, actual, tolerance) \
  { \
    EXPECT_NEAR(expected(0), actual(0), tolerance); \
    EXPECT_NEAR(expected(1), actual(1), tolerance); \
    EXPECT_NEAR(expected(2), actual(2), tolerance); \
  }

#define EXPECT_TUPLE2_NEAR(expected, actual, tolerance) \
  { \
    EXPECT_NEAR(expected(0), actual(0), tolerance); \
    EXPECT_NEAR(expected(1), actual(1), tolerance); \
  }

const double kEpsilon = 1e-10;
const double kHardcodedError = 1e-5;
const double kCoeffErr = 5e-2;
// Use a lower sample count than the default so the tests complete faster.
const int kTestSampleCount = 5000;

// Use a very small image since computing the diffuse irradiance explicitly
// is an O(N^4) operation on the resolution.
const int kImageWidth = 32;
const int kImageHeight = 16;

const double kIrradianceError = 0.01;

// Because the test is limited to an environment of 32x16, the error is high
// but most irradiance values are on the order of 2 or 3 so relatively this is
// reasonable. This error factor is also dependent on the SH order used to
// represent the cosine lobe in irradiance calculations. Band 2 does add some
// ringing caused by the clamping function applied to the lobe.
const double kEnvMapIrradianceError = 0.08;

// Clamp the first argument to be greater than or equal to the second
// and less than or equal to the third.
double Clamp(double val, double min, double max) {
  if (val < min) {
    val = min;
  }
  if (val > max) {
    val = max;
  }
  return val;
}

// Return true if the first value is within epsilon of the second value.
bool NearByMargin(double actual, double expected) {
  double diff = actual - expected;
  if (diff < 0.0) {
    diff = -diff;
  }
  return diff < 1e-16;
}

void ExpectMatrixNear(const Eigen::MatrixXd& expected,
                      const Eigen::MatrixXd& actual, double tolerance) {
  EXPECT_EQ(expected.rows(), actual.rows());
  EXPECT_EQ(expected.cols(), actual.cols());

  for (int i = 0; i < expected.rows(); i++) {
    for (int j = 0; j < expected.cols(); j++) {
      EXPECT_NEAR(expected(i, j), actual(i, j), tolerance);
    }
  }
}

void GenerateTestEnvironment(Image* env_map) {
  for (int y = 0; y < kImageHeight; y++) {
    for (int x = 0; x < kImageWidth; x++) {
      float red = x < kImageWidth / 2 ? 1.0 : 0.0;
      float green = x >= kImageWidth / 2 ? 1.0 : 0.0;
      float blue = y > kImageHeight / 2 ? 1.0 : 0.0;

      env_map->SetPixel(x, y, Eigen::Array3f(red, green, blue));
    }
  }
}

void ComputeExplicitDiffuseIrradiance(const Image& env_map,
                                      Image* diffuse) {
  double pixel_area = 2 * M_PI / env_map.width() * M_PI / env_map.height();
  for (int y = 0; y < kImageHeight; y++) {
    for (int x = 0; x < kImageWidth; x++) {
      Eigen::Vector3d normal = ToVector((x + 0.5) * 2 * M_PI / kImageWidth,
                                        (y + 0.5) * M_PI / kImageHeight);
      Eigen::Array3f irradiance(0.0, 0.0, 0.0);
      for (int ey = 0; ey < env_map.height(); ey++) {
        double theta = (ey + 0.5) * M_PI / env_map.height();
        double sa = pixel_area * sin(theta);
        for (int ex = 0; ex < env_map.width(); ex++) {
          Eigen::Vector3d light = ToVector(
              (ex + 0.5) * 2 * M_PI / env_map.width(),
              (ey + 0.5) * M_PI / env_map.height());
          irradiance += sa * Clamp(light.dot(normal), 0.0, 1.0) *
              env_map.GetPixel(ex, ey);
        }
      }
      diffuse->SetPixel(x, y, irradiance);
    }
  }
}

}  // namespace

TEST(SphericalHarmonicsTest, ProjectFunction) {
  // The expected coefficients used to define the analytic spherical function
  const std::vector<double> coeffs = {-1.028, 0.779, -0.275, 0.601, -0.256,
                                      1.891, -1.658, -0.370, -0.772};

  // Project and compare the fitted coefficients, which should be near identical
  // to the initial coefficients
  SphericalFunction func = [&] (double phi, double theta) {
    return EvalSHSum(2, coeffs, phi, theta); };
  std::unique_ptr<std::vector<double>> fitted = ProjectFunction(
      2, func, kTestSampleCount);
  ASSERT_TRUE(fitted != nullptr);

  for (int i = 0; i < 9; i++) {
    EXPECT_NEAR(coeffs[i], (*fitted)[i], kCoeffErr);
  }
}

TEST(SphericalHarmonicsTest, ProjectSparseSamples) {
  // These are the expected coefficients that define the sparse samples of
  // the underyling spherical function
  const std::vector<double> coeffs = {-0.591, -0.713, 0.191, 1.206, -0.587,
                                      -0.051, 1.543, -0.818, 1.482};

  // Generate sparse samples
  std::vector<Eigen::Vector3d> sample_dirs;
  std::vector<double> sample_vals;
  for (int t = 0; t < 6; t++) {
    double theta = t * M_PI / 6.0;
    for (int p = 0; p < 8; p++) {
      double phi = p * 2.0 * M_PI / 8.0;
      Eigen::Vector3d dir = ToVector(phi, theta);
      double value = EvalSHSum(2, coeffs, phi, theta);
      sample_dirs.push_back(dir);
      sample_vals.push_back(value);
    }
  }

  // Compute the sparse fit and given that the samples were drawn from the
  // spherical basis functions this should be a pretty ideal match
  std::unique_ptr<std::vector<double>> fitted = ProjectSparseSamples(
      2, sample_dirs, sample_vals);
  ASSERT_TRUE(fitted != nullptr);

  for (int i = 0; i < 9; i++) {
    EXPECT_NEAR(coeffs[i], (*fitted)[i], kCoeffErr);
  }
}

TEST(SphericalHarmonicsTest, ProjectEnvironment) {
  // These are the expected coefficients that define the environment map
  // passed into Project()
  const std::vector<double> c_red = {-1.028, 0.779, -0.275, 0.601, -0.256,
                                     1.891, -1.658, -0.370, -0.772};
  const std::vector<double> c_green = {-0.591, -0.713, 0.191, 1.206, -0.587,
                                       -0.051, 1.543, -0.818, 1.482};
  const std::vector<double> c_blue = {-1.119, 0.559, 0.433, -0.680, -1.815,
                                      -0.915, 1.345, 1.572, -0.622};

  // Generate an environment map based off of c_red, c_green, and c_blue
  DefaultImage env_map(64, 32);  // This does not need to be a large map for the test
  for (int t = 0; t < env_map.height(); t++) {
    double theta = (t + 0.5) * M_PI / env_map.height();
    for (int p = 0; p < env_map.width(); p++) {
      double phi = (p + 0.5) * 2.0 * M_PI / env_map.width();
      env_map.SetPixel(p, t, Eigen::Array3f(EvalSHSum(2, c_red, phi, theta),
                                            EvalSHSum(2, c_green, phi, theta),
                                            EvalSHSum(2, c_blue, phi, theta)));
    }
  }

  // Fit the environment to spherical functions. Given that we formed it from
  // the spherical basis we should get a very near perfect fit.
  std::unique_ptr<std::vector<Eigen::Array3f>> c_fit = ProjectEnvironment(
      2, env_map);
  ASSERT_TRUE(c_fit != nullptr);

  for (int i = 0; i < 9; i++) {
    Eigen::Array3f fitted = (*c_fit)[i];
    EXPECT_NEAR(c_red[i], fitted(0), kCoeffErr);
    EXPECT_NEAR(c_green[i], fitted(1), kCoeffErr);
    EXPECT_NEAR(c_blue[i], fitted(2), kCoeffErr);
  }
}

TEST(SphericalHarmonicsTest, EvalSHSum) {
  const std::vector<double> coeffs = {-1.119, 0.559, 0.433, -0.680, -1.815,
                                      -0.915, 1.345, 1.572, -0.622};
  double expected = 0.0;
  for (int l = 0; l <= 2; l++) {
    for (int m = -l; m <= l; m++) {
      expected += coeffs[GetIndex(l, m)] * EvalSH(l, m, M_PI / 4, M_PI / 4);
    }
  }

  EXPECT_EQ(expected, EvalSHSum(2, coeffs, M_PI / 4, M_PI / 4));
}

TEST(SphericalHarmonicsTest, EvalSHSumArray3f) {
  const std::vector<double> coeffs_0 = {-1.119, 0.559, 0.433, -0.680, -1.815,
                                        -0.915, 1.345, 1.572, -0.622};
  const std::vector<double> coeffs_1 = {-0.591, -0.713, 0.191, 1.206, -0.587,
                                        -0.051, 1.543, -0.818, 1.482};
  const std::vector<double> coeffs_2 = {-1.119, 0.559, 0.433, -0.680, -1.815,
                                        -0.915, 1.345, 1.572, -0.622};
  std::vector<Eigen::Array3f> coeffs;
  for (unsigned int i = 0; i < coeffs_0.size(); i++) {
    coeffs.push_back(Eigen::Array3f(coeffs_0[i], coeffs_1[i], coeffs_2[i]));
  }

  Eigen::Array3f expected(0.0, 0.0, 0.0);
  for (int l = 0; l <= 2; l++) {
    for (int m = -l; m <= l; m++) {
      expected += coeffs[GetIndex(l, m)] * EvalSH(l, m, M_PI / 4, M_PI / 4);
    }
  }

  Eigen::Array3f actual = EvalSHSum(2, coeffs, M_PI / 4, M_PI / 4);
  EXPECT_TUPLE3_NEAR(expected, actual, kEpsilon);
}

TEST(SphericalHarmonicsTest, GetIndex) {
  // Indices are arranged from low band to high degree, and from low order
  // to high order within a band.
  EXPECT_EQ(0, GetIndex(0, 0));
  EXPECT_EQ(1, GetIndex(1, -1));
  EXPECT_EQ(2, GetIndex(1, 0));
  EXPECT_EQ(3, GetIndex(1, 1));
  EXPECT_EQ(4, GetIndex(2, -2));
  EXPECT_EQ(5, GetIndex(2, -1));
  EXPECT_EQ(6, GetIndex(2, 0));
  EXPECT_EQ(7, GetIndex(2, 1));
  EXPECT_EQ(8, GetIndex(2, 2));
}

TEST(SphericalHarmonicsTest, GetCoefficientCount) {
  // For up to order n SH representation, there are (n+1)^2 coefficients.
  EXPECT_EQ(1, GetCoefficientCount(0));
  EXPECT_EQ(9, GetCoefficientCount(2));
  EXPECT_EQ(16, GetCoefficientCount(3));
}

TEST(SphericalHarmonicsTest, ToVector) {
  // Compare spherical coordinates with their known direction vectors.
  EXPECT_TUPLE3_NEAR(Eigen::Vector3d(1, 0, 0), ToVector(0.0, M_PI / 2), 
                     kEpsilon);
  EXPECT_TUPLE3_NEAR(Eigen::Vector3d(0, 1, 0), ToVector(M_PI / 2, M_PI / 2), 
                     kEpsilon);
  EXPECT_TUPLE3_NEAR(Eigen::Vector3d(0, 0, 1), ToVector(0.0, 0.0), kEpsilon);
  EXPECT_TUPLE3_NEAR(Eigen::Vector3d(0.5, 0.5, sqrt(0.5)),
                     ToVector(M_PI / 4, M_PI / 4), kEpsilon);
  EXPECT_TUPLE3_NEAR(Eigen::Vector3d(0.5, 0.5, -sqrt(0.5)),
                     ToVector(M_PI / 4, 3 * M_PI / 4), kEpsilon);
  EXPECT_TUPLE3_NEAR(Eigen::Vector3d(-0.5, 0.5, -sqrt(0.5)),
                     ToVector(3 * M_PI / 4, 3 * M_PI / 4), kEpsilon);
  EXPECT_TUPLE3_NEAR(Eigen::Vector3d(0.5, -0.5, -sqrt(0.5)),
                     ToVector(-M_PI / 4, 3 * M_PI / 4), kEpsilon);
}

TEST(SphericalHarmonicsTest, ToSphericalCoords) {
  // Compare vectors with their known spherical coordinates.
  double phi, theta;
  ToSphericalCoords(Eigen::Vector3d(1, 0, 0), &phi, &theta);
  EXPECT_EQ(0.0, phi);
  EXPECT_EQ(M_PI / 2, theta);
  ToSphericalCoords(Eigen::Vector3d(0, 1, 0), &phi, &theta);
  EXPECT_EQ(M_PI / 2, phi);
  EXPECT_EQ(M_PI / 2, theta);
  ToSphericalCoords(Eigen::Vector3d(0, 0, 1), &phi, &theta);
  EXPECT_EQ(0.0, phi);
  EXPECT_EQ(0.0, theta);
  ToSphericalCoords(Eigen::Vector3d(0.5, 0.5, sqrt(0.5)), &phi, &theta);
  EXPECT_EQ(M_PI / 4, phi);
  EXPECT_EQ(M_PI / 4, theta);
  ToSphericalCoords(Eigen::Vector3d(0.5, 0.5, -sqrt(0.5)), &phi, &theta);
  EXPECT_EQ(M_PI / 4, phi);
  EXPECT_EQ(3 * M_PI / 4, theta);
  ToSphericalCoords(Eigen::Vector3d(-0.5, 0.5, -sqrt(0.5)), &phi, &theta);
  EXPECT_EQ(3 * M_PI / 4, phi);
  EXPECT_EQ(3 * M_PI / 4, theta);
  ToSphericalCoords(Eigen::Vector3d(0.5, -0.5, -sqrt(0.5)), &phi, &theta);
  EXPECT_EQ(-M_PI / 4, phi);
  EXPECT_EQ(3 * M_PI / 4, theta);
}

TEST(SphericalHarmonicsTest, EvalSHSlow) {
  // Compare the general SH implementation to the closed form functions for
  // several bands, from: http://en.wikipedia.org/wiki/Table_of_spherical_harmonics#Real_spherical_harmonics
  // It's assumed that if the implementation matches these for this subset, the
  // probability it's correct overall is high.
  //
  // Note that for all cases |m|=1 below, we negate compared to what Wikipedia
  // lists. After careful review, it seems they do not include the (-1)^m term
  // (the Condon-Shortley phase) in their calculations.
  const double phi = M_PI / 4;
  const double theta = M_PI / 3;
  const Eigen::Vector3d d = ToVector(phi, theta);

  // l = 0
  EXPECT_NEAR(0.5 * sqrt(1 / M_PI), EvalSHSlow(0, 0, phi, theta), kEpsilon);

  // l = 1, m = -1
  EXPECT_NEAR(-sqrt(3 / (4 * M_PI)) * d.y(), EvalSHSlow(1, -1, phi, theta),
              kEpsilon);
  // l = 1, m = 0
  EXPECT_NEAR(sqrt(3 / (4 * M_PI)) * d.z(), EvalSHSlow(1, 0, phi, theta),
              kEpsilon);
  // l = 1, m = 1
  EXPECT_NEAR(-sqrt(3 / (4 * M_PI)) * d.x(), EvalSHSlow(1, 1, phi, theta),
              kEpsilon);

  // l = 2, m = -2
  EXPECT_NEAR(0.5 * sqrt(15 / M_PI) * d.x() * d.y(),
              EvalSHSlow(2, -2, phi, theta), kEpsilon);
  // l = 2, m = -1
  EXPECT_NEAR(-0.5 * sqrt(15 / M_PI) * d.y() * d.z(),
              EvalSHSlow(2, -1, phi, theta), kEpsilon);
  // l = 2, m = 0
  EXPECT_NEAR(0.25 * sqrt(5 / M_PI) *
              (-d.x() * d.x() - d.y() * d.y() + 2 * d.z() * d.z()),
              EvalSHSlow(2, 0, phi, theta), kEpsilon);
  // l = 2, m = 1
  EXPECT_NEAR(-0.5 * sqrt(15 / M_PI) * d.z() * d.x(),
              EvalSHSlow(2, 1, phi, theta), kEpsilon);
  // l = 2, m = 2
  EXPECT_NEAR(0.25 * sqrt(15 / M_PI) * (d.x() * d.x() - d.y() * d.y()),
              EvalSHSlow(2, 2, phi, theta), kEpsilon);
}

TEST(SphericalHarmonicsTest, EvalSHHardcoded) {
  // Arbitrary coordinates
  const double phi = 0.4296;
  const double theta = 1.73234;
  const Eigen::Vector3d d = ToVector(phi, theta);

  for (int l = 0; l <= 4; l++) {
    for (int m = -l; m <= l; m++) {
      double expected = EvalSHSlow(l, m, phi, theta);
      EXPECT_NEAR(expected, EvalSH(l, m, phi, theta), kHardcodedError)
          << "Spherical coord failure for l, m = (" << l << ", " << m << ")";
      EXPECT_NEAR(expected, EvalSH(l, m, d), kHardcodedError)
          << "Vector failure for l, m = (" << l << ", " << m << ")";
    }
  }
}

TEST(SphericalHarmonicsDeathTest, EvalSHBadInputs) {
  const double phi = M_PI / 4;
  const double theta = M_PI / 3;
  const Eigen::Vector3d d = ToVector(phi, theta);

  // l < 0
  EXPECT_DEATH(EvalSH(-1, 0, phi, theta), "l must be at least 0.");
  EXPECT_DEATH(EvalSH(-1, 0, d), "l must be at least 0.");

  // m > l
  EXPECT_DEATH(EvalSH(1, 2, phi, theta), "m must be between -l and l.");
  EXPECT_DEATH(EvalSH(1, 2, d), "m must be between -l and l.");

  // m < -l
  EXPECT_DEATH(EvalSH(1, -2, phi, theta), "m must be between -l and l.");
  EXPECT_DEATH(EvalSH(1, -2, phi, theta), "m must be between -l and l.");
}

TEST(SphericalHarmonicsDeathTest, ProjectFunctionBadInputs) {
  const std::vector<double> coeffs = {-1.028};

  SphericalFunction func = [&] (double phi, double theta) {
    return EvalSHSum(0, coeffs, phi, theta); };

  // order < 0
  EXPECT_DEATH(ProjectFunction(-1, func, kTestSampleCount),
               "Order must be at least zero.");
  // sample_count <= 0
  EXPECT_DEATH(ProjectFunction(2, func, 0),
               "Sample count must be at least one.");
  EXPECT_DEATH(ProjectFunction(2, func, -1),
               "Sample count must be at least one.");
}

TEST(SphericalHarmonicsDeathTest, ProjectEnvironmentBadInputs) {
  DefaultImage env(64, 32);

  // order < 0
  EXPECT_DEATH(ProjectEnvironment(-1, env), "Order must be at least zero.");
}

TEST(SphericalHarmonicsDeathTest, ProjectSparseSamplesBadInputs) {
  // These are the expected coefficients that define the sparse samples of
  // the underyling spherical function
  const std::vector<double> coeffs = {-0.591};

  // Generate sparse samples
  std::vector<Eigen::Vector3d> sample_dirs;
  std::vector<double> sample_vals;
  for (int t = 0; t < 6; t++) {
    double theta = t * M_PI / 6.0;
    for (int p = 0; p < 8; p++) {
      double phi = p * 2.0 * M_PI / 8.0;
      Eigen::Vector3d dir = ToVector(phi, theta);
      double value = EvalSHSum<double>(0, coeffs, phi, theta);
      sample_dirs.push_back(dir);
      sample_vals.push_back(value);
    }
  }

  // order < 0
  EXPECT_DEATH(ProjectSparseSamples(-1, sample_dirs, sample_vals),
               "Order must be at least zero.");

  // unequal directions and values
  sample_vals.push_back(0.4);
  EXPECT_DEATH(ProjectSparseSamples(2, sample_dirs, sample_vals),
               "Directions and values must have the same size.");
}

TEST(SphericalHarmonicsDeathTest, EvalSHSumBadInputs) {
  // These are the expected coefficients that define the sparse samples of
  // the underyling spherical function
  const std::vector<double> coeffs = {-0.591, -0.713, 0.191, 1.206, -0.587,
                                      -0.051, 1.543, -0.818, 1.482};
  EXPECT_DEATH(EvalSHSum(3, coeffs, M_PI / 4, M_PI / 4),
               "Incorrect number of coefficients provided.");
}

TEST(SphericalHarmonicsDeathTest, ToSphericalCoordsBadInputs) {
  double phi, theta;
  EXPECT_DEATH(ToSphericalCoords(Eigen::Vector3d(2.0, 0.0, 0.4), &phi, &theta),
               "");
}

TEST(SphericalHarmonicsRotationTest, ClosedFormZAxisRotation) {
  // The band-level rotation matrices for a rotation about the z-axis are
  // relatively simple so we can compute them closed form and make sure the
  // recursive general approach works properly.
  // This closed form comes from [1].
  double alpha = M_PI / 4.0;
  Eigen::Quaterniond rz(Eigen::AngleAxisd(alpha, Eigen::Vector3d::UnitZ()));
  std::unique_ptr<Rotation> rz_sh(Rotation::Create(3, rz));

  // order 0
  Eigen::MatrixXd r0(1, 1);
  r0 << 1.0;
  ExpectMatrixNear(r0, rz_sh->band_rotation(0), kEpsilon);

  // order 1
  Eigen::MatrixXd r1(3, 3);
  r1 << cos(alpha), 0, sin(alpha),
        0, 1, 0,
        -sin(alpha), 0, cos(alpha);
  ExpectMatrixNear(r1, rz_sh->band_rotation(1), kEpsilon);

  // order 2
  Eigen::MatrixXd r2(5, 5);
  r2 << cos(2 * alpha), 0, 0, 0, sin(2 * alpha),
        0, cos(alpha), 0, sin(alpha), 0,
        0, 0, 1, 0, 0,
        0, -sin(alpha), 0, cos(alpha), 0,
        -sin(2 * alpha), 0, 0, 0, cos(2 * alpha);
  ExpectMatrixNear(r2, rz_sh->band_rotation(2), kEpsilon);

  // order 3
  Eigen::MatrixXd r3(7, 7);
  r3 << cos(3 * alpha), 0, 0, 0, 0, 0, sin(3 * alpha),
        0, cos(2 * alpha), 0, 0, 0, sin(2 * alpha), 0,
        0, 0, cos(alpha), 0, sin(alpha), 0, 0,
        0, 0, 0, 1, 0, 0, 0,
        0, 0, -sin(alpha), 0, cos(alpha), 0, 0,
        0, -sin(2 * alpha), 0, 0, 0, cos(2 * alpha), 0,
        -sin(3 * alpha), 0, 0, 0, 0, 0, cos(3 * alpha);
  ExpectMatrixNear(r3, rz_sh->band_rotation(3), kEpsilon);
}

TEST(SphericalHarmonicsRotationTest, ClosedFormBands) {
  // Use an arbitrary rotation
  Eigen::Quaterniond r(Eigen::AngleAxisd(
      0.423, Eigen::Vector3d(0.234, -0.642, 0.829).normalized()));
  Eigen::Matrix3d r_mat = r.toRotationMatrix();

  // Create rotation for band 1 and 2
  std::unique_ptr<Rotation> sh_rot(Rotation::Create(2, r));

  // For l = 1, the transformation matrix for the coefficients is relatively
  // easy to derive. If R is the rotation matrix, the elements of the transform
  // can be described as: Mij = integral_over_sphere Yi(R * s)Yj(s) ds.
  // For l = 1, we have:
  //   Y0(s) = -0.5sqrt(3/pi)s.y
  //   Y1(s) = 0.5sqrt(3/pi)s.z
  //   Y2(s) = -0.5sqrt(3/pi)s.x
  // Note that these Yi include the Condon-Shortely phase. The expectent matrix
  // M is equal to:
  //   [ R11  -R12   R10
  //    -R21   R22  -R20
  //     R01  -R02   R00 ]
  // In [1]'s Appendix summarizing [4], this is given without the negative signs
  // and is a simple permutation, but that is because [4] does not include the
  // Condon-Shortely phase in their definition of the SH basis functions.
  Eigen::Matrix3d band_1 = sh_rot->band_rotation(1);

  EXPECT_DOUBLE_EQ(r_mat(1, 1), band_1(0, 0));
  EXPECT_DOUBLE_EQ(-r_mat(1, 2), band_1(0, 1));
  EXPECT_DOUBLE_EQ(r_mat(1, 0), band_1(0, 2));
  EXPECT_DOUBLE_EQ(-r_mat(2, 1), band_1(1, 0));
  EXPECT_DOUBLE_EQ(r_mat(2, 2), band_1(1, 1));
  EXPECT_DOUBLE_EQ(-r_mat(2, 0), band_1(1, 2));
  EXPECT_DOUBLE_EQ(r_mat(0, 1), band_1(2, 0));
  EXPECT_DOUBLE_EQ(-r_mat(0, 2), band_1(2, 1));
  EXPECT_DOUBLE_EQ(r_mat(0, 0), band_1(2, 2));

  // The l = 2 band transformation is significantly more complex in terms of R,
  // and a CAS program was used to arrive at these equations (plus a fair
  // amount of simplification by hand afterwards).
  Eigen::MatrixXd band_2 = sh_rot->band_rotation(2);

  EXPECT_NEAR(r_mat(0, 0) * r_mat(1, 1) + r_mat(0, 1) * r_mat(1, 0),
              band_2(0, 0), kEpsilon);
  EXPECT_NEAR(-r_mat(0, 1) * r_mat(1, 2) - r_mat(0, 2) * r_mat(1, 1),
              band_2(0, 1), kEpsilon);
  EXPECT_NEAR(-sqrt(3) / 3 * (r_mat(0, 0) * r_mat(1, 0) +
                              r_mat(0, 1) * r_mat(1, 1) -
                              2 * r_mat(0, 2) * r_mat(1, 2)),
              band_2(0, 2), kEpsilon);
  EXPECT_NEAR(-r_mat(0, 0) * r_mat(1, 2) - r_mat(0, 2) * r_mat(1, 0),
              band_2(0, 3), kEpsilon);
  EXPECT_NEAR(r_mat(0, 0) * r_mat(1, 0) - r_mat(0, 1) * r_mat(1, 1),
              band_2(0, 4), kEpsilon);

  EXPECT_NEAR(-r_mat(1, 0) * r_mat(2, 1) - r_mat(1, 1) * r_mat(2, 0),
              band_2(1, 0), kEpsilon);
  EXPECT_NEAR(r_mat(1, 1) * r_mat(2, 2) + r_mat(1, 2) * r_mat(2, 1),
              band_2(1, 1), kEpsilon);
  EXPECT_NEAR(sqrt(3) / 3* (r_mat(1, 0) * r_mat(2, 0) +
                            r_mat(1, 1) * r_mat(2, 1) -
                            2 * r_mat(1, 2) * r_mat(2, 2)),
              band_2(1, 2), kEpsilon);
  EXPECT_NEAR(r_mat(1, 0) * r_mat(2, 2) + r_mat(1, 2) * r_mat(2, 0),
              band_2(1, 3), kEpsilon);
  EXPECT_NEAR(-r_mat(1, 0) * r_mat(2, 0) + r_mat(1, 1) * r_mat(2, 1),
              band_2(1, 4), kEpsilon);

  EXPECT_NEAR(-sqrt(3) / 3 * (r_mat(0, 0) * r_mat(0, 1) +
                              r_mat(1, 0) * r_mat(1, 1) -
                              2 * r_mat(2, 0) * r_mat(2, 1)),
              band_2(2, 0), kEpsilon);
  EXPECT_NEAR(sqrt(3) / 3 * (r_mat(0, 1) * r_mat(0, 2) +
                             r_mat(1, 1) * r_mat(1, 2) -
                             2 * r_mat(2, 1) * r_mat(2, 2)),
              band_2(2, 1), kEpsilon);
  EXPECT_NEAR(-0.5 * (1 - 3 * r_mat(2, 2) * r_mat(2, 2)),
              band_2(2, 2), kEpsilon);
  EXPECT_NEAR(sqrt(3) / 3 * (r_mat(0, 0) * r_mat(0, 2) +
                             r_mat(1, 0) * r_mat(1, 2) -
                             2 * r_mat(2, 0) * r_mat(2, 2)),
              band_2(2, 3), kEpsilon);
  EXPECT_NEAR(sqrt(3) / 6 * (-r_mat(0, 0) * r_mat(0, 0) +
                             r_mat(0, 1) * r_mat(0, 1) -
                             r_mat(1, 0) * r_mat(1, 0) +
                             r_mat(1, 1) * r_mat(1, 1) +
                             2 * r_mat(2, 0) * r_mat(2, 0) -
                             2 * r_mat(2, 1) * r_mat(2, 1)),
              band_2(2, 4), kEpsilon);

  EXPECT_NEAR(-r_mat(0, 0) * r_mat(2, 1) - r_mat(0, 1) * r_mat(2, 0),
              band_2(3, 0), kEpsilon);
  EXPECT_NEAR(r_mat(0, 1) * r_mat(2, 2) + r_mat(0, 2) * r_mat(2, 1),
              band_2(3, 1), kEpsilon);
  EXPECT_NEAR(sqrt(3) / 3 * (r_mat(0, 0) * r_mat(2, 0) +
                             r_mat(0, 1) * r_mat(2, 1) -
                             2 * r_mat(0, 2) * r_mat(2, 2)),
              band_2(3, 2), kEpsilon);
  EXPECT_NEAR(r_mat(0, 0) * r_mat(2, 2) + r_mat(0, 2) * r_mat(2, 0),
              band_2(3, 3), kEpsilon);
  EXPECT_NEAR(-r_mat(0, 0) * r_mat(2, 0) + r_mat(0, 1) * r_mat(2, 1),
              band_2(3, 4), kEpsilon);

  EXPECT_NEAR(r_mat(0, 0) * r_mat(0, 1) - r_mat(1, 0) * r_mat(1, 1),
              band_2(4, 0), kEpsilon);
  EXPECT_NEAR(-r_mat(0, 1) * r_mat(0, 2) + r_mat(1, 1) * r_mat(1, 2),
              band_2(4, 1), kEpsilon);
  EXPECT_NEAR(sqrt(3) / 6 * (-r_mat(0, 0) * r_mat(0, 0) -
                             r_mat(0, 1) * r_mat(0, 1) +
                             r_mat(1, 0) * r_mat(1, 0) +
                             r_mat(1, 1) * r_mat(1, 1) +
                             2 * r_mat(0, 2) * r_mat(0, 2) -
                             2 * r_mat(1, 2) * r_mat(1, 2)),
              band_2(4, 2), kEpsilon);
  EXPECT_NEAR(-r_mat(0, 0) * r_mat(0, 2) + r_mat(1, 0) * r_mat(1, 2),
              band_2(4, 3), kEpsilon);
  EXPECT_NEAR(0.5 * (r_mat(0, 0) * r_mat(0, 0) -
                     r_mat(0, 1) * r_mat(0, 1) -
                     r_mat(1, 0) * r_mat(1, 0) +
                     r_mat(1, 1) * r_mat(1, 1)),
              band_2(4, 4), kEpsilon);
}

TEST(SphericalHarmonicsRotationTest, CreateFromSHRotation) {
  Eigen::Quaterniond r(Eigen::AngleAxisd(M_PI / 4.0, Eigen::Vector3d::UnitY()));
  std::unique_ptr<Rotation> low_band(Rotation::Create(3, r));
  std::unique_ptr<Rotation> high_band(Rotation::Create(5, r));
  std::unique_ptr<Rotation> from_low(Rotation::Create(5, *low_band));

  for (int l = 0; l <= 5; l++) {
    ExpectMatrixNear(high_band->band_rotation(l),
                     from_low->band_rotation(l), kEpsilon);
  }
}

TEST(SphericalHarmonicsRotationTest, RotateSymmetricFunction) {
  SphericalFunction function = [] (double phi, double theta) {
    Eigen::Vector3d d = ToVector(phi, theta);
    return Eigen::Vector3d::UnitZ().dot(d);
  };
  std::unique_ptr<std::vector<double>> coeff = ProjectFunction(
      3, function, kTestSampleCount);

  // Rotation about the z-axis, but the function used is rotationally symmetric
  // about the z-axis so the rotated coefficients should have little change.
  for (double angle = 0.0; angle < 2.0 * M_PI; angle += M_PI / 8) {
    Eigen::Quaterniond r1(Eigen::AngleAxisd(angle, Eigen::Vector3d::UnitZ()));
    std::unique_ptr<Rotation> r1_sh(Rotation::Create(3, r1));

    std::vector<double> r1_coeff;
    r1_sh->Apply(*coeff, &r1_coeff);

    // Compare the rotated coefficients to the coefficients fitted to the
    // rotated source function.
    for (int i = 0; i < 16; i++) {
      // r1 was a rotation about the z-axis, so even though the rotation isn't
      // the identity, the transformed coefficients should be equal to the
      // original coefficients.
      EXPECT_NEAR((*coeff)[i], r1_coeff[i], kCoeffErr);
    }
  }

  // Rotate about more arbitrary angles to test a simple function's rotation.
  Eigen::Vector3d axis = Eigen::Vector3d(0.234, -0.642, 0.829).normalized();
  std::vector<double> rotated_coeff;
  for (double angle = 0.0; angle < 2.0 * M_PI; angle += M_PI / 8) {
    Eigen::Quaterniond rotation(Eigen::AngleAxisd(angle, axis));
    Eigen::Quaterniond r_inv = rotation.inverse();
    std::unique_ptr<Rotation> r_sh(Rotation::Create(3, rotation));

    SphericalFunction rotated_function = [&] (double p, double t) {
      const Eigen::Vector3d n(0, 0, 1);
      Eigen::Vector3d d = r_inv * ToVector(p, t);
      return n.dot(d);
    };

    std::unique_ptr<std::vector<double>> expected_coeff = ProjectFunction(
        3, rotated_function, kTestSampleCount);
    r_sh->Apply(*coeff, &rotated_coeff);

    for (int i = 0; i < 16; i++) {
      EXPECT_NEAR((*expected_coeff)[i], rotated_coeff[i], kCoeffErr);
    }
  }
}

TEST(SphericalHarmonicsRotationTest, RotateComplexFunction) {
  const std::vector<double> coeff = {-1.028, 0.779, -0.275, 0.601, -0.256,
                                     1.891, -1.658, -0.370, -0.772,
                                     -0.591, -0.713, 0.191, 1.206, -0.587,
                                     -0.051, 1.543, -0.370, -0.772,
                                     -0.591, -0.713, 0.191, 1.206, -0.587,
                                     -0.051, 1.543};

  Eigen::Vector3d axis = Eigen::Vector3d(-.43, 0.19, 0.634).normalized();
  std::vector<double> rotated_coeff;
  for (double angle = 0.0; angle < 2.0 * M_PI; angle += M_PI / 8) {
    Eigen::Quaterniond rotation(Eigen::AngleAxisd(angle, axis));
    Eigen::Quaterniond r_inv = rotation.inverse();
    std::unique_ptr<Rotation> r_sh(Rotation::Create(4, rotation));

    SphericalFunction rotated_function = [&] (double p, double t) {
      ToSphericalCoords(r_inv * ToVector(p, t), &p, &t);
      return EvalSHSum(4, coeff, p, t);
    };
    std::unique_ptr<std::vector<double>> expected_coeff = ProjectFunction(
        4, rotated_function, kTestSampleCount);
    r_sh->Apply(coeff, &rotated_coeff);

    for (int i = 0; i < 25; i++) {
      EXPECT_NEAR((*expected_coeff)[i], rotated_coeff[i], kCoeffErr);
    }
  }
}

TEST(SphericalHarmonicsRotationTest, RotateInPlace) {
  // Rotate the function about the y axis by pi/4, which is no longer an
  // identity for the SH coefficients
  const Eigen::Quaterniond r(Eigen::AngleAxisd(
      M_PI / 4.0, Eigen::Vector3d::UnitY()));
  const Eigen::Quaterniond r_inv = r.inverse();
  std::unique_ptr<Rotation> r_sh(Rotation::Create(3, r));

  SphericalFunction function = [] (double phi, double theta) {
    Eigen::Vector3d d = ToVector(phi, theta);
    return Clamp(Eigen::Vector3d::UnitZ().dot(d), 0.0, 1.0);
  };
  SphericalFunction rotated_function = [&] (double phi, double theta) {
    Eigen::Vector3d d = r_inv * ToVector(phi, theta);
    return Clamp(Eigen::Vector3d::UnitZ().dot(d), 0.0, 1.0);
  };

  std::unique_ptr<std::vector<double>> coeff = ProjectFunction(
      3, function, kTestSampleCount);
  std::unique_ptr<std::vector<double>> rotated_coeff = ProjectFunction(
      3, rotated_function, kTestSampleCount);

  r_sh->Apply(*coeff, coeff.get());

  // Compare the rotated coefficients to the coefficients fitted to the
  // rotated source function
  for (int i = 0; i < 16; i++) {
    EXPECT_NEAR((*rotated_coeff)[i], (*coeff)[i], kCoeffErr);
  }
}

TEST(SphericalHarmonicsRotationTest, RotateArray3f) {
  // The coefficients for red, green, and blue channels
  const std::vector<double> c_red = {-1.028, 0.779, -0.275, 0.601, -0.256,
                                     1.891, -1.658, -0.370, -0.772};
  const std::vector<double> c_green = {-0.591, -0.713, 0.191, 1.206, -0.587,
                                       -0.051, 1.543, -0.818, 1.482};
  const std::vector<double> c_blue = {-1.119, 0.559, 0.433, -0.680, -1.815,
                                      -0.915, 1.345, 1.572, -0.622};

  // Combined as an Array3f
  std::vector<Eigen::Array3f> combined;
  for (unsigned int i = 0; i < c_red.size(); i++) {
    combined.push_back(Eigen::Array3f(c_red[i], c_green[i], c_blue[i]));
  }


  // Rotate the function about the y axis by pi/4, which is no longer an
  // identity for the SH coefficients
  const Eigen::Quaterniond r(Eigen::AngleAxisd(
      M_PI / 4.0, Eigen::Vector3d::UnitY()));
  std::unique_ptr<Rotation> r_sh(Rotation::Create(2, r));

  std::vector<double> rotated_r;
  std::vector<double> rotated_g;
  std::vector<double> rotated_b;
  std::vector<Eigen::Array3f> rotated_combined;

  r_sh->Apply(c_red, &rotated_r);
  r_sh->Apply(c_green, &rotated_g);
  r_sh->Apply(c_blue, &rotated_b);
  r_sh->Apply(combined, &rotated_combined);

  for (unsigned int i = 0; i < c_red.size(); i++) {
    EXPECT_FLOAT_EQ(rotated_r[i], rotated_combined[i](0));
    EXPECT_FLOAT_EQ(rotated_g[i], rotated_combined[i](1));
    EXPECT_FLOAT_EQ(rotated_b[i], rotated_combined[i](2));
  }
}

TEST(SphericalHarmonicsRotationTest, RotateArray3fInPlace) {
  // The coefficients for red, green, and blue channels
  const std::vector<double> c_red = {-1.028, 0.779, -0.275, 0.601, -0.256,
                                     1.891, -1.658, -0.370, -0.772};
  const std::vector<double> c_green = {-0.591, -0.713, 0.191, 1.206, -0.587,
                                       -0.051, 1.543, -0.818, 1.482};
  const std::vector<double> c_blue = {-1.119, 0.559, 0.433, -0.680, -1.815,
                                      -0.915, 1.345, 1.572, -0.622};

  // Combined as an Array3f
  std::vector<Eigen::Array3f> combined;
  for (unsigned int i = 0; i < c_red.size(); i++) {
    combined.push_back(Eigen::Array3f(c_red[i], c_green[i], c_blue[i]));
  }


  // Rotate the function about the y axis by pi/4, which is no longer an
  // identity for the SH coefficients
  const Eigen::Quaterniond r(Eigen::AngleAxisd(
      M_PI / 4.0, Eigen::Vector3d::UnitY()));
  std::unique_ptr<Rotation> r_sh(Rotation::Create(2, r));

  std::vector<double> rotated_r;
  std::vector<double> rotated_g;
  std::vector<double> rotated_b;

  r_sh->Apply(c_red, &rotated_r);
  r_sh->Apply(c_green, &rotated_g);
  r_sh->Apply(c_blue, &rotated_b);
  r_sh->Apply(combined, &combined);

  for (unsigned int i = 0; i < c_red.size(); i++) {
    EXPECT_FLOAT_EQ(rotated_r[i], combined[i](0));
    EXPECT_FLOAT_EQ(rotated_g[i], combined[i](1));
    EXPECT_FLOAT_EQ(rotated_b[i], combined[i](2));
  }
}

TEST(SphericalHarmonicsRotationDeathTest, CreateFromMatrixBadInputs) {
  Eigen::Quaterniond good_r(Eigen::AngleAxisd(
      M_PI / 4.0, Eigen::Vector3d::UnitY()));
  Eigen::Quaterniond bad_r(0.0, 0.0, 0.1, 0.3);

  EXPECT_DEATH(Rotation::Create(-1, good_r), "Order must be at least 0.");
  EXPECT_DEATH(Rotation::Create(2, bad_r), "Rotation must be normalized.");
}

TEST(SphericalHarmonicsRotationDeathTest, CreateFromSHBadInputs) {
  Eigen::Quaterniond good_r(Eigen::AngleAxisd(
      M_PI / 4.0, Eigen::Vector3d::UnitY()));
  std::unique_ptr<Rotation> good_sh(Rotation::Create(2, good_r));

  EXPECT_DEATH(Rotation::Create(-1, *good_sh), "Order must be at least 0.");
}

TEST(SphericalHarmonicsRotationDeathTest, RotateBadInputs) {
  Eigen::Quaterniond r(Eigen::AngleAxisd(M_PI / 4.0, Eigen::Vector3d::UnitY()));
  std::unique_ptr<Rotation> sh(Rotation::Create(2, r));

  const std::vector<double> bad_coeffs = {-0.459, 0.3242};

  std::vector<double> result;
  EXPECT_DEATH(sh->Apply(bad_coeffs, &result),
               "Incorrect number of coefficients provided.");
}

TEST(SphericalHarmonicsTest, ImageCoordsToSphericalCoordsTest) {
  // Half-pixel angular increments.
  double pixel_phi = M_PI / kImageWidth;
  double pixel_theta = 0.5 * M_PI / kImageHeight;

  EXPECT_NEAR(pixel_phi,  ImageXToPhi(0, kImageWidth), kEpsilon);
  EXPECT_NEAR(pixel_theta, ImageYToTheta(0, kImageHeight), kEpsilon);

  EXPECT_NEAR(M_PI + pixel_phi, ImageXToPhi(kImageWidth / 2, kImageWidth),
              kEpsilon);
  EXPECT_NEAR(M_PI / 2 + pixel_theta,
              ImageYToTheta(kImageHeight / 2, kImageHeight), kEpsilon);

  EXPECT_NEAR(2 * M_PI - pixel_phi, ImageXToPhi(kImageWidth - 1, kImageWidth),
              kEpsilon);
  EXPECT_NEAR(M_PI - pixel_theta,
              ImageYToTheta(kImageHeight - 1, kImageHeight), kEpsilon);

  // Out of bounds pixels on either side of the image range
  EXPECT_NEAR(-pixel_phi, ImageXToPhi(-1, kImageWidth), kEpsilon);
  EXPECT_NEAR(-pixel_theta, ImageYToTheta(-1, kImageHeight), kEpsilon);

  EXPECT_NEAR(2 * M_PI + pixel_phi, ImageXToPhi(kImageWidth, kImageWidth),
              kEpsilon);
  EXPECT_NEAR(M_PI + pixel_theta, ImageYToTheta(kImageHeight, kImageHeight),
              kEpsilon);
}

TEST(SphericalHarmonicsTest, SphericalCoordsToImageCoordsTest) {
  EXPECT_TUPLE2_NEAR(Eigen::Vector2d(0.0, 0.0),
                     ToImageCoords(0.0, 0.0, kImageWidth, kImageHeight), 
                     kEpsilon);

  EXPECT_TUPLE2_NEAR(Eigen::Vector2d(kImageWidth / 2.0, kImageHeight / 2.0),
                     ToImageCoords(M_PI, M_PI / 2.0, kImageWidth, kImageHeight),
                     kEpsilon);

  EXPECT_TUPLE2_NEAR(Eigen::Vector2d(kImageWidth, kImageHeight),
                     ToImageCoords(2 * M_PI - kEpsilon * 1e-3,
                                   M_PI, kImageWidth, kImageHeight), kEpsilon);

  // Out of the normal phi, theta ranges

  // Negative rotation half a pixel past 0 in the xy plane and a negative
  // rotation half a pixel past 0 in the z axis. The equivalent in-range angles
  // are a half pixel rotation along the z axis and a rotation just shy of
  // 180 in the xy plane (full rotation to bring it into range, and a 180
  // adjust for the z axis).
  EXPECT_TUPLE2_NEAR(Eigen::Vector2d(kImageWidth / 2 - 0.5, 0.5),
                     ToImageCoords(-M_PI / kImageWidth, 
                                   -0.5 * M_PI / kImageHeight,
                                   kImageWidth, kImageHeight), kEpsilon);

  // A half pixel past one full rotation in the xy plane and the through the
  // z axis. The equivalent in-range angles are a half pixel past 180 degrees
  // in the xy plane and half a pixel shy of 180 along the z axis.
  EXPECT_TUPLE2_NEAR(Eigen::Vector2d(kImageWidth / 2 + 0.5, kImageHeight - 0.5),
                     ToImageCoords(
                         (kImageWidth + 0.5) * 2 * M_PI / kImageWidth,
                         (kImageHeight + 0.5) * M_PI / kImageHeight,
                         kImageWidth, kImageHeight), kEpsilon);
}

TEST(SphericalHarmonicsTest, RenderDiffuseIrradianceTest) {
  // Use a simple function where it's convenient to analytically determine the
  // diffuse irradiance. If the light is a constant for the positive z
  // hemisphere, and the query normal is the positive z axis, then the diffuse
  // irradiance is the integral of cos(theta) over that hemisphere, which is
  // equal to:
  //   int(int(light * cos(theta) * sin(theta), 0, 2pi, phi), 0, pi/2, theta) ->
  //   2pi * light * int(cos(theta) * sin(theta), 0, pi/2, theta) ->
  //   2pi * light * (-0.5cos(theta)^2|0,pi/2) ->
  //   pi * light * (-cos(pi/2)+cos(0)) ->
  //   pi * light
  // Set light = 1 for simplicity.

  std::unique_ptr<std::vector<double>> env =
      ProjectFunction(2, [] (double phi, double theta) {
        return theta > M_PI / 2 ? 0.0 : 1.0;
      }, 2500);
  // Convert it to RGB
  std::vector<Eigen::Array3f> env_rgb(9);
  for (int i = 0; i < 9; i++) {
    env_rgb[i] = Eigen::Array3f((*env)[i], (*env)[i], (*env)[i]);
  }

  Eigen::Array3f diffuse_irradiance = RenderDiffuseIrradiance(
      env_rgb, Eigen::Vector3d::UnitZ());
  EXPECT_TUPLE3_NEAR(Eigen::Array3f(M_PI, M_PI, M_PI), diffuse_irradiance, 
                     kIrradianceError);
}

TEST(SphericalHarmonicsTest, RenderDiffuseIrradianceMapTest) {
  DefaultImage env_map(kImageWidth, kImageHeight);
  DefaultImage expected_diffuse(kImageWidth, kImageHeight);

  GenerateTestEnvironment(&env_map);
  ComputeExplicitDiffuseIrradiance(env_map, &expected_diffuse);

  DefaultImage rendered_diffuse(kImageWidth, kImageHeight);
  RenderDiffuseIrradianceMap(env_map, &rendered_diffuse);

  // Compare the two diffuse irradiance maps
  for (int y = 0; y < expected_diffuse.height(); y++) {
    for (int x = 0; x < expected_diffuse.width(); x++) {
      Eigen::Array3f expected = expected_diffuse.GetPixel(x, y);
      Eigen::Array3f rendered = rendered_diffuse.GetPixel(x, y);

      EXPECT_TUPLE3_NEAR(expected, rendered, kEnvMapIrradianceError);
    }
  }
}

TEST(SphericalHarmonicsTest, RenderDiffuseIrradianceMoreCoefficientsTest) {
  // Coefficients of the expected size (9).
  std::vector<Eigen::Array3f> coeffs9;
  for (int i = 0; i < 9; i++) {
    float c = static_cast<float>(i);
    coeffs9.push_back({c, c, c});
  }
  // Coefficients equal to coeffs9 for the first 9 and then 3 additional
  // coefficients that should be ignored.
  std::vector<Eigen::Array3f> coeffs12;
  for (int i = 0; i < 12; i++) {
    float c = static_cast<float>(i);
    coeffs12.push_back({c, c, c});
  }
  EXPECT_TUPLE3_NEAR(
      RenderDiffuseIrradiance(coeffs9, Eigen::Vector3d::UnitZ()),
      RenderDiffuseIrradiance(coeffs12, Eigen::Vector3d::UnitZ()), kEpsilon);
}

TEST(SphericalHarmonicsTest, RenderDiffuseIrradianceFewCoefficientsTest) {
  std::vector<Eigen::Array3f> coeffs3;
  for (int i = 0; i < 3; i++) {
    float c = static_cast<float>(i);
    coeffs3.push_back({c, c, c});
  }
  // Coefficients of the expected size (9), padded with zeros so the first 3
  // are equal to coeffs3.
  std::vector<Eigen::Array3f> coeffs9 = coeffs3;
  for (unsigned int i = coeffs3.size(); i < 9; i++) {
    coeffs9.push_back({0.0, 0.0, 0.0});
  }

  EXPECT_TUPLE3_NEAR(RenderDiffuseIrradiance(coeffs9, Eigen::Vector3d::UnitZ()),
                     RenderDiffuseIrradiance(coeffs3, Eigen::Vector3d::UnitZ()),
                     kEpsilon);
}

TEST(SphericalHarmonicUtilsTest, RenderDiffuseIrradianceNoCoefficientsTest) {
  EXPECT_TUPLE3_NEAR(Eigen::Array3f(0.0, 0.0, 0.0),
                     RenderDiffuseIrradiance({}, Eigen::Vector3d::UnitZ()), 
                     kEpsilon);
}

}  // namespace sh
