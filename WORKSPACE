workspace(name = "spherical_harmonics")

load("@bazel_tools//tools/build_defs/repo:http.bzl", "http_archive")

http_archive(
  name = "eigen3",
  url = "https://bitbucket.org/eigen/eigen/get/3.3.5.zip",
  strip_prefix = "eigen-eigen-b3f3d4950030",
  sha256 = "35fa84bc23114b9d37c4597745f8b4e03354a5077579fdba597019f595a602b6",
  build_file = "@//third_party:eigen3.BUILD",
)

http_archive(
  name = "gtest",
  url = "https://github.com/google/googletest/archive/release-1.8.0.zip",
  strip_prefix = "googletest-release-1.8.0",
  sha256 = "f3ed3b58511efd272eb074a3a6d6fb79d7c2e6a0e374323d1e6bcbcc1ef141bf",
  build_file = "@//third_party:gtest.BUILD",
)
