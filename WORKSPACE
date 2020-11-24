workspace(name = "spherical_harmonics")

load("@bazel_tools//tools/build_defs/repo:http.bzl", "http_archive")

http_archive(
  name = "eigen3",
  url = "https://gitlab.com/libeigen/eigen/-/archive/3.3.5/eigen-3.3.5.zip",
  strip_prefix = "eigen-3.3.5",
  sha256 = "0e7aeece6c8874146c2a4addc437eebdf1ec4026680270f00e76705c8186f0b5",
  build_file = "@//third_party:eigen3.BUILD",
)

http_archive(
  name = "gtest",
  url = "https://github.com/google/googletest/archive/release-1.8.0.zip",
  strip_prefix = "googletest-release-1.8.0",
  sha256 = "f3ed3b58511efd272eb074a3a6d6fb79d7c2e6a0e374323d1e6bcbcc1ef141bf",
  build_file = "@//third_party:gtest.BUILD",
)
