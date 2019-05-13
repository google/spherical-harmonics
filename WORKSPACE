load("@bazel_tools//tools/build_defs/repo:http.bzl", "http_archive", "http_file")

http_archive(
  name = "http_eigen3",
  url = "https://bitbucket.org/eigen/eigen/get/3.3.5.zip",
  sha256 = "35fa84bc23114b9d37c4597745f8b4e03354a5077579fdba597019f595a602b6",
  build_file = "third_party/eigen3.BUILD",
)

http_archive(
  name = "http_gtest",
  url = "https://github.com/google/googletest/archive/release-1.8.0.zip",
  sha256 = "f3ed3b58511efd272eb074a3a6d6fb79d7c2e6a0e374323d1e6bcbcc1ef141bf",
  build_file = "third_party/gtest.BUILD",
)

bind(
  name = "eigen3",
  actual = "@http_eigen3//:eigen3",
)

bind(
  name = "gtest",
  actual = "@http_gtest//:main",
)
