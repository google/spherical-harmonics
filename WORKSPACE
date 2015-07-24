new_http_archive(
  name = "http-eigen3",
  url = "https://bitbucket.org/eigen/eigen/get/3.2.5.zip",
  sha256 = "73eeb8231acf1b1c498ca0f9d3b94eeb25f2594c0fc3f11c4bb823d019ec0c01",
  build_file = "third_party/eigen3.BUILD",
)

new_http_archive(
  name = "http-gtest",
  url = "https://googletest.googlecode.com/files/gtest-1.7.0.zip",
  sha256 = "247ca18dd83f53deb1328be17e4b1be31514cedfc1e3424f672bf11fd7e0d60d",
  build_file = "third_party/gtest.BUILD",
)

bind(
  name = "eigen3",
  actual = "@http-eigen3//:eigen3",
)

bind(
  name = "gtest",
  actual = "@http-gtest//:main",
)
