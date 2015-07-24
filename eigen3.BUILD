
licenses(["restricted"])  # MPL2, portions GPL v3, LGPL v3, BSD-like

exports_files(["LICENSE"])

package(
    default_hdrs_check = "strict",
)

# All Eigen files, including those under a more restrictive license.
eigen_header_files = glob(
    "Eigen/**",
    exclude = ["Eigen/**/CMakeLists.txt"],
)

cc_library(
    name = "eigen3",
    hdrs = eigen_header_files,
    defines = ["EIGEN_MPL2_ONLY"],
    includes = ["."],
    licenses = ["notice"],
    visibility = ["//visibility:public"],
    alwayslink = 1,
)

