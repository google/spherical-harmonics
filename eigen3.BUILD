licenses(["restricted"])  # MPL2, portions GPL v3, LGPL v3, BSD-like

exports_files(["LICENSE"])

# Name of the eigen3 source after being unpacked from its zip
eigen_project_dir = "eigen-eigen-bdd17ee3b1b3/"

# All Eigen files, including those under a more restrictive license.
eigen_header_files = glob(
    include = [eigen_project_dir + "Eigen/**"],
    exclude = [eigen_project_dir + "Eigen/**/CMakeLists.txt"],
)

cc_library(
    name = "eigen3",
    hdrs = eigen_header_files,
    defines = ["EIGEN_MPL2_ONLY"],
    includes = [eigen_project_dir],
    licenses = ["notice"],
    visibility = ["//visibility:public"],
    alwayslink = 1,
)

