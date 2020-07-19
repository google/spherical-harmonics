package(
    default_visibility = ["//visibility:public"],
)

cc_library(
    name = "eigen3",
    hdrs = glob(
        include = ["Eigen/**"],
        exclude = ["Eigen/**/CMakeLists.txt"],
    ),
    defines = ["EIGEN_MPL2_ONLY", "EIGEN_NO_DEBUG"],
    alwayslink = 1,
)