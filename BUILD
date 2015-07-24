cc_library(
    name = "image",
    hdrs = ["image.h"],
    deps = [
        "//external:eigen3",
    ],
    visibility = ["//visibility:public"],
    alwayslink = 1
)

cc_library(
    name = "spherical_harmonics",
    srcs = ["spherical_harmonics.cc"],
    hdrs = ["spherical_harmonics.h"],
    deps = [
        ":image",
        "//external:eigen3",
    ],
    visibility = ["//visibility:public"],
)

cc_library(
    name = "default_image",
    srcs = ["default_image.cc"],
    hdrs = ["default_image.h"],
    deps = [
        ":image",
    ]
)

cc_test(
    name = "spherical_harmonics_test",
    srcs = [
        "spherical_harmonics_test.cc",
    ],
    deps = [
        ":default_image",
        ":spherical_harmonics",
        "//external:gtest",
    ],
    linkopts = ["-lm"],
)

