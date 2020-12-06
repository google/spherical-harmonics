config_setting(
    name = "windows",
    constraint_values = [
        "@bazel_tools//platforms:windows",
    ],
)

cc_library(
    name = "image",
    hdrs = ["image.h"],
    deps = [
        "@eigen3//:eigen3",
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
        "@eigen3//:eigen3",
    ],
    defines = select({
        ":windows": ["_USE_MATH_DEFINES",],
        "//conditions:default": [],
    }),
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
    size = "small",
    srcs = [
        "spherical_harmonics_test.cc",
    ],
    deps = [
        ":default_image",
        ":spherical_harmonics",
        "@gtest//:main",
    ],
    linkopts = ["-lm"],
)

