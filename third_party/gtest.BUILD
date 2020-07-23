cc_library(
    name = "main",
    srcs = glob(
        include = ["googletest/src/*.h",
                   "googletest/include/gtest/internal/**/*.h",
                   "googletest/src/*.cc"],
        exclude = ["googletest/src/gtest-all.cc"]
    ),
    hdrs = glob(["googletest/include/gtest/*.h"]),
    includes = ["googletest/",
                "googletest/include"],
    linkopts = ["-pthread"],
    visibility = ["//visibility:public"],
)
