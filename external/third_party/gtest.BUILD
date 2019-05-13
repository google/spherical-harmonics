# Base directory for project after unpacking
googletest_project_dir = "googletest-release-1.8.0/googletest"

cc_library(
    name = "main",
    srcs = glob(
        [googletest_project_dir + "/src/*.h",
         googletest_project_dir + "/include/gtest/internal/**/*.h",
         googletest_project_dir + "/src/*.cc"],
        exclude = [googletest_project_dir + "/src/gtest-all.cc"]
    ),
    hdrs = glob([googletest_project_dir + "/include/gtest/*.h"]),
    includes = [
        googletest_project_dir + "/",
        googletest_project_dir + "/include"
    ],
    linkopts = ["-pthread"],
    visibility = ["//visibility:public"],
)
