licenses(["notice"])

cc_library(
    name = "TreeDecomp",
    hdrs = ["TreeDecomp.h"],
    deps = [
        "//benchmarks/KCore/JulienneDBS17:KCore",
        "//gbbs",
        "//TreeDecomp/common:common"
    ],
)

cc_binary(
    name = "TreeDecomp_main",
    srcs = ["TreeDecomp.cc"],
    deps = [":TreeDecomp"],
)

package(
    default_visibility = ["//visibility:public"],
)
