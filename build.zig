const std = @import("std");

pub fn build(b: *std.Build) void {
    const optimize = b.standardOptimizeOption(.{});
    const target = b.standardTargetOptions(.{});

    // Define the nlopt library
    const lib = b.addLibrary(.{
        .name = "nlopt",
        .linkage = .static, // Default to static, will be conditional later
        .root_module = b.createModule(.{
            .target = target,
            .optimize = optimize,
        }),
    });

    // 2. Configuration Header Generation (nlopt_config.h)
    // Adapted from ik/build.zig's generated_config_h
    const generated_nlopt_config_h = b.addWriteFiles();
    _ = generated_nlopt_config_h.add("nlopt_config.h",
        \\/* Bugfix version number. */
        \\#define BUGFIX_VERSION 0 
        \\
        \\/* Define to enable extra debugging code. */
        \\#undef DEBUG
        \\#undef HAVE_BSDGETTIMEOFDAY
        \\
        \\#define HAVE_COPYSIGN // Assuming true for modern systems, will be conditional later
        \\#define HAVE_FPCLASSIFY // Assuming true for modern systems, will be conditional later
        \\#define HAVE_GETOPT_H // Assuming true for modern systems, will be conditional later
        \\
        \\#define HAVE_GETOPT // Assuming true for modern systems, will be conditional later
        \\#define HAVE_GETPID // Assuming true for modern systems, will be conditional later
        \\#undef HAVE_GETTID_SYSCALL
        \\#define HAVE_GETTIMEOFDAY // Assuming true for modern systems, will be conditional later
        \\#define HAVE_ISINF // Assuming true for modern systems, will be conditional later
        \\#define HAVE_ISNAN // Assuming true for modern systems, will be conditional later
        \\#define HAVE_LIBM // Assuming true for modern systems, will be conditional later
        \\#define HAVE_QSORT_R // Assuming true for modern systems, will be conditional later
        \\#define HAVE_STDINT_H // Assuming true for modern systems, will be conditional later
        \\#define HAVE_SYS_TIME_H // Assuming true for modern systems, will be conditional later
        \\#define HAVE_TIME // Assuming true for modern systems, will be conditional later
        \\#define HAVE_UINT32_T // Assuming true for modern systems, will be conditional later
        \\#define HAVE_UNISTD_H // Assuming true for modern systems, will be conditional later
//        \\#define HAVE_STRNLEN 1
        \\#undef LT_OBJDIR
        \\#define MAJOR_VERSION 2
        \\#define MINOR_VERSION 10
        \\#undef PACKAGE
        \\#undef PACKAGE_BUGREPORT
        \\#undef PACKAGE_NAME
        \\#undef PACKAGE_STRING
        \\#undef PACKAGE_TARNAME
        \\#undef PACKAGE_URL
        \\#undef PACKAGE_VERSION
        \\#undef REPLACEMENT_HUGE_VAL
        \\#define SIZEOF_UNSIGNED_INT 4 // Placeholder, will be determined by Zig
        \\#define SIZEOF_UNSIGNED_LONG 8 // Placeholder, will be determined by Zig
        \\#undef STDC_HEADERS
        \\#define THREADLOCAL // Placeholder, will be determined by Zig
        \\#define TIME_WITH_SYS_TIME // Assuming true for modern systems, will be conditional later
        \\#undef VERSION
        \\#define NLOPT_CXX // Assuming true for now, will be conditional later
        \\#undef const
        \\#ifndef __cplusplus
        \\#undef inline
        \\#endif
        \\
    );
    lib.addIncludePath(generated_nlopt_config_h.getDirectory());

    // 3. Source Files (NLOPT_SOURCES) - Initial C sources
    // This will be a large list, starting with the C files.
    lib.root_module.addCSourceFiles(.{
        .files = &.{ 
            "src/algs/direct/DIRect.c", "src/algs/direct/direct_wrap.c", "src/algs/direct/DIRserial.c", "src/algs/direct/DIRsubrout.c",
            "src/algs/cdirect/cdirect.c", "src/algs/cdirect/hybrid.c",
            "src/algs/praxis/praxis.c",
            "src/algs/crs/crs.c",
            "src/algs/mlsl/mlsl.c",
            "src/algs/mma/mma.c", "src/algs/mma/ccsa_quadratic.c",
            "src/algs/cobyla/cobyla.c",
            "src/algs/newuoa/newuoa.c",
            "src/algs/neldermead/nldrmd.c", "src/algs/neldermead/sbplx.c",
            "src/algs/auglag/auglag.c",
            "src/algs/bobyqa/bobyqa.c",
            "src/algs/isres/isres.c",
            "src/algs/slsqp/slsqp.c",
            "src/algs/esch/esch.c",
            "src/api/general.c", "src/api/options.c", "src/api/optimize.c", "src/api/deprecated.c", "src/api/f77api.c",
            "src/util/mt19937ar.c", "src/util/sobolseq.c", "src/util/timer.c", "src/util/stop.c", "src/util/redblack.c", "src/util/qsort_r.c", "src/util/rescale.c",
        },
        .flags = &.{ 
            "-std=c11",
            "-D_POSIX_C_SOURCE=200809L",
            // Add other common C flags here as needed, from CMakeLists.txt or ik/build.zig
        },
    });

    lib.root_module.addCSourceFiles(.{
        .files = &.{ 
            "src/algs/ags/ags.cc", "src/algs/ags/evolvent.cc", "src/algs/ags/local_optimizer.cc", "src/algs/ags/solver.cc",
            "src/algs/stogo/global.cc", "src/algs/stogo/linalg.cc", "src/algs/stogo/local.cc", "src/algs/stogo/stogo.cc", "src/algs/stogo/tools.cc",
        },
        .flags = &.{ 
            "-std=c++11",
        },
    });

    // 4. Include Paths - Initial paths
    lib.addIncludePath(b.path("src/api")); // For nlopt.h
    lib.addIncludePath(b.path("src/util")); // For nlopt-util.h
    lib.addIncludePath(b.path("src/algs/direct"));
    lib.addIncludePath(b.path("src/algs/cdirect"));
    lib.addIncludePath(b.path("src/algs/praxis"));
    lib.addIncludePath(b.path("src/algs/crs"));
    lib.addIncludePath(b.path("src/algs/mlsl"));
    lib.addIncludePath(b.path("src/algs/mma"));
    lib.addIncludePath(b.path("src/algs/cobyla"));
    lib.addIncludePath(b.path("src/algs/newuoa"));
    lib.addIncludePath(b.path("src/algs/neldermead"));
    lib.addIncludePath(b.path("src/algs/auglag"));
    lib.addIncludePath(b.path("src/algs/bobyqa"));
    lib.addIncludePath(b.path("src/algs/isres"));
    lib.addIncludePath(b.path("src/algs/slsqp"));
    lib.addIncludePath(b.path("src/algs/esch"));
    lib.addIncludePath(b.path("src/algs/stogo"));
    lib.addIncludePath(b.path("src/algs/ags"));
    // 6. Linker Libraries
    lib.linkLibC(); // Link against the C standard library
    lib.linkLibCpp(); // Link against the C++ standard library

    // Install the library artifact
    b.installArtifact(lib);

    // --- Create and export the lib_zik Zig module (C bindings) ---
    const lib_znlopt_mit_module = b.createModule(.{
        .root_source_file = b.path("zig/c.zig"),
        .target = target,
        .optimize = optimize,
    });

    lib_znlopt_mit_module.addIncludePath(b.path("src/api"));

    lib_znlopt_mit_module.linkLibrary(lib);

    b.modules.put("lib_znlopt_mit", lib_znlopt_mit_module) catch @panic("failed to register lib_zik module");

    // --- Build an example executable --- 
    const example_exe = b.addExecutable(.{
        .name = "znlopt_example",
        .root_module = b.createModule(.{
            .root_source_file = b.path("zig/example.zig"),
            .target = target,
            .optimize = optimize,
            .imports = &.{ 
                .{ .name = "lib_znlopt_mit", .module = lib_znlopt_mit_module },
            },
        }),
    });

    b.installArtifact(example_exe);

    const run_example_step = b.addRunArtifact(example_exe);
    const run_step = b.step("run", "Run the nlopt example");
    run_step.dependOn(&run_example_step.step);
}

    // Placeholder for Luksan sources (conditional NLOPT_LUKSAN)
    // TODO: Implement conditional addition of Luksan sources and flags

