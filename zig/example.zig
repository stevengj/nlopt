const std = @import("std");
const nlopt = @import("lib_znlopt_mit");

const MyFuncData = struct {
    a: f64,
    b: f64,
};

// Objective function: f(x) = (a - x[0])^2 + b * (x[1] - x[0]^2)^2
// COBYLA is derivative-free, so we don't provide a gradient.
fn rosenbrock_df(n: u32, x: [*c]const f64, grad: [*c]f64, func_data: ?*anyopaque) callconv(.c) f64 {
    _ = n;
    _ = grad; // COBYLA is derivative-free, so grad is ignored
    const data: *MyFuncData = @ptrCast(@alignCast(func_data.?));

    const a = data.a;
    const b = data.b;
    const t1 = a - x[0];
    const t2 = x[1] - x[0] * x[0];
    return t1 * t1 + b * t2 * t2;
}

pub fn main() !void {
    const n: u32 = 2; // Dimensionality of the problem
    var x: [2]f64 = .{ 0.0, 0.0 }; // Initial guess
    var min_f: f64 = 0.0; // Minimum objective value

    var data = MyFuncData{ .a = 1.0, .b = 100.0 };

    // Create an NLopt optimizer using COBYLA (Local, Derivative-Free)
    const opt = nlopt.nlopt.nlopt_create(nlopt.nlopt.NLOPT_LN_COBYLA, n);
    if (opt == null) {
        std.debug.print("NLopt failed to create optimizer\n", .{});
        return error.NloptError;
    }

    const set_obj_res = nlopt.nlopt.nlopt_set_min_objective(opt, rosenbrock_df, &data);
    if (set_obj_res != nlopt.nlopt.NLOPT_SUCCESS) {
        std.debug.print("NLopt failed to set objective: {}\n", .{set_obj_res});
        nlopt.nlopt.nlopt_destroy(opt);
        return error.NloptError;
    }

    var lb = [_]f64{ -1.5, -1.5 };
    const set_lb_res = nlopt.nlopt.nlopt_set_lower_bounds(opt, &lb[0]);
    if (set_lb_res != nlopt.nlopt.NLOPT_SUCCESS) {
        std.debug.print("NLopt failed to set lower bounds: {}\n", .{set_lb_res});
        nlopt.nlopt.nlopt_destroy(opt);
        return error.NloptError;
    }

    var ub = [_]f64{ 1.5, 1.5 };
    const set_ub_res = nlopt.nlopt.nlopt_set_upper_bounds(opt, &ub[0]);
    if (set_ub_res != nlopt.nlopt.NLOPT_SUCCESS) {
        std.debug.print("NLopt failed to set upper bounds: {}\n", .{set_ub_res});
        nlopt.nlopt.nlopt_destroy(opt);
        return error.NloptError;
    }

    // Set stopping criteria
    const set_ftol_res = nlopt.nlopt.nlopt_set_ftol_abs(opt, 1e-6);
    if (set_ftol_res != nlopt.nlopt.NLOPT_SUCCESS) {
        std.debug.print("NLopt failed to set ftol_abs: {}\n", .{set_ftol_res});
        nlopt.nlopt.nlopt_destroy(opt);
        return error.NloptError;
    }

    const result = nlopt.nlopt.nlopt_optimize(opt, &x[0], &min_f);

    if (result < 0) {
        std.debug.print("NLopt optimization failed with result: {}\n", .{result});
    } else {
        std.debug.print("NLopt optimization successful!\n", .{});
        std.debug.print("found minimum at f({d}, {d}) = {d}\n", .{ x[0], x[1], min_f });
    }

    nlopt.nlopt.nlopt_destroy(opt);
}

const NloptError = error{};
