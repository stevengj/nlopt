const std = @import("std");
pub const nlopt = @cImport({
    @cInclude("nlopt.h");
});