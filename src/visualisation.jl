"""
Placeholder plotting hooks.  Keep plotting separate from solvers so that headless
CI or benchmarking runs do not pull in heavy graphical dependencies.
"""
function plot_field end

function plot_residuals end
