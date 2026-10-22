function report=run_direct_components()
%RUN_DIRECT_COMPONENTS CFL SSPRK3 and direct polar spatial regressions.
time=test_direct_cfl();space=test_spatial_shared();mesh=test_uniform_core_mesh();poisson=test_bspline_poisson();
report=struct('kind','ipm_c_direct_cfl_components_v3','passed',time.passed&&space.passed&&mesh.passed&&poisson.passed, ...
    'time',time,'spatial',space,'mesh',mesh,'poisson',poisson,'implicitSolverUsed',false,'stepRetryUsed',false);
assert(report.passed,'ipm_c:DirectComponents');
end
