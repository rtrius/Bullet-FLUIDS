
@echo off

premake4 --file=stringifyKernel.lua --kernelfile="../Demos/FluidSphDemo/BulletMultiThreaded/SphSolverOpenCL/fluidSph.cl" --headerfile="../Demos/FluidSphDemo/BulletMultiThreaded/SphSolverOpenCL/fluidSphCL.h" --stringname="fluidSphCL" stringify

pause