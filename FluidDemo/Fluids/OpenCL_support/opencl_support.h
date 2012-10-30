/* opencl_support.h
*/
#ifndef OPENCL_SUPPORT_H_INCLUDED
#define OPENCL_SUPPORT_H_INCLUDED

#include <CL/cl.h>


cl_program compileProgramOpenCL(cl_context context, cl_device_id device, const char* programPath);



#endif