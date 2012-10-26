/* opencl_support.h
*/
#ifndef OPENCL_SUPPORT_H_INCLUDED
#define OPENCL_SUPPORT_H_INCLUDED

#include <CL/cl.h>

cl_program compileProgramOpenCL(cl_context context, cl_device_id device, const char *programPath);

class OpenCLConfig
{
	static const cl_uint MAX_PLATFORMS = 16;		//Arbitrary value
	static const cl_uint MAX_DEVICES = 16;			//Arbitrary value

	cl_platform_id m_platformId;
	
public:	
	cl_device_id m_device;
	
	cl_context m_context;
	cl_command_queue m_commandQueue;
	
	OpenCLConfig();
	void initialize();
	void deactivate();

private:
	void initialize_stage1_platform();
	void initialize_stage2_device();
	void initialize_stage3_context_and_queue();
};

#endif