/* opencl_support.h
*/
#ifndef OPENCL_SUPPORT_H_INCLUDED
#define OPENCL_SUPPORT_H_INCLUDED

#include <string>

#include <CL/cl.h>


std::string load_text_file(const char *path);

const char* get_cl_device_type(cl_device_type device_type);
const char* get_cl_error_string(cl_int error_code);

#define CHECK_CL_ERROR(ERROR_CODE) check_cl_error(ERROR_CODE, __FILE__, __LINE__);
void check_cl_error(cl_int error, const char *pFile = 0, int line = 0);

cl_program compileProgramOpenCL(cl_context context, cl_device_id device, const char *programPath);

class OpenCLBuffer
{
	unsigned int m_size;
	cl_mem m_buffer;
	
public:	
	OpenCLBuffer() : m_size(0), m_buffer(0) {}

	unsigned int getSize() const { return m_size; }
	cl_mem getBuffer() { return m_buffer; }
	
	void allocate(cl_context context, unsigned int size);
	void deallocate();
	
	void writeToBuffer(cl_command_queue command_queue, const void *source, unsigned int size);
	void readFromBuffer(cl_command_queue command_queue, void *target, unsigned int size);
};

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