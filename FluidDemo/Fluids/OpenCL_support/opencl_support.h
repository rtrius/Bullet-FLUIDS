/* opencl_support.h
*/
#ifndef OPENCL_SUPPORT_H_INCLUDED
#define OPENCL_SUPPORT_H_INCLUDED

#include <string>

#include <CL/cl.h>

// from <CL/cl.h> (OpenCL 1.1):
// 	typedef struct _cl_platform_id *    cl_platform_id;
// 	typedef struct _cl_device_id *      cl_device_id;
// 	typedef struct _cl_context *        cl_context;
// 	typedef struct _cl_command_queue *  cl_command_queue;
// 	typedef struct _cl_mem *            cl_mem;
// 	typedef struct _cl_program *        cl_program;
// 	typedef struct _cl_kernel *         cl_kernel;
// 	typedef struct _cl_event *          cl_event;
// 	typedef struct _cl_sampler *        cl_sampler;
//
// --Since the above typedefs are all pointers,
// it is likely that 0 may be used to represent an invalid value.
const cl_platform_id INVALID_PLATFORM_ID = 0;	
const cl_device_id INVALID_DEVICE_ID = 0;	
const cl_context INVALID_CONTEXT = 0;	
const cl_command_queue INVALID_COMMAND_QUEUE = 0;	
const cl_program INVALID_PROGRAM = 0;	
const cl_kernel INVALID_KERNEL = 0;	
const cl_mem INVALID_BUFFER = 0;




std::string load_text_file(const char *path);

const char* get_cl_device_type(cl_device_type device_type);
const char* get_cl_error_string(cl_int error_code);

#define CHECK_CL_ERROR(ERROR_CODE) check_cl_error(ERROR_CODE, __FILE__, __LINE__);
void check_cl_error(cl_int error, const char *pFile = 0, int line = 0);


class OpenCLBuffer
{
	unsigned int m_size;
	cl_mem m_buffer;
	
public:	
	OpenCLBuffer() : m_size(0), m_buffer(INVALID_BUFFER) {}

	unsigned int getSize() const { return m_size; }
	void* getAddress() { return static_cast<void*>(&m_buffer); }
	
	void allocate(cl_context context, unsigned int size);
	void deallocate();
	
	void writeToBuffer(cl_command_queue command_queue, const void *source, unsigned int size);
	void readFromBuffer(cl_command_queue command_queue, void *target, unsigned int size);
};

#endif