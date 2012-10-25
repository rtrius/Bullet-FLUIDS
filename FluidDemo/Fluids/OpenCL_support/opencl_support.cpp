/* opencl_support.cpp
*/

#include "opencl_support.h"

#include <iostream>
#include <fstream>

#include "LinearMath/btAlignedObjectArray.h"

std::string load_text_file(const char *path)
{	
	std::string out = "";

	std::ifstream file(path, std::ios_base::binary);
	if( file.good() )
	{
		std::cout << "Loaded file: " << path << '\n';
		while( file.good() )
		{
			int read_char = file.get();
			if(read_char != EOF) out += static_cast<char>(read_char); 	//Last value returned by file.get() is 'EOF'
		}
	}
	else std::cout << "Failed to load file: " << path << '\n';

	file.close();

	return out;
}

const char* get_cl_device_type(cl_device_type device_type)
{
	switch(device_type)
	{
		case CL_DEVICE_TYPE_CPU:
			return "CL_DEVICE_TYPE_CPU";
		case CL_DEVICE_TYPE_GPU:
			return "CL_DEVICE_TYPE_GPU";
		case CL_DEVICE_TYPE_ACCELERATOR:
			return "CL_DEVICE_TYPE_ACCELERATOR";
		case CL_DEVICE_TYPE_DEFAULT:
			return "CL_DEVICE_TYPE_DEFAULT";
	}
	
	return "Unknown or invalid OpenCL device type.";
}
const char* get_cl_error_string(cl_int error_code)
{
	switch(error_code)
	{
		case CL_SUCCESS:
			return "CL_SUCCESS";
		case CL_DEVICE_NOT_FOUND:
			return "CL_DEVICE_NOT_FOUND";
		case CL_DEVICE_NOT_AVAILABLE:
			return "CL_DEVICE_NOT_AVAILABLE";
		case CL_COMPILER_NOT_AVAILABLE:
			return "CL_COMPILER_NOT_AVAILABLE";
		case CL_MEM_OBJECT_ALLOCATION_FAILURE:
			return "CL_MEM_OBJECT_ALLOCATION_FAILURE";
		case CL_OUT_OF_RESOURCES:
			return "CL_OUT_OF_RESOURCES";
		case CL_OUT_OF_HOST_MEMORY:
			return "CL_OUT_OF_HOST_MEMORY";
		case CL_PROFILING_INFO_NOT_AVAILABLE:
			return "CL_PROFILING_INFO_NOT_AVAILABLE";
		case CL_MEM_COPY_OVERLAP:
			return "CL_MEM_COPY_OVERLAP";
		case CL_IMAGE_FORMAT_MISMATCH:
			return "CL_IMAGE_FORMAT_MISMATCH";
		case CL_IMAGE_FORMAT_NOT_SUPPORTED:
			return "CL_IMAGE_FORMAT_NOT_SUPPORTED";
		case CL_BUILD_PROGRAM_FAILURE:
			return "CL_BUILD_PROGRAM_FAILURE";
		case CL_MAP_FAILURE:
			return "CL_MAP_FAILURE";
		case CL_MISALIGNED_SUB_BUFFER_OFFSET:
			return "CL_MISALIGNED_SUB_BUFFER_OFFSET";
		case CL_EXEC_STATUS_ERROR_FOR_EVENTS_IN_WAIT_LIST:
			return "CL_EXEC_STATUS_ERROR_FOR_EVENTS_IN_WAIT_LIST";

		case CL_INVALID_VALUE:
			return "CL_INVALID_VALUE";
		case CL_INVALID_DEVICE_TYPE:
			return "CL_INVALID_DEVICE_TYPE";
		case CL_INVALID_PLATFORM:
			return "CL_INVALID_PLATFORM";
		case CL_INVALID_DEVICE:
			return "CL_INVALID_DEVICE";
		case CL_INVALID_CONTEXT:
			return "CL_INVALID_CONTEXT";
		case CL_INVALID_QUEUE_PROPERTIES:
			return "CL_INVALID_QUEUE_PROPERTIES";
		case CL_INVALID_COMMAND_QUEUE:
			return "CL_INVALID_COMMAND_QUEUE";
		case CL_INVALID_HOST_PTR:
			return "CL_INVALID_HOST_PTR";
		case CL_INVALID_MEM_OBJECT:
			return "CL_INVALID_MEM_OBJECT";
		case CL_INVALID_IMAGE_FORMAT_DESCRIPTOR:
			return "CL_INVALID_IMAGE_FORMAT_DESCRIPTOR";
		case CL_INVALID_IMAGE_SIZE:
			return "CL_INVALID_IMAGE_SIZE";
		case CL_INVALID_SAMPLER:
			return "CL_INVALID_SAMPLER";
		case CL_INVALID_BINARY:
			return "CL_INVALID_BINARY";
		case CL_INVALID_BUILD_OPTIONS:
			return "CL_INVALID_BUILD_OPTIONS";
		case CL_INVALID_PROGRAM:
			return "CL_INVALID_PROGRAM";
		case CL_INVALID_PROGRAM_EXECUTABLE:
			return "CL_INVALID_PROGRAM_EXECUTABLE";
		case CL_INVALID_KERNEL_NAME:
			return "CL_INVALID_KERNEL_NAME";
		case CL_INVALID_KERNEL_DEFINITION:
			return "CL_INVALID_KERNEL_DEFINITION";
		case CL_INVALID_KERNEL:
			return "CL_INVALID_KERNEL";
		case CL_INVALID_ARG_INDEX:
			return "CL_INVALID_ARG_INDEX";
		case CL_INVALID_ARG_VALUE:
			return "CL_INVALID_ARG_VALUE";
		case CL_INVALID_ARG_SIZE:
			return "CL_INVALID_ARG_SIZE";
		case CL_INVALID_KERNEL_ARGS:
			return "CL_INVALID_KERNEL_ARGS";
		case CL_INVALID_WORK_DIMENSION:
			return "CL_INVALID_WORK_DIMENSION";
		case CL_INVALID_WORK_GROUP_SIZE:
			return "CL_INVALID_WORK_GROUP_SIZE";
		case CL_INVALID_WORK_ITEM_SIZE:
			return "CL_INVALID_WORK_ITEM_SIZE";
		case CL_INVALID_GLOBAL_OFFSET:
			return "CL_INVALID_GLOBAL_OFFSET";
		case CL_INVALID_EVENT_WAIT_LIST:
			return "CL_INVALID_EVENT_WAIT_LIST";
		case CL_INVALID_EVENT:
			return "CL_INVALID_EVENT";
		case CL_INVALID_OPERATION:
			return "CL_INVALID_OPERATION";
		case CL_INVALID_GL_OBJECT:
			return "CL_INVALID_GL_OBJECT";
		case CL_INVALID_BUFFER_SIZE:
			return "CL_INVALID_BUFFER_SIZE";
		case CL_INVALID_MIP_LEVEL:
			return "CL_INVALID_MIP_LEVEL";
		case CL_INVALID_GLOBAL_WORK_SIZE:
			return "CL_INVALID_GLOBAL_WORK_SIZE";
		case CL_INVALID_PROPERTY:
			return "CL_INVALID_PROPERTY";
	}
	
	return "Unknown or invalid OpenCL error code.";
}

void check_cl_error(cl_int error, const char *pFile, int line) 
{ 
	if(error != CL_SUCCESS)
	{
		if(pFile)
			printf("OpenCL error - %s(%d): %s(%d)\n", pFile, line, get_cl_error_string(error), error);
		else
			printf("OpenCL error: %s)%d)\n", get_cl_error_string(error), error);
	}
}

cl_program compileProgramOpenCL(cl_context context, cl_device_id device, const char *programPath)
{
	cl_program program = 0;

	std::string programText = load_text_file(programPath);
	const char *programData = programText.c_str();
	size_t programLength = programText.length();
	
	//Program
	cl_int error_code;
	program = clCreateProgramWithSource(context, 1, const_cast<const char**>(&programData), NULL, &error_code);
	CHECK_CL_ERROR(error_code);

	error_code = clBuildProgram(program, 1, &device, NULL, NULL, NULL);
	if(error_code == CL_SUCCESS) printf("%s compiled successfully.\n", programPath);
	CHECK_CL_ERROR(error_code);
	
	//if(error_code != CL_SUCCESS)
	{
		const int MAX_STRING_LENGTH = 65536;
		btAlignedObjectArray<char> string;
		string.resize(MAX_STRING_LENGTH);
		char *stringStart = &string[0];
	
		error_code = clGetDeviceInfo(device, CL_DEVICE_NAME, MAX_STRING_LENGTH, stringStart, NULL);
		CHECK_CL_ERROR(error_code);
		printf("for CL_DEVICE_NAME: %s\n", stringStart);
		
		error_code = clGetProgramBuildInfo( program, device, CL_PROGRAM_BUILD_LOG, 
											MAX_STRING_LENGTH, stringStart, NULL );
		CHECK_CL_ERROR(error_code);
		
		printf("----------CL Program Build Log - start----------------------\n");
		printf("%s\n", stringStart);
		printf("----------CL Program Build Log - end------------------------\n");
		printf("\n");
	}
	
	
	return program;
}

//#include <CL/cl_ext.h>	//for CL_MEM_USE_PERSISTENT_MEM_AMD -- improves writeToOpenCL() performance
void OpenCLBuffer::allocate(cl_context context, unsigned int size)
{
	cl_int error_code;

	if(!m_buffer)
	{
		//m_buffer = clCreateBuffer(context, CL_MEM_READ_WRITE | CL_MEM_USE_PERSISTENT_MEM_AMD, size, NULL, &error_code);
		m_buffer = clCreateBuffer(context, CL_MEM_READ_WRITE, size, NULL, &error_code);
		CHECK_CL_ERROR(error_code);
		
		m_size = size;
	}
	else
	{
		deallocate();
		allocate(context, size);
	}
}
void OpenCLBuffer::deallocate()
{
	cl_int error_code;

	if(m_buffer)
	{
		error_code = clReleaseMemObject(m_buffer);
		CHECK_CL_ERROR(error_code);
	
		m_buffer = 0;
		m_size = 0;
	}
}

void OpenCLBuffer::writeToBuffer(cl_command_queue command_queue, const void *source, unsigned int size)
{
	btAssert(m_buffer);

	cl_int error_code;

	const cl_bool BLOCK_WRITES = CL_FALSE;
	error_code = clEnqueueWriteBuffer(command_queue, m_buffer, BLOCK_WRITES, 0, size, source, 0, NULL, NULL);
	CHECK_CL_ERROR(error_code);
}

void OpenCLBuffer::readFromBuffer(cl_command_queue command_queue, void *target, unsigned int size)
{	
	btAssert(m_buffer);
	
	cl_int error_code;
	
	const cl_bool BLOCK_READS = CL_FALSE;
	error_code = clEnqueueReadBuffer(command_queue, m_buffer, BLOCK_READS, 0, size, target, 0, NULL, NULL);
	CHECK_CL_ERROR(error_code);
}
	
	
OpenCLConfig::OpenCLConfig()
{
	m_platformId = 0;
	m_context = 0;
	
	m_device = 0;
	m_commandQueue = 0;
}
void OpenCLConfig::initialize()
{
	initialize_stage1_platform();
	initialize_stage2_device();
	initialize_stage3_context_and_queue();
}
void OpenCLConfig::deactivate()
{
	cl_int error_code;
	
	//Command queues
	error_code = clReleaseCommandQueue(m_commandQueue);
	CHECK_CL_ERROR(error_code);
	
	//Context
	error_code = clReleaseContext(m_context);
	CHECK_CL_ERROR(error_code);
	
	//
	m_commandQueue = 0;
	m_context = 0;
}

void OpenCLConfig::initialize_stage1_platform()
{	
	cl_int error_code;
	const size_t MAX_STRING_LENGTH = 1024;
	char string[MAX_STRING_LENGTH];
	
	//Select platform
	cl_uint num_platforms;
	cl_platform_id platforms[MAX_PLATFORMS];
	error_code = clGetPlatformIDs(MAX_PLATFORMS, platforms, &num_platforms);
	CHECK_CL_ERROR(error_code);
	
	printf("Platforms available: \n");
	for(cl_uint i = 0; i < num_platforms; ++i)
	{
		error_code = clGetPlatformInfo(platforms[i], CL_PLATFORM_NAME, MAX_STRING_LENGTH, string, NULL);
		CHECK_CL_ERROR(error_code);
		printf("CL_PLATFORM_NAME: %s\n", string);
		
		error_code = clGetPlatformInfo(platforms[i], CL_PLATFORM_VENDOR, MAX_STRING_LENGTH, string, NULL);
		CHECK_CL_ERROR(error_code);
		printf("CL_PLATFORM_VENDOR: %s\n", string);
		
		error_code = clGetPlatformInfo(platforms[i], CL_PLATFORM_VERSION, MAX_STRING_LENGTH, string, NULL);
		CHECK_CL_ERROR(error_code);
		printf("CL_PLATFORM_VERSION: %s\n", string);
		
		error_code = clGetPlatformInfo(platforms[i], CL_PLATFORM_PROFILE, MAX_STRING_LENGTH, string, NULL);
		CHECK_CL_ERROR(error_code);
		printf("CL_PLATFORM_PROFILE: %s\n", string);
		
		error_code = clGetPlatformInfo(platforms[i], CL_PLATFORM_EXTENSIONS, MAX_STRING_LENGTH, string, NULL);
		CHECK_CL_ERROR(error_code);
		printf("CL_PLATFORM_EXTENSIONS: %s\n", string);
		
		//Select any platform with a GPU device
		cl_uint num_gpu_devices;
		cl_device_id devices[MAX_DEVICES];	//Replacing 'devices' with 'NULL' in clGetDeviceIDs() results in (error_code != CL_SUCCESS)
		error_code = clGetDeviceIDs(platforms[i], CL_DEVICE_TYPE_GPU, MAX_DEVICES, devices, &num_gpu_devices);
		CHECK_CL_ERROR(error_code);	
		
		if( !m_platformId && num_gpu_devices > 0) 
		{
			m_platformId = platforms[i];
			printf("-----Above platform selected.\n");
		}
	}
	printf("\n");
}
void OpenCLConfig::initialize_stage2_device()
{
	cl_int error_code;
	const size_t MAX_STRING_LENGTH = 1024;
	char string[MAX_STRING_LENGTH];
	
	//Select device
	cl_uint num_devices;
	cl_device_id devices[MAX_DEVICES];
	if(m_platformId)
	{
		//Get devices
		error_code = clGetDeviceIDs(m_platformId, CL_DEVICE_TYPE_GPU, MAX_DEVICES, devices, &num_devices);
		CHECK_CL_ERROR(error_code);	
		
		printf("num_devices: %d\n", num_devices);
		for(cl_uint i = 0; i < num_devices; ++i)
		{
			error_code = clGetDeviceInfo(devices[i], CL_DEVICE_NAME, MAX_STRING_LENGTH, string, NULL);
			CHECK_CL_ERROR(error_code);
			printf("CL_DEVICE_NAME: %s\n", string);
			
			error_code = clGetDeviceInfo(devices[i], CL_DEVICE_PLATFORM, MAX_STRING_LENGTH, string, NULL);
			CHECK_CL_ERROR(error_code);
			printf( "CL_DEVICE_PLATFORM: %d\n", *reinterpret_cast<const int*>(string) );

			error_code = clGetDeviceInfo(devices[i], CL_DEVICE_VERSION, MAX_STRING_LENGTH, string, NULL);
			CHECK_CL_ERROR(error_code);
			printf("CL_DEVICE_VERSION: %s\n", string);	
			
			error_code = clGetDeviceInfo(devices[i], CL_DEVICE_OPENCL_C_VERSION, MAX_STRING_LENGTH, string, NULL);
			CHECK_CL_ERROR(error_code);
			printf("CL_DEVICE_OPENCL_C_VERSION: %s\n", string);	
			
			error_code = clGetDeviceInfo(devices[i], CL_DRIVER_VERSION, MAX_STRING_LENGTH, string, NULL);
			CHECK_CL_ERROR(error_code);
			printf("CL_DRIVER_VERSION: %s\n", string);	
				
			error_code = clGetDeviceInfo(devices[i], CL_DEVICE_PROFILE, MAX_STRING_LENGTH, string, NULL);
			CHECK_CL_ERROR(error_code);
			printf("CL_DEVICE_PROFILE: %s\n", string);
			
			error_code = clGetDeviceInfo(devices[i], CL_DEVICE_AVAILABLE, MAX_STRING_LENGTH, string, NULL);
			CHECK_CL_ERROR(error_code);
			printf( "CL_DEVICE_AVAILABLE: %d\n", *reinterpret_cast<const cl_bool*>(string) );	
			
			error_code = clGetDeviceInfo(devices[i], CL_DEVICE_COMPILER_AVAILABLE, MAX_STRING_LENGTH, string, NULL);
			CHECK_CL_ERROR(error_code);
			printf( "CL_DEVICE_COMPILER_AVAILABLE: %d\n", *reinterpret_cast<const cl_bool*>(string) );
			
			error_code = clGetDeviceInfo(devices[i], CL_DEVICE_ENDIAN_LITTLE, MAX_STRING_LENGTH, string, NULL);
			CHECK_CL_ERROR(error_code);
			printf( "CL_DEVICE_ENDIAN_LITTLE: %d\n", *reinterpret_cast<const cl_bool*>(string) );
			if(*reinterpret_cast<const cl_bool*>(string) != CL_TRUE) printf(" Warning: device does not use little endian encoding.\n");
			
			error_code = clGetDeviceInfo(devices[i], CL_DEVICE_TYPE, MAX_STRING_LENGTH, string, NULL);
			CHECK_CL_ERROR(error_code);
			printf( "CL_DEVICE_TYPE: %s\n", get_cl_device_type(*reinterpret_cast<cl_device_type*>(string)) );
			
			//Select the first available GPU device
			if(!m_device)
			{
				m_device = devices[i];
				printf("-----Above device selected.\n");
			}
			
			printf("\n");
		}
		printf("\n");
	}
	else
	{
		printf("initialize_stage2_device() error: invalid m_platformId. \n");
	}
}
void OpenCLConfig::initialize_stage3_context_and_queue()
{
	cl_int error_code;
	
	if(m_device)
	{
		//Create context
		cl_context_properties context_properties[3] = { CL_CONTEXT_PLATFORM, cl_context_properties(m_platformId), 0 };
		m_context = clCreateContext(context_properties, 1, &m_device, NULL, NULL, &error_code);
		CHECK_CL_ERROR(error_code);	
		
		//Create command queue
		const cl_command_queue_properties COMMAND_QUEUE_PROPERTIES = 0;		//CL_QUEUE_OUT_OF_ORDER_EXEC_MODE_ENABLE;
		m_commandQueue = clCreateCommandQueue(m_context, m_device, COMMAND_QUEUE_PROPERTIES, &error_code);
		CHECK_CL_ERROR(error_code);
	}
	else
	{
		printf("initialize_stage3_context_and_queue() error: invalid m_device. \n");
	}
}
