/* opencl_support.cpp
*/

#include "opencl_support.h"


#include <iostream>
#include <fstream>
#include <string>

#include "LinearMath/btAlignedObjectArray.h"
#include "btExperimentsOpenCL/btOpenCLInclude.h"

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

cl_program compileProgramOpenCL(cl_context context, cl_device_id device, const char *programPath)
{
	cl_program program = 0;

	std::string programText = load_text_file(programPath);
	const char *programData = programText.c_str();
	size_t programLength = programText.length();
	
	//Program
	cl_int error_code;
	program = clCreateProgramWithSource(context, 1, const_cast<const char**>(&programData), NULL, &error_code);
	oclCHECKERROR(error_code, CL_SUCCESS);

	error_code = clBuildProgram(program, 1, &device, NULL, NULL, NULL);
	if(error_code == CL_SUCCESS) printf("%s compiled successfully.\n", programPath);
	oclCHECKERROR(error_code, CL_SUCCESS);
	
	//if(error_code != CL_SUCCESS)
	{
		const int MAX_STRING_LENGTH = 65536;
		btAlignedObjectArray<char> string;
		string.resize(MAX_STRING_LENGTH);
		char *stringStart = &string[0];
	
		error_code = clGetDeviceInfo(device, CL_DEVICE_NAME, MAX_STRING_LENGTH, stringStart, NULL);
		oclCHECKERROR(error_code, CL_SUCCESS);
		printf("for CL_DEVICE_NAME: %s\n", stringStart);
		
		error_code = clGetProgramBuildInfo( program, device, CL_PROGRAM_BUILD_LOG, 
											MAX_STRING_LENGTH, stringStart, NULL );
		oclCHECKERROR(error_code, CL_SUCCESS);
		
		printf("----------CL Program Build Log - start----------------------\n");
		printf("%s\n", stringStart);
		printf("----------CL Program Build Log - end------------------------\n");
		printf("\n");
	}
	
	
	return program;
}
