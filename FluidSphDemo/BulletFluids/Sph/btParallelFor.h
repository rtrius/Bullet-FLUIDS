/*
Bullet-FLUIDS 
Copyright (c) 2012 Jackson Lee

This software is provided 'as-is', without any express or implied warranty.
In no event will the authors be held liable for any damages arising from the use of this software.
Permission is granted to anyone to use this software for any purpose, 
including commercial applications, and to alter it and redistribute it freely, 
subject to the following restrictions:

1. The origin of this software must not be misrepresented; you must not claim that you wrote the original software. 
   If you use this software in a product, an acknowledgment in the product documentation would be appreciated but is not required.
2. Altered source versions must be plainly marked as such, and must not be misrepresented as being the original software.
3. This notice may not be removed or altered from any source distribution.
*/
#ifndef BT_PARALLEL_FOR_H
#define BT_PARALLEL_FOR_H

#include "BulletMultiThreaded/SpuCollisionTaskProcess.h"

#ifdef _WIN32
	#include "BulletMultiThreaded/Win32ThreadSupport.h"
	typedef Win32ThreadFunc ThreadFunc;
	typedef Win32lsMemorySetupFunc lsMemorySetupFunc;
#elif defined (USE_PTHREADS)
	#include "BulletMultiThreaded/PosixThreadSupport.h"
	typedef PosixThreadFunc ThreadFunc;
	typedef PosixlsMemorySetupFunc lsMemorySetupFunc;
#else
	//Other platforms run the parallel code sequentially (until pthread support or other parallel implementation is added)
	#include "BulletMultiThreaded/SequentialThreadSupport.h"
	typedef SequentialThreadFunc ThreadFunc;
	typedef SequentiallsMemorySetupFunc lsMemorySetupFunc;
#endif

//void example_ThreadFunc(void* userPtr, void* lsMemory) {}
//void* example_lsMemorySetupFunc() { return 0; }
inline btThreadSupportInterface* createThreadInterface(const char* uniqueThreadName, ThreadFunc userThreadFunc, lsMemorySetupFunc lsMemoryFunc,
													   int numThreads = 1, int threadStackSize = 65535)
{
	btThreadSupportInterface* interface = 0;
	
	#ifdef _WIN32
		Win32ThreadSupport::Win32ThreadConstructionInfo TCI(uniqueThreadName, userThreadFunc, lsMemoryFunc, numThreads, threadStackSize);
		void* ptr = btAlignedAlloc( sizeof(Win32ThreadSupport), 16 );
		interface = new(ptr) Win32ThreadSupport(TCI);
	#elif defined (USE_PTHREADS)
		PosixThreadSupport::ThreadConstructionInfo TCI(uniqueThreadName, userThreadFunc, lsMemoryFunc, numThreads, threadStackSize);
		void* ptr = btAlignedAlloc( sizeof(PosixThreadSupport), 16 );
		interface = new(ptr) PosixThreadSupport(TCI);
	#else
		SequentialThreadSupport::SequentialThreadConstructionInfo TCI(uniqueThreadName, userThreadFunc, lsMemoryFunc);
		void* ptr = btAlignedAlloc( sizeof(SequentialThreadSupport), 16 );
		interface = new(ptr) SequentialThreadSupport(TCI);
	#endif
	
	return interface;
}
inline void destroyThreadInterface(btThreadSupportInterface* interface)
{
	interface->~btThreadSupportInterface();
	btAlignedFree(interface);
}



typedef void (*btParallelForFunction) (void* parameters, int index);

///Simple implementation of 'parallel for' using the BulletMultiThreaded library.
class btParallelFor
{
	btThreadSupportInterface* m_threadInterface;
	btCriticalSection* m_mutex;
	
public:
	btParallelFor(const char* uniqueName, unsigned int numThreads = 1)
	{
		m_threadInterface = createThreadInterface(uniqueName, btParallelFor::mainFunction, btParallelFor::memorySetupFunc, numThreads);
		m_mutex = m_threadInterface->createCriticalSection();
	}
	~btParallelFor()
	{
		m_threadInterface->deleteCriticalSection(m_mutex);
		destroyThreadInterface(m_threadInterface);
	}
	
	///Executes the following loop in parallel:
	///@code
	///for(int i = firstIndex; i <= lastIndex; ++i) function(parameters, i);
	///@endcode
	///@param grainSize Number of indicies for a single thread to execute before syncronization.
	///Overly high values lead to an uneven distribution of work among threads, 
	///while extremely low values increase the syncronization overhead between threads.
	void execute(btParallelForFunction function, void* parameters, int firstIndex, int lastIndex, int grainSize)
	{
		if( !(firstIndex <= lastIndex) || grainSize <= 0 ) 
		{
			btAssert(0);
			return;
		}
	
		btParallelFor::IndexPool P(m_mutex, function, parameters, firstIndex, lastIndex, grainSize);
		
		//Defined in SpuCollisionTaskProcess.h; if this value is not used, either an assert will occur or nothing will happen.
		const int UICOMMAND = CMD_GATHER_AND_PROCESS_PAIRLIST;	
		
		//Execute btParallelFor::mainFunction() with userPtr == &P on all threads
		for(int i = 0; i < m_threadInterface->getNumTasks(); ++i)
			m_threadInterface->sendRequest(UICOMMAND, reinterpret_cast<ppu_address_t>(&P), i);
		
		//
		for(int i = 0; i < m_threadInterface->getNumTasks(); ++i) 
		{
			unsigned int taskId, status;
			m_threadInterface->waitForResponse(&taskId, &status);
		}
	}
	
private:
	static void* memorySetupFunc() { return 0; }
	static void mainFunction(void* userPtr, void* lsMemory) 
	{
		btParallelFor::IndexPool* P = static_cast<btParallelFor::IndexPool*>(userPtr);
		
		int firstIndex, lastIndex;
		for(bool indiciesValid = P->getIndexRange(firstIndex, lastIndex); indiciesValid; 
				indiciesValid = P->getIndexRange(firstIndex, lastIndex) )
		{
			for(int i = firstIndex; i <= lastIndex; ++i) P->m_function(P->m_parameters, i);
		}
	}
	
	///Contains a pool of indicies that multiple threads may draw from simultaneously
	struct IndexPool
	{
		btCriticalSection* m_mutex;

		btParallelForFunction m_function;
		void* m_parameters;
		
		int m_firstIndex;
		int m_lastIndex;
		int m_grainSize;
		
		IndexPool(btCriticalSection* mutex, btParallelForFunction function, void* parameters, int firstIndex, int lastIndex, int grainSize)
		: m_mutex(mutex), m_function(function), m_parameters(parameters), m_firstIndex(firstIndex), m_lastIndex(lastIndex), m_grainSize(grainSize) {}

		///Returns true if the index range in [out_firstIndex, out_lastIndex] is valid
		bool getIndexRange(int& out_firstIndex, int& out_lastIndex)
		{
			m_mutex->lock();
				int currentFirstIndex = m_firstIndex;	//	should use 'atomic fetch and add' here 
				m_firstIndex += m_grainSize;
			m_mutex->unlock();
			
			int potentialLastIndex = currentFirstIndex + m_grainSize - 1;
			if(potentialLastIndex <= m_lastIndex)
			{
				out_firstIndex = currentFirstIndex;
				out_lastIndex = potentialLastIndex;
				return true;
			}
			else if(currentFirstIndex <= m_lastIndex)
			{
				out_firstIndex = currentFirstIndex;
				out_lastIndex = m_lastIndex;
				return true;
			}
			
			//Set indicies such that this will not execute: for(int i = firstIndex; i <= lastIndex; ++i)
			out_firstIndex = -1;
			out_lastIndex = -2;
			return false;
		}
	};
};


#endif