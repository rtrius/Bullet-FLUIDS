/* BulletMultiThreadedSupport.h
*/
#ifndef BULLET_MULTI_THREADED_SUPPORT_H_INCLUDED
#define BULLET_MULTI_THREADED_SUPPORT_H_INCLUDED

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
inline btThreadSupportInterface* createThreadInterface(const char *uniqueThreadName, ThreadFunc userThreadFunc, lsMemorySetupFunc lsMemoryFunc,
													   int numThreads = 1, int threadStackSize = 65535)
{
	btThreadSupportInterface *interface = 0;

	#ifdef _WIN32
		Win32ThreadSupport::Win32ThreadConstructionInfo TCI(uniqueThreadName, userThreadFunc, lsMemoryFunc, numThreads, threadStackSize);
		interface = new Win32ThreadSupport(TCI);
	#elif defined (USE_PTHREADS)
		PosixThreadSupport::ThreadConstructionInfo TCI(uniqueThreadName, userThreadFunc, lsMemoryFunc, numThreads, threadStackSize);
		interface = new PosixThreadSupport(TCI);
	#else
		SequentialThreadSupport::SequentialThreadConstructionInfo TCI(uniqueThreadName, userThreadFunc, lsMemoryFunc);
		interface = new SequentialThreadSupport(TCI);
	#endif
	
	return interface;
}
//void destroyThreadInterface(btThreadSupportInterface *interface) { delete interface; }




struct IndexRange
{
	int m_firstIndex;
	int m_lastIndex;
	IndexRange(int firstIndex, int lastIndex) : m_firstIndex(firstIndex), m_lastIndex(lastIndex) {}
};

typedef void (*ParallelForFunction) (void *parameters, int index);
struct IndexPool
{
	btCriticalSection *m_mutex;

	ParallelForFunction m_function;
	void *m_parameters;
	
	int m_firstIndex;
	int m_lastIndex;
	int m_grainSize;
	
	IndexPool(btCriticalSection *mutex, ParallelForFunction function, void *parameters, int firstIndex, int lastIndex, int grainSize)
	: m_mutex(mutex), m_function(function), m_parameters(parameters), m_firstIndex(firstIndex), m_lastIndex(lastIndex), m_grainSize(grainSize) {}

	IndexRange getIndexRange()
	{
		m_mutex->lock();
			int currentFirstIndex = m_firstIndex;	//	should use 'atomic fetch and add' here 
			m_firstIndex += m_grainSize;
		m_mutex->unlock();
		
		IndexRange result(-1, -1);
		
		int potentialLastIndex = currentFirstIndex + m_grainSize - 1;
		if(potentialLastIndex <= m_lastIndex)
		{
			result.m_firstIndex = currentFirstIndex;
			result.m_lastIndex = potentialLastIndex;
		}
		else if(currentFirstIndex <= m_lastIndex)
		{
			result.m_firstIndex = currentFirstIndex;
			result.m_lastIndex = m_lastIndex;
		}
		
		return result;
	}
};

void* ParallelFor_MemorySetupFunc() { return 0; }
void ParallelFor_MainFunction(void* userPtr, void* lsMemory) 
{
	IndexPool *P = static_cast<IndexPool*>(userPtr);
	
	for( IndexRange IR = P->getIndexRange(); IR.m_firstIndex != -1; IR = P->getIndexRange() )
		for(int i = IR.m_firstIndex; i <= IR.m_lastIndex; ++i)
			P->m_function(P->m_parameters, i);
}

///@brief Simple implementation of 'parallel for' using the BulletMultiThreaded library.
class ParallelFor
{
	btThreadSupportInterface *m_threadInterface;
	btCriticalSection *m_mutex;
	
public:
	ParallelFor(const char *uniqueName, unsigned int numThreads = 1) 
	{
		m_threadInterface = createThreadInterface(uniqueName, ParallelFor_MainFunction, ParallelFor_MemorySetupFunc, numThreads);
		m_mutex = m_threadInterface->createCriticalSection();
	}
	~ParallelFor() 
	{
		//m_mutex is destroyed in m_threadInterface destructor
		delete m_threadInterface;
	}
	
	void execute(ParallelForFunction function, void *parameters, int firstIndex, int lastIndex, int grainSize)
	{
		IndexPool P(m_mutex, function, parameters, firstIndex, lastIndex, grainSize);
		
		//
		const int UICOMMAND = CMD_GATHER_AND_PROCESS_PAIRLIST;	//defined in SpuCollisionTaskProcess.h
		for(int i = 0; i < m_threadInterface->getNumTasks(); ++i)
			m_threadInterface->sendRequest(UICOMMAND, reinterpret_cast<ppu_address_t>(&P), i);
		
		//
		for(int i = 0; i < m_threadInterface->getNumTasks(); ++i) 
		{
			unsigned int taskId, status;
			m_threadInterface->waitForResponse(&taskId, &status);
		}
	}
};


#endif