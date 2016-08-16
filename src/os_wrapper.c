#include "os_wrapper.h"

#ifdef _WIN32
#include <windows.h>
#endif // #ifdef _WIN32

int os_get_cpu_core_count( void )
{
#ifdef _WIN32
	SYSTEM_INFO si;
	GetSystemInfo(&si);
	int core_num = si.dwNumberOfProcessors;
	return core_num;
#else
	return 4;
#endif

}




