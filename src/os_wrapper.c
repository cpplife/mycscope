#include "os_wrapper.h"

#include <fcntl.h>
#ifdef _WIN32
#include <windows.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <io.h>
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




int os_mmap( const char* file_path, mmap_info_t* out_info )
{
	int fd;
    struct stat statbuf;
	int len;
	char* buf;

#ifdef _WIN32
    fd = _open(file_path, _O_RDONLY);
#else
    fd = open(file_path, O_RDONLY);
#endif
	if ( fd < 0 ) return -1;
	if ( fstat( fd, &statbuf ) != 0 ) {
		if ( fd != -1 ) close( fd );
		return -1;
	}

	len = (int)statbuf.st_size;
	if ( out_info->size == 0 ) {
		if ( fd != -1 ) close( fd );
		return 1;
	}
#ifdef _WIN32
    {
        HANDLE hmmap = CreateFileMapping( (HANDLE)_get_osfhandle(fd), 0, PAGE_READONLY, 0, len, NULL);
        buf = (char *)MapViewOfFile(hmmap, FILE_SHARE_READ, 0, 0, len);
        if (hmmap != NULL) CloseHandle(hmmap);
    }
    if (buf == NULL) {
		if ( fd != -1 ) close( fd );
		return -1;
    }
#else
    buf = mmap(0, len, PROT_READ, MAP_SHARED, fd, 0);
    if (buf == MAP_FAILED) {
		if ( fd != -1 ) close( fd );
		return -1;
    }
#endif

	out_info->file_handle = fd;
	out_info->buffer      = buf;
	out_info->size        = len;

	return 0;
}

void os_munmap( mmap_info_t* info )
{
	if ( info->buffer != NULL ) {
#ifdef _WIN32
		UnmapViewOfFile( info->buffer );
#else
		munmap( info->buffer, info->size );
#endif
	}
	if ( info->file_handle != -1 ) {
		close( info->file_handle );
	}
}
