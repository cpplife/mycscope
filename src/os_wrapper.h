

int os_get_cpu_core_count( void );

typedef struct {
	int   file_handle;
	char* buffer;
	int   size;
} mmap_info_t;

int os_mmap( const char* file_path, mmap_info_t* out_info );
void os_munmap( mmap_info_t* info );
