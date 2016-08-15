#include <stdio.h>

typedef struct {
	int id;
} thread_info_t;

typedef struct match_info_s {
	char*  filename;
	int    line_number;
	size_t offset; 
	struct match_info_s* next;
} match_info_t;

typedef struct file_info_s {
	char* filename;
	struct file_info_s* next;
} file_info_t;

typedef struct worker_info_s {
	file_info_t*  file_list;
	file_info_t*  file_last;
	match_info_t* match_list;
	match_info_t* match_last;
} worker_info_t;

extern worker_info_t* bm_search_info;

int bm_search(char *file, FILE *output, char *format, char* pattern);

int  bm_search_worker_init( FILE* output, char* fmt, char* pat, int file_count );
void bm_search_worker_deinit( void );
int  bm_search_worker_add( int index, char* file );
int  bm_search_worker_run( void );
int  bm_search_and_output_worker_run( void );

int  bm_search_and_output_global_worker_run( char** file_name_list, int file_count );
