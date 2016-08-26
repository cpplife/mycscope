#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <pthread.h>

#include "global.h"
#include "constants.h"
#include "bm_search.h"
#include "library.h"
#include "os_wrapper.h"

#define ALPHABET_LEN 256
#define NOT_FOUND patlen
#define max(a, b) ((a < b) ? b : a)

// delta1 table: delta1[c] contains the distance between the last
// character of pat and the rightmost occurrence of c in pat.
// If c does not occur in pat, then delta1[c] = patlen.
// If c is at string[i] and c != pat[patlen-1], we can
// safely shift i over by delta1[c], which is the minimum distance
// needed to shift pat forward to get string[i] lined up 
// with some character in pat.
// this algorithm runs in alphabet_len+patlen time.
void make_delta1(int *delta1, uint8_t *pat, int32_t patlen) {
    int i;
    for (i=0; i < ALPHABET_LEN; i++) {
        delta1[i] = NOT_FOUND;
    }
    for (i=0; i < patlen-1; i++) {
        delta1[pat[i]] = patlen-1 - i;
    }
}

// true if the suffix of word starting from word[pos] is a prefix 
// of word
int is_prefix(uint8_t *word, int wordlen, int pos) {
    int i;
    int suffixlen = wordlen - pos;
    // could also use the strncmp() library function here
    for (i = 0; i < suffixlen; i++) {
        if (word[i] != word[pos+i]) {
            return 0;
        }
    }
    return 1;
}

// length of the longest suffix of word ending on word[pos].
// suffix_length("dddbcabc", 8, 4) = 2
int suffix_length(uint8_t *word, int wordlen, int pos) {
    int i;
    // increment suffix length i to the first mismatch or beginning
    // of the word
    for (i = 0; (word[pos-i] == word[wordlen-1-i]) && (i < pos); i++);
    return i;
}

// delta2 table: given a mismatch at pat[pos], we want to align 
// with the next possible full match could be based on what we
// know about pat[pos+1] to pat[patlen-1].
//
// In case 1:
// pat[pos+1] to pat[patlen-1] does not occur elsewhere in pat,
// the next plausible match starts at or after the mismatch.
// If, within the substring pat[pos+1 .. patlen-1], lies a prefix
// of pat, the next plausible match is here (if there are multiple
// prefixes in the substring, pick the longest). Otherwise, the
// next plausible match starts past the character aligned with 
// pat[patlen-1].
// 
// In case 2:
// pat[pos+1] to pat[patlen-1] does occur elsewhere in pat. The
// mismatch tells us that we are not looking at the end of a match.
// We may, however, be looking at the middle of a match.
// 
// The first loop, which takes care of case 1, is analogous to
// the KMP table, adapted for a 'backwards' scan order with the
// additional restriction that the substrings it considers as 
// potential prefixes are all suffixes. In the worst case scenario
// pat consists of the same letter repeated, so every suffix is
// a prefix. This loop alone is not sufficient, however:
// Suppose that pat is "ABYXCDBYX", and text is ".....ABYXCDEYX".
// We will match X, Y, and find B != E. There is no prefix of pat
// in the suffix "YX", so the first loop tells us to skip forward
// by 9 characters.
// Although superficially similar to the KMP table, the KMP table
// relies on information about the beginning of the partial match
// that the BM algorithm does not have.
//
// The second loop addresses case 2. Since suffix_length may not be
// unique, we want to take the minimum value, which will tell us
// how far away the closest potential match is.
void make_delta2(int *delta2, uint8_t *pat, int32_t patlen) {
    int p;
    int last_prefix_index = patlen-1;

    // first loop
    for (p=patlen-1; p>=0; p--) {
        if (is_prefix(pat, patlen, p+1)) {
            last_prefix_index = p+1;
        }
        delta2[p] = last_prefix_index + (patlen-1 - p);
    }

    // second loop
    for (p=0; p < patlen-1; p++) {
        int slen = suffix_length(pat, patlen, p);
        if (pat[p - slen] != pat[patlen-1 - slen]) {
            delta2[patlen-1 - slen] = patlen-1 - p + slen;
        }
    }
}

static struct {
	/* parameters for output.*/
	FILE* output;
	char* format;
	/* search file count */
	int file_count;
	/* search pattern */
	char* pat;
	int* delta1;
	int* delta2;
	worker_info_t* workers;
	/* mutex for print*/
	pthread_mutex_t output_lock;

	/* mutex for global worker. */
	pthread_mutex_t global_worker_lock;
	/* cond for global worker. */
	pthread_cond_t global_worker_cond;
	/* global worker for all threads. */
	worker_info_t global_worker;
	/* flag for checking if all files have been added into global worker */
	int finish_of_adding_files;

} bm_search_data;


typedef struct {
	int* delta1;
	int* delta2;
	uint8_t* text;
	uint32_t text_len;
	uint8_t* pat;
	uint32_t pat_len;
	int line_number;
} bm_search_state_t;

static void bm_search_init( bm_search_state_t* state, 
		uint8_t *string, uint32_t stringlen, uint8_t *pat, uint32_t patlen )
{
	state->delta1 = (int*)malloc( ALPHABET_LEN * sizeof(int) );
	state->delta2 = (int*)malloc( patlen * sizeof(int) );

    make_delta1(state->delta1, pat, patlen);
    make_delta2(state->delta2, pat, patlen);

	state->text = string;
	state->text_len = stringlen;
	state->pat = pat;
	state->pat_len = patlen;

	state->line_number = 1;
}

static void bm_search_deinit( bm_search_state_t* state )
{
	free( state->delta1 );
	free( state->delta2 );
}

static void 
bm_search_set( bm_search_state_t* state, int* delta1, int* delta2,
		uint8_t *string, uint32_t stringlen, uint8_t *pat, uint32_t patlen )
{
	state->delta1 = delta1;
	state->delta2 = delta2;

	state->text = string;
	state->text_len = stringlen;
	state->pat = pat;
	state->pat_len = patlen;

	state->line_number = 1;
}

static uint8_t* boyer_moore( bm_search_state_t* state ) {
    int i;
    int* delta1 = state->delta1;
    int* delta2 = state->delta2;
	uint8_t* string = state->text;
	uint32_t stringlen = state->text_len;
	uint8_t* pat = state->pat;
	uint32_t patlen = state->pat_len;
	int step;

    // The empty pattern must be considered specially
    if (patlen == 0) {
        return string;
    }

    i = patlen-1;
	step = i;
    while (i < stringlen) {
		int k = 0;
		for ( ; k < step; ++k ) {
			if (string[i - k ] == '\n') ++state->line_number;
		}
        int j = patlen-1;
        while (j >= 0 && (string[i] == pat[j])) {
            --i;
            --j;
        }
        if (j < 0) {
            return (string + i+1);
        }

        step = max(delta1[string[i]], delta2[j]);
		i += step;
    }
    return NULL;
}

static int get_line_number( uint8_t* buffer, int buffer_offset )
{
	int line = 1;
	int i = 0;
	while ( i < buffer_offset ) {
		if ( buffer[i] == '\n' ) ++line;
		++i;
	}
	return line;
}

static int is_printable_char( uint8_t c )
{
	return c >= 20 && c <= 126;
}

static int is_line_ending_char( char c )
{
	return c == '\n' || c == '\r';
}

int bm_search(char *file, FILE *output, char *format, char* pat)
{
	FILE* fptr;
	int line_number = 1;

	fptr = fopen(file, "rb");
    if (fptr == NULL) return(-1);
	fseek( fptr, 0, SEEK_END );
	size_t size = ftell( fptr );
	uint8_t* buffer = (uint8_t*)mymalloc( size );
	if ( buffer != NULL ) {
		fseek( fptr, 0, SEEK_SET );
		int rv = fread( buffer, 1, size, fptr );
		if ( rv == size ) {
			int pat_len = (int)strlen( pat );
			int buffer_offset = 0;


			bm_search_state_t state;
			bm_search_init( &state, buffer, size, (uint8_t*)pat, pat_len );
			uint8_t* target = boyer_moore( &state );
			while ( target != NULL ) {
				line_number = state.line_number;
				fprintf( output, format, file, line_number );
				const char* p = (const char*)target;
				while ( *p != '\n' && (uint8_t*)p >= buffer ) {
					--p;
				}
				++p;
				while ( *p != '\n' && ((uint8_t*)p - buffer) < size  ) {
					if ( is_printable_char( *p ) ) {
						putc( *p, output );
					}
					++p;
				}
				putc( '\n', output );
				buffer_offset = (target - buffer ) + pat_len;
				target = NULL;
				if ( buffer_offset < size ) {
					state.text = buffer + buffer_offset;
					state.text_len = size - buffer_offset;
					target = boyer_moore( &state );
				}
			}
			bm_search_deinit( &state );
		}

		free( buffer );
	}

	fclose(fptr);
	return 0;
}


static void bm_search_set_match( match_info_t* match,
		char* filename, int line_number, size_t offset )
{
	match->filename    = filename;
	match->line_number = line_number;
	match->offset      = offset;
	match->next        = NULL;
}

static match_info_t* 
bm_search_match(char *file, char* pat, int* delta1, int* delta2, match_info_t** last )
{
	FILE* fptr;
	size_t size;
	uint8_t* buffer, *target;
	int pat_len, buffer_offset;
	match_info_t* match_list, *match_last;

	match_list = NULL; match_last = NULL;

	fptr = fopen(file, "rb");
    if (fptr == NULL) return NULL;
	fseek( fptr, 0, SEEK_END );
	size = ftell( fptr );
	buffer = (uint8_t*)mymalloc( size );
	if ( buffer != NULL ) {
		fseek( fptr, 0, SEEK_SET );
		if ( fread( buffer, 1, size, fptr ) == size ) {
			pat_len = (int)strlen( pat );
			buffer_offset = 0;

			bm_search_state_t state;
			bm_search_set( &state, delta1, delta2, buffer, size, (uint8_t*)pat, pat_len ); 
			target = boyer_moore( &state );
			while ( target != NULL ) {

				match_info_t* m = malloc( sizeof( match_info_t ) );
				bm_search_set_match( m, file, state.line_number, target - buffer );
				if ( match_list == NULL ) {
					match_list = m;
					match_last = m;
				}
				else {
					match_last->next = m;
					match_last = m;
				}

				buffer_offset = (target - buffer ) + pat_len;
				target = NULL;
				if ( buffer_offset < size ) {
					state.text = buffer + buffer_offset;
					state.text_len = size - buffer_offset;
					target = boyer_moore( &state );
				}
			}
		}

		free( buffer );
	}

	fclose(fptr);
	if ( last != NULL ) *last = match_last;
	return match_list;
}

static void* bm_search_worker( void* p )
{
	int i = *(int*)p;
	file_info_t* f;
	match_info_t* m, *l;
	worker_info_t* w = &bm_search_data.workers[i];

	f = w->file_list;
	while ( f != NULL ) {
		m = bm_search_match( f->filename, bm_search_data.pat, 
				bm_search_data.delta1, bm_search_data.delta2, &l );
		if ( w->match_list == NULL ) {
			w->match_list = m;
			w->match_last = l;
		}
		else if ( m != NULL ) {
			w->match_last->next = m;
			w->match_last = l;
		}
		f = f->next;
	}
}

static void
bm_search_output_match( match_info_t* m, FILE* output, char* fmt, char* buffer, size_t buffer_len )
{
	int offset;
	static const int PREFIX_MAX_LEN = 64;
	static const int POSTFIX_MAX_LEN = 64;
	int rv;
	const char* p;
	int char_count;

	pthread_mutex_lock( &bm_search_data.output_lock );

	fprintf( output, fmt, m->filename, m->line_number );

	offset = m->offset;
	{
		p = (char*)buffer + offset;
		char_count = PREFIX_MAX_LEN;
		while ( p >= buffer && !is_line_ending_char( *p ) && char_count > 0) {
			--p;
			--char_count;
		}
		++p;
		char_count = PREFIX_MAX_LEN + POSTFIX_MAX_LEN;
		while ( !is_line_ending_char( *p ) && (p - buffer ) < buffer_len && char_count > 0 ) {
			if ( is_printable_char( *p ) ) {
				putc( *p, output );
				--char_count;
			}
			++p;
		}
		putc( '\n', output );
	}
	pthread_mutex_unlock( &bm_search_data.output_lock );
}

static void 
bm_search_and_output_match(char *file, char* pat, int* delta1, int* delta2 )
{
	FILE* fptr;
	size_t size;
	uint8_t* buffer, *target;
	int pat_len, buffer_offset;
	match_info_t* match_list, *match_last, *m;

	match_list = NULL; match_last = NULL;

	fptr = fopen(file, "rb");
    if (fptr == NULL) return;
	fseek( fptr, 0, SEEK_END );
	size = ftell( fptr );
	buffer = (uint8_t*)mymalloc( size );
	if ( buffer != NULL ) {
		fseek( fptr, 0, SEEK_SET );
		if ( fread( buffer, 1, size, fptr ) == size ) {
			pat_len = (int)strlen( pat );
			buffer_offset = 0;

			bm_search_state_t state;
			bm_search_set( &state, delta1, delta2, buffer, size, (uint8_t*)pat, pat_len ); 
			target = boyer_moore( &state );
			while ( target != NULL ) {

				match_info_t* m = malloc( sizeof( match_info_t ) );
				bm_search_set_match( m, file, state.line_number, target - buffer );
				if ( match_list == NULL ) {
					match_list = m;
					match_last = m;
				}
				else {
					match_last->next = m;
					match_last = m;
				}

				buffer_offset = (target - buffer ) + pat_len;
				target = NULL;
				if ( buffer_offset < size ) {
					state.text = buffer + buffer_offset;
					state.text_len = size - buffer_offset;
					target = boyer_moore( &state );
				}
			}
		}
		/* output the match result */
		m = match_list;
		while ( m != NULL ) {
			bm_search_output_match( m, bm_search_data.output, bm_search_data.format, buffer, size );
			m = m->next;
		}

		free( buffer );
	}
	/* clear the match list */
	m = match_list;
	while ( m != NULL ) {
		match_info_t* temp = m;
		m = m->next;
		free( temp );
	}

	fclose(fptr);
	return;
}

static void 
bm_search_and_output_match_with_mmap(char *file, char* pat, int* delta1, int* delta2 )
{
	mmap_info_t mmap_info;
	size_t size;
	uint8_t* buffer, *target;
	int pat_len, buffer_offset;
	match_info_t* match_list, *match_last, *m;

	match_list = NULL; match_last = NULL;

	if ( os_mmap( file, &mmap_info ) != 0 ) return;

	size = (size_t)mmap_info.size;
	buffer = (uint8_t*)mmap_info.buffer;

	if ( buffer != NULL ) {
		pat_len = (int)strlen( pat );
		buffer_offset = 0;

		bm_search_state_t state;
		bm_search_set( &state, delta1, delta2, buffer, size, (uint8_t*)pat, pat_len ); 
		target = boyer_moore( &state );
		while ( target != NULL ) {

			match_info_t* m = malloc( sizeof( match_info_t ) );
			bm_search_set_match( m, file, state.line_number, target - buffer );
			if ( match_list == NULL ) {
				match_list = m;
				match_last = m;
			}
			else {
				match_last->next = m;
				match_last = m;
			}

			buffer_offset = (target - buffer ) + pat_len;
			target = NULL;
			if ( buffer_offset < size ) {
				state.text = buffer + buffer_offset;
				state.text_len = size - buffer_offset;
				target = boyer_moore( &state );
			}
		}

		/* output the match result */
		m = match_list;
		while ( m != NULL ) {
			bm_search_output_match( m, bm_search_data.output, bm_search_data.format, buffer, size );
			m = m->next;
		}
	}
	/* clear the match list */
	m = match_list;
	while ( m != NULL ) {
		match_info_t* temp = m;
		m = m->next;
		free( temp );
	}

	os_munmap( &mmap_info );
	return;
}

static void* bm_search_and_output_worker( void* p )
{
	int i = *(int*)p;
	file_info_t* f;
	worker_info_t* w = &bm_search_data.workers[i];
	char path[PATHLEN + 1];
	char* filename;
	filename = filepath( f->filename, path, PATHLEN + 1 );

	f = w->file_list;
	while ( f != NULL ) {
		bm_search_and_output_match( filename, bm_search_data.pat, 
				bm_search_data.delta1, bm_search_data.delta2 );
		f = f->next;
	}
}

extern int thread_worker_count;

int bm_search_worker_init( FILE* output, char* fmt, char* pat, int file_count )
{
	int i;
	int patlen;
	if ( thread_worker_count < 2 ) return -1;

	bm_search_data.output     = output;
	bm_search_data.format     = fmt;
	bm_search_data.file_count = file_count;
	bm_search_data.pat        = pat;

	patlen = strlen( pat );

	bm_search_data.delta1 = (int*)malloc( ALPHABET_LEN * sizeof(int) );
	bm_search_data.delta2 = (int*)malloc( patlen * sizeof(int) );

    make_delta1(bm_search_data.delta1, pat, patlen);
    make_delta2(bm_search_data.delta2, pat, patlen);

	bm_search_data.workers = malloc( sizeof(worker_info_t) * thread_worker_count );
	for ( i = 0; i < thread_worker_count; ++i ) {
		bm_search_data.workers[i].file_list  = NULL;
		bm_search_data.workers[i].file_last  = NULL;
		bm_search_data.workers[i].match_list = NULL;
		bm_search_data.workers[i].match_last = NULL;
	}

	 pthread_mutex_init( &bm_search_data.output_lock, NULL );
	 pthread_mutex_init( &bm_search_data.global_worker_lock, NULL );
	 pthread_cond_init( &bm_search_data.global_worker_cond, NULL );

	 bm_search_data.global_worker.file_list = NULL;
	 bm_search_data.global_worker.file_last = NULL;

	 bm_search_data.finish_of_adding_files = 0;

	return 0;
}

void bm_search_worker_deinit( void )
{
	file_info_t* f, *temp_f;
	match_info_t* m, *temp_m;
	int i;
	for ( i = 0; i < thread_worker_count; ++i ) {
		f = bm_search_data.workers[i].file_list;
		while ( f != NULL ) {
			temp_f = f->next;
			free( f );
			f = temp_f;
		}
		m = bm_search_data.workers[i].match_list;
		while ( m != NULL ) {
			temp_m = m->next;
			free( m );
			m = temp_m;
		}
	}

	f = bm_search_data.global_worker.file_list;
	while ( f != NULL ) {
		temp_f = f->next;
		free( f );
		f = temp_f;
	}

	pthread_cond_destroy( &bm_search_data.global_worker_cond );
	pthread_mutex_destroy( &bm_search_data.global_worker_lock );
	pthread_mutex_destroy( &bm_search_data.output_lock );

	free( bm_search_data.workers );
	free( bm_search_data.delta1 );
	free( bm_search_data.delta2 );

	bm_search_data.workers = NULL;
	bm_search_data.delta1 = NULL;
	bm_search_data.delta2 = NULL;

	bm_search_data.output = NULL;
	bm_search_data.pat    = NULL;

	bm_search_data.global_worker.file_list = NULL;
	bm_search_data.global_worker.file_last = NULL;

	bm_search_data.finish_of_adding_files = 0;
}

int bm_search_worker_add( int index, char* file )
{
	if ( index < 0 || index >= bm_search_data.file_count ) return -1;

	int file_count_for_on_worker = bm_search_data.file_count / thread_worker_count;

	int id = index / file_count_for_on_worker;
	if ( id > thread_worker_count - 1 ) id = thread_worker_count - 1;


	file_info_t* f = malloc( sizeof( file_info_t ) );
	f->filename = file;
	f->next = NULL;
	if ( bm_search_data.workers[id].file_list == NULL ) {
		bm_search_data.workers[id].file_list = f;
		bm_search_data.workers[id].file_last = f;
	}
	else {
		bm_search_data.workers[id].file_last->next = f;
		bm_search_data.workers[id].file_last = f;
	}

	return 0;
}

static void bm_search_print_match( match_info_t* m, FILE* output, char* fmt )
{
	FILE* fptr;
	int offset;
	static const int BUF_LEN = 128;
	char buf[BUF_LEN];
	int rv;
	const char* p;

	fprintf( output, fmt, m->filename, m->line_number );

	fptr = fopen(m->filename, "rb");
	if ( fptr == NULL ) return;

	offset = m->offset - BUF_LEN / 2;
	if (offset < 0 ) offset = 0;
	fseek( fptr, offset, SEEK_SET );
	rv = fread( buf, 1, BUF_LEN, fptr );
	if ( rv > 0 ) {
		p = buf + BUF_LEN / 2;
		while ( *p != '\n' && p >= buf ) {
			--p;
		}
		++p;
		while ( *p != '\n' && (p - buf ) < BUF_LEN  ) {
			if ( is_printable_char( *p ) ) {
				putc( *p, output );
			}
			++p;
		}
		putc( '\n', output );
	}
	fclose( fptr );
}

void bm_search_print_out( void )
{
	int i;
	match_info_t* m;

	for ( i = 0; i < thread_worker_count; ++i ) {
		m = bm_search_data.workers[i].match_list;
		while ( m != NULL ) {
			bm_search_print_match( m, bm_search_data.output, bm_search_data.format );
			m = m->next;
		}
	}
}


int bm_search_worker_run( void )
{
	int* thread_id;
	pthread_t* thread;
	int i, rv;

	thread = malloc( sizeof( pthread_t ) * thread_worker_count );
	thread_id = malloc( sizeof(int)*thread_worker_count );

	for ( i = 0; i < thread_worker_count; ++i ) {
		thread_id[i] = i;
		rv = pthread_create( &thread[i], NULL, bm_search_worker,  &thread_id[i] );
		if ( rv != 0 ) {
			return -1;
		}
	}

	/* block until all threads finish. */
	for ( i = 0; i < thread_worker_count; ++i ) {
		pthread_join( thread[i], NULL );
	}

	/* print out the result */
	bm_search_print_out();

	return 0;
}

int bm_search_and_output_worker_run( void )
{
	int* thread_id;
	pthread_t* thread;
	int i, rv;

	thread = malloc( sizeof( pthread_t ) * thread_worker_count );
	thread_id = malloc( sizeof(int)*thread_worker_count );

	for ( i = 0; i < thread_worker_count; ++i ) {
		thread_id[i] = i;
		rv = pthread_create( &thread[i], NULL, bm_search_and_output_worker,  &thread_id[i] );
		if ( rv != 0 ) {
			return -1;
		}
	}

	/* block until all threads finish. */
	for ( i = 0; i < thread_worker_count; ++i ) {
		pthread_join( thread[i], NULL );
	}

	return 0;
}

/*****************************************************************************/

static void* bm_search_and_output_from_global_worker( void* p )
{
	int i = *(int*)p;
	(void)i;
	file_info_t* f;
	char path[PATHLEN + 1];
	char* filename;

	while ( 1 ) {
		pthread_mutex_lock( &bm_search_data.global_worker_lock );

		while ( bm_search_data.global_worker.file_list == NULL ) {
			if ( bm_search_data.finish_of_adding_files ) {
				pthread_mutex_unlock( &bm_search_data.global_worker_lock );
				pthread_exit( NULL );
			}
			pthread_cond_wait( &bm_search_data.global_worker_cond, &bm_search_data.global_worker_lock );
		}
		
		f = bm_search_data.global_worker.file_list;
		bm_search_data.global_worker.file_list = bm_search_data.global_worker.file_list->next;
		if ( bm_search_data.global_worker.file_list == NULL ) {
			bm_search_data.global_worker.file_last = NULL;
		}
		pthread_mutex_unlock( &bm_search_data.global_worker_lock );

		filename = filepath( f->filename, path, PATHLEN + 1 );
		bm_search_and_output_match( 
				filename, bm_search_data.pat, 
				bm_search_data.delta1, bm_search_data.delta2 );

		free( f );
	}
}

static void bm_search_worker_add_into_global( char* file )
{
	file_info_t* f = malloc( sizeof( file_info_t ) );
	f->filename = file;
	f->next = NULL;

	pthread_mutex_lock( &bm_search_data.global_worker_lock );

	if ( bm_search_data.global_worker.file_list == NULL ) {
		bm_search_data.global_worker.file_list = f;
		bm_search_data.global_worker.file_last = f;
	}
	else {
		bm_search_data.global_worker.file_last->next = f;
		bm_search_data.global_worker.file_last = f;
	}

	pthread_cond_signal( &bm_search_data.global_worker_cond );
	pthread_mutex_unlock( &bm_search_data.global_worker_lock );
}


int bm_search_and_output_global_worker_run( char** file_name_list, int file_count )
{
	int* thread_id;
	pthread_t* thread;
	int i, rv;

	if ( file_count <= 0 ) return -1;

	
	bm_search_data.finish_of_adding_files = 0;

	/* create work thread */
	thread = malloc( sizeof( pthread_t ) * thread_worker_count );
	thread_id = malloc( sizeof(int)*thread_worker_count );

	for ( i = 0; i < thread_worker_count; ++i ) {
		thread_id[i] = i;
		rv = pthread_create( &thread[i], NULL, bm_search_and_output_from_global_worker,  &thread_id[i] );
		if ( rv != 0 ) {
			return -1;
		}
	}
	/* add files into global list */
	for ( i = 0; i < file_count; ++i ) {
		bm_search_worker_add_into_global( file_name_list[i] );
	}

	/* set the flag that finish global list.*/
	pthread_mutex_lock( &bm_search_data.global_worker_lock );
	bm_search_data.finish_of_adding_files = 1;
	pthread_cond_broadcast( &bm_search_data.global_worker_cond );
	pthread_mutex_unlock( &bm_search_data.global_worker_lock );

	/* block until all threads finish. */
	for ( i = 0; i < thread_worker_count; ++i ) {
		pthread_join( thread[i], NULL );
	}

	return 0;
}

