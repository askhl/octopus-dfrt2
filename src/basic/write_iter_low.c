/*
 Copyright (C) 2002 M. Marques, A. Castro, A. Rubio, G. Bertsch

 This program is free software; you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation; either version 2, or (at your option)
 any later version.

 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.

 You should have received a copy of the GNU General Public License
 along with this program; if not, write to the Free Software
 Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA
 02111-1307, USA.

 $Id: write_iter_low.c 6722 2010-06-13 12:44:43Z acastro $
*/

#include <stdio.h>
#include <stdlib.h>
#include <errno.h>
#include <string.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <unistd.h>

#include <config.h>
#include "string_f.h"

#define max(a, b) ((a) > (b) ? (a) : (b))
#define min(a, b) ((a) < (b) ? (a) : (b))

#define CHUNK 1024

typedef struct {
	char *filename, *buf;
	size_t iter, pos, size;
	double dt;
} write_iter;

/* internal functions */
static void write_iter_realloc(write_iter *w, int s)
{
	if(w->pos+s+1 <= w->size) return;
	w->size += CHUNK;
	w->buf = (char *)realloc(w->buf, w->size);
}

static void write_iter_string_work(write_iter *w, char *s)
{
	int l;

	l = strlen(s);
	write_iter_realloc(w, l);
	strcpy(w->buf+w->pos, s);
	w->pos += l;
}

static char* _str_center(char *s1, int l2)
{
	char *s2;
	int i, j, l1;
	l1 = strlen(s1);

	s2 = malloc(l2+1);
	for(i=0; i<(l2-l1)/2; i++) s2[i] = ' ';
	for(j=0; j<l1 && i<l2; j++, i++) s2[i] = s1[j];
	for(; i<l2; i++) s2[i] = ' ';
	s2[l2] = '\0';

	return s2;
}

void write_iter_header_work(write_iter *w, char *s)
{
	char *c;
	c = _str_center(s, 20);
	write_iter_string_work(w, c);
	free(c);
}

/* Functions called from FORTRAN */
void FC_FUNC_(write_iter_init, WRITE_ITER_INIT)
#ifdef SINGLE_PRECISION
		 (void **v, int *i, float  *d, STR_F_TYPE fname STR_ARG1)
#else
		 (void **v, int *i, double *d, STR_F_TYPE fname STR_ARG1)
#endif
{
	write_iter *w;
	w = (write_iter *)malloc(sizeof(write_iter));
	TO_C_STR1(fname, w->filename);
	w->buf = NULL;
	w->pos = w->size = 0;
	w->iter = *i;
	w->dt = *d;

	*v = w;
}

void FC_FUNC_(write_iter_clear, WRITE_ITER_CLEAR)(void **v)
{
	write_iter *w=*v;

	if(creat(w->filename, 0666) == -1){
		fprintf(stderr, "Could not create file '%s' (%s)", w->filename, strerror(errno));
		exit(1);
	}
}

void FC_FUNC_(write_iter_flush, WRITE_ITER_FLUSH)(void **v)
{
	int fd;
	write_iter *w=*v;
	if(!w->buf) return;

	fd = open(w->filename, O_WRONLY | O_CREAT | O_APPEND, 0666);
	if(fd == -1){
		fprintf(stderr, "Could not open file '%s' (%s)", w->filename, strerror(errno));
		exit(1);
	}
	
	write(fd, w->buf, w->pos);
	close(fd);

	w->pos = 0;
}

void FC_FUNC_(write_iter_end, WRITE_ITER_END)
		 (void **v)
{
	write_iter *w=*v;

	FC_FUNC_(write_iter_flush, WRITE_ITER_FLUSH)(v);

	free(w->filename);
	if(w->buf) free(w->buf);
	free(w);
}

void FC_FUNC_(write_iter_start, WRITE_ITER_START)
		 (void **v)
{
	write_iter *w=(write_iter *)*v;
	write_iter_realloc(w, 8+20);
	sprintf(w->buf+w->pos, "%8u%20.12e", (unsigned)(w->iter), w->iter*w->dt); w->pos+=8+20;
	w->iter++;
}

void FC_FUNC_(write_iter_double_1, WRITE_ITER_DOUBLE_1)
		 (void **v, double *d, int *no)
{
	write_iter *w=(write_iter *)*v;
	int i;

	write_iter_realloc(w, (*no)*20);
	for(i=0; i<*no; i++){
		sprintf(w->buf+w->pos, "%20.12e", d[i]); 
		w->pos += 20;
	}
}

void FC_FUNC_(write_iter_float_1, WRITE_ITER_FLOAT_1)
		 (void **v, float *d, int *no)
{
	write_iter *w=(write_iter *)*v;
	int i;

	write_iter_realloc(w, (*no)*20);
	for(i=0; i<*no; i++){
		sprintf(w->buf+w->pos, "%20.12e", d[i]); 
		w->pos += 20;
	}
}

void FC_FUNC_(write_iter_double_n, WRITE_ITER_DOUBLE_N)
		 (void **v, double *d, int *no)
{
	write_iter *w=(write_iter *)*v;
	int i;

	write_iter_realloc(w, (*no)*20);
	for(i=0; i<*no; i++){
		sprintf(w->buf+w->pos, "%20.12e", d[i]); 
		w->pos += 20;
	}
}

void FC_FUNC_(write_iter_float_n, WRITE_ITER_FLOAT_N)
		 (void **v, float *d, int *no)
{
	write_iter *w=(write_iter *)*v;
	int i;

	write_iter_realloc(w, (*no)*20);
	for(i=0; i<*no; i++){
		sprintf(w->buf+w->pos, "%20.12e", d[i]); 
		w->pos += 20;
	}
}

void FC_FUNC_(write_iter_int_1, WRITE_ITER_INT_1)
		 (void **v, int *d, int *no)
{
	write_iter *w=(write_iter *)*v;
	int i;

	write_iter_realloc(w, (*no)*8);
	for(i=0; i<*no; i++){
		sprintf(w->buf+w->pos, "%8d", d[i]); 
		w->pos += 8;
	}
}

void FC_FUNC_(write_iter_int_n, WRITE_ITER_INT_N)
		 (void **v, int *d, int *no)
{
	write_iter *w=(write_iter *)*v;
	int i;

	write_iter_realloc(w, (*no)*8);
	for(i=0; i<*no; i++){
		sprintf(w->buf+w->pos, "%8d", d[i]); 
		w->pos += 8;
	}
}

void FC_FUNC_(write_iter_string, WRITE_ITER_STRING)
		 (void **v, STR_F_TYPE s STR_ARG1)
{
  char *c;

  TO_C_STR1(s, c);
  write_iter_string_work((write_iter *)*v, c);
  free(c);
}

void FC_FUNC_(write_iter_header_start, WRITE_ITER_HEADER_START)
		 (void **v)
{
	write_iter *w=(write_iter *)*v;

	write_iter_string_work(w, "# Iter    ");
	write_iter_header_work(w, "t");
}

void FC_FUNC_(write_iter_header, WRITE_ITER_HEADER)
		 (void **v, STR_F_TYPE s STR_ARG1)
{
  char *c;
  TO_C_STR1(s, c);
  write_iter_header_work((write_iter *)*v, c);
  free(c);
}

void FC_FUNC_(write_iter_nl, WRITE_ITER_NL)
		 (void **v)
{
	write_iter_string_work((write_iter *)*v, "\n");
}
