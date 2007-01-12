#include <stdio.h>
#include <stdlib.h>
#include <fcntl.h>
#include <unistd.h>
#include <sys/types.h>
#include <sys/mman.h>
#define RANDOM 0

int main(int argc, char** argv){
  if (argc < 4){
    fprintf(stderr, "Usage: test <fasta-file> <page-number> <iterations>\n");
    exit(1);
  }
  int fd, pagesize;

  if ( (fd = open(argv[1], O_RDONLY)) == -1){
    fprintf(stderr, "Failed to open %s!\n", argv[1]);
    exit(1);
  }

  int page_number = atoi(argv[2]);
  int iterations = atoi(argv[3]);
  pagesize = getpagesize();

  int offset = 0;
  char *data = mmap((caddr_t)0, page_number * pagesize, PROT_READ, MAP_SHARED, fd, offset);
  srand(10000);


  char *array = malloc(5 * sizeof(char));
  int max_offset = pagesize * page_number;
  printf("max = %i\n", max_offset);
  int idx;
  char temp;
  for (idx=0; idx < iterations; idx++){
    int offset;
    if (RANDOM==1){
      double a = rand()/((double)RAND_MAX+1);
      offset = (int)(max_offset * a);
    } else{
      offset = (int)((float)idx / iterations * max_offset);
    }
    temp = (char)data[offset];
    //printf("%c", temp);
    array[idx%5] = temp;
  }


}
