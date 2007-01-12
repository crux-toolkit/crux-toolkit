#include <stdio.h>
#include <stdlib.h>
#include <fcntl.h>
#include <unistd.h>
#include <sys/types.h>
#include <sys/mman.h>

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
  /*int bufsize = sizeof(char) * pagesize;
  c = (char *) malloc(bufsize);
  int retval = read(fd, c, pagesize * sizeof(char));
   printf("%s\n", c);*/

  // data = mmap((caddr_t)0, pagesize, PROT_READ, MAP_SHARED, fd, pagesize);

  int offset = 0;
  srand(10000);


  char *array = malloc(5 * sizeof(char));
  int max_offset = pagesize * page_number;
  printf("max = %i\n", max_offset);
  int idx;
  char temp; 
  for (idx=0; idx < iterations; idx++){
    double a = rand()/((double)RAND_MAX+1);
    int offset = (int)(max_offset * a);
    lseek(fd, offset, SEEK_SET);
    read(fd, &temp, 1);
    array[idx%5] = temp;
  }


}
