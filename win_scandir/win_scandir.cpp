// win_scandir.cpp : Defines the entry point for the console application.
//

#include <errno.h>
#include <stdio.h>
#include <vector>
#include "IndexFile.h"
#include "WinCrux.h"
#include "windirent.h"

using namespace std;

int is_text_file_name(const struct dirent *entry){

  const char* filename = entry->d_name; //w/o const gcc warning
  const char* suffix = "txt";

  int name_length = strlen(filename);
  int suffix_length = strlen(suffix);

  if( suffix_length > name_length){
    return 0;
  }

  // point to the last bit of filename
  filename += name_length - suffix_length;
  int matches = strncmp(filename, suffix, suffix_length); // 0 if matches

  if( matches == 0 ) {
    return 1;
  }//else

  return 0;
}

int win_scandir(
	const char *dirp,
	struct dirent ***namelist,
    int (*filter)(const struct dirent *),
    int (*compar)(const struct dirent **, const struct dirent **)
){
	errno = 0;
	
	DIR *directory = opendir(dirp);
	if (!directory) {
		return -1;
	}

	struct dirent *entry;
	vector<struct dirent *> *entries = new vector<struct dirent *>;
	while(entry = readdir(directory)) {
		if (!filter || filter(entry)) {
			struct dirent *valid_entry = new struct dirent;
			copy_dirent(entry, valid_entry);
			entries->push_back(valid_entry);
		}
	}
	if (errno) {
		return -1;
	}

	*namelist = static_cast<struct dirent **>(malloc(sizeof(struct dirent *) * entries->size()));
	copy(entries->begin(), entries->end(), *namelist);

	return entries->size();
}

int win_alphasort(const void *a, const void *b) {
	return 0;
}


int main(int argc, char* argv[])
{
	struct dirent **dirlist;

	int result = win_scandir(".", &dirlist, is_text_file_name, NULL);
	if (result < 0) {
		perror("win_scandir failed");
	}

	for(int i = 0; i < result; ++i) {
	  printf("Name: %s\n", dirlist[i]->d_name);
	}

	return 0;
}

