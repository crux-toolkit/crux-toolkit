/*
 decoyFASTA.cpp - a tool for building decoy databases.  See usage() below for details.

***************************************************************************************
Copyright (C) Insilicos LLC 2008 All Rights Reserved.

* This library is free software; you can redistribute it and/or
* modify it under the terms of the GNU Lesser General Public
* License as published by the Free Software Foundation; either
* version 2.1 of the License, or (at your option) any later version.
*
* This library is distributed in the hope that it will be useful,
* but WITHOUT ANY WARRANTY; without even the implied warranty of
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
* Lesser General Public License for more details.
*
* You should have received a copy of the GNU Lesser General Public
* License along with this library; if not, write to the Free Software
* Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
*
****************************************************************************************/

#include <stdio.h>

void usage(const char *prog) {
	printf("%s: a tool for creating decoy search databases.  In its default operation, opens "
        "a FASTA formatted database file and writes it to a new file with the addition of a "
		" decoy reverse sequence entry for each protein.\n",prog);
	printf("usage: %s [options] <input_file> <output_file> [<filter_file>]\n",prog);
	puts("optional filter file is a FASTA file, only entries found in both filter and input are "
		"processed (actually it can just be file of FASTA headers, it looks for the leading > and ignores anything else)");
	puts("options:");
	puts("-no_orig");  
	puts("   do not copy the original proteins from the input file to the output file");  
	puts("-no_reverse");  
	puts("   do not reverse the sequences of the new decoy proteins being written to the output file");  
	puts("[-t[freq] <decoy_tag> [-t[freq] <2nd_decoy_tag> [-t[freq] <3rd...>]...]]");
	puts("   By default, decoy sequences have the tag \"decoy_\" added to the protein name, this can be customized with the -t switch.");
	puts("Use multiple -t switches to have decoy names randomly assigned from the provided tags.");
	puts("Use the [freq] option on -t to dictate relative frequency of decoy names.");
	puts("A negative [freq] value can be used to cause a portion of the input file to be skipped, which is useful for adjusting the ratio of decoys to actual proteins in your output database.");
	printf("example: %s -t5 decoy1_ -t3 decoy2_ -t-1 mydb.fasta decoydb.fasta\n will convert 5/9 (5+3+1=9) of the inputs to decoy1_, 3/9 to decoy2_, and skip 1/9.\n",prog);
}

#include <stdlib.h>
#include <string.h>
#include <string>
#include <vector>
#include <sstream>
#include <ctype.h>
#include <errno.h>

#ifndef min
#define min(a,b) ((a<b)?a:b)
#endif

class decoyTag {
public:
	decoyTag(const char *nameTag,double w=1.0) : name(nameTag),weight(w) {
		count = 0.0;
	}
	const char *name;
	double count;
	double weight; // use can specify relative frequency
};
static void write_protein(const char *description, const char *sequence, FILE *out);
static char description[128000]; // hack for unlikely crazy big description

int main(int argc, char *argv[]) {
	// allow multiple tags
	std::vector<decoyTag> decoyTags;
	srand(0); // we'll pick randomly from list of tags (although its usually length 1)
	bool bCopyOrig = true;
	bool bReverseDecoys = true;

	// check options for custom tag, etc
	for (int a=1;a<argc;a++) {
		if (*(argv[a])=='-') {
			int count = 1;
			switch (*(argv[a]+1)) {
			case 't': 
				{
					// any weighting?
					double weight = atof(argv[a]+2);
					if (weight < 0) {
						// negative weighting - we'll skip some
						decoyTags.push_back(decoyTag(NULL,-weight));
					} else {
						// custom decoy tag
						decoyTags.push_back(decoyTag(argv[a+1],weight));
						count++; // this is two args
					}
				}
				break;
			default:
				if (!strcmp(argv[a]+1,"no_orig")) {
					bCopyOrig = false;
				} else if (!strcmp(argv[a]+1,"no_reverse")) {
					bReverseDecoys = false;
				} else {
					printf("unknown arg %s - quit\n",argv[a]);
					exit(1);
				}
				break;
			}
			// eat this arg
			int t=a;
			for (int n=t+count;n<argc;) {
				argv[t++] = argv[n++];
			}
			argc -= count;
			a--; // don't advance - we moved the arglist down
		}
	} // end options processing
	if (!decoyTags.size()) {
		decoyTags.push_back(decoyTag("decoy_")); // the default
	}

	if (argc < 3) {
		const char *prog = strrchr(argv[0],'/');
		if (!prog) {
			prog = strrchr(argv[0],'\\');
		}
		if (!prog) {
			prog = argv[0];
		} else {
			prog++;
		}
		usage(prog);
		exit(1);
	}
	FILE *in = fopen(argv[1],"rb");
	if (!in) {
		printf("failed to load input file %s (\"%s\"), quitting\n",argv[1],strerror(errno));
		exit(1);
	}
	FILE *out = fopen(argv[2],"wb");
	if (!out) {
		printf("failed to open output file %s (\"%s\"), quitting\n",argv[2],strerror(errno));
		fclose(in);
		exit(1);
	}
	FILE *tmp = bCopyOrig?tmpfile():out; // stash the decoy sequences here until the end
	std::vector<std::string> filters;
	if (argc==4) {
		FILE *ff=fopen(argv[3],"rb");
		if (!ff) {
			printf("failed to load filter file %s(\"%s\"), quitting\n",argv[3],strerror(errno));
			exit(1);
		}     
		while (fgets(description,sizeof(description),ff)) {
			if ('>'==description[0]) {
				char *bar = strchr(description,'|');
				if (bar) {
					*bar = 0; // ignore annotations
				}
				filters.push_back(std::string(description+1));
			}
		}
		fclose(ff);
	}
	// deal with relative frequency of tags
	double totalWeight = 0;
	double totalDecoys = 0;
	bool gotWeights = false;
	for (int tg=decoyTags.size();tg--;) {
		totalWeight += decoyTags[tg].weight;
		gotWeights |= (decoyTags[tg].weight!=0.0);
	}
	if (gotWeights) {
		if (totalWeight <= 0) {
			printf("relative decoy frequencies total to 0 or less, quitting\n");
			exit(1);
		}
	} else {
		totalWeight = decoyTags.size();
	}

	for (int tg=decoyTags.size();tg--;) {
		if (gotWeights) {
			decoyTags[tg].weight/=totalWeight;
		} else {
			decoyTags[tg].weight = 1.0/totalWeight;
		}
	}
	int max_sequence=32000;
	char *sequence=(char *)malloc(max_sequence);
	bool first = true;
	while (fgets(description+!first,sizeof(description)-!first,in)) {
		first = false;
		char peek;
		sequence[0]=0;
		while ((peek=fgetc(in))!=EOF) {
			if ('>'==peek) {
				break;
			}
			if ('\n'==peek) {  // skip empty lines
				continue;
			}
			int sequence_len = strlen(sequence);
			if (sequence_len >= (max_sequence-1)) {
				max_sequence = max_sequence + max_sequence/2;
				sequence=(char *)realloc(sequence,max_sequence);
			}
			sequence[sequence_len++] = peek;
			char *fgot=fgets(sequence+sequence_len,max_sequence-sequence_len,in);
		}
		int len = strlen(description)-1; // ignore that trailing \n
		if (description[len-1]=='\r') {
			len--;
		}
		description[len]=0;

		// filter match?
		bool match;
		if (filters.size()) {
			match = false;
			for (int i=filters.size();i--;) {
				if (strstr(description,filters[i].c_str())) {
					match = true;
					break;
				}
			}
		} else {
			match = true;
		}
		for (int t=decoyTags.size();t--;) {
			if (match && decoyTags[t].name && strstr(description,decoyTags[t].name)) {
				match = false; // don't match on already-processed entries
				break;
			}
		}

		if (match) {

			if (bCopyOrig) {// write the original
				write_protein(description,sequence,out);
				if (ferror(out)) {
					perror("error writing output file");
					exit(1);
				}
			}
			// write the decoy
			totalDecoys++;
			std::string decoy_description(">");
			int tag;
			while (true) { // handle multi tag relative frequency
				tag = rand()%decoyTags.size();
				if (decoyTags[tag].count/totalDecoys <= decoyTags[tag].weight) {
					decoyTags[tag].count++; // we'll use another of these
					break;
				}
			}
			if (decoyTags[tag].name) { // NULL name means skip this one
				decoy_description += decoyTags[tag].name;
				std::stringstream count;
				count << (int)(decoyTags[tag].count);
				decoy_description += count.str();

				if (bReverseDecoys) {
					// write the decoy description
					if (!fwrite(decoy_description.c_str(),decoy_description.length(),1,tmp)) {
						perror("write error");
						exit(1);
					}
					fputc('\n',tmp);
					// write the decoy reversed sequence
					len=0;
					for (int sequence_len = strlen(sequence);sequence_len--;) {
						if (isalpha(sequence[sequence_len])) {
							fputc(sequence[sequence_len],tmp);
							if (++len>=80) {
								fputc('\n',tmp);
								len=0;
							}
						}
					}
					if(len) {
						fputc('\n',tmp);
					}
				} else { // write the decoy description and unreversed sequence
					write_protein(decoy_description.c_str(),sequence,out);
				}

				if (ferror(tmp)) {
					perror("error writing output file");
					exit(1);
				}
			} // end if we're writing this decoy
		} // end if match
	}
	fclose(in);

	if (bCopyOrig) {  // we wrote decoys to tmp
		// now tack the decoys onto the end of the output file
		rewind(tmp);
		size_t nread;
		while ((nread=fread(sequence,1,max_sequence,tmp))) {
			if ((nread!=fwrite(sequence,1,nread,out)) || ferror(out)) {
				perror("error writing output file");
				exit(1);
			}
		}
		fclose(tmp); 
	}
	fclose(out);
	free(sequence);
	return 0;
}

static void write_protein(const char *description, const char *sequence, FILE *out) {
	if (!fwrite(description,strlen(description),1,out)) {
		perror("write error");
		exit(1);
	}
	fputc('\n',out);
	int linelen=0;
	int seqlen = strlen(sequence);
	for (int orig_len = 0;orig_len < seqlen; orig_len++) {
	        if (isalpha(sequence[orig_len])) {
			fputc(sequence[orig_len],out);
			if (++linelen>=80) {
				fputc('\n',out);
				linelen=0;
			}
		} else if (linelen && ('\n'==sequence[orig_len])) {
			// may as well retain formatting for originals
			fputc('\n',out);
			linelen=0;
		}
	}
	if(linelen) {
		fputc('\n',out);
	}
}
