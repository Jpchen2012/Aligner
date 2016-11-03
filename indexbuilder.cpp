/* 
Qalinger an ultrafast short read aligner
Copyright (C) 2016 BGI-shenzhen
Jianping Chen <chenjianping@genomics.cn>
Qaligner is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
Qaligner is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
You should have received a copy of the GNU General Public License along with this program. If not, see <http://www.gnu.org/licenses/>.
*/

/*This part of code is a reference index building program for the reference sequence smaller than 4G */

#include <stdio.h>
#include <stdlib.h>
#include <math.h> 
#include <string.h> 
#include <malloc.h>
#include <sys/time.h>
#include <stdint.h>

uint32_t c_to_i(char* ch) {
	switch (*ch) {
	case 'A':
		return 0;
	case 'a':
		return 0;
	case 'C':
		return 1;
	case 'c':
		return 1;
	case 'G':
		return 2;
	case 'g':
		return 2;
	case 'T':
		return 3;
	case 't':
		return 3;
	default:
		return 4;
	}
}


int main(int argv, char* argc[])
{
		struct timeval begin;
		struct timeval end;
		gettimeofday(&begin, NULL);
		
		int i = 0;
		int nline = 0;
		uint32_t j = 0;
		uint32_t p = 0;
		uint32_t in = 0; //number of nsegments
		uint32_t nchr = 0;
		uint32_t totnpos = 0;
		int ncon = 0; // length of continuous n
		int maxnchr = 10000;
		int maxnns = 40000;
		uint32_t* npos = (uint32_t*)calloc(maxnns, sizeof(uint32_t));
		uint32_t* lchr = (uint32_t*)calloc(maxnchr, sizeof(uint32_t));
		uint32_t* totlchr = (uint32_t*)calloc(maxnchr, sizeof(uint32_t));
		const uint64_t maxnseqint = 150000000; //150M 4G/32
		uint32_t maxnseeds = 268435456; //4**14
		uint32_t maxnpos = 2000000000; //4G/2
		uint64_t* seqint = (uint64_t*) calloc(maxnseqint, sizeof(uint64_t));
		uint32_t* seedind = (uint32_t*) calloc(maxnseeds, sizeof(uint32_t));
		uint16_t* seednum = (uint16_t*) calloc(maxnseeds, sizeof(uint32_t));
		uint32_t* seedindshift = (uint32_t*) calloc(maxnseeds, sizeof(uint32_t));
		uint32_t* seedpos = (uint32_t*) calloc(maxnpos, sizeof(uint32_t));
		uint32_t seedseq = 0;
		uint32_t mask = 0Xfffffff;
		int iseed = 0;
		uint32_t iseq = 0;
		char refpath[5000];
		char chrpath[5000];
		char infopath[5000];
		char pospath[5000];
		char snpath[5000];
		char sipath[5000];
		
		if(argv <= 1){
			fprintf(stderr, "usage: \n");
			fprintf(stderr, "\t indexbuilder reference.fa \n");
			fprintf(stderr, "The above command will write index to the directory where reference.fa is located. \n");
			fprintf(stderr, "If you want to choose a different path to store index, enter the following command: \n");
			fprintf(stderr, "\t indexbuilder reference.fa index/ref \n");
			exit(0);
		}
		
		if(argv ==2 ){
			sprintf(chrpath, "%s.chrs", argc[1]);
			sprintf(infopath, "%s.info", argc[1]);
			sprintf(pospath, "%s.spos", argc[1]);
			sprintf(snpath, "%s.sn", argc[1]);
			sprintf(sipath, "%s.si", argc[1]);
			sprintf(refpath, "%s.cref", argc[1]);
		}
		else if(argv == 3){
			sprintf(chrpath, "%s.chrs", argc[2]);
			sprintf(infopath, "%s.info", argc[2]);
			sprintf(pospath, "%s.spos", argc[2]);
			sprintf(snpath, "%s.sn", argc[2]);
			sprintf(sipath, "%s.si", argc[2]);
			sprintf(refpath, "%s.cref", argc[2]);
		}
		
		FILE *f = fopen(argc[1], "r");
		if(f == NULL){
			fprintf(stderr, " read reference file failed \n");
			exit(0);
		}

		char* seq;
		char seqid[1024];
		size_t seqlen = 0;
		ssize_t linelength;
		bool includen = true;
		
		FILE *fchr = fopen(chrpath, "a");
		while((linelength = getline(&seq, &seqlen, f)) != -1){
				nline++;
				if(*seq == '>'){
						nchr++;
						sscanf(seq, "%s", seqid);
						fprintf(fchr, "%s\n", (seqid+1));
				}
				else if(nchr > 0){
						for(i = 0; i < (linelength-1); i++){
								lchr[nchr-1] ++;
								p = j/32;
								j++;
								seedseq = mask & ((seedseq << 2)|(c_to_i(seq + i) & 3));
								*(seqint+p) = ((*(seqint + p))<< 2)|(c_to_i(seq + i) & 3);
							
								if((*(seq+i) == 'N')||(*(seq+i) == 'n')){
										iseed = 0;
								}
								else{
										iseed ++;
										if(iseed == 14){
												if(j%2 == 0){
													if(seednum[seedseq] < 65535)
														seednum[seedseq] ++;
												}
												iseed --;
										}
								}

								
								if((*(seq+i) == 'N')||(*(seq+i) == 'n')){
										ncon ++;
										if(ncon == 10){
												npos[in] = j - 9;
												in++;
										}
								}
								else{
										if(ncon >= 10){
											npos[in] = j - 1;
											in ++;
										}
										ncon = 0;
								}
								

						}
				}
		}
		fclose(f);
		fclose(fchr);
		
		if(nchr == 0){
			fprintf(stderr, " read reference file failed \n");
			exit(0);
		}
		
		if(j%32 != 0)
			*(seqint+p) = (*(seqint + p)) << (2*(32-j%32));
		totlchr[0] = lchr[0];
		for(i = 1; i < nchr; i++)
			totlchr[i] = totlchr[i-1] + lchr[i];
		
		int nmulseed = 0;
		int ntotseeds = 0;
		int nrepseed = 0;
		maxnseeds = 268435456;
		uint64_t sntot = 0;
		for(i = 0; i < maxnseeds; i++){
				if(seednum[i] > 0)
					ntotseeds ++;
				
				if(seednum[i] > 32){
					nrepseed ++;
				}

				if(seednum[i] > 65530){
						seednum[i] = 0;
						nmulseed ++;
				}
				else{
					sntot += seednum[i];
				}
				if(i < (maxnseeds-1))
					seedind[i+1] = seedind[i]+seednum[i];
				totnpos += seednum[i];
		}
		
		
		f = fopen(argc[1],"r");
		
		nline = 0;
		j = 0;
		seedseq = 0;
		iseed = 0;
		while((linelength = getline(&seq, &seqlen, f)) != -1){
				nline++;
				if(*seq != '>'){
						for(i = 0; i < (linelength-1); i++){
								
								j++;
								seedseq = mask & ((seedseq << 2)|(c_to_i(seq + i) & 3));
					
								if((*(seq+i) == 'N')||(*(seq+i) == 'n')){
										iseed = 0;
								}
								else{
										iseed ++;
										if(iseed == 14){
												if(j%2 == 0){
														if(seednum[seedseq] > 0){
															seedpos[seedind[seedseq]+seedindshift[seedseq]] = j - 13;
															seedindshift[seedseq]++;
														}
												}
												iseed --;
										}
								}
						}
				}
		}
		fclose(f);

		FILE *fpos = fopen(pospath,"wb+");
		if(fpos == NULL){
				printf(" create file failed \n");
		}
		else{
				fwrite(seedpos, sizeof(uint32_t), totnpos, fpos);
		}
		fclose(fpos);

		
		FILE *fs = fopen(snpath,"wb+");
		if(fs == NULL){
			printf(" create file failed \n");
		}
		else{
				fwrite(seednum, sizeof(uint16_t), maxnseeds, fs);
		}
		fclose(fs);
		
		FILE *fi = fopen(sipath,"wb+");
		if(fi == NULL){
			printf(" create file failed \n");
		}
		else{
				fwrite(seedind, sizeof(uint32_t), maxnseeds, fi);
		}
		fclose(fi);

		p++;
		
		FILE *fn = fopen(infopath,"wb+");
		if(fn == NULL){
			printf(" create file failed \n");
		}
		else{
			fwrite(&nchr, sizeof(uint32_t), 1, fn);
			fwrite(&in, sizeof(uint32_t), 1, fn);
			fwrite(&p, sizeof(uint32_t), 1, fn);
			fwrite(&totnpos, sizeof(uint32_t), 1, fn);
			fwrite(lchr, sizeof(uint32_t), nchr, fn);
			fwrite(totlchr, sizeof(uint32_t), nchr, fn);
			fwrite(npos, sizeof(uint32_t), in, fn);
		}
		fclose(fn);

		FILE *fc = fopen(refpath,"wb+");
		if(fc == NULL){
			printf(" create file failed \n");
		}
		else{
			fwrite(seqint, sizeof(uint64_t), p, fc);
		}
		fclose(fc);

		free(seqint);
		free(npos);
		free(lchr);
		free(totlchr);
		free(seedpos);
		free(seednum);
		free(seedind);
		free(seedindshift);
		gettimeofday(&end, NULL);
		printf("search index cost = %fsec\n",end.tv_sec - begin.tv_sec + float(end.tv_usec - begin.tv_usec) / 1000000);
		return 0;
}

