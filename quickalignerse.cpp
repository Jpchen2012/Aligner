/* Qalinger an ultrafast short read aligner
Copyright (C) 2016 BGI-shenzhen
Jianping Chen <chenjianping@genomics.cn>
Qaligner is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
Qaligner is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
You should have received a copy of the GNU General Public License along with this program. If not, see <http://www.gnu.org/licenses/>.
*/

/* This is a mapping program for single end short reads data produced by Hiseq or CPAS platform. */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <malloc.h>
#include <sys/time.h>
#include <stdint.h>
#include <pthread.h>
#include <zlib.h>

using namespace std;

int subsize, buffersize;	
int nthreads;	
int isid;
uint64_t totnfq = 0;
uint64_t totnprint = 0;
uint64_t totmapped = 0;
uint64_t totunimapped = 0;
uint64_t maxreflength = 1024;												
char chrlist[10000][1000];
char idname[200];
char *pret = NULL;
bool isspace = false;
bool isslash = false;
bool isend = false;

typedef struct read{
	char *sid; 
	char *seq; 
	char *qual;
}reads;

typedef struct result{
	int flag;
	int jc;
	uint32_t pos;
	int mapscore;
	int nmis;
	char cigar[20];
	bool iscomp;
	bool isprint;
}res;

typedef struct farg{
	reads* seqbuffer;
	res* resbuffer;
	int ip;
	int bufsize;
	uint32_t nmapped;
	uint32_t nunimapped;
	bool iscomp;
	bool isprint;
	uint32_t nprint;
	uint64_t* seqint;
	uint16_t* seednum;
	uint32_t* seedind;
	uint32_t* seedpos;
	uint32_t* chrpos;
	uint32_t nchr;
}args;

typedef struct datarg{
	args* arglist;
}dargs;

int match_num[65536];
const uint32_t mask[9] = {0, 0xf, 0xff, 0xfff, 0xffff, 0xfffff, 0xffffff, 0xfffffff, 0xffffffff};
const uint64_t maskl[33] = {0, 0x3, 0xf, 0x3f, 0xff, 0x3ff, 0xfff, 0x3fff, 0xffff, 0x3ffff, 0xfffff, 0x3fffff, 0xffffff, 
	0x3ffffff, 0xfffffff, 0x3fffffff, 0xffffffff, 0x3ffffffff, 0xfffffffff, 0x3fffffffff, 
	0xffffffffff, 0x3ffffffffff, 0xfffffffffff, 0x3fffffffffff, 0xffffffffffff,
	0x3ffffffffffff, 0xfffffffffffff, 0x3fffffffffffff, 0xffffffffffffff,
	0x3ffffffffffffff, 0xfffffffffffffff, 0x3fffffffffffffff, 0xffffffffffffffff};

int32_t match_32(int32_t n) {
	if (n == 0)
		return n;
	n = (n | (n << 1)) & 0x0AAAAAAA;
	n = (n & 0x55555555) + ((n >> 1) & 0x55555555);
	n = (n & 0x33333333) + ((n >> 2) & 0x33333333);
	n = (n & 0x0f0f0f0f) + ((n >> 4) & 0x0f0f0f0f);
	n = (n & 0x00ff00ff) + ((n >> 8) & 0x00ff00ff);
	n = (n & 0x0000ffff) + ((n >> 16) & 0x0000ffff);
	return n;
}

void initial_match_num() {
	int i = 0;
	while (i < 65536) {
		int match = match_32(i) & 0xf;
		int misFlagBit = 0;
		if (match > 4) {
			bool is_bp = false;
			int bp = 0;
			int ep = 0;
			int k = 0;
			int j = 7;
			while (j >= 0) {
				if (((i >> (j * 2)) & 3) > 0) {
					k++;
					if (k < 3)
						misFlagBit = misFlagBit | ((7 - j) << (4 * (4 - k)));
					else if (match - k < 2)
						misFlagBit = misFlagBit
							| ((7 - j) << (4 * (match - k)));
				}
				j--;
			}
		} else if (match <= 4 && match > 0) {
			int k = 0;
			int j = 7;
			while (j >= 0) {
				if (((i >> (j * 2)) & 3) > 0) {
					k++;
					misFlagBit =
						k == match ?
						((misFlagBit | ((7 - j) << (4 * (4 - k))))
						 | ((7 - j) & 0xf)) :
						(misFlagBit | ((7 - j) << (4 * (4 - k))));
				}
				j--;
			}
		}
		match_num[i] = (match << 16) | (misFlagBit & 0xffff);
		//match_num[i] = (match << 16);
		i++;
	}
}

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

char rev(char* ch){
	switch(*ch){
		case 'A':
			return 'T';
		case 'a':
			return 'T';
		case 'C':
			return 'G';
		case 'c':
			return 'G';
		case 'G':
			return 'C';
		case 'g':
			return 'C';
		case 'T':
			return 'A';
		case 't':
			return 'A';
		default:
			return 'N';
	}
}

char upper(char* ch){
	switch(*ch){
		case 'A':
			return 'A';
		case 'a':
			return 'A';
		case 'C':
			return 'C';
		case 'c':
			return 'C';
		case 'G':
			return 'G';
		case 'g':
			return 'G';
		case 'T':
			return 'T';
		case 't':
			return 'T';
		default:
			return 'N';
	}
}

char i_to_c(int a) {
	if (a == 0)
		return 'A';
	else if (a == 1)
		return 'C';
	else if (a == 2)
		return 'G';
	else if (a == 3)
		return 'T';
	else
		return 'N';
}

/*shifting sequences to see whether there are gaps */
int Shiftseq(uint64_t lastseqint, uint32_t ankpos, int direct, int sreadlen, int seffl, uint64_t *seqint){
	uint64_t targetint = 0;
	int gshift = 8;
	int igshift = 0;
	int shiftmis = 0;
	int minsm = 10;
	int modpos = ankpos%32;
	if(modpos == 0)
		modpos = 32;

	if((modpos >= (gshift+1))&&(modpos <= (33-8-gshift)))  
		targetint = seqint[ankpos/32] >> 2*(33-8-gshift-modpos);
	else if(modpos < (gshift+1))
		targetint = (maskl[8+gshift-1+modpos]&(seqint[ankpos/32]>>2*(33-8-gshift-modpos))) | (~maskl[8-1+gshift+modpos]&(seqint[ankpos/32-1]<<2*(8-1+gshift+modpos)));
	else
		targetint = (~maskl[modpos-33+8+gshift]&(seqint[(ankpos-1)/32] << 2*(modpos-33+8+gshift))) | (maskl[32-modpos+(33-8-gshift)]&(seqint[(ankpos-1)/32+1] >> 2*(32-modpos+33-8-gshift)));
	minsm = 5;
	igshift = -20;
	for(int ig = 0; ig <= 2*gshift; ig++){
		shiftmis = match_num[(lastseqint^(targetint>>2*ig))&maskl[8]] >> 16;

		if(shiftmis <= 2){
			if(shiftmis < minsm){
				igshift = ig;
				minsm = shiftmis;
			}
			else if((shiftmis == minsm)&&(abs(ig-gshift) < abs(igshift-gshift))){
				igshift = ig;
			}
		}
	}
	if(minsm > 2){
		if((sreadlen-seffl) <= 6){
			for(int ig = gshift-6; ig <= gshift+6; ig++){
				if(direct == 1)
					shiftmis = match_num[(lastseqint^(targetint>>2*ig))&maskl[4]] >> 16;
				else
					shiftmis = match_num[((lastseqint^(targetint>>2*ig))>>8)&maskl[4]] >> 16;
				if(shiftmis <= 1){
					if(shiftmis < minsm){
						igshift = ig;
						minsm = shiftmis;
					}
					else if((shiftmis == minsm)&&(abs(ig-gshift)<abs(igshift-gshift)))
						igshift = ig;
				}
			}
		}
	}

	return (igshift-gshift);
}

/* search mapping results for all reference positions of a seed */ 
void Gapalign(int ipos, int readlen, char* fqid, char* fqseq, uint64_t* fqseqint, uint32_t* fqseedint, uint64_t* seqint, uint16_t* seednum, uint32_t* seedind, uint32_t* seedpos, uint32_t* chrpos, uint32_t nchr, int* nmin, uint32_t* bestpos, int* nbest, int* gaplen, int* forwlen, int* backlen, char* cigar){
	int i = 0;
	int j = 0;
	int k = 0;
	int m = 0;
	int p = 1;
	int nmis = 0;
	int found = 0;
	int forw[100];
	int back[100];
	uint32_t spos, epos, ankpos;
	uint64_t lastseqint;
	int sp = 0;
	int ep = 0;
	int nmpos = 100;
	int normis = 15;
	int bytemis = 100;
	int shift = 0;
	int begmis = 10;
	int effmis = 0;
	int termis = 0;
	int termisb = 0;
	int nmisth = 0;
	int effl = 0;
	int oldeffl = 0;
	uint32_t gapl, delta;
	int kgap = 0;
	int searchindel = 0;
	int shiftsearch = 0;
	int igshift = 0;
	int nminmis = 0;
	int gapth = 0;
	uint32_t nshiftmis = 0;
	uint32_t nshiftmis2 = 0;
	uint32_t nshiftshit = 0;
	uint32_t nshiftgood = 0;
	uint64_t seqdot[100];
	int jg, kg, nforw, nback, kmax, imisp, ntotmis, deltap, realgap;
	
	for(i = 0; i < seednum[fqseedint[ipos]]; i++){
		if((seedpos[seedind[fqseedint[ipos]]+i] > (ipos + 32+readlen/2))&&(seedpos[seedind[fqseedint[ipos]]+i] < (maxreflength + ipos - 3*readlen/2 - 36))){
		spos = seedpos[seedind[fqseedint[ipos]]+i]-ipos;
		epos = spos;

		if((spos%32) == 0)
			shift = 62;
		else
			shift = 2*((spos%32) - 1);
		nmis = 0;
		bytemis = 0;
		begmis = 0;
		termis = 0;
		termisb = 0;
		effl = 0;
		for(j = 0; (j <= (readlen-1)/32)&&(nmis <= *nmin); j++){
			if((spos%32) == 1)
				seqdot[j] = fqseqint[j]^seqint[spos/32+j];
			else
				seqdot[j] = fqseqint[j]^(((seqint[(spos-1)/32+j]<<shift)&(~maskl[shift/2]))|((seqint[(spos-1)/32+j+1]>>(64-shift))& maskl[shift/2]));
			if(j >= readlen/32)
				seqdot[j] &= ~maskl[32-readlen%32];
			bytemis = 0;
			for(p = 0; (p < 4) && (nmis <= *nmin); p++){
				effmis = nmis - bytemis;
				normis = nmis;
				if(j < readlen/32)
					effl += 8;
				else{
					if((p+1)*8 <= readlen%32)
						effl += 8;
					else if( 8*p < readlen%32)
						effl += (readlen%32 - 8*p);
				}
				bytemis = match_num[(seqdot[j]>>(16*(3-p)))&maskl[8]] >> 16;
				nmis += (bytemis > 3) ? 1000:bytemis;
				if(j ==0){
					if(p < 3)
						begmis += bytemis;
					if(p < 4)
						termis += bytemis;
				}
				if((readlen-effl) < 33)
					termisb += bytemis;
			}
		}
		oldeffl = effl;
		shiftsearch = 0;
		if(nmis <= *nmin){
			if(readlen > 63){
				if((nmis > 1)&&(*nmin > 0)){
					if((termisb > 1)&&(termisb*10/nmis > 6)){
						if(termisb == 2){
							if(nmis == 2)
								shiftsearch = 1;
						}
						else
							shiftsearch = 1;
					}
					if((termis > 1)&&(termis*10/nmis > 6)){
						if(termis == 2){
							if(nmis == 2)
								shiftsearch = -1;
						}
						else
							shiftsearch = -1;
					}
				}
			}
			if((*nbest == 0)||(nmis < *nmin)){
				*nmin = nmis;
				*bestpos = spos;
				*nbest = 1;
				*gaplen = 0;
				sprintf(cigar, "%dM", readlen);
			}
			else if(*bestpos != spos){
				*nbest = 2;
			}
		}
		else if((*nmin > gapth)&&(readlen > 47)){
			searchindel = 0;
			if(effl >= 28){
				if(begmis <= nmisth)
					searchindel = 1;
				else{
					if(termis <= (nmisth+1))
						searchindel = 1;
					else{
						if(((effl-16)>= 28)&&(effmis <= (effl-16)*(6+2*nmisth)/100))
							searchindel = 1;
						else if(((effl-8)>= 28)&&(normis <= (effl-8)*(6+2*nmisth)/100))
							searchindel = 1;
					}
				}
				if(searchindel == 1){
					found = 0;
					for(j = (readlen - 14); (j >= effl-12)&&(found < 1); j--){
						sp = 0; 
						ep = seednum[fqseedint[j]]-1;
						gapl = 2500;
						while(ep >= sp){
							k = (sp+ep)/2;
							delta = (seedpos[seedind[fqseedint[j]]+k] >= (spos+j))? (seedpos[seedind[fqseedint[j]]+k]-spos-j): (spos+j - seedpos[seedind[fqseedint[j]]+k]);
							if((delta < gapl)&&(seedpos[seedind[fqseedint[j]]+k] > j)){
								gapl = delta;
								kgap = k;
							}
							
							if(seedpos[seedind[fqseedint[j]]+k] > (spos+j)){
								ep = k-1;
							}
							else if (seedpos[seedind[fqseedint[j]]+k] < (spos+j)){
								sp = k+1;
							}
							else 
								ep = sp -1;
						} 
						if(gapl < readlen/4){
							found = 1;
							if(gapl != 0)
								epos = seedpos[seedind[fqseedint[j]]+kgap] - j;
						}
					}
					if((found == 0)&&((readlen-effl) < 29))
						shiftsearch = 1;
				}
			}
			if((searchindel == 0)&&(effl <= 16)){
				effl = 0;
				bytemis = 0;
				nmis = 0;
				begmis = 0;
				termis = 0;
				termisb = 0;
				for(j = (readlen-1)/32; (j >= 0)&&(nmis <= *nmin); j--){
					if((spos%32) == 1)
						seqdot[j] = fqseqint[j]^seqint[spos/32+j];
					else
						seqdot[j] = fqseqint[j]^(((seqint[(spos-1)/32+j]<<shift)&(~maskl[shift/2]))|((seqint[(spos-1)/32+j+1]>>(64-shift))& maskl[shift/2]));
					if(j >= readlen/32)
						seqdot[j] &= ~maskl[32-readlen%32];
					bytemis = 0;
					for(p = 3; (p >= 0) && (nmis <= *nmin); p--){
						effmis = nmis - bytemis;
						normis = nmis;
						if(j < readlen/32)
							effl += 8;
						else{
							if((p+1)*8 <= readlen%32)
								effl += 8;
							else if( 8*p < readlen%32)
								effl += (readlen%32 - 8*p);
						}
						bytemis = match_num[(seqdot[j]>>(16*(3-p)))&maskl[8]] >> 16;
						nmis += (bytemis > 3) ? 1000:bytemis;
						if(effl <= 28)
							begmis += bytemis;
						if(effl <= 39)
							termis += bytemis;
						if((readlen-effl) < 33)
							termisb += bytemis;
					}
				}
				if(effl >= 28){
					if((begmis <= nmisth)||(termis <= (nmisth+1)))
						searchindel = 1;
					else{
						if(((effl-16)>= 28)&&(effmis <= (effl-16)*(6+2*nmisth)/100))
							searchindel = 1;
						else if(((effl-8)>= 28)&&(normis <= (effl-8)*(6+2*nmisth)/100))
							searchindel = 1;
					}
					if(searchindel == 1){
						found = 0;
						for(j = 0; (j < (readlen-effl))&&(found < 1); j++){
							sp = 0;
							ep = seednum[fqseedint[j]]-1;
							gapl = 2500;
							while(ep >= sp){
								k = (sp+ep)/2;
								delta = (seedpos[seedind[fqseedint[j]]+k] >= (spos+j))? (seedpos[seedind[fqseedint[j]]+k]-spos-j): (spos+j - seedpos[seedind[fqseedint[j]]+k]);
								if((delta < gapl)&&(seedpos[seedind[fqseedint[j]]+k] > j)){
									gapl = delta;
									kgap = k;
								}
								if(seedpos[seedind[fqseedint[j]]+k] > (spos+j)){
									ep = k-1;
								}
								else if (seedpos[seedind[fqseedint[j]]+k] < (spos+j)){
									sp = k+1;
								}
								else
									ep = sp -1;
							}
							if(gapl < readlen/4){
								found = 1;
								if(gapl != 0)
									spos = seedpos[seedind[fqseedint[j]]+kgap] - j;
							}
						}

						if((found == 0)&&((readlen-effl) < 29))
							shiftsearch = -1;
					}
				}
			}
		}

		if((shiftsearch != 0)&&(readlen > 47)){
			lastseqint = 0;
			ankpos = spos+readlen-8;
			igshift = -30;
			if(shiftsearch == 1){
				if(readlen % 32 == 0)
					lastseqint = fqseqint[readlen/32-1];
				else
					lastseqint = (maskl[readlen%32]&(fqseqint[readlen/32]>>(2*(32-readlen%32))))
						| (~maskl[readlen%32]&(fqseqint[readlen/32-1] << (2*(readlen%32))));
				lastseqint &= maskl[8];
				igshift = Shiftseq(lastseqint, ankpos, shiftsearch, readlen, effl, seqint);
				if(igshift > -9){
					epos = spos - igshift;
				}
			}
			else{
				lastseqint = fqseqint[0] >> 48;
				lastseqint &= maskl[8];
				igshift = Shiftseq(lastseqint, spos, shiftsearch, readlen, effl, seqint);
				if(igshift > -9){
					spos = epos -igshift;
				}
			}
		}
		if(spos != epos){
			for(j = 0; j < 20; j++){
				forw[j] = -1;
				back[j] = -1;
			}
			nforw = 0;
			nback = 0;
			kmax = 0;
			imisp = 0;
			ntotmis = 100;
			kg = -1;
			jg = -1;
			deltap = (epos > spos) ? 0:(spos - epos);
			realgap = (epos > spos) ? (epos-spos) : (spos-epos);
			if((realgap > 0)&&(realgap <= readlen/4)){
			nmis = 0;
			if((spos%32) == 0)
				shift = 62;
			else
				shift = 2*((spos%32) - 1);
			for(j = 0; (j <= (readlen-1)/32)&&(nmis <= *nmin); j++){
				if((spos%32) == 1)
					seqdot[j] = fqseqint[j]^seqint[spos/32+j];
				else
					seqdot[j] = fqseqint[j]^(((seqint[(spos-1)/32+j]<<shift)&(~maskl[shift/2]))|((seqint[(spos-1)/32+j+1]>>(64-shift))& maskl[shift/2]));
				if(j >= readlen/32)
					seqdot[j] &= ~maskl[32-readlen%32];
				for(p = 0; (p < 4) && (nmis <= *nmin); p++){
					bytemis = match_num[(seqdot[j]>>(16*(3-p)))&maskl[8]] >> 16;
					nmis += (bytemis > 3)? 1000:bytemis;
					kmax = bytemis;
					if(kmax > 4){
						kmax = 2;
					}
					for(k = 0; k < kmax; k++){
						forw[imisp] = (match_num[(seqdot[j]>>(16*(3-p)))&maskl[8]]>>(12-4*k))&0XF;
						forw[imisp] += 32*j + 8*p ;
						nforw ++;
						imisp ++;
					}
				}
			}
			imisp = 0;
			nmis = 0;
			if((epos%32) == 0)
				shift = 62;
			else
				shift = 2*((epos%32) - 1);
			for(j = (readlen-1)/32; (j >= 0)&&(nmis <= *nmin); j--){
				if((epos%32) == 1)
					seqdot[j] = fqseqint[j]^seqint[epos/32+j];
				else
					seqdot[j] = fqseqint[j]^(((seqint[(epos-1)/32+j]<<shift)&(~maskl[shift/2]))|((seqint[(epos-1)/32+j+1]>>(64-shift))& maskl[shift/2]));
				if(j >= readlen/32)
					seqdot[j] &= ~maskl[32-readlen%32];
				for(p = 3; (p >= 0) && (nmis <= *nmin); p--){
					bytemis = match_num[(seqdot[j]>>(16*(3-p)))&maskl[8]] >> 16;
					if(bytemis > 3)
						nmis = 1000;
					else
						nmis += bytemis;
					kmax = (bytemis > 4) ? 4:bytemis;
					for(k = (kmax-1); k >= 0; k--){
						back[imisp] = (match_num[(seqdot[j]>>(16*(3-p)))&maskl[8]]>>(12-4*k))&0XF;
						back[imisp] += j*32 + p*8;
						nback ++;
						imisp ++;
					}
					if(bytemis > 4)
						nback -= 2;
				}
			}
			for(j = 0; j < nforw; j++){
				for(k = 0; k < nback; k++){
					if(back[k] < (forw[j]+deltap)){
						if((j+k) <= ntotmis){
							if(((j+k) == ntotmis)&&(k > kg)){
								kg = k;
								jg = j;
							}
							if((j+k) < ntotmis){
								ntotmis = j+k;
								kg = k;
								jg = j;
							}
						}
					}
				}
			}
			if((jg > -1)&&(kg > -1)){
				if((forw[jg] <= jg*2)||(forw[jg] <= realgap/2)){
					if(kg > 3)
						kg = 3;
					deltap = 1+back[kg];
					realgap = deltap;
					ntotmis = kg;
				}
				else if(((readlen-back[kg]-1) <= kg*2)||((readlen-back[kg]-1) <= realgap/2)){
					back[kg] = readlen -1;
					if(jg > 3)
						jg = 3;
					realgap = readlen - forw[jg];
					deltap = realgap;
					ntotmis = jg;
				}
			}
			if((ntotmis + 1 + realgap/10) <= *nmin){ 
				if((*nbest == 0)||((ntotmis + 1 + realgap/10) < *nmin)){ 
					*nmin = ntotmis + 1 + realgap/10;
					*bestpos = spos;
					*nbest = 1;
					nmis = *nmin;
					if(back[kg] <= (deltap+1)){
						sprintf(cigar, "%dS%dM", (back[kg]+1), (readlen-1-back[kg]));
						*bestpos = epos + back[kg] + 1;
						*gaplen = back[kg]+1;
						*forwlen = 0;
						*backlen = readlen-1-back[kg];
					}
					else if (back[kg] < (readlen-1)){
						if(deltap == 0){
							*gaplen = -1*realgap;
							sprintf(cigar, "%dM%dD%dM", (back[kg]+1-deltap) , realgap, (readlen-1-back[kg]));
						}
						else{
							*gaplen = realgap;
							sprintf(cigar, "%dM%dI%dM", (back[kg]+1-deltap), realgap, (readlen-1-back[kg]));
						}
						
							*forwlen = back[kg]+1-deltap;
							*backlen = readlen-1-back[kg];
					}
					else{
						if(deltap == 0){
							sprintf(cigar, "%dM", readlen);
							*gaplen = 0;
						}
						else{
							sprintf(cigar, "%dM%dS", (readlen-realgap), realgap);
							*gaplen = realgap;
							*forwlen = readlen - realgap;
							*backlen = 0;
						}
					}
				}
				else if(abs(*bestpos-spos) > readlen/4){
					*nbest = 2;
				}
			}
			}
		}
		if(*nmin > 3){
			if((nmis > *nmin)&&(readlen > 47)){
				if((effmis <= 2)&&(effl >= readlen*8/10)){
					if(oldeffl >= 36){
						if(normis <= 2){
							if(oldeffl == readlen){
								nmis = normis;
								effl = readlen - readlen%8;
							}
							else{ 
								nmis = normis + (readlen - oldeffl + 8)/10;
								effl = oldeffl - 8;
							}
						}
						else{
							if(oldeffl == readlen){
								nmis = effmis + (readlen%8 + 8)/10;
								effl = readlen - 8 - readlen%8;
							}
							else{
								nmis = effmis + (readlen - oldeffl + 16)/10;
								effl = oldeffl - 16;
							}
						}
						nmis ++;
						if(nmis < *nmin){
							*nmin = nmis;
							sprintf(cigar, "%dM%dS", effl, (readlen -effl));
							*bestpos = seedpos[seedind[fqseedint[ipos]]+i]-ipos;
							*nbest = 1;
							*gaplen = readlen-effl;
							*forwlen = effl;
							*backlen = 0;
						}
					}
					else{
						if(normis <= 2){
							if(effl == readlen)
								nmis = normis;
							else
								nmis = normis + (readlen - effl + 8)/10;
							effl -= 8;
						}
						else{
							if(effl == readlen)
								nmis = effmis + 1;
							else
								nmis = effmis + (readlen - effl + 16)/10;
							effl -= 16;
						}
						nmis ++;
						if(nmis < *nmin){
							*nmin = nmis;
							sprintf(cigar, "%dS%dM", (readlen-effl), effl);
							*bestpos = seedpos[seedind[fqseedint[ipos]]+i]-ipos+readlen-effl;
							*nbest = 1;
							*gaplen = readlen-effl;
							*forwlen = 0;
							*backlen = effl;
						}
					}
				}
			}
		}
		}
	}

}

/* search for mapping results for single end reads in buffer */
void* aln_se(void *marg){
	args * data = (args*) marg;
	int maxseqlen = 500;
	int readlen = 500;
	char fqid[maxseqlen];
	char fqseq[maxseqlen];
	char revfqseq[maxseqlen];
	char fqual[maxseqlen];
	char revfqual[maxseqlen];
	uint64_t fqseqint [20];
	uint64_t revfqseqint [20];
	uint64_t seqdot [20];
	uint64_t revseqdot [20]; 
	uint32_t fqseedint [500];
	uint32_t revfqseedint [500];
	uint32_t nmapped = 0;
	uint32_t nunimapped = 0;
	uint32_t ntotfq = 0;
	uint32_t nmpos = 100;
	int iseed = 0;
	int ifq = 0;
	int nmin = readlen/10+1;
	int oldnmin = nmin;
	int nbest = 0;
	int direction = 0;
	int score = 0;
	int minpos1 = 100000;
	int minpos2 = 100000;
	int minpos3 = 100000;
	int minpos4 = 100000;
	uint32_t bestpos = 0;
	char cigar[20];

	int ipos[40];
	bool iscon = false;
	bool isunmap = false;

	int i, j, jc, p, nsfq, isbg, gaplen, forwlen, backlen, found, numn;
	for(p = (data->ip); p < (data->ip+data->bufsize); p++){
		(data->resbuffer+p)->flag = 4;
		gaplen = 0;
		forwlen = 0;
		backlen = 0;
		numn = 0;
		//(data->resbuffer+p)->iscomp = false;
		//(data->resbuffer+p)->isprint = false;
		iseed = 0;
		for(i = 0; i < 20; i++){
			ipos[i] = -1;
		}
		readlen = strlen((data->seqbuffer+p)->seq);
		strcpy(fqid, (data->seqbuffer+p)->sid);
		strcpy(fqseq, (data->seqbuffer+p)->seq);
		strcpy(fqual, (data->seqbuffer+p)->qual);
		for(i = 0; i < (readlen/32+1); i++){
			fqseqint[i] = 0;
			revfqseqint[i] = 0;
		}
		for(i = 0; i < readlen; i++){
			if((fqseq[i] == 'N')||(fqseq[i] == 'n')){
				numn ++;
			}
			ifq = i/32;																	
			fqseqint[ifq] = (fqseqint[ifq] << 2)|(c_to_i(fqseq + i) & 3);
			if(i >= 14)
				iseed ++;
			if(iseed == 0)
				fqseedint[iseed] = mask[7] & ((fqseedint[iseed] << 2)|(c_to_i(fqseq + i) & 3));
			else
				fqseedint[iseed] = mask[7] & ((fqseedint[iseed-1] << 2)|(c_to_i(fqseq + i) & 3));
			revfqseq[i] = rev(fqseq+readlen-1-i);
			revfqual[i] = *(fqual+readlen-1-i);
		}
		
	if(numn < readlen/8){
		if(readlen%32 != 0)
			fqseqint[readlen/32] = fqseqint[readlen/32] << (2*(32-readlen%32));

		iseed = 0;
		for(i = 0; i < readlen; i++){
			ifq = i/32;								 
			revfqseqint[ifq] = (revfqseqint[ifq] << 2)|(c_to_i(revfqseq + i) & 3);
			if(i >= 14)
				iseed ++;
			if(iseed == 0)
				revfqseedint[iseed] = mask[7] & ((revfqseedint[iseed] << 2)|(c_to_i(revfqseq + i) & 3));
			else
				revfqseedint[iseed] = mask[7] & ((revfqseedint[iseed-1] << 2)|(c_to_i(revfqseq + i) & 3));
		}
		if(readlen%32 != 0)
			revfqseqint[readlen/32] = revfqseqint[readlen/32] << (2*(32-readlen%32));
		minpos1 = 8000000;
		minpos2 = 8000000;
		for(i = 0; i <= (readlen-14); i++){
			if(data->seednum[fqseedint[i]] > 0){
				if(i % 2 == 0){
					if(data->seednum[fqseedint[i]] < minpos1){
						minpos1 = data->seednum[fqseedint[i]];
						ipos[0] = i;
					}
				}
				else{
					if(data->seednum[fqseedint[i]] < minpos2){
						minpos2 = data->seednum[fqseedint[i]];
						ipos[2] = i;
					}
				}
			}
		}

		minpos1 = 8000000;
		minpos2 = 8000000;
		for(i = 0; i <= (readlen-14); i++){
			if(data->seednum[revfqseedint[i]] > 0){
				if(i % 2 == 0){
					if(data->seednum[revfqseedint[i]] < minpos1){
						minpos1 = data->seednum[revfqseedint[i]];
						ipos[1] = i;
					}
				}
				else{
					if(data->seednum[revfqseedint[i]] < minpos2){
						minpos2 = data->seednum[revfqseedint[i]];
						ipos[3] = i;
					}
				}
			}
		}
		nmin = readlen/10+1;
		*cigar = '*';
		nbest = 0;
		for(i = 0; i < 4; i++){
			if(ipos[i] >= 0){
				oldnmin = nmin;
				if(i%2 == 0){
					Gapalign(ipos[i], readlen, fqid, fqseq, fqseqint, fqseedint, data->seqint, data->seednum, data->seedind, data->seedpos, data->chrpos, data->nchr, &nmin, &bestpos, &nbest, &gaplen, &forwlen, &backlen, cigar);
				}
				else{
					Gapalign(ipos[i], readlen, fqid, revfqseq, revfqseqint, revfqseedint, data->seqint, data->seednum, data->seedind, data->seedpos, data->chrpos, data->nchr, &nmin, &bestpos, &nbest, &gaplen, &forwlen, &backlen, cigar);
				}
				if(nmin < oldnmin){
					if(i%2 == 0)
						direction = 1;
					else
						direction = -1;
					oldnmin = nmin;
				}
			}
		}
		
		iscon = true;
		if(nmin == 0)
			iscon = false;
		else if((nbest > 1)&&(nmin == 1))
			iscon = false;

		if(iscon){
			oldnmin = nmin;
			nsfq = (readlen > 75) ? 4 : 3;
			if(readlen < 45)
				nsfq = 2;
			isbg = (readlen > nsfq*14) ? (readlen/nsfq-14) : 0;
			for(i = 0; i < nsfq; i++){
				minpos1 = 1000;
				minpos2 = 1000;
				minpos3 = 1000;
				minpos4 = 1000;
				for(j = i*readlen/nsfq; j <= ((i+1)*readlen/nsfq -14); j++){
					if((data->seednum[fqseedint[j]] > 0)&&(data->seednum[fqseedint[j]] < 500)){
						if(j%2 == 0){
							if(data->seednum[fqseedint[j]] < minpos1){
								minpos1 = data->seednum[fqseedint[j]];
								ipos[4+4*i] = j;
							}
						}
						else{
							if(data->seednum[fqseedint[j]] < minpos2){
								minpos2 = data->seednum[fqseedint[j]];
								ipos[6+4*i] = j;
							}
						}
					}
					if((data->seednum[revfqseedint[j]] > 0)&&(data->seednum[revfqseedint[j]] < 500)){
						if(j%2 == 0){
							if(data->seednum[revfqseedint[j]] < minpos3){
								minpos3 = data->seednum[revfqseedint[j]];
								ipos[5+4*i] = j;
							}
						}
						else{
							if(data->seednum[revfqseedint[j]] < minpos4){
								minpos4 = data->seednum[revfqseedint[j]];
								ipos[7+4*i] = j;
							}
						}
					}
				}
				if((minpos1 > 500)||(minpos2 > 500)){
					for(j = (isbg+i*readlen/nsfq); j <= ((i+1)*readlen/nsfq-8); j++){
						if((data->seednum[fqseedint[j]] > 0)&&(data->seednum[fqseedint[j]] < 500)){
							if((j%2 == 0)&&(minpos1 > 500)){
								minpos1 = data->seednum[fqseedint[j]];
								ipos[4+4*i] = j;
							}
							if((j%2 == 1)&&(minpos2 > 500)){
								minpos2 = data->seednum[fqseedint[j]];
								ipos[6+4*i] = j;
							}
						}
					}
				}
				if((minpos3 > 500)||(minpos4 > 500)){
					for(j = (isbg+i*readlen/nsfq); j <= ((i+1)*readlen/nsfq-8); j++){
						if((data->seednum[revfqseedint[j]] > 0)&&(data->seednum[revfqseedint[j]] < 500)){
							if((j%2 == 0)&&(minpos3 > 500)){
								minpos3 = data->seednum[revfqseedint[j]];
								ipos[5+4*i] = j;
							}
							if((j%2 == 1)&&(minpos4 > 500)){
								minpos4 = data->seednum[revfqseedint[j]];
								ipos[7+4*i] = j;
							}
						}
					}
				}

				if(ipos[4+4*i] == ipos[0]){
					ipos[4+4*i] = -1;
				}
				if(ipos[5+4*i] == ipos[1]){
					ipos[5+4*i] = -1;
				}
				if(ipos[6+4*i] == ipos[2]){
					ipos[6+4*i] = -1;
				}
				if(ipos[7+4*i] == ipos[3]){
					ipos[7+4*i] = -1;
				}
			}
			
			/*
			for(i = 0; i < (nsfq-1); i++){
				if((ipos[4+4*i]+ipos[6+4*i]+ipos[5+4*i]+ipos[7+4*i]) > (ipos[8+4*i]+ipos[10+4*i]+ipos[9+4*i]+ipos[11+4*i])){
					ipos[0] = ipos[4+4*i];
					ipos[2] = ipos[6+4*i];
					ipos[4+4*i] = ipos[8+4*i];
					ipos[6+4*i] = ipos[10+4*i];
					ipos[8+4*i] = ipos[4+4*i];
					ipos[10+4*i] = ipos[6+4*i];
					ipos[1] = ipos[5+4*i];
					ipos[3] = ipos[7+4*i];
					ipos[5+4*i] = ipos[9+4*i];
					ipos[7+4*i] = ipos[11+4*i];
					ipos[9+4*i] = ipos[5+4*i];
					ipos[11+4*i] = ipos[7+4*i];
				}
			}
			*/
			found = 0;
			for(i = 4; (i < (nsfq+1)*4)&&(found < 1); i++){
				if(ipos[i] >= 0){
					if(i%2 == 0){
						Gapalign(ipos[i], readlen, fqid, fqseq, fqseqint, fqseedint, data->seqint, data->seednum, data->seedind, data->seedpos, data->chrpos, data->nchr, &nmin, &bestpos, &nbest, &gaplen, &forwlen, &backlen, cigar);
					}
					else{
						Gapalign(ipos[i], readlen, fqid, revfqseq, revfqseqint, revfqseedint, data->seqint, data->seednum, data->seedind, data->seedpos, data->chrpos, data->nchr, &nmin, &bestpos, &nbest, &gaplen, &forwlen, &backlen, cigar);
					}
					if(nmin < oldnmin){
						if(i%2 == 0)
							direction = 1;
						else
							direction = -1;
						oldnmin = nmin;
					}
				
					if(nmin < 3){
						if(nmin == 1){
							if((nbest > 1)||(i >= 7))
								found = 1;
						}
						else if(i >= 11)
							found = 1;
					}
				}
			}
		}

		if(nmin <= readlen/10){
			nmapped++;
			if(direction == 1)
				(data->resbuffer+p)->flag = 0;
			else{
				(data->resbuffer+p)->flag = 16;
				strcpy((data->seqbuffer+p)->seq, revfqseq);
				strcpy((data->seqbuffer+p)->qual, revfqual);
			}
			isunmap = false;
			(data->resbuffer+p)->pos = bestpos;
			for(jc = 0; jc < data->nchr; jc++){
				if(bestpos < data->chrpos[jc]){
					if(jc > 0){
						(data->resbuffer+p)->jc = jc;
						(data->resbuffer+p)->pos -= data->chrpos[jc-1];
					}
					else
						(data->resbuffer+p)->jc = 0;
					if((bestpos + 5*readlen/4+1) > data->chrpos[jc]){
						if(gaplen == 0){
							if((bestpos + readlen*3/4-1) > data->chrpos[jc]){
								if((bestpos + readlen/4) > data->chrpos[jc]){
									sprintf(cigar, "%dS%dM", (data->chrpos[jc]+1-bestpos), (bestpos+ readlen -1-data->chrpos[jc]));
									(data->resbuffer+p)->jc = jc +1;
									(data->resbuffer+p)->pos = 1;
								}
								else
									isunmap = true;
							}
							else if((bestpos + readlen -1) > data->chrpos[jc]){
								sprintf(cigar, "%dM%dS", (data->chrpos[jc]+1-bestpos), (bestpos+ readlen -1-data->chrpos[jc]));
							}
						}
						else if(forwlen == 0){
							if((backlen+bestpos-1) > data->chrpos[jc]){
								if((backlen+bestpos-1-readlen/4) < data->chrpos[jc]){
									sprintf(cigar, "%dS%dM%dS", (readlen-backlen), (data->chrpos[jc]+1-bestpos), (backlen+bestpos-1-data->chrpos[jc]));
								}
								else if((bestpos+readlen/4) > data->chrpos[jc]){
									sprintf(cigar, "%dS%dM", (data->chrpos[jc]+1-bestpos+gaplen), (readlen+bestpos-1-gaplen-data->chrpos[jc]));
									(data->resbuffer+p)->jc = jc +1;
									(data->resbuffer+p)->pos = 1;
								}
								else
									isunmap = true;
							}
						}
						else if(backlen == 0){
							if((forwlen+bestpos-1) > data->chrpos[jc]){
								if((forwlen+bestpos-1-readlen/4) < data->chrpos[jc]){
									sprintf(cigar, "%dM%dS", (data->chrpos[jc]+1-bestpos), (readlen-1+bestpos-data->chrpos[jc]));
								}
								else
									isunmap = true;
							}
						}
						else if(gaplen< 0){
							if((bestpos+readlen-1-gaplen) > data->chrpos[jc]){
								if((bestpos+readlen-1-gaplen-readlen/4) < data->chrpos[jc]){
									if(backlen > (bestpos+readlen+3-gaplen-data->chrpos[jc]))
										sprintf(cigar, "%dM%dD%dM%dS", forwlen, -1*gaplen, (backlen-bestpos-readlen+1+gaplen+data->chrpos[jc]), (bestpos+readlen-1-gaplen-data->chrpos[jc]));
									else if((bestpos+forwlen-1) <= data->chrpos[jc])
										sprintf(cigar, "%dM%dS", forwlen, (readlen-forwlen));
									else if((bestpos+3*readlen/4) < data->chrpos[jc])
										sprintf(cigar, "%dM%dS", (data->chrpos[jc]+1-bestpos), (readlen-1+bestpos-data->chrpos[jc]));
									else
										isunmap = true;
								}
								else
									isunmap = true;
							}
						}
						else{
							if((bestpos+readlen-1-gaplen) > data->chrpos[jc]){
								if((bestpos+readlen-1-gaplen-readlen/4) < data->chrpos[jc]){
									if(backlen > (bestpos+readlen+3-gaplen-data->chrpos[jc]))
										sprintf(cigar, "%dM%dI%dM%dS", forwlen, gaplen, (backlen-bestpos-readlen+1+gaplen+data->chrpos[jc]), (bestpos+readlen-1-gaplen-data->chrpos[jc]));
									else if((bestpos+forwlen-1) <= data->chrpos[jc])
										sprintf(cigar, "%dM%dS", forwlen, (readlen-forwlen));
									else if((bestpos+3*readlen/4) < data->chrpos[jc])
										sprintf(cigar, "%dM%dS", (data->chrpos[jc]+1-bestpos), (readlen-1+bestpos-data->chrpos[jc]));
									else
										isunmap = true;
								}
								else
									isunmap = true;
							}
						}
						if(isunmap){
							(data->resbuffer+p)->flag = 4;
							nbest = 0;
							score = 0;
							direction = 0;
							nmin = readlen;
							nmapped --;
						}
					}
					break;
				}
				else if(bestpos == data->chrpos[jc]){
					(data->resbuffer+p)->pos = 1;
					(data->resbuffer+p)->jc = jc + 1;
					if(gaplen == 0){
						sprintf(cigar, "1S%dM", (readlen-1));
					}
					else if(forwlen == 0){
						sprintf(cigar, "%dS%dM", (readlen+1-backlen), (backlen-1));
					}
					else if(backlen == 0){
						sprintf(cigar, "1S%dM%dS", (forwlen-1), (readlen-forwlen));
					}
					else if(gaplen < 0){
						if(forwlen > 3)
							sprintf(cigar, "1S%dM%dD%dM", (forwlen-1), -1*gaplen, backlen);
						else{
							sprintf(cigar, "%dS%dM", (readlen-backlen), backlen);
							(data->resbuffer+p)->pos = forwlen-gaplen;
						}
					}
					else{
						if(forwlen > 3)
							sprintf(cigar, "1S%dM%dI%dM", (forwlen-1), gaplen, backlen);
						else{
							sprintf(cigar, "%dS%dM", (readlen-backlen), backlen);
							(data->resbuffer+p)->pos = forwlen;
						}
					}	
					break;
				}
			}
														
			if(nbest == 1){
				if(nmin < 3)
					score = readlen;
				else
					score = readlen*readlen/(nmin*nmin*16);
				if(score > 60)
					score = 60;
			}
			else if(nbest > 1)
				score = 1;
			else
				score = 0;
			(data->resbuffer+p)-> mapscore = score;
			if(score > 1)
				nunimapped ++;
			
			if(*cigar == '*')
				sprintf((data->resbuffer+p)->cigar, "%dM", readlen);
			else
				strcpy((data->resbuffer+p)->cigar, cigar);
			//printf("nmin %d gaplen %d  %d %d \n", nmin, gaplen, forwlen, backlen);
			if(gaplen == 0)
				(data->resbuffer+p)->nmis = nmin;
			else{
				if((forwlen != 0)&&(backlen != 0)){
					(data->resbuffer+p)->nmis = nmin -1 - abs(gaplen)/10 + abs(gaplen);
				}
				else{
					(data->resbuffer+p)->nmis = nmin -1 - abs(gaplen)/10;
				}
			}
			if((data->resbuffer+p)->nmis < 0)
				(data->resbuffer+p)->nmis = nmin;
		}
	}
		(data->resbuffer+p)->iscomp = true;
	}
	data->nmapped = nmapped;
	data->nunimapped = nunimapped;
	data->iscomp = 1;
}

void* dataprint(void* arg){
	dargs* data = (dargs*) arg;
	int i, j, jc;
	
	char psid[256];
	
	while(!isend){
		for(j = 0; j < nthreads; j++){
			if(((data->arglist+j)->resbuffer+(data->arglist+j)->ip)->iscomp && (!(data->arglist+j)->isprint)){	
				for(i = (data->arglist+j)->ip ; i < ((data->arglist+j)->ip+(data->arglist+j)->bufsize); i++){
					if(((data->arglist+j)->resbuffer+i)->iscomp && (!((data->arglist+j)->resbuffer+i)->isprint)){
							if(isspace){
								sscanf(((data->arglist+j)->seqbuffer+i)->sid, "%s", psid);
								printf("%s\t", (psid+1));
							}
							else if(!isslash){
								printf("%s\t", (((data->arglist+j)->seqbuffer+i)->sid+1));
							}
							else{
								printf("%*.*s\t", 2, (strlen(((data->arglist+j)->seqbuffer+i)->sid)-3), (((data->arglist+j)->seqbuffer+i)->sid+1));
							}
								
							
								if(((data->arglist+j)->resbuffer+i)->flag == 4)
									printf("4\t*\t0\t0\t*\t*\t0\t0\t%s\t%s", ((data->arglist+j)->seqbuffer+i)->seq, ((data->arglist+j)->seqbuffer+i)->qual);
								else{
									printf("%d\t", ((data->arglist+j)->resbuffer+i)->flag);
									printf("%s%lu\t", chrlist[((data->arglist+j)->resbuffer+i)->jc], ((data->arglist+j)->resbuffer+i)->pos);
									
									printf("%d\t", ((data->arglist+j)->resbuffer+i)->mapscore);
									printf("%s\t", ((data->arglist+j)->resbuffer+i)->cigar);
									printf("*\t0\t0\t%s\t%s\t", ((data->arglist+j)->seqbuffer+i)->seq, ((data->arglist+j)->seqbuffer+i)->qual);
									printf("NM:i:%d", ((data->arglist+j)->resbuffer+i)->nmis);
									//printf("AS:i:%d", -1*(((data->arglist+j)->resbuffer+i)->nmis));
								}
								if(isid == 1){
									printf("\tRG:Z:%s\n", idname);
								}
								else
									printf("\n");
								((data->arglist+j)->resbuffer+i)->isprint = true;
								(data->arglist+j)->nprint ++;
					}
				}
				if((data->arglist+j)->nprint == (data->arglist+j)-> bufsize){
					(data->arglist+j)->isprint = true;
					totmapped += (data->arglist+j)->nmapped;
					totunimapped += (data->arglist+j)->nunimapped;
					totnprint += (data->arglist+j)->nprint;
				}
			}
		}
		if(pret == NULL){
			if(totnprint == totnfq){
				isend = true;
			}
		}
	}						
}
			
int main(int argv, char* argc[])
{
	struct timeval begin;
	struct timeval end;
	

	int nline = 0;
	int i = 0;
	int j = 0;
	int k = 0;
	int m = 0;
	int p = 1;
	int nmis = 0;
	int found = 0;
	initial_match_num();

	//int maxnchr = 10000;
	//int maxchrlen = 1000;
	//char chrlist[maxnchr][maxchrlen];
	uint32_t nchr = 93;
	uint32_t nns = 780; // n N seg * 2
	uint32_t reflength = 3137161264;
	uint32_t nseqint = 98036290;
	uint32_t maxnseeds = 268435456;
	uint32_t totpos = 1437553481;
	subsize = 50000;
	nthreads = 1;
	buffersize = subsize*nthreads;
	size_t seqlen = 50;
	ssize_t linelength;
	int readlen = 0;
	int sidlen = 0;
	//int infq = -1;
	char firstseq[500];
	char firstqual[500];
	char firstid[500];
	char plus[500];
	char plname[200];
	char lbname[200];
	char smname[200];
	int issm = 0;
	int ispl = 0;
	int islb = 0;
	//int totmapped = 0;
	uint32_t iseq = 0;
	int buflen = 0;

	char refpath[5000];
	char chrpath[5000];
	char infopath[5000];
	char pospath[5000];
	char snpath[5000];
	char sipath[5000];

	if(argv > 1){
		sprintf(chrpath, "%s.chrs", argc[1]);
		sprintf(infopath, "%s.info", argc[1]);
		sprintf(pospath, "%s.spos", argc[1]);
		sprintf(snpath, "%s.sn", argc[1]);
		sprintf(sipath, "%s.si", argc[1]);
		sprintf(refpath, "%s.cref", argc[1]);
	}
	else{
		fprintf(stderr, "usage: \n");
		fprintf(stderr, "saligner index reads1.fq -n nthreads -id sample_id -sm sample_name -pl platform_name -lb library_name > results.sam \n");
		exit(0);
	}

	if(argv > 4){
		for(i = 3; i < argv; i++){
			if(strcmp(argc[i], "-n") == 0)
				nthreads = atoi(argc[i+1]);				
			if(strcmp(argc[i], "-b") == 0)
				subsize = atoi(argc[i+1]);				
			if(strcmp(argc[i], "-id") == 0){
				strcpy(idname, argc[i+1]);	
				isid = 1;
			}
			if(strcmp(argc[i], "-sm") == 0){
				strcpy(smname, argc[i+1]);
				issm = 1;
			}
			if(strcmp(argc[i], "-pl") == 0){
				strcpy(plname, argc[i+1]);
				ispl = 1;
			}
			if(strcmp(argc[i], "-lb") == 0){
				strcpy(lbname, argc[i+1]);
				islb = 1;
			}
		}
		buffersize = subsize*nthreads;
	}

	FILE *fn = fopen(infopath,"r");
	if(fn == NULL){
		fprintf(stderr, " read info file failed \n");
		exit(0);
	}
	else{
		fread(&nchr, sizeof(uint32_t), 1, fn);
		fread(&nns, sizeof(uint32_t), 1, fn);
		fread(&nseqint, sizeof(uint32_t), 1, fn);
		fread(&totpos, sizeof(uint32_t), 1, fn);		
	}
	maxreflength = nseqint*32;
	
	char* seq = (char*)calloc(5000, sizeof(char));
	uint64_t* seqint = (uint64_t*) calloc(nseqint, sizeof(uint64_t));
	uint16_t* seednum = (uint16_t*) calloc(maxnseeds, sizeof(uint16_t));
	uint32_t* seedind = (uint32_t*) calloc(maxnseeds, sizeof(uint32_t));
	uint32_t* seedpos = (uint32_t*) calloc(totpos, sizeof(uint32_t));
	uint32_t* chrlen = (uint32_t*) calloc(nchr, sizeof(uint32_t));
	uint32_t* chrpos = (uint32_t*) calloc(nchr, sizeof(uint32_t));
	uint32_t* npos = (uint32_t*) calloc(nns, sizeof(uint32_t)); 

	if(fn != NULL){
		fread(chrlen, sizeof(uint32_t), nchr, fn);
		fread(chrpos, sizeof(uint32_t), nchr, fn);
		fread(npos, sizeof(uint32_t), nns, fn);
	}
	fclose(fn);

	FILE *fchr = fopen(chrpath,"r");

	i = 0;
	if(fchr == NULL){
		printf(" read chrs file failed \n");
		exit(0);
	}
	else{
		while((linelength = getline(&seq, &seqlen, fchr)) != -1){
			strcpy(chrlist[i], seq);
			chrlist[i][strlen(chrlist[i])-1] = '\t';
			i++;
		}
	}
	fclose(fchr);

	FILE *fseqint = fopen(refpath,"r");
	if(fseqint == NULL){
		printf(" read seq int file failed \n");
		exit(0);
	}
	else
		fread(seqint, sizeof(uint64_t), nseqint, fseqint);

	fclose(fseqint);

	FILE *fseednum = fopen(snpath,"r");
	if(fseednum == NULL){
		printf(" read seed num file failed \n");
		exit(0);
	}
	else
		fread(seednum, sizeof(uint16_t), maxnseeds, fseednum);

	fclose(fseednum);

	FILE *fseedind = fopen(sipath,"r");
	if(fseedind == NULL){
		printf(" read seed ind file failed \n");
		exit(0);
	}
	else
		fread(seedind, sizeof(uint32_t), maxnseeds, fseedind);

	fclose(fseedind);

	FILE *fseedpos = fopen(pospath,"r");
	if(fseedpos == NULL){
		printf(" read seed pos file failed \n");
		exit(0);
	}
	else
		fread(seedpos, sizeof(uint32_t), totpos, fseedpos);

	fclose(fseedpos);

	gettimeofday(&begin, NULL);
	
	FILE *fq = NULL;
	gzFile gzfp;
	bool isgz = false;
	int ret = 0;
	unsigned gzbufsize = 256*1024;
	
	if (strcmp(argc[2] + strlen(argc[2]) - 3, ".gz") == 0) {
		isgz = true;
		gzfp = gzopen(argc[2], "r");
		gzbuffer(gzfp, gzbufsize);
		if(!gzfp) {
			fprintf(stderr, "read fastq file err \n");
			exit(0);
		}
		if (!gzgets(gzfp, firstid, 512)) {
			fprintf(stderr, "read fastq file err \n");
			exit(0);
		}
		firstid[strlen(firstid) - 1] = '\0';

		gzgets(gzfp, firstseq, 512);
		firstseq[strlen(firstseq) - 1] = '\0';
		gzgets(gzfp, plus, 512);
		plus[strlen(plus) - 1] = '\0';
		gzgets(gzfp, firstqual, 512);
		firstqual[strlen(firstqual) - 1] = '\0';
	} else {
		fq = fopen(argc[2],"r");

		if(fq == NULL){
			fprintf(stderr, "read fastq file err \n");
			exit(0);
		}

		if((pret = fgets(firstid, 512, fq)) == NULL){
			fprintf(stderr, "read fastq file err \n");
			exit(0);
		}
		firstid[strlen(firstid) - 1] = '\0';

		fgets(firstseq, 512, fq);
		firstseq[strlen(firstseq) - 1] = '\0';
		fgets(plus, 512, fq);
		plus[strlen(plus) - 1] = '\0';
		fgets(firstqual, 512, fq);
		firstqual[strlen(firstqual) - 1] = '\0';
	}
	readlen = strlen(firstseq);
	sidlen = strlen(firstid);
	iseq = 1;
	totnfq = 1;
	if(readlen > 17){
		subsize = 10000*50/readlen;
		buffersize = subsize*nthreads;
		for(i = 0; (i < sidlen)&&(!isspace); i++){
			if((firstid[i] == ' ')||(firstid[i] == '\t')){
				isspace = true;
			}
		}
		if(sidlen > 2){
			if(firstid[sidlen-2] == '/')
				isslash = true;
		}
	}
	else{
		fprintf(stderr, "read fastq file err \n");
		exit(0);
	}
	for(i = 0; i < nchr; i++){
		printf("@SQ\tSN:%sLN:%lu\n", chrlist[i],chrlen[i]);
	}
	if(isid == 1){
		printf("@RG\tID:%s", idname);
		if(issm == 1)
			printf("\tSM:%s", smname);
		if(ispl == 1)
			printf("\tPL:%s", plname);
		if(islb == 1)
			printf("\tLB:%s", lbname);
		printf("\n");
	}

	reads* seqbuffer = (reads*) calloc(buffersize, sizeof(reads));
	res* resbuffer = (res*) calloc(buffersize, sizeof(res));
	pthread_t tid[nthreads+1];
	ret = -1;
	
	int ithread = 0;
	int fthreads = 0;
	int jc = 0;
	bool isfinish = false;
	bool ispc = false;

	for(i = 0; (readlen > 0)&&(i < buffersize); i++){
		(seqbuffer+i)->sid = (char*)calloc((sidlen+50), sizeof(char));
		(seqbuffer+i)->seq = (char*)calloc((readlen+1), sizeof(char));
		(seqbuffer+i)->qual = (char*)calloc((readlen+1), sizeof(char));
		(resbuffer+i)->flag = 4;
		(resbuffer+i)->iscomp = false;
		(resbuffer+i)->isprint = false;
	}


	args* arglist = (args*)calloc(nthreads, sizeof(args));
	for(i = 0; i < nthreads; i++){
		(arglist+i)->seqbuffer = seqbuffer;
		(arglist+i)->resbuffer = resbuffer;
		(arglist+i)->ip = i*subsize;
		(arglist+i)->bufsize = subsize;
		(arglist+i)->seqint = seqint;
		(arglist+i)->seednum = seednum;
		(arglist+i)->seedind = seedind;
		(arglist+i)->seedpos = seedpos;
		(arglist+i)->chrpos = chrpos;
		(arglist+i)->nchr = nchr;
		(arglist+i)->iscomp = 0;
		(arglist+i)->isprint = 1;
	}
	
	dargs dataarg;
	dataarg.arglist = arglist;
		
	strcpy(seqbuffer->sid, firstid);
	strcpy(seqbuffer->seq, firstseq);
	strcpy(seqbuffer->qual, firstqual);
	while(!isend){
		if(ithread != -1){
		if(isgz){
			if ((pret = gzgets(gzfp, (seqbuffer+iseq)->sid, 512)) != NULL) {
				(seqbuffer+iseq)->sid[strlen((seqbuffer+iseq)->sid) - 1] = '\0';
				gzgets(gzfp, (seqbuffer+iseq)->seq, 512);
				(seqbuffer+iseq)->seq[strlen((seqbuffer+iseq)->seq) - 1] = '\0';
				gzgets(gzfp, plus, 512);
				plus[strlen(plus) - 1] = '\0';
				gzgets(gzfp, (seqbuffer+iseq)->qual, 512);
				(seqbuffer+iseq)->qual[strlen((seqbuffer+iseq)->qual) - 1] = '\0';
				iseq++;
				totnfq++;
			}
		}
		else{
			if((pret = fgets((seqbuffer+iseq)->sid, 512, fq)) != NULL){
				(seqbuffer+iseq)->sid[strlen((seqbuffer+iseq)->sid) - 1] = '\0';
				fgets((seqbuffer+iseq)->seq, 512, fq);
				(seqbuffer+iseq)->seq[strlen((seqbuffer+iseq)->seq) - 1] = '\0';
				fgets(plus, 512, fq);
				plus[strlen(plus) - 1] = '\0';
				fgets((seqbuffer+iseq)->qual, 512, fq);
				(seqbuffer+iseq)->qual[strlen((seqbuffer+iseq)->qual) - 1] = '\0';
				iseq++;
				totnfq++;
			}
		}
		}
		if((iseq%subsize == 0)||(pret == NULL)){
			if(ithread != -1){
			(arglist+ithread)->ip = ithread*subsize;
			(arglist+ithread)->bufsize = iseq - ithread*subsize;
			(arglist+ithread)->iscomp = 0;
			(arglist+ithread)->isprint = 0;
			(arglist+ithread)->nprint = 0;
			for(i = (arglist+ithread)->ip; i < ((arglist+ithread)->ip + (arglist+ithread)-> bufsize); i++){
				(resbuffer+i)->iscomp = false;
				(resbuffer+i)->isprint = false;
			}
			ret = pthread_create((tid+ithread), NULL, aln_se, (arglist+ithread));
			if(ret != 0){
				printf("Thread Create Error\n");
				exit(0);
			}
			ithread = -1;
			}
			/*
			for(j = 0; j < nthreads; j++){
				if((arglist+j)->iscomp){
					pthread_join(tid[j], NULL);
					//totmapped += (arglist+j)->nmapped;
					(arglist+j)->iscomp = false;
				}
			}
			*/			
			if(pret != NULL){
				for(j = 0; (j < nthreads)&&(ithread == -1); j++){
						if((arglist+j)->isprint){
								ithread = j;
								iseq = ithread*subsize;
						}
						if((arglist+j)->iscomp){
							pthread_join(tid[j], NULL);
							(arglist+j)->iscomp = false;
						}
				}
			}
		
				if(!ispc){
					ret = pthread_create((tid+nthreads), NULL, dataprint, &dataarg); 
					if(ret != 0){
						fprintf(stderr, "Thread Create Error\n");
						exit(0);
					}
					else{
						ispc = true;
					}
				}
			
			if(pret == NULL){
				for(j = 0; (j < nthreads); j++){
					if((arglist+j)->iscomp){
						pthread_join(tid[j], NULL);
						//totmapped += (arglist+j)->nmapped;
						(arglist+j)->iscomp = false;
					}
				}
				if(isend)
					pthread_join(tid[nthreads], NULL);
			}
		}
	}

	if(isgz){
		gzclose(gzfp);
	} 
	else{
		fclose(fq);
	}
	fprintf(stderr, "number of total seqs\t%d\n", totnfq);
	fprintf(stderr, "number of total mapped seqs\t%d\t%.2f\%\n", totmapped, 100*totmapped/(totnfq+0.00001));
	fprintf(stderr, "number of uniquely mapped seqs\t%d\t%.2f\%\n", totunimapped, 100*totunimapped/(totnfq+0.00001));

	for(i = 0; i < buffersize; i++){
			if(strlen((seqbuffer+i)->seq) > 0){
				free((seqbuffer+i)->sid);
				free((seqbuffer+i)->seq);
				free((seqbuffer+i)->qual);
			}
	}
	
	free(seqbuffer);
	free(resbuffer);
	free(arglist);
	free(seq);
	free(seqint);
	free(seednum);
	free(seedind);
	free(seedpos);
	free(npos);
	free(chrlen);
	free(chrpos);

	gettimeofday(&end, NULL);
	fprintf(stderr,"search cost = %fsec\n",end.tv_sec - begin.tv_sec + float(end.tv_usec - begin.tv_usec) / 1000000);

	return 0;
}

