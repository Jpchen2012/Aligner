/* 
Qalinger an ultrafast short read aligner
Copyright (C) 2016 BGI-shenzhen
Jianping Chen <chenjianping@genomics.cn>
Qaligner is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
Qaligner is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
You should have received a copy of the GNU General Public License along with this program. If not, see <http://www.gnu.org/licenses/>.
*/

/* This is a mapping program for pair ends short reads data produced by Hiseq or CPAS platform. */

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

int infq;
int subsize, buffersize;	
int nthreads;	
int isid;
int maxnmappos = 99;
uint64_t totnfq = 0;
uint64_t totnprint = 0;
uint64_t totmapped = 0;
uint64_t totunimapped = 0;
uint64_t totlunimapped = 0;
uint64_t totrunimapped = 0;
uint64_t totlmapped = 0;
uint64_t totrmapped = 0;
uint64_t maxreflength = 1024;												
char chrlist[10000][1000];
char idname[200];
char *pret = NULL;
bool isspace = false;
bool isslash = false;

typedef struct read{
			char *sid; 
			char *seq; 
			char *qual;
}reads;

typedef struct result{
		int flag;
		int jc;
		uint32_t pos;
		int deltapos;
		int mapscore;
		int nmis;
		int nipos;
		char cigar[20];
		bool iscomp;
}res;

typedef struct farg{
		reads* seqbuffer;
		res* resbuffer;
		reads* seqbufferb;
		res* resbufferb;
		int ip;
		int bufsize;
		int nmapped;
		int nunimapped;
		int nlmapped;
		int nrmapped;
		int nlunimapped;
		int nrunimapped;
		bool iscomp;
		int nprint;
		bool isprint;
		int isfw;
		uint64_t* seqint;
		uint16_t* seednum;
		uint32_t* seedind;
		uint32_t* seedpos;
		uint32_t* chrpos;
		uint32_t nchr;
}args;

typedef struct datarg{
	args* arglist;
	int readlen;
	int readlenb;
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
				if((shiftmis+abs(ig-gshift)/3) < (minsm+abs(igshift-gshift)/3)){
					igshift = ig;
					minsm = shiftmis;
				}
			}
		}
		if((minsm > 2)||(abs(igshift-gshift) > 6)){
			if((sreadlen-seffl) <= 6){
				for(int ig = gshift-6; ig <= gshift+6; ig++){
						if(direct == 1)
							shiftmis = match_num[(lastseqint^(targetint>>2*ig))&maskl[4]] >> 16;
						else
							shiftmis = match_num[((lastseqint^(targetint>>2*ig))>>8)&maskl[4]] >> 16;
						if(shiftmis <= 1){
							if((shiftmis <= minsm)&&(abs(ig-gshift) < abs(igshift-gshift))){
								igshift = ig;
								minsm = shiftmis;
							}
						}
				}
			}
		}
																						
		return (igshift-gshift);
}

/* search mapping results for a single reference position */ 
void Spalign(uint32_t ispos, int readlen, uint64_t* fqseqint, uint32_t* fqseedint, uint64_t* seqint, uint16_t* seednum, uint32_t* seedind, uint32_t* seedpos, uint32_t* chrpos, uint32_t nchr, int* nmin, uint32_t* bestpos, int* gaplen, int* forwlen, int* backlen, char* cigar){
		int i = 0;
		int j = 0;
		int k = 0;
		int m = 0;
		int p = 1;
		int nmis = 0;
		int found = 0;
		int forw[100];
		int back[100];
		uint32_t spos, epos, delta, gapl;
		int sp = 0;
		int ep = 0;
		int normis = 15;
		int bytemis = 100;
		int shift = 0;
		int begmis = 10;
		int effmis = 0;
		int termis = 0;
		int termisb = 0;
		int nmisth = 2;
		int effl = 0;
		int oldeffl = 0;
		int kgap = 0;
		int searchindel = 0;
		int shiftsearch = 0;
		int igshift = 0;
		int match = 2;
		int mis = 5;
		int gapo = 5;
		int gape = 1;
		int gapmis = 10;
		int nforw, nback, kmax, imisp, ntotmis, kg, jg, deltap, realgap, mscore, dclip, iclip;
		uint32_t nshiftmis = 0;
		uint32_t nshiftmis2 = 0;
		uint32_t nshiftshit = 0;
		uint32_t nshiftgood = 0;
		uint64_t seqdot[20];
												spos = ispos;
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
														if(nmis < *nmin){												
															*nmin = nmis;
															*bestpos = spos;
															*gaplen = 0;
															sprintf(cigar, "%dM", readlen);
														}														
												}
												else if((*nmin > 2)&&(readlen > 47)){
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
																		if(gapl == 0)
																			found = 0;
																		
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
																				if(gapl == 0)
																					found = 0;
																				if((found == 0)&&((readlen-effl) < 29))
																					shiftsearch = -1;
																		}
																}
														}
												}

												if((shiftsearch != 0)&&(readlen > 47)){
														uint64_t lastseqint = 0;
														uint32_t ankpos = spos+readlen-8;
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
												if((spos != epos )&&(spos < maxreflength)){
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
													if((realgap < readlen/4)&&(realgap > 0)){	
														nmis = 0;
														mscore = 0;
														iclip = -1;
														dclip = 0;
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
																				if(((forw[imisp]-imisp)*match-imisp*mis) > mscore){
																					mscore = (forw[imisp]-imisp)*match-imisp*mis;
																					iclip = imisp;
																					dclip = 1;
																				}
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
																				if(((readlen-back[imisp]-1-imisp)*match-imisp*mis) > mscore){
																					mscore = (readlen-back[imisp]-1-imisp)*match-imisp*mis;
																					iclip = imisp;
																					dclip = -1;
																				}
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
														mscore -= 2;
														if((jg > -1)&&(kg > -1)){
															if((dclip != 0)&&(((readlen-deltap-ntotmis)*match-ntotmis*mis-gapo-(realgap-1)*gape) > mscore)){
																iclip = -1;
																dclip = 0;
															}
															else if(ntotmis < 1){
																if((forw[jg] > (2+realgap/2))&&(forw[jg] < (readlen-(2+deltap+realgap/2)))){
																	if((back[kg] > (2+deltap+realgap/2))&&(back[kg] < (readlen-(2+realgap/2)))){
																		iclip = -1;
																		dclip = 0;
																	}
																}
															}
															
									
															if(dclip == 0){
																if((ntotmis + 1 + realgap/gapmis) < *nmin){ 
																	*nmin = ntotmis + 1 + realgap/gapmis;
																	*bestpos = spos;																								
																	nmis = *nmin;
																	if(back[kg] <= (deltap+1)){
																		sprintf(cigar, "%dS%dM", (back[kg]+1), (readlen-1-back[kg]));
																		*bestpos = epos + back[kg] + 1;
																		*gaplen = back[kg]+1;
																		*forwlen = 0;
																		*backlen = readlen-1-back[kg];
																	}
																	else if((readlen-1-back[kg]) > 0){
																		if(deltap == 0){
																			*gaplen = -1*realgap;
																			*forwlen = back[kg]+1;
																			*backlen = readlen - 1 -back[kg];
																			sprintf(cigar, "%dM%dD%dM", (back[kg]+1) , realgap, (readlen-1-back[kg]));
																		}
																		else{
																			*gaplen = realgap;
																			*forwlen = back[kg]+1-deltap;
																			*backlen = readlen - 1 -back[kg];
																			sprintf(cigar, "%dM%dI%dM", (back[kg]+1-deltap), realgap, (readlen-1-back[kg]));
																		}
																	}
																	else{
																		if(deltap == 0){
																			sprintf(cigar, "%dM", readlen);
																			*gaplen = 0;
																		}
																		else{
																			*gaplen = realgap;
																			*forwlen = readlen-realgap;
																			*backlen = 0;
																			sprintf(cigar, "%dM%dS", (readlen-realgap), realgap);
																		}
																	}
																}
															}
														}
														if((dclip != 0)&&(mscore > readlen*3/5*match)){
															if(dclip == 1){
																mscore = 0;
																j = 0;
																for(k = 0; (k <= iclip)&&(k <= readlen/30); k++){
																	if(((forw[k]-k)*match-k*mis*3) > mscore){
																		mscore = (forw[k]-k)*match-k*mis*3;
																		j = k;
																	}
																}
																if(forw[j] > readlen/2){
																	iclip = j;
																if((iclip + 1 + (readlen-forw[iclip])*10/readlen) < *nmin){
																	*nmin = iclip + 1 + (readlen-forw[iclip])*10/readlen;
																	*gaplen = readlen - forw[iclip];
																	*forwlen = forw[iclip];
																	*backlen = 0;
																	sprintf(cigar, "%dM%dS", *forwlen, *gaplen);
																	*bestpos = spos;
																	nmis = *nmin;
																}
																}
															}
															else{
																mscore = 0;
																j = 0;
																for(k = 0; (k <= iclip)&&(k <= readlen/30); k++){
																	if(((readlen-back[k]-k)*match-k*mis*3) > mscore){
																		mscore = (readlen-back[k]-k)*match-k*mis*3;
																		j = k;
																	}
																}
																if(back[j] < readlen/2){
																	iclip = j;
																if((iclip + 1 + (1+back[iclip])*10/readlen) < *nmin){
																	*nmin = iclip + 1 + (1+back[iclip])*10/readlen;
																	sprintf(cigar, "%dS%dM", (back[iclip]+1), (readlen-1-back[iclip]));
																		*bestpos = epos + back[iclip] + 1;
																		*gaplen = back[iclip]+1;
																		*forwlen = 0;
																		*backlen = readlen-1-back[iclip];
																		nmis = *nmin;
																}
																}
															}
														}
													}
												}
												
												if((nmis > *nmin)&&(readlen > 47)){
													if((effmis <= 2)&&(effl >= readlen*78/100)){
															if(oldeffl >= 30){
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
										
																		sprintf(cigar, "%dM%dS", effl, (readlen -effl));
																		*bestpos = ispos;
																		*nmin = nmis;
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
																			*bestpos = ispos+readlen-effl;
																			*gaplen = readlen-effl;
																			*forwlen = 0;
																			*backlen = effl;
																	}
																	
															}
													}
												}
												
									
		
}

/* search mapping results for all reference positions of a seed */ 
void Gapalign(int ipos, int readlen, uint64_t* fqseqint, uint32_t* fqseedint, uint64_t* seqint, uint16_t* seednum, uint32_t* seedind, uint32_t* seedpos, uint32_t* chrpos, uint32_t nchr, int* nmin, uint32_t* bestpos, int* inbest, uint32_t* poslist, int* nmappos, int* nmislist, int* direct, int direction, int* gaplen, int* forwlen, int* backlen, int* isgap, char* cigar){
		int i = 0;
		int j = 0;
		int k = 0;
		int m = 0;
		int p = 1;
		int nmis = 0;
		int found = 0;
		int forw[100];
		int back[100];
		uint32_t spos, epos, gapl, delta;
		int sp = 0;
		int ep = 0;
		int normis = 15;
		int bytemis = 100;
		int shift = 0;
		int begmis = 10;
		int effmis = 0;
		int termis = 0;
		int termisb = 0;
		int nmisth = 0;
		int gapth = 1;
		int effl = 0;
		int oldeffl = 0;
		int kgap = 0;
		int searchindel = 0;
		int shiftsearch = 0;
		int igshift = 0;
		int match = 2;
		int mis = 5;
		int gapo = 5;
		int gape = 1;
		int gapmis = 10;
		int nforw, nback, kmax, imisp, ntotmis, kg, jg, deltap, realgap, mscore, iclip, dclip;
		int forw2, back2, nforwmis2, nbackmis2, n2bytemis; 
		
		bool newres = false;
		bool ismis = false;
		uint32_t nshiftmis = 0;
		uint32_t nshiftmis2 = 0;
		uint32_t nshiftshit = 0;
		uint32_t nshiftgood = 0;
		uint64_t seqdot[20];
						for(i = 0; (i < seednum[fqseedint[ipos]])&&(*nmin >= 0); i++){
								if((seedpos[seedind[fqseedint[ipos]]+i] > (ipos+32))&&(seedpos[seedind[fqseedint[ipos]]+i] < maxreflength)){
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
												forw2 = 0;
												back2 = 0;
												n2bytemis = 0;
												nbackmis2 = readlen/10;
												nforwmis2 = readlen/10;
												ismis = false;
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
																		if(p < 2)
																			n2bytemis += bytemis;
																		if(p < 3)
																			begmis += bytemis;
																		if(p < 4)
																			termis += bytemis;
																}
																if((readlen-effl) < 33)
																	termisb += bytemis;
																
																if((nmis < readlen/25)&&(bytemis < 3)){
																	if((bytemis < 2)||(forw2 < 33)){
																		forw2 = effl;
																		nforwmis2 = nmis;
																	}
																}
																
														}
												}
												oldeffl = effl;
												shiftsearch = 0;
												if((forw2 < readlen*3/5)||(forw2 < 32)){
													forw2 = 0;
												}
												else if((nforwmis2 + 1 + (readlen-forw2)/10) >= nmis){
													forw2 = 0;
												}
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
														
														if((nmis > readlen/30)&&(termis > nmis*50/readlen)){
														if(((nmis - termis) < readlen/30)&&(forw2 == 0)){
															if(nbackmis2 >= (termis-1)){
																nbackmis2 = nmis - n2bytemis;
																back2 = readlen - 16;
															}
															else if(begmis == termis){
																back2 = readlen - 24;
																nbackmis2 = nmis - termis;
															}
															else{
																back2 = readlen - 32;
																nbackmis2 = nmis - termis;
															}
														}
														}
														if(((nbackmis2 + 1 + (readlen-back2)/10) >= nmis)||(back2 < readlen*3/5)){
															back2 = 0;
														}
														
														if((*nmappos == 0)||(nmis < *nmin)){
															if((*nmin > (nmis+1))&&(*nmin > 2)){
																*nmappos = 0;
															}
															if(*nmappos >= maxnmappos)
																*nmappos = maxnmappos-1;
															*inbest = *nmappos;
															*nmin = nmis;
															*bestpos = spos;
															poslist[*nmappos] = spos;
															nmislist[*nmappos] = nmis;
															direct[*nmappos] = direction;
															gaplen[*nmappos] = 0;
															(*nmappos) += 1;
															ismis = true;
															sprintf(cigar, "%dM", readlen);												
														}
														else if((poslist[*nmappos-1] != spos)&&(*bestpos != spos)){
															if((nmis < (*nmin + 2))&&(*nmappos < maxnmappos)){
																poslist[*nmappos] = spos;
																nmislist[*nmappos] = nmis;
																direct[*nmappos] = direction;
																gaplen[*nmappos] = 0;
																(*nmappos) += 1;
																ismis = true;
															}
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
																		if(gapl == 0)
																			found = 0;
																				
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
																				if((nmis < readlen/25)&&(bytemis < 3)){
																					if((bytemis < 2)||(back2 < 33)){
																						back2 = effl;
																						nbackmis2 = nmis;
																					}
																				}
																		}
																}
																if(((nbackmis2 + 1 + (readlen-back2)/10) > *nmin)||(back2 < readlen*3/5)){
																	back2 = 0;
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
																				if(gapl == 0)
																					found = 0;
																				
																				if((found == 0)&&((readlen-effl) < 29))
																					shiftsearch = -1;
																		}
																}
														}
												}

												if((shiftsearch != 0)&&(readlen > 47)){
														if((spos > (readlen/2+36))&&(spos < maxreflength)){
															uint64_t lastseqint = 0;
															uint32_t ankpos = spos+readlen-8;
															igshift = -30;
															if(shiftsearch == 1){
																	if(readlen % 32 == 0)
																			lastseqint = fqseqint[readlen/32-1];
																	else
																			lastseqint = (maskl[readlen%32]&(fqseqint[readlen/32]>>(2*(32-readlen%32))))
																							| (~maskl[readlen%32]&(fqseqint[readlen/32-1] << (2*(readlen%32))));
																	lastseqint &= maskl[8];
																	igshift = Shiftseq(lastseqint, ankpos, shiftsearch, readlen, effl, seqint);
																	if(igshift > -9)
																			epos = spos - igshift;
															}
															else{
																	lastseqint = fqseqint[0] >> 48;
																	lastseqint &= maskl[8];
																	igshift = Shiftseq(lastseqint, spos, shiftsearch, readlen, effl, seqint);
																	if(igshift > -9)
																			spos = epos -igshift;
															}
													
														}
												}
												if((spos != epos )&&(spos < maxreflength)){
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
													if((realgap < readlen/4)&&(realgap > 0)){
														nmis = 0;
														mscore = 0;
														iclip = -1;
														dclip = 0;
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
																				if(((forw[imisp]-imisp)*match-imisp*mis) > mscore){
																					mscore = (forw[imisp]-imisp)*match-imisp*mis;
																					iclip = imisp;
																					dclip = 1;
																				}
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
																				if(((readlen-back[imisp]-1-imisp)*match-imisp*mis) > mscore){
																					mscore = (readlen-back[imisp]-1-imisp)*match-imisp*mis;
																					iclip = imisp;
																					dclip = -1;
																				}
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
														
														mscore -= 2;
														if((jg > -1)&&(kg > -1)){
														
														if((dclip != 0)&&(((readlen-deltap-ntotmis)*match-ntotmis*mis-gapo-(realgap-1)*gape) > mscore)){
																iclip = -1;
																dclip = 0;
														}
														else if(ntotmis < 1){
																if(forw[jg] > (back[kg]+deltap+1)){
																		iclip = -1;
																		dclip = 0;
																}
														}
														
														if(dclip == 0){
															if((ntotmis + 1 + realgap/gapmis) <= *nmin){		
																newres = false;
																if((ntotmis + 1 + realgap/gapmis) == *nmin){
																	if(*nmappos == 0)
																		newres = true;
																	else if(!ismis){
																		if((abs(*bestpos-spos) > readlen/3) && (abs(poslist[*nmappos-1]-spos) > readlen/3)){
																			newres = true;
																		}
																	}
																}
																else
																	newres = true;
																	
																if(newres){
																	if(ismis)
																			(*nmappos) -=  1;
																	else if(*nmappos >= maxnmappos){
																			*nmappos = maxnmappos-1;
																	}
																	if((*nmappos == 0)||((ntotmis + 1 + realgap/gapmis) < *nmin)){
																			if(*nmin > 2)
																				*nmappos = 0;
																		*nmin = ntotmis + 1 + realgap/gapmis;
																		*bestpos = spos;
																		*inbest = *nmappos;
																	}
																	
																if(*nmappos < maxnmappos){
																	poslist[*nmappos] = spos;
																	nmislist[*nmappos] = *nmin;
																	direct[*nmappos] = direction;																	
																	nmis = *nmin;
																	if(back[kg] <= (deltap+1)){
																		if(*inbest == *nmappos){
																			sprintf(cigar, "%dS%dM", (back[kg]+1), (readlen-1-back[kg]));
																			*bestpos = epos + back[kg] + 1;
																		}
																		poslist[*nmappos] = *bestpos;
																		gaplen[*nmappos] = (back[kg]+1);
																		forwlen[*nmappos] = 0;
																		backlen[*nmappos] = readlen-(back[kg]+1);
																	}
																	else if(back[kg] < (readlen -1)){
																		if(deltap == 0){
																			gaplen[*nmappos] = -1*realgap;
																			forwlen[*nmappos] = back[kg]+1;
																			backlen[*nmappos] = readlen-1-back[kg];
																			if(*inbest == *nmappos)
																				sprintf(cigar, "%dM%dD%dM", (back[kg]+1) , realgap, (readlen-1-back[kg]));
																		}
																		else{
																			gaplen[*nmappos] = realgap;
																			forwlen[*nmappos] = back[kg]+1-deltap;
																			backlen[*nmappos] = readlen-1-back[kg];
																			if(*inbest == *nmappos)
																				sprintf(cigar, "%dM%dI%dM", (back[kg]+1-deltap), realgap, (readlen-1-back[kg]));
																		}
																	}
																	else{
																		if(deltap == 0){
																			if(*inbest == *nmappos)
																				sprintf(cigar, "%dM", readlen);
																			gaplen[*nmappos] = 0;
																		}
																		else{
																			if(*inbest == *nmappos)
																				sprintf(cigar, "%dM%dS", (readlen-realgap), realgap);
																			gaplen[*nmappos] = realgap;
																			forwlen[*nmappos] = readlen - realgap;
																			backlen[*nmappos] = 0;
																		}
																	}
																	(*nmappos) += 1;
																}
																*isgap = 1;
																}
															}
														}
														}
													
														if((dclip != 0)&&(mscore > readlen*3/5*match)){
															if(dclip == 1){
																mscore = 0;
																j = 0;
																for(k = 0; (k <= iclip)&&(k <= readlen/30); k++){
																	if(((forw[k]-k)*match-k*mis*3) > mscore){
																		mscore = (forw[k]-k)*match-k*mis*3;
																		j = k;
																	}
																}
																if(forw[j] > readlen/2){
																	iclip = j;
																if((iclip + 1 + (readlen-forw[iclip])*10/readlen) < *nmin){
																	
																	if(ismis)
																			(*nmappos) -=  1;
																	else{
																		if(*nmin > 2)
																			*nmappos = 0;
																		if(*nmappos >= maxnmappos)
																			*nmappos = maxnmappos-1;
																	}
																	*nmin = iclip + 1 + (readlen-forw[iclip])*10/readlen;
																	*bestpos = spos;
																	*inbest = *nmappos;
																	poslist[*nmappos] = spos;
																	nmislist[*nmappos] = *nmin;
																	direct[*nmappos] = direction;																	
																	nmis = *nmin;
																	
																	gaplen[*nmappos] = readlen - forw[iclip];
																	forwlen[*nmappos] = forw[iclip];
																	backlen[*nmappos] = 0;
																	sprintf(cigar, "%dM%dS", forwlen[*nmappos], gaplen[*nmappos]);
																	(*nmappos) += 1;
																}
																}
															}
															else{
																mscore = 0;
																j = 0;
																for(k = 0; (k <= iclip)&&(k <= readlen/30); k++){
																	if(((readlen-back[k]-k)*match-k*mis*3) > mscore){
																		mscore = (readlen-back[k]-k)*match-k*mis*3;
																		j = k;
																	}
																}
																if(back[j] < readlen/2){
																	iclip = j;
																if((iclip + 1 + (1+back[iclip])*10/readlen) < *nmin){
																	if(ismis)
																			(*nmappos) -=  1;
																	else{
																		if(*nmin > 2)
																			*nmappos = 0;
																		if(*nmappos >= maxnmappos)
																			*nmappos = maxnmappos-1;
																	}
																	*nmin = iclip + 1 + (1+back[iclip])*10/readlen;
																	*inbest = *nmappos;
																	nmislist[*nmappos] = *nmin;
																	direct[*nmappos] = direction;																	
																	nmis = *nmin;
																	*bestpos = epos + back[iclip] + 1;
																	poslist[*nmappos] = *bestpos;
																	gaplen[*nmappos] = (back[iclip]+1);
																	forwlen[*nmappos] = 0;
																	backlen[*nmappos] = readlen-(back[iclip]+1);
																	sprintf(cigar, "%dS%dM", (back[iclip]+1), (readlen-1-back[iclip]));
																	(*nmappos) += 1;
																}
																}
															}
															
														}
													}
													//}
												}
												
											if((*nmin > readlen/30)&&(*nmin > 2)){
												
												if(forw2 > readlen/2){
													nmis = nforwmis2 + 1 + (readlen-forw2)/10;
													if(nmis < *nmin){
																		sprintf(cigar, "%dM%dS", forw2, (readlen - forw2));
																		*bestpos = seedpos[seedind[fqseedint[ipos]]+i]-ipos;
																		if(ismis)
																			(*nmappos) -=  1;
																		else{
																			if(*nmin > 3)
																				*nmappos = 0;
																			if(*nmappos >= maxnmappos)
																				*nmappos = maxnmappos-1;
																		}
																		*inbest = *nmappos;
																		*nmin = nmis;
																		poslist[*nmappos] = *bestpos;
																		nmislist[*nmappos] = *nmin;
																		direct[*nmappos] = direction;
																		gaplen[*nmappos] = readlen -forw2;
																		forwlen[*nmappos] = forw2;
																		backlen[*nmappos] = 0;
																		(*nmappos) += 1;
																		*isgap = 1;
													}
												}
												else if (back2 > readlen/2){
													nmis = nbackmis2 + 1 + (readlen-back2)/10;
													if(nmis < *nmin){
																			sprintf(cigar, "%dS%dM", (readlen-back2), back2);
																			*bestpos = seedpos[seedind[fqseedint[ipos]]+i]-ipos+readlen-back2;
																		if(ismis)
																			(*nmappos) -=  1;
																		else{
																			if(*nmin > 3)
																				*nmappos = 0;
																			if(*nmappos >= maxnmappos)
																				*nmappos = maxnmappos-1;
																		}
																			*nmin = nmis;
																			*inbest = *nmappos;
																			poslist[*nmappos] = *bestpos;
																			nmislist[*nmappos] = *nmin;
																			direct[*nmappos] = direction;
																			gaplen[*nmappos] = readlen -back2;
																			forwlen[*nmappos] = 0;
																			backlen[*nmappos] = back2;
																			(*nmappos) += 1;
																			*isgap = 1;
													}
												}
												
						
													else if((effmis <= 2)&&(effl >= readlen*75/100)){
															
															if(oldeffl >= 30){
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
										
																		sprintf(cigar, "%dM%dS", effl, (readlen -effl));
																		*bestpos = seedpos[seedind[fqseedint[ipos]]+i]-ipos;
																		if(ismis)
																			(*nmappos) -=  1;
																		else{
																			if(*nmin > 3)
																				*nmappos = 0;
																			if(*nmappos >= maxnmappos)
																				*nmappos = maxnmappos-1;
																		}
																		*inbest = *nmappos;
																		*nmin = nmis;
																		poslist[*nmappos] = *bestpos;
																		nmislist[*nmappos] = *nmin;
																		direct[*nmappos] = direction;
																		gaplen[*nmappos] = readlen -effl;
																		forwlen[*nmappos] = effl;
																		backlen[*nmappos] = 0;
																		(*nmappos) += 1;
																		*isgap = 1;
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
																			
																			sprintf(cigar, "%dS%dM", (readlen-effl), effl);
																			*bestpos = seedpos[seedind[fqseedint[ipos]]+i]-ipos+readlen-effl;
																		if(ismis)
																			(*nmappos) -=  1;
																		else{
																			if(*nmin > 3)
																				*nmappos = 0;
																			if(*nmappos >= maxnmappos)
																				*nmappos = maxnmappos-1;
																		}
																			*nmin = nmis;
																			*inbest = *nmappos;
																			poslist[*nmappos] = *bestpos;
																			nmislist[*nmappos] = *nmin;
																			direct[*nmappos] = direction;
																			gaplen[*nmappos] = readlen -effl;
																			forwlen[*nmappos] = 0;
																			backlen[*nmappos] = effl;
																			(*nmappos) += 1;
																			*isgap = 1;
																	}
																	
															}
													}
											
											}
											
								}			
						}				
}

/* search for mapping results for pairs in buffer */
void* aln_pe(void *marg){
		args * data = (args*) marg;
				int maxseqlen = 500;
				int readlen = 500;
				char fqid[maxseqlen];
				char fqseq[maxseqlen];
				char revfqseq[maxseqlen];
				char fqual[maxseqlen];
				char revfqual[maxseqlen];
				int readlenb = 500;
				char fqidb[maxseqlen];
				char fqseqb[maxseqlen];
				char revfqseqb[maxseqlen];
				char fqualb[maxseqlen];
				char revfqualb[maxseqlen];
				uint64_t fqseqint [20];
				uint64_t revfqseqint [20];
				uint64_t seqdot [20];
				uint64_t revseqdot [20]; 
				uint32_t fqseedint [500];
				uint32_t revfqseedint [500];
				uint64_t fqseqintb [20];
				uint64_t revfqseqintb [20];
				uint64_t seqdotb [20];
				uint64_t revseqdotb [20]; 
				uint32_t fqseedintb [500];
				uint32_t revfqseedintb [500];
				uint64_t nmapped = 0;
				uint64_t nunimapped = 0;
				uint64_t nlmapped = 0;
				uint64_t nrmapped = 0;
				uint64_t nlunimapped = 0;
				uint64_t nrunimapped = 0;
				int iseed = 0;
				int ifq = 0;
				int nmin = readlen*10/100+1;
				int nminb = readlen*10/100+1;
				int oldnmin = nmin;
				int nbest = 0;
				int direction = 0;
				int directiona = 0;
				int directionb = 0;
				int score, scoreb;
				int minpos1 = 100000;
				int minpos2 = 100000;
				int minpos3 = 100000;
				int minpos4 = 100000;
				uint32_t bestpos = 0;
				int maxnpos = maxnmappos+1;
				uint32_t mappos[maxnpos];
				int mapmis[maxnpos];
				int direct[maxnpos];
				int nmappos, oldnmappos;
				char cigar[20];
				uint32_t bestposb = 0;
				uint32_t ispos = 0;
				
				uint32_t mapposb[maxnpos];
				int mapmisb[maxnpos];
				int directb[maxnpos];
				int gaplen[maxnpos];
				int forwlen[maxnpos];
				int backlen[maxnpos];
				int gaplenb[maxnpos];
				int forwlenb[maxnpos];
				int backlenb[maxnpos];
				int nmapposb, oldnmapposb;
				int imappos, imapposb, misab, ismapa, ismapb;
				char cigarb[20];
				char newcigar[20];
				char newcigarb[20];			
				int ipos[40], iposb[40];
				int nsfq, isbg, isgap, isgapb, niposa, niposb;
				int newnmin, newnminb;
				int inbest = -1;
				int inbestb = -1;
				uint32_t pairdis, newpos, newposb, delta, gapl;

				int i, j, k, p, l, m, sp, ep, kgap, found, numn;
				int pairmis = 2;
				int step = 6;
				bool searchpair = false;
				pairdis = 1000;
				
				for(p = (data->ip); p < (data->ip+data->bufsize); p++){
								iseed = 0;
								nmappos = 0;
								oldnmappos = 0;
								oldnmapposb = 0;
								directiona = 0;
								directionb = 0;
								isgap = 0;
								isgapb = 0;
								numn = 0;
								
								for(i = 0; i < 20; i++){
										ipos[i] = -1;
								}
								readlen = strlen((data->seqbuffer+p)->seq);
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
					if(numn < readlen/10){
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
						oldnmin = readlen*10/100+1;
						nmin = oldnmin;
						*cigar = '*';
						bestpos = 0;
						isgap = 0;
						for(i = 0; i < 4; i++){					
								if(ipos[i] >= 0){
										if(i%2 == 0){
												direction = 0;
												Gapalign(ipos[i], readlen, fqseqint, fqseedint, data->seqint, data->seednum, data->seedind, data->seedpos, data->chrpos, data->nchr, &nmin, &bestpos, &inbest, mappos, &nmappos, mapmis, direct, direction, gaplen, forwlen, backlen, &isgap, cigar);
						
										}
										else{
												direction = 1;
												Gapalign(ipos[i], readlen, revfqseqint, revfqseedint, data->seqint, data->seednum, data->seedind, data->seedpos, data->chrpos, data->nchr, &nmin, &bestpos, &inbest, mappos, &nmappos, mapmis, direct, direction, gaplen, forwlen, backlen, &isgap, cigar);
										}
										if(nmin < oldnmin){
												
												if(i%2 == 0)
														directiona = 0;
												else
														directiona = 1;
												
												oldnmin = nmin;
										}
										else if(nmin == oldnmin){
												if((nmappos > 0)&&(oldnmappos == 0)){
													if(i%2 == 0)
														directiona = 0;
													else
														directiona = 1;
												}
										}
										oldnmappos = nmappos;
								}
						}
						niposa = 0;
						if(nmin > 0){
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

								found = 0;
								for(i = 4; (i < (nsfq+1)*4)&&(found < 1); i++){			
									if(ipos[i] >= 0){
										niposa ++;
											if(i%2 == 0){
													direction = 0;
													Gapalign(ipos[i], readlen, fqseqint, fqseedint, data->seqint, data->seednum, data->seedind, data->seedpos, data->chrpos, data->nchr, &nmin, &bestpos, &inbest, mappos, &nmappos, mapmis, direct, direction, gaplen, forwlen, backlen, &isgap, cigar);
											}
											else{
													direction = 1;
													Gapalign(ipos[i], readlen, revfqseqint, revfqseedint, data->seqint, data->seednum, data->seedind, data->seedpos, data->chrpos, data->nchr, &nmin, &bestpos, &inbest, mappos, &nmappos, mapmis, direct, direction, gaplen, forwlen, backlen, &isgap, cigar);
											}
											if(nmin < oldnmin){
													
													if(i%2 == 0)
														directiona = 0;
													else
														directiona = 1;
													oldnmin = nmin;
											}
											else if(nmin == oldnmin){
												if((nmappos > 0)&&(oldnmappos == 0)){
													if(i%2 == 0)
														directiona = 0;
													else
														directiona = 1;
												}
											}
											oldnmappos = nmappos;
									
									
									if((nmin < 3)&&(i%4 == 3)){
											if(nmin == 1)
													found = 1;
											else if(i >= 11)
													found = 1;
									}
									}
								}
						}
					}
								iseed = 0;
								nmapposb = 0;
								numn = 0;
								niposb = 0;
								for(i = 0; i < 20; i++){
										iposb[i] = -1;
								}
								readlenb = strlen((data->seqbufferb+p)->seq);
								strcpy(fqseqb, (data->seqbufferb+p)->seq);
								strcpy(fqualb, (data->seqbufferb+p)->qual);
								for(i = 0; i < (readlenb/32+1); i++){
										fqseqintb[i] = 0;
										revfqseqintb[i] = 0;
								}
								for(i = 0; i < readlenb; i++){
										if((fqseqb[i] == 'N')||(fqseqb[i] == 'n')){
											numn ++;
										}
										ifq = i/32;																	
										fqseqintb[ifq] = (fqseqintb[ifq] << 2)|(c_to_i(fqseqb + i) & 3);
										if(i >= 14)
											iseed ++;
										if(iseed == 0)
											fqseedintb[iseed] = mask[7] & ((fqseedintb[iseed] << 2)|(c_to_i(fqseqb + i) & 3));
										else
											fqseedintb[iseed] = mask[7] & ((fqseedintb[iseed-1] << 2)|(c_to_i(fqseqb + i) & 3));
										revfqseqb[i] = rev(fqseqb+readlenb-1-i);
										revfqualb[i] = *(fqualb+readlenb-1-i);
								}
				if(numn < readlenb/10){
								if(readlenb%32 != 0)
										fqseqintb[readlenb/32] = fqseqintb[readlenb/32] << (2*(32-readlenb%32));								
								iseed = 0;
								for(i = 0; i < readlenb; i++){
										ifq = i/32;								 
										revfqseqintb[ifq] = (revfqseqintb[ifq] << 2)|(c_to_i(revfqseqb + i) & 3);
										if(i >= 14)
											iseed ++;
										if(iseed == 0)
											revfqseedintb[iseed] = mask[7] & ((revfqseedintb[iseed] << 2)|(c_to_i(revfqseqb + i) & 3));
										else
											revfqseedintb[iseed] = mask[7] & ((revfqseedintb[iseed-1] << 2)|(c_to_i(revfqseqb + i) & 3));
								}
								if(readlenb%32 != 0)
										revfqseqintb[readlenb/32] = revfqseqintb[readlenb/32] << (2*(32-readlenb%32));
								minpos1 = 8000000;
								minpos2 = 8000000;
								for(i = 0; i <= (readlenb-14); i++){
										if(data->seednum[fqseedintb[i]] > 0){
												if(i % 2 == 0){
													if(data->seednum[fqseedintb[i]] < minpos1){
														minpos1 = data->seednum[fqseedintb[i]];
														iposb[0] = i;
													}
												}
												else{
													if(data->seednum[fqseedintb[i]] < minpos2){
														minpos2 = data->seednum[fqseedintb[i]];
														iposb[2] = i;
													}
												}
										}
								}

						minpos1 = 8000000;
						minpos2 = 8000000;
						for(i = 0; i <= (readlenb-14); i++){
								if(data->seednum[revfqseedintb[i]] > 0){
										if(i % 2 == 0){
												if(data->seednum[revfqseedintb[i]] < minpos1){
														minpos1 = data->seednum[revfqseedintb[i]];
														iposb[1] = i;
												}
										}
										else{
												if(data->seednum[revfqseedintb[i]] < minpos2){
														minpos2 = data->seednum[revfqseedintb[i]];
														iposb[3] = i;
												}
										}
								}
						}
						oldnmin = readlenb*10/100+1;
						nminb = oldnmin;
						*cigarb = '*';
						bestposb = 0;
						isgapb = 0;
						for(i = 0; i < 4; i++){										
								if(iposb[i] >= 0){
										if(i%2 == 0){
												direction = 0;												
												Gapalign(iposb[i], readlenb,  fqseqintb, fqseedintb, data->seqint, data->seednum, data->seedind, data->seedpos, data->chrpos, data->nchr, &nminb, &bestposb, &inbestb, mapposb, &nmapposb, mapmisb, directb, direction, gaplenb, forwlenb, backlenb, &isgapb,cigarb);
										}
										else{
												direction = 1;											
												Gapalign(iposb[i], readlenb, revfqseqintb, revfqseedintb, data->seqint, data->seednum, data->seedind, data->seedpos, data->chrpos, data->nchr, &nminb, &bestposb, &inbestb, mapposb, &nmapposb, mapmisb, directb, direction, gaplenb, forwlenb, backlenb, &isgapb, cigarb);
										}
										if(nminb < oldnmin){
												
												if(i%2 == 0)
														directionb = 0;
												else
														directionb = 1;
												oldnmin = nminb;
										}
										else if(nminb == oldnmin){
												if((nmapposb > 0)&&(oldnmapposb == 0)){
													if(i%2 == 0)
														directionb = 0;
													else
														directionb = 1;
												}
										}
										oldnmapposb = nmapposb;
								}
						}

						if(nminb > 0){
								nsfq = (readlenb > 75) ? 4 : 3;
								if(readlenb < 45)
										nsfq = 2;
								isbg = (readlenb > nsfq*14) ? (readlenb/nsfq-14) : 0;
				
								for(i = 0; i < nsfq; i++){
										minpos1 = 1000;
										minpos2 = 1000;
										minpos3 = 1000;
										minpos4 = 1000;
										for(j = i*readlenb/nsfq; j <= ((i+1)*readlenb/nsfq -14); j++){
												if((data->seednum[fqseedintb[j]] > 0)&&(data->seednum[fqseedintb[j]] < 500)){
														if(j%2 == 0){
																if(data->seednum[fqseedintb[j]] < minpos1){
																		minpos1 = data->seednum[fqseedintb[j]];
																		iposb[4+4*i] = j;
																}
														}
														else{
																if(data->seednum[fqseedintb[j]] < minpos2){
																		minpos2 = data->seednum[fqseedintb[j]];
																		iposb[6+4*i] = j;
																}
														}
												}
												if((data->seednum[revfqseedintb[j]] > 0)&&(data->seednum[revfqseedintb[j]] < 500)){
														if(j%2 == 0){
																if(data->seednum[revfqseedintb[j]] < minpos3){
																		minpos3 = data->seednum[revfqseedintb[j]];
																		iposb[5+4*i] = j;
																}
														}
														else{
																if(data->seednum[revfqseedintb[j]] < minpos4){
																		minpos4 = data->seednum[revfqseedintb[j]];
																		iposb[7+4*i] = j;
																}
														}
												}
										}
										if((minpos1 > 500)||(minpos2 > 500)){
												for(j = (isbg+i*readlenb/nsfq); j <= ((i+1)*readlenb/nsfq-8); j++){
														if((data->seednum[fqseedintb[j]] > 0)&&(data->seednum[fqseedintb[j]] < 500)){
																if((j%2 == 0)&&(minpos1 > 500)){
																		minpos1 = data->seednum[fqseedintb[j]];
																		iposb[4+4*i] = j;
																}
																if((j%2 == 1)&&(minpos2 > 500)){
																		minpos2 = data->seednum[fqseedintb[j]];
																		iposb[6+4*i] = j;
																}
														}
												}
										}
										if((minpos3 > 500)||(minpos4 > 500)){
												for(j = (isbg+i*readlenb/nsfq); j <= ((i+1)*readlenb/nsfq-8); j++){
														if((data->seednum[revfqseedintb[j]] > 0)&&(data->seednum[revfqseedintb[j]] < 500)){
																if((j%2 == 0)&&(minpos3 > 500)){
																		minpos3 = data->seednum[revfqseedintb[j]];
																		iposb[5+4*i] = j;
																}
																if((j%2 == 1)&&(minpos4 > 500)){
																		minpos4 = data->seednum[revfqseedintb[j]];
																		iposb[7+4*i] = j;
																}
														}
												}
										}

										if(iposb[4+4*i] == iposb[0]){
											iposb[4+4*i] = -1;
										}
										if(iposb[5+4*i] == iposb[1]){
											iposb[5+4*i] = -1;
										}
										if(iposb[6+4*i] == iposb[2]){
											iposb[6+4*i] = -1;
										}
										if(iposb[7+4*i] == iposb[3]){
											iposb[7+4*i] = -1;
										}
								}
								found = 0;
								for(i = 4; (i < (nsfq+1)*4)&&(found < 1); i++){									
									if(iposb[i] >= 0){
										niposb ++;
											if(i%2 == 0){
													direction = 0;
													Gapalign(iposb[i], readlenb, fqseqintb, fqseedintb, data->seqint, data->seednum, data->seedind, data->seedpos, data->chrpos, data->nchr, &nminb, &bestposb, &inbestb, mapposb, &nmapposb, mapmisb, directb, direction, gaplenb, forwlenb, backlenb,  &isgapb, cigarb);
											}
											else{
													direction = 1;
													Gapalign(iposb[i], readlenb, revfqseqintb, revfqseedintb, data->seqint, data->seednum, data->seedind, data->seedpos, data->chrpos, data->nchr, &nminb, &bestposb, &inbestb, mapposb, &nmapposb, mapmisb, directb, direction, gaplenb, forwlenb, backlenb, &isgapb, cigarb);
											}
											if(nminb < oldnmin){
													
													if(i%2 == 0)
														directionb = 0;
													else
														directionb = 1;
													
													oldnmin = nminb;
											}
											else if(nminb == oldnmin){
												if((nmapposb > 0)&&(oldnmapposb == 0)){
													if(i%2 == 0)
														directionb = 0;
													else
														directionb = 1;
												}
											}
											oldnmapposb = nmapposb;
									
									if((nminb < 3)&&(i%4 == 3)){
											if(nminb == 1)
													found = 1;
											else if(i >= 11)
													found = 1;
									}
									}
								}
						}
				}						
						found = 0;
						nbest = 0;
						misab = readlen/10+readlenb/10;
						imappos = -1;
						imapposb = -1;
						int gapfound = 0;
						for(l = 0; l < nmappos; l++){	
							for(m = 0; m < nmapposb; m++){
								delta = (mappos[l] > mapposb[m]) ? (mappos[l] - mapposb[m]): (mapposb[m] - mappos[l]);
								if((delta < pairdis)&&((mapmis[l]+mapmisb[m]) <= misab)){
									if((nbest == 0)||((mapmis[l]+mapmisb[m]) < misab)){
									nbest = 1;
									
									if(((mapmis[l]+mapmisb[m]-misab+1) == 0)&&(misab > 5)){
										if((mappos[l] != mappos[imappos])&&(mapposb[m] != mapposb[imapposb]))
											nbest = 2;
									}
									
									imappos = l;
									imapposb = m;
									misab = mapmis[l]+mapmisb[m];
									}
									else if((mappos[l] != mappos[imappos])&&(mapposb[m] != mapposb[imapposb]))
										nbest ++;
									found = 1;
									
								}
								
							}
						}
						
						
						if(found == 1){
							if((nmappos == maxnmappos)&&(nmapposb == maxnmappos))
								nbest = 2;
							directiona = direct[imappos];
							directionb = directb[imapposb];
							inbest = imappos;
							inbestb = imapposb;
							if(nmappos > 1){
								if(gaplen[imappos] == 0){
									sprintf(cigar, "%dM", readlen);
								}
								else if(backlen[imappos] == 0){
									sprintf(cigar, "%dM%dS", forwlen[imappos], gaplen[imappos]);
								}
								else if(forwlen[imappos] == 0){
									sprintf(cigar, "%dS%dM", gaplen[imappos], backlen[imappos]);
								}
								else if(gaplen[imappos] < 0){
									sprintf(cigar, "%dM%dD%dM", forwlen[imappos], abs(gaplen[imappos]), backlen[imappos]);
								}
								else if(gaplen[imappos] > 0){
									sprintf(cigar, "%dM%dI%dM", forwlen[imappos], gaplen[imappos], backlen[imappos]);
								}			
							}
							if(nmapposb > 1){
								if(gaplenb[imapposb] == 0){
									sprintf(cigarb, "%dM", readlenb);
								}
								else if(backlenb[imapposb] == 0){
									sprintf(cigarb, "%dM%dS", forwlenb[imapposb], gaplenb[imapposb]);
								}
								else if(forwlenb[imapposb] == 0){
									sprintf(cigarb, "%dS%dM", gaplenb[imapposb], backlenb[imapposb]);
								}
								else if(gaplenb[imapposb] < 0){
									sprintf(cigarb, "%dM%dD%dM", forwlenb[imapposb], abs(gaplenb[imapposb]), backlenb[imapposb]);
								}
								else if(gaplenb[imapposb] > 0){
									sprintf(cigarb, "%dM%dI%dM", forwlenb[imapposb], gaplenb[imapposb], backlenb[imapposb]);
								}				
							}
						}
						searchpair = false;
						if(found == 0)
							searchpair = true;
						else{
							if((nbest == 1)&&((nmappos-maxnmappos)*(nmapposb-maxnmappos) == 0))
								searchpair = true;
						}
						if(searchpair){
							oldnmin = readlenb/10+1;
							newpos = 0;
							newposb = 0;
								if(readlenb > 100)
									step = 6*(readlenb/50);
								for(l = 0; (l < nmappos)&&(nmappos < maxnmappos); l++){
									if(l != imappos){
										newnminb = oldnmin;
									for(i = 0; (i < (readlenb-14))&&(newnminb >= oldnmin); i++){
										for(j = 0; (j < 2)&&(newnminb >= oldnmin); j++){
												gapfound = 0;
												sp = 0;
												if(j%2 == 0)
													ep = data->seednum[fqseedintb[i]]-1;
												else
													ep = data->seednum[revfqseedintb[i]]-1;
												gapl = 2500;
												while(ep >= sp){
													k = (sp+ep)/2;
													if(j%2 == 0)
														ispos = data->seedpos[data->seedind[fqseedintb[i]]+k];
													else
														ispos = data->seedpos[data->seedind[revfqseedintb[i]]+k];
													if(ispos > i){
														ispos -= i;
														delta = (ispos > mappos[l]) ? (ispos - mappos[l]) : (mappos[l]-ispos);
														if(delta < gapl){
															gapl = delta;
															kgap = k;
														}
														
													}
													else{
														ispos = 0;
													}
													if(ispos > mappos[l])
														ep = k-1;
													else if (ispos < mappos[l])
														sp = k+1;
													else
														ep = sp -1;
												}
													if(gapl < pairdis){
														if(j%2 == 0)
															ispos = data->seedpos[data->seedind[fqseedintb[i]]+kgap]-i;
														else
															ispos = data->seedpos[data->seedind[revfqseedintb[i]]+kgap]-i;
														
														if((ispos > (36+readlenb/2))&&(ispos < maxreflength)){
															gapfound = 1;
														}
													}										
											if(gapfound == 1){
												if(j%2 == 0)
													Spalign(ispos, readlenb, fqseqintb, fqseedintb, data->seqint, data->seednum, data->seedind, data->seedpos, data->chrpos, data->nchr, &newnminb, &newposb, (gaplenb+maxnpos-1), (forwlenb+maxnpos-1), (backlenb+maxnpos-1), newcigarb);
												else
													Spalign(ispos, readlenb, revfqseqintb, revfqseedintb, data->seqint, data->seednum, data->seedind, data->seedpos, data->chrpos, data->nchr, &newnminb, &newposb, (gaplenb+maxnpos-1), (forwlenb+maxnpos-1), (backlenb+maxnpos-1), newcigarb);
												
												if((newnminb < oldnmin)&&((newnminb+mapmis[l]) <= misab)){
													if(newnminb < nminb){
														nminb = newnminb;
														bestposb = newposb;
													}
													if(nbest > 0){
														if((newnminb+mapmis[l]) == misab){
															if((newposb != mapposb[imapposb])&&(mappos[l] != mappos[imappos]))
																nbest ++;
														}
													}
													if((nbest == 0)||((newnminb+mapmis[l]) < misab)){
														
														if(((newnminb+mapmis[l]+1-misab) == 0)&&(misab > 5))
															nbest = 2;
														else
															nbest = 1;
														
														misab = newnminb+mapmis[l];
														strcpy(cigarb, newcigarb);
														found = 1;
														imappos = l;
														imapposb = maxnmappos-1;
														mapposb[imapposb] = newposb;
														mapmisb[imapposb] = newnminb;
														directiona = direct[l];
														if(j%2 == 0)
															directionb = 0;
														else
															directionb = 1;
															if(gaplen[imappos] == 0)
																sprintf(cigar, "%dM", readlen);
															else if(backlen[imappos] == 0)
																sprintf(cigar, "%dM%dS", forwlen[imappos], gaplen[imappos]);
															else if(forwlen[imappos] == 0)
																sprintf(cigar, "%dS%dM", gaplen[imappos], backlen[imappos]);
															else if(gaplen[imappos] < 0)
																sprintf(cigar, "%dM%dD%dM", forwlen[imappos], abs(gaplen[imappos]), backlen[imappos]);
															else if(gaplen[imappos] > 0)
																sprintf(cigar, "%dM%dI%dM", forwlen[imappos], gaplen[imappos], backlen[imappos]);
													}
												}
												
											}
										}
										if(i%2 == 1)
											i += step;
									}
									}
								}
					
								if(readlen > 100)
									step = 6*(readlen/50);
								oldnmin = readlen/10+1;
								for(l = 0; (l < nmapposb)&&(nmapposb < maxnmappos); l++){
									
									if(l != imapposb){
										newnmin = oldnmin;
									for(i = 0; (i < (readlen-14))&&(newnmin >= oldnmin); i++){
										for(j = 0; (j < 2)&&(newnmin >= oldnmin); j++){
												gapfound = 0;
												sp = 0;
												if(j%2 == 0)
													ep = data->seednum[fqseedint[i]]-1;
												else
													ep = data->seednum[revfqseedint[i]]-1;
												gapl = 2500;
												while(ep >= sp){
													k = (sp+ep)/2;
													if(j%2 == 0)
														ispos = data->seedpos[data->seedind[fqseedint[i]]+k];
													else
														ispos = data->seedpos[data->seedind[revfqseedint[i]]+k];
													if(ispos > i){
														ispos -= i;
														delta = (ispos > mapposb[l]) ? (ispos - mapposb[l]) : (mapposb[l]-ispos);
														if(delta < gapl){
															gapl = delta;
															kgap = k;
														}
													}
													else{
														ispos = 0;
													}
													if(ispos > mapposb[l])
														ep = k-1;
													else if (ispos < mapposb[l])
														sp = k+1;
													else
														ep = sp -1;
												}
													if(gapl < pairdis){
														if(j%2 == 0)
															ispos = data->seedpos[data->seedind[fqseedint[i]]+kgap]-i;
														else
															ispos = data->seedpos[data->seedind[revfqseedint[i]]+kgap]-i;
														
														if((ispos > (36+readlen/2))&&(ispos < maxreflength)){
															gapfound = 1;
														}
													}
											
											if(gapfound == 1){
												if(j%2 == 0)
													Spalign(ispos, readlen, fqseqint, fqseedint, data->seqint, data->seednum, data->seedind, data->seedpos, data->chrpos, data->nchr, &newnmin, &newpos, (gaplen+maxnpos-1), (forwlen+maxnpos-1), (backlen+maxnpos-1), newcigar);
												else
													Spalign(ispos, readlen, revfqseqint, revfqseedint, data->seqint, data->seednum, data->seedind, data->seedpos, data->chrpos, data->nchr, &newnmin, &newpos, (gaplen+maxnpos-1), (forwlen+maxnpos-1), (backlen+maxnpos-1), newcigar);
											
												if((newnmin < oldnmin)&&((newnmin+mapmisb[l]) <= misab)){
													if(newnmin < nmin){
														nmin = newnmin;
														bestpos = newpos;
													}
													if(nbest > 0){
														if((newnmin+mapmisb[l]) == misab){
															if((newpos != mappos[imappos])&&(mapposb[l] != mapposb[imapposb]))
																nbest ++;
														}
													}
													if((nbest == 0)||((newnmin+mapmisb[l]) < misab)){
														
														if(((misab-newnmin-mapmisb[l]-1) == 0)&&(misab > 5))
															nbest = 2;
														else
															nbest = 1;
														
														misab = newnmin+mapmisb[l];
														//nbest = 1;
														strcpy(cigar, newcigar);
														found = 1;
														imappos = maxnmappos - 1;
														imapposb = l;
														mappos[imappos] = newpos;
														mapmis[imappos] = newnmin;
														directionb = directb[l];
														if(j%2 == 0)
															directiona = 0;
														else
															directiona = 1;	
															if(gaplenb[imapposb] == 0){
																sprintf(cigarb, "%dM", readlenb);
															}
															else if(backlenb[imapposb] == 0){
																sprintf(cigarb, "%dM%dS", forwlenb[imapposb], gaplenb[imapposb]);
															}
															else if(forwlenb[imapposb] == 0){
																sprintf(cigarb, "%dS%dM", gaplenb[imapposb], backlenb[imapposb]);
															}
															else if(gaplenb[imapposb] < 0){
																sprintf(cigarb, "%dM%dD%dM", forwlenb[imapposb], abs(gaplenb[imapposb]), backlenb[imapposb]);
															}
															else if(gaplenb[imapposb] > 0){
																sprintf(cigarb, "%dM%dI%dM", forwlenb[imapposb], gaplenb[imapposb], backlenb[imapposb]);
															}	
													}
												}
											}
										}
										if(i%2 == 1)
											i += step;
									}
									}
								}
					
							if(found == 1){
								if(nmappos == 0)
									nmappos = 1;
								if(nmapposb == 0)
									nmapposb = 1;
								inbest = imappos;
								inbestb = imapposb;
							}
						}
						
						ismapa = 1;
						ismapb = 1;
						score = 0;
						scoreb = 0;
						
						bool cigarerr = false;
						if(nmappos > 0){
						if(*cigar == '*'){
							cigarerr = true;
							nmappos = 0;
							score = 0;
							directiona = 0;
						}
						}
						if(nmapposb > 0){
						if(*cigarb == '*'){
							cigarerr = true;
							nmapposb = 0;
							scoreb = 0;
							directionb = 0;
						}
						}
						if(cigarerr){
							found = 0;
							nbest = 0;
						}
						
									if(nbest == 0){
										if(nmappos == 1){
											if(nmin < 3)
												score = readlen;
											else
												score = readlen*readlen/(nmin*nmin*16);
										}
										else if(nmappos > 1)
											score = 1;
									}
									else if(nbest > 1)
										score = 1;
									else{
										if(mapmis[imappos] < 3)
											score = readlen;
										else
											score = readlen*readlen/(mapmis[imappos]*mapmis[imappos]*16);
									}
									if(nbest == 0){
										if(nmapposb == 1){
											if(nminb < 3)
												scoreb = readlenb;
											else
												scoreb = readlenb*readlenb/(nminb*nminb*16);
										}
										else if(nmapposb > 1)
											scoreb = 1;
									}
									else if(nbest > 1)
										scoreb = 1;
									else{
										if(mapmisb[imapposb] < 3)
											scoreb = readlenb;
										else
											scoreb = readlenb*readlenb/(mapmisb[imapposb]*mapmisb[imapposb]*16);
									}
									
									if(score > 60)
										score = 60;
									if(scoreb > 60)
										scoreb = 60;
									
									if(nbest == 1){
										
										if(score < scoreb)
											score = (score*4+scoreb)/5;
										else
											scoreb = (score+scoreb*4)/5;
									}
									
									
												int jca = -1;
												int jcb = -1;
												uint32_t posa = 0;
												uint32_t posb = 0;
												bool isunmap = false;
						if(nmappos > 0){
									ismapa = 0;
									nlmapped++;
									if(directiona == 1){
										strcpy((data->seqbuffer+p)->seq, revfqseq);
										strcpy((data->seqbuffer+p)->qual, revfqual);
									}
									
										posa = mappos[inbest];
									isunmap = false;
									for(int jc = 0; jc < data->nchr; jc++){
											if(posa < data->chrpos[jc]){
													jca = jc;
												if(jc > 0)
													(data->resbuffer+p)->pos = posa -  data->chrpos[jc-1];
												else
													(data->resbuffer+p)->pos = posa;
												if((posa + 5*readlen/4+1) > data->chrpos[jc]){
													if(gaplen[inbest] == 0){
														if((posa + readlen*3/4-1) > data->chrpos[jc]){
															if((posa + readlen/4) > data->chrpos[jc]){
																sprintf(cigar, "%dS%dM", (data->chrpos[jc]+1-posa), (posa+ readlen -1-data->chrpos[jc]));
																jca = jc +1;
																(data->resbuffer+p)->pos = 1;
															}
															else
																isunmap = true;
														}
														else if((posa + readlen -1) > data->chrpos[jc]){
															sprintf(cigar, "%dM%dS", (data->chrpos[jc]+1-posa), (posa + readlen -1-data->chrpos[jc]));
														}
													}
													else if(forwlen[inbest] == 0){
														if((backlen[inbest]+posa-1) > data->chrpos[jc]){
															if((backlen[inbest]+posa-1-readlen/4) < data->chrpos[jc]){
																sprintf(cigar, "%dS%dM%dS", (readlen-backlen[inbest]), (data->chrpos[jc]+1-posa), (backlen[inbest]+posa-1-data->chrpos[jc]));
															}
															else if((posa+readlen/4) > data->chrpos[jc]){
																sprintf(cigar, "%dS%dM", (data->chrpos[jc]+1-posa+gaplen[inbest]), (readlen+posa-1-gaplen[inbest]-data->chrpos[jc]));
																jca = jc +1;
																(data->resbuffer+p)->pos = 1;
															}
															else
																isunmap = true;
														}
													}
													else if(backlen[inbest] == 0){
														if((forwlen[inbest]+posa-1) > data->chrpos[jc]){
															if((forwlen[inbest]+posa-1-readlen/4) < data->chrpos[jc]){
																sprintf(cigar, "%dM%dS", (data->chrpos[jc]+1-posa), (readlen-1+posa-data->chrpos[jc]));
															}
															else
																isunmap = true;
														}
													}
													else if(gaplen[inbest] < 0){
														if((posa+readlen-1-gaplen[inbest]) > data->chrpos[jc]){
															if((posa+readlen-1-gaplen[inbest]-readlen/4) < data->chrpos[jc]){
																if(backlen[inbest] > (posa+readlen+3-gaplen[inbest]-data->chrpos[jc]))
																	sprintf(cigar, "%dM%dD%dM%dS", forwlen[inbest], -1*gaplen[inbest], (backlen[inbest]-posa-readlen+1+gaplen[inbest]+data->chrpos[jc]), (posa+readlen-1-gaplen[inbest]-data->chrpos[jc]));
																else if((posa+forwlen[inbest]-1) <= data->chrpos[jc])
																	sprintf(cigar, "%dM%dS", forwlen[inbest], (readlen-forwlen[inbest]));
																else if((posa+3*readlen/4) < data->chrpos[jc])
																	sprintf(cigar, "%dM%dS", (data->chrpos[jc]+1-posa), (readlen-1+posa-data->chrpos[jc]));
																else
																	isunmap = true;
															}
															else
																isunmap = true;
														}
													}
													else{
														if((posa+readlen-1-gaplen[inbest]) > data->chrpos[jc]){
															if((posa+readlen-1-gaplen[inbest]-readlen/4) < data->chrpos[jc]){
																if(backlen[inbest] > (posa+readlen+3-gaplen[inbest]-data->chrpos[jc]))
																	sprintf(cigar, "%dM%dI%dM%dS", forwlen[inbest], gaplen[inbest], (backlen[inbest]-posa-readlen+1+gaplen[inbest]+data->chrpos[jc]), (posa+readlen-1-gaplen[inbest]-data->chrpos[jc]));
																else if((posa+forwlen[inbest]-1) <= data->chrpos[jc])
																	sprintf(cigar, "%dM%dS", forwlen[inbest], (readlen-forwlen[inbest]));
																else if((posa+3*readlen/4) < data->chrpos[jc])
																	sprintf(cigar, "%dM%dS", (data->chrpos[jc]+1-posa), (readlen-1+posa-data->chrpos[jc]));
																else
																	isunmap = true;
															}
															else
																isunmap = true;
														}
													}
													
													if(isunmap){
															nmappos = 0;
															ismapa = 1;
															score = 0;
															directiona = 0;
															jca = -1;
															nlmapped --;
													}
													
												}
												break;
											}
											else if(posa == data->chrpos[jc]){
												(data->resbuffer+p)->pos = 1;
												jca = jc + 1;
												if(gaplen[inbest] == 0){
													sprintf(cigar, "1S%dM", (readlen-1));
												}
												else if(forwlen[inbest] == 0){
													sprintf(cigar, "%dS%dM", (readlen+1-backlen[inbest]), (backlen[inbest]-1));
												}
												else if(backlen[inbest] == 0){
													sprintf(cigar, "1S%dM%dS", (forwlen[inbest]-1), (readlen-forwlen[inbest]));
												}
												else if(gaplen[inbest] < 0){
													if(forwlen[inbest] > 3)
														sprintf(cigar, "1S%dM%dD%dM", (forwlen[inbest]-1), -1*gaplen[inbest], backlen[inbest]);
													else{
														sprintf(cigar, "%dS%dM", (readlen-backlen[inbest]), backlen[inbest]);
														(data->resbuffer+p)->pos = forwlen[inbest]-gaplen[inbest];
													}
												}
												else{
													if(forwlen[inbest] > 3)
														sprintf(cigar, "1S%dM%dI%dM", (forwlen[inbest]-1), gaplen[inbest], backlen[inbest]);
													else{
														sprintf(cigar, "%dS%dM", (readlen-backlen[inbest]), backlen[inbest]);
														(data->resbuffer+p)->pos = forwlen[inbest];
													}
												}	
												break;
											}
									}
									if(gaplen[inbest] == 0){
										(data->resbuffer+p)->nmis = mapmis[inbest];
									}
									else if((forwlen[inbest] != 0)&&(backlen[inbest] != 0)){
										(data->resbuffer+p)->nmis = mapmis[inbest] - 1 - abs(gaplen[inbest])/10 + abs(gaplen[inbest]);
									}
									else{
										(data->resbuffer+p)->nmis = mapmis[inbest] - 1 - abs(gaplen[inbest])/10;
									}
									if((data->resbuffer+p)->nmis < 0)
										(data->resbuffer+p)->nmis = mapmis[inbest];
									(data->resbuffer+p)-> mapscore = score;
								
									strcpy((data->resbuffer+p)->cigar, cigar);
						}
												
						if(nmapposb > 0){
									ismapb = 0;
									nrmapped++;
									if(directionb == 1){
										strcpy((data->seqbufferb+p)->seq, revfqseqb);
										strcpy((data->seqbufferb+p)->qual, revfqualb);
									}
									
										posb = mapposb[inbestb];
									isunmap = false;
									for(int jc = 0; jc < data->nchr; jc++){
											if(posb < data->chrpos[jc]){
													jcb = jc;
												if(jc > 0)
													(data->resbufferb+p)->pos = posb -  data->chrpos[jc-1];
												else
													(data->resbufferb+p)->pos = posb;
												
												if((posb + 5*readlenb/4+1) > data->chrpos[jc]){
													if(gaplenb[inbestb] == 0){
														if((posb + readlenb*3/4-1) > data->chrpos[jc]){
															if((posb + readlenb/4) > data->chrpos[jc]){
																sprintf(cigarb, "%dS%dM", (data->chrpos[jc]+1-posb), (posb+ readlenb -1-data->chrpos[jc]));
																jcb = jc +1;
																(data->resbufferb+p)->pos = 1;
															}
															else
															isunmap = true;
														}
														else if((posb + readlenb -1) > data->chrpos[jc]){
															sprintf(cigarb, "%dM%dS", (data->chrpos[jc]+1-posb), (posb + readlenb -1-data->chrpos[jc]));
														}
													}
													else if(forwlenb[inbestb] == 0){
														if((backlenb[inbestb]+posb-1) > data->chrpos[jc]){
															if((backlenb[inbestb]+posb-1-readlenb/4) < data->chrpos[jc]){
																sprintf(cigarb, "%dS%dM%dS", (readlenb-backlenb[inbestb]), (data->chrpos[jc]+1-posb), (backlenb[inbestb]+posb-1-data->chrpos[jc]));
															}
															else if((posb+readlenb/4) > data->chrpos[jc]){
																sprintf(cigarb, "%dS%dM", (data->chrpos[jc]+1-posb+gaplenb[inbestb]), (readlenb+posb-1-gaplenb[inbestb]-data->chrpos[jc]));
																jcb = jc +1;
																(data->resbufferb+p)->pos = 1;
															}
															else
																isunmap = true;
														}
													}
													else if(backlenb[inbestb] == 0){
														if((forwlenb[inbestb]+posb-1) > data->chrpos[jc]){
															if((forwlenb[inbestb]+posb-1-readlenb/4) < data->chrpos[jc]){
																sprintf(cigarb, "%dM%dS", (data->chrpos[jc]+1-posb), (readlenb-1+posb-data->chrpos[jc]));
															}
															else
																isunmap = true;
														}
													}
													else if(gaplenb[inbestb] < 0){
														if((posb+readlenb-1-gaplenb[inbestb]) > data->chrpos[jc]){
															if((posb+readlenb-1-gaplenb[inbestb]-readlenb/4) < data->chrpos[jc]){
																if(backlenb[inbestb] > (posb+readlenb+3-gaplenb[inbestb]-data->chrpos[jc]))
																	sprintf(cigarb, "%dM%dD%dM%dS", forwlenb[inbestb], -1*gaplenb[inbestb], (backlenb[inbestb]-posb-readlenb+1+gaplenb[inbestb]+data->chrpos[jc]), (posb+readlenb-1-gaplenb[inbestb]-data->chrpos[jc]));
																else if((posb+forwlenb[inbestb]-1) <= data->chrpos[jc])
																	sprintf(cigarb, "%dM%dS", forwlenb[inbestb], (readlenb-forwlenb[inbestb]));
																else if((posb+3*readlenb/4) < data->chrpos[jc])
																	sprintf(cigarb, "%dM%dS", (data->chrpos[jc]+1-posb), (readlenb-1+posb-data->chrpos[jc]));
																else
																	isunmap = true;
															}
															else
																isunmap = true;
														}
													}
													else{
														if((posb+readlenb-1-gaplenb[inbestb]) > data->chrpos[jc]){
															if((posb+readlenb-1-gaplenb[inbestb]-readlenb/4) < data->chrpos[jc]){
																if(backlenb[inbestb] > (posb+readlenb+3-gaplenb[inbestb]-data->chrpos[jc]))
																	sprintf(cigarb, "%dM%dI%dM%dS", forwlenb[inbestb], gaplenb[inbestb], (backlenb[inbestb]-posb-readlenb+1+gaplenb[inbestb]+data->chrpos[jc]), (posb+readlenb-1-gaplenb[inbestb]-data->chrpos[jc]));
																else if((posb+forwlenb[inbestb]-1) <= data->chrpos[jc])
																	sprintf(cigarb, "%dM%dS", forwlenb[inbestb], (readlenb-forwlenb[inbestb]));
																else if((posb+3*readlenb/4) < data->chrpos[jc])
																	sprintf(cigarb, "%dM%dS", (data->chrpos[jc]+1-posb), (readlenb-1+posb-data->chrpos[jc]));
																else
																	isunmap = true;
															}
															else
																isunmap = true;
														}
													}
													
													if(isunmap){
															nmapposb = 0;
															ismapb = 1;
															scoreb = 0;
															directionb = 0;
															jcb = -1;
															nrmapped --;
													}
													
												}
												break;
											}
											else if(posb == data->chrpos[jc]){
												(data->resbufferb+p)->pos = 1;
												jcb = jc + 1;
												if(gaplenb[inbestb] == 0){
													sprintf(cigarb, "1S%dM", (readlenb-1));
												}
												else if(forwlenb[inbestb] == 0){
													sprintf(cigarb, "%dS%dM", (readlenb+1-backlenb[inbestb]), (backlenb[inbestb]-1));
												}
												else if(backlenb[inbestb] == 0){
													sprintf(cigarb, "1S%dM%dS", (forwlenb[inbestb]-1), (readlenb-forwlenb[inbestb]));
												}
												else if(gaplenb[inbestb] < 0){
													if(forwlenb[inbestb] > 3)
														sprintf(cigarb, "1S%dM%dD%dM", (forwlenb[inbestb]-1), -1*gaplenb[inbestb], backlenb[inbestb]);
													else{
														sprintf(cigarb, "%dS%dM", (readlenb-backlenb[inbestb]), backlenb[inbestb]);
														(data->resbufferb+p)->pos = forwlenb[inbestb]-gaplenb[inbestb];
													}
												}
												else{
													if(forwlenb[inbestb] > 3)
														sprintf(cigarb, "1S%dM%dI%dM", (forwlenb[inbestb]-1), gaplenb[inbestb], backlenb[inbestb]);
													else{
														sprintf(cigarb, "%dS%dM", (readlenb-backlenb[inbestb]), backlenb[inbestb]);
														(data->resbufferb+p)->pos = forwlenb[inbestb];
													}
												}	
												
												break;
											}
									}	
									
									if(gaplenb[inbestb] == 0){
										(data->resbufferb+p)->nmis = mapmisb[inbestb];
									}
									else if((forwlenb[inbestb] != 0)&&(backlenb[inbestb] != 0)){
										(data->resbufferb+p)->nmis = mapmisb[inbestb] - 1 - abs(gaplenb[inbestb])/10 + abs(gaplenb[inbestb]);
									}
									else{
										(data->resbufferb+p)->nmis = mapmisb[inbestb] - 1 - abs(gaplenb[inbestb])/10;
									}
									if((data->resbufferb+p)->nmis < 0)
										(data->resbufferb+p)->nmis = mapmisb[inbestb];
									(data->resbufferb+p)-> mapscore = scoreb;
									strcpy((data->resbufferb+p)->cigar, cigarb);
						}
						if((nmappos == 0)||(nmapposb == 0)){
							found = 0;
						}
						
						if((found == 1)&&(jca != jcb)){
								found = 0;
						}
						if((nmapposb > 0)&&((data->resbufferb+p)->pos == 0)){
							found = 0;
							ismapb = 1;
							directionb = 0;
							
						}
						if((nmappos > 0)&&((data->resbuffer+p)->pos == 0)){
							found = 0;
							ismapa = 1;
							directiona = 0;
							
						}
						
						if(found == 1){
							if(directiona == directionb){
								if(data->isfw == 0)
									found = 0;
							}
							else{
								if(data->isfw != 0)
									found = 0;
							}
						}
						
						if(found == 1){
							nmapped ++;
							if((score > 1)&&(scoreb > 1))
								nunimapped ++;
						}
						
						if(score > 1)
							nlunimapped ++;
						if(scoreb > 1)
							nrunimapped ++;
						
						(data->resbuffer+p)->nipos = niposa;
						(data->resbufferb+p)->nipos = niposb;
						(data->resbuffer+p)->flag = 1 + 2*found + ismapa*4 + ismapb*8 + 16*directiona + 32*directionb + 64;
						(data->resbufferb+p)->flag = 1 + 2*found + ismapb*4 + ismapa*8 + 16*directionb + 32*directiona + 128;
												
												
												(data->resbuffer+p)->jc = jca;
												(data->resbufferb+p)->jc = jcb;
												posa = (data->resbuffer+p)->pos;
												posb = (data->resbufferb+p)->pos;
												
												if((jca == jcb)&&(jcb > -1)){
													posa = (data->resbuffer+p)->pos;
													posb = (data->resbufferb+p)->pos;
													if(directiona == 0){
														if(directionb == 0){
															(data->resbuffer+p)->deltapos = posb - posa;
															if((data->resbuffer+p)->deltapos >= 0)
																(data->resbuffer+p)->deltapos += 1;
															else
																(data->resbuffer+p)->deltapos -= 1;
														}
														else{
															(data->resbuffer+p)->deltapos = posb + readlenb - posa;
														}
													}
													else{
														if(directionb == 0){
															(data->resbuffer+p)->deltapos = posb - posa - readlen;
														}
														else{
															(data->resbuffer+p)->deltapos = posb + readlenb - posa - readlen;
															if((data->resbuffer+p)->deltapos >= 0)
																(data->resbuffer+p)->deltapos += 1;
															else
																(data->resbuffer+p)->deltapos -= 1;
														}
													}
													if((data->resbuffer+p)->deltapos == 0)
														(data->resbuffer+p)->deltapos += 1;
													
													(data->resbufferb+p)->deltapos = -1*(data->resbuffer+p)->deltapos;
												}
						(data->resbuffer+p)->iscomp = true;
				}
				data->nmapped = nmapped;
				data->nunimapped = nunimapped;
				data->nlmapped = nlmapped;
				data->nrmapped = nrmapped;
				data->nlunimapped = nlunimapped;
				data->nrunimapped = nrunimapped;
				data->iscomp = true;
				
}

void* dataprint(void* arg){
	dargs* data = (dargs*) arg;
	int i, j, jca, jcb, posa, posb, deltaposa, deltaposb;
	bool isend = false;
	char psid[256];
	
	while(!isend){
		for(j = 0; j < nthreads; j++){
			if(((data->arglist+j)->resbuffer+(data->arglist+j)->ip)->iscomp && (!(data->arglist+j)->isprint)){	
				for(i = (data->arglist+j)->ip ; i < ((data->arglist+j)->ip+(data->arglist+j)->bufsize); i++){
					if(((data->arglist+j)->resbuffer+i)->iscomp && (!((data->arglist+j)->resbufferb+i)->iscomp)){	
							jca = ((data->arglist+j)->resbuffer+i)->jc;
							jcb = ((data->arglist+j)->resbufferb+i)->jc;
							posa = ((data->arglist+j)->resbuffer+i)->pos;
							posb = ((data->arglist+j)->resbufferb+i)->pos;
							deltaposa = ((data->arglist+j)->resbuffer+i)->deltapos;
							deltaposb = ((data->arglist+j)->resbufferb+i)->deltapos;
							if(isspace){
								sscanf(((data->arglist+j)->seqbuffer+i)->sid, "%s", psid);
								printf("%s\t", (psid+1));
							}
							else if(!isslash){
								printf("%s\t", (((data->arglist+j)->seqbuffer+i)->sid+1));
							}
							else
								printf("%*.*s\t", 2, (strlen(((data->arglist+j)->seqbuffer+i)->sid)-3), (((data->arglist+j)->seqbuffer+i)->sid+1));
							if((((data->arglist+j)->resbuffer+i)->flag/4)%2 == 1){
								if(jcb == -1)
									printf("%d\t*\t0\t0\t*\t*\t0\t0\t%s\t%s", ((data->arglist+j)->resbuffer+i)->flag, ((data->arglist+j)->seqbuffer+i)->seq, ((data->arglist+j)->seqbuffer+i)->qual);
								else
									printf("%d\t%s%lu\t0\t*\t=\t%lu\t0\t%s\t%s", ((data->arglist+j)->resbuffer+i)->flag, chrlist[jcb], posb, posb, ((data->arglist+j)->seqbuffer+i)->seq, ((data->arglist+j)->seqbuffer+i)->qual);
							}
							else{
									printf("%d\t%s%lu\t%d\t%s\t", ((data->arglist+j)->resbuffer+i)->flag, chrlist[jca], posa, ((data->arglist+j)-
