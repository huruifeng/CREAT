/*
 * common.h
 *
 *  Created on: Dec 12, 2016
 *      Author: Ruifeng HU
 */

#ifndef COMMON_H_
#define COMMON_H_

#include <dirent.h>
#include <unistd.h>

#include <stdio.h>
#include <sys/time.h>
#include <pthread.h>

#include <fstream>
#include <string.h>
#include <sys/stat.h>
#include <errno.h>
#include <stdint.h>
#include <ctime>
#include <time.h>
#include <iostream>
#include <memory>
#include <sstream>
#include <algorithm>
#include <string>
#include <vector>
#include <list>
#include <set>
#include <map>
#include <ext/hash_map>
#include <stdlib.h>
#include <math.h>

using namespace std;
using namespace __gnu_cxx;

typedef set<int> SETINT;
typedef set<string> SETSTR;
typedef set<char> SETCHAR;


typedef list<int> LISTINT;
typedef list<string> LISTSTR;
typedef list<char> LISTCHAR;

typedef vector<int> VECTORINT;
typedef vector<string> VECTORSTR;
typedef vector<char> VECTORCHAR;

typedef map<int,string> MAPINT_STR;
typedef map<string,int> MAPSTR_INT;
typedef map<string,unsigned long> MAPSTR_UNLONG;
typedef map<char,int> MAPCHR_INT;
typedef map<string,string> MAPSTR_STR;

//////////////////////////////////////////////////////
extern MAPSTR_STR pe1_name_1st;
extern MAPSTR_STR pe1_name_2nd;
extern MAPSTR_STR pe1_name_3rd;
extern MAPSTR_STR pe1_name_4th;

extern MAPSTR_STR pe2_name_1st;
extern MAPSTR_STR pe2_name_2nd;
extern MAPSTR_STR pe2_name_3rd;
extern MAPSTR_STR pe2_name_4th;

extern int tag_type;

extern string tag_tail_1;
extern string tag_tail_2;
///////////////////////////////////////////////////////
/* run.cfg variable*/
extern int cfg_running_mode;//0
extern string cfg_reads_type;// paired-end
extern string cfg_file_type;//fastq
extern string cfg_SE_path;//data/SE_file.fastq
extern string cfg_PE_path_1;// data/PE_1.fastq
extern string cfg_PE_path_2;//data/PE_2.fastq
extern int cfg_Thresold_Value_Left;//-1
extern int cfg_Thresold_Value_Right;//-1
extern int cfg_reads_count;
extern int cfg_kmer;//65
extern int cfg_bin_size;// 1
extern string cfg_Nthread;//4
extern string cfg_spades_max_mem;//-1
extern string cfg_spades_options;//-k 33,53,63,73 --careful

extern VECTORSTR reads_in_thread[64];
extern map<string,VECTORSTR> kmers_of_read[64];
extern MAPSTR_UNLONG kmer_count[64];
extern MAPSTR_UNLONG Reads_MaxKmerCount[64];

extern MAPSTR_UNLONG kmer_count_all;
extern MAPSTR_UNLONG Reads_MaxKmerCount_all;


// Functions

string LTrim(const string& str) ;

string RTrim(const string& str) ;

string Trim(const string& str) ;

string get_time();
void logPrint(string logstr,int time_flag=1);
void logPrint_int(int);

void viewPercent(float);

MAPSTR_STR Readcfgfile(string cfg_filename);

int IsFileExist(const char* path);

vector<string> split(const string src, char separator);

void read_fasta(string file_path,MAPSTR_STR &);
void read_fastq(string,MAPSTR_STR &,MAPSTR_STR &,MAPSTR_STR &,MAPSTR_STR &,string &);
vector<VECTORSTR> pair_reads_filter(vector<VECTORSTR>,vector<VECTORSTR>);

void* exec_KmerSplit(void*);
void* exec_KmerCount(void*);

#endif /* COMMON_H_ */
