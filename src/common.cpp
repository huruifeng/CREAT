/*
 * common.cpp
 *
 *  Created on: Dec 12, 2016
 *      Author: Rufeng Hu
 */

#include "common.h"
/* Global variable defination;*/
MAPSTR_STR pe1_name_1st;
MAPSTR_STR pe1_name_2nd;//common_name:sequence
MAPSTR_STR pe1_name_3rd;
MAPSTR_STR pe1_name_4th;

MAPSTR_STR pe2_name_1st;
MAPSTR_STR pe2_name_2nd;//common_name:sequence
MAPSTR_STR pe2_name_3rd;
MAPSTR_STR pe2_name_4th;

int tag_type = 1; //1:space; 2:#

string tag_tail_1;
string tag_tail_2;
/* Global  variable end;*/

/* run.cfg variable*/
int cfg_running_mode;//0
string cfg_reads_type;// paired-end
string cfg_file_type;//fastq
string cfg_SE_path;//data/SE_file.fastq
string cfg_PE_path_1;// data/PE_1.fastq
string cfg_PE_path_2;//data/PE_2.fastq
int cfg_Thresold_Value_Left;//-1
int cfg_Thresold_Value_Right;//-1
int cfg_reads_count; //1000
int cfg_kmer;//65
int cfg_bin_size;// 3
string cfg_Nthread;//4
string cfg_spades_max_mem;//-1
string cfg_spades_options;//-k 33,53,63,73 --careful

VECTORSTR reads_in_thread[64];
map<string,VECTORSTR> kmers_of_read[64];
MAPSTR_UNLONG kmer_count[64];
MAPSTR_UNLONG Reads_MaxKmerCount[64];

MAPSTR_UNLONG kmer_count_all;
MAPSTR_UNLONG Reads_MaxKmerCount_all;
/////////////////////////////////////////////////////////
string LTrim(const string& str)
 {
 return str.substr(str.find_first_not_of(" \n\r\t"));
 }

string RTrim(const string& str)
 {
 return str.substr(0,str.find_last_not_of(" \n\r\t")+1);
 }

string Trim(const string& str)
 {
 return LTrim(RTrim(str));
 }

void logPrint(string logstr,int time_flag){
	ofstream logFile("temp/PEAK_logs.txt",ios::app);
	if(time_flag==1)
		logFile<<get_time();
	logFile<<logstr;
	logFile.close();
}

void logPrint_int(int logstr){
	ofstream logFile("temp/PEAK_logs.txt",ios::app);
	logFile<<logstr;
	logFile.close();
}


vector<string> split(const string& str, string sep)
{
    if (str.empty())
    {
    	cout<<"******************************************"<<endl;
    	cout<<"!!!ERROR::split_string, String is EMPTY !!!"<<endl;
    	cout<<"******************************************"<<endl;
    	cout<<get_time()<<" Pipeline Stopped and Exited !"<<endl;

    	logPrint(" !!!ERROR::split_string, String is EMPTY !!!\n");
    	logPrint(" Pipeline Stopped and Exited !\n");

    	exit(1);
    }

    if(sep==""){
    	sep = " ";
    }

    vector<string> ret_;
    string tmp;
    string::size_type pos_begin = str.find_first_not_of(sep);
    string::size_type comma_pos = 0;

    while (pos_begin != string::npos)
    {
        comma_pos = str.find(sep, pos_begin);
        if (comma_pos != string::npos)
        {
            tmp = str.substr(pos_begin, comma_pos - pos_begin);
            pos_begin = comma_pos + sep.length();
        }
        else
        {
            tmp = str.substr(pos_begin);
            pos_begin = comma_pos;
        }

        if (!tmp.empty())
        {
            ret_.push_back(tmp);
            tmp.clear();
        }
    }
    return ret_;
}



MAPSTR_STR Readcfgfile(string cfg_filename){
	MAPSTR_STR results;
	ifstream cfg;
	VECTORSTR ls;
	string in_text,temp_line;
	cfg.open(cfg_filename.c_str());
	if(cfg.is_open()){
		while(getline(cfg,temp_line)){
			 temp_line=RTrim(temp_line);
			 //cout<<"AA"<<temp_line<<"AA"<<endl; //for test
			 if((temp_line[0]=='#') or temp_line==""){
				continue;
			 }
			 else{
				 ls.clear();
				 ls = split(temp_line,"=");
				 results[Trim(ls[0])]=Trim(ls[1]);
			 }
		}
	}
	else{
		cout<<"######################################\n";
		cout<<"ERROR:*.cfg file open goes wrong !!!\n";
		cout<<"######################################\n";
		cout<<get_time()<<" Pipeline Stopped and Exited !"<<endl;
		logPrint(" ERROR:*.cfg file open goes wrong !!!\n");
		logPrint(" Pipeline Stopped and Exited !\n");
		exit(1);
	}
	cfg.close();
	for(MAPSTR_STR::iterator Itr=results.begin();Itr!=results.end();++Itr){
			if(Itr->first=="reads_type"){
				cfg_reads_type = Itr->second;
			}
			else if(Itr->first=="running_mode"){
				cfg_running_mode = atoi(Itr->second.c_str());
			}
			else if(Itr->first=="file_type"){
				cfg_file_type = Itr->second;
			}
			else if(Itr->first=="SE_path"){
				cfg_SE_path = Itr->second;
			}
			else if(Itr->first=="PE_path_1"){
				cfg_PE_path_1 = Itr->second;
			}
			else if(Itr->first=="PE_path_2"){
				 cfg_PE_path_2 = Itr->second;
			}
			else if(Itr->first=="Thresold_Value_Left"){
				cfg_Thresold_Value_Left = atoi(Itr->second.c_str());
			}
			else if(Itr->first=="Thresold_Value_Right"){
				cfg_Thresold_Value_Right = atoi(Itr->second.c_str());
			}
			else if(Itr->first=="reads_count"){
				cfg_reads_count = atoi(Itr->second.c_str());
			}
			else if(Itr->first=="kmer"){
				cfg_kmer = atoi(Itr->second.c_str());
			}
			else if(Itr->first=="bin_size"){
				cfg_bin_size=atoi(Itr->second.c_str());
			}
			else if(Itr->first=="Nthread"){
				cfg_Nthread = Itr->second;
			}
			else if(Itr->first=="spades_max_mem"){
				cfg_spades_max_mem = Itr->second;
			}
			else if(Itr->first=="spades_options"){
				cfg_spades_options = Itr->second;
			}
			else{
				cout<<"#######################################################"<<endl;
				cout<<"ERROR:*.cfg, There is no variable '"<<Itr->first<<"'!"<<endl;
				cout<<"Please do not change the default variable name !"<<endl;
				cout<<"#######################################################"<<endl;
				logPrint(" ERROR:*.cfg, There is no variable '");
				logPrint(Itr->first,0);
				logPrint("'!\n",0);
				logPrint(" Please do not change the default variable name !\n");
				exit(1);
			}
		}
	return results;
}



string get_time(){
	string str;
	time_t t0 = time(NULL);
	tm *t = localtime(&t0);
	ostringstream oss_h,oss_m,oss_s;
	oss_h << t->tm_hour;
	oss_m << t->tm_min;
	oss_s << t->tm_sec;
	str = oss_h.str()+":"+oss_m.str()+":"+oss_s.str()+" ";
	return str;

}

void viewPercent(float x){
	printf("  Complete Percent: %.2f%% \r",x);
	fflush(stdout);
}




int IsFileExist(const char* path)
{
    return !access(path, F_OK);
}

void read_fastq(string file_path,MAPSTR_STR &pe_name_1st,MAPSTR_STR &pe_name_2nd,MAPSTR_STR &pe_name_3rd,MAPSTR_STR &pe_name_4th,string &tag_tail){

	ifstream fp_in;

	string temp_line,temp_name,temp_read,common_name;

	MAPSTR_STR	name_read;

	fp_in.open(file_path.c_str());
	if(!fp_in){
		cout<<"********************************************************************************"<<endl;
		cout<<"!!!ERROR::open fastq file, File Open Failed ! Check the File path and names !!!"<<endl;
		cout<<"********************************************************************************"<<endl;
		cout<<get_time()<<" Pipeline Stopped and Exited !"<<endl;
		logPrint(" !!!ERROR::open fastq file, File Open Failed ! Check the File path and names !!!\n");
		logPrint(" Pipeline Stopped and Exited !\n");
		exit(1);
	}

	bool set_tag_type = false;
	size_t tag_pos;
	while(getline(fp_in,temp_line)){
		if (temp_line[0] != '@'){
			cout<<"*************************************************************************"<<endl;
			cout<<"!!!ERROR::fastq file format, File format error ! Check the File content !"<<endl;
			cout<<"*************************************************************************"<<endl;
			cout<<get_time()<<" Pipeline Stopped and Exited !"<<endl;
			logPrint(" !!!ERROR::fastq file format, File format error ! Check the File content !!!\n");
			logPrint(" Pipeline Stopped and Exited !\n");
			exit(1);
		}else{
			//cout<<Trim(temp_line)<<endl; // For test
			temp_name = Trim(temp_line);
			if(!set_tag_type){
				tag_pos = temp_name.find_first_of(" ");
				if(tag_pos!=string::npos){
					tag_type = 1;
					tag_tail = temp_name.substr(tag_pos);
				}
				else{
					tag_type = 2;
					tag_pos = temp_name.find_first_of("#");
					tag_tail = temp_name.substr(tag_pos);
				}
				set_tag_type=true;
			}
			//cout<<tag_tail<<endl; //For test
			//cout<<temp_name<<endl; //For test
			common_name = temp_name.substr(0,tag_pos);
			pe_name_1st[common_name] = temp_name;
			//cout<<temp_name.substr(0,tag_pos)<<endl; //For test

		}
		getline(fp_in,temp_line);
		pe_name_2nd[common_name] = Trim(temp_line);
		//cout<<Trim(temp_line)<<endl; //For test

		getline(fp_in,temp_line);
		pe_name_3rd[common_name] = Trim(temp_line);
		//cout<<Trim(temp_line)<<endl; //For test

		getline(fp_in,temp_line);
		pe_name_4th[common_name] = Trim(temp_line);
		//cout<<Trim(temp_line)<<endl; //For test

	}//end of while
	fp_in.close();
}

void read_fasta(string file_path,MAPSTR_STR &pe_name_2nd){

	ifstream fp_in;

	string temp_line,temp_name,temp_read="None";

	fp_in.open(file_path.c_str());

	if(!fp_in){
		cout<<"*************************************************************************"<<endl;
		cout<<"!!!ERROR::open file, File Open Failed ! Check the File path and names !!!"<<endl;
		cout<<"*************************************************************************"<<endl;
		cout<<get_time()<<" Pipeline Stopped and Exited !"<<endl;
		logPrint(" !!!ERROR::open file, File Open Failed ! Check the File path and names !!!\n");
		logPrint(" Pipeline Stopped and Exited !\n");
		exit(1);
	}

	while(getline(fp_in,temp_line)){
		if (temp_line[0] == '>'){
			if(temp_read!="None"){
				pe_name_2nd[temp_name] = temp_read;
				temp_read = "";
			}
			temp_name = Trim(temp_line);
		}
		else if(temp_line[0] != '>'){
			temp_read += Trim(temp_line);
		}
		else{

		}

	}//end of while
	//Save the last one read
	pe_name_2nd[temp_name] = temp_read;

	fp_in.close();
}

vector<VECTORSTR>pair_reads_filter(vector<VECTORSTR> pe1_names_reads,vector<VECTORSTR> pe2_names_reads){
	cout<<get_time()<<" Paired-Reads Filtering..."<<endl;
	logPrint(" Paired-Reads Filtering...\n");

	vector<VECTORSTR> result;

	MAPSTR_INT names;
	MAPSTR_INT::iterator name_iter;

	unsigned long i = 0;

	VECTORSTR pe1_names,pe2_names,pe1_reads,pe2_reads;
	pe1_names = pe1_names_reads[0];
	pe2_names = pe2_names_reads[0];
	pe1_reads = pe1_names_reads[1];
	pe2_reads = pe2_names_reads[1];

	string temp_name;
	unsigned long len_pe1_names = pe1_names.size();
	//cout<<len_pe1_names<<endl;//for test
	for(i=0;i<len_pe1_names;i++){
		names[split(split(pe1_names[i]," ")[0],"#")[0]]=i;
		//cout<<split(split(pe1_names[i]," ")[0],"#")[0]<<endl;// For test
	}

	VECTORSTR outpe_1_names,outpe_1_seqs,outpe_2_names,outpe_2_seqs;
	string this_name;
	unsigned long len_pe2_names = pe2_names.size();
	unsigned long id_num;
	for(i=0;i<len_pe2_names;i++){
		this_name = split(split(pe2_names[i]," ")[0],"#")[0];
		name_iter = names.find(this_name);
		if(name_iter != names.end()){
			id_num = names[this_name];
			outpe_1_names.push_back(pe1_names[id_num]);
			outpe_1_seqs.push_back(pe1_reads[id_num]);
			outpe_2_names.push_back(pe2_names[i]);
		    outpe_2_seqs.push_back(pe2_reads[i]);

		    //for test
		    /**
		    cout<<pe1_names[id_num]<<endl;
		    cout<<pe1_reads[id_num]<<endl;
		    cout<<pe2_names[i]<<endl;
		    cout<<pe2_reads[i]<<endl;
		    **/
		}
	}
	result.push_back(outpe_1_names);
	result.push_back(outpe_1_seqs);
	result.push_back(outpe_2_names);
	result.push_back(outpe_2_seqs);
	return result;
}


/////////////////////////////////////////////////////////////////////////////////////
void kmerSplit_in_thread(VECTORSTR temp_name,map<string,VECTORSTR> &kmers_of_read_temp,MAPSTR_UNLONG &kmer_count_temp){

	int kmer_len = cfg_kmer;
	//cout<<"kmer_len:"<<kmer_len<<endl; //for test

	kmer_count_temp.clear();
	MAPSTR_UNLONG::iterator kmer_it;

	VECTORSTR kmer_list;
	VECTORSTR::iterator kmer_list_it;

	map<string,VECTORSTR>::iterator kmers_of_read_it;

	string sequence,kmer_seq;

	int start_i,seq_len;

	//cout<<len_names<<endl; //for test
	VECTORSTR::iterator vecstr_it;

	for(vecstr_it=temp_name.begin();vecstr_it!=temp_name.end();vecstr_it++){
		//get kmers
		sequence = pe1_name_2nd[*vecstr_it];
		//cout<<sequence<<endl; //for test
		seq_len = sequence.length();
		//cout<<"Seq_length:"<<seq_len<<endl; //for test
		kmer_list.clear();

		for(start_i = 0;start_i < seq_len-kmer_len+1;start_i ++){
			kmer_seq = sequence.substr(start_i,kmer_len);
			kmer_list.push_back(kmer_seq);

			kmer_it = kmer_count_temp.find(kmer_seq);
			if(kmer_it != kmer_count_temp.end()){
				kmer_count_temp[kmer_seq] += 1;
			}else{
				kmer_count_temp[kmer_seq] = 1;
			}
		}
		kmers_of_read_temp[*vecstr_it] = kmer_list;
		//cout<<"kmer_list_size:"<<kmer_list.size()<<endl; //for test
	}//end of for::mapss_it
}


void* exec_KmerSplit(void* arg){
	int x = *(int *)arg;
	cout<<"K-merSplitting in thread:"<<x<<endl;
	kmerSplit_in_thread(reads_in_thread[x],kmers_of_read[x],kmer_count[x]);
	cout<<"K-merSplitting in thread:"<<x<<" Completed!"<<endl;

	pthread_exit((void*)0);
}


void kmerCount_in_thread(MAPSTR_UNLONG &Reads_MaxKmerCount_temp,map<string,VECTORSTR> kmers_of_read_temp,MAPSTR_UNLONG kmer_count_temp){
	//TODO
	unsigned long max_num;
	VECTORSTR kmer_list;
	VECTORSTR::iterator kmer_list_it;
	MAPSTR_UNLONG::iterator kmer_num_it;
	map<string,VECTORSTR>::iterator kmers_of_read_it;
	for(kmers_of_read_it=kmers_of_read_temp.begin();kmers_of_read_it!=kmers_of_read_temp.end();kmers_of_read_it++){
		kmer_list=kmers_of_read_it->second;
		max_num = 0;
		for(kmer_list_it=kmer_list.begin();kmer_list_it!=kmer_list.begin();kmer_list_it++){
			if(kmer_count_all[*kmer_list_it ]> max_num){
				max_num = kmer_count_all[*kmer_list_it];
			}

		}
		Reads_MaxKmerCount_temp[kmers_of_read_it->first]=max_num;
	}
}


void* exec_KmerCount(void* arg){
	int x = *(int *)arg;
	cout<<"K-merCounting in thread:"<<x<<endl;
	kmerCount_in_thread(Reads_MaxKmerCount[x],kmers_of_read[x],kmer_count_all);
	cout<<"K-merCounting in thread:"<<x<<" Completed!"<<endl;

	pthread_exit((void*)0);
}




