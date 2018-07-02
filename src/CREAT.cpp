//===========================================================
// Name        : CREAT.cpp
// Author      : Ruifeng Hu
// Version     : v17
// Copyright   : copyright @ Ruifeng Hu, IMPLAD
// Description : C++, Ansi-style: assembling complete plastid genome
//               without reference from genome skimming data
//===========================================================

#include "CREAT.h"
using namespace std;

string version="V17.08.31_2";
void printLogo(){
	cout<<"##############################################################\n";
	cout<<"                Welcome to use CREAT(Ver:"<<version<<")    \n";
	cout<<"##############################################################\n";
	cout<<"      *****     *******    ********       **     ************ \n";
	cout<<"    ****  **    ********   ********      ****    ************ \n";
	cout<<"   **      **   **     **  **           **  **        **      \n";
	cout<<"  **            **     **  **          **    **       **      \n";
	cout<<"  **            ********   *******    **********      **      \n";
	cout<<"  **            *******    *******   ************     **      \n";
	cout<<"  **            **  **     **        **        **     **      \n";
	cout<<"   **      **   **   **    **        **        **     **      \n";
	cout<<"    ****  **    **    **   ********  **        **     **      \n";
	cout<<"      *****     **     **  ********  **        **     **      \n";
	cout<<"##############################################################\n";
	cout<<"\nChloroplast Reads Enrichment and Assembly Toolkits(CREAT)\n";
	cout<<"Assembling Complete Chloroplast Genome without Reference from Genome Skimming Data\n\n";
}

int main(int argc, char* argv[]) {

	printLogo();
	string pathname = "data/run.cfg";
	//cout<<argc<<endl;
	int n=0;
	while(n<argc){
		cout<<argv[n++]<<" ";
	}
	cout<<endl;

	if (argc <= 1){
		pathname = "data/run.cfg";
	}else if(strcmp(argv[1],"-h")==0){
		cout<<"Usage:"<<endl;
		cout<<"    ./CREAT [dir/run.cfg]"<<endl;
		cout<<"     The default path to *.cfg is 'data/run.cfg'."<<endl;
		cout<<endl;
		return 0;
	}
	else{
		pathname = argv[1];
	}
	cout<<pathname<<endl;
	MAPSTR_STR cfgs = Readcfgfile(pathname);

	//Create the folder;
	string clean_temp_cmd;
	if (IsFileExist("temp")){
		clean_temp_cmd = "rm -rf temp";
		system(clean_temp_cmd.c_str());
	}
	mkdir("temp",0777);

	if (IsFileExist("output")){
		clean_temp_cmd = "rm -rf output";
		system(clean_temp_cmd.c_str());
	}
	mkdir("output",0777);

	logPrint("log file: temp/CREAT_logs.txt\n");
	logPrint("Running Options:\n",0);
	logPrint(cfg_PE_path_1,0);
	logPrint("\n",0);
	logPrint(cfg_PE_path_2,0);
	logPrint("\n",0);
	logPrint("File Type:",0);
	logPrint(cfg_file_type,0);
	logPrint("\n",0);
	logPrint("Reads Type:",0);
	logPrint(cfg_reads_type,0);
	logPrint("\n",0);

	cout<<"Running Options:\n";
	cout<<"Data Path:"<<cfg_PE_path_1<<"; "<<cfg_PE_path_2<<endl;
	cout<<"File Type:"<<cfg_file_type<<endl;
	cout<<"Reads Type:"<<cfg_reads_type<<endl;
	cout<<"Threads Num:"<<cfg_Nthread<<endl;



	cout<<"\nStart the Job ? (Y/N):";
	char Y = getchar();
	//char Y = 'Y';
	if(Y=='Y'||Y=='y'){
		cout<<get_time()<<" Pipeline Starting..."<<endl;
		logPrint(" Pipeline Starting...\n");
	}
	else{
		cout<<get_time()<<" Pipeline Exited Manually !"<<endl;
		logPrint(" Pipeline Exited Manually !\n");
		exit(1);
	}

	/////////////////////////////////////////////////
	string file_1_path,file_2_path;

	if(cfg_reads_type=="paired-end"){
		file_1_path=cfg_PE_path_1;
		file_2_path=cfg_PE_path_2;
	}
	else if(cfg_reads_type=="single-end"){
		file_1_path=cfg_SE_path;
	}
	// Define data path
	//string file_1_path="data/test1.fastq";
	//string file_2_path="data/test2.fastq";



	if(cfg_file_type == "fastq" || cfg_file_type == "fq"){
		// reads from file 1
		cout<<get_time()<<" File Reading 1..."<<endl;
		logPrint(" File Reading 1...\n");
		cout<<cfg_PE_path_1<<endl;
		read_fastq(file_1_path,pe1_name_1st,pe1_name_2nd,pe1_name_3rd,pe1_name_4th,tag_tail_1);

		if(cfg_reads_type=="paired-end"){
			// reads from file 2
			cout<<get_time()<<" File Reading 2..."<<endl;
			logPrint(" File Reading 2...\n");
			read_fastq(file_2_path,pe2_name_1st,pe2_name_2nd,pe2_name_3rd,pe2_name_4th,tag_tail_2);
			cout<<cfg_PE_path_2<<endl;
		}
	}//end of if

	else if(cfg_file_type == "fasta" || cfg_file_type == "fa"){
		// reads from file 1
		read_fasta(file_1_path,pe1_name_2nd);

		if(cfg_reads_type=="paired-end"){
			// reads from file 2
			read_fasta(file_2_path,pe2_name_2nd);
		}

	}else{
		cout<<"****************************************************"<<endl;
		cout<<"!!!ERROR::file_type, Wrong File Type Setting::"<<cfg_file_type<<endl;
		cout<<"****************************************************"<<endl;
		logPrint(" !!!ERROR::file_type, Wrong File Type Setting:");
		logPrint(cfg_file_type);
		logPrint("\n",0);
	}
	cout<<get_time()<<" File Reading Completed !"<<endl;
	logPrint(" File Reading Completed !\n");

	cout<<get_time()<<" K-mer Splitting..."<<endl;
	logPrint(" K-mer Splitting...\n");
	unsigned long num_per_thread,n_thread = atoi(cfg_Nthread.c_str()),i,k;
	num_per_thread = 1 + pe1_name_1st.size()/n_thread;

	cout<<"Reads Sum:"<<pe1_name_1st.size()<<endl;

	MAPSTR_STR::iterator mss_it;
	i = 0;k=0;
	for(mss_it=pe1_name_1st.begin();mss_it!=pe1_name_1st.end();mss_it++){
		i++;
		if(i>num_per_thread){
			k += 1;
			i = 1;
		}
		reads_in_thread[k].push_back(mss_it->first);
	}


	int ret;
	vector<pthread_t> vecThreadId;
	pthread_t threadId;
	for(i = 0;i < n_thread;i++){
		ret= pthread_create(&threadId,NULL, exec_KmerSplit, (void *)&i);
		sleep(1);
		//cout<<"Sleeping...!"<<i<<endl; // for test
		if(ret!=0){
			cout<<"ERROR:kmerSplit_Thread_Creat: "<<i<<", error,ErrNo:"<<ret<<endl;
		}
		else{
			vecThreadId.push_back(threadId);
		}

	}

	for(vector<pthread_t>::iterator it = vecThreadId.begin(); it != vecThreadId.end(); ++it){
		pthread_join(*it, NULL);
	}
	cout<<get_time()<<" K-mer Splitting Completed !"<<endl;
	logPrint(" K-mer Splitting Completed !\n");
	cout<<"Reads Number in per-Thread:\n";
	for(i=0;i<n_thread;i++){
		cout<<"    "<<i<<"-->"<<kmers_of_read[i].size()<<endl;
	}

//////////////////////////////////////////////////////////
	cout<<get_time()<<" K-mer Counting......"<<endl;
	logPrint(" K-mer Counting......\n");

	MAPSTR_UNLONG::iterator msul_it,kmer_it;
	kmer_count_all.clear();
	for(i = 0;i < n_thread;i++){
		// k-mer Summarizing
		cout<<"Summarizing in Thread:"<<i<<endl;//for test
		for(msul_it=kmer_count[i].begin();msul_it!=kmer_count[i].end();msul_it++){
			kmer_it = kmer_count_all.find(msul_it->first);
			if(kmer_it != kmer_count_all.end()){
				kmer_count_all[msul_it->first] += msul_it->second;
			}else{
				kmer_count_all[msul_it->first] = msul_it->second;
			}

		}//kmer_count[i]
		//kmer_count[i].clear();
	}//n_thread
	cout<<"kmerCounting in Thread..."<<endl;
	vecThreadId.clear();
	for(i = 0;i < n_thread;i++){
		ret= pthread_create(&threadId,NULL, exec_KmerCount, (void *)&i);
		sleep(1);
		//cout<<"Sleeping...!"<<i<<endl; // for test
		if(ret!=0){
			cout<<"ERROR:kmerCount_Thread_Creat: "<<i<<", error,ErrNo:"<<ret<<endl;
		}
		else{
			vecThreadId.push_back(threadId);
		}

	}

	for(vector<pthread_t>::iterator it = vecThreadId.begin(); it != vecThreadId.end(); ++it){
		pthread_join(*it, NULL);
	}

	cout<<"Wriring in file..."<<endl;
	ofstream out_file;
	out_file.open("output/reads_kmerMaxCount.txt");
	map<long,long> reads_num;
	map<long,long>::iterator mll_it;
	for(i = 0;i < n_thread;i++){
		// k-mer Counting
		for(msul_it=Reads_MaxKmerCount[i].begin();msul_it!=Reads_MaxKmerCount[i].end();msul_it++){
			Reads_MaxKmerCount_all[msul_it->first] = msul_it->second;

			mll_it = reads_num.find(msul_it->second);
			if(mll_it != reads_num.end()){
				reads_num[msul_it->second] += 1;
			}else{
				reads_num[msul_it->second] = 1;
			}

			out_file<<msul_it->first<<"\t"<<msul_it->second<<"\n";
		}
		//Reads_MaxKmerCount[i].clear();
	}//n_thread
	out_file.close();
	cout<<"Grouping..."<<endl;
	out_file.open("output/reads_groupCount.txt");
	for(mll_it = reads_num.begin();mll_it != reads_num.end();mll_it++){
		out_file<<mll_it->first<<"\t"<<mll_it->second<<"\n";
	}
	out_file.close();

	cout<<get_time()<<" K-mer Counting Complete!"<<endl;
	logPrint(" K-mer Counting Complete!\n");

	exit(0);

//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

	if(cfg_running_mode==0){
		cout<<get_time()<<" Filtered Reads have been write into temp folder !\n";
		cout<<get_time()<<" Genome Assembly will be started !\n";
		logPrint(" Filtered Reads have been write into temp folder !\n");
		logPrint(" Genome Assembly will be started !\n");
	}else{
		cout<<get_time()<<" Filtered Reads have been write into temp folder !\n";
		cout<<get_time()<<" The reads can be assembly by other programs !\n";
		cout<<get_time()<<" All Work have been Done !"<<endl;
		cout<<get_time()<<" CREAT Pipeline Finished !"<<endl;
		cout<<get_time()<<" Thank you for using CREAT !"<<endl;
		logPrint(" Filtered Reads have been write into temp folder !\n The reads can be assembly by other programs !\n");
		logPrint(" All Work have been Done !\n CREAT Pipeline Finished !\n Thank you for using CREAT !\n");
		return 0;
	}
	//////////////////////////////////////////////////////////
	cout<<get_time()<<" Genome Assembly starting..."<<endl;
	logPrint(" Genome Assembly starting...\n");
	cout<<get_time()<<" Genome Assembling..."<<endl;
	logPrint(" Genome Assembling...\n");
	string command_option = "";
	if(atoi(cfg_Nthread.c_str())>0 ){
		command_option = " -t " + cfg_Nthread+" ";
	}else{
		command_option = " -t 4 ";
	}

	if(atoi(cfg_spades_max_mem.c_str())>0 ){
		command_option = " -m " + cfg_spades_max_mem+" ";
	}else{
		command_option = " -m 64 ";
	}

	if(cfg_file_type == "fasta" || cfg_file_type == "fa"){
		 command_option=command_option + cfg_spades_options +" --only-assembler ";
	 }
	 if(cfg_file_type == "fastq" || cfg_file_type == "fq"){
	 	command_option= command_option + cfg_spades_options + " ";
	 }
	string command_reads;
	if(cfg_reads_type=="paired-end"){
		if(cfg_file_type == "fastq" || cfg_file_type == "fq")
			command_reads= "-1 temp/paired_reads_1_filtered.fq -2 temp/paired_reads_2_filtered.fq";
		if(cfg_file_type == "fasta" || cfg_file_type == "fa")
			command_reads= "-1 temp/paired_reads_1_filtered.fa -2 temp/paired_reads_2_filtered.fa";
	}
	if(cfg_reads_type=="single-end"){
		if(cfg_file_type == "fastq" || cfg_file_type == "fq")
			command_reads= "-s temp/single_reads_filtered.fq";
		if(cfg_file_type == "fasta" || cfg_file_type == "fa")
			command_reads= "-s temp/single_reads_filtered.fa";
	}

	string SPAdes_Base = getcwd(NULL,0);
	//cout<<SPAdes_Base<<endl;//for test
	string command_str ="python "+SPAdes_Base+"/SPAdes/bin/spades.py "+command_option+command_reads+" -o output";
	logPrint(command_str);
	logPrint("\n",0);
	//cout<<command_str<<endl;//for test
	system(command_str.c_str());

	cout<<get_time()<<" Genome Assembly completed !"<<endl;
	logPrint(" Genome Assembly completed !\n");

	cout<<get_time()<<" All Work have been Done !"<<endl;
	logPrint(" All Work have been Done !\n");
	cout<<get_time()<<" CREAT Pipeline Finished !"<<endl;
	logPrint(" CREAT Pipeline Finished !\n");

	cout<<" Thank you for using CREAT !"<<endl;
	logPrint(" Thank you for using CREAT !\n");

	return 0;
}
