##########################################################
import sys
def read_fasta(fasta_dir):
    fasta_file = open(fasta_dir, 'r')
    names_seqs = {}
    names_lens = {}
    seq = ''
    name = ''
    
    for this_line in fasta_file:
        if this_line.startswith('>'):
            if seq != '':
                names_seqs[name]=seq
                names_lens[name] = len(seq)
                seq = ''
            name = this_line.strip()
        else:
            seq = seq + this_line.strip()
                
    fasta_file.close()
    names_seqs[name]=seq ##Store the last record
    names_lens[name] = len(seq)
    return [names_seqs,names_lens]

############################################################
def read_fastq(fastq_path):
    fastq_file = open(fastq_path, 'r')
    names_seqs = {}
    read_tag = 0
    names_lens = {}
    for seq_record in fastq_file:
        if seq_record[0]=="@" and read_tag == 0:
            name = seq_record.strip()
            read_tag = 1
            continue
        if read_tag == 1:
            read_tag = 0
            seq = seq_record.strip()
            names_seqs[name] = seq
            names_lens[name] = len(seq)
    fastq_file.close()
    return [names_seqs,names_lens]
##########################################################
###############################################################
def printUsage():
    print "Usage:"
    print "    python seqSelect.py -f input_file [-option value,...][-o output_dir]\n"
    print "Note:"
    print "The input file could be the result file of SPAdes(scaffolds.fsata)."\
          "This script can select the top N sequences,and write them in one file(-t)."\
          "or, separate the top N sequences into n files(-s)."
    print "\nOptions:"
    for option_i in option_dict:
        print '    '+option_i + "    " + option_dict[option_i]
###############################################################

## init
file_path = ""
s_n = 5
t_n = 5
output = ""

option_dict = {"-f":"The input file path(Peired-end files: The path of one of the files).",\
               "-s":"select the top n longest sequences and write into n files with file name seq_<n>.fa(Default:5).",\
               "-t":"select the top n longest sequences and write into one file with file name top_<n>.fa(Default:5).",\
               "-h":"Print the usage and option list.",\
               "-o":"The output directory",\
               "-v":"Print the version of PEAK."}

if len(sys.argv) >= 2:
    option_list = sys.argv[1:]
    print option_list
    running_dict = {}
    if len(option_list)==1:
        if option_list[0] =="-h":
            printUsage()
            sys.exit(1)
        elif option_list[0] =="-v":
            print "PEAK V17.01.24"
            sys.exit(1)
        else:
            print "ERROR!!!"
            printUsage()
            sys.exit(1)
    for option_i in range(0,len(option_list),2):
        running_dict[option_list[option_i]] = option_list[option_i+1]
    for running_i in running_dict:
        if running_i not in option_dict:
            print "Option "+ running_i +" is invalid. It is not one of the Options."
            sys.exit(1)
        else:
            ##TODO
            
            if running_i =="-f":
                file_path = running_dict[running_i]
            if running_i =="-s":
                s_n = int(running_dict[running_i])
            if running_i =="-t":
                t_n = int(running_dict[running_i])
            if running_i =="-o":
                output = running_dict[running_i]
               
else:
    printUsage()
    sys.exit(1)

kmer_seq_count={}
read_tag = 0

reads = {}

if file_path == "":
    print "ERROR:Input file is not specified."
    sys.exit(1)
else:
    file_type = file_path.split('.')[-1]

if file_type == "fq" or file_type == "fastq":
    reads_length = read_fastq(file_path)
elif file_type == "fa" or file_type == "fasta":
    reads_length = read_fasta(file_path)
else:
    fp_file = open(file_path, "r")
    fp_line = fp_file.readline()
    fp_file.close()
    if fp_line.strip()[0] == ">":
        reads_length = read_fasta(file_path)
    elif fp_line.strip()[0] == "@":
        reads_length = read_fasta(file_path)
    else:
        print "ERROR: File type or File format is wrong!"
        sys.exit(1)

reads = reads_length[0]
reads_len = reads_length[1]

##
i = 0
reads_sorted = sorted(reads_len.iteritems(),key=lambda d:d[1], reverse = True)
if "-s" in running_dict:
    for read_i in reads_sorted:
        if output != "":
            output= output.strip("").strip("\\").strip("/").strip()
            fp = open(output+"/seq_"+str(i+1)+".fa",'w')
        else:
            fp = open("seq_"+str(i+1)+".fa",'w')
        fp.write(read_i[0]+"\n")
        sequence = reads[read_i[0]]
        set_len = 60
        while True:
            if len(sequence) >=set_len:
                sub_seq = sequence[0:set_len]
                fp.write(sub_seq+"\n")
                sequence=sequence[set_len:]
            else:
                sub_seq = sequence
                fp.write(sub_seq+"\n")
                break    
        fp.close()
        i = i + 1
        if i == s_n:
            break
i = 0
if "-t" in running_dict:
    fp = open("top"+str(t_n)+".fa",'w')
    if output != "":
        output= output.strip("").strip("\\").strip("/").strip()
        fp = open(output+"/top"+str(t_n)+".fa",'w')
    else:
        fp = open("top"+str(t_n)+".fa",'w')
    for read_i in reads_sorted:
        print read_i[0]
        fp.write(read_i[0]+"\n")
        sequence = reads[read_i[0]]
        set_len = 60
        while True:
            if len(sequence) > set_len:
                sub_seq = sequence[0:set_len]
                fp.write(sub_seq+"\n")
                sequence=sequence[set_len:]
            else:
                sub_seq = sequence
                fp.write(sub_seq+"\n")
                break    
        
        i = i + 1
        if i == t_n:
            break
    fp.close()

