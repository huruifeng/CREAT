import sys
def viewBar(x):
    output = sys.stdout
    output.write('\rComplete percent: %.2f%%' % x)
    output.flush()

##########################################################
def read_fasta(fasta_dir):
    fasta_file = open(fasta_dir, 'r')
    names_seqs = {}
    seq = ''
    name = ''
    
    for this_line in fasta_file:
        if this_line.startswith('>'):
            if seq != '':
                names_seqs[name]=seq
                seq = ''
            name = this_line.strip()
        else:
            seq = seq + this_line.strip()
                
        this_line = fasta_file.readline()
    fasta_file.close()
    names_seqs[name]=seq ##Store the last record
    return names_seqs

############################################################
def read_fastq(fastq_path):
    fastq_file = open(fastq_path, 'r')
    names_seqs = {}
    read_tag = 0 
    for seq_record in fastq_file:
        if seq_record[0]=="@" and read_tag == 0:
            name = seq_record.strip()
            read_tag = 1
            continue
        if read_tag == 1:
            read_tag = 0
            seq = seq_record.strip()
            names_seqs[name] = seq
    fastq_file.close()
    return names_seqs
###############################################################
def printUsage():
    print "Usage:"
    print "    python kmerAnalysis.py -f input_file [-option value,...][-o output_file]\n"
    print "Note:"
    print "The input file could be Single-end or Peired-end reads file in fasta or fastq format."\
          "The output file contains two columns. The first column is the CoverageDepth distribution values of"\
          " all kmers. The second column is the kmer types Counts_Number corresponding to the CoverageDepth value."\
          " e.g. There are <Counts_Number> kmers appear <CoverageDepth> times in file <input>."
    print "\nOptions:"
    for option_i in option_dict:
        print '    '+option_i + "    " + option_dict[option_i]
###############################################################

## init
file_path = ""
kmer = 60
bin_size = 1
output = ""
thresold = 500
option_dict = {"-f":"The input file path(Peired-end files: The path of one of the files).",\
               "-o":"The output file path and name(Default:'result_kmer<k>bin<b>.txt",\
               "-k":"kmer size(Default:60).",\
               "-b":"bin size for kmer grouping(Default:1).",\
               "-t":"The reads in a group with a number of members that is less than this thresold value will be discarded(Default:1000).",\
               "-h":"Print the usage and option list.",\
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
            if running_i =="-o":
                output = running_dict[running_i]
            if running_i =="-k":
                kmer = int(running_dict[running_i])
            if running_i =="-t":
                thresold = int(running_dict[running_i])
            if running_i =="-b":
                bin_size = int(running_dict[running_i])
               
else:
    printUsage()
    sys.exit(1)


kmer_seq_count={}
read_tag = 0

reads = {}
print "Read sequence file..."
if file_path == "":
    print "ERROR:Input file is not specified."
    sys.exit(1)
else:
    file_type = file_path.split('.')[-1]

if file_type == "fq" or file_type == "fastq":
    reads = read_fastq(file_path)
elif file_type == "fa" or file_type == "fasta":
    reads = read_fasta(file_path)
else:
    fp_file = open(file_path, "r")
    fp_line = fp_file.readline()
    fp_file.close()
    if fp_line.strip()[0] == ">":
        reads = read_fasta(file_path)
    elif fp_line.strip()[0] == "@":
        reads = read_fastq(file_path)
    else:
        print "ERROR: File type or File format is wrong!"
        sys.exit(1)
print "Read file complete!"
##Obtain kmers
last_percent = 0.00
i = 0
len_count = len(reads)
print "kmer Counting..."
for read_i in reads:
    i+=1
    current_percent = round(float(i)/len_count,4)*100
    if current_percent != last_percent:
        viewBar(current_percent)
        last_percent = current_percent
    start_i = 0
    sequence = reads[read_i]
    seq_len = len(sequence)
    for start_i in range(seq_len+1-kmer):
        kmer_seq = sequence[start_i:start_i + kmer]
        if kmer_seq in kmer_seq_count:
            kmer_seq_count[kmer_seq] += 1
        else:
            kmer_seq_count[kmer_seq] = 1
kmer_count = []
print "\nkmer Counting complete !"
print "kmer sorting..."
for kmer_seq in kmer_seq_count:
    kmer_count.append(kmer_seq_count[kmer_seq])

kmer_count.sort()
#print kmer_count
min_count = kmer_count[0]
max_count = kmer_count[-1]

len_count = len(kmer_count)
 
#print min_count,max_count
print "kmer sorting complete"

print "K-mer grouping..."
print "####k-mer:"+str(kmer)+" ####"

print "####Bin_Size:"+str(bin_size)+" ####"
static = {}
k = min_count
i = 0
last_percent = 0.00
while k <= max_count and i < len_count:
    current_percent = round(float(i)/len_count,4)*100
    if current_percent != last_percent:
        viewBar(current_percent)
        last_percent = current_percent
    if kmer_count[i] >= k and kmer_count[i] < k+bin_size:
        if k in static:
            static[k] += 1
        else:
            static[k] = 1
        i = i + 1
    else:
        k = k+bin_size
        static[k] = 0

print "\nK-mer grouping completed !"


print "File Writing..."
if output =="":
    fp_r = open('result_kmer'+str(kmer)+"_bin"+str(bin_size)+".txt",'w')
else:
    fp_r = open(output,'w')
fp_r.write("kmerCoverageDepth\tKmer_Counts\n")
for r_i in static:
    fp_r.write(str(r_i)+"\t"+str(static[r_i])+"\n")

fp_r.close()
print "File Writing Completed !"
print "Calculating the thresold values..."
thresold_left = 0
static_sorted = sorted(static.iteritems(),key=lambda d:d[0])
for k_i in range(len(static_sorted)-1):
    if static_sorted[k_i][1]<static_sorted[k_i+1][1] and thresold_left == 0:
        thresold_left = static_sorted[k_i][0]
    if static_sorted[k_i][1]<thresold and static_sorted[k_i+1][1]<thresold:
        thresold_right = static_sorted[k_i][0]
        break
if thresold_left<=30 or thresold_left>=160 or thresold_left>=thresold_right:
    print "The thresold values can not be calculated, Please specify them manully by refer to the result file!"
else:
    print "Recommend value:["+str(thresold_left)+":"+str(thresold_right)+"]"
print "All Work have been Done !"

