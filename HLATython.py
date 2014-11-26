import sys,os,glob,re,shutil
import math
import ConfigParser

# Open configuration file

config = ConfigParser.ConfigParser()
config.readfp(open(r'sys.argv[0]'))

# Note that the dictionary of loci should be modified by the experimental run Fastq file name.

dict_loci = {"A":  {"pre":["A"],    "exon_num":"8", "primer_start":"2", "primer_end":"4", "paired":config.get('Mapping', 'paired_end'), "filter_cutoff":"10"}, 
             "B":  {"pre":["B"],    "exon_num":"7", "primer_start":"2", "primer_end":"4", "paired":config.get('Mapping', 'paired_end'), "filter_cutoff":"10"}, 
             "C":  {"pre":["C"],    "exon_num":"8", "primer_start":"2", "primer_end":"4", "paired":config.get('Mapping', 'paired_end'), "filter_cutoff":"10"},
             "DPA":{"pre":["DPA1"], "exon_num":"4", "primer_start":"2", "primer_end":"3", "paired":config.get('Mapping', 'paired_end'), "filter_cutoff":"10"},
             "DPB":{"pre":["DPB1"], "exon_num":"5", "primer_start":"2", "primer_end":"3", "paired":config.get('Mapping', 'paired_end'), "filter_cutoff":"10"},
             "DQA":{"pre":["DQA1"], "exon_num":"4", "primer_start":"2", "primer_end":"3", "paired":config.get('Mapping', 'paired_end'), "filter_cutoff":"10"},
             "DQB":{"pre":["DQB1"], "exon_num":"6", "primer_start":"2", "primer_end":"3", "paired":config.get('Mapping', 'paired_end'), "filter_cutoff":"10"},
             "DRB1":{"pre":["DRB1"], "exon_num":"6", "primer_start":"2", "primer_end":"3", "paired":config.get('Mapping', 'paired_end'), "filter_cutoff":"10"}
            }

data_dir   = config.get('EnvironmentPath', 'data_dir')
proj_dir   = config.get('EnvironmentPath', 'proj_dir')
scproj_dir = config.get('EnvironmentPath', 'scproj_dir')
source_dir = config.get('EnvironmentPath', 'mapping')

RUN_ID	   = config.get('Mapping', 'RUN_ID')
equipment  = "@"

if(not os.path.exists(proj_dir+"/qsublog")):
    os.mkdir(proj_dir+"/qsublog")
 
dict_sample = {}
list_loci   = []
list_file1  = []
list_file   = []
sample_loci_pair_set = set()
sample_loci_pair_list = []


os.chdir(data_dir)
file_list = glob.glob("*.fastq")


for file in file_list:
    for loci in dict_loci:
        m = re.search('\w+(?='+loci+'_)' , file)
        if(loci == "A" and m): 
            if(("DP" not in file) and ("DQ" not in file)):
                if(not m.group(0) in dict_sample):
                    dict_sample[m.group(0)] = {}
                if(not loci in dict_sample[m.group(0)]):
                    dict_sample[m.group(0)][loci] = {"first":"", "second":""}
                if("_R1_" in file):
                    dict_sample[m.group(0)][loci]["first"]  = file
                if("_R2_" in file):
                    dict_sample[m.group(0)][loci]["second"] = file
                sample_loci_pair_set.add(m.group(0) + "_" + loci) 
        if(loci == "B" and m):
            if(("DP" not in file) and ("DQ" not in file) and ("DR" not in file)):
                if(not m.group(0) in dict_sample):
                    dict_sample[m.group(0)] = {}
                if(not loci in dict_sample[m.group(0)]):
                    dict_sample[m.group(0)][loci] = {"first":"", "second":""}
                if("_R1_" in file):
                    dict_sample[m.group(0)][loci]["first"]  = file
                if("_R2_" in file):
                    dict_sample[m.group(0)][loci]["second"] = file
                sample_loci_pair_set.add(m.group(0) + "_" + loci)
        if(loci != "A" and loci != "B" and m):
            if(not m.group(0) in dict_sample):
                dict_sample[m.group(0)] = {}
            if(not loci in dict_sample[m.group(0)]):
                dict_sample[m.group(0)][loci] = {"first":"", "second":""}
            if("_R1_" in file):
                    dict_sample[m.group(0)][loci]["first"]  = file
            if("_R2_" in file):
                    dict_sample[m.group(0)][loci]["second"] = file
            sample_loci_pair_set.add(m.group(0) + "_" + loci)

sample_list = dict_sample.keys()
sample_loci_pair_list = list(sample_loci_pair_set)

print("Total sample Loci pairs: ")
print(len(sample_loci_pair_list))
print("\n")

for turn in range(0,192):
    
    job_name = "JOB" + str(turn + 1) + ".TXT"
    run_name = str(turn + 1) + "_" + RUN_ID
    
    JOB_FILE = open(proj_dir + "/" + job_name, "w")
    QSUB_PARALLEL_FILE = open(proj_dir + "/PARALLEL" + str(turn + 1) + ".sh", "w")
    
    # Write in the parallel job submission file:
    
    QSUB_PARALLEL_FILE.write("#!" + config.get('RequiredPackages', 'bash') + "\n\n")
    QSUB_PARALLEL_FILE.write("#PBS -N " + run_name + "\n")
    QSUB_PARALLEL_FILE.write("#PBS -o " + proj_dir + "/" + run_name + ".out\n")
    QSUB_PARALLEL_FILE.write("#PBS -e " + proj_dir + "/" + run_name + ".err\n")
    QSUB_PARALLEL_FILE.write("#PBS -l walltime=200:00:00\n")
    QSUB_PARALLEL_FILE.write("#PBS -q default\n")
    QSUB_PARALLEL_FILE.write("#PBS -l mem=2400MB\n")
    QSUB_PARALLEL_FILE.write("#PBS -m abe\n")
    QSUB_PARALLEL_FILE.write("#PBS -l nodes=1:ppn=2\n\n")
    QSUB_PARALLEL_FILE.write(config.get('RequiredPackages', 'parallel') + " -j 0 -a " + proj_dir + "/" + job_name + "\n\n")
    QSUB_PARALLEL_FILE.write("exit\n")
    
    for ii in range(2*turn, min(len(sample_loci_pair_list),2*(turn+1))):  
        pair_list = re.split("_", sample_loci_pair_list[ii])
        sample = pair_list[0]
        loci   = pair_list[1]
        
        if(not os.path.exists(proj_dir+"/"+sample)):
            os.mkdir(proj_dir+"/"+sample)

        subloci = dict_loci[loci]["pre"][0]
        
        if(not os.path.exists(proj_dir+"/"+sample+"/"+subloci+"_nuc")):
            os.mkdir(proj_dir+"/"+sample+"/"+subloci+"_nuc")
        for gfile in glob.glob(proj_dir+"/"+sample+"/"+subloci+"_nuc/*.sh"):
            os.remove(gfile)
        for gfile in glob.glob(proj_dir+"/"+sample+"/"+subloci+"_nuc/*.py"):
            os.remove(gfile)
        for gfile in glob.glob(source_dir+"/*.sh"):
            shutil.copy(gfile,proj_dir+"/"+sample+"/"+subloci+"_nuc/")
        for gfile in glob.glob(source_dir+"/*.py"):
            shutil.copy(gfile,proj_dir+"/"+sample+"/"+subloci+"_nuc/")
        pbs_list = []
        pbs_list.append("#!" + config.get('RequiredPackages', 'bash'))
        pbs_list.append("map_single=1")
        pbs_list.append("select_single=1")

        pbs_list.append(config.get('RequiredPackages', 'bw2_build'))
        pbs_list.append(config.get('RequiredPackages', 'bw2_align'))
        pbs_list.append(config.get('RequiredPackages', 'samtools'))
        pbs_list.append(config.get('RequiredPackages', 'python'))

        pbs_list.append("fq_first="+data_dir+"/"+dict_sample[sample][loci]["first"])
        pbs_list.append("fq_second="+data_dir+"/"+dict_sample[sample][loci]["second"])
        pbs_list.append("bw2_index_pre="+subloci+"_nuc")
        pbs_list.append("id_file=${bw2_index_dir}/${bw2_index_pre}/${bw2_index_pre}_id.txt")
        pbs_list.append("seq_file=${bw2_index_dir}/${bw2_index_pre}/${bw2_index_pre}_seq.txt")

        pbs_list.append("PROJ=" + proj_dir + "/" + sample + "/" + subloci + "_nuc")
      	pbs_list.append("SCPROJ=" + scproj_dir + "/" + RUN_ID + "/" + sample + "/" + subloci + "_nuc")

        pbs_list.append(config.get('cDNAReferenceDatabase', 'bw2_index_dir'))
        pbs_list.append("exon_bd=" + config.get('cDNAReferenceDatabase', 'ExonBoundary') + "/" + subloci + "_nuc_EX.txt")

        pbs_list.append("filter_cutoff="+dict_loci[loci]["filter_cutoff"])
        pbs_list.append("minimum_mapped_reads=5")
        pbs_list.append("primer_start="+dict_loci[loci]["primer_start"])
        pbs_list.append("primer_end="+dict_loci[loci]["primer_end"])
        pbs_list.append("exon_num="+dict_loci[loci]["exon_num"])
        pbs_list.append("paired="+dict_loci[loci]["paired"])
        pbs_list.append("debug=0")
        pbs_list.append("EQUIPMENT="+equipment)

        with open(source_dir+"/HLA_parts.sh","r") as fi:
            for line in fi:
                pbs_list.append(line.strip("\n"))
        with open(proj_dir+"/"+sample+"/"+subloci+"_nuc/"+"pipeline_HLA.sh", "w") as fo:
            for line in pbs_list:
                fo.write(line+"\n")
        JOB_FILE.write("cd " + proj_dir + "/" + sample + "/" + subloci + "_nuc;")
        JOB_FILE.write("source pipeline_HLA.sh > run.log 2>&1 \n")

    ## close files
    JOB_FILE.close()	
    QSUB_PARALLEL_FILE.close()

QSUB_SUBMISSION_FILE = open(proj_dir + "/QSUB_SUBMIT.sh", "w")

for i in range(0,192):
    
    job_name = "JOB" + str(i + 1) + ".TXT"
    JOB_TXT_FILE = proj_dir + "/" + job_name
    QSUB_FILE = proj_dir + "/PARALLEL" + str(i + 1) + ".sh"
    
    if (os.stat(JOB_TXT_FILE).st_size != 0):
        QSUB_SUBMISSION_FILE.write("qsub " + QSUB_FILE + "\n")
        
QSUB_SUBMISSION_FILE.close()






