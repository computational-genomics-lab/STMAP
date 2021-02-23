#! /usr/bin/env python3
import os
import re
import sys
import shutil
from collections import OrderedDict
import optparse
from sys import exit
import psutil
import subprocess
import pandas as pd

totalcpus = psutil.cpu_count()
threads = (totalcpus-1)
mem = psutil.virtual_memory()
maxMemory= int((mem.available/1073741824) -1)


	
required="dataDir  emailAddress".split()

parser = optparse.OptionParser()
parser.add_option('-i', '--inputDir',
				  help="Path to Input Data Directory, Default 'raw_data'",
				  dest='dataDir',
				  default=(os.path.abspath(os.path.join((os.getcwd()),"raw_data")))
				  )

parser.add_option('-p', '--projectName',
				  help="Name of the project directory, Default 'metagenome_project'",
				  dest='projectName',
				  default=(os.path.abspath(os.path.join((os.getcwd()),"metagenome_project")))
				  )

parser.add_option('-d', '--domain',
				  type='choice',
				  choices=['prokaryote', 'eukaryote'],
				  help="Domain of the Organism. Choose from 'prokaryote' or 'eukaryote' ",
		  		  default='prokaryote',
				  dest='domain'
				  )

parser.add_option('-o', '--symLinkDir',
				  help="[Optional] Path to symbolic link Data Directory, Default 'raw_data_symlink'",
				  dest='symLinkDir',
				  default=(os.path.abspath(os.path.join((os.getcwd()),"raw_data_symlink")))
				  )

parser.add_option('-s', '--schedulerPort',
				  help="[Optional Parameter] Scheduler Port Number. default =int[8888] ",
				  type="int",
				  default=8082
				 )

parser.add_option('-e', '--emailAddress',
				  help="Provide your email address =email[abc@xyz.com]",
				  dest='emailAddress'				  
				 )

parser.add_option('-t', '--threads',
				  help="[Optional Parameter, ] Number of threads. Default = (total threads -1)",
				  dest='threads',
				  type="int",
				  default = threads)

parser.add_option('-x', '--maxMemory',
				  help="[Optional Parameter] Maximum Memory in GB. Default = (available memory in GB -1)",
				  dest='maxMemory',
				  type="int",
				  default=maxMemory
				 )

options,args = parser.parse_args()

for r in required:
	if options.__dict__[r] is None:
		parser.error("parameter %s required" %r)


option_dict = vars(options)

dataDir = option_dict.get('dataDir')
symLinkDir = option_dict.get('symLinkDir')
projectDir=option_dict.get('projectName')
projectName=os.path.basename(os.path.dirname(projectDir))
domain=option_dict.get('domain')
adapter=os.path.join(os.getcwd(),"tasks","utility",'adapters.fasta.gz')
email=option_dict.get('emailAddress')
port = int(option_dict.get('schedulerPort'))
cpu = int(option_dict.get('threads'))
memory = int(option_dict.get('maxMemory'))


def createFolder(directory):
	try:
		if not os.path.exists(directory):
			os.makedirs(directory)
	except OSError:
		print ('Error: Creating directory. ' + directory)


def run_cmd(cmd):
	p = subprocess.Popen(cmd, bufsize=-1,
				 shell=True,
				 universal_newlines=True,
				 stdout=subprocess.PIPE,
				 executable='/bin/bash')
	output = p.communicate()[0]
	return output



current_folder=os.path.join(os.getcwd())
symLinkDir=os.path.abspath(os.path.join(symLinkDir))
adapter=os.path.abspath(os.path.join(os.getcwd(),"tasks","utility",'adapters.fasta.gz'))
paired_end_read_dir=os.path.abspath(os.path.join(symLinkDir, "pe" ))
pac_read_dir=os.path.abspath(os.path.join(symLinkDir, "pac" ))
ont_read_dir=os.path.abspath(os.path.join(symLinkDir, "ont" ))
projectDir=os.path.abspath(os.path.join(projectDir))
inputDir=os.path.abspath(os.path.join(dataDir))
excludeDir=os.path.abspath(os.path.join(symLinkDir, "ex" ))


createFolder(excludeDir)
createFolder(paired_end_read_dir)
createFolder(pac_read_dir)
createFolder(ont_read_dir)
createFolder(projectDir)
createFolder(symLinkDir)
createFolder("config")


if os.path.isdir(inputDir):
	files = [f for f in os.listdir(inputDir) if os.path.isfile(os.path.join(inputDir, f))]
	keys = []
	fileList = re.compile(r'^(.+?).(fastq|fq|fastq\.gz|fq\.gz)?$')
	for file in files:
		if fileList.search(file):
			keys.append(file)

	dicts = OrderedDict ()
		#keys = [f for f in os.listdir(".") if os.path.isfile(os.path.join(".", f))]
	for i in keys:
		accepted_values="pe ont pac ex".split()
		val = str(input("Enter Data Type of {data}: \tchoose from [pe:paired-end, ont:nanopore, pac:pacbio ex:exclude]: ".format(data=i)))
		if val in accepted_values:
			dicts[i] = val
		else:
				print(f'{val} is not a valid option. \tchoose from [pe, ont, pac, ex]: ')
				val = str(input("Enter Data Type of {data}: \tchoose from [pe,ont, pac, ex]: ".format(data=i)))

	for key, val in dicts.items():
		if not os.path.exists(os.path.join(symLinkDir, val)):
			os.mkdir(os.path.join(symLinkDir, val))

		##ln -nsf method
	for key, val in dicts.items():
		dest = (os.path.join(symLinkDir,val,key))
		src = (os.path.join(inputDir,key))
		source = os.path.abspath(src)
		destination = os.path.abspath(dest)
		escape="\'"
		

		link_cmd = "ln -nsf "
		create_symlink = "{link_cmd} {source} {destination}".format(link_cmd=link_cmd,source=source,destination=destination)
		print("****** NOW RUNNING COMMAND ******: " + create_symlink)
		print (run_cmd(create_symlink))

		###########################################
def paired_end_samples(pe_dir):
	pe_read_list=os.listdir(paired_end_read_dir)
	sample_list=[]
	for read in pe_read_list:
		pe_allowed_extn=["_R1.fq","_R1.fastq","_R1.fq.gz","_R1.fastq.gz"]
		if any (read.endswith(ext) for ext in pe_allowed_extn):
			sample_name=read.split("_R1.f",1)[0]
			sample_list.append(sample_name)
			file_extn=read.split('.',1)[1]

	with open ((os.path.join(os.getcwd(),"config",'pe_samples.lst')),'w') as file:
		for sample in sample_list:
			file.write("%s\n" % sample)
	file.close()
	return file_extn




def paired_end_group(pe_dir):
	if ((len(os.listdir(paired_end_read_dir))!=0)):
		pe_read_list=os.listdir(paired_end_read_dir)
		sample_list=[]
		for read in pe_read_list:
			pe_allowed_extn=["_R1.fq","_R1.fastq","_R1.fq.gz","_R1.fastq.gz"]
			if any (read.endswith(ext) for ext in pe_allowed_extn):
				sample_name=read.split("_R1.f",1)[0]
				sample_list.append(sample_name)

		print ('''

			For Genome resolved metagenomics you may provide a single condition associated with each sample
		
			Example:
		
			samples         conditions
			sample_s1_rep1  control
			sample_s1_rep2  control
			sample_s2_rep1  control
			sample_s2_rep2  control

			Note: (1) Condition name is case sensitive.


			For Enrichment analysis you have to enter the biological conditions associated with each sample
		
			Example:
		
			samples         conditions
			sample_s1_rep1  control
			sample_s1_rep2  control
			sample_s2_rep1  treated
			sample_s2_rep2  treated

			Note: (1) For Enrichment Analysis, user must provide two conditions. condition name is case sensitive.
			  	  (2) Each sample must have atleast one replicate


			 ''')
						
		dictionary=OrderedDict()
		for sample in range(len(sample_list)):
			condition=input('Enter Condition for %s: ' %(sample_list[sample]))
			dictionary.update({sample_list[sample]:condition})
			#print (dictionary)

			df = pd.DataFrame()
			df['lable'] = dictionary.keys()
			df['samples'] = dictionary.keys()
			df['conditions'] = dictionary.values()
			sampledf = df[["samples"]]
			groupdf = df[["samples","conditions"]]


		if (len(df.conditions.unique()) >2):
			print ('\n\n.................NOTE.................... \n')
			sys.exit("You have provided more than two conditions. please check the values you entered")

		#elif ((df.groupby('conditions').size().min()) == 1):
		if (len(df.conditions.unique()) ==1):
			sampledf.to_csv("config" + "/" + "samples.txt", index=False, header=False)
			groupdf.to_csv("config" + "/" + "metagenome_group.tsv", sep='\t', index=False, header=False)
			print ('\n\n.................NOTE.................... \nYou have provided single condition for each sample. Enrichment analysis Can not be performed\n')
			print('Config folder generated')
			
		else:
			sampledf.to_csv("config" + "/" + "samples.txt", index=False, header=False)
			groupdf.to_csv("config" + "/" + "metagenome_group.tsv", sep='\t', index=False, header=False)
			groupdf.to_csv("config" + "/" + "metagenome_condition.tsv", sep='\t', index=False)
			print('Config folder generated')
		
#################################################################################

def pacbio_samples(pb_dir):
	raw_pb_read_list=os.listdir(pac_read_dir)
	sample_list=[]
	for read in raw_pb_read_list:
		raw_pb_allowed_extn=[".fq",".fastq",".fq.gz",".fastq.gz"]
				
		if any (read.endswith(ext) for ext in raw_pb_allowed_extn):
					
			sample_name=read.split(".",1)[0]
			sample_list.append(sample_name)

			file_extn=read.split('.',1)[1]


	with open ((os.path.join(os.getcwd(),"config",'pac_samples.lst')),'w') as file:
		for sample in sample_list:
			file.write("%s\n" % sample)
	file.close()

	return file_extn

################################################################################
def ont_samples(ont_raw_dir):
	corr_ont_read_list=os.listdir(ont_read_dir)
	sample_list=[]
	for read in corr_ont_read_list:
		corr_ont_allowed_extn=[".fq",".fastq",".fq.gz",".fastq.gz"]
				
		if any (read.endswith(ext) for ext in corr_ont_allowed_extn):
					
			sample_name=read.split(".",1)[0]
			sample_list.append(sample_name)
			file_extn=read.split('.',1)[1]


	with open ((os.path.join(os.getcwd(),"config",'ont_samples.lst')),'w') as file:
		for sample in sample_list:
			file.write("%s\n" % sample)
	file.close()

	return file_extn

#Get Read Extension


if ((len(os.listdir(paired_end_read_dir))!=0)):
	paired_end_read_suffix=paired_end_samples(paired_end_read_dir)
	paired_end_group(paired_end_read_dir)


if ((len(os.listdir(ont_read_dir))!=0)):
	ont_read_suffix=ont_samples(ont_read_dir)

if ((len(os.listdir(pac_read_dir))!=0)):
	pac_read_suffix=pacbio_samples(pac_read_dir)





with open('luigi.cfg', 'w') as config:
	config.write('[core]\n')
	config.write('default-scheduler-port:{port}\n'.format(port=port))
	config.write('error-email={email}\n\n'.format(email=email))
	config.write('[GlobalParameter]\n')
	config.write('projectName={projectName}\n'.format(projectName=projectName))
	config.write('projectDir={projectDir}/\n'.format(projectDir=projectDir))
	config.write('domain={domain}\n'.format(domain=domain))
	config.write('adapter={adapter}\n'.format(adapter=adapter))
	
	#PE READ
	if ((len(os.listdir(paired_end_read_dir))!=0) and (len(os.listdir(ont_read_dir))==0) and (len(os.listdir(pac_read_dir))==0)):
		paired_end_read_suffix=paired_end_samples(paired_end_read_dir)
		config.write('pe_read_dir={paired_end_read_dir}/\n'.format(paired_end_read_dir=paired_end_read_dir))
		config.write('pe_read_suffix={paired_end_read_suffix}\n'.format(paired_end_read_suffix=paired_end_read_suffix))
		config.write('seq_platforms=pe\n')
		config.write('pac_read_dir=NA\n')		
		config.write('pac_read_suffix=NA\n')
		config.write('ont_read_dir=NA\n')
		config.write('ont_read_suffix=NA\n')


	if ((len(os.listdir(ont_read_dir))!=0)  and (len(os.listdir(paired_end_read_dir))==0) and (len(os.listdir(pac_read_dir))==0)):
		ont_read_suffix=ont_samples(ont_read_dir)
		config.write('ont_read_dir={ont_read_dir}/\n'.format(ont_read_dir=ont_read_dir))
		config.write('ont_read_suffix={ont_read_suffix}\n'.format(ont_read_suffix=ont_read_suffix))
		config.write('seq_platforms=ont\n')
		config.write('pac_read_dir=NA\n')		
		config.write('pac_read_suffix=NA\n')
		config.write('pe_read_dir=NA\n')
		config.write('pe_read_suffix=NA\n')



	if ((len(os.listdir(pac_read_dir))!=0) and (len(os.listdir(ont_read_dir))==0) and (len(os.listdir(paired_end_read_dir))==0)):
		pac_read_suffix=pacbio_samples(pac_read_dir)

		config.write('pac_read_dir={pac_read_dir}/\n'.format(pac_read_dir=pac_read_dir))
		config.write('pac_read_suffix={pac_read_suffix}\n'.format(pac_read_suffix=pac_read_suffix))
		config.write('seq_platforms=pac\n')

		config.write('pe_read_dir=NA\n')		
		config.write('pe_read_suffix=NA\n')
		config.write('ont_read_dir=NA\n')
		config.write('ont_read_suffix=NA\n')


	if ((len(os.listdir(paired_end_read_dir))!=0) and (len(os.listdir(pac_read_dir))!=0) and (len(os.listdir(ont_read_dir))==0)):
		paired_end_read_suffix=paired_end_samples(paired_end_read_dir)
		pac_read_suffix=pacbio_samples(pac_read_dir)

		config.write('pe_read_dir={paired_end_read_dir}/\n'.format(paired_end_read_dir=paired_end_read_dir))
		config.write('pac_read_dir={pacbio_read_dir}/\n'.format(pacbio_read_dir=pacbio_read_dir))
		config.write('pe_read_suffix={paired_end_read_suffix}\n'.format(paired_end_read_suffix=paired_end_read_suffix))
		config.write('pac_read_suffix={pac_read_suffix}\n'.format(pac_read_suffix=pac_read_suffix))
		config.write('seq_platforms=pe-pac\n')

		config.write('ont_read_dir=NA\n')
		config.write('ont_read_suffix=NA\n')


	if ((len(os.listdir(paired_end_read_dir))!=0) and (len(os.listdir(ont_read_dir))!=0) and (len(os.listdir(pac_read_dir))==0)):
		paired_end_read_suffix=paired_end_samples(paired_end_read_dir)
		ont_read_suffix=ont_samples(ont_read_dir)

		config.write('pe_read_dir={paired_end_read_dir}/\n'.format(paired_end_read_dir=paired_end_read_dir))
		config.write('ont_read_dir={ont_read_dir}/\n'.format(ont_read_dir=ont_read_dir))
		config.write('pe_read_suffix={paired_end_read_suffix}\n'.format(paired_end_read_suffix=paired_end_read_suffix))
		config.write('ont_read_suffix={ont_read_suffix}\n'.format(ont_read_suffix=ont_read_suffix))
		config.write('seq_platforms=pe-ont\n')

		config.write('pac_read_dir=NA\n')		
		config.write('pac_read_suffix=NA\n')


	config.write('threads={cpu}\n'.format(cpu=cpu)) 
	config.write('maxMemory={memory}\n'.format(memory=memory))
	config.close()

	print("the luigi config file generated")
