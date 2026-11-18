import os
from multiprocessing import Process, Queue
import subprocess
import time
import datetime
import numpy as np
import matplotlib.pyplot as plt

startTime = time.time()
elapTime = time.time()

cores = 24#set number of cores
readBatchSize = 10000000#Limits RAM usage
minmal_rc = 1000#plate threshold

input_dir = '/mnt/appl/Erlantz/Run5_221122/Plate_input/Gastruloids/2_Cell_Annotation_input/'#where the input is located, also make a 'Ouput' folder
output_dir = '/mnt/appl/Erlantz/Run5_221122/Plate_input/Gastruloids/3_Process_input_mouse/'
publicFolder = '/mnt/DevBio_data/Software/'#general location

def getTime(elapTime,elap=True):
	timeString = datetime.date.today().strftime('%d')+'-'+datetime.date.today().strftime('%m')+'-'+datetime.date.today().strftime('%Y')+'_'+datetime.datetime.now().strftime('%H')+'-'+datetime.datetime.now().strftime('%M')+'-'+datetime.datetime.now().strftime('%S')+'.'+datetime.datetime.now().strftime('%f')
	if elap:
		elapsed = time.time()-elapTime
		#global elapTime
		elapTime = time.time()
		return timeString+' '+str(datetime.timedelta(seconds=elapsed)),elapTime
	else:
		return timeString,elapTime

comp = {'T':'A','A':'T','C':'G','G':'C'}
def rev_comp(s):
	rev_comp_s = ''
	for c in s[::-1]:
		rev_comp_s+=comp[c]
	return rev_comp_s

adapters = []
prime_dime = {}
header = True
with open('Plate.txt','r') as f:#Plate sequences
	for line in f:
		if header:
			header = False; continue
		split = line.strip().split('\t')
		adapters.append(split[2])
		prime_dime[split[2]] = rev_comp(split[2])
adapters_set = set(adapters)
adap_len = len(adapters[0])
off_set = 3
print('plate count: ',len(adapters_set))
print('plate len: ',adap_len)
print('Off set: ',off_set,'N'*off_set)

#seq_adapter = ['TGGAATTCTCG','TCGGACTGT','CGAGAATTCCA','ACAGTCCGA']#
seq_adapter = ['CAAGCAGAAGA','TGGAATTCTCG']#

def tag_read(reads,data,start,n,queue):
	#print(str(start+1)+'/'+str(n))
	plate_dic = {}

	total = len(reads)
	write_reads = ''
	i = start
	process_count = 0
	found = 0;	missed = 0;	dime_count = 0;	seq_trim = 0
	while i < total:
		if reads[i]=='':#[0]
			break
		split = reads[i].split('\n')
		seq = split[1]
		plate_umi = seq[:off_set]
		plate_id = seq[off_set:off_set+adap_len]
		if plate_id in adapters_set:
			found+=1
			seq = seq[off_set+adap_len:]
			quality_s =  split[3][off_set+adap_len:]
			dime_pos = seq.find(prime_dime[plate_id])
			if dime_pos != -1:
				dime_count+=1
			else:
				seq_adap_found = -1
				for s_a in seq_adapter:
					seq_adap_pos = seq.find(s_a)
					if seq_adap_pos > seq_adap_found and seq_adap_found == -1:
						seq_adap_found = seq_adap_pos
				if seq_adap_found != -1:
					seq = seq[:seq_adap_found]
					quality_s = quality_s[:seq_adap_found]
					seq_trim+=1
				fastq_read = split[0].split(' ')[0]+'_'+plate_umi+'#'+plate_id+'\n'+seq+'\n'+'+\n'+quality_s+'\n'
				write_reads+=fastq_read				
				if plate_id not in plate_dic:
					plate_dic[plate_id] = [0,0]
				plate_dic[plate_id][0]+=1
				plate_dic[plate_id][1]+=len(quality_s)
		else:
			missed+=1
		process_count+=1
		i+=n
	#print(str(start+1)+' done')
	#print(found,missed,dime_count,seq_trim)
	queue.put([process_count,write_reads,plate_dic,found,missed,dime_count,seq_trim])

def merge_count_dic(main_dic, part_dic):
	for key in part_dic.keys():
		if key not in main_dic:#dictionary.has_key(key):
			main_dic[key] = [0,0]
		main_dic[key][0]+=part_dic[key][0]
		main_dic[key][1]+=part_dic[key][1]
	return main_dic

def Parallel_tagging(reads,data):
	queue = Queue()
	procs = []
	#print('.:start threads:.')
	for core in range(cores):
		p = Process(target=tag_read, args=(reads,data,core,cores,queue))
		procs.append(p)
		p.start()
	main= {}
	write_reads = ''
	process_count = 0
	found = 0;	missed = 0;	dime_count = 0;	seq_trim = 0
	for p in range(cores):
		part = queue.get()
		if part == ():
			exit(1)
		merge_count_dic(main, part[2])
		write_reads+=part[1]
		process_count+=part[0]
		found += part[3];missed += part[4]; dime_count += part[5];	seq_trim += part[6]
	for p in procs:
		p.join()
	print('.:threads done:.')
	return [process_count,write_reads,main,found,missed,dime_count,seq_trim]

def batchTag(f_name):
	# Load reads in batches
	reads =  ['']*readBatchSize#[['','','+\n',''] for x in range(readBatchSize)]
	batchI = 0
	j= 0
	ReadCount = 0
	main= {}
	w_a = 'w'
	process_count =0
	write_name = output_dir+sample_name+'_annotated.fastq'
	found = 0;	missed = 0;	dime_count = 0;	seq_trim = 0
	print(write_name)
	with open(input_dir+f_name, 'r') as f:
		for i,line in enumerate(f):
			if j < readBatchSize*(batchI+1): # in next batch
				reads[j-readBatchSize*batchI]+=line#.strip()[i%4]
				if i%4  == 3:
					ReadCount+=1
			else:
				#Process in parallel
				print('batch: '+str(batchI+1))
				#print(reads[0])
				print(ReadCount,len(reads))
				#plate_dic = tag_read(reads,[],0,1)
				part = Parallel_tagging(reads,[])
				merge_count_dic(main, part[2])
				if batchI > 0:
					w_a = 'a'
				with open(write_name,w_a) as of:
					of.write(part[1])
				process_count+=part[0]
				found += part[3];missed += part[4]; dime_count += part[5];	seq_trim += part[6]
				#Setup next batch
				reads = ['']*readBatchSize#[['','','+\n',''] for x in range(readBatchSize)]
				reads[0]=line#[i%4]
				#ReadCount+=1
				batchI+=1
			if i%4  == 3:
				j+=1
		if reads[0]!='':#[0]
			part = Parallel_tagging(reads,[])
			merge_count_dic(main, part[2])
			with open(write_name,w_a) as of:
				of.write(part[1])
			process_count+=part[0]
			found += part[3];missed += part[4]; dime_count += part[5];	seq_trim += part[6]
	print(process_count,ReadCount,process_count == ReadCount)
	if process_count != ReadCount:
		print('\t*\t Reads in and out incorrect!\t'+str(ReadCount)+'\t'+str(process_count))
	#	exit(1)
	#return total[0:3]
	print('Total\t'+str(found+missed))
	print('found\t'+str(found)+'\t'+str(round((found/(found+missed)*100),2)))
	print('missed\t'+str(missed)+'\t'+str(round((missed/(found+missed)*100),2)))
	print('primer-dimer\t'+str(dime_count)+'\t'+str(round((dime_count/(found+missed)*100),2)))
	print('seq - trim\t'+str(seq_trim)+'\t'+str(round((seq_trim/(found+missed)*100),2)))
	return main

timeString,elapTime = getTime(elapTime);print(timeString)
for f_name in os.listdir(input_dir):
	sample_name = f_name.split('_')[0]+'_'+f_name.split('_')[1]
	print(f_name+'\n'+sample_name)
	
	print('.:Find/Write combi:.') #_annotated.fastq
	plate_dic = batchTag(f_name)
	timeString,elapTime = getTime(elapTime);print(timeString)

	print('.:Find/Write combi:.')
	combi_count = 0
	low_combi = 0
	array = []
	array_avg_len = []
	names = []
	name =f_name.split('.')[0]
	for key in plate_dic.keys():
		names.append(key)
		array.append(plate_dic[key][0])
		array_avg_len.append(round(float(plate_dic[key][1])/plate_dic[key][0],2))

	with open(output_dir+sample_name+'_cell_stats.tsv','w') as of_all:
		of_all.write('plate\treadcount\tavg_len\n')
		for i in np.array(array).argsort()[::-1]:
			of_all.write(names[i]+'\t'+str(array[i])+'\t'+str(array_avg_len[i])+'\n')
	
	plt.title(name)
	plot_array = sorted(set(array))[::-1]
	plt.plot(plot_array,'b')
	plt.plot([minmal_rc]*len(plot_array),'r')
	plt.xlabel('unique count')
	plt.ylabel('read count')
	#plt.plot(sorted(array)[::-1],'b')
	#plt.show()
	plt.savefig(output_dir+sample_name+'.png',dpi=300)
	plt.clf()

	timeString,elapTime = getTime(elapTime);print(timeString)

	
timeString = datetime.date.today().strftime('%d')+'-'+datetime.date.today().strftime('%m')+'-'+datetime.date.today().strftime('%Y')+'_'+datetime.datetime.now().strftime('%H')+':'+datetime.datetime.now().strftime('%M')+':'+datetime.datetime.now().strftime('%S')
print(' End: '+timeString)
timeSpend = time.time() - startTime
timeElapsed = 'Elapsed Time: ' + str(int(timeSpend/60/60)) + ':' + str(int(timeSpend/60) - int(timeSpend/60/60)*60) + ':' + str(int(timeSpend) - int(timeSpend/60)*60)
print(timeElapsed)
