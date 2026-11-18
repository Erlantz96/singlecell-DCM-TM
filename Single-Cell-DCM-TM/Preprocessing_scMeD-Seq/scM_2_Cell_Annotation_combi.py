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
minmal_rc = 1000#combi threshold

#input_dir = '/mnt/appl/Erlantz/scNMT/Pool_37/Chase_combi_enterocytes/2_Cell_Annotation_input_mouse/'#where the input is located, also make a 'Ouput' folder
#output_dir = '/mnt/appl/Erlantz/scNMT/Pool_37/Chase_combi_enterocytes/3_Process_input_mouse/'
#input_dir = '/mnt/appl/Erlantz/Run5_221122/Combi_input/2_Cell_Annotation_input/'#where the input is located, also make a 'Ouput' folder
#output_dir = '/mnt/appl/Erlantz/Run5_221122/Combi_input/3_Process_input/'
input_dir = '/mnt/appl/Erlantz/Run8_scMed_SI_NPC_Kidney/2_Cell_Annotation_input/'#where the input is located, also make a 'Ouput' folder
output_dir = '/mnt/appl/Erlantz/Run8_scMed_SI_NPC_Kidney/3_Process_input/'
publicFolder = '/mnt/DevBio_data/Software/'#general location

tbc = []
even = []
odd = []
header = True
with open('Combi.txt','r') as f:#Combi sequences
	for line in f:
		if header:
			header = False
			continue
		split = line.strip().split('\t')
		tbc.append(split[0])
		even.append(split[1])
		odd.append(split[2])

even_start = 'TGACTTG'
odd_start =  'GACAACT'
tbc_len = len(tbc[0])
even_len = len(even[0])
odd_len = len(odd[0])
print('combi len: ',tbc_len,even_len,odd_len)

def getTime(elapTime,elap=True):
	timeString = datetime.date.today().strftime('%d')+'-'+datetime.date.today().strftime('%m')+'-'+datetime.date.today().strftime('%Y')+'_'+datetime.datetime.now().strftime('%H')+'-'+datetime.datetime.now().strftime('%M')+'-'+datetime.datetime.now().strftime('%S')+'.'+datetime.datetime.now().strftime('%f')
	if elap:
		elapsed = time.time()-elapTime
		#global elapTime
		elapTime = time.time()
		return timeString+' '+str(datetime.timedelta(seconds=elapsed)),elapTime
	else:
		return timeString,elapTime

def tag_read(reads,data,start,n,queue):
	#print(str(start+1)+'/'+str(n))
	combi_dic = {}

	total = len(reads)
	write_reads = ''
	i = start
	process_count = 0
	while i < total:
		if reads[i]=='':#[0]
			break
		split = reads[i].split('\n')
		seq = split[1]
		found_even_s = seq.find(even_start)
		found_odd_s = seq.find(odd_start)
		tbc_s = seq[1:found_even_s]
		even_s = seq[found_even_s:found_odd_s]
		even_ID = seq[found_even_s+7:found_odd_s]
		odd_s = seq[found_odd_s:found_odd_s+odd_len]
		odd_ID = seq[found_odd_s+7:found_odd_s+odd_len]
		seq_start = found_odd_s+odd_len
		genome_s =  seq[seq_start:]

		#print(tbc_s+'\t'+str(tbc_s in tbc))
		#print(even_s+'\t'+str(even_s in even))
		#print('.'*7+even_ID)
		#print(odd_s+'\t'+str(odd_s in odd))
		#print('.'*7+odd_ID)
		combi_id = tbc_s+'_'+even_ID+'_'+odd_ID
		#print(len(tbc_s),found_even_s)
		if tbc_s in tbc and even_s in even and odd_s in odd:
			quality_s =  split[3][seq_start:]+'\n'
			fastq_read = split[0].split(' ')[0]+'#'+combi_id+'\n'+genome_s+'\n'+'+\n'+quality_s
			write_reads+=fastq_read
			if combi_id not in combi_dic:
				combi_dic[combi_id] = [0,0]
			combi_dic[combi_id][0]+=1
			combi_dic[combi_id][1]+=len(quality_s)
		process_count+=1
		i+=n
	#print(str(start+1)+' done')
	queue.put([process_count,write_reads,combi_dic])

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
	for p in range(cores):
		part = queue.get()
		if part == ():
			exit(1)
		merge_count_dic(main, part[2])
		write_reads+=part[1]
		process_count+=part[0]
	for p in procs:
		p.join()
	print('.:threads done:.')
	return [process_count,write_reads,main]

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
				#combi_dic = tag_read(reads,[],0,1)
				part = Parallel_tagging(reads,[])
				merge_count_dic(main, part[2])
				if batchI > 0:
					w_a = 'a'
				with open(write_name,w_a) as of:
					of.write(part[1])
				process_count+=part[0]
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
	print(process_count,ReadCount,process_count == ReadCount)
	if process_count != ReadCount:
		print('\t*\t Reads in and out incorrect!\t'+str(ReadCount)+'\t'+str(process_count))
	#	exit(1)
	#return total[0:3]
	return main

timeString,elapTime = getTime(elapTime);print(timeString)
for f_name in os.listdir(input_dir):
	sample_name = f_name.split('_')[0]+'_'+f_name.split('_')[1]
	print(f_name+'\n'+sample_name)
	
	print('.:Find/Write combi:.') #_annotated.fastq
	combi_dic = batchTag(f_name)
	timeString,elapTime = getTime(elapTime);print(timeString)

	print('.:Find/Write combi:.')
	combi_count = 0
	low_combi = 0
	array = []
	array_avg_len = []
	names = []
	name =f_name.split('.')[0]
	for key in combi_dic.keys():
		names.append(key)
		array.append(combi_dic[key][0])
		array_avg_len.append(round(float(combi_dic[key][1])/combi_dic[key][0],2))

	with open(output_dir+sample_name+'_cell_stats.tsv','w') as of_all:
		of_all.write('combi\treadcount\tavg_len\n')
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
