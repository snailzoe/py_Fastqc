#!/usr/bin/env python3
import os
import sys
import re
import argparse
import math
from itertools import product
from collections import Counter
import numpy as np
from scipy.stats import norm
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.mlab as mlab
import matplotlib.collections as collections
import matplotlib.ticker as ticker
import seaborn as sns

__author__='snail~'
__mail__= 'snail1993@outlook.com'
__doc__='base stat'

def check_outdir(outdir):
	file_outdir = os.path.abspath(outdir)
	if not os.path.exists(file_outdir):
		os.makedirs(file_outdir)
	return(file_outdir)

def base_stat(seq):
	input_file = open (seq, 'r')
	length, qual = [], {}
	line_num, GC_num, base_num, read_num = 0, 0, 0, 0
	for line in input_file:
		line_num += 1
		if line_num %4 == 2:
			#统计reads数
			read_num += 1
			GC_num += line.count("C") + line.count("G")
			base_num += len(line.strip())
			length.append(len(line.strip()))
		#统计碱基质量，字典嵌套列表
		if line_num %4 == 0:
			qual[read_num] = []
			for q in line.strip():
				qual[read_num].append(ord(q)-33)
	#统计GC含量
	GC_percent = math.floor(GC_num/base_num*100)
	read_length = np.unique(length)[0]
#	print(GC_percent)
	return(read_num, GC_percent, read_length, qual, base_num)

def show_stat (filename, read_num, GC_percent, read_length, output):
	outfile1 = open(os.path.join(output, "base.txt"), 'w')
	col_name = ["Total Sequences", "File type", "Encoding", "Total Sequences", "Sequences flagged as poor quality", "Sequence length", "%GC"]
	value = [os.path.basename(filename), "Conventional base calls", "Sanger / Illumina 1.9", read_num, 0, read_length, GC_percent]
	#左对齐或居中显示
	print("\033[0;37;44mMeasure\033[0m".center(35), "\033[0;37;44mValue\033[0m".center(55), file = outfile1)
	for i in range(len(col_name)):
		print((col_name[i]).ljust(35), value[i], file = outfile1)

def base_qual_show(read_length, qual, read_num, output):
	#按碱基位置统计质量
	##字典形式
	base_qual = {}
	for base in range(1,read_length+1):
		base_qual[base] = []
		for r in range(1,read_num+1):
			base_qual[base].append(qual[r][base-1])
	#统计每个位点质量均值
	base_qual_mean = []
	for key in base_qual.keys():
		base_qual_mean.append(np.mean(base_qual[key]))
	##整理数据，列：base_pos; quality，行：碱基数
	boxplot_data = pd.DataFrame(base_qual)
	plt.figure(figsize = (20, 16), facecolor = "white")
	#sym=''不显示异常点，whiskerprops须颜色，boxprops箱体线以及箱内部颜色，medianprops众数线，capprops极值线，whis极值范围
	boxplot_data.boxplot(patch_artist=True,
						whiskerprops = {'color': "black"}, boxprops = {'color':'black', 'facecolor':'yellow'},medianprops = {'color':'red'}, capprops = {'color': "black"}, whis = [10,90], 
						sym='')
	#添加均值线
	plt.plot(base_qual.keys(),base_qual_mean)
	#纵轴分隔单位2
	plt.gca().yaxis.set_major_locator(ticker.MultipleLocator(2))
	plt.ylim(0,41)
	#设置背景颜色，axhspan横向，axvspan纵向
	plt.axhspan(0, 20, facecolor = ('#fdc1c5'))
	plt.axhspan(20, 28, facecolor = ('#fbeeac'))
	plt.axhspan(28, 41, facecolor = ('#bcf5a6'))
	for i in range(1, read_length, 2):
		plt.axvspan(i-0.5, i+0.5, facecolor = ('whitesmoke'), alpha=0.3)
	plt.xlabel("Position in read (bp)", fontsize = 15)
	plt.title("Quality scores across all bases (Sanger / Illumina 1.9 encoding)", fontsize = 15)
	#网格线
	plt.grid(b=None)
	outfile2 = os.path.join(output, "base_qual.png")
	plt.savefig(outfile2)
#	plt.show()
	plt.close()

def reads_qual_show(qual, output):
	reads_qual_mean = []
	#计算每条reads的质量值
	for key in qual.keys():
		reads_qual_mean.append(math.floor(np.mean(qual[key])))
	#统计质量值的reads数
	reads_qual_show = {}
	for s in reads_qual_mean :
		if s in reads_qual_show:
			reads_qual_show[s] += 1
		else:
			reads_qual_show[s] = 1
	#字典按照键排序
	plt.figure(figsize = (20, 16), facecolor = "white")
	plt.plot(dict(sorted(reads_qual_show.items())).keys(),dict(sorted(reads_qual_show.items())).values(),color='red')
	plt.gca().set_xlim([8,40])
	plt.gca().yaxis.set_major_locator(ticker.MultipleLocator(50))
	plt.gca().xaxis.set_major_locator(ticker.MultipleLocator(1))
	plt.gca().set_ylim([0,320])
	plt.xlabel("Mean Sequence Quality (Phred Score)", fontsize = 15)
	for i in range(1, 40, 2):
		plt.axvspan(i-0.5, i+0.5, facecolor = ('lightgrey'), alpha=0.5)
	#加注释框
	plt.text(35, 310, "Average Quality per read", fontsize = 18, color = "red", bbox = {'facecolor': 'white'})
	plt.title("Quality score distribution over all sequences", fontsize = 15)
	#去边框
	plt.gca().spines['top'].set_visible(False)
	plt.gca().spines['right'].set_visible(False)
#	plt.show()
	outfile3 = os.path.join(output, "reads_qual.png")
	plt.savefig(outfile3)
	plt.close()

def base_GC(seq, read_length):
	input_file = open (seq, 'r')
	reads, GC, GC_content =  [], {}, []
	line_num = 0
	for line in input_file:
		line_num += 1
		if line_num %4 == 2:
			#统计每条reads中的GC含量，返还为GC列表
			GC_content.append(math.floor((line.count("C") + line.count("G")) / read_length * 100))
			reads.append(list(line.strip()))
	#reads为长度为read_num的列表，每个子列表为每条reads，统计每个位置上所有reads的ATCG数量
	for p in range(0, read_length):
		#调取列表套列表中每个子列表中的第一位；
		GC[p+1] = []
		GC[p+1].append(list(list(zip(*reads))[p]).count("A"))
		GC[p+1].append(list(list(zip(*reads))[p]).count("T"))
		GC[p+1].append(list(list(zip(*reads))[p]).count("C"))
		GC[p+1].append(list(list(zip(*reads))[p]).count("G"))
		#统计每条reads中N的含量，返还为N_count列表
		GC[p+1].append(list(list(zip(*reads))[p]).count("N"))
	return(GC, GC_content, reads)

def Dup(seq, read_num):
	input_file = open (seq, 'r')
	line_num = 0
	Readsline, Dup, Dup_show = [], {}, {}
	for line in input_file:
		line_num += 1
		if line_num %4 == 2:
			Readsline.append(line.strip())
	#Dup{seq:count}
	for r in Readsline:
		Dup[r] = Readsline.count(r)
	#temp{Dup level:reads num}
	temp = Counter(list(Dup.values()))
	dupnum = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 50, 100, 500, 1000, 5000, 10000]
	for i in temp.keys():
		if i in dupnum:
			Dup_show[i] = math.floor(temp[i] / read_num *100)
		elif i >= 10 and i < 50:
			Dup_show[10] = math.floor(temp[i] / read_num *100)
		elif i >= 50 and i < 100:
			Dup_show[50] = math.floor(temp[i] / read_num *100)
		elif i >= 100 and i < 500:
			Dup_show[100] = math.floor(temp[i] / read_num *100)
		elif i >= 500 and i < 1000:
			Dup_show[500] = math.floor(temp[i] / read_num *100)
		elif i >= 1000 and i < 5000:
			Dup_show[1000] = math.floor(temp[i] / read_num *100)
		elif i >= 5000 and i < 10000:
			Dup_show[5000] = math.floor(temp[i] / read_num *100)
		elif i >= 10000:
			Dup_show[10000] = math.floor(temp[i] / read_num *100)
	for c in dupnum:
		if c not in Dup_show.keys():
			Dup_show[c] = 0
	return(Dup_show, Readsline)

def Adapter(Readsline, read_length):
#	SOLID = "CAGCGAGCCCCAGCAGCCA"
#	Universal = "CAGCGAGCCCCAGCAGCCA"
	Universal = "AGATCGGAAGAG"
	small3 = "TGGAATTCTCGG"
	small5 = "GATCGTCGGACT"
	Transposase = "CTGTCTCTTATA"
	SOLID = "CGCCTTGGCCGT"
	A1, A2, A3, A4, A5 = {}, {}, {}, {}, {}
	a1, a2, a3, a4, a5 = [], [], [], [], []
	for R in Readsline:
		#不存在返回-1，存在返回比对上的第一个字符的位置
		a1.append(R.find(Universal))
		a2.append(R.find(small3))
		a3.append(R.find(small5))
		a4.append(R.find(Transposase))
		a5.append(R.find(SOLID))
	#统计每个位置上出现adapter的频数
	for a in range(1, read_length-10):
		A1[a] = (np.array(a1)+1).tolist().count(a)
		A2[a] = (np.array(a2)+1).tolist().count(a)
		A3[a] = (np.array(a3)+1).tolist().count(a)
		A4[a] = (np.array(a4)+1).tolist().count(a)
		A5[a] = (np.array(a5)+1).tolist().count(a)
#	print(A1, A2, A3, A4, A5)
	return(A1, A2, A3, A4, A5)

def Adapter_show(A1, A2, A3, A4, A5, read_length, output):
	plt.figure(figsize = (20, 16), facecolor = "white")
	plt.plot(A1.keys(), A1.values(), color = "red", label = "Illumina Universal Adapter")
	plt.plot(A2.keys(), A2.values(), color = "blue", label = "Illumina Small RNA 3' Adapter")
	plt.plot(A3.keys(), A3.values(), color = "green", label = "Illumina Small RNA 5' Adapter")
	plt.plot(A4.keys(), A4.values(), color = "black", label = "Nextera Transposase Sequence")
	plt.plot(A5.keys(), A5.values(), color = "deeppink", label = "SOLID Small RNA Adapter")
	plt.gca().set_ylim([0, 100])
	plt.gca().yaxis.set_major_locator(ticker.MultipleLocator(10))
	plt.gca().xaxis.set_major_locator(ticker.MultipleLocator(1))
	for i in range(2, read_length-10, 2):
		plt.axvspan(i-0.5, i+0.5, facecolor = ('lightgrey'), alpha=0.5)
	plt.gca().spines['top'].set_visible(False)
	plt.gca().spines['right'].set_visible(False)
	plt.grid(axis = "y")
	plt.title("% Adapter", fontsize = 15)
	plt.xlabel("Position in read (bp)", fontsize = 15)
	plt.legend()
#	plt.show()
	outfile9 = os.path.join(output, "Adapter.png")
	plt.savefig(outfile9)
	plt.close()

def Dup_level_show(Dup_show, read_num, output):
	per = Dup_show[1]
	plt.figure(figsize = (20, 16), facecolor = "white")
	plt.plot(list(range(1,17,1)), list(Dup_show.values()), color = "red", label = "% Deduplicated sequences")
	plt.gca().set_ylim([0, 100])
	plt.gca().yaxis.set_major_locator(ticker.MultipleLocator(10))
	plt.gca().xaxis.set_major_locator(ticker.MultipleLocator(1))
	x_ticks_labels = ["1", "2", "3", "4", "5", "6", "7", "8", "9", ">10", ">50", ">100", ">500", ">1k", ">5k", ">10k"]
	plt.xticks(range(1,17,1), x_ticks_labels)
	for i in range(2, 17, 2):
		plt.axvspan(i-0.5, i+0.5, facecolor = ('lightgrey'), alpha=0.5)
	plt.gca().spines['top'].set_visible(False)
	plt.gca().spines['right'].set_visible(False)
	plt.grid(axis = "y")
	plt.title(f"Percent of seqs remaining if deduplicated {per}%", fontsize = 15)
	plt.xlabel("Sequence Duplication Level", fontsize = 15)
	plt.legend()
#	plt.show()
	outfile8 = os.path.join(output, "Dup_Level.png")
	plt.savefig(outfile8)
	plt.close()

def base_show(GC, read_num, read_length, output):
	A, T, C, G = [], [], [], []
	for key, value in GC.items():
		A.append(value[0]/read_num*100)
		T.append(value[1]/read_num*100)
		C.append(value[2]/read_num*100)
		G.append(value[3]/read_num*100)
#	print(A)
	plt.figure(figsize = (20, 16), facecolor = "white")
	plt.plot(GC.keys(), A, color='green', label = "%A")
	plt.plot(GC.keys(), T, color='red', label = "%T")
	plt.plot(GC.keys(), C, color='blue', label = "%C")
	plt.plot(GC.keys(), G, color='black', label = "%G")
	plt.gca().set_ylim([0, 100])
	plt.gca().set_xlim([0, read_length])
	plt.gca().xaxis.set_major_locator(ticker.MultipleLocator(1))
	plt.gca().yaxis.set_major_locator(ticker.MultipleLocator(10))
	plt.title("Sequence content across all bases", fontsize = 15)
	plt.xlabel("Position in read (bp)", fontsize = 15)
	for i in range(1, read_length+1, 2):
		plt.axvspan(i, i+1, facecolor = ('lightgrey'), alpha=0.5)
	#去边框
	plt.gca().spines['top'].set_visible(False)
	plt.gca().spines['right'].set_visible(False)
	plt.grid(axis = "y")
	#http://www.matplotlib.org.cn/gallery/text_labels_and_annotations/rainbow_text.html 需要尝试
	plt.legend()
#	plt.show()
	outfile4 = os.path.join(output, "ATCG_content.png")
	plt.savefig(outfile4)
	plt.close()

def N_count_show(GC, read_num, read_length, output):
	N = []
	for key, value in GC.items():
		N.append(value[4]/read_num*100)
	plt.figure(figsize = (20, 16), facecolor = "white")
	plt.plot(GC.keys(), N, color='red', label = "%N", marker='o', markersize=2)
	plt.gca().set_ylim([0, 100])
	plt.gca().set_xlim([0.5, read_length+0.5])
	plt.gca().xaxis.set_major_locator(ticker.MultipleLocator(1))
	plt.gca().yaxis.set_major_locator(ticker.MultipleLocator(10))
	plt.title("N content across all bases", fontsize = 15)
	plt.xlabel("Position in read (bp)", fontsize = 15)
	for i in range(2, read_length+1, 2):
		plt.axvspan(i-0.5, i+0.5, facecolor = ('lightgrey'), alpha=0.5)
	plt.gca().spines['top'].set_visible(False)
	plt.gca().spines['right'].set_visible(False)
	plt.grid(axis = "y")
	plt.text(read_length+1, 96, "%N", fontsize = 18, color = "red", bbox = {'facecolor': 'white'})
#	plt.show()
	outfile6 = os.path.join(output, "N_content.png")
	plt.savefig(outfile6)
	plt.close()

def GC_content_show(GC_content, GC_percent, base_num, output):
	#整理成字典，键是GC百分含量，值是reads数，在小于最小值和大于最大值处填补0
	reads_GC = {}
	std = np.std(GC_content)
	for key in sorted(GC_content):
		reads_GC[key] = reads_GC.get(key, 0) + 1
	for i in range(1,101):
		if i not in reads_GC.keys() and (i < min(GC_content) or i > max(GC_content)):
			reads_GC[i] = 0
	reads_GC_show= dict(sorted(reads_GC.items(), key=lambda x: x[0]))
#	print(reads_GC_show)
	plt.figure(figsize = (20, 16), facecolor = "white")
	plt.plot(reads_GC_show.keys(),reads_GC_show.values(), label = "GC count per read", color = "red")
	#标准分布
	plt.plot(np.linspace(0,100,100), 1400*norm.pdf(np.linspace(0,100,100),loc = GC_percent, scale = std), color = "blue", label = "Theoretical Distribution")
#	n, bins, patches = plt.hist(GC_content,100)
#	y = ((1 / (np.sqrt(2 * np.pi) * std)) * np.exp(-0.5 * (1 / std * (bins - GC_percent))**2))
#	plt.plot(bins, 1000*y)
	plt.gca().set_ylim([0,max(reads_GC_show.values())+1])
	plt.gca().set_xlim([0,100])
	plt.gca().yaxis.set_major_locator(ticker.MultipleLocator(10))
	plt.gca().xaxis.set_major_locator(ticker.MultipleLocator(2))
	plt.grid(axis = "y")
	plt.title("GC distribution over all sequences", fontsize = 15)
	plt.xlabel("MeanGC content (%)", fontsize =15)
	for i in range(1, 100, 2):
		plt.axvspan(i-0.5, i+0.5, facecolor = ('lightgrey'), alpha=0.5)
	plt.gca().spines['top'].set_visible(False)
	plt.gca().spines['right'].set_visible(False)
	plt.legend(loc=1)
#	plt.show()
	outfile5 = os.path.join(output, "GC_content.png")
	plt.savefig(outfile5)
	plt.close()

def Length_show(read_length, read_num, output):
	X = [read_length-1, read_length, read_length+1]
	Y = [0, read_num, 0]
	plt.figure(figsize = (20, 16), facecolor = "white")
	plt.plot(X, Y, color = "red")
	plt.gca().set_ylim([0,read_num])
	plt.gca().yaxis.set_major_locator(ticker.MultipleLocator(100))
	plt.gca().xaxis.set_major_locator(ticker.MultipleLocator(1))
	plt.grid(axis = "y")
	plt.title("Distribution of sequence lengths over all sequences", fontsize =15)
	plt.xlabel("Sequence Length (bp)", fontsize = 15)
	plt.axvspan(read_length-0.5, read_length+0.5, facecolor = ('lightgrey'), alpha=0.5)
	plt.gca().spines['top'].set_visible(False)
	plt.gca().spines['right'].set_visible(False)
	plt.text(read_length+1, 978, "Sequence Length", fontsize = 15, color = "red", bbox = {'facecolor': 'white'})
#	plt.show()
	outfile7 = os.path.join(output, "Sequence_Length.png")
	plt.savefig(outfile7)
	plt.close()

def main():
	parser=argparse.ArgumentParser(description=__doc__,formatter_class=argparse.RawDescriptionHelpFormatter)
	parser.add_argument('-s','--seq',help='fastq file',required=True)
	parser.add_argument('-o','--outdir',help='output_dir',required=True)
#	parser.add_argument('-o1','--outfile1',help='output filename of Basic Statistics',required=True)
#	parser.add_argument('-o2','--outfile2',help='Per base sequence quality',required=True)
#	parser.add_argument('-o3','--outfile3',help='Per sequence quality scores',required=True)
#	parser.add_argument('-o4','--outfile4',help='Per base sequence content',required=True)
#	parser.add_argument('-o5','--outfile5',help='Per sequence GC content',required=True)
#	parser.add_argument('-o6','--outfile6',help='Per base N content',required=True)
#	parser.add_argument('-o7','--outfile7',help='Sequence Length Distribution',required=True)
#	parser.add_argument('-o8','--outfile8',help='Sequence Duplication Levels',required=True)
#	parser.add_argument('-o9','--outfile9',help='Adapter Content',required=True)
	args=parser.parse_args()

	output = check_outdir(args.outdir)
	read_num, GC_percent, read_length, qual, base_num = base_stat(args.seq)
	show_stat(args.seq, read_num, GC_percent, read_length, output)
	base_qual_show(read_length, qual, read_num, output)
	reads_qual_show(qual, output)
	GC, GC_content, reads = base_GC(args.seq, read_length)
	base_show(GC, read_num, read_length, output)
	GC_content_show(GC_content, GC_percent, base_num, output)
	N_count_show(GC, read_num, read_length, output)
	Length_show(read_length, read_num, output)
	Dup_show, Readsline = Dup(args.seq, read_num)
	Dup_level_show(Dup_show, read_num, output)
	A1, A2, A3, A4, A5 = Adapter(Readsline, read_length)
	Adapter_show(A1, A2, A3, A4, A5, read_length, output)

if __name__=="__main__":
	main()
