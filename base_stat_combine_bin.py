#!/usr/bin/env python3
import os
import sys
import re
import argparse
from itertools import product
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.collections as collections
import matplotlib.ticker as ticker

bindir = os.path.abspath(os.path.dirname(__file__))
__author__='snail~'
__mail__= 'snail1993@outlook.com'
__doc__='base stat'

def base_stat(seq,outname):
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
#	print(qual)
	GC_percent = round(GC_num/base_num*100)
	return(read_num, GC_percent, length, qual)

def base_qual_show(length, qual, read_num):
	#按碱基位置统计质量
	##字典形式
	base_qual = {}
	for base in range(1,length[0]+1):
		base_qual[base] = []
		for r in range(1,read_num+1):
			base_qual[base].append(qual[r][base-1])
	#绘图
	##整理数据，X：base_pos; Y: quality
	base_qual_motified = {}
	i = 1
	if len(list(base_qual.keys())) > 15:
		while i < len(list(base_qual.keys())):
			if i < 10:
				base_qual_motified[i] = base_qual[i]
				i += 1
			else:
				if i+3 < (len(list(base_qual.keys())) +1):
					temp = base_qual[i] + base_qual[i+1] + base_qual[i+2]
					base_qual_motified[i] = temp
					i += 3
				else:
#					K = str(i) + '-' + str(len(list(base_qual.keys())))
					n = len(list(base_qual.keys())) - i
					temp = []
					for m in range(1, n, 1):
						temp += base_qual[m]
					base_qual_motified[i] = temp
					break

	#统计每个位点质量均值
	base_qual_mean = []
	for key in base_qual_motified.keys():
		base_qual_mean.append(np.mean(base_qual_motified[key]))
	print(base_qual_motified.keys(), base_qual_mean)
	plt.figure(figsize = (20, 16), facecolor = "white")
	#sym=''不显示异常点，whiskerprops须颜色，boxprops箱体线以及箱内部颜色，medianprops众数线，capprops极值线，whis极值范围
	plt.boxplot(base_qual_motified.values(), base_qual_motified.keys(), patch_artist = True, 
						whiskerprops = {'color': "black"}, boxprops = {'color':'black', 'facecolor':'yellow'}, medianprops = {'color':'red'}, capprops = {'color': "black"}, whis = [10,90], 
						sym='')
	#加均值线
	plt.plot(range(1,len(base_qual_motified.keys())+1),base_qual_mean)
	#纵轴分隔单位2
	plt.gca().yaxis.set_major_locator(ticker.MultipleLocator(2))
	plt.ylim(0,41)
	plt.xticks(range(1,len(base_qual_motified.keys())+1), list(base_qual_motified.keys()))
	#设置背景颜色，axhspan横向，axvspan纵向
	plt.axhspan(0, 20, facecolor = ('#fdc1c5'))
	plt.axhspan(20, 28, facecolor = ('#fbeeac'))
	plt.axhspan(28, 41, facecolor = ('#bcf5a6'))
	for i in range(1, len(base_qual_motified.keys()), 2):
		plt.axvspan(i-0.5, i+0.5, facecolor = ('whitesmoke'), alpha=0.3)
	plt.xlabel("Position in read (bp)")
	plt.title("Quality scores across all bases (Sanger / Illumina 1.9 encoding)")
	plt.show()

def main():
	parser=argparse.ArgumentParser(description=__doc__,formatter_class=argparse.RawDescriptionHelpFormatter)
	parser.add_argument('-s','--seq',help='fastq file',required=True)
	parser.add_argument('-o','--outfile',help='output filename',required=True)
	args=parser.parse_args()

	read_num, GC_percent, length, qual = base_stat(args.seq, args.outfile)
	base_qual_show(length, qual, read_num)

if __name__=="__main__":
	main()
