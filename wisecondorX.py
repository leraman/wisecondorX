#!/usr/bin/env python
# Copyright (C) 2016 VU University Medical Center Amsterdam
# Author: Roy Straver (github.com/rstraver)
#
# This file is part of WISECONDOR
# WISECONDOR is distributed under the following license:
# Attribution-NonCommercial-ShareAlike, CC BY-NC-SA (https://creativecommons.org/licenses/by-nc-sa/4.0/legalcode)
# This license is governed by Dutch law and this license is subject to the exclusive jurisdiction of the courts of the Netherlands.

import argparse
import sys
from scipy.stats import norm
import json
import os
import math
from wisetools import *

def toolConvert(args):
	converted, qual_info = convertBam(args.infile, binsize=args.binsize, minShift=args.retdist, threshold=args.retthres)
	np.savez_compressed(args.outfile,
						arguments=vars(args),
						runtime=getRuntime(),
						sample=converted,
						quality=qual_info)
	print 'Conversion finished'


def toolNewref(args):

	if args.gender != "M" and args.gender != "F":
		print "ERROR: Unknown gender"
		print "Select \"M\" or \"F\""
		exit(1)

	splitPath = list(os.path.split(args.outfile))
	if splitPath[-1][-4:] == '.npz':
		splitPath[-1] = splitPath[-1][:-4]
	basePath = os.path.join(splitPath[0], splitPath[1])

	# Add single thread information used for parallel processing functions
	args.prepfile = basePath + "_prep.npz"
	args.partfile = basePath + "_part"
	args.parts = args.cpus

	toolNewrefPrep(args)

	# Use multiple cores if requested
	if args.cpus != 1:
		import concurrent.futures
		import copy
		with concurrent.futures.ProcessPoolExecutor(max_workers=args.cpus) as executor:
			for part in xrange(1, args.parts + 1):
				if not os.path.isfile(args.partfile + "_" + str(part) + ".npz"):
					thisArgs = copy.copy(args)
					thisArgs.part = [part, args.parts]
					executor.submit(toolNewrefPart, thisArgs)
			executor.shutdown(wait=True)
	else:
		for part in xrange(1, args.parts + 1):
			if not os.path.isfile(args.partfile + "_" + str(part) + ".npz"):
				args.part = [part, args.parts]
				toolNewrefPart(args)

	# Put it all together
	toolNewrefPost(args)

	# Remove parallel processing temp data
	os.remove(args.prepfile)
	for part in xrange(1, args.parts + 1):
		os.remove(args.partfile + '_' + str(part) + '.npz')


def toolNewrefPrep(args):
	samples = []
	nreads = []
	binsizes = set()
	for infile in args.infiles:  # glob.glob(args.infiles):
		print 'Loading:', infile,
		npzdata = np.load(infile)
		sample = npzdata['sample'].item()
		if args.gender == "M":
			sample = normalizeRefMale(sample)
		print ' \tbinsize:', int(npzdata['arguments'].item()['binsize'])
		samples.append(scaleSample(sample, npzdata['arguments'].item()['binsize'], args.binsize))
		nreads.append(sum([sum(sample[x]) for x in sample.keys()]))
		binsizes.add(npzdata['arguments'].item()['binsize'])

	if args.binsize is None and len(binsizes) != 1:
		print 'ERROR: There appears to be a mismatch in binsizes in your dataset:', binsizes
		print 'Either remove the offending sample or use -binsize to scale all samples'
		exit(1)

	binsize = args.binsize
	if args.binsize is None:
		binsize = binsizes.pop()

	maskedData, chromosomeBins, mask = toNumpyArray(samples, args.gender)
	del samples
	maskedChromBins = [sum(mask[sum(chromosomeBins[:i]):sum(chromosomeBins[:i]) + x]) for i, x in
						enumerate(chromosomeBins)]
	maskedChromBinSums = [sum(maskedChromBins[:x + 1]) for x in range(len(maskedChromBins))]
	correctedData, pca = trainPCA(maskedData)
	np.savez_compressed(args.prepfile,
						arguments=vars(args),
						runtime=getRuntime(),
						binsize=binsize,
						chromosomeBins=chromosomeBins,
						maskedData=maskedData,
						mask=mask,
						maskedChromBins=maskedChromBins,
						maskedChromBinSums=maskedChromBinSums,
						correctedData=correctedData,
						pca_components=pca.components_,
						pca_mean=pca.mean_)


def toolNewrefPart(args):
	if args.part[0] > args.part[1]:
		print 'ERROR: Part should be smaller or equal to total parts:', args.part[0], '>', args.part[
			1], 'is wrong'
		exit(1)
	if args.part[0] < 0:
		print 'ERROR: Part should be at least zero:', args.part[0], '<', 0, 'is wrong'
		exit(1)

	npzdata = np.load(args.prepfile)
	correctedData = npzdata['correctedData']
	maskedChromBins = npzdata['maskedChromBins']
	maskedChromBinSums = npzdata['maskedChromBinSums']

	indexes, distances = getReference(correctedData, maskedChromBins, maskedChromBinSums, args.gender,
										selectRefAmount=args.refsize, part=args.part[0], splitParts=args.part[1])

	np.savez_compressed(args.partfile + '_' + str(args.part[0]) + '.npz',
						arguments=vars(args),
						runtime=getRuntime(),
						indexes=indexes,
						distances=distances)


def toolNewrefPost(args):
	# Load prep file data
	npzdata = np.load(args.prepfile)
	maskedChromBins = npzdata['maskedChromBins']
	chromosomeBins = npzdata['chromosomeBins']
	mask = npzdata['mask']
	pca_components = npzdata['pca_components']
	pca_mean = npzdata['pca_mean']
	binsize = npzdata['binsize'].item()

	# Load and combine part file data
	bigIndexes = []
	bigDistances = []
	for part in xrange(1, args.parts + 1):  # glob.glob(args.infiles):
		infile = args.partfile + '_' + str(part) + '.npz'
		print 'Loading:', infile
		npzdata = np.load(infile)
		bigIndexes.extend(npzdata['indexes'])
		bigDistances.extend(npzdata['distances'])

		print part, npzdata['indexes'].shape

	indexes = np.array(bigIndexes)
	distances = np.array(bigDistances)

	np.savez_compressed(args.outfile,
						arguments=vars(args),
						runtime=getRuntime(),
						binsize=binsize,
						indexes=indexes,
						distances=distances,
						chromosome_sizes=chromosomeBins,
						mask=mask,
						masked_sizes=maskedChromBins,
						pca_components=pca_components,
						pca_mean=pca_mean,
						gender=args.gender)


# Most of this "tool" is actually basic functionality and should be in wisetools.py instead (and still is in WCX!)
def toolTest(args):

	time_at_start = datetime.datetime.now()
	start_time = [time_at_start.year, time_at_start.month, time_at_start.day,
				  time_at_start.hour, time_at_start.minute, time_at_start.second]

	if not args.json and not args.txt and not args.plot:
		print "ERROR: No output format selected"
		print "Select at least one of the supported output formats (-json, -txt, -plot)"
		exit(1)

	# Reference data handling
	mask_list = []
	referenceFile = np.load(args.reference)

	gender = referenceFile['gender']

	binsize = referenceFile['binsize'].item()
	indexes = referenceFile['indexes']
	distances = referenceFile['distances']
	chromosome_sizes = referenceFile['chromosome_sizes']
	mask = referenceFile['mask']
	mask_list.append(mask)
	masked_sizes = referenceFile['masked_sizes']
	maskedChromBinSums = [sum(masked_sizes[:x + 1]) for x in range(len(masked_sizes))]

	pca_mean = referenceFile['pca_mean']
	pca_components = referenceFile['pca_components']

	del referenceFile

	# Test sample data handling
	sampleFile = np.load(args.infile)
	sample = sampleFile['sample'].item()

	nreads = sum([sum(sample[x]) for x in sample.keys()])
	sampleBinSize = sampleFile['arguments'].item()['binsize']
	sample = scaleSample(sample, sampleBinSize, binsize)

	testData = toNumpyRefFormat(sample, chromosome_sizes, mask, gender)
	testData = applyPCA(testData, pca_mean, pca_components)
	optimalCutoff, cutOffMask = getOptimalCutoff(distances, args.maskrepeats)

	num_tests = sum(masked_sizes)
	z_threshold = norm.ppf(1 - 1. / (num_tests * 0.5))

	print 'Per bin z-score threshold for first testing cycles:', z_threshold

	testCopy = np.copy(testData)
	resultsZ, resultsR, refSizes, stdDevAvg = repeatTest(testCopy, indexes, distances, masked_sizes, maskedChromBinSums,
														 optimalCutoff, z_threshold, 5)
	print 'ASDES:', stdDevAvg, '\nAASDEF:', stdDevAvg * z_threshold
	# Get rid of infinite values caused by having no reference bins or only zeros in the reference
	infinite_mask = (refSizes >= args.minrefbins)
	mask_list.append(infinite_mask)
	cleanedR = resultsR[infinite_mask]
	cleanedZ = resultsZ[infinite_mask]
	cleanedBinSums = [np_sum(infinite_mask[:val]) for val in maskedChromBinSums]
	cleanedBins = [cleanedBinSums[i] - cleanedBinSums[i - 1] for i in range(1, len(cleanedBinSums))]
	cleanedBins.insert(0, cleanedBinSums[0])

	resultsZ = []
	resultsR = []
	inflatedZ = inflateArrayMulti(cleanedZ, mask_list)
	inflatedR = inflateArrayMulti(np.log2(cleanedR) - np.median(np.log2(cleanedR)), mask_list)
	for chrom in xrange(len(chromosome_sizes)):
		chromData = inflatedZ[sum(chromosome_sizes[:chrom]):sum(chromosome_sizes[:chrom + 1])]
		resultsZ.append(chromData)
		chromData = inflatedR[sum(chromosome_sizes[:chrom]):sum(chromosome_sizes[:chrom + 1])]
		resultsR.append(chromData)

	# Apply blacklist
	if args.blacklist != None:
		applyBlacklist(args, binsize, resultsR, resultsZ, sample, gender)

	# Make R interpretable
	resultsR = [x.tolist() for x in resultsR]
	nchrs = len(resultsR)
	for c in range(nchrs):
		for i, rR in enumerate(resultsR[c]):
			if not np.isfinite(rR):
				resultsR[c][i] = 0 # 0 -> not interesting
	resultsZ = [x.tolist() for x in resultsZ]
	for c in range(nchrs):
		for i, rR in enumerate(resultsZ[c]):
			if not np.isfinite(rR):
				resultsZ[c][i] = 0 # 0 -> not interesting

	cbs_calls = getBMScore(args, resultsR, gender)

	time_at_end = datetime.datetime.now()
	end_time = [time_at_end.year, time_at_end.month, time_at_end.day,
				time_at_end.hour, time_at_end.minute, time_at_end.second]

	json_out = {'binsize': binsize,
				'results_r': resultsR,
				'results_z': resultsZ,
				'threshold_z': z_threshold.tolist(),
				'asdef': stdDevAvg.tolist(),
				'aasdef': (stdDevAvg * z_threshold).tolist(),
				'nreads': nreads,
				'runtime': [start_time, end_time],
				'cbs_calls': cbs_calls}

	json_file = open(args.outid + ".json", "w")
	json.dump(json_out, json_file)
	json_file.close()

	# Save txt: optional
	if args.txt:
		generateTxtOuts(args, binsize, json_out)

	# Plot: optional
	if args.plot:
		plot_script = str(os.path.dirname(os.path.realpath(__file__))) + "/plotter.R"
		if gender == "M":
			sexchrom = "XY"
		else:
			sexchrom = "X"
		os.popen("Rscript \"" + plot_script + "\" --infile \"" +args.outid + ".json\" --outdir \"" +
				 args.outid + ".plots\"" + " --sexchroms " + sexchrom)

	if not args.json:
		os.remove(args.outid + ".json")

	print("Done!")
	exit(0)


def generateTxtOuts(args, binsize, json_out):
	output_dir = os.path.dirname(os.path.abspath(args.outid))

	bed_file = open(output_dir + "/" + args.outid + ".bed","w")
	bed_file.write("track name=\"" + args.outid + "\" description=\"rscore\"\n")
	resultsR = json_out["results_r"]
	for chr_i in range(len(resultsR)):
		chr = str(chr_i + 1)
		if chr == "23":
			chr = "X"
		if chr == "24":
			chr = "Y"
		feat = 1
		for feat_i in range(len(resultsR[chr_i])):
			r = resultsR[chr_i][feat_i]
			feat_str = chr + ":" + str(feat) + "-" + str(feat + binsize - 1)
			it = [chr, feat - 1, feat + binsize - 1, feat_str, r, "+"]
			it = [str(x) for x in it]
			bed_file.write("\t".join(it) + "\n")
			feat += binsize
	bed_file.close()

	segments_file = open(output_dir + "/" + args.outid + "_segments.txt","w")
	for segment in json_out["cbs_calls"]:
		chr = str(int(segment[0]))
		if chr == "23":
			chr = "X"
		if chr == "24":
			chr = "Y"
		it = [chr, int(segment[1] * binsize + 1), int((segment[2] + 1) * binsize), segment[4]]
		it = [str(x) for x in it]
		segments_file.write("\t".join(it) + "\n")
	segments_file.close()

	statistics_file = open(output_dir + "/statistics.out","w")
	chrom_scores = []
	for chr_i in range(len(resultsR)):
		chr = str(chr_i + 1)
		if chr == "23":
			chr = "X"
		if chr == "24":
			chr = "Y"
		chrom_score = np.mean(resultsR[chr_i])
		statistics_file.write("\"chr" + str(chr) + "\" " + str(chrom_score) + "\n")
		chrom_scores.append(chrom_score)

	statistics_file.write("\"stdev\" " + str(np.std(chrom_scores)))
	statistics_file.close()


def applyBlacklist(args, binsize, resultsR, resultsZ, sample, gender):
	blacklist = {}

	for line in open(args.blacklist):
		bchr, bstart, bstop = line.strip().split("\t")
		bchr = bchr[3:]
		if bchr not in blacklist.keys():
			blacklist[bchr] = []
		blacklist[bchr].append([int(int(bstart) / binsize), int(int(bstop) / binsize) + 1])

	for chr in blacklist.keys():
		for s_s in blacklist[chr]:
			if chr == "X":
				chr = "23"
			if chr == "Y":
				chr = "24"
			for pos in range(s_s[0], s_s[1]):
				if gender == "F" and chr == "24":
					continue
				resultsR[int(chr) - 1][pos] = 0
				resultsZ[int(chr) - 1][pos] = 0
				sample[chr][pos] = 0


def getBMScore(args, resultsR, gender):
	json_cbs_temp_dir = os.path.abspath(args.outid + "_CBS_tmp")
	json_cbs_file = open(json_cbs_temp_dir + "_01.json", "w")
	json.dump({"results_r": resultsR}, json_cbs_file)
	json_cbs_file.close()
	CBS_script = str(os.path.dirname(os.path.realpath(__file__))) + "/CBS.R"

	if gender == "M":
		sexchrom = "XY"
	else:
		sexchrom = "X"

	os.popen(
		"Rscript \"" + CBS_script + "\" --infile \"" + json_cbs_temp_dir + "_01.json\" --outfile \"" +
		json_cbs_temp_dir + "_02.json\"" + " --sexchroms " + sexchrom)
	os.remove(json_cbs_temp_dir + "_01.json")
	cbs_data = json.load(open(json_cbs_temp_dir + "_02.json"))[1:]
	cbs_data = [[float(y.encode("utf-8")) for y in x] for x in cbs_data]
	os.remove(json_cbs_temp_dir + "_02.json")

	BM_scores = []
	for cbs_call_index in range(len(cbs_data[0])):
		chr_i = int(cbs_data[0][cbs_call_index]) - 1
		start = int(cbs_data[1][cbs_call_index]) - 1
		end = int(cbs_data[2][cbs_call_index])  # no - 1! (closed interval in python)

		BM_score = np.median(resultsR[chr_i][start:end])
		if math.isnan(BM_score):
			BM_score = 0.0
		BM_scores.append(BM_score)

	# Save results

	cbs_calls = []
	for cbs_call_index in range(len(cbs_data[0])):
		cbs_calls.append(
			[cbs_data[0][cbs_call_index], cbs_data[1][cbs_call_index] - 1, cbs_data[2][cbs_call_index] - 1,
			 BM_scores[cbs_call_index], cbs_data[4][cbs_call_index]])
	return cbs_calls


def reformat(args):
	original = np.load(args.infile)
	sample_data = original['sample']
	updated_sample_data = dict()
	chrs = [str(x) for x in range(1, 23)]
	chrs += ["X", "Y"]
	i = 1
	for chr in chrs:
		updated_sample_data[str(i)] = sample_data.item()[chr]
		i += 1

	np.savez_compressed(args.outid, arguments=original['arguments'], runtime=original['runtime'],
						sample=updated_sample_data, quality=original['quality'])

def get_gender(args):
	npzfile = np.load(args.infile)
	sample = npzfile['sample'].item()
	nonY = float(sum([np.sum(sample[str(chr)]) for chr in range(1,24)]))
	Y = float(np.sum(sample["24"]))
	permilleY = Y / (nonY + Y) * 1000.0
	if permilleY > args.cutoff:
		print("Gender = M")
	else:
		print("Gender = F")
	exit(1)


def main():

	parser = argparse.ArgumentParser(
		description="wisecondorX")
	subparsers = parser.add_subparsers()

	# File conversion
	parser_convert = subparsers.add_parser('convert',
		description='Convert and filter a bam file to an npz')
	parser_convert.add_argument('infile',
		type=str,
		help='Bam input file for conversion')
	parser_convert.add_argument('outfile',
		type=str,
		help='File to write binned information to')
	parser_convert.add_argument('-binsize',
		type=float, default=5e3,
		help='Size per bin in bp')
	parser_convert.add_argument('-retdist',
		type=int, default=4,
		help='Maximum amount of base pairs difference between sequential reads '
			 'to consider them part of the same tower')
	parser_convert.add_argument('-retthres',
		type=int, default=4,
		help='Threshold for when a group of reads is considered a tower and will be removed')
	parser_convert.set_defaults(func=toolConvert)

	# Reformat
	parser_reformat = subparsers.add_parser('reformat',
		description='Reformat an original .npz to a wisecondorX .npz')
	parser_reformat.add_argument('infile',
		type=str, help='.npz input file for conversion')
	parser_reformat.add_argument('outfile',
		type=str, help='resulting .npz output file')
	parser_reformat.set_defaults(func=reformat)

	# Find gender
	parser_gender = subparsers.add_parser('gender',
        description='Predict the gender of a sample')
	parser_gender.add_argument('infile',
		type=str, help='.npz input file, resulting from .bam conversion')
	parser_gender.add_argument('-cutoff',
		type=float, default=2.5, help='Y-read permille cut-off. Below is female, above is male')
	parser_gender.set_defaults(func=get_gender)

	# New reference creation
	parser_newref = subparsers.add_parser('newref',
		description='Create a new reference using healthy reference samples')
	parser_newref.add_argument('infiles',
		type=str, nargs='+',
		help='Path to all reference data files (i.e. ./reference/*.npz)')
	parser_newref.add_argument('outfile',
		type=str,
		help='Path and filename for the reference output (i.e. ./reference/myref.npz)')
	parser_newref.add_argument('-refsize',
		type=int, default=300,
		help='Amount of reference locations per target')
	parser_newref.add_argument('-binsize',
		type=int, default=1e5,
		help='Try to scale samples to this binsize, multiples of existing binsize only')
	parser_newref.add_argument('-gender',
		type=str, default="F",
		help='Gender of the reference data files')
	parser_newref.add_argument('-cpus',
		type=int, default=1,
		help='EXPERIMENTAL: Use multiple cores to find reference bins')
	parser_newref.set_defaults(func=toolNewref)

	# Find CNAs
	parser_test = subparsers.add_parser('predict',
		description='Find Copy Number Aberrations')
	parser_test.add_argument('infile',
		type=str,
		help='Sample.npz of which the CNAs will be predicted')
	parser_test.add_argument('reference',
		type=str,
		help='Reference as previously created')
	parser_test.add_argument('outid',
		type=str,
		help='Basename (w/o extension) of output files')
	parser_test.add_argument('-minrefbins',
		type=int, default=150,
		help='Minimum amount of sensible ref bins per target bin')
	parser_test.add_argument('-maskrepeats',
		type=int, default=5,
		help='Regions with distances > mean + sd * 3 will be masked. Number of masking cycles')
	parser_test.add_argument('-blacklist',type=str, default=None,
		help='Blacklist that masks regions in output, structure of header-less'
		     'file chr...(/t)startpos(/t)endpos(/n)')
	parser_test.add_argument('-json',
		action="store_true",
		help='Outputs .json file, containing all generated information')
	parser_test.add_argument('-txt',
		action="store_true",
		help='Outputs tab-delimited text files, containing the most important information')
	parser_test.add_argument('-plot',
		action="store_true",
		help='Outputs .png plots')
	parser_test.set_defaults(func=toolTest)


	args = parser.parse_args(sys.argv[1:])
	printArgs(args)
	args.func(args)


if __name__ == '__main__':
	import warnings
	warnings.filterwarnings("ignore")
	main()