# Copyright (C) 2016 VU University Medical Center Amsterdam
# Author: Roy Straver (github.com/rstraver)
#
# This file is part of WISECONDOR
# WISECONDOR is distributed under the following license:
# Attribution-NonCommercial-ShareAlike, CC BY-NC-SA (https://creativecommons.org/licenses/by-nc-sa/4.0/legalcode)
# This license is governed by Dutch law and this license is subject to the exclusive jurisdiction of the courts of the Netherlands.

import matplotlib
matplotlib.use('agg')
import time
import bisect
from sklearn.decomposition import PCA
from sklearn.utils.extmath import fast_dot
import pysam
from triarray import *
import getpass
import datetime
import socket

# Get rid of some useless warnings, I know there are emtpy slices
import warnings
warnings.filterwarnings('ignore', 'Mean of empty slice')
warnings.filterwarnings('ignore', 'Degrees of freedom <= 0 for slice')

curTime = datetime.datetime.now()

np.seterr('ignore')
np_sum = np.sum
np_pow = np.power
np_max = np.argmax
np_mean = np.mean
np_median = np.median
np_std = np.std
np_abs = np.abs
np_sqrt = np.sqrt

find_pos = bisect.bisect


def getRuntime():
	runtime = dict()
	runtime['datetime']=curTime
	runtime['hostname']=socket.gethostname()
	runtime['username']=getpass.getuser()
	return runtime

def printArgs(args):
	argdict=vars(args)
	print 'tool =', str(argdict['func']).split()[1][4:]
	for arg in sorted(argdict.keys()):
		if arg != 'func':
			print arg,'=',argdict[arg]


def loadCytoBands(cytoFile):
	cytoDict = dict()
	curChrom = None
	curCyto = None
	with open(cytoFile, 'r') as cytoData:
		for line in cytoData:
			splitLine = line.split()
			if splitLine[0][3:] != curChrom:
				if curChrom != None:
					cytoDict[curChrom] = curCyto
				curCyto = []
				curChrom = splitLine[0][3:]
			curCyto.append(splitLine[1:])
	return cytoDict


def trainPCA(refData,pcacomp=3):
	tData = refData.T
	pca = PCA(n_components=pcacomp)
	pca.fit(tData)
	PCA(copy=True, whiten=False)
	transformed = pca.transform(tData)
	inversed = pca.inverse_transform(transformed)
	corrected = tData / inversed

	return corrected.T, pca


def applyPCA(sampleData, mean, components):
	pca = PCA(n_components=components.shape[0])
	pca.components_ = components
	pca.mean_ = mean

	transform = pca.transform(np.array([sampleData]))

	reconstructed = fast_dot(transform, pca.components_) + pca.mean_
	reconstructed = reconstructed[0]
	return sampleData / reconstructed


def convertBam(bamfile, binsize=100000, minShift=4, threshold=6, mapq=1, demandPair=False):
	# Prepare the list of chromosomes
	chromosomes = dict()
	for chromosome in range(1, 25):
		chromosomes[str(chromosome)] = None

	# Flush the current stack of reads
	def flush(readBuff, counts):
		stairSize = len(readBuff)
		if stairSize <= threshold or threshold < 0:
			for read in readBuff:
				location = read.pos / binsize
				counts[int(location)] += 1
				#if location >= len(counts):
				#	print read

	sam_file = pysam.AlignmentFile(bamfile, "rb")
	reads_seen = 0
	reads_kept = 0
	reads_mapq = 0
	reads_rmdup = 0
	reads_pairf = 0
	larp = -1 # LAst Read Position...
	larp2 = -1

	for index,chrom in enumerate(sam_file.references):

		chromName = chrom
		if chromName[:3].lower() == 'chr':
			chromName = chromName[3:]
		if chromName not in chromosomes and chromName != "X" and chromName != "Y":
			continue

		print chrom,'length:', sam_file.lengths[index], 'bins:', int(sam_file.lengths[index] / float(binsize) + 1)
		counts = np.zeros(int(sam_file.lengths[index] / float(binsize) + 1), dtype=np.int32)

		readBuff = []
		sam_iter = sam_file.fetch(chrom)

		if chromName == 'X':
			chromName = '23'
		if chromName == 'Y':
			chromName = '24'

		prevRead = sam_iter.next()

		# Split paths here, for-loop was heavily slowed down by if-statements otherwise
		if demandPair:
			for read in sam_iter:
				if ((int(read.pos) - int(prevRead.pos)) > minShift):
					flush(readBuff, counts)
					readBuff = []
				# Normal ndups will be appended here

				if not (read.is_proper_pair and read.is_read1):
					reads_pairf += 1
					continue

				if larp == read.pos and larp2 == read.next_reference_start:
						reads_rmdup += 1
				else:
					if read.mapping_quality >= mapq:
						readBuff.append(read)
						prevRead = read
					else:
						reads_mapq += 1

				larp2 = read.next_reference_start

				reads_seen += 1
				larp = read.pos
		else:
			for read in sam_iter:
				if ((int(read.pos) - int(prevRead.pos)) > minShift):
					flush(readBuff, counts)
					readBuff = []
				# Normal ndups will be appended here

				if larp == read.pos:
						reads_rmdup += 1
				else:
					if read.mapping_quality >= mapq:
						readBuff.append(read)
						prevRead = read
					else:
						reads_mapq += 1

				reads_seen += 1
				larp = read.pos

		# Flush after we're done
		flush(readBuff, counts)
		chromosomes[chromName] = counts
		reads_kept += sum(counts)

	#print reads_seen,reads_kept
	qual_info = {'mapped':sam_file.mapped,
				 'unmapped':sam_file.unmapped,
				 'no_coordinate':sam_file.nocoordinate,
				 'filter_rmdup':reads_rmdup,
				 'filter_mapq':reads_mapq,
				 'pre_retro':reads_seen,
				 'post_retro':reads_kept,
				 'pair_fail':reads_pairf}
	return chromosomes, qual_info


def scaleSample(sample, fromSize, toSize):
	if fromSize == toSize or toSize == None:
		return sample

	if toSize == 0 or fromSize == 0 or toSize < fromSize or toSize % fromSize > 0:
		print 'ERROR: Impossible binsize scaling requested:', fromSize, 'to', toSize
		exit(1)

	returnSample = dict()
	scale = toSize/fromSize
	for chrom in sample:
		chromData = sample[chrom]
		newLen = int(np.ceil(len(chromData)/float(scale)))
		scaledChrom = np.zeros(newLen, dtype=np.int32)
		for i in range(newLen):
			scaledChrom[i] = np_sum(chromData[int(i*scale):int(i*scale+scale)])
			returnSample[chrom] = scaledChrom
	return returnSample


def toNumpyArray(samples, gender):
	byChrom = []
	chromBins = []
	sampleCount = len(samples)

	if gender == "F":
		chromosomes = range(1, 24)
	else:
		chromosomes = range(1, 25)

	for chromosome in chromosomes:
		maxLen = max([sample[str(chromosome)].shape[0] for sample in samples])
		thisChrom = np.zeros((maxLen, sampleCount), dtype=float)
		chromBins.append(maxLen)
		i = 0
		for sample in samples:
			thisChrom[:, i] = sample[str(chromosome)]
			i += 1
		byChrom.append(thisChrom)
	allData = np.concatenate(byChrom, axis=0)

	sumPerSample = np_sum(allData, 0)
	allData = allData / sumPerSample

	print 'Applying nonzero mask on the data:', allData.shape,
	sumPerBin = np_sum(allData, 1)
	mask = sumPerBin > 0
	maskedData = allData[mask, :]
	print 'becomes',maskedData.shape

	return maskedData, chromBins, mask


def toNumpyRefFormat(sample, chromBins, mask, gender):
	byChrom = []

	if gender == "M":
		chrs = range(1, 25)
	else:
		chrs = range(1, 24)

	for chromosome in chrs:
		thisChrom = np.zeros(chromBins[chromosome - 1], dtype=float)
		minLen = min(chromBins[chromosome - 1], len(sample[str(chromosome)]))
		thisChrom[:minLen] = sample[str(chromosome)][:minLen]
		byChrom.append(thisChrom)
	allData = np.concatenate(byChrom, axis=0)
	allData = allData / np_sum(allData)
	maskedData = allData[mask]

	return maskedData


def inflateArray(array, mask):
	temp = np.zeros(mask.shape[0])
	j = 0
	for i, val in enumerate(mask):
		if val:
			temp[i] = array[j]
			j += 1
	return temp


def inflateArrayMulti(array, mask_list):
	temp = array
	for mask in reversed(mask_list):
		temp = inflateArray(temp, mask)
	return temp


def getRefForBins(amount, start, end, sampleData, otherData, knit_length):
	refIndexes = np.zeros((end - start, amount), dtype=np.int32)
	refDistances = np.ones((end - start, amount))
	for thisBin in xrange(start, end):
		thisMask = np_sum(np_pow(otherData - sampleData[thisBin, :], 2), 1)
		thisIndexes = [-1 for i in xrange(amount)]
		thisDistances = [1e10 for i in xrange(amount)]
		removeIndex = thisIndexes.pop
		removeDist = thisDistances.pop
		insertIndex = thisIndexes.insert
		insertDist = thisDistances.insert
		curMax = 1e10
		totLen = len(thisMask) + (end - start)
		for i, binVal in enumerate(thisMask):
			if (end - start) + i >= totLen - knit_length: # skip XY referenced bins
				continue
			if binVal < curMax:
				pos = find_pos(thisDistances, binVal)
				removeIndex(-1)
				removeDist(-1)
				insertIndex(pos, i)
				insertDist(pos, binVal)
				curMax = thisDistances[-1]
		refIndexes[thisBin - start, :] = thisIndexes
		refDistances[thisBin - start, :] = thisDistances
	return refIndexes, refDistances


def getOptimalCutoff(reference, repeats):
	optimalCutoff = float("inf")
	mask = np.zeros(reference.shape)
	for i in range(0, repeats):
		mask = reference < optimalCutoff
		average = np.average(reference[mask])
		stddev = np.std(reference[mask])
		optimalCutoff = average + 3 * stddev
	return optimalCutoff, mask


# Returns: Chromosome index, startBinNumber, endBinNumber
def splitByChrom(start, end, chromosomeBinSums):
	areas = []
	tmp = [0, start, 0]
	for i, val in enumerate(chromosomeBinSums):
		tmp[0] = i
		if val >= end:
			break
		if start < val < end:
			tmp[2] = val
			areas.append(tmp)
			tmp = [i, val, 0]
		tmp[1] = val
	tmp[2] = end
	areas.append(tmp)
	return areas


# Returns: Start and end bin numbers this instance should work on
def getPart(partnum, outof, bincount):
	startBin = int(bincount / float(outof) * partnum)
	endBin = int(bincount / float(outof) * (partnum + 1))
	return startBin, endBin

def normalizeRefMale(sample):
	sexChrs = ["23", "24"]
	for sexChr in sexChrs:
		sample[sexChr] = sample[sexChr] * 2
	return sample

def getReference(correctedData, chromosomeBins, chromosomeBinSums, gender, selectRefAmount=100, part=1, splitParts=1):
	timeStartSelection = time.time()
	bigIndexes = []
	bigDistances = []

	bincount = chromosomeBinSums[-1]

	startNum, endNum = getPart(part - 1, splitParts, bincount)
	print 'Working on thread', part, 'of', splitParts, 'meaning bins', startNum, 'up to', endNum
	regions = splitByChrom(startNum, endNum, chromosomeBinSums)

	for region in regions:
		chrom = region[0]
		start = region[1]
		end = region[2]

		if startNum > start:
			start = startNum
		if endNum < end:
			end = endNum

		print "Thread", part, '| Actual Chromosome area', chromosomeBinSums[chrom] - chromosomeBins[chrom], chromosomeBinSums[chrom],\
			"| chr" + str(chrom+1)
		chromData = np.concatenate((correctedData[:chromosomeBinSums[chrom] - chromosomeBins[chrom], :],
									correctedData[chromosomeBinSums[chrom]:, :]))

		X_length = chromosomeBinSums[22] - (chromosomeBinSums[22] - chromosomeBins[22]) # index 22 -> chrX

		if gender == "M":
			Y_length = chromosomeBinSums[23] - (chromosomeBinSums[23] - chromosomeBins[23])
			if chrom == 22:
				knit_length = Y_length
			elif chrom == 23:
				knit_length = X_length
			else:
				knit_length = X_length + Y_length
		else:
			if chrom == 22:
				knit_length = 0
			else:
				knit_length = X_length

		partIndexes, partDistances = getRefForBins(selectRefAmount, start, end, correctedData, chromData, knit_length)
		bigIndexes.extend(partIndexes)
		bigDistances.extend(partDistances)

		print part, 'Time spent:', int(time.time() - timeStartSelection), 'seconds'

	indexArray = np.array(bigIndexes)
	distanceArray = np.array(bigDistances)

	return indexArray, distanceArray

def trySample(testData, testCopy, indexes, distances, chromosomeBins, chromosomeBinSums, cutoff):
	bincount = chromosomeBinSums[-1]

	resultsZ = np.zeros(bincount)
	resultsR = np.zeros(bincount)
	refSizes = np.zeros(bincount)
	stdDevSum = 0.
	stdDevNum = 0
	i = 0

	for chrom in xrange(len(chromosomeBins)):
		start = chromosomeBinSums[chrom] - chromosomeBins[chrom]
		end = chromosomeBinSums[chrom]
		chromData = np.concatenate(
			(testCopy[:chromosomeBinSums[chrom] - chromosomeBins[chrom]], testCopy[chromosomeBinSums[chrom]:]))

		for index in indexes[start:end]:
			refData = chromData[index[distances[i] < cutoff]]
			refData = refData[refData >= 0]  # Previously found aberrations may be marked by negative values
			refMean = np_mean(refData)
			refStdDev = np_std(refData)
			if not np.isnan(refStdDev):
				stdDevSum += refStdDev
				stdDevNum += 1
			resultsZ[i] = (testData[i] - refMean) / refStdDev
			resultsR[i] = testData[i] / refMean
			refSizes[i] = refData.shape[0]
			i += 1
	return resultsZ, resultsR, refSizes, stdDevSum/stdDevNum


def repeatTest(testData, indexes, distances, chromosomeBins, chromosomeBinSums, cutoff, threshold, repeats):
	timeStartTest = time.time()
	resultsZ = None
	resultsR = None
	testCopy = np.copy(testData)
	for i in xrange(repeats):
		resultsZ, resultsR, refSizes, stdDevAvg = trySample(testData, testCopy, indexes, distances, chromosomeBins,
												  chromosomeBinSums, cutoff)
		testCopy[np_abs(resultsZ) >= threshold] = -1
	print 'Time spent on obtaining z-scores:', int(time.time() - timeStartTest), 'seconds'
	return resultsZ, resultsR, refSizes, stdDevAvg