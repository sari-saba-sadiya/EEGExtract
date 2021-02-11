import bisect
import numpy as np
import pandas as pd
import pywt
from scipy import stats, signal, integrate
from dit.other import tsallis_entropy
import dit
import librosa
import statsmodels.api as sm
import itertools
from pyinform import mutualinfo
from statsmodels import tsa
from sklearn.metrics import mutual_info_score
import numpy as np
from scipy import signal,integrate
from sklearn.metrics.cluster import normalized_mutual_info_score as normed_mutual_info 

################################################
#	Auxiliary Functions
################################################

##########
# Filter the eegData, midpass filter 
#	eegData: 3D np array [chans x ms x epochs] 
def filt_data(eegData, lowcut, highcut, fs, order=7):
    nyq = 0.5 * fs
    low = lowcut / nyq
    high = highcut / nyq
    b, a = signal.butter(order, [low, high], btype='band')
    filt_eegData = signal.lfilter(b, a, eegData, axis = 1)
    return filt_eegData

#########
# remove short bursts / spikes 
def fcnRemoveShortEvents(z,n):
    for chan in range(z.shape[0]):
        # check for too-short suppressions
        ct=0
        i0=1
        i1=1 
        for i in range(2,len(z[chan,:])):
            if z[chan,i]==z[chan,i-1]:
                ct=ct+1
                i1=i
            else:
                if ct<n:
                    z[chan,i0:i1] = 0
                    z[chan,i1] = 0 #nasty little bug
                ct=0
                i0=i
                i1=i
        if z[chan,0] == 1 and z[chan,1] == 0:
            z[chan,0] = 0
    return z

##########
# Find interval of consistent values in binary 1D numpy array
def get_intervals(A,B,endIdx=500):
    # This function gives you intervals (a1,b1), (a2,b3) for every a in A=[a1,a2,a3,..]
    # and the smallest element in b that is larger than a.
    intervals = []
    for ii,A_idx_lst in enumerate(A):
        B_idx_lst = [bisect.bisect_left(B[ii], idx) for idx in A_idx_lst]
        chan_intervals = []
        for jj,idx_l in enumerate(B_idx_lst):
            if idx_l == len(B[ii]):
                chan_intervals.append((A_idx_lst[jj],endIdx))
            else:
                chan_intervals.append((A_idx_lst[jj],B[ii][idx_l]))
        intervals.append(chan_intervals)
        # previous code already takes care of the [] possibility
        #if B_idx_lst == []:
        #    intervals.append([])
    return intervals

##########
# Detect bursts and supressions in eeg data
def burst_supression_detection(x,fs,suppression_threshold = 10):
	'''
	# DETECT EMG ARTIFACTS.
	nyq = 0.5 * fs
	low = low / nyq
	high = high / nyq
	be, ae = signal.butter(order, [low, high], btype='band')
	'''
	# CALCULATE ENVELOPE
	e = abs(signal.hilbert(x,axis=1));
	# same as smooth(e,Fs/4) in MATLAB, apply 1/2 second smoothing
	ME = np.array([np.convolve(el,np.ones(int(fs/4))/(fs/4),'same') for el in e.tolist()])
	e = ME
	# DETECT SUPRESSIONS
	# apply threshold -- 10uv
	z = (ME<suppression_threshold)
	# remove too-short suppression segments
	z = fcnRemoveShortEvents(z,fs/2)
	# remove too-short burst segments
	b = fcnRemoveShortEvents(1-z,fs/2)
	z = 1-b
	went_high = [np.where(np.array(chD[:-1]) < np.array(chD[1:]))[0].tolist() for chD in z.tolist()]
	went_low = [np.where(np.array(chD[:-1]) > np.array(chD[1:]))[0].tolist() for chD in z.tolist()]

	bursts = get_intervals(went_high,went_low)
	supressions = get_intervals(went_low,went_high)

	return bursts,supressions

##########
# Coherence in the Delta Band
def CoherenceDelta(eegData, i, j, fs=100):
    nfft=eegData.shape[1]
    f, Cxy = signal.coherence(eegData[i,:,:], eegData[j,:,:], fs=fs, nfft=nfft, axis=0)#, window=np.hanning(nfft))
    out = np.mean(Cxy[np.all([f >= 0.5, f<=4], axis=0)], axis=0)
    return out

##########
# correlation across channels
def PhaseLagIndex(eegData, i, j):
    hxi = ss.hilbert(eegData[i,:,:])
    hxj = ss.hilbert(eegData[j,:,:])
    # calculating the INSTANTANEOUS PHASE
    inst_phasei = np.arctan(np.angle(hxi))
    inst_phasej = np.arctan(np.angle(hxj))

    out = np.abs(np.mean(np.sign(inst_phasej - inst_phasei), axis=0))
    return out

##########
# Cross Correlation
def crossCorrelation(eegData, i, j):
    out = np.zeros(eegData.shape[2])
    for epoch in range(eegData.shape[2]):
        ccor = np.correlate(eegData[i,:,epoch], eegData[j,:,epoch], mode="full")
        absccor = np.abs(ccor)
        out[epoch] = (np.max(absccor) - np.mean(absccor)) / np.std(absccor)
    return out

##########
# Auxilary Cross-correlation Lag
def corrCorrLagAux(eegData,ii,jj,Fs=100):
    out = np.zeros(eegData.shape[2])
    lagCorr = []
    for lag in range(0,eegData.shape[1],int(0.2*Fs)):
        tmp = eegData.copy()
        tmp[jj,:,:] = np.roll(tmp[jj,:,:], lag, axis=0)
        lagCorr.append(CrossCorrelation(tmp, ii, jj, Fs))
    return np.argmax(lagCorr,axis=0)

################################################
#	bandpower Functions
################################################

##########
# compute the bandpower (area under segment (from fband[0] to fband[1] in Hz)
# of curve in freqency domain) of data, at sampling frequency of Fs (100 ussually)
def bandpower(data, fs, fband):
    freqs, powers = periodogram(data, fs)
    idx_min = np.argmax(freqs > fband[0]) - 1
    idx_max = np.argmax(freqs > fband[1]) - 1
    idx_delta = np.zeros(dtype=bool, shape=freqs.shape)
    idx_delta[idx_min:idx_max] = True
    bpower = simps(powers[idx_delta], freqs[idx_delta])
    return bpower

##########
# computes the same thing as vecbandpower but with a loop
def pfvecbandpower(data, fs, fband):
    bpowers = np.zeros((data.shape[0], data.shape[2]))
    for i in range(data.shape[0]):
        freqs, powers = periodogram(data[i, :, :], fs, axis=0)
        idx_min = np.argmax(freqs > fband[0]) - 1
        idx_max = np.argmax(freqs > fband[1]) - 1
        idx_delta = np.zeros(dtype=bool, shape=freqs.shape)
        idx_delta[idx_min:idx_max] = True
       
        bpower = simps(powers[idx_delta, :], freqs[idx_delta], axis=0)
        bpowers[i, :] = bpower

    return bpowers

################################################
#	Complexity features
################################################    

##########
# Extract the Shannon Entropy
# threshold the signal and make it discrete, normalize it and then compute entropy
def shannonEntropy(eegData, bin_min, bin_max, binWidth):
    H = np.zeros((eegData.shape[0], eegData.shape[2]))
    for chan in range(H.shape[0]):
        for epoch in range(H.shape[1]):
            counts, binCenters = np.histogram(eegData[chan,:,epoch], bins=np.arange(bin_min+1, bin_max, binWidth))
            nz = counts > 0
            prob = counts[nz] / np.sum(counts[nz])
            H[chan, epoch] = -np.dot(prob, np.log2(prob/binWidth))
    return H
    
##########
# Extract the tsalis Entropy
def tsalisEntropy(eegData, bin_min, bin_max, binWidth, orders = [1]):
    H = [np.zeros((eegData.shape[0], eegData.shape[2]))]*len(orders)
    for chan in range(H[0].shape[0]):
        for epoch in range(H[0].shape[1]):
            counts, bins = np.histogram(eegData[chan,:,epoch], bins=np.arange(-200+1, 200, 2))
            dist = dit.Distribution([str(bc).zfill(5) for bc in bins[:-1]],counts/sum(counts))
            for ii,order in enumerate(orders):
                H[ii][chan,epoch] = tsallis_entropy(dist,order)
    return H

##########
# Cepstrum Coefficients (n=2)
def mfcc(eegData,fs,order=2):
    H = np.zeros((eegData.shape[0], eegData.shape[2],order))
    for chan in range(H.shape[0]):
        for epoch in range(H.shape[1]):
            H[chan, epoch, : ] = librosa.feature.mfcc(np.asfortranarray(eegData[chan,:,epoch]), sr=fs)[0:order].T
    return H

##########
# Lyapunov exponent
def lyapunov(eegData):
    return np.mean(np.log(np.abs(np.gradient(eegData,axis=1))),axis=1)
    
##########
# Fractal Embedding Dimension
# From pyrem: packadge for sleep scoring from EEG data
# https://github.com/gilestrolab/pyrem/blob/master/src/pyrem/univariate.py
def hFD(a, k_max): #Higuchi FD
    L = []
    x = []
    N = len(a)

    for k in range(1,k_max):
        Lk = 0
        for m in range(0,k):
            #we pregenerate all idxs
            idxs = np.arange(1,int(np.floor((N-m)/k)),dtype=np.int32)
            Lmk = np.sum(np.abs(a[m+idxs*k] - a[m+k*(idxs-1)]))
            Lmk = (Lmk*(N - 1)/(((N - m)/ k)* k)) / k
            Lk += Lmk

        L.append(np.log(Lk/(m+1)))
        x.append([np.log(1.0/ k), 1])

    (p, r1, r2, s)=np.linalg.lstsq(x, L)
    return p[0]
    
##########
# Hjorth Mobility
# Hjorth Complexity
# variance = mean(signal^2) iff mean(signal)=0
# which it is be because I normalized the signal
# Assuming signals have mean 0
# Mobility = sqrt( mean(dx^2) / mean(x^2) )
def hjorthParameters(xV):
    dxV = np.diff(xV, axis=1)
    ddxV = np.diff(dxV, axis=1)

    mx2 = np.mean(np.square(xV), axis=1)
    mdx2 = np.mean(np.square(dxV), axis=1)
    mddx2 = np.mean(np.square(ddxV), axis=1)

    mob = mdx2 / mx2
    complexity = np.sqrt((mddx2 / mdx2) / mob)
    mobility = np.sqrt(mob)

    # PLEASE NOTE that Mohammad did NOT ACTUALLY use hjorth complexity,
    # in the matlab code for hjorth complexity subtraction by mob not division was used 
    return mobility, complexity

##########
# false nearest neighbor descriptor
def falseNearestNeighbor(eegData, fast=True):
    # Average Mutual Information
    # There exist good arguments that if the time delayed mutual
    # information exhibits a marked minimum at a certain value of tex2html_wrap_inline6553,
    # then this is a good candidate for a reasonable time delay.
    npts = 1000   # not sure about this?
    maxdims = 50
    max_delay = 2 # max_delay = 200  # TODO: need to use 200, but also need to speed this up
    distance_thresh = 0.5
    
    out = np.zeros((eegData.shape[0], eegData.shape[2]))
    for chan in range(eegData.shape[0]):
        for epoch in range(eegData.shape[2]):
            if fast:
                out[chan, epoch] = 0
            else:
                cur_eegData = eegData[chan, :, epoch]
                lagidx = 0  # we are looking for the index of the lag that makes the signal maximally uncorrelated to the original
                # # minNMI = 1  # normed_mutual_info is from 1 (perfectly correlated) to 0 (not at all correlated) 
                # # for lag in range(1, max_delay):
                # #     x = cur_eegData[:-lag]
                # #     xlag = cur_eegData[lag:]
                # #     # convert float data into histogram bins
                # #     nbins = int(np.floor(1 + np.log2(len(x)) + 0.5))
                # #     x_discrete = np.histogram(x, bins=nbins)[0]
                # #     xlag_discrete = np.histogram(xlag, bins=nbins)[0]
                # #     cNMI = normed_mutual_info(x_discrete, xlag_discrete)
                # #     if cNMI < minNMI:
                # #         minNMI = cNMI
                # #         lagidx = lag
                # nearest neighbors part
                knn = int(max(2, 6*lagidx))  # heuristic (number of nearest neighbors to look up)
                m = 1 # lagidx + 1

                # y is the embedded version of the signal
                y = np.zeros((maxdims+1, npts))
                for d in range(maxdims+1):
                    tmp = cur_eegData[d*m:d*m + npts]
                    y[d, :tmp.shape[0]] = tmp
                
                nnd = np.ones((npts, maxdims))
                nnz = np.zeros((npts, maxdims))
                
                # see where it tends to settle
                for d in range(1, maxdims):
                    for k in range(0, npts):
                        # get the distances to all points in the window (distance given embedding dimension)
                        dists = []
                        for nextpt in range(1, knn+1):
                            if k+nextpt < npts:
                                dists.append(np.linalg.norm(y[:d, k] - y[:d, k+nextpt]))
                        if len(dists) > 0:
                            minIdx = np.argmin(dists)
                            if dists[minIdx] == 0:
                                dists[minIdx] = 0.0000001  # essentially 0 just silence the error
                            nnd[k, d-1] = dists[minIdx]
                            nnz[k, d-1] = np.abs( y[d+1, k] - y[d+1, minIdx+1+k] )
                # aggregate results
                mindim = np.mean(nnz/nnd > distance_thresh, axis=0) < 0.1
                # get the index of the first occurence of the value true
                # (a 1 in the binary representation of true and false)
                out[chan, epoch] = np.argmax(mindim)
        
    return out 

##########
# ARMA coefficients
def arma(eegData,order=2):
    H = np.zeros((eegData.shape[0], eegData.shape[2],order))
    for chan in range(H.shape[0]):
        for epoch in range(H.shape[1]):
            arma_mod = sm.tsa.ARMA(eegData[chan,:,epoch], order=(order,order))
            arma_res = arma_mod.fit(trend='nc', disp=-1)
            H[chan, epoch, : ] = arma_res.arparams
    return H

################################################
#	Continuity features
################################################  

##########
# median frequency
def medianFreq(eegData,fs):
    H = np.zeros((eegData.shape[0], eegData.shape[2]))
    for chan in range(H.shape[0]):
        freqs, powers = signal.periodogram(eegData[chan, :, :], fs, axis=0)
        H[chan,:] = freqs[np.argsort(powers,axis=0)[len(powers)//2]]
    return H

##########
# calculate band power
def bandPower(eegData, lowcut, highcut, fs):
	eegData_band = filt_data(eegData, lowcut, highcut, fs, order=7)
	freqs, powers = signal.periodogram(eegData_band, fs, axis=1)
	bandPwr = np.mean(powers,axis=1)
	return bandPwr

##########
# numberOfSpikes    
def spikeNum(eegData,minNumSamples=7,stdAway = 3):
    H = np.zeros((eegData.shape[0], eegData.shape[2]))
    for chan in range(H.shape[0]):
        for epoch in range(H.shape[1]):
            mean = np.mean(eegData[chan, :, epoch])
            std = np.std(eegData[chan,:,epoch],axis=1)
            H[chan,epoch] = len(signal.find_peaks(abs(eegData[chan,:,epoch]-mean), 3*std,epoch,width=7)[0])
    return H

##########    
# Standard Deviation
def eegStd(eegData):
	std_res = np.std(eegData,axis=1)
	return std_res

##########
# α/δ Ratio
def eegRatio(eegData,fs):
	# alpha (8–12 Hz)
	eegData_alpha = filt_data(eegData, 8, 12, fs)
	# delta (0.5–4 Hz)
	eegData_delta = filt_data(eegData, 0.5, 4, fs)
	# calculate the power
	powers_alpha = bandPower(eegData, 8, 12, fs)
	powers_delta = bandPower(eegData, 0.5, 4, fs)
	ratio_res = np.sum(powers_alpha,axis=0) / np.sum(powers_delta,axis=0)
	return np.expand_dims(x, axis=0)

###########
# Regularity (burst-suppression)
# Regularity of eeg
# filter with a window of 0.5 seconds to create a nonnegative smooth signal.
# In this technique, we first squared the signal and applied a moving-average
# The window length of the moving average was set at 0.5 seconds.
def eegRegularity(eegData, Fs=100):
    in_x = np.square(eegData)  # square signal
    num_wts = Fs//2  # find the filter length in samples - we want 0.5 seconds.
    q = signal.lfilter(np.ones(num_wts) / num_wts, 1, in_x, axis=1)
    q = -np.sort(-q, axis=1) # descending sort on smooth signal
    N = q.shape[1]
    u2 = np.square(np.arange(1, N+1))
    # COMPUTE THE Regularity
    # dot each 5min epoch with the quadratic data points and then normalize by the size of the dotted things    
    reg = np.sqrt( np.einsum('ijk,j->ik', q, u2) / (np.sum(q, axis=1)*(N**2)/3) )
    return reg

###########
# Voltage < (5μ, 10μ, 20μ)
def eegVoltage(eegData,voltage=20):
	eegFilt = eegData.copy()
	eegFilt[abs(eegFilt) > voltage] = np.nan
	volt_res = np.nanmean(eegFilt,axis=1)
	return volt_res

##########
# Diffuse Slowing
# look for diffuse slowing (bandpower max from frequency domain integral)
# repeated integration of a huge tensor is really expensive
def diffuseSlowing(eegData, Fs=100, fast=True):
    maxBP = np.zeros((eegData.shape[0], eegData.shape[2]))
    idx = np.zeros((eegData.shape[0], eegData.shape[2]))
    if fast:
        return idx
    for j in range(1, Fs//2):
        print("BP", j)
        cbp = bandpower(eegData, Fs, [j-1, j])
        biggerCIdx = cbp > maxBP
        idx[biggerCIdx] = j
        maxBP[biggerCIdx] = cbp[biggerCIdx]
    return (idx < 8)

##########
# Spikes
def spikeNum(eegData,minNumSamples=7,stdAway = 3):
    H = np.zeros((eegData.shape[0], eegData.shape[2]))
    for chan in range(H.shape[0]):
        for epoch in range(H.shape[1]):
            mean = np.mean(eegData[chan, :, epoch])
            std = np.std(eegData[chan,:,epoch])
            H[chan,epoch] = len(signal.find_peaks(abs(eegData[chan,:,epoch]-mean), 3*std,epoch,width=7)[0])
    return H

##########
# Delta Burst after spike
def burstAfterSpike(eegData,eegData_subband,minNumSamples=7,stdAway = 3):
    H = np.zeros((eegData.shape[0], eegData.shape[2]))
    for chan in range(H.shape[0]):
        for epoch in range(H.shape[1]):
            preBurst = 0
            postBurst = 0
            mean = np.mean(eegData[chan, :, epoch])
            std = np.std(eegData[chan,:,epoch])
            idxList = signal.find_peaks(abs(eegData[chan,:,epoch]-mean), stdAway*std,epoch,width=minNumSamples)[0]
            for idx in idxList:
                preBurst += np.mean(eegData_subband[chan,idx-7:idx-1,epoch])
                postBurst += np.mean(eegData_subband[chan,idx+1:idx+7,epoch])
            H[chan,epoch] = postBurst - preBurst
    return H

##########
# Sharp spike
def shortSpikeNum(eegData,minNumSamples=7,stdAway = 3):
    H = np.zeros((eegData.shape[0], eegData.shape[2]))
    for chan in range(H.shape[0]):
        for epoch in range(H.shape[1]):
            mean = np.mean(eegData[chan, :, epoch])
            std = np.std(eegData[chan,:,epoch])
            longSpikes = set(signal.find_peaks(abs(eegData[chan,:,epoch]-mean), 3*std,epoch,width=7)[0])
            shortSpikes = set(signal.find_peaks(abs(eegData[chan,:,epoch]-mean), 3*std,epoch,width=1)[0])
            H[chan,epoch] = len(shortSpikes.difference(longSpikes))
    return H

##########
# Number of Bursts
def numBursts(eegData,fs):
	bursts = []
	supressions = []
	for epoch in range(eegData.shape[2]):
		epochBurst,epochSupressions = burst_supression_detection(eegData[:,:,epoch],fs,suppression_threshold=10)#,low=30,high=49)
		bursts.append(epochBurst)
		supressions.append(epochSupressions)
	# Number of Bursts
	numBursts_res = np.zeros((eegData.shape[0], eegData.shape[2]))
	for chan in range(numBursts_res.shape[0]):
		for epoch in range(numBursts_res.shape[1]):
			numBursts_res[chan,epoch] = len(bursts[epoch][chan])
	return numBursts_res
	
##########
# Burst length μ and σ
def burstLengthStats(eegData,fs):
	bursts = []
	supressions = []
	for epoch in range(eegData.shape[2]):
		epochBurst,epochSupressions = burst_supression_detection(eegData[:,:,epoch],fs,suppression_threshold=10)#,low=30,high=49)
		bursts.append(epochBurst)
		supressions.append(epochSupressions)
	# Number of Bursts
	burstMean_res = np.zeros((eegData.shape[0], eegData.shape[2]))
	burstStd_res = np.zeros((eegData.shape[0], eegData.shape[2]))
	for chan in range(burstMean_res.shape[0]):
		for epoch in range(burstMean_res.shape[1]):
			burstMean_res[chan,epoch] = np.mean([burst[1]-burst[0] for burst in bursts[epoch][chan]])
			burstStd_res[chan,epoch] = np.std([burst[1]-burst[0] for burst in bursts[epoch][chan]])
	burstMean_res = np.nan_to_num(burstMean_res)
	burstStd_res = np.nan_to_num(burstStd_res)
	return burstMean_res,burstStd_res

##########
# Burst band powers (δ, α, θ, β, γ)
def burstBandPowers(eegData, lowcut, highcut, fs, order=7):
	band_burst_powers = np.zeros((eegData.shape[0], eegData.shape[2]))
	bursts = []
	supressions = []
	for epoch in range(eegData.shape[2]):
		epochBurst,epochSupressions = burst_supression_detection(eegData[:,:,epoch],fs,suppression_threshold=10)#,low=30,high=49)
		bursts.append(epochBurst)
		supressions.append(epochSupressions)
	eegData_band = filt_data(eegData, lowcut, highcut, fs, order=7)
	for epoch,epochBursts in enumerate(bursts):
		for chan,chanBursts in enumerate(epochBursts):
			epochPowers = []  
			for burst in chanBursts:
				if burst[1] == eegData.shape[1]:
					burstData =  eegData_band[:,burst[0]:,epoch]
				else:
					burstData =  eegData_band[:,burst[0]:burst[1],epoch]
				freqs, powers = signal.periodogram(burstData, fs, axis=1)
				epochPowers.append(np.mean(powers,axis=1))
			band_burst_powers[chan,epoch] = np.mean(epochPowers)	
	return band_burst_powers

##########
# Number of Suppressions
def numSuppressions(eegData,fs,suppression_threshold=10):
	bursts = []
	supressions = []
	for epoch in range(eegData.shape[2]):
		epochBurst,epochSupressions = burst_supression_detection(eegData[:,:,epoch],fs,suppression_threshold=suppression_threshold)#,low=30,high=49)
		bursts.append(epochBurst)
		supressions.append(epochSupressions)
	numSupprs_res = np.zeros((eegData.shape[0], eegData.shape[2]))
	for chan in range(numSupprs_res.shape[0]):
		for epoch in range(numSupprs_res.shape[1]):
			numSupprs_res[chan,epoch] = len(supressions[epoch][chan])
	return numSupprs_res

##########
# Suppression length μ and σ
def suppressionLengthStats(eegData,fs,suppression_threshold=10):
	bursts = []
	supressions = []
	for epoch in range(eegData.shape[2]):
		epochBurst,epochSupressions = burst_supression_detection(eegData[:,:,epoch],fs,suppression_threshold=suppression_threshold)#,low=30,high=49)
		bursts.append(epochBurst)
		supressions.append(epochSupressions)
	supressionMean_res = np.zeros((eegData.shape[0], eegData.shape[2]))
	supressionStd_res = np.zeros((eegData.shape[0], eegData.shape[2]))
	for chan in range(supressionMean_res.shape[0]):
		for epoch in range(supressionMean_res.shape[1]):
			supressionMean_res[chan,epoch] = np.mean([suppr[1]-suppr[0] for suppr in supressions[epoch][chan]])
			supressionStd_res[chan,epoch] = np.std([suppr[1]-suppr[0] for suppr in supressions[epoch][chan]])
	supressionMean_res = np.nan_to_num(supressionMean_res)
	supressionStd_res = np.nan_to_num(supressionStd_res)
	return supressionMean_res, supressionStd_res

################################################
#	Connectivity features
################################################

##########
# Coherence - δ
def coherence(eegData,fs):
	coh_res = []
	for ii, jj in itertools.combinations(range(eegData.shape[0]), 2):
		coh_res.append(CoherenceDelta(eegData, ii, jj, fs=fs))
	coh_res = np.array(coh_res)
	return coh_res

##########
# Mutual information
def calculate2Chan_MI(eegData,ii,jj,bin_min=-200, bin_max=200, binWidth=2):
    H = np.zeros(eegData.shape[2])
    bins = np.arange(bin_min+1, bin_max, binWidth)
    for epoch in range(eegData.shape[2]):
        c_xy = np.histogram2d(eegData[ii,:,epoch],eegData[jj,:,epoch],bins)[0]
        H[epoch] = mutual_info_score(None, None, contingency=c_xy)
    return H

##########
# Granger causality
def calcGrangerCausality(eegData,ii,jj):
    H = np.zeros(eegData.shape[2])
    for epoch in range(eegData.shape[2]):
        X = np.vstack([eegData[ii,:,epoch],eegData[jj,:,epoch]]).T
        H[epoch] = tsa.stattools.grangercausalitytests(X, 1, addconst=True, verbose=False)[1][0]['ssr_ftest'][0]
    return H

##########
# phase Lag Index
def phaseLagIndex(eegData, i, j):
    hxi = ss.hilbert(eegData[i,:,:])
    hxj = ss.hilbert(eegData[j,:,:])
    # calculating the INSTANTANEOUS PHASE
    inst_phasei = np.arctan(np.angle(hxi))
    inst_phasej = np.arctan(np.angle(hxj))

    out = np.abs(np.mean(np.sign(inst_phasej - inst_phasei), axis=0))
    return out

##########
# Cross-correlation Magnitude
def crossCorrMag(eegData,ii,jj):
	crossCorr_res = []
	for ii, jj in itertools.combinations(range(eegData.shape[0]), 2):
		crossCorr_res.append(crossCorrelation(eegData, ii, jj))
	crossCorr_res = np.array(crossCorr_res)
	return crossCorr_res

##########
# Cross-correlation Lag
def corrCorrLag(eegData,ii,jj,fs=100):
	crossCorrLag_res = []
	for ii, jj in itertools.combinations(range(eegData.shape[0]), 2):
		crossCorrLag_res.append(corrCorrLag(eegData, ii, jj, fs))
	crossCorrLag_res = np.array(crossCorrLag_res)
	return crossCorrLag_res
