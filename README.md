# EEGExtract
[Sari Saba-Sadiya](https://cse.msu.edu/~sadiyasa/)<sup>1,2</sup>,
[Eric Chantland]()<sup>2</sup>,
[Taosheng Liu](https://npal.psy.msu.edu/)<sup>2</sup>,
[Tuka Alhanai](https://talhanai.xyz/)<sup>3</sup>,
[Mohammad Ghassemi](https://ghassemi.xyz/)<sup>1</sup><br>
<sub>
<sup>1</sup> Human Augmentation and Artificial Intelligence lab, Michigan State University, Department of Computer Science<br>
<sup>2</sup> Neuroimaging of Perception and Attention Lab, Michigan State University, Department of Psychology<br>
<sup>3</sup> Computer Human Intelligence Lab, New York University Abu Dhabi, Department of Electrical and Computer Engineering<br>
</sub>


A python package for extracting EEG features. First developed for the paper ["Unsupervised EEG Artifact Detection and Correction"](https://www.frontiersin.org/articles/10.3389/fdgth.2020.608920/abstract), published in Frontiers in Digital Health, Special issue on Machine Learning in Clinical Decision-Making. [Press here](https://www.frontiersin.org/articles/10.3389/fdgth.2020.608920/bibTex) for a BibTex citation (or scroll to the bottom of this page).

To the best of our knowledge EEGExtract is the most comprehensive library for EEG feature extraction currently available. This library is actively maintained, __please open an issue if you believe adding a specific feature will be of benefit for the community!__


## Setup
1. Make sure that you have the required packages listed in `requirements.txt`. Use `pip install -r requirements.txt` if unsure. 
2. Simply download and place the `EEGExtract.py` file in the same folder as your repo. You can then use `import EEGExtract as eeg`.

## License

<a href="https://www.gnu.org/licenses/gpl-3.0.en.html">
<img src="https://cdn2.iconfinder.com/data/icons/business-and-finance-311/32/Business_and_Finance_scale_weight_weights_scales_balance-256.png" alt="GPL-3.0 License" height="16" style="vertical-align:middle"> GPL-3.0 License </img></a>

Free to use and modify, but must cite the original publication below.

## The Features

| _Signal Descriptor_                       | _Brief Description_ | _Function_ |
| --------------- | --------------- |  --------------- |
| __Complexity Features__                                 | _degree of randomness or irregularity_  |
| Shannon Entropy        | additive measure of signal stochasticity    | shannonEntropy |
| Tsalis Entropy (n=10)  | non-additive measure of signal stochasticity  | tsalisEntropy |
| Information Quantity  (δ,α,θ,β,γ)  | entropy of a wavelet decomposed signal         | filt_data + shannonEntropy |
| Cepstrum Coefficients (n=2)                        | rate of change in signal spectral band power   | mfcc |
| Lyapunov Exponent                                | separation between signals with similar trajectories   | lyapunov|
| Fractal Embedding Dimension                      | how signal properties change with scale | hFD |
| Hjorth Mobility                                  | mean signal frequency   | hjorthParameters |
| Hjorth Complexity                                | rate of change in mean signal frequency   | hjorthParameters |
| False Nearest Neighbor                           | signal continuity and smoothness | falseNearestNeighbor |
| ARMA Coefficients (n=2)                       | autoregressive coefficient of signal at (t-1) and (t-2)  | arma |
| __Continuity Features__                               | _clinically grounded signal characteristics_  |
| Median Frequency                               |   the median spectral power    | medianFreq |
| δ band Power                             |  spectral power in the 0-3Hz range  | bandPower |
| α band Power                             |  spectral power in the 4-7Hz range   | bandPower |
| θ band Power                             |  spectral power in the 8-15Hz range   | bandPower |
| β band Power                              |  spectral power in the 16-31Hz range   | bandPower |
| γ band Power                             |  spectral power above 32Hz | bandPower |
| Median Frequency                               |   median spectral power   | medianFreq |
| Standard Deviation                              | average difference between signal value and it's mean   | eegStd |
| α/δ Ratio                          | ratio of the power spectral density in α and δ bands | eegRatio |
| Regularity (burst-suppression)                 | measure of signal stationarity / spectral consistency  | eegRegularity |
| Voltage < (5μ, 10μ, 20μ)          |  low signal amplitude | eegVoltage |
| Normal EEG                                |        Peak spectral power textgreater= 8Hz   | |
| Diffuse Slowing                           |        indicator of peak power spectral density less than 8Hz   | diffuseSlowing |
| Spikes                                    |        signal amplitude exceeds μ by 3σ for 70 ms or less   | spikeNum |
| Delta Burst after spike                   |        Increased δ after spike, relative to δ before spike | burstAfterSpike |
| Sharp spike                               |        spikes lasting less than 70 ms   | shortSpikeNum |
| Number of Bursts                          |        number of amplitude bursts  | numBursts |
| Burst length μ and σ             |        statistical properties of bursts | burstLengthStats |
| Burst band powers (δ,α,θ,β,γ)   | spectral power of bursts  | burstBandPowers |
| Number of Suppressions                            | segments with contiguous amplitude suppression  | numSuppressions |
| Suppression length μ and σ             | statistical properties of suppressions | suppressionLengthStats |
| __Connectivity Features__                             |  _interactions between EEG electrode pairs_  |
| Coherence - δ                            | correlation in in 0-4 Hz power between signals    | filt_data + coherence |
| Coherence - All                                 | correlation in overall power between signals | coherence |
| Mutual Information                               | measure of dependence | calculate2Chan_MI |
| Granger causality - All                          | measure of causality | calcGrangerCausality |
| Phase Lag Index                                  | association between the instantaneous phase of signals | phaseLagIndex |
| Cross-correlation Magnitude                      | maximum correlation between two signals | crossCorrMag |
| Crosscorrelation - Lag                           | time-delay that maximizes correlation between signals | corrCorrLagAux |

Additionally, `EEGExtract` also contains implementations for a number of auxiliary functions

| _Function_  | params |  _Brief Description_ |
| --------------- | --------------- |  --------------- |
| filt_data | eegData, lowcut, highcut, fs, order=7 | midpass filter between lowcut and highcut |
| fcnRemoveShortEvents | z,n| z=[chan x samples ], n is threshold |
| get_intervals | A,B,endIdx | Find interval of consistent values in binary 1D numpy array |

## Important Note
The feature extractor is an independent section that can be used with any artifact correction method (recently there have been quite a few including some notable example [1,2]). If you are interested in the specific setup that was used in the paper, as well as a link to the data, please visit [the following repository](https://github.com/sari-saba-sadiya/EEG-Artifact-Correction-Via-Completion).

[1] S. Phadikar, N. Sinha, and R. Ghosh, “Automatic eeg eyeblink artifact identification and removal technique using independent component analysis in combination with support vector machines and denoising autoencoder”

[2] B. Somers, T. Francart, and A. Bertrand, “A generic eeg artifact removal algorithm based on the multi-channel wiener filter.”



## Cite
```
@article{saba2020unsupervised,
  title={Unsupervised EEG Artifact Detection and Correction},
  author={Saba-Sadiya, Sari and Chantland, Eric and Alhanai, Tuka and Liu, Taosheng and Ghassemi, Mohammad Mahdi},
  journal={Frontiers in Digital Health},
  volume={2},
  pages={57},
  year={2020},
  publisher={Frontiers}
}
```

