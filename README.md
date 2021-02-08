# EEGExtract
[Sari Saba-Sadiya](https://cse.msu.edu/~sadiyasa/)<sup>1,2</sup>,
[Eric Chantland]()<sup>1</sup>,
[Taushang Liu](https://npal.psy.msu.edu/)<sup>1</sup>,
[Tuka Alhanai](https://talhanai.xyz/)<sup>3</sup>,
[Mohammad Ghassemi](https://ghassemi.xyz/)<sup>1</sup><br>
<sub>
<sup>1</sup> Human Augmentation and Artificial Intelligence lab, Michigan State University, Department of Computer Science<br>
<sup>2</sup> Neuroimaging of Perception and Attention Lab, Michigan State University, Department of Psychology<br>
<sup>3</sup> Computer Human Intelligence Lab, New York University Abu Dhabi, Department of Electrical and Computer Engineering<br>
</sub>


A pyhton package for extracting EEG features. First developed for the paper ["Unsupervised EEG Artifact Detection and Correction"](https://www.frontiersin.org/articles/10.3389/fdgth.2020.608920/abstract), published in Frontiers in Digital Health, Special issue on Machine Learning in Clinical Decision-Making.


## License

<a href="https://www.gnu.org/licenses/gpl-3.0.en.html">
<img src="https://cdn2.iconfinder.com/data/icons/business-and-finance-311/32/Business_and_Finance_scale_weight_weights_scales_balance-256.png" alt="GPL-3.0 License" height="16" style="vertical-align:middle"> GPL-3.0 License </img></a>

Free to use and modify, but must cite the original publication below.

## The Features
| Signal Descriptor                       | Brief Description|
| --------------- | --------------- |
| __Complexity Features__                                 | degree of randomness or irregularity  |
| Shannon Entropy        | additive measure of signal stochasticity    |
| Tsalis Entropy (n=10)  | non-additive measure of signal stochasticity  |
| Information Quantity  (δ,α,θ,β,γ)  | entropy of a wavelet decomposed signal         |
| Cepstrum Coefficients (n=2)                        | rate of change in signal spectral band power   |
| Lyapunov Exponent                                | separation between signals with similar trajectories   |
| Fractal Embedding Dimension                      | how signal properties change with scale |
| Hjorth Mobility                                  | mean signal frequency   |
| Hjorth Complexity                                | rate of change in mean signal frequency   |
| False Nearest Neighbor                           | signal continuity and smoothness |
| ARMA Coefficients (n=2)                       | autoregressive coefficient of signal at (t-1) and (t-2)  |
| __Continuity Features__                               | clinically grounded signal characteristics |
| Median Frequency                               |   the median spectral power    |
| δ band Power                             |  spectral power in the 0-3Hz range  |
| α band Power                             |  spectral power in the 4-7Hz range   |
| θ band Power                             |  spectral power in the 8-15Hz range   |
| β band Power                              |  spectral power in the 16-31Hz range   |
| γ band Power                             |  spectral power above 32Hz |
| Median Frequency                               |   median spectral power   |
| Standard Deviation                              | average difference between signal value and it's mean   |
| α/δ Ratio                          | ratio of the power spectral density in α and δ bands |
| Regularity (burst-suppression)                 | measure of signal stationarity / spectral consistency  |
| Voltage < (5μ, 10μ, 20μ)          |  low signal amplitude |
| Normal EEG                                |        Peak spectral power textgreater= 8Hz   |
| Diffuse Slowing                           |        indicator of peak power spectral density less than 8Hz   |
| Spikes                                    |        signal amplitude exceeds μ by 3σ for 70 ms or less   |
| Delta Burst after spike                   |        Increased δ after spike, relative to δ before spike |
| Sharp spike                               |        spikes lasting less than 70 ms   |
| Number of Bursts                          |        number of amplitude bursts  |
| Burst length μ and σ             |        statistical properties of bursts |
| Burst band powers (δ,α,θ,β,γ)   | spectral power of bursts  |
| Number of Suppressions                            | segments with contiguous amplitude suppression  |
| Suppression length μ and σ             | statistical properties of suppressions |
| __Connectivity Features__                             |             interactions between EEG electrode pairs  |
| Coherence - δ                            | correlation in in 0-4 Hz power between signals    |
| Coherence - All                                 | correlation in overall power between signals |
| Mutual Information                               | measure of dependence |
| Granger causality - All                          | measure of causality |
| Phase Lag Index                                  | association between the instantaneous phase of signals |
| Cross-correlation Magnitude                      | maximum correlation between two signals |
| Crosscorrelation - Lag                           | time-delay that maximizes correlation between signals |

## Important Note
the feature extractor is an independent section that can be used with any artifact correction method (recently there have been quite a few including some notable example [1,2]). If you are interested in the specific setup that was used in the paper, as well as a link to the data, please visit [the following repository](https://github.com/sari-saba-sadiya/EEG-Artifact-Correction-Via-Completion).

[1] S. Phadikar, N. Sinha, and R. Ghosh, “Automatic eeg eyeblink artefactidentification and removal technique using independent componentanalysis in combination with support vector machines and denoisingautoencoder”

[2] B. Somers, T. Francart, and A. Bertrand, “A generic eeg artifactremoval algorithm based on the multi-channel wiener filter.”



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

