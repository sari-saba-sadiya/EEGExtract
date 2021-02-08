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

