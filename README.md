## Centripetal propagation of calcium signals in astrocytes

Custom code to compute **centripetal propagation** patterns (delay maps) from **astrocytic calcium** recordings using pixel-wise **correlation functions**.

<p align="center"><img src="https://github.com/HelmchenLabSoftware/Centripetal_propagation_astrocytes/blob/main/Github%20centripetal%20propagation.png"  width="65%"></p>

This image illustrates mean fluorescence (left) and computed delay maps (right). Delay maps are organized such that positive delays are close to astrocytic somata and negative delays far from somata (middle). These delay maps illustrate that astrocytic calcium signals propagate from distal processes to somatic compartments for individual astrocytes.

### How the algorithm works

The code computes the delay with respect to a reference time trace for each pixel of a 3D movie. This enables to extract typical spatio-temporal delays with respect to a mean time trace across the movie. It works best for either long or relatively noise-free movies. The details of the algorithm are described in Fig. 6 and the associated supplementary figures in **[this preprint](https://www.biorxiv.org/content/10.1101/2022.08.16.504030v1.full)**.

### How to use the algorithm

- The code includes extensive comments and should be mostly self-explanatory.

- The code is not optimized for performance in order to keep it as simple as possible.

- The scripts `Demo_analysis_Matlab.m` and `Demo_analysis_Python.py` use a calcium imaging recording as a tif-file and produce a delay map as shown here.

- The code is provided as scripts in *Matlab* (tested with Matlab R2020b) and *Python 3* (tested with various Python 3 versions) and with an example dataset (small 80x80 FOV excerpt covering a single astrocyte). The Python version requires the standard packages `numpy`, `glob`, `tifffile`, `skimage`, and `matplotlib` for visualization.

- The small FOV excerpt was chosen to minimize the size of the Github repository. It is possible to run the code on much larger files.

<p align="center"><img src="https://github.com/HelmchenLabSoftware/Centripetal_propagation_astrocytes/blob/main/Example_astrocyte_delay_map_Github.png"  width="55%"></p>

### When will centripetal propagation patterns be visible?

- Centripetal propagation requires either noise-less or rather long recordings. We denoised our recordings and used 10-20 min recordings. Patterns might however also be visible in short recordings, depending on the recorded calcium activity and behaviors.

- Centripetal propagation requires behavioral events that lead to arousal, e.g., running episodes or negative valence-stimuli like air puffs.

- Centripetal propagation requires neuromodulation as control factor. Centripetal propagation will normally not be observed under anesthesia or in slices.

- We do not know whether astrocytes in brain regions other than hippocampus and cortex will also exhibit centripetal propagation.

### How to cite this work

If you use these scripts in its original or in a  modified version, please cite the following publication as reference: 

> Rupprecht, P., Lewis, C. M., & Helmchen, F. *Centripetal integration of past events by hippocampal astrocytes.* bioRxiv, 2022-08 (2022). [https://www.biorxiv.org/content/10.1101/2022.08.16.504030v1](https://www.biorxiv.org/content/10.1101/2022.08.16.504030v1)

### Questions?

Any questions should be addressed to [Peter Rupprecht](mailto:p.t.r.rupprecht+centripetal+propagation@gmail.com). Alternatively, file an issue on our [Github page](https://github.com/HelmchenLabSoftware/Centripetal_propagation_astrocytes/issues).
