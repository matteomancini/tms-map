# tms-map - Significant electrical sources mapping of TMS-EEG data

This code processes the electrical sources estimated from electroencephalography (EEG) data acquired during transcranial magnetic stimulation (TMS) to obtain a map of the significant sources, compute descriptive measures of their distribution and scattering, and finally compute the perturbation complexity index as a measure of complexity.

The general approach and the related measures are described in the following journal publications:

- Casali et al. (2010), "General indices to characterize the electrical response of the cerebral cortex to TMS", NeuroImage [link](https://www.sciencedirect.com/science/article/pii/S1053811909010052?via%3Dihub);
- Casali et al. (2013), "A Theoretically Based Index of Consciousness Independent of Sensory Processing and Behavior", Sci Transl Med [link](https://www.science.org/doi/10.1126/scitranslmed.3006294);

## Setup

After cloning the repository, you'll need a recent Python interpreted (tested on `3.12.1`) and to install the required packages listed in `requirements.txt`:

    $ pip install -r requirements.txt

To use the C implementation of the Lempel-Ziv complexity index, you'll need to compile the related source file:

	$ cc -fPIC -shared -o lzc.so lzc.c

## Use

The use of the computing and plotting functions is illustrated in the Jupyter Notebook `processing.ipynb`. The data required are the source-reconstructed data that can be estimated with software packages such as [BrainStorm](https://neuroimage.usc.edu/brainstorm) and [MNE-Python](https://mne.tools/stable/index.html).

## Acknowledgements

The development of this code was supported by a Starting Grant awarded in 2022 from the Italian National Institute of Health to Matteo Mancini. 

## License

This code is licensed under the MIT license, which you can find in
the `MIT-LICENSE.txt` file.
