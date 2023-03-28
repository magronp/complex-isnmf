# Complex ISNMF

Here, you will find the code related to complex ISNMF

If you use any of the things existing in this repository, please cite the [corresponding paper](https://arxiv.org/abs/1802.03156). 

You can also find an online demo with sound examples related to this work on the [companion website](https://magronp.github.io/demos/taslp19_cisnmf.html).

### Data setup

To reproduce the experiments conducted in our paper, you will need to download the [Dexmixing Secret Database (DSD100)](http://www.sisec17.audiolabs-erlangen.de) and to place its content in the `dataset/DSD100` folder.

If you use this dataset, you will end up with the proper directory structure and file names, as used in the `functions/get_data_DSD` function.

If you want to use a different dataset, then you have two options: 
- either you format your file names and directory structure to match the one from DSD100;
- or you modify the file reading function `get_data_DSD` to suit your needs.


### Scripts

The script to reproduce the experiments are placed in the `scripts` folder. They will notably record audio files in the `audio_files` folder, and some metrics (SDR, SIR and SAR) in the `metrics` folder.
