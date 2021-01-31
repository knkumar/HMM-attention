- [HMM for describing attention from cursor movements](#hmm-for-describing-attention-from-cursor-movements)
  * [Data](#data)
  * [Code](#code)
    + [Transform data](#transform-data)
    + [Generate fits](#generate-fits)
    + [Fix fits based on operational constraints of attention](#fix-fits-based-on-operational-constraints-of-attention)
    + [Summarize fits to observe an individual's attention](#summarize-fits-to-observe-an-individual-s-attention)

# HMM for describing attention from cursor movements

The code here can be used to quantify aspects of human attention in an experimental paradigm described [here](https://cogsci.mindmodeling.org/2018/papers/0139/0139.pdf). This experimental paradigm is an improvement from an [earlier paradigm](https://www.jstor.org/stable/10.5406/amerjpsyc.128.2.0253?seq=1) to obtain a continuous measure of attention. A continuous measure of attention has the benefit of assessing an individual's selective and sustained attention. More importantly, the vast amounts of data provided by the paradigm allows sophisticated analysis to assess individual differences. 

This repository describes the use of HMM modeling in 4 stages:

1. Transform data for use in HMM (in folder pre)
2. Generate fits from HMM (in folder model)
3. Fix fits from HMM based on knowledge of attention (in folder fix)
4. Summarize fits to generate necessary observations about an individual's attention. 

## Data

The *data* for one subject is stored in folder **data**, available in file *S4-05-18-2017.mat*. 

The data from the [experimental paradigm](https://cogsci.mindmodeling.org/2018/papers/0139/0139.pdf) is organized in a structure as below:

block_data

- Session **One**
  - Block **Red1**
    - mouse
      - Trial R31 (cursor data for trial used in the model)
      - ...
  - Block **Red2**
  - ...
- Session **Two**
  - Same structure as session 1

The naming convention for the trial is as follows:

1. R31 - Red target presented at position 3 on first trial

   >  `R <Target R:Red S:Size Q:square>`
>
   > `3 <Position on screen>`
>
   > `1 <Trial number>`

2. RS25 - Size foil presented when searching for a red target at position 2 on fifth trial

   >  `R <Target R:Red S:Size Q:square>`
>
   > `S <Foil R:Red S:Size Q:square>`
>
   > `2 <Position on screen>`
>
   > `5 <Trial number>`

## Code

The main driver program is *att_hmm.m*.

The following line saves a fit mat file in the same folder as the data by using the ```fitmodel``` function.

```matlab
% sName: directory structure with filename for data
% Resolution: screen resolution for experimental paradigm
% iterations: number of iterations for optimization
fitmodel(sName, Resolution, iterations);
```

### Transform data

The first stage in fitting a model involved the preparation of data. It uses the following files from the *pre* folder.

- `getDataAnalysis.m`

  The function starts by loading the law data *block_data* and generates *mouse_data* used to calculate `hand` and `targets` data structures.

  - `hand` structure contains data for each period, block, trial and XY coordinates.
  - `targets` structure contains state and target flag.

- `generateMouseData.m`

  The function uses raw data to calculate normalized mouse data (from *normalizeData.m*) based on the modified resolution to downsample raw data. 

### Generate fits

- `optimHMM.m`

  The function is used to find optimal parameters for the first and last 10000 data points using fminsearch.

- `analyze10.m`

### Fix fits based on operational constraints of attention

### Summarize fits to observe an individual's attention

