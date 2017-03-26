# ribomethseq
Implements the RiboMethSeq scoring analysis as described in [Birkedal et al.](http://onlinelibrary.wiley.com/doi/10.1002/ange.201408362/full)


## Installation
Tested in Python 3.4.

First, make sure pip is up to date:

`pip3 install pip --upgrade`

Then just install directly from github:

`pip3 install https://github.com/caethan/ribomethseq/zipball/master`

Alternatively, you can clone the repo and install locally:

`python setup.py install`

`pip install -r requirements.txt`

## Use
Installation will add a script called `calculate-ribomethseq-from-bam` into your path.  Run it as follows:

`calculate-ribomethseq-from-bam --input-bam-file YOUR_PAIRED_READ_FILE.bam --output-dir YOUR_OUTPUT_DIRECTORY`

The current version works _only on paired-end reads_ - it will not work on data from single-end sequencing experiments.
This is something that could be added in the future if there is demand.

The output directory will end up with a large number of Wiggle-formatted .wig files and a single .csv file.

There will be a .wig file for each chromosome in the input data with:

* 5' end counts
* 3' end counts
* total end counts (sum of previous two)   
* scores A, B, and C, as described in the Birkedal paper
    
These files are suitable for loading in to IGV, for example, while the .csv file is suitable for loading into Excel.
