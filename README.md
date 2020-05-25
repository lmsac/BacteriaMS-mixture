# BacteriaMS-mixture
MALDI-MS based bacterial mixture identification. 


## Requirement
* [R] (https://www.r-project.org/). As an alternative, the latest version of [Microsoft R Open] (https://mran.microsoft.com/open/) should be fine.

* [RStudio] (https://www.rstudio.com/) is recommended but optional.


## Database
BacteriaMS-mixture uses peak lists for bacteria identification. A peak list is stored in a space-separated text (.txt) file including the *m/z* and *intensity* of all the peaks from a mass spectrum of a strain of bacteria with the first column as *m/z* and the second column as *intensity*, e.g.:

```
4376.517617 6369.468
4545.36698 14923.287008
4616.504007 3381.736378
4778.424564 6757.753407
4830.992604 10117.37512
……
9899.339424 4964.88939
10159.44407 1130.018489
10323.81304 5016.852179
10707.91308 1477.481354
11240.44102 1693.14288
```

Some example peak lists are provided in this repo.


## Usage
If you have everything installed, you can run identification for a sample spectrum as follows:

1. Run [main.R] (main.R) to load the functions.

2. Load the reference peak lists.
```R
setwd(REF_PATH) 
# REF_PATH is the path of the directory that contains the reference peak list files.

reference.peaklists = local({
  filenames = list.files(pattern = '.txt')
  peaklists = lapply(filenames, read.table)
  names(peaklists) = sub('.txt', '', filenames)
  lapply(peaklists, normalize.intensity)
})
```

3. Load the sample peak list.
```R
sample.peaklist = read.table(SAMPLE_PATH)
# SAMPLE_PATH is the path of the sample peak list file.
sample.peaklist = normalize.intensity(sample.peaklist)
```


4. Search database.
```R
result = search.datebase.mixture(
  sample.peaklist,
  reference.peaklists,
  tolerance = 2000, # Unit: ppm
  mix.number = 4 # Max component count
)

result[['cosine']][1, ] # Print the top match.
```
See also: [identify.R](scripts/identify.R).


5. Search database with jackknife assessment.
```R
result = search.datebase.mixture.jackknife(  
  sample.peaklist,
  reference.peaklists,
  tolerance = 2000, # Unit: ppm
  mix.number = 4, # Max component count
  jackknife.ratio = 1/3, 
  jackknife.number = 100
)

conf = get.confidence.score(result)
conf # Print the species ranked by confidence score.
```
See also: [jackknife.R](scripts/jackknife.R) for more information.


Notes
See [preprocess.R](scripts/preprocess.R) for baseline subtraction and peak extraction.


## Publications
1. Yang, Y., Lin, Y., Chen, Z., Gong, T., Yang, P., Girault, H., Liu, B., Qiao, L. Bacterial whole cell typing by mass spectra pattern matching with bootstrapping assessment. *Anal Chem* **89**, 12556–12561 (2017). https://doi.org/10.1021/acs.analchem.7b03820.
2. Yang, Y., Lin, Y., Qiao, L. Direct MALDI-TOF MS identification of bacterial mixtures. *Anal Chem* **90**, 10400–10408 (2018). https://doi.org/10.1021/acs.analchem.8b02258.

## License
BacteriaMS-mixture is distributed under a BSD license. See the LICENSE file for details.


## Contacts
Please report any problems directly to the github issue tracker. Also, you can send feedback to liang_qiao@fudan.edu.cn  
