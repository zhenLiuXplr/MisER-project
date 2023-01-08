# MisSER-project
An annotation based method to find and fix small exons missed alignment defects in Nanopore long reads.

<img src="examples/pictures/ex1.png" width=800>
The example plot above show the misaligned exon in the middle is find by compaired the reference annotation BED and realigned to the correct position. When we use [minimap2](https://github.com/lh3/minimap2) to align ONT cDNA reads, it is easy to miss the exon with small size, because of the difficulty to find an exact match anchor in these region.

### INSTALLATION

`pip install MisSER`

Upgrade to a newer version using:
`pip install MisSER --upgrade`

The package is written for python3

### INPUT

MisSER requires four essential arguments.
- BAM file
- genome reference fasta (with index fai in the same directory)
- transcript annotation in BED12 format
- the target file for the output exon-missed regions information 

### OUTPUT
- regions which have missed small exons
- realigned BAM file which are fixed 

### USAGE
```
MisSER [-h] [-v] [-c N] [-s N] [-d float] [-f N] [--strandSpecific]
              [--allTranscripts] [--fixFlankLen] [--debugMode] [--setTag]
              [-o file | --onlyRegion]
              inBam genomeFasta annotBed outRegion

positional arguments:
  inBam                 Input original bam file.
  genomeFasta           Reference genome fasta file, with fai index under the same directory.
  annotBed              Annotated transcripts file in BED12 format.
  outRegion             Output Region file, regions contain missed small exons.

optional arguments:
  -h, --help            show this help message and exit
  -v, --version         Print version and exit.
  -c, --coreNum N       The number of cpu cores we used. (default: 1)
  -s, --exonSizeThd N   The threshold of exons size, ignore exons with size > exonSizeThd. (default: 80)
  -d, --deltaRatioThd float
                        The threshold of absolute delta ratio, ignore abs(delta ratio) > deltaRatioThd.
                        (default: 0.5)
  -f, --flankLen N      The extended length on the both sides of realign region. (default: 20)
  --strandSpecific      Only compare reads consistent with annotated strand. (default: False)
  --allTranscripts      Return all possible missed exons among overlapping transcripts. (default: False)
  --fixFlankLen         Flank length will not extend to cover the adjacent indels. (default: False)
  --debugMode           Won't stop when meet an error in one read. (default: False)
  --setTag              Set fr tags on the realigned reads (default: False)
  -o, --outBam file     Output realigned bam file. (default: None)
  --onlyRegion          Only return the Region file without realign process. (default: False)

```

### NOTES
<img src="examples/pictures/illustrator.png" width=400>
Misaligned exons are absent from annotated position and connected to the neighbour exon as extra protruding parts with high error rate. MisSER will first select all introns on reads which overlap with annotated exons and set borders on the suspected misaligned regions. Then MisSER tries to realign the read sequence in region and compare the alignment score before and after realignment. If alignment score improves, the region will be assigned as misaligned region.

- **Delta length** As the picture show above, *Delta length* is defined as the extra bases between reads flank exons and reference flank exons. Intuitively, it is the bases count in the extra protruding parts. We find the reads introns if the intron overlaps with the annotated exon, and then check whether the *Delta length* close to the *exon length*. It is easy to think that *Delta length* should be close to the *Exon length* if they come from the annotated exon.
- **Flank length** In order to define the range of the realign region. we set the region start = min(annoted splice site, read splice site) - *Flank length* and the region end = max(annoted splice site, read splice site) + *Flank length*
- **Delta ratio** *Delta ratio* = (*Delta length* - *Exon length*) / *Exon length*. We choose *Exon length* <= `--exonSizeThd` and *Delta ratio* <= `--deltaRatioThd` to filter out the false positives and record the regions as suspected exon-missed region. The regions information is output in `outRegion`.
- `--allTranscripts` One Read intron may be overlapped with exons on multiple transcripts. As default, we only return the annotated transcript with minimum *Delta raio* as the correct reference. You can set `--allTranscripts` to keep all possible transcripts reference. In the realign process, we only keep one realign result among different transcript, the one with highest realignment score will be kept.

We use [parasail](https://github.com/jeffdaily/parasail-python) to do the global pairwise alignment between read sequence and annotated exon sequence in the realign region. We only keep the realign result if realigned score > original score. 
ATTENTION: the original score does not refer to the AS field in BAM if provided. We calculate the realigned score and origial score based on the realignment and original alignment in the realign region, and the score is equal to the matched bases count - edit distance. As different alignment tools may have different score system, we do not change the AS of NM field in BAM if provided.
- `--onlyRegion` Although we default to return the realigned result in `-o, --outBam file`, you can set the `--onlyRegion` to skip the realign process (although the realign process is not the bottleneck at present).
- `--strandSpecific` This argument is used if the original reads is stranded RNA-seq. We will only try the strand of mapping read to find the overlapped exons.

### EXAMPLE USAGE
```bash
git clone https://github.com/zhenLiuExplr/fixalign-project
cd examples
MisSER ex.bam ex.fa ex_annotation.bed ex_realign_region -o realn.bam
```
Be careful that chromosome in annotBed, genomeFasta and inBam should have same naming style (All in UCSC style like "chr1" or in Ensembl style like "1"). Inconsistent naming style will lead to failed judgement.

### ACKNOWLEDGMENTS/CONTRIBUTORS
- Zhen Liu for building and maintance
- Wu Wei and Chenchen Zhu for advising

### CONTRIBUTING
Welcome for all suggestions, bug reports, feature request and contributions. You can leave an [issue](https://github.com/zhenLiuExplr/MisSER-project/issues) or open a pull request.

### CITATION
If you use this tool, please consider citing our [publication]()