
# Filter reads
Set of tools to select alignments based on:
- length
- first nucleotide of sequence read

This is particularly useful for the analysis of _C. elegans_ small RNA analysis which have fairly specific properties. 

The main script is `filterReads/filterSmallRNAclasses.py`. Another script, `filterReads/summarizeNucleotideByReadLenght.py` will determine the frequency of each nucleotide per read length from any bam file. 

Use case and tutorial to come.

## Dependencies

- pysam (working with 0.8.1)


```bash
rm data/a.sam
for b in /home/adomingu/imb-kettinggr/adomingues/projects/imb_ketting_2014_14_almeida_smallRNA_celegans/results/small_RNA_classes/N2_adult_r1_untreated.*.bam; do
    echo $b
    samtools view $b | head >> data/a.sam
done

samtools view -H $b > data/header.sam
cat data/header.sam data/a.sam | samtools sort - | samtools view -bS - > data/a.bam
```

This creates a file in which we have exactly:

| length | base | count |
|:-------|:-----|------:|
| 21     | T    | 10    |
| 22     | G    | 10    |
| 26     | G    | 10    |


Let's see if our script is working:

## Testing

```bash
cd imb-kettinggr/adomingues/projects/filterReads/

genes="/fsimb/groups/imb-kettinggr/genomes/Caenorhabditis_elegans/Ensembl/WBcel235/Annotation/Genes/genes.gtf"

python filterReads/filterSmallRNAclasses.py -i data/a.bam -o stdout -c 21U | intersectBed -abam stdin -b $genes | samtools view 

python filterReads/summarizeNucleotideByReadLenght.py -i data/a.bam -o res.out
cat res.out 
```

> 26      G       10
> 21      T       10
> 22      G       10

And would this work with a proper file?

```bash
bsub -q short -n1 -app Reserve5G -J summary python filterReads/summarizeNucleotideByReadLenght.py -i /home/adomingu/imb-kettinggr/adomingues/projects/imb_ketting_2014_14_almeida_smallRNA_celegans/results/mapped/multimapped/N2_adult_r1_untreated.bam -o res.out
cat res.out 
```
> 21      T       490832
> 23      T       318399
> 22      T       289332
> 24      T       102213
> 22      G       71496
> 23      A       61258
> 26      G       55746
> 22      A       51805
> 20      T       46093
> 24      C       45768
> 23      C       39615

And does this match what we have seen previously?

```bash
samtools view -c /home/adomingu/imb-kettinggr/adomingues/projects/imb_ketting_2014_14_almeida_smallRNA_celegans/results/small_RNA_classes/N2_adult_r1_untreated.21U.bam
samtools view -c /home/adomingu/imb-kettinggr/adomingues/projects/imb_ketting_2014_14_almeida_smallRNA_celegans/results/small_RNA_classes/N2_adult_r1_untreated.22G.bam
```
>490832
>71496

Yes it does. We have a working script. It is also very memory efficient for these small files:

> ------------------------------------------------------------
> # LSBATCH: User input
> python filterReads/summarizeNucleotideByReadLenght.py -i /home/adomingu/imb-kettinggr/adomingues/projects/imb_ketting_2014_14_almeida_smallRNA_celegans/results/mapped/multimapped/N2_adult_r1_untreated.bam -o res.out
> ------------------------------------------------------------
> 
> Successfully completed.
> 
> Resource usage summary:
> 
>     CPU time :               21.08 sec.
>     Max Memory :             2.99 MB
>     Average Memory :         2.99 MB
>     Total Requested Memory : 5120.00 MB
>     Delta Memory :           5117.01 MB
>     (Delta: the difference between total requested memory and actual max usage.)
>     Max Swap :               84 MB
> 
>     Max Processes :          1
>     Max Threads :            1


## filter with bed file

Create example bed which should contain the `21U` and `22G` reads:

```bash
echo -e "I\t1\t1000
I\t3500\t4000
I\t4500\t5000" > data/a.bed
cat data/a.bed
```


## Some errors

Writing to `stdout` was working but reporting:

>close failed in file object destructor:
sys.excepthook is missing
lost sys.stderr

Turns it was an issue with writing to `stderr`. Solved with the function:

```python
def eprint(*args, **kwargs):
    """
    helper function to print to stderr
    source:
    https://stackoverflow.com/questions/5574702/how-to-print-to-stderr-in-python
    """
    print(*args, file=sys.stderr, **kwargs)
```

and replacing `print` with `eprint`:

```python
    eprint("Number of %s reads: %d" % (args.RNAclass, piRNA_reads))
    eprint("Number of reads filtered out: %d" % other_reads)
```


From discussions with René, we should be more flexible with the definition of small RNA classes and play less importance to the first nucleotide - there could be perfectly good 22G without the G. What we should do is to filter according to:
1. length range
2. intersect with relevant features, for instance, sense to 21U genes, or antisense to annotated genes

We can then plot the length distribution of reads along with the first nucleotide bias. Let's test:

```bash
cd imb-kettinggr/adomingues/projects/filterReads/

python filterReads/filterSmallRNAclasses.py -i data/a.bam -o stdout -m 21 -M 26 --nuc T | samtools view -c

```
