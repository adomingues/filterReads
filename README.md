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

```bash
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

