# Thing!
Code that takes the SAM/BAM files output by genome alignment programs and
outputs TCC matrix like those output by [**kallisto**](https://pachterlab.github.io/kallisto/).
So far has only been used on Mac OS X and Linux.

## Installation
Binaries may eventually be provided. For now, they will have to be compiled from
the source code.

### Dependencies

### Making the executable
The code is compiled with the g++ compiler. Install it if you don't have it.

Clone repository:
```
    $ git clone https://github.com/laureneliu/tcc-from-alignment ~/clone/path
```

GFF and SAM/BAM I/O uses parts of the [SeqAn](https://www.seqan.de/) library.
Follow instructions [here](http://seqan.readthedocs.io/en/master/Infrastructure/Use/Install.html)
to download. Or, change into ~/clone/path/tcc-from-alignment and run:

```
    $ wget http://packages.seqan.de/seqan-library/seqan-library-2.4.0.tar.xz
    $ tar -xf seqan-library-2.4.0.tar.xz
    $ mv seqan-library-2.4.0/include ./
```

Depending on where the seqan folder is and what it is named, you may have to go
into the Makefile and change the `SEQAN_PATH` variable at the top of the
Makefile to lead to the seqan folder.

Go to `~/clone/path/tcc-from-alignment`, wherever that may be. Make two new
directories, `obj/` and `bin/` with `mkdir obj/ bin/`.

In the source directory, run `make`. The executable is `~/clone/path/bin/main`.

## Running the program
### Basic workflow
Given some a dataset of RNA-seq reads contained in some FASTQ files:
1. Use a genome alignment program such as [hisat2](https://ccb.jhu.edu/software/hisat2/manual.shtml)
or [STAR](https://github.com/alexdobin/STAR/) to align the reads. Your output
should be a SAM or BAM file.
1. Use this program to read the SAM/BAM file and output the appropriate TCC
matrix, contained in .ec, .tsv, and .cells files.

### Basic command line
```
Usage: ~/clone/path/bin/main [options]* -g <GTF> -S <SAM> [-o <output>]
```

`<GTF>` is the file of annotated sequences to use. Note that it must be in the
[ensembl GTF format](https://uswest.ensembl.org/info/website/upload/gff.html) or
it cannot be read properly. If there are multiple GTFs you wish to use, you may
input them as a comma-separated list of values.

`<SAM>` is the SAM/BAM file output by the genome aligner. Again, it must be in
the [specified format](https://samtools.github.io/hts-specs/SAMv1.pdf). Either
of the two example aligners listed above (hisat2, STAR) should give a
correctly-formated file. As with the GTFs, this option accepts a comma-separated
lis of values.

`<output>` is the name/directory of your output files. Appropriate file
extensions will be added to the name your provide. Default is out.ec, out.tsv,
and out.cells. See below for brief description of what these files look like.

### Other options
* **-p, --threads <int>** The number of threads to use. Defaults to 1.

* **-q** By default, the program is highly verbose and gives
updates on how far it has progressed, along with various warnings. Use this
option to suppress such output. _Note_: this will not suppress any error
messages in case the program must exit.

* **-t, --transcriptome <fa>** The transcriptome file. Only useful with
kallisto. kallisto's output zero-indexes the transcripts. Using this option will
give the output the same indexing. If there are more transcripts in the GTF than
in the transcriptome, the program will continue where kallisto leaves off.

* **-e, --ec <ec>** A kallisto (or other) ec file. Order the TCCs in the same
order as this ec file. Useful for direct comparison of program output with
kallisto pseudo output. However, it might be more useful to look at kallisto
quant output.

* **--full-matrix** Output a full matrix instead of the default sparse matrix.
See below for more information on these formats.

* **-u, --unmatched <SAM>** Also output a SAM file containing all reads that
didn't align to any transcripts. This excludes those reads that didn't align
anywhere on the genome. Currently outputs a bad header.

## File formats
This section will only describe the ec, tsv, and cell files and add some
restrictions to SAM/BAM formatting. For more information on
the [GTF](https://uswest.ensembl.org/info/website/upload/gff.html)
 and [SAM/BAM](https://samtools.github.io/hts-specs/SAMv1.pdf)
file formats, you should look at their documentation.

### .ec files
This file contains two tab-separated columns. The first column is the
zero-indexed row number (the line number). The second contains the equivalence
class. From [the paper on kallisto](https://www.nature.com/articles/nbt.3519):

> An equivalence class for a read is a multi-set of transcripts associated with
> the read; ideally it represents the transcripts a read could have originated
> from and provides a sufficient statistic for quantification.

### .cells files
This is a list of the SAM/BAM files used. They will be listed in the order that
they were input.

### .tsv files
This file may be in either the sparse matrix format (default) or full matrix
format. You may change the output format using the `--full-matrix`` option.

#### Full matrix
This file contains _n + 1_ tab-separated columns, where _n_ is the number of
SAM/BAM files you input. The first column contains the zero-indexed row number
(the line number). The _(i, j)_ th entry, where _j >= 1_, is the number of reads
file _j_ (one-index the SAM/BAM files) that had the equivalence class of row _i_
described in the ec file. A short example:
```
0 0 0 1
1 2 0 1
2 0 1 1
3 3 0 0
```
This tsv files says that three reads in the first SAM/BAM file aligned to the
equivalence class described by row three. Say the ec file says this:
```
0   0
1   1
2   2
3   1,3
```
Then, three reads in the first SAM/BAM file aligned to equivalence class 1,3.

Note that you can look at the .cells file to see which SAM/BAM file this was.

#### Sparse matrix
This file contains three tab-separated columns. Take a full matrix. The first
column of this file gives the zero-indexed row number of an entry. The second
gives the -1-indexed column number (i.e. it one-indexes, but does not count the
first column containing row numbers). The third is the entry itself. It only
shows those entries that are nonzero. Take again the sample tsv file from above:
```
0 0 0 1
1 2 0 1
2 0 1 1
3 3 0 0
```
The sparse matrix looks like this:
```
1   0   2
3   0   3
2   1   1
0   3   1
1   3   1
2   3   1
```

### SAM/BAM restrictions
For paired-end reads, the program will only count as "correct" a pair of
alignments  with flags: 0x01, 0x02, 0x10 for one alignment and 0x20 for the
other, and 0x40 for one alignment and 0x80 for the other.

Also, all paired reads must have at least two alignments, one for each end,
regardless of whether it maps. hisat2 will give a SAM file readable by the
program, but STAR requires the additional option `--outSAMunmapped Within
KeepPairs`. Consult the [manual](https://github.com/alexdobin/STAR/blob/master/doc/STARmanual.pdf)
for more information about the available options.

## Things that still need to be done
- [ ] Better documentation! It's, like, several months out of date.
