# bam2tcc
Code that takes the SAM/BAM files output by genome alignment programs and
outputs TCC matrix like those output by [**kallisto**](https://pachterlab.github.io/kallisto/).
So far has only been used/tested on Mac OS X and Linux.

## Installation.
Binaries may eventually be provided. For now, they will have to be compiled from
the source code. To compile `bam2tcc` you will need a c++14 compatible compiler such as gcc 5.0 or Clang 3.4 or later

### Dependencies.
[CMake](https://cmake.org/) for compiling and [SeqAn](https://www.seqan.de/) for
SAM/BAM I/O. To download SeqAn, follow instructions
[here](http://seqan.readthedocs.io/en/master/Infrastructure/Use/Install.html).

For compiling `bam2tcc` is is enough to download the [library package](http://packages.seqan.de/seqan-library/seqan-library-2.4.0.zip) and unzip to a directory.

### Making the executable.
Clone repository:
```
    $ git clone https://github.com/laureneliu/bam2tcc ~/clone/path
```

Make a new directory called `build` and change into it:
```
    $ mkdir build && cd build
```

Run CMake (`cmake ..`) from the `build` directory. 
Depending on where/how you installed SeqAn, you may need to append some extra
options, `-DCMAKE_PREFIX_PATH and -DCMAKE_INCLUDE_PATH` so that CMake can find
the relevant SeqAn folder. See [here](https://seqan.readthedocs.io/en/master/Infrastructure/Use/FindSeqAnCMake.html#install-seqan-into-user-defined-prefix-or-clone-from-github). If you downloaded the SeqAn library it should be run with

```cmake -DCMAKE_PREFIX_PATH=download-path/seqan-library-2.4.0/share/cmake     -DSEQAN_INCLUDE_PATH=download-path/seqan-library-2.4.0/include ..```

The executable is `build/src/bam2tcc`.

## Running the program.
### Basic workflow.
Given a set of RNA-seq reads contained in some FASTQ files:
1. Use a genome alignment program such as [HISAT2](https://ccb.jhu.edu/software/hisat2/manual.shtml)
or [STAR](https://github.com/alexdobin/STAR/) to align the reads. Your output
should be a SAM or BAM file.
1. Use samtools to order the output by genomic coordinate if necessary. The
command is `samtools sort -T [tempPrefix] -@ [nthreads] -o [outfile] infile`.
[Here](http://www.htslib.org/doc/samtools.html) for more information.
1. Use this program to read the SAM/BAM file and output the appropriate TCC
matrix, contained in .ec, .tsv, and .cells files.

### Basic command line.
```
Usage: /path/bam2tcc/build/src/bam2tcc [options]* -g <GFF> -S <SAM/BA> [-o <output>]
```

`<GFF>` is the file of annotated sequences to use. Note that it must be in the
[ensembl GFF format](https://uswest.ensembl.org/info/website/upload/gff.html) or
it cannot be read properly. If there are multiple GFFs you wish to use, you may
input them as a comma-separated list of values. Entries should be grouped by
chromosome, then by transcript. The entry of type "transcript" must come before
the entries of type "exon" for one given transcript. Exons should be in order of
increasing genomic start coordinate if the transcript is on the forward strand,
and in order of decreasing genomic start coordinate if the transcript is on the
reverse strand. The program will check the GFF for this format before using it.
GFFs downloaded from ensembl should be in the correct format.

`<SAM/BAM>` is the SAM/BAM file output by the genome aligner. Again, it must be in
the [specified format](https://samtools.github.io/hts-specs/SAMv1.pdf). Either
of the two example aligners listed above (HISAT2, STAR) should give a
correctly-formated file. As with the GTFs, this option accepts a comma-separated
list of values. Alignments should be sorted by chromosome, then by genomic start
coordinate. samtools provides a way to sort alignments in this way.

`<output>` is the name/directory of your output files. Appropriate file
extensions will be added to the name your provide. Default is matrix.ec,
matrix.tsv, and matrix.cells. See below for brief description of what these
files look like.

### Other options.
* **-U** Indicate that the reads are unpaired. Only necessary for running on
kallisto --genomebam output.i

* **-k** Indicate that all input SAM/BAM files are kallisto genomebam files.
Required if input file is in the BAM format, or if it does not specify the
@PG:ID as kallisto in the header.

* **-R** Indicate that all input SAM/BAM files are RapMap `lightweight` files.
Required if input file is in the BAM format, or if it does not speciy the @PG:ID
as RapMap in the header.

* **-t, --transcriptome <fa>** The transcriptome file. Only useful with
kallisto. kallisto's output zero-indexes the transcripts. Using this option will
give the output the same indexing. If there are more transcripts in the GTF than
in the transcriptome, the program will continue where kallisto leaves off.

* **-e, --ec <ec>** A kallisto (or other) ec file. Order the TCCs in the same
order as this ec file. Useful for direct comparison of program output with
kallisto pseudo output. However, it might be more useful to look at kallisto
quant output.

* **-p, --threads <int>** The number of threads to use. Defaults to 1.

* **--full-matrix** Output a full matrix instead of the default sparse matrix.
See below for more information on these formats.

* **-u, --unmatched <SAM>** Also output a SAM file containing all reads that
didn't align to any transcripts. This excludes those reads that didn't align
anywhere on the genome. Currently outputs a bad header.

* **--check-gff** Only check GFF format.

### Alternative compilation options.
There exists a preprocessor macro, READ_DIST, that gives the option to output
the names of those reads that map to some transcript (i.e., whose TCCs are
nonempty). The macros is in Mapper.hpp in the src folder.
kallisto --genomebam gives information on how specific reads map, so it is
possible to determine on which reads kallisto and genomic aligners differ.
In the kallisto-generated file, all reads where at least one alignment does not
have the "unmapped flag" (0x04) are considered to have mapped.

## File formats.
This section will only describe the ec, tsv, and cell files and add some
restrictions to SAM/BAM formatting. For more information on
the [GFF](https://uswest.ensembl.org/info/website/upload/gff.html)
 and [SAM/BAM](https://samtools.github.io/hts-specs/SAMv1.pdf)
file formats, you should look at their documentation.

### .ec files.
This file contains two tab-separated columns. The first column is the
zero-indexed row number (the line number). The second contains the equivalence
class. From [the paper on kallisto](https://www.nature.com/articles/nbt.3519):

> An equivalence class for a read is a multi-set of transcripts associated with
> the read; ideally it represents the transcripts a read could have originated
> from and provides a sufficient statistic for quantification.

### .cells files.
This is a list of the SAM/BAM files used. They will be listed in the order that
they were input.

### .tsv files.
This file may be in either the sparse matrix format (default) or full matrix
format. You may change the output format using the `--full-matrix`` option.

#### Full matrix.
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

#### Sparse matrix.
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

### SAM/BAM restrictions.
For paired-end reads, the program will only count as "correct" a pair of
alignments  with flags: 0x01, 0x02, 0x10 for one alignment and 0x20 for the
other, and 0x40 for one alignment and 0x80 for the other.

Also, all paired reads must have at least two alignments, one for each end,
regardless of whether it maps. HISAT2 will give a SAM file readable by the
program, but STAR requires the additional option `--outSAMunmapped Within
KeepPairs`. Consult the [manual](https://github.com/alexdobin/STAR/blob/master/doc/STARmanual.pdf)
for more information about the available options.

## Various poorly-document utilities.
There exists another executable, `/path/bam2tcc/build/src/debug`. It provides various useful
functions, but does not go into detail on what they do. Also, any mistakes
in input format of either the files or of the commandline arguments will result
in a potentially generic error (e.g. a segfault). However, none of the functions
are particularly complicated so a brief look at the source code should be all
the explanation necessary. Run the executable with no arguments for brief and
utterly uninformative descriptions of all functions and their arguments.

Of particular use are:
* "sparse" (`-x`), which turns a full matrix from ONE INPUT
FILE into a sparse matrix

* "Unmapped cat" (`-o`), which takes in an input
SAM (not a BAM), presumably the file containing unmapped reads output by
`bam2tcc --unmapped`, and determines which reads were not properly aligned to
the genome, and therefore did not map to any transcripts, and which reads did
properly align somewhere in the genome, but in an area which did not correspond
to any transcripts in the transcriptome. `sameQName` should be 1 if all segments
of a read are named the same way (as opposed to being named, for example, read/1
and read/2 to indicate that the first is the first mate, and the second the
second mate). `genomebam` should be 1 if the input SAM alignments were
originally generated by kallisto --genomebam.

* "EQ to TCC" (`-t`), which takes as input the eq_classes.txt file generated by
Salmon --dumpEQ and outputs a TCC in kallisto's format.
