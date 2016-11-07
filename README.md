# Aligner

Qaligner is an ultrafast short reads mapper. It indexes reference sequences using k-mer hash tables, and quickly searches mapping results using an efficient alignment algorithm. 

1. Hardware requirements:

	64bit CPU
	20G memory

2. Software environment:

	linux operating system, zlib library version 1.2.8

3. Input Data Format:

	FASTQ, read length 20-500, optimal for reads between 50-150bp. Right now, Qalinger can only handle sequencing data from ILLUMINA or CPaS. ABI proton or Pacbio data are not acceptable. 

4. Reference size limit:

	Right now, the total length of reference seqences should not exceed 4G. 

5. complie
	a.  g++ indexbuilder.cpp -o indexbuilder （reference indexing）
	b.  g++ quickalignerpe.cpp -o paligner -lpthread -lz（pair ends aligner）
	c.  g++ quickalignerse.cpp -o saligner -lpthread -lz（single end aligner）

5. Usage

	Index:
		./indexbuilder	reference.fa index 
	For example,
		 ./indexbuilder /home/ref/hg19.fa /home/index/hg19

	/home/index directory should be created before indexing.

	Mapping:

	a) pair reads 

	./paligner index reads1.fq reads2.fq -n nthreads -fr 0 -id sample_id -sm sample_name -pl platform_name -lb library_name > results.sam

	-n, number of threads

	-fr, direction of reads, the default value is 0, which means two arms of reads are belong to different directions, one being forward, the other being reverse. When reads are belong to the same direction, this parameter should be set to 1.

	-id, id of the sample

	-sm, name of the sample

	-pl, name of plateform, ILLUMINA or CPaS

	-lb, name of the library.

	-id, -sm, -pl, -lb are parameters for read group header line, the same as -R parameter of BWA.

	For example,

	./paligner /home/index/hg19 /home/data/reads1.fq /home/data/reads2.fq -n 40 -fr 0 -id normal -sm human -pl ILLUMINA -lb SRR54768 > results.sam

	b)single end

	./saligner index reads1.fq -n nthreads -id sample_id -sm sample_name -pl platform_name -lb library_name > results.sam

	All parameters have the same meaning as those of paligner. 

	For exameple,

	./saligner /home/index/hg19 /home/data/reads1.fq -n 40 -id normal -sm human -pl ILLUMINA -lb SRR54768 > results.sam

6. Results 

	Unique or multple mapping: when mapping score is equal to 1, the reads can be mapped to multiple positions.If the score is larger than 1, the mapping is unique.
