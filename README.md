# Mako
Mako is a model-free ultra fast genome structural variants detection tool. It also provides a machine learning based method to classify SVs based on their mutational signature sequential features as well as providing scoring method for complex SVs.

# Install and run

Mako requires Java JDK (>=v1.6) and Python (2.7) to run. 

#### Dependency
* htsjdk(https://github.com/samtools/htsjdk): A Java API for high-throughput sequencing data (HTS) formats.
* Numpy: used in Python script for BAM parameter estimation.

#### Usage
```sh
$ git clone https://github.com/jiadong324/Morphling.git
```

Run bamConfig.py to get BAM statistics, including read length, libraray average and standard deviation of insert size. To run it, you have to specify how much standard deviation (-X) away are considered as discordant insert size and number of read-pairs (-N) you would like to use for the estimation.
```sh
Run bamConfig.py
$ samtools view your.bam | python bamConfig.py -X 3 -N 30000
```
Further you have to create a Mako configuration file, including BAM statistics, absolute path to your BAM file and working directory. 
```
Example of Mako config file
readlen:126
mean:550
stdev:100
bam: /path/to/sample.bam
workDir:/path/to/output
```

Get help info and run SV discovery.
```
$ java -jar /path/to/Morphling/dist/Mako.jar
```

Run mode one: run with BAM file, you can either output all SuperItems to file or keep it in the memory. It is suggested to keep them in file, so you don't need to go through the BAM file agian for next run with different parameters. The command line used to run the program with default parameter settings:
```
$ java -jar /path/to/Mako.jar fa=file.fa bamCfg=bam.cfg
```
Run mode two: run on SuperItem files directly without BAM file

```
$ java -jar /path/to/Mako.jar fa=file.fa itemOut=item.txt
```

#### Output format

The SV output file contains predicted SV position on the genome. Additional information includes SupType, Pattern, Region (genome region spanned by pattern), weights, ratio (allele fraction of each Super-Item), orientation (orientation of reads in Super-Item). A single SV can be supported by more than one evidence, more evidence indicates more confident calls.
* SupType=ARP_Span: indicates SV is combined by two patterns that is able to link together through read-pair. Each pattern of the SV might be a breakpoint. Number of read pairs support such relation is provided.
* SupType=ARP_Self: a pattern is self-linked through read pairs. We estimate potential breakpoint based on abnormal read pairs. Number, quality and weight of these supporting read pairs is provided.
* SupType=Split: indicates SV is discovered based on split alignment. We provide additional information, such as number of split read support, split read mapping quality.
* SupType=Cross: indicates SV is discovered based on local sequence cross links. Additional information includes number of reads support the cross, the maximum cross matched sequence length.
* SupType=Realign: for region with multiple clipped Super-Items, we usually do realignment, this helps discover INDELS and small SVs. Information includes minus and plus strand support read is provided.
* SupType=OEM: one-end-unmapped reads formed cluster may indicate potential insertion breakpoint near OEM Super-Item. This is not a very confident evidence, but we report such abnormal.

# Classification and score CSV

Please go https://github.com/ttbond/makov_classifier to download Mako prediction module. 

# Contact
If you have questions or encouter problems, please feel free to contact: jiadong324@gmail.com, ccxtbut@gmail.com.

License
----



[//]: # (These are reference links used in the body of this note and get stripped out when the markdown processor does its job. There is no need to format nicely because it shouldn't be seen. Thanks SO - http://stackoverflow.com/questions/4823468/store-comments-in-markdown-syntax)


   [dill]: <https://github.com/joemccann/dillinger>
   [git-repo-url]: <https://github.com/joemccann/dillinger.git>
   [john gruber]: <http://daringfireball.net>
   [df1]: <http://daringfireball.net/projects/markdown/>
   [markdown-it]: <https://github.com/markdown-it/markdown-it>
   [Ace Editor]: <http://ace.ajax.org>
   [node.js]: <http://nodejs.org>
   [Twitter Bootstrap]: <http://twitter.github.com/bootstrap/>
   [jQuery]: <http://jquery.com>
   [@tjholowaychuk]: <http://twitter.com/tjholowaychuk>
   [express]: <http://expressjs.com>
   [AngularJS]: <http://angularjs.org>
   [Gulp]: <http://gulpjs.com>

   [PlDb]: <https://github.com/joemccann/dillinger/tree/master/plugins/dropbox/README.md>
   [PlGh]: <https://github.com/joemccann/dillinger/tree/master/plugins/github/README.md>
   [PlGd]: <https://github.com/joemccann/dillinger/tree/master/plugins/googledrive/README.md>
   [PlOd]: <https://github.com/joemccann/dillinger/tree/master/plugins/onedrive/README.md>
   [PlMe]: <https://github.com/joemccann/dillinger/tree/master/plugins/medium/README.md>
   [PlGa]: <https://github.com/RahulHP/dillinger/blob/master/plugins/googleanalytics/README.md>