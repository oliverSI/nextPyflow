# nextPyflow
NextPyflow is a Python module that helps you build complex bioinformatics workflow.

**Note: This project is still in alpha stage.**

## Example
Example pipeline for mapping sequence files against a reference genome.
```
import nextPyflow
from nextPyflow.utils import localfile

class bwa_index(nextPyflow.Task):
    def parameter(self, ref_fasta):
        self.docker = 'kathrinklee/bwa:latest'
    def requires(self):
        return [localfile(self.ref_fasta)]
    def run(self):
        return '''                                                                                                                  
        <%! import os %>                                                                                                            
        bwa index ${os.path.basename(self_.ref_fasta)}                                                                              
        '''
        
class bwa(nextPyflow.Task):
    def parameter(self, ref_fasta, fastq1, fastq2):
        self.docker = 'kathrinklee/bwa:latest'
    def requires(self):
        return [localfile(self.ref_fasta),
                localfile(self.fastq1),
                localfile(self.fastq2),
                bwa_index(self.ref_fasta)]
    def run(self):
        return '''                                                                                                                  
        <%! import os %>                                                                                                            
        bwa mem -M \                                                                                                                
        ${os.path.basename(self_.ref_fasta)} \                                                                                      
        ${os.path.basename(self_.fastq1)} \                                                                                         
        ${os.path.basename(self_.fastq2)} \                                                                                         
        > aln-pe.sam                                                                                                                
        '''
        
class bam_sort(nextPyflow.Task):
    def parameter(self, ref_fasta, fastq1, fastq2):
        self.docker = 'comics/samtools:latest'
    def requires(self):
        return [bwa(self.ref_fasta, self.fastq1, self.fastq2)]
    def run(self):
        return '''                                                                                                                  
        samtools sort \                                                                                                             
        -o aln-pe-sorted.bam \                                                                                                      
        -O BAM \                                                                                                                    
        aln-pe.sam                                                                                                                  
        samtools index aln-pe-sorted.bam                                                                                            
        '''
```

Pipelines should be run like this
```
workflow = nextPyflow.Workflow(max_core=1)
workflow.run(bam_sort('/path/to/ref_fasta', '/path/to/fastq1', '/path/to/fastq2'))
```
