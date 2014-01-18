A simple, stateful VCF (http://www.1000genomes.org/node/101) parser for python.

Emphasizes ease of use over performance.

Useful for reading and modifying VCFs.

```python
import sys
from vcfreader import VCFReader

vcfr = VCFReader('/path/to/my_vcf.vcf')

# show which samples are in this file
print ', '.join(vcfr.samples)

while vcfr.nextentry():
    # let me see the genotype of the first sample in this VCF
    gt = vcfr.getgenotype(vcfr.samples[0])
    print ','.join(gt)

    # let me see the AF at this variant
    print vcfr.INFO.get('AF')

    # let's set a custom annotation in the INFO
    vcfr.INFO['ANNOT'] = "HELLOWORLD"

    # now if we write out this variant with outputentry, the update to INFO will be reflected
    vcfr.outputentry()
