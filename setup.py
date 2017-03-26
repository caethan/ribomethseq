from distutils.core import setup

setup(
    name='RiboMethSeq',
    version='0.1.0',
    author='Brett N. Olsen',
    author_email='brett.olsen@gmail.com',
    packages=['ribomethseq'],
    scripts=['bin/calculate-ribomethseq-from-bam', 'bin/bed-to-wig', 'bin/split-bed'],
    license='LICENSE',
    description='Implementation of the RiboMethSeq analysis algorithms',
    requires=[
        "numpy",
        "pysam",
    ],
)
