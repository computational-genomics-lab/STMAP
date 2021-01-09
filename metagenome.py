#!/usr/bin/env python3
import luigi
import sys

from tasks.metagenome.readCleaning.rawReadQC import readqc
from tasks.metagenome.readCleaning.rawReadQC import rawReadsQC
from tasks.metagenome.readCleaning.preProcessReads import cleanFastq
from tasks.metagenome.readCleaning.preProcessReads import cleanReads
from tasks.metagenome.readCleaning.preProcessReads import filtlong


from tasks.metagenome.utility.luigiconfig import configureProject

from tasks.metagenome import metaAssembly
from tasks.metagenome import genome_binning
from tasks.metagenome import bin_refinement
from tasks.metagenome import genome_dereplicate
from tasks.metagenome import genome_annotation
from tasks.metagenome import genome_enrichment
from tasks.metagenome import metagenome_profiling


if __name__ == '__main__':

    luigi.run()
