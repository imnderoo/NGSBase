*.bed can be obtained from UCSC browser. But, make sure to replace coordinate with manual custom coordinates.
*.genelist can be entered manually
*.exoncoords can be obtained from UCSC browser using selected columns from the refGene table. The last column, cds.name, however needs to be manually added by querying the NM number on NCBI

Note. 

1. Manifest used by MiSeqReporter is 1-based. In-House GATK BED File should be 0-based.
    What this means is that the start position that is -1 compared to MiSeqReporter manifests. The MiSeqReporter should be +-15 from 

