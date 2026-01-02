/quobyte/bmhenngrp/from-lssc0/projects/NCR_scRNAseq/data/batch*.txt # these are the paths to the 143 that have RNA seq. I have to ensure that its the actual 138 instead since we still use them to demultiplex with vireo as we specif only partial sample genotypes are available

the you run /quobyte/bmhenngrp/from-lssc0/projects/NCR_scRNAseq/scripts/make_batch_vcfs/make_batch_vcfs.sh to make a vcf for each batch then we ude this vcf dor vireo. cellsnp lite has been run already

then i used this  make batch vcf 
thse are the data lists i used in the script above /quobyte/bmhenngrp/from-lssc0/projects/NCR_scRNAseq/data/batch5clean.txt

no cellsnp lite has alrady been run

so time for vireo - RUN THESE /quobyte/bmhenngrp/from-lssc0/projects/NCR_scRNAseq/scripts/demultiplexing/vireo/vireonewb2.sh
