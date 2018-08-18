# #!/bin/bash
#
# FILES=/data/khandekara2/cancer_WGBS/raw_data/*_WGBS.bed
# SUFFIX=sorted
# for f in $FILES
# do
#    if [[$f == MALY_tumor_*]];
#    then
#    bedtools sort -i $f > temp.bed
#    mv temp.bed $f$SUFFIX
#    fi
#
# done
import os
import pybedtools
for file in os.listdir('/data/khandekara2/cancer_WGBS/raw_data'):
    if file.startswith('MALY_tumor_') and file.endswith('_WGBS.bed'):
        a = pybedtools.BedTool(file)
        a.sort().saveas(file + '.sorted')
import os
import pybedtools
for file in os.listdir('/data/khandekara2/imputation'):
    if file.endswith('.sorted'):
        a = pybedtools.BedTool(file)
        b = pybedtools.BedTool('/data/khandekara2/bed_CpGs/cds_canonical.bed')
        a.intersect(b).saveas(file + '.cds')
