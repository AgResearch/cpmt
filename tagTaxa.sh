#!/bin/sh
#---------------------------------------------------------------------------- 
#Based on Tassel ver5.2.38
#https://bitbucket.org/tasseladmin/tassel-5-source/wiki/Tassel5GBSv2Pipeline
#precond:  see the link above for deploying raw data and keyfile
#postcond: tag sequences and counts of each samples
#-----------------------------------------------------------------------------
gbs="~/tools/tasseladmin-tassel-5-standalone-d937dbc2338c/run_pipeline.pl"
proj="cpmt"
keyfile="keyfile_cpmt.txt"
enzyme="ApeKI"
dbname=$proj".db"
output=$proj"TagTaxa.txt"

$gbs -fork1 -GBSSeqToTagDBPlugin -e $enzyme -i ./01_RawSequence/ -k $keyfile -db $dbname -endPlugin -runfork1
$gbs -fork1 -TagExportToFastqPlugin -db $dbname -o tagsForAlign.fq.gz -c 5 -endPlugin -runfork1
$gbs -fork1 -GetTagTaxaDistFromDBPlugin -db $dbname -o output -endPlugin -runfork1
