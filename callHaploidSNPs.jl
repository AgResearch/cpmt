#-------------------------------------------------
#Tag alignemt and snp calling
#precond:  fasta files -- tag sequences
#postcond: SNP output in vcf
#Julia 1.1.0; samtools v1.9
#-------------------------------------------------
faFiles = filter(a->endswith(a, ".fa"), readdir())
snames = map(a->replace(a, ".fa"=>""), faFiles)
bamFiles=snames.*".bam"
REF="~/ref/LpPlastids.fa"

tic=time()
for i in 1:length(snames)
  id=string(i)
  sname1=snames[i]
  fa1=faFiles[i]
  bam1=bamFiles[i]
  align=pipeline(`bwa mem -t 5 -R "@RG\tID:$id\tSM:$sname1" $REF $fa1`,
               `  samtools view -hbS`,
               `  samtools sort --threads 2 -o $bam1`)

  run(align)
end
time()-tic                

#indexing
for j in 1:length(bamFiles)
  bam1=bamFiles[j]
  run(`samtools index -@ 5 $bam1`)
end

#OR switch and work on shell
;ls -1 *.bam > bamfileList        
;REF="~/ref/LpPlastids.fa"
;bcftools mpileup -C 2 -f $REF -b bamfileList > cpmt.mpileup  
;bcftools call --ploidy 1 -vc cpmt.mpileup > cpmt.raw.vcf
