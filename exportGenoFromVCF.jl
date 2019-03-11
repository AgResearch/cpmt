#-------------------------------------------
#export genotype data from VCF
#precond:  VCF file from bcftools
#postcond: export SNP data in a text file
#tested in Julia v1.1.0
#-------------------------------------------
using GeneticVariation, DelimitedFiles

function exportGenoFromVCF(vcffile)
  fh=open(vcffile, "r")
  vcfReader=VCF.Reader(fh)
  sampleID=header(vcfReader).sampleID
  nsample=length(sampleID)
  G=fill("", nsample)
  snpID=String[]           
  qual=Float64[]          

  seekstart(fh)
  vcfReader=VCF.Reader(fh)
  for rec in vcfReader                            
   snp1=VCF.genotype(rec, 1:nsample, "GT")
   G=hcat(G, snp1)    
   push!(snpID, string(VCF.chrom(rec), "_", VCF.pos(rec)))
   push!(qual, VCF.qual(rec))
  end
  close(vcfReader)

  #prepare genotype data matrix to export
  G=G[:,2:end]
  nSNP=size(G,2)
  dataHeader=vcat("sample", snpID)
  res=hcat(sampleID, G)
  tmp=reshape(dataHeader, 1, length(dataHeader))
  res2=vcat(tmp, res)

  outfile=replace(vcffile, ".vcf"=>".txt")
  writedlm(outfile, res2)
  println("genotype data saved as: ", abspath(outfile))
end #>>>

vcffile="pop2_cpmt.raw.vcf"
@time exportGenoFromVCF(vcffile)



