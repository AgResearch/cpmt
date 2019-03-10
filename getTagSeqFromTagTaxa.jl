#------------------------------------------------------------
#extract tag sequences and save as fasta files 
#precond:  tagTaxa text file -- output file from tagTaxa.sh
#postcond: tag sequences saved in individual files  
#tested in Julia v1.1.0
#------------------------------------------------------------
using JuliaDB, BioSequences
 
@time d=loadtable("cpmt_TagTaxa.txt"; delim='\t')
@time save(d, "./cpmt_tagtaxa.juliadb")
tagTaxa_cpmt=load("cpmt_tagtaxa.juliadb")
snames=colnames(tagTaxa_cpmt)             
snames=snames[2:end]                     

#;mkdir ./tmp
#;cd tmp

tic=time()
for sidx in 1:length(snames)
  tagcount1=select(tagTaxa_cpmt, snames[sidx])
  sum(tagcount1)
  k=findall(a->a!=0, tagcount1)
  length(k)
  countPerTag=tagcount1[k]
  tag1=select(tagTaxa_cpmt, :Tag)[k]
  #describe(countPerTag)

  #save to fasta file
  outfile=string(snames[sidx], ".fa")
  writer = FASTA.Writer(open(outfile, "w"))
  for i in 1:length(countPerTag)
    n=countPerTag[i]
    if(n==1)
      rec = FASTA.Record(string(i), DNASequence(tag1[i]))
      write(writer, rec)
    else
      #more than one counts, no zeros as filtered
      for j in 1:n
        id=string(i)*"-"*string(j)
        rec = FASTA.Record(id, DNASequence(tag1[i]))
        write(writer, rec)
    end 
    end
   end
   close(writer)
end                                
time()-tic




