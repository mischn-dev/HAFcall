module HAFcall


import DataFrames, Statistics, CSV, ProgressMeter, Crayons, Gadfly, Compose, StatsBase

using DataFrames, Statistics, CSV, ProgressMeter, Crayons, Gadfly, Compose, StatsBase

export  writeDict, readDict, AFcall, locate, haplotyping, contig, stacker, annotate


"""
SAVE created files from the HAFcall pipeline \n
\n
This function should help to write the dictionary generated from Allelefreqcalling() to multiple text seperated files in a folder of choice \n
csv or text file as format can be choosen by "text" or "csv"
# Example
writeDict("C:/Users/me/Desktop/", file, "csv")
"""
function writeDict(disk, ram, how="text")
                             if how == "text"
                                           for i in keys(ram)
                                                            println(i)
                                                            if typeof(ram[i]) == DataFrame
                                                            CSV.write(string(disk,string(i),".txt"),ram[i];delim = '\t')
                                                            else
                                                              writedlm(string(disk,string(i),".txt"),ram[i], '\t')
                                                            end
                                            end
                             elseif how == "csv"
                                     for i in keys(ram)
                                                      println(i)
                                                      if typeof(ram[i]) == DataFrame
                                                      CSV.write(string(disk,string(i),".csv"),ram[i];delim = ',')
                                                    else
                                                      writedlm(string(disk,string(i),".txt"),ram[i], ',')
                                                     end
                                     end

                             end
 end

"""
READ files created by HAFcall pipline into julia environment \n
\n
Text can be read that has been saved by "writeDict()" in the exact same shape as is has been created by "Allelefreqcalling()" \n
The function also allows only to read in certain files identified by the second input - "A" will read all files containing an "A" in there name \n

# Example \n

readDict("C:/Users/me/Desktop/folder", "txt")
"""
function readDict(loc, ok = "txt")
                        p = Dict()
                        if ok != "std"
                                            for i in readdir(loc)[occursin.(ok,readdir(loc))]
                                                                 println(i)
                                                                 r = string(loc,"/", i)
                                                                 r2 = CSV.read(r, DataFrame)

                                                                 p[Symbol(split(i,".")[1])] = r2
                                            end
                                            return p
                        else
                                            for i in readdir(loc)[.!(occursin.(r"f[0-9]", readdir(loc)))]
                                                                 println(i)
                                                                 r = string(loc,"/", i)
                                                                 r2 = CSV.read(r, DataFrame)

                                                                 p[Symbol(split(i,".")[1])] = r2
                                            end
                                            return p

                        end
        end


"""
ALLELE calling from an VCF file \n
\n
This is the starting Function where the SNP information is extracted from the VCF file created by a varaint caller \n
it works fine for both SNPs and INDELs and can handle infinite numbers of Samples in one file
An OVERVIEW file is required as well , where the name and some other informations are listed by which the files can be merged / linked / compared \n
this has to fit to the order of the VCF file !! \n
the style of the OVERVIEW file should follow this style: \n
Name   \t    Env    \t   Generation \n
E1  \t        0  \t        0      \t            # this is the first parent \n
E2    \t      0    \t      0       \t          # this is the second parent \n
F3env1R1 \t   1     \t     3        \t         # followed by the progenies of different environments and generations ( the Name can be chosen as personal preference) \n
F3env1R2  \t  1    \t      3 \n
F3env2R1  \t  2    \t      3 \n
F3env2R2  \t  2   \t       3 \n

 - other settings to make in the Allelefreqcalling function: \n
 minreaddepth => minimum read depth - default 1 \n
 minqual => minimum quality - default 30 \n
 maxqual => maximum quality - default 1000 \n
 posE1 => the row in the OVERVIEW file that the first parent is written \n
 posE2 => the row in the OVERVIEW file that the second parent is written \n

AFcall(vcffile, genolist, minreaddepth=1, minqual=30, maxqual=1000,posE1=1, posE2=2) \n

# Example \n
overviewfile = "info.txt" \n
freq  = AFcall(SNPfile.vcf, overviewfile, 3,30,1000,1,2)
"""
function AFcall(vcffile, genolist, minreaddepth=1, minqual=30, maxqual=1000,posE1=1, posE2=2)

                                                                                                              vcf2 = @time CSV.read(vcffile,DataFrame; delim="\t", comment="#", header=false)
                                                                                                              glist = CSV.read(genolist,DataFrame; delim="\t")
                                                                                                              nam = glist[!,1]
                                                                                                              pl = ["Chr", "Pos","id", "Ref", "Alt", "Qual", "Filter", "Info", "Format"]
                                                                                                              # set the name of the snp file
                                                                                                              append!(pl, nam)
                                                                                                              rename!(vcf2,Symbol.(pl))


                                                                                                            println("Data is loaded successfully")
                                                                                                            frequency = Dict()

                                                                                                            ## the file for linking to the individual genotype dicts is processed according to the needs of CuArrays - make all strings become integers
                                                                                                            frequency[:General] = Matrix(vcf2[!,[1,2,4,5,6]])

                                                                                                            # transfer the information to a seperate library for each genotype
                                                                                                            # split the individual colum for each library into 3 columns giving info about the call type, REf readcount and ALT readcount
                                                                                                            Threads.@threads  for j in nam
                                                                                                                         p = vcf2[!,Symbol(j)]
                                                                                                                         k = String[]
                                                                                                                         m = Int64[]
                                                                                                                         n = Int64[]

                                                                                                                        for i in 1:size(p,1)
                                                                                                                                             push!(k, string(split(p[i],":")[1]))
                                                                                                                                             push!(m,parse(Int64,string(split(split(p[i],":")[3],",")[1])))
                                                                                                                                             push!(n,parse(Int64,string(split(split(p[i],":")[3],",")[2])))
                                                                                                                         end

                                                                                                                         frequency[Symbol(j)] = hcat(k,m,n)

                                                                                                            end
                                                                                                             println("seperation in Dictionaries completed")

                                                                                                                            ############## experimental step
                                                                                                                            # this step should introduce an correction of the snp calls generated - it is reported that the reference base has a higher change to be alligend than the alt base has
                                                                                                                            # therefore a correction will be introduced where the over all Loc calculated freq is estimated and the corrected to 0.5
                                                                                                                            println("normalisation of the frequency calls is calculated")
                                                                                                                            re = repeat([0.5],inner=length(nam)-2)
                                                                                                                          Threads.@threads  for i in 3:length(nam)

                                                                                                                                                 re[i-2] = sum(frequency[Symbol(nam[i])][:,3])/sum(frequency[Symbol(nam[i])][:,2] .+ frequency[Symbol(nam[i])][:,3])


                                                                                                                                               end
                                                                                                                                               frequency[:correction_alt] = vcat(fill(0,2),re)


                                                                                                    # remove homogen snps in between the parents
                                                                                                                                                println("Process filtering..")
                                                                                                                                                # highlight the snps where the parents have low coverage - convert these to a "./."
                                                                                                                          Threads.@threads     for n in nam[[posE1, posE2]]
                                                                                                                                                                            # calulate the readdepth
                                                                                                                                                                            p = frequency[Symbol(n)][:,2] .+ frequency[Symbol(n)][:,3]
                                                                                                                                                                            frequency[Symbol(n)] = hcat(frequency[Symbol(n)],p)
                                                                                                                                                                            frequency[Symbol(n)] = Array(frequency[Symbol(n)])

                                                                                                                                                                            for i in 1:size(frequency[Symbol(n)],1)
                                                                                                                                                                                                                 if frequency[Symbol(n)][i,4] < minreaddepth
                                                                                                                                                                                                                 frequency[Symbol(n)][i, 1] = "./."
                                                                                                                                                                                                                 end
                                                                                                                                                                             end
                                                                                                                                               end
                                                                                                                                                # check which snps are equal in between the parents and HETEROZYGOT - label them to remove these from ALL Dicts
                                                                                                                                                sel = (.&(frequency[Symbol(nam[posE1])][:,1] .!= "./.", frequency[Symbol(nam[posE1])][:,1] .!= frequency[Symbol(nam[posE2])][:,1], frequency[Symbol(nam[posE2])][:,1] .!= "./.", frequency[Symbol(nam[posE1])][:,1] .!= "0/1", frequency[Symbol(nam[posE2])][:,1] .!= "0/1"))
                                                                                                                         Threads.@threads     for j in nam; frequency[Symbol(j)] = frequency[Symbol(j)][sel,:]; end

                                                                                                                                                frequency[:General] = frequency[:General][sel,:] # the same for the General file
                                                                                                                                                sel=nothing


                                                                                                                                            # filter the Quality of the SNP calls - everything lower qual will be removed
                                                                                                                          Threads.@threads  for i in nam
                                                                                                                                                                        frequency[Symbol(i)] = frequency[Symbol(i)][.&(frequency[:General][:,5] .> minqual, frequency[:General][:,5] .< maxqual),:]
                                                                                                                                                         end

                                                                                                                                                frequency[:General] = frequency[:General][.&(frequency[:General][:,5] .> minqual, frequency[:General][:,5] .< maxqual),:] # the same for the General file


                                                                                                                                              #### get rid of the snps with more than two alleles
                                                                                                                                               # create a new dict to store information from all dict allel calls in one table
                                                                                                                                               frequency[:Alleles] = frequency[:General][:,[1,2]]
                                                                                                                                               for i in 1:size(nam,1); frequency[:Alleles] = hcat(frequency[:Alleles],frequency[Symbol(nam[i])][:,1]); end

                                                                                                                                             frequency[:Alleles] = frequency[:Alleles][:,3:size(frequency[:Alleles],2)]
                                                                                                                                             # save information of the number of different alleles in one value per column
                                                                                                                                             kind = Bool[]
                                                                                                                                             for i in 1:size(frequency[:Alleles],1); append!(kind,all(.&(frequency[:Alleles][i,:] .!= "1/2", frequency[:Alleles][i,:] .!= "0/2", frequency[:Alleles][i,:] .!= "2/2", frequency[:Alleles][i,:] .!= "2/3", frequency[:Alleles][i,:] .!= "1/3", frequency[:Alleles][i,:] .!= "0/3"))
                                                                                                                                                     ) ; end
                                                                                                                                             frequency[:Alleles] = kind

                                                                                                                    Threads.@threads        for i in nam
                                                                                                                                                         frequency[Symbol(i)] = frequency[Symbol(i)][frequency[:Alleles],:]
                                                                                                                                               end
                                                                                                                                             # remove the information also for the "General" file
                                                                                                                                             frequency[:General] = frequency[:General][frequency[:Alleles],:]

                                                                                                                                             delete!(frequency, :Alleles)
                                                                                                                                             println("Loci with unnormal Allel pattern removed")


                                                                                                                                            freq = frequency

                                                                                                                                            ### when the coverage is lower 5, convert :Alleles to ./. - calculation of readdepth necessary
                                                                                                                                            # parents should have another consideration that is set as standard, so that only the pools are effected by the change of the minreaddepth setting
                                                                                                                                            # the minimum for parents will be set to 5 statically
                                                                                                                                            nam_pool = nam[.!(.|((nam .==nam[posE1]), (nam .== nam[posE2])))]
                                                                                                                                            nam_parents = nam[.|((nam .==nam[posE1]), (nam .== nam[posE2]))]
                                                                                                                                             # progenies
                                                                                                                                          Threads.@threads for n in nam_pool
                                                                                                                                                                            # calulate the readdepth
                                                                                                                                                                            p = freq[Symbol(n)][:,2].+ freq[Symbol(n)][:,3]
                                                                                                                                                                            freq[Symbol(n)] = hcat(freq[Symbol(n)],p)
                                                                                                                                                                            freq[Symbol(n)] = Array(freq[Symbol(n)])

                                                                                                                                                             end


                                                                                                                                             # calculate the allel freq for each progeny
                                                                                                                                          println("prepare empty frames for Allel freq calculation")

                                                                                                                    Threads.@threads      for i in 3:length(nam); freq[Symbol(nam[i])]= hcat(freq[Symbol(nam[i])], Float32.(repeat(2:2, size(freq[Symbol(nam[i])],1))), Float32.(repeat(2:2, size(freq[Symbol(nam[i])],1)))); end
                                                                                                                                                            println("calculate Allel frequency for each progeny")
                                                                                                                                                            Threads.@threads for i in 3:length(nam)



                                                                                                                                                                  println(nam[i])
                                                                                                                                                                  freq[Symbol(nam[i])] = hcat(freq[:General][:,1:2], freq[Symbol(nam[i])])
                                                                                                                                                                  a1 = freq[Symbol(nam[i])][freq[Symbol(nam[posE1])][:,1].=="0/0",:]
                                                                                                                                                                  a1[:,7] .= a1[:,4] ./ a1[:,6]
                                                                                                                                                                  a1[:,8] .= a1[:,5] ./ a1[:,6]
                                                                                                                                                                  #a1 = hcat(a1,fill("g", 6917))
                                                                                                                                                                  a2 = freq[Symbol(nam[i])][freq[Symbol(nam[posE1])][:,1].=="1/1",:]
                                                                                                                                                                  a2[:,7] .= a2[:,5] ./ a2[:,6]
                                                                                                                                                                  a2[:,8] .= a2[:,4] ./ a2[:,6]
                                                                                                                                                                  #a2 = hcat(a2,fill("i", 3422))
                                                                                                                                                                  a3 = DataFrame(vcat(a1,a2))
                                                                                                                                                                  sort!(a3, [:x1,:x2])
                                                                                                                                                                  freq[Symbol(nam[i])] = Matrix(a3[:,3:8])

                                                                                                                                                               end


                                                                                                                                             ## calculate the allelfreq for both parents
                                                                                                                                            Threads.@threads for i in 3:length(nam)
                                                                                                                                                                      x = Symbol(nam[i])
                                                                                                                                                                      # subset by > -1 to get rid of the NaN values
                                                                                                                                                                      FREQ = sum((freq[Symbol(nam[i])][freq[Symbol(nam[i])][:,5].>-1,6] .* freq[Symbol(nam[i])][freq[Symbol(nam[i])][:,5].>-1,4]) ./ sum(freq[Symbol(nam[i])][freq[Symbol(nam[i])][:,5].>-1,4]))

                                                                                                                                                                     print(Crayon(foreground = :blue, bold = false),"Allelfreq over all loci of $x: $FREQ\n")

                                                                                                                                                                    end

                                                                                                                                                     #freq[Symbol(nam[i])][Symbol("Allelfreq",nam[posE1])][j]
                                                                                                                                                     pl = ["Chr", "Pos", "Ref", "Alt", "Qual", "Alleles", "Refcount", "Altcount", "Readcount"]
                                                                                                                                                     pl2 =  [string("Allelfreq_",nam[posE1]), string("Allelfreq_",nam[posE2])]


                                                                                                                                                     # make ever array a dataframe
                                                                                                                                                     Threads.@threads  for i in nam
                                                                                                                                                                                     freq[Symbol(i)] = hcat(freq[:General],freq[Symbol(i)])
                                                                                                                                                                                     freq[Symbol(i)] = DataFrame(freq[Symbol(i)])
                                                                                                                                                                                     if size(freq[Symbol(i)],2) == 9
                                                                                                                                                                                                                     rename!(freq[Symbol(i)], Symbol.(pl))
                                                                                                                                                                                     elseif size(freq[Symbol(i)],2) == 11
                                                                                                                                                                                                                     rename!(freq[Symbol(i)] ,Symbol.(vcat(pl,pl2)))

                                                                                                                                                                                    end
                                                                                                                                                                     end
                                                                                                                                                        # add the position information of the parents in the genolist file to the dictoinary - will be used in following steps to adress the correct name of founders
                                                                                                                                                        freq[:Parents] = nam[[posE1,posE2]]

                                                                                                                                                            return freq

                                                                                                                                                    end

"""
LINK SNP / INDEL information to an Object on the GENOME (like a GENE) \n
\n
relies on the output of function Allelfreqcalling \n
the  function is designed to extract the information of related genes or Marker to the detected SNPS with its Names, functions and Go terms. \n
The output id written to an other dict entry in the Allelfreqcalling output file, so that all information is stored in one single dictionary \n
Relatives can be Genes or Markers and extention of their bounds can be performed by the extention setting (default = true, extention will be performed) \n

The reference file style can have to different types \n
if merged is specified as false (default) \n
Gene_ID     chromosome:start-stop     [adiditonal like :  description       GO-IDs      InterPro-IDs] \n
Gene1         1A:32300-40332                          Ring Finger protein   GO:0000021, GO:00000120  IPR000184, IPR00000250 \n

 - if you use a marker instead of genes, create a table like this:
 Gene_ID    chromosome:start-stop     [adiditonal like :  description       GO-IDs      InterPro-IDs]
 Marker_1   1A:323001:323003            [..]

if merged is specified as true \n
- the only major difference is the already splited "chromosome:start-stop"   column. The rest is additional \n
Marker      Chr      start    end     Pos_genetic   Pos_phys    Gene_ID   chromosome:start-stop   description   GO-IDs  ..... \n
M1          1A        25440   32200       0.24        29250       Gene1     1A:25440-32200          HMA          GO:000014 \n


locate(info, freqfile, genolist,  presplit=false,  extention=true) \n
\n
- presplit - false => reference file looks like default type  \n
           - true  => reference file is already split for the position \n
- extension => shall the Objects be wided up or shall they keep their start and end position ? - true if extention should be performed (recommended) , false if not \n
  -> including the extention will give better results for the haplotype estimate \n
- gap => how far should the Object be extended ? e.g., how far a gene should be extended into the intergenic region (default = 0.45 - covering 90% of the intergenic region for up- and downstream gene) \n
 -> do not use more than 0.5 -> the regions will overlap (not recommended, but can be done) \n
\n
locate(info, freqfile, genolist, presplit=false,  extention=true, gap = 0.45) \n
\n
In case an Error occurs (like "key 5 not found"), just try to rerun the code again
\n
\n
# Example \n

freq  = Allelefreqcalling(SNPfile.vcf, overviewfile, 3,30,1000,1,2) \n
locate(refernce_file, freq, overviewfile)
"""
function locate(info, freqfile, genolist, presplit=false,  extention=true, gap = 0.45)
                                                  # gff file might not even be needed
                                                  println("read file")
                                                  label = CSV.read(genolist,DataFrame;  delim="\t")[1,1]
                                                  infofile = CSV.read(info,DataFrame; delim="\t",header=1)

                                # two possible ways can be gone - genetic map or classic physical map
                                if presplit == true
                                                    print(Crayon(foreground = :yellow, bold = false),"Extension is not calculated and Position should already be splited in a start and end position - should already exist \n")
                                                    ni = names(infofile)
                                                    inf = Matrix(infofile)
                                                    np = names(freqfile[Symbol(label)])
                                                    pp = Matrix(freqfile[Symbol(label)])
                                                    qerf = Dict()
                        Threads.@threads             for i in unique(pp[:,1])
                                                                                # subset
                                                                                println(i)
                                                                                ofni = inf[inf[:,2].== i ,:]
                                                                                qerf[Threads.threadid()] = pp[pp[:,1].== i , 1:2]
                                                                                # add the columns from the infotable to the snp table
                                                                                qerf[Threads.threadid()] = hcat(qerf[Threads.threadid()], Array{Union{Nothing, Any}}(nothing, size(qerf[Threads.threadid()],1), size(ni,1)))
                                                                                # create a "trash" file containing some value to set wenn no gene is matching
                                                                                pl = [".",".", 0, 0, 0.0, 0, ".", ".", ".", ".", ".", ".", "."]
                                                                                # write pl to all rows
                                                                                for n in 1:size(inf,2); qerf[Threads.threadid()][:,n+2] .= pl[n]; end
                                                                                # in the first iteration, give the column names to the Dict of Allelfreqcalling
                                                                                if haskey(freqfile, :Location)==false; freqfile[:Location] = Array{Union{Nothing, Any}}(nothing,0, 15); end
                                                                                #if haskey(freqfile, :Location)==false; freqfile[:Location] = DataFrame(describe(qerf)[:eltype], names(qerf), 0); end
                                                                                # do the check and write the information to an new DataFrame with position information of the snp and the gene both

                                                                                 for j in 1:size(qerf[Threads.threadid()],1), k in 1:size(ofni,1)

                                                                                                                  if qerf[Threads.threadid()][j,2] > ofni[k,3] && qerf[Threads.threadid()][j,2] < ofni[k,4]
                                                                                                                           for u in 1:size(inf,2)
                                                                                                                                           qerf[Threads.threadid()][j,u+2] = ofni[k,u]
                                                                                                                           end
                                                                                                                           continue
                                                                                                                  end
                                                                                              end
                                                   end
                                                   println("done")
                                                                 # write the information from the loop to the Dict entry
                                                                 for i in keys(qerf); freqfile[:Location] = vcat(freqfile[:Location],qerf[i]); end
                                                                 # convert array to a dataframe
                                                                 freqfile[:Location] = DataFrame(freqfile[:Location])
                                                                 rename!(freqfile[:Location], vcat(np[1:2], ni),makeunique=true)
                                                                 # rearrange the file according to the llcation file coming from the physical positions
                                                                 freqfile[:Location] = freqfile[:Location][:,[1,2,4,5,6,3,7,8,9,10,11,12,13,14,15]]
                                                                 println("Allocation of SNPs to the Genes successfully completed")

                                else # physical map
                                                  print(Crayon(foreground = :yellow, bold = false),"Position is splitted and the extention will be calculated \n")
                                                  # the chromosom, start and end have to be splited up
                                                  # two new split functions have been created just for this problem don here
                                                  function split3(file,sep1=":",sep2="-")
                                                                  es = vcat(split.(file,sep1)...)
                                                                  es2 = es[collect(2:2:length(es))]
                                                                  es21 = vcat(split.(es2,sep2)...)[collect(1:2:length(es))]
                                                                  es22 = vcat(split.(es2,sep2)...)[collect(2:2:length(es))]
                                                                  em = DataFrame(Array{Union{Nothing, Any}}(nothing, length(es2), 3))
                                                                  em[!,1] = es[collect(1:2:length(es))]
                                                                  em[!,2] = es21
                                                                  em[!,3] = es22
                                                                  return em
                                                  end

                                                  # split the position column into three columns
                                                  print(Crayon(foreground = :white, bold = false),"Split Position Column \n")
                                                  inq = split3(infofile[!,2],":","-")
                                                  inq[!,:x2] = parse.(Int,inq[!,:x2]) # convert to an integer from string
                                                  inq[!,:x3] = parse.(Int,inq[!,:x3])



                                                  #= extend the gene boundaries, due to the fact that with GBS data most of the SNPs will be located outside of genes
                                                     the rule for extention follows the same rule as has been applied to the extention of wheat with different extention sizes
                                                  =#
                                                  if extention == true
                                                           print(Crayon(foreground = :white, bold = false),"Extention of Gene Bounds is Performed \n")
                                                          inq[!, :start] .= 1
                                                          inq[!, :end] .= 1
                                                          Threads.@threads for i in 1:size(inq,1)
                                                                                  if i == 1; inq[!,:start][i] = inq[!,:x2][i]; end
                                                                                  if i == size(inq,1)
                                                                                                      inq[!,:end][i] = inq[!,:x3][i]
                                                                                  else
                                                                                          w = inq[!,:x2][i+1] - inq[!,:x3][i]
                                                                                          if w > 1000
                                                                                                  inq[!,:start][i+1] = floor(inq[!,:x2][i+1] - w * 0.45) # convert an float value to an integer
                                                                                                  inq[!,:end][i] = floor(inq[!,:x3][i] + w * 0.45)
                                                                                          else
                                                                                                  inq[!,:end][i] = inq[!,:x3][i]
                                                                                                  inq[!,:start][i+1] = inq[!,:x2][i+1]
                                                                                          end
                                                                                  end
                                                                          end
                                                                     print(Crayon(foreground = :green, bold = false),"Extention performed \n")
                                                                  else
                                                                          inq[!,:start] .= inq[!,:x2]
                                                                          inq[!,:end] .= inq[!,:x3]
                                                                  end
                                                  # test if there are overlappings - the overlaps results from the original ref data, extention of overlaps is due to the code above not intended
                                                  # overlaps might cause problems when trying to locate the snp correctly to a unique gene
                                           ### OVERLAPS will be kept if the following 3 '#' are existing #######
                                           #         inq[:check] = vcat(false,inq[2:size(inq)[1],:start] .> inq[1:size(inq)[1]-1,:end])
                                                  #link the intermed file with the infofile again
                                                  infofile = hcat(inq, infofile)
                                                  # if its a check == false and a low confidence gene, then remove it from data set
                                                  select!(infofile, Not(:x2))
                                                  select!(infofile, Not(:x3))

                                                  # link the SNPs to the gene - it only has to be done once, due to the fact that all libraries of SNPs reffer to the same loci
                                                  # saved in an additional dictionary
                                                  #= procedure:
                                                  subset the chromosoms, so that no mismatch with an other chromosom is possible
                                                  =#
                                                    ni = names(infofile)
                                                    inf = Matrix(infofile)
                                                    np = names(freqfile[Symbol(label)])
                                                    pp = Matrix(freqfile[Symbol(label)])
                                                    qerf = Dict()


                    print(Crayon(foreground = :blue, bold = false),"Allocate Polymorphsims to haplotype \n")
                      Threads.@threads              for i in unique(infofile[:,1])
                                                                                # subset
                                                                                println(i)
                                                                                ofni = inf[inf[:,1].== i ,:]
                                                                                qerf[Threads.threadid()] = pp[pp[:,1].== i , 1:2]
                                                                                # add the columns from the infotable to the snp table
                                                                                qerf[Threads.threadid()] = hcat(  qerf[Threads.threadid()], Array{Union{Nothing, Any}}(nothing, size(qerf[Threads.threadid()],1), size(ni,1)))
                                                                                # create a "trash" file containing some value to set wenn no gene is matching
                                                                                pl = [".", 0, 0, ".",".",".",".",".",".","." ]
                                                                                # write pl to all rows
                                                                                for n in 1:size(inf,2); qerf[Threads.threadid()][:,n+2] .= pl[n]; end
                                                                                # in the first iteration, give the column names to the Dict of Allelfreqcalling
                                                                                if haskey(freqfile, :Location)==false; freqfile[:Location] = Array{Union{Nothing, Any}}(nothing,0, 12); end
                                                                                #if haskey(freqfile, :Location)==false; freqfile[:Location] = DataFrame(describe(qerf)[:eltype], names(qerf), 0); end
                                                                                # do the check and write the information to an new DataFrame with position information of the snp and the gene both
                                                                                 for j in 1:size(qerf[Threads.threadid()],1), k in 1:size(ofni,1)

                                                                                                                  if   qerf[Threads.threadid()][j,2] > ofni[k,2] &&   qerf[Threads.threadid()][j,2] < ofni[k,3]
                                                                                                                           for u in 1:size(inf,2)
                                                                                                                                             qerf[Threads.threadid()][j,u+2] = ofni[k,u]
                                                                                                                           end
                                                                                                                           continue
                                                                                                                  end
                                                                                end
                                                   end

                                                                 # write the information from the loop to the Dict entry
                                                                 for i in keys(qerf); freqfile[:Location] = vcat(freqfile[:Location],qerf[i]); end
                                                                 # convert array to a dataframe
                                                                 freqfile[:Location] = DataFrame(freqfile[:Location])
                                                                 rename!(freqfile[:Location], vcat(np[1:2], ni))
                                                                 print(Crayon(foreground = :green, bold = false),"Allocation of SNPs to the Genes successfully completed \n")

                                          end # genetic of physical map

                                          print(Crayon(foreground = :white, bold = false),"Location file is sorted.. \n")
                                          sort!(freqfile[:Location], [:Chr, :Pos])

                                     end # genloc function


######## HAPLOTYPING ############
"""
HAPLOTYPE Frequency calculation based on locate information \n
\n
create haplotypes of the snps using the gene bound information \n
 haplotype have been proven to give a higher consistency in the calling of the correct allele frequency \n
 for the haplotype calculation 2 ways can be applyed, the first is refering to the readdepth given, so higher rd will contribute more to the Frequency (weighted = true - default, recommended) \n
 second way weights all snps equally, not depending if the are unweighted in readdepth (weighted = false)\n
\n
 - gene = haplotyping by Genes (Column has to have the name "Gene_ID") \n
 - Marker = haplotyping by Marker (Column has to have the name "Marker") \n
 - GO = haplotyping by go terms (Go Enrichment) (column has to have the name "GO-IDs") \n
 - PF = haplotyping by PF terms (column has to have the name "PFAM-IDs") \n

 haplotyping(freqfile,genolist,haplotype="gene", weighted=true) \n

#Example \n

freq  = Allelefreqcalling(SNPfile.vcf, overviewfile, 3,30,1000,1,2) \n
locate(refernce_file, freq, overviewfile)  \n
haplotyping(freq, overviewfile, "gene")
"""
function haplotyping(freqfile,genolist,haplotype="gene", weighted=true)
                                                            nam = CSV.read(genolist,DataFrame;  delim='\t')[!,:Name]

                                                            # get a unique list of all GO terms without "none" and "." ; lstrip - remove space from front if exists; vcat(split) - get a vector kind list from the GO terms in the Location subdict
                                                            if haplotype == "gene"
                                                                                    hap = unique(freqfile[:Location][!,:Gene_ID])
                                                                                    search = Symbol("Gene_ID")
                                                            elseif haplotype == "GO"
                                                                                    hap = lstrip.(unique(vcat(split.(freqfile[:Location][!,Symbol("GO-IDs")],',')...))[occursin.("GO",unique(vcat(split.(freqfile[:Location][!,Symbol("GO-IDs")],',')...)))])
                                                                                    search = Symbol("GO-IDs")
                                                            elseif haplotype == "Marker"
                                                                                    hap = unique(freqfile[:Location][!,:Marker])
                                                                                    search = Symbol("Marker")
                                                            elseif haplotype == "PF"
                                                                                    hap = lstrip.(unique(vcat(split.(freqfile[:Location][!,Symbol("PFAM-IDs") ],',')...))[occursin.("PF",unique(vcat(split.(freqfile[:Location][!,Symbol("PFAM-IDs") ],',')...)))])
                                                                                    search = Symbol("PFAM-IDs")
                                                            else
                                                            print(Crayon(foreground = :red, bold = false),"wrong haplotype specified, chose 'gene', 'GO' or 'PF' ")

                                                            end

                                                            # create a new dictionary for each progeny
                                                            for k in 3:length(nam); freqfile[Symbol(string(nam[k], "_Haplotypes"))] = DataFrame(Gene_ID=String[], Allelfreq_P2=Float64[], Allelfreq_P1=Float64[], SNPcount=Int64[], Readcount=Int64[], MinFreq=Float64[], MaxFreq=Float64[]); end
                                                                                 # the haplotype file information coming with info about the gene, function and FRequency of this particular haplotype
                                            println("Make Haplotypes using the Bounds calculated with Genloc..")
                                        ###########################################################


                                        # this part only preformes the haplotyping by gene
                                        if haplotype == "gene" || haplotype == "Marker"
                                           for j in 3:length(nam)
                                                              dat = hcat(freqfile[Symbol(nam[j])],freqfile[:Location][!,Symbol(search)])
                                                                                 function callsnps(dat,k)
                                                                                                         m = Int32[]
                                                                                                         for i in 1:nrow(dat)
                                                                                                                             if dat[i,:x1] != k
                                                                                                                             append!(m,i)
                                                                                                                             break
                                                                                                                             end
                                                                                                                         end
                                                                                                       if  !(isempty(m))
                                                                                                          dat2 = dat[1:m[1]-1,:]
                                                                                                          deleterows!(dat, 1:m[1]-1) # remove the gene snps that have already been worked on
                                                                                                        else
                                                                                                          dat2 = dat[1:nrow(dat),:]
                                                                                                        end

                                                                                                     return dat2
                                                                                end

                                                                                dat = dat[dat[!,:x1].!= ".",:] # remove not annotated snps from the list, highlighted by a "." in the x1 gene column
                                                                                sort!(dat, [:Chr, :x1])
                                                                                kap = unique(dat[!,:x1])
                                                                                println(nam[j])
                                                                                @showprogress for k in kap
                                                                                                          qwe = callsnps(dat,k)
                                                                                                          qwe = qwe[.!(isnan.(qwe[:,:11])),:]
                                                                                                          if isempty(qwe) == false
                                                                                                                 if size(qwe,1) > 1
                                                                                                                             if weighted == true
                                                                                                                                     p =  sum(qwe[!,:Readcount] .* qwe[!,Symbol("Allelfreq_$(freqfile[:Parents][2])")])/sum(qwe[!,:Readcount])
                                                                                                                                     q = 1 - sum(qwe[!,:Readcount] .* qwe[!,Symbol("Allelfreq_$(freqfile[:Parents][2])")])/sum(qwe[!,:Readcount]) # P1 freq
                                                                                                                                     wi = sum(qwe[!,:Readcount]) # readcount of the haplotype snps
                                                                                                                                     mi = minimum(qwe[!,Symbol("Allelfreq_$(freqfile[:Parents][2])")]) # min Allelfreq P2 of all markers for P2
                                                                                                                                     mx = maximum(qwe[!,Symbol("Allelfreq_$(freqfile[:Parents][2])")]) # max Allelfreq of all markers for P2
                                                                                                                                     push!(freqfile[Symbol(string(nam[j], "_Haplotypes"))],[k,p,q,size(qwe,1),wi, mi, mx])
                                                                                                                                     #append!(freqfile[Symbol(string(nam[j], "_Haplotypes"))],asd)
                                                                                                                             else # use the none weighted version where each snp contributes equally
                                                                                                                                     p = sum(qwe[!,Symbol("Allelfreq_$(freqfile[:Parents][2])")]) ./ size(qwe,1) # P2 freq
                                                                                                                                     q = 1 - sum(qwe[!,Symbol("Allelfreq_$(freqfile[:Parents][2])")]) ./ size(qwe,1) # P1 freq
                                                                                                                                     wi = sum(qwe[!,:Readcount]) # readcount of the haplotype snps
                                                                                                                                     mi = minimum(qwe[!,Symbol("Allelfreq_$(freqfile[:Parents][2])")]) # min Allelfreq P2 of all markers for P2
                                                                                                                                     mx = maximum(qwe[!,Symbol("Allelfreq_$(freqfile[:Parents][2])")]) # max Allelfreq of all markers for P2
                                                                                                                                     push!(freqfile[Symbol(string(nam[j], "_Haplotypes"))],[k,p,q,size(qwe,1),wi, mi,mx])
                                                                                                                            end
                                                                                                                 else
                                                                                                                     # calculate the correction of the Freq Calling based on the correction factor that has been saved in the freqfile dict "correction"

                                                                                                                     # make sure the Alt Alleles will be corrected - it is not directly constructable if P1 or P2 are the Ref Base call
                                                                                                                     if qwe[!,Symbol("Allelfreq_$(freqfile[:Parents][2])")] == qwe[!,:Altcount] ./ qwe[!,:Readcount]
                                                                                                                                             p = qwe[!,Symbol("Allelfreq_$(freqfile[:Parents][2])")]
                                                                                                                                             q = 1 - qwe[!,Symbol("Allelfreq_$(freqfile[:Parents][2])")][1]
                                                                                                                                             wi = sum(qwe[!,:Readcount]) # readcount of the haplotype snps
                                                                                                                                             mi = 0.0
                                                                                                                                             mx = 0.0
                                                                                                                                             push!(freqfile[Symbol(string(nam[j], "_Haplotypes"))],[k,p[1],q[1],size(qwe,1),wi, mi, mx])
                                                                                                                     else
                                                                                                                                             p = qwe[!,Symbol("Allelfreq_$(freqfile[:Parents][2])")]
                                                                                                                                             q = 1 - qwe[!,Symbol("Allelfreq_$(freqfile[:Parents][2])")][1]
                                                                                                                                             wi = sum(qwe[!,:Readcount]) # readcount of the haplotype snps
                                                                                                                                             mi = 0.0
                                                                                                                                             mx = 0.0
                                                                                                                                             push!(freqfile[Symbol(string(nam[j], "_Haplotypes"))],[k,p[1],q[1],size(qwe,1),wi,mi,mx])
                                                                                                                     end
                                                                                                                 end
                                                                                                             end

                                                                                                      end # of kap

                                                                                    end # of nam

                                                            else
                                                            Threads.@threads       for j in 3:length(nam)
                                                                                          for i in hap
                                                                                                               # the alogrhytm is following the list of genes, so all entries in the list and then alloaction the information to the single pools
                                                                                                               qwe = freqfile[Symbol(nam[j])][freqfile[:Location][search].== i,:] # get all snps in the tested gene
                                                                                                               filter!(row -> !isnan(row[:Allelfreq_P2]),qwe) # if the Freq contains any NaN, remove these
                                                                                                               # calculate the frequency based on which method has been choosen
                                                                                                                    if isempty(qwe) == false
                                                                                                                         if size(qwe,1) > 1
                                                                                                                                     if weighted == true
                                                                                                                                             p =  sum(qwe[!,:Readcount] .* qwe[!,Symbol("Allelfreq_$(freqfile[:Parents][2])")])/sum(qwe[!,:Readcount])
                                                                                                                                             q = 1 - sum(qwe[!,:Readcount] .* qwe[!,Symbol("Allelfreq_$(freqfile[:Parents][2])")])/sum(qwe[!,:Readcount]) # P1 freq
                                                                                                                                             wi = sum(qwe[!,:Readcount]) # readcount of the haplotype snps
                                                                                                                                             mi = minimum(qwe[!,Symbol("Allelfreq_$(freqfile[:Parents][2])")]) # min Allelfreq P2 of all markers for P2
                                                                                                                                             mx = maximum(qwe[!,Symbol("Allelfreq_$(freqfile[:Parents][2])")]) # max Allelfreq of all markers for P2
                                                                                                                                             push!(freqfile[Symbol(string(nam[j], "_Haplotypes"))],[i,p,q,size(qwe,1),wi, mi, mx])
                                                                                                                                             #append!(freqfile[Symbol(string(nam[j], "_Haplotypes"))],asd)
                                                                                                                                     else # use the none weighted version where each snp contributes equally
                                                                                                                                             p = sum(qwe[!,Symbol("Allelfreq_$(freqfile[:Parents][2])")]) ./ size(qwe,1) # P2 freq
                                                                                                                                             q = 1 - sum(qwe[!,Symbol("Allelfreq_$(freqfile[:Parents][2])")]) ./ size(qwe,1) # P1 freq
                                                                                                                                             wi = sum(qwe[!,:Readcount]) # readcount of the haplotype snps
                                                                                                                                             mi = minimum(qwe[!,Symbol("Allelfreq_$(freqfile[:Parents][2])")]) # min Allelfreq P2 of all markers for P2
                                                                                                                                             mx = maximum(qwe[!,Symbol("Allelfreq_$(freqfile[:Parents][2])")]) # max Allelfreq of all markers for P2
                                                                                                                                             push!(freqfile[Symbol(string(nam[j], "_Haplotypes"))],[i,p,q,size(qwe,1),wi, mi,mx])
                                                                                                                                    end
                                                                                                                         else
                                                                                                                             # calculate the correction of the Freq Calling based on the correction factor that has been saved in the freqfile dict "correction"

                                                                                                                             # make sure the Alt Alleles will be corrected - it is not directly constructable if P1 or P2 are the Ref Base call
                                                                                                                             if qwe[!,Symbol("Allelfreq_$(freqfile[:Parents][2])")] == qwe[!,:Altcount] ./ qwe[!,:Readcount]
                                                                                                                                                     p = qwe[!,Symbol("Allelfreq_$(freqfile[:Parents][2])")]
                                                                                                                                                     q = 1 - qwe[!,Symbol("Allelfreq_$(freqfile[:Parents][2])")][1]
                                                                                                                                                     wi = sum(qwe[!,:Readcount]) # readcount of the haplotype snps
                                                                                                                                                     mi = 0.0
                                                                                                                                                     mx = 0.0
                                                                                                                                                     push!(freqfile[Symbol(string(nam[j], "_Haplotypes"))],[i,p[1],q[1],size(qwe,1),wi, mi, mx])
                                                                                                                             else
                                                                                                                                                     p = qwe[!,Symbol("Allelfreq_$(freqfile[:Parents][2])")]
                                                                                                                                                     q = 1 - qwe[!,Symbol("Allelfreq_$(freqfile[:Parents][2])")][1]
                                                                                                                                                     wi = sum(qwe[!,:Readcount]) # readcount of the haplotype snps
                                                                                                                                                     mi = 0.0
                                                                                                                                                     mx = 0.0
                                                                                                                                                     push!(freqfile[Symbol(string(nam[j], "_Haplotypes"))],[i,p[1],q[1],size(qwe,1),wi,mi,mx])
                                                                                                                             end
                                                                                                                         end
                                                                                                                      end
                                                                                                        end # of i loop
                                                                                            end # of j loop

                                                            end # of haplotype selection algorythm

                                                            # give an information how many Genes have 1,2,3,4 or more than 4 snps linked to an haplotype
                                                            println("Calculate overview over all Levels of Haplotypes")
                                                                    # build an empty table where the information will be written to
                                                                    haplotypecheck = DataFrame(SNPperGene=["all","1","2","3","4",">4"])
                                                          Threads.@threads  for j in 3:length(nam)
                                                                                  haplotypecheck[Symbol(nam[j])] = [size(freqfile[Symbol(string(nam[j], "_Haplotypes"))],1),size(freqfile[Symbol(string(nam[j], "_Haplotypes"))][freqfile[Symbol(string(nam[j], "_Haplotypes"))][!,:SNPcount] .==1,:],1),
                                                                                  size(freqfile[Symbol(string(nam[j], "_Haplotypes"))][freqfile[Symbol(string(nam[j], "_Haplotypes"))][!,:SNPcount] .==2,:],1),
                                                                                  size(freqfile[Symbol(string(nam[j], "_Haplotypes"))][freqfile[Symbol(string(nam[j], "_Haplotypes"))][!,:SNPcount] .==3,:],1),
                                                                                  size(freqfile[Symbol(string(nam[j], "_Haplotypes"))][freqfile[Symbol(string(nam[j], "_Haplotypes"))][!,:SNPcount] .==4,:],1),
                                                                                  size(freqfile[Symbol(string(nam[j], "_Haplotypes"))][freqfile[Symbol(string(nam[j], "_Haplotypes"))][!,:SNPcount] .>4,:],1)]

                                                                          end

                                                                          # change the name to the once give in the genolist
                                                                          for k in 3:length(nam); rename!(freqfile[Symbol(string(nam[k], "_Haplotypes"))], :Allelfreq_P2 => "Allelefreq_$(freqfile[:Parents][2])", :Allelfreq_P1 => "Allelefreq_$(freqfile[:Parents][1])"); end

                                                                           print(Crayon(foreground = :green, bold = false),"Haplotyping Finished \n")
                                                                          return haplotypecheck
                                     end # of function haplotyping

######## Contig mean ###########
# create summary stats and mean Frequency values for Contigs of a defined size - homolog to haplotyping, but with fixed windows
"""
if no information on the gene or marker position can be retained from any source, the haplotypes can be creaded by generating contigs according to a certain size. \n
  inbreeding species can use bigger contigs, while outbreeding species should have smaller contig size \n
  default contig size = 10 mio bp \n

# Example \n

freq  = Allelefreqcalling(SNPfile.vcf, overviewfile, 3,30,1000,1,2) \n
contig(freq, overviewfile, 10000)
"""
function contig(freqfile, genolist, contigsize = 10000000)
                                                              nam = CSV.read(genolist,DataFrame;  delim='\t')[!,:Name]
                                                              n = contigsize
                                                              # calculate the maximum count of steps
                                                              steps = Int64(round(maximum(freqfile[:Location][!,:Pos]) / contigsize))

                                                              freqfile[:Location][!,:Contigs] .= "."
                                                              freqfile[:Location][!,:Cstart] .= 0
                                                              freqfile[:Location][!,:Cend] .= 0


                                                              println("prepare the Location file - create the contig bounds")
                                                              Threads.@threads for i in unique(freqfile[:Location][!,:Chr])
                                                                                                      set = freqfile[:Location][freqfile[:Location][!,:Chr] .== i,:]
                                                                                                      set[!,:start] .= 0
                                                                                                      set[!,:end] .= 0

                                                                                                      for j in 1:steps
                                                                                                                     set[.&(set[!,:Pos].> n* (j-1), set[!,:Pos].< n*j),:Contigs] .= "$(i)Contig$j"
                                                                                                                     set[.&(set[!,:Pos].> n* (j-1), set[!,:Pos].< n*j),:start] .= n * (j-1)
                                                                                                                     set[.&(set[!,:Pos].> n* (j-1), set[!,:Pos].< n*j),:end] .= n * j

                                                                                                      end
                                                                                                      freqfile[:Location][freqfile[:Location][!,:Chr] .== i,:Contigs] .= set[!,:Contigs]
                                                                                                      freqfile[:Location][freqfile[:Location][!,:Chr] .== i,:Cstart] .= set[!,:start]
                                                                                                      freqfile[:Location][freqfile[:Location][!,:Chr] .== i,:Cend] .= set[!,:end]

                                                                  end

                                                              # merge by Contigs
                                                              for k in 3:length(nam); freqfile[Symbol(string(nam[k], "_Contigs"))] = DataFrame(Gene_ID=String[], Allelfreq_P2=Float64[], Allelfreq_P1=Float64[], SNPcount=Int64[], Readcount=Int64[], MinFreq=Float64[], MaxFreq=Float64[]); end

                                                                println("calculate the Contig Haplotypes and save it to a new Dict entry")
                                                                Threads.@threads for j in 3:length(nam)
                                                                                     for k in unique(freqfile[:Location][!,:Contigs])

                                                                                                            qwe = freqfile[Symbol(nam[j])][freqfile[:Location][!,:Contigs].==k,:]
                                                                                                            qwe = qwe[.!(isnan.(qwe[:,:11])),:]
                                                                                                            if isempty(qwe) == false
                                                                                                                                    p =  sum(qwe[!,:Readcount] .* qwe[!,Symbol("Allelfreq_$(freqfile[:Parents][2])")])/sum(qwe[!,:Readcount])
                                                                                                                                    q = 1 - sum(qwe[!,:Readcount] .* qwe[!,Symbol("Allelfreq_$(freqfile[:Parents][2])")])/sum(qwe[!,:Readcount]) # P1 freq
                                                                                                                                    wi = sum(qwe[!,:Readcount]) # readcount of the haplotype snps
                                                                                                                                    mi = minimum(qwe[!,Symbol("Allelfreq_$(freqfile[:Parents][2])")]) # min Allelfreq P2 of all markers for P2
                                                                                                                                    mx = maximum(qwe[!,Symbol("Allelfreq_$(freqfile[:Parents][2])")]) # max Allelfreq of all markers for P2
                                                                                                                                    push!(freqfile[Symbol(string(nam[j], "_Contigs"))],[k,p,q,size(qwe,1),wi, mi, mx])
                                                                                                            end
                                                                                                    end
                                                                end
                                                                # change the name to the once give in the genolist
                                                                for k in 3:length(nam); rename!(freqfile[Symbol(string(nam[k], "_Contigs"))], :Allelfreq_P2 => "Allelefreq_$(freqfile[:Parents][2])", :Allelfreq_P1 => "Allelefreq_$(freqfile[:Parents][1])"); end

                                                                 print(Crayon(foreground = :green, bold = false),"Haplotyping Finished \n")

            end # of contig function


"""
merge the replicates together for \n
check the Dict() and check how the naming is done  - can be either _Haplotypes (default) or _Contigs \n
you can merge replicates if 2 or 3 replicate sets are present - for more replicates please check the code and adjust by yourself \n

# Example \n

freq  = Allelefreqcalling(SNPfile.vcf, overviewfile, 3,30,1000,1,2) \n
locate(refernce_file, freq, overviewfile) \n
haplotyping(freq, overviewfile, "gene") \n
stacker(freq, overviewfile, "_Haplotypes") \n
"""
function stacker(freq, genolist, s = "_Haplotypes")
                                                     println("merge the replicates for $s together")
                                                     set = collect(keys(freq))[occursin.(s, String.(collect(keys(freq))))]
                                                     gl = CSV.read(genolist, DataFrame; delim="\t")

                                                     ev = unique(gl[!,:Env])[2:size(unique(gl[!,:Env]),1)]
                                                     gen = unique(gl[!,:Generation])[2:size(unique(gl[!,:Generation]),1)]
                                                     # define the sets that should be matched togehter
                                                     for i in ev , j in gen
                                                                         qay = Symbol.(string.(gl[.&(gl[!,:Env].==i, gl[!,:Generation].==j),:Name], s))
                                                                         println(qay)
                                                                         # calculate the mean allelfreq and the sum of the readcount
                                                                         if size(qay,1) == 2
                                                                             # delete the haplotypes that have no replicate in the other set
                                                                             freq[qay[1]] = freq[qay[1]][ .!(isnothing.(indexin(freq[qay[1]][!,:Gene_ID], freq[qay[2]][!,:Gene_ID]))),:]
                                                                             freq[qay[2]] = freq[qay[2]][ .!(isnothing.(indexin(freq[qay[2]][!,:Gene_ID], freq[qay[1]][!,:Gene_ID]))),:]

                                                                                             a1 = (freq[qay[1]][!,Symbol("Allelefreq_$(freq[:Parents][2])")] .* freq[qay[1]][!,:Readcount] + freq[qay[2]][!,Symbol("Allelefreq_$(freq[:Parents][2])")] .* freq[qay[2]][!,:Readcount]) ./ (freq[qay[1]][!,:Readcount] .+ freq[qay[2]][!,:Readcount] ) # Allelfreq ISR
                                                                                             a2 = (freq[qay[1]][!,Symbol("Allelefreq_$(freq[:Parents][1])")] .* freq[qay[1]][!,:Readcount] + freq[qay[2]][!,Symbol("Allelefreq_$(freq[:Parents][1])")] .* freq[qay[2]][!,:Readcount]) ./ (freq[qay[1]][!,:Readcount] .+ freq[qay[2]][!,:Readcount] ) # Allelfreq Golf
                                                                                             a3 = maximum.(vcat.(freq[qay[1]][!,:SNPcount],freq[qay[2]][!,:SNPcount]))
                                                                                             a4 = freq[qay[1]][!,:Readcount] .+ freq[qay[2]][!,:Readcount]
                                                                                             a5 = minimum.(vcat.(freq[qay[1]][!,:MinFreq],freq[qay[2]][!,:MinFreq] ))
                                                                                             a6 = maximum.(vcat.(freq[qay[1]][!,:MaxFreq],freq[qay[2]][!,:MaxFreq] ))
                                                                        elseif size(qay,1) == 3
                                                                                    set = vcat(freq[qay[1]][!,:Gene_ID], freq[qay[2]][!,:Gene_ID],freq[qay[3]][!,:Gene_ID])
                                                                                    set = countmap(set)
                                                                                    # get all haplotypes that are present in all 3
                                                                                    df = String[]
                                                                                    for i in String.(keys(set))
                                                                                                                  if set[i] == 3; push!(df, i); end
                                                                                    end
                                                                                    freq[qay[1]] = freq[qay[1]][ .!(isnothing.(indexin(freq[qay[1]][!,:Gene_ID], df))),:]
                                                                                    freq[qay[2]] = freq[qay[2]][ .!(isnothing.(indexin(freq[qay[2]][!,:Gene_ID], df))),:]
                                                                                    freq[qay[3]] = freq[qay[3]][ .!(isnothing.(indexin(freq[qay[3]][!,:Gene_ID], df))),:]



                                                                                          a1 = (freq[qay[1]][!,Symbol("Allelefreq_$(freq[:Parents][2])")] .* freq[qay[1]][!,:Readcount] + freq[qay[2]][!,Symbol("Allelefreq_$(freq[:Parents][2])")] .* freq[qay[2]][!,:Readcount] + freq[qay[3]][!,Symbol("Allelefreq_$(freq[:Parents][2])")] .* freq[qay[3]][!,:Readcount]) ./ (freq[qay[1]][!,:Readcount] .+ freq[qay[2]][!,:Readcount] .+ freq[qay[3]][!,:Readcount]) # Allelfreq P2
                                                                                          a2 = (freq[qay[1]][!,Symbol("Allelefreq_$(freq[:Parents][1])")] .* freq[qay[1]][!,:Readcount] + freq[qay[2]][!,Symbol("Allelefreq_$(freq[:Parents][1])")] .* freq[qay[2]][!,:Readcount] + freq[qay[3]][!,Symbol("Allelefreq_$(freq[:Parents][1])")] .* freq[qay[3]][!,:Readcount]) ./ (freq[qay[1]][!,:Readcount] .+ freq[qay[2]][!,:Readcount] .+ freq[qay[3]][!,:Readcount]) # Allelfreq P1
                                                                                          a3 = maximum.(vcat.(freq[qay[1]][!,:SNPcount],freq[qay[2]][!,:SNPcount], freq[qay[3]][!,:SNPcount]))
                                                                                          a4 = freq[qay[1]][!,:Readcount] .+ freq[qay[2]][!,:Readcount] .+ freq[qay[3]][!,:Readcount]
                                                                                          a5 = minimum.(vcat.(freq[qay[1]][!,:MinFreq],freq[qay[2]][!,:MinFreq], freq[qay[3]][!,:MinFreq] ))
                                                                                          a6 = maximum.(vcat.(freq[qay[1]][!,:MaxFreq],freq[qay[2]][!,:MaxFreq], freq[qay[3]][!,:MaxFreq] ))
                                                                         else
                                                                              println("No replicates to merge")
                                                                              continue

                                                                         end

                                                                         a = DataFrame(hcat(freq[qay[1]][!,:Gene_ID], a1,a2,a3,a4,a5,a6))
                                                                         rename!(a, names(freq[qay[1]]))
                                                                         # push to the Dictionary
                                                                         saver = Symbol(string("F",j,"E",i,s))
                                                                         freq[saver] = a
                                                     end
      end
        # stacker(freq, genolist)

"""
add the position information to the Haplotype \n
\n

col = "Gene_ID" |  col = ["Gene_ID", "Marker"] | col = ["Gene_ID", "GO-IDs"] | col = ["Gene_ID", "PFAM-IDs"] | col = ["Gene_ID", "InterPro-IDs"] \n
col = "Gene_ID" => if both the freq[:Location] and freq[:Haplotype] column have the same name code \n
 col = ["Gene_ID", "Marker"] => if the column to merge within freq[:Location] is named "Marker" \n
# Example: \n

freq  = Allelefreqcalling(SNPfile.vcf, overviewfile, 3,30,1000,1,2) \n
locate(refernce_file, freq, overviewfile) \n
haplotyping(freq, overviewfile, "gene") \n
annotate( freq[:F2E2_Haplotypes], freq, "Gene_ID")

"""
function annotate(file, freqfile, col = "Gene_ID")
                                      q3 = unique(freqfile[:Location][!,3:12])
                                      if typeof(col) == String
                                                                            file = join(file, q3, on=Symbol(col))
                                      elseif typeof(col) == Array{String,1}
                                                                            file = join(file, q3, on=[Symbol(col[1]) => Symbol(col[2])], makeunique=true)
                                      else
                                      println("""The columns to merge with are not correctly specified, use: \n col = "Gene_ID" |  col = ["Gene_ID", "Marker"] | col = ["Gene_ID", "GO-IDs"] | col = ["Gene_ID", "PFAM-IDs"] | col = ["Gene_ID", "InterPro-IDs"]""")
                                      end

                                      return file

     end
        # annotate( freq[:F23E2_Haplotypes], freq)



end # module
