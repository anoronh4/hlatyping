hlascan_bin="/juno/work/ccs/noronhaa/tools/hla_scan_r_v2.1.4"
db="/juno/work/ccs/noronhaa/hlascan/db/HLA-ALL.IMGT"
samtools_="/opt/common/CentOS_7/samtools/samtools-1.9/bin/samtools"
bed_folder="/juno/work/ccs/noronhaa/hlascan/gen_hlascan_WES25TN/filtered_bed"

outDir = params.outDir
bamMapping = Channel.fromPath(params.bamMapping)
	.splitCsv(sep:'\t', header:true)
	.map{ row -> 
		[row.SAMPLE, file(row.BAM), file(row.BAI)]
	}.unique()
	.into{inBam_Ch; inBam4Polysolver3; inBam4Polysolver4; inBam4Bam2Fastq; inBam4HLALA; inBam4SOAP}

process filterBam {
tag {sampleid}

input:
set sampleid, file(bam), file(bai) from inBam_Ch
file(bed_folder) from Channel.value([file(bed_folder)])
output:
set sampleid, file("*filtered.bam"),file("*filtered.bam.bai") into filteredBam_Ch1, filteredBam_Ch2, filteredBam_Ch3
// set sampleid, file("${sampleid}.1.fastq"), file("${sampleid}.2.fastq") into filteredFastq

//when: 1==0

script:
bed=bed_folder + "/filtered.bed"
"""
$samtools_ view -b -hL $bed $bam > ${sampleid}.filtered.bam
$samtools_ index ${sampleid}.filtered.bam
"""
}

Channel.from("HLA-A","HLA-B","HLA-C","HLA-DMA","HLA-DMB","HLA-DOA","HLA-DOB","HLA-DPA1","HLA-DPB1","HLA-DQA1","HLA-DQB1","HLA-DRA","HLA-DRB1","HLA-DRB5").into{ GOI_Ch ; GOI_Ch2 }

process filterGene_and_hlascan {
tag { sampleid + "@" + gene }
publishDir "${outDir}/filter1/${sampleid}/", mode: params.publishDirMode

input:
each gene from GOI_Ch
set sampleid, file(bam), file(bai) from filteredBam_Ch1
file(bed_folder) from Channel.value([file(bed_folder)])

output:
file("${sampleid}.${gene}.results.txt") into hlascan_output

script:
bed=bed_folder + "/filtered.${gene}.bed"
"""
$samtools_ view -b -hL $bed $bam > ${sampleid}.filtered.${gene}.bam
$samtools_ index ${sampleid}.filtered.${gene}.bam

${hlascan_bin} \\
	-b ${sampleid}.filtered.${gene}.bam \\
	-d ${db} \\
	-v 37 -t 2 \\
	-g $gene \\
	> ${sampleid}.${gene}.results.txt | true

"""
}

process filterGene_and_hlascan2 {
tag { sampleid + "@" + gene }
publishDir "${outDir}/filter2/${sampleid}/", mode: params.publishDirMode

input:
each gene from GOI_Ch2
set sampleid, file(bam), file(bai) from filteredBam_Ch2
file(bed_folder) from Channel.value([file(bed_folder)])

output:
file("${sampleid}.${gene}.results.txt") into hlascan_output2 

script:
bed=bed_folder + "/filtered.${gene}.bed"
"""
${hlascan_bin} \\
	-b $bam \\
	-d ${db} \\
	-v 37 -t 2 \\
	-g $gene \\
	> ${sampleid}.${gene}.results.txt | true
"""
}


process calculate_coverage {
tag {sampleid}
publishDir "${outDir}/coverage/", mode: params.publishDirMode

input:
set sampleid, file(bam), file(bai) from filteredBam_Ch3
file(bed) from Channel.value([file("/juno/work/ccs/noronhaa/hlascan/filtered_regions/annot/filtered.bed")])

output: 
file("${sampleid}.coverage")

script:
"""
samtools bedcov ${bed} ${bam} > ${sampleid}.coverage 
"""
}

process polysolver_v4 {
tag {sampleid}
container = "sachet/polysolver:v4"

cpus = { 1 + (1 * task.attempt) } 
memory = 8.GB 

publishDir "${outDir}/polysolver_v4/", mode: params.publishDirMode

input: 
set sampleid, file(bam), file(bai) from inBam4Polysolver4

output:
file("${outputPrefix}.hla.txt") 

script:
outputPrefix = "${sampleid}"
outputDir = "."
tmpDir = "${outputDir}-nf-scratch"
genome_ = "hg19" 
"""
cp /home/polysolver/scripts/shell_call_hla_type .

sed -i "171s/TMP_DIR=.*/TMP_DIR=${tmpDir}/" shell_call_hla_type 
bash shell_call_hla_type \\
	${bam} \\
	Unknown \\
	1 \\
	${genome_} \\
	STDFQ \\
	0 \\
	${outputDir}

mv winners.hla.txt ${outputPrefix}.hla.txt
"""
}

process polysolver_v3 {
tag {sampleid}
container = "/opt/common/CentOS_6-dev/polysolver/v3/polysolver.img"
scratch = true

cpus = { 1 + (1 * task.attempt) }
memory = 8.GB 

publishDir "${outDir}/polysolver_v3/", mode: params.publishDirMode

input: 
set sampleid, file(bam), file(bai) from inBam4Polysolver3

output:
file("${outputPrefix}.hla.txt") 

script:
outputPrefix = "${sampleid}"
outputDir = "."
tmpDir = "${outputDir}-nf-scratch"
genome_ = "hg19" 
bin_path = "/usr/local/libexec/polysolver/scripts/shell_call_hla_type"
"""
export PSHOME=/usr/local/libexec/polysolver
export SAMTOOLS_DIR=/usr/local/libexec/samtools
export JAVA_DIR=/usr/bin
export NOVOALIGN_DIR=\$PSHOME/binaries
export NUM_THREADS=4


bash ${bin_path} ${bam} Unknown 1 ${genome_} STDFQ 0 ${outputDir}
mv winners.hla.txt ${outputPrefix}.hla.txt
"""
}

process Bam2Fastq {
tag {sampleid}
scratch = true
container = "broadinstitute/gatk:4.1.9.0"

input:
set val(sampleid), file(bam), file(bai) from inBam4Bam2Fastq

output:
//set val(sampleid), file("${sampleid}.1.fastq"), file("${sampleid}.2.fastq") into inFq4HlaHd
set val(sampleid), file("sample.hla.1.fastq"), file("sample.hla.2.fastq") into inFq4HlaHd, inFq4OptiType, inFq4Kourami

script:
"""
#samtools view -bh -F 268 ${bam} 6:28,477,797-33,448,354 | samtools sort -n - > sorted.bam
#samtools fixmate -O BAM sorted.bam fixed.bam 
#samtools view -bh -f 1 fixed.bam > sample.mhc.bam
#samtools view -bh -F 260 -f 8 ${bam} 6:28,477,797-33,448,354 > sample.mhcR1.bam
#samtools view -bh -F 264 -f 4 ${bam} 6:28,477,797-33,448,354 > sample.mhcR2.bam
#samtools view -b -f 12 ${bam} > sample.unmap.bam
#samtools merge sample.merge.nsort.bam sample.unmap.bam sample.mhc.bam sample.mhcR1.bam sample.mhcR2.bam
#samtools merge sample.merge.nsort.bam sample.unmap.bam sample.mhc.bam
#samtools sort sample.merge.nsort.bam > sample.merge.bam
samtools  view -@4 -h -b -u -f 4 ${bam} > unmapped_bam
samtools  view -@4 -h -b -u ${bam} 6:28,477,797-33,448,354 > mhc_mapped_bam 
samtools merge -@4 -u unsorted_bam unmapped_bam mhc_mapped_bam
samtools sort -@4 -n unsorted_bam > sorted_bam
samtools fastq -@2 -N -1 sample.hlatmp.1.fastq -2 sample.hlatmp.2.fastq -s /dev/null -0 /dev/null sorted_bam
#gatk SamToFastq I=sample.merge.bam F=sample.hlatmp.1.fastq F2=sample.hlatmp.2.fastq
cat sample.hlatmp.1.fastq |awk '{if(NR%4 == 1){O=\$0;gsub("/1"," 1",O);print O}else{print \$0}}' > sample.hla.1.fastq
cat sample.hlatmp.2.fastq |awk '{if(NR%4 == 1){O=\$0;gsub("/2"," 2",O);print O}else{print \$0}}' > sample.hla.2.fastq

#samtools view -@ ${task.cpus} -h -F 4 -f 0x40 -b1 ${bam} 6:28,477,797-33,448,354 > mapped.optitype.1.bam
#samtools view -@ ${task.cpus} -h -F 4 -f 0x80 -b1 ${bam} 6:28,477,797-33,448,354 > mapped.optitype.2.bam
"""

}

process HLAHD {
tag {sampleid}
container = "cmopipeline/hlahd:1.4"
scratch = true

cpus = { 3 * task.attempt }
memory = 10.GB

publishDir "${outDir}/hla-hd/", mode: params.publishDirMode

input:
set val(sampleid), file(fq1), file(fq2) from inFq4HlaHd 

output:
file("${sampleid}_final.result.txt")

script:
"""
install_dir=/hlahd.1.4.0
if [[ \$( ulimit -n ) -lt 1024 ]] ; then ulimit -n 1024 ;fi
bash \$install_dir/bin/hlahd.sh -t ${task.cpus} -m 100 -f \$install_dir/freq_data \
  ${fq1} ${fq2} \
  \$install_dir/HLA_gene.split.txt \
  \$install_dir/dictionary ${sampleid} .

cp ${sampleid}/result/${sampleid}_final.result.txt . 
"""
}

process OptiType {

tag {sampleid}
//container = "fred2/optitype:release-v1.3.1"
//container = "fred2/optitype:latest"
container = "quay.io/biocontainers/optitype:1.3.5--0"
scratch = true

cpus = { 2 * task.attempt }
memory = 8.GB

publishDir "${outDir}/optitype/", mode: params.publishDirMode

input:
set val(sampleid), file(fq1), file(fq2) from inFq4OptiType 

output:
file "$sampleid"

script:
"""
echo "[mapping]
razers3=\$(which razers3)
threads=16

[ilp]
solver=glpk
threads=1

[behavior]
deletebam=true
unpaired_weight=0
use_discordant=false
" > config.ini

OptiTypePipeline.py -i ${fq1} ${fq2} -c config.ini --dna --prefix $sampleid --outdir ${sampleid} 
"""
}

process HLALA {

tag { "${sampleid}" }
publishDir "${outDir}/hlala/", mode: params.publishDirMode

container = "cmopipeline/hlala:0.0.1-test"
cpus = { 4 * task.attempt }
memory = 8.GB


input:
  set sampleid, file(bam), file(bai) from inBam4HLALA 
  file(graphdir) from Channel.value([file("/juno/work/ccs/noronhaa/hlascan/gen_hlascan_repo/hlascan_nextflow/small_test/HLA-LA/graphs/")])

output:
  file("${sampleid}.hlala.tsv")

script:
  """
  export LC_ALL=C
  HLA-LA.pl \\
    --BAM ${bam} \\
    --graph PRG_MHC_GRCh38_withIMGT \\
    --sampleID sample \\
    --maxThreads ${task.cpus * 2} \\
    --customGraphDir ${graphdir} \\
    --workingDir ./

  cp sample/hla/R1_bestguess_G.txt ${sampleid}.hlala.tsv 
  """
}

process SOAPHLA {

tag { "${sampleid}" }
publishDir "${outDir}/soap/", mode: params.publishDirMode

input:
set sampleid, file(bam), file(bai) from inBam4SOAP

output:
file("output.tsv")

when: 1== 0

script:
"""
MHC_autopipeline -i $bam -od output/ -v hg19
"""
}

process Kourami {

tag { "${sampleid}" }
publishDir "${outDir}/kourami/", mode: params.publishDirMode

container = "cmopipeline/kourami:0.0.1"

input:
set sampleid, file(fq1), file(fq2) from inFq4Kourami 
file(resourcedir) from Channel.value([file("/juno/work/ccs/noronhaa/tools/kourami/")])

output:
file("${sampleid}.result")

script:
"""
bwa mem -t ${task.cpus *2 } ${resourcedir}/db/All_FINAL_with_Decoy.fa.gz ${fq1} ${fq2} | samtools view -Sb - > ${sampleid}.kourami.bam
java -jar /opt/kourami-0.9.6/target/Kourami.jar -d ${resourcedir}/db/ -o ${sampleid} ${sampleid}.kourami.bam
"""
}
