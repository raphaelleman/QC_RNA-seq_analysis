########## some useful functions
in_error=0 # will be 1 if a file or cmd not exist
echo_on_stderr () {
    (>&2 echo -e "$*")
}

test_command_if_exist () {
    command -v $* >/dev/null 2>&1 && echo 0 || { echo 1; }
}

test_file_if_exist () {
    if [ ! -f $* ]; then
        echo_on_stderr "File $* not found! Will abort."
        echo 1
    else
        echo 0
    fi
}

########## Help message
messageHelp="Usage: $0 [options] <bam files>\n
    \n
    -B, --bam /path/to/bam/\n\t\trepository of the BAM files\n
    -O, --output /path/to/output/\n\trepository of the output files\n
    -p, --paired\n\tIf data are paired\n
    --samtools /path/to/samtools\n\tpath to samtools executable\n
    --htseq /path/to/htseq\n\tpath to htseq count executable (python script 'htseq-count')\n
    --dexseq /path/to/dexseq script\n\tpath to the dexseq python script\n
    --refGene /path/to/geneCoordonate.gtf\n\tpath to the gene coordinate in GTF format\n
    --refExon /path/to/exonCoordonate.gtf\n\tpath to the exon coordinate in GTF format\n
    -h, --help\n\tPrint this help message and exit"

## exit if not enough arguments
if [ $# -lt 1 ]; then
    echo -e $messageHelp
    exit
fi

paired="no"

while [[ $# -gt 0 ]]; do
   key=$1
   case $key in

       -O|--output)
       out_path="`readlink -v -f $2`"
       shift 2 # shift past argument and past value
       ;;

       -p|--paired)
       paired="yes"
       shift # shift past argument
       ;;

       --htseq)
       htseq="$2"
       shift 2 # shift past argument and past value
       ;;

       --samtools)
       samtools="$2"
       shift 2 # shift past argument and past value
       ;;

       --dexseq)
       dexseq="$2"
       shift 2 # shift past argument and past value
       ;;

       --refGene)
       refGene="`readlink -v -f $2`"
       shift 2 # shift past argument and past value
       ;;

       --refExon)
       refExon="`readlink -v -f $2`"
       shift 2 # shift past argument and past value
       ;;

       -B|--bam)
       bam_path="`readlink -v -f $2`"
       shift 2 # shift past argument and past value
       ;;

       -h|--help)
       echo -e "${messageHelp}"
       exit # shift past argument and past value
       ;;

       *)  # unknown option
       POSITIONAL+=("$key") # save it in an array for later
       echo -e "    Unknown option ${key}"
       shift # shift past argument
      ;;
   esac
done

for i in samtools dexseq htseq; do
    if [[ -z ${!i} || $(test_command_if_exist ${!i}) -ne 0 ]]; then
        in_error=1
        echo_on_stderr "require ${i} but it's not installed. Will abort."
    else
        echo "${!i} OK."
    fi
done

for i in refExon refGene; do
    if [[ -z ${!i} || $(test_file_if_exist "${!i}") -ne 0 ]]; then
        in_error=1
        echo_on_stderr "${i} not found! Will abort."
    else
        echo "${i} = ${!i}"
    fi
done

for i in bam_path out_path; do
    if [[ ! -d ${!i} ]]; then
        in_error=1
        echo_on_stderr "${i} not found! Will abort."
    else
        echo "${i} = ${!i}"
    fi
done

if [ $in_error -eq 1 ]; then
    echo -e "=> Aborting."
    exit
fi

echo "Your option are:"
echo -e "BAM directory: ${bam_path}"
echo -e "output directory: ${out_path}"
echo -e "htseq count: ${htseq}"
echo -e "SAMtools: ${samtools}"
echo -e "DEXseq: ${dexseq}"
echo -e "Gene coordinates: ${refGene}"
echo -e "Exon coordinates: ${refExon}"

###################################################
# read count for the genes by HTSEq-count v.0.6.1 #
####################################################

##Parameters##
# stranded : "reverse"
# a : "10" minimal quality to count read
# m : "union"
echo '### sort and convert in sam file'
cd ${bam_path}

mkdir ${out_path}/htseq
HtseqPath="${out_path}/htseq"

for file in *.bam; do echo "treatment of $file..."; ${samtools} sort -n $file -o ${file%.bam}_sortN.bam; done
for file in *_sortN.bam; do echo "treatment of $file..."; ${samtools} view -h $file > ${HtseqPath}/${file%_sortN.bam}.sam; done

echo '### HTseq count analysis'

cd $HtseqPath

for file in *.sam; do echo "treatment of $file..."; ${htseq} --stranded=reverse -a 10 -m union $file ${refGene} > ${file%.sam}.count; done

################################################
# read count for the exon by DEXseq v.1.12.1 #
################################################

##Parameters##
# s : "reverse" ou "yes"
echo '### DEXseq count analysis'

mkdir ${out_path}/Dexseq
DexSeqPath="${out_path}/Dexseq"

#comptage Exon
for file in *.sam; do echo "treatment of $file..."; python ${dexseq} -p ${paired} -s reverse ${refExon} $file ${DexSeqPath}/${file%.sam}.count; done
