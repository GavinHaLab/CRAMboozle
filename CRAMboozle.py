# CRAMboozle.py - de-identify sequencing data in BAM or CRAM format
# Originally BAMboozle.py by Christoph Ziegenhain / christoph.ziegenhain@ki.se
# Modified to CRAMboozle.py to handle both BAM and CRAM files as input and output CRAM files by default
# Author: Robert Patton / rpatton@fredhutch.org
# Last update: 10-06-2025

import os
import pysam
import argparse
import multiprocessing as mp
import itertools

def get_cpu_count():
    """
    Get the number of available CPUs for processing.
    
    Returns:
        int: Number of available CPUs
    """
    try:
        # Get the number of CPUs available to the current process
        cpu_count = len(os.sched_getaffinity(0))
    except AttributeError:
        # os.sched_getaffinity is not available on all systems (e.g., macOS)
        cpu_count = os.cpu_count()
    
    # Fallback to multiprocessing if os.cpu_count() returns None
    if cpu_count is None:
        cpu_count = mp.cpu_count()
    
    return cpu_count

def get_file_format(filepath):
    """
    Detect if input file is BAM or CRAM based on file extension.
    
    Args:
        filepath (str): Path to the input file
        
    Returns:
        tuple: (format_string, mode_string) where format is 'bam' or 'cram' 
               and mode is 'rb' for BAM or 'rc' for CRAM
    """
    if filepath.lower().endswith('.cram'):
        return 'cram', 'rc'
    elif filepath.lower().endswith('.bam'):
        return 'bam', 'rb'
    else:
        raise ValueError(f"Unsupported file format. File must have .bam or .cram extension: {filepath}")

def get_output_format(output_path, default_format='cram'):
    """
    Determine output format based on file extension or use default.
    
    Args:
        output_path (str): Path to the output file
        default_format (str): Default format if extension is ambiguous ('cram' or 'bam')
        
    Returns:
        tuple: (format_string, mode_string, extension) where format is 'bam' or 'cram',
               mode is 'wb'/'wc', and extension is the proper file extension
    """
    if output_path.lower().endswith('.cram'):
        return 'cram', 'wc', '.cram'
    elif output_path.lower().endswith('.bam'):
        return 'bam', 'wb', '.bam'
    else:
        # Use default format and append appropriate extension if not specified
        if default_format == 'cram':
            return 'cram', 'wc', '.cram'
        else:
            return 'bam', 'wb', '.bam'

def makeAlignmentHeader(args, v):
    """
    Create alignment file header from input BAM or CRAM file.
    
    Args:
        args: Command line arguments containing file paths and options
        v: Version string for the program
        
    Returns:
        dict: Header dictionary for the output alignment file
    """
    # Detect input format and open appropriately
    input_format, input_mode = get_file_format(args.input)
    
    if input_format == 'cram':
        # CRAM files require reference for reading
        alignment_file = pysam.AlignmentFile(args.input, input_mode, reference_filename=args.fa)
    else:
        alignment_file = pysam.AlignmentFile(args.input, input_mode)
    
    hdr = alignment_file.header.to_dict()
    alignment_file.close()

    cmdlinecall = 'CRAMboozle --input '+args.input+' --out '+args.out+' --fa '+args.fa+' --p '+str(args.p)
    if args.strict:
        cmdlinecall = cmdlinecall+' --strict'
    if args.keepunmapped:
        cmdlinecall = cmdlinecall+' --keepunmapped'
    if args.keepsecondary:
        cmdlinecall = cmdlinecall+' --keepsecondary'

    pg = {'ID': 'CRAMboozle', 'PN': 'CRAMboozle',
          'CL': cmdlinecall, 'VN': v}

    if 'PG' in hdr:
        pglines = hdr['PG']
        pglines.append(pg)
    else:
        pglines = [pg]
    hdr['PG'] = pglines

    return hdr

def idx_alignment_file(alignment_file_path, threads):
    """
    Index a BAM or CRAM file using pysam.
    
    Args:
        alignment_file_path (str): Path to the BAM or CRAM file to index
        threads (int): Number of threads to use for indexing
        
    Returns:
        str: Path to the indexed file (may be different if sorting was required)
    """
    threads = str(threads)
    file_format, _ = get_file_format(alignment_file_path)
    
    try:
        pysam.index("-@"+threads, alignment_file_path)
    except:
        outcome = 'idxerror'
    else:
        outcome = 'idxsuccess'
        
    if outcome == 'idxerror':
        print(f"Indexing failed, trying to sort {file_format.upper()} file...")
        input_file = alignment_file_path
        
        if file_format == 'cram':
            sorted_file = alignment_file_path + ".sorted.cram" 
        else:
            sorted_file = alignment_file_path + ".sorted.bam"
            
        pysam.sort("-@"+threads, "-o", sorted_file, input_file)
        print(f"Indexing {file_format.upper()} file...")
        pysam.index(sorted_file)
        alignment_file_path = sorted_file
        
    return alignment_file_path

def collect_alignment_chunks(inpath, chrs, outpath, unmapped, reference_fasta=None):
    """
    Collect temporary alignment file chunks into final output file.
    
    Args:
        inpath (str): Path to input file (used to generate tmp file names)
        chrs (list): List of chromosome names
        outpath (str): Path to final output file
        unmapped (bool): Whether unmapped reads were processed
        reference_fasta (str): Path to reference FASTA (required for CRAM output)
    """
    # Determine output format to use appropriate extension for temp files
    output_format, _, _ = get_output_format(outpath)
    temp_ext = '.cram' if output_format == 'cram' else '.bam'
    
    allpaths = [inpath+".tmp."+c+temp_ext for c in chrs]
    if unmapped:
        allpaths = [inpath+".tmp."+c+temp_ext for c in chrs[:-1]]
        allpaths.append(inpath+".tmp."+"unmapped"+temp_ext)
    
    # Use appropriate pysam.cat arguments based on output format
    if output_format == 'cram' and reference_fasta:
        # For CRAM output, specify reference genome
        cat_args = ['-o', outpath, '--reference', reference_fasta] + allpaths
    else:
        # For BAM output or when no reference specified
        cat_args = ['-o', outpath] + allpaths
    
    pysam.cat(*cat_args)
    x = [os.remove(f) for f in allpaths]

def remove_tag(read, rtag):
    all_tags = read.get_tags()
    to_keep = [t[0] != rtag for t in all_tags]
    kept_tags = [tag for tag, keep in zip(all_tags, to_keep) if keep]
    read.set_tags(kept_tags)
    return read

def count_ref_consuming_bases(cigartuples):
    bases = 0
    for cig in cigartuples:
        if cig[0] in [0,2,7,8]:
            bases = bases+cig[1]
    return bases

def clean_alignment_file(inpath, threads, fastapath, chr, strict, keepunmapped, keepsecondary, anonheader, output_format='cram'):
    """
    Clean alignment data for a specific chromosome from BAM or CRAM file.
    
    Args:
        inpath (str): Path to input alignment file
        threads (int): Number of threads to use
        fastapath (str): Path to reference FASTA file
        chr (str): Chromosome name to process ('*' for unmapped)
        strict (bool): Whether to apply strict cleaning
        keepunmapped (bool): Whether to keep unmapped reads
        keepsecondary (bool): Whether to keep secondary alignments
        anonheader (dict): Header for output file
        output_format (str): Output format ('cram' or 'bam') for temporary files
    """
    fa = pysam.FastaFile(fastapath)

    if chr == '*':
        chrlabel = 'unmapped'
    else:
        chrlabel = chr

    # Detect input format
    input_format, input_mode = get_file_format(inpath)
    
    # Use specified output format for temporary files to ensure proper concatenation
    temp_ext = '.cram' if output_format == 'cram' else '.bam'
    temp_mode = 'wc' if output_format == 'cram' else 'wb'

    #open in/out files
    outpath = inpath+".tmp."+chrlabel+temp_ext
    tmppath = inpath+".tmp."+chrlabel+".tmp"+temp_ext
    
    # Open input file based on its format
    if input_format == 'cram':
        inp = pysam.AlignmentFile(inpath, input_mode, reference_filename=fastapath, threads=threads)
    else:
        inp = pysam.AlignmentFile(inpath, input_mode, threads=threads)
    
    # Open output file based on desired output format
    if output_format == 'cram':
        out = pysam.AlignmentFile(tmppath, temp_mode, header=anonheader, reference_filename=fastapath, threads=threads)
    else:
        out = pysam.AlignmentFile(tmppath, temp_mode, header=anonheader, threads=threads)
    for read in inp.fetch(chr):
        # deal with unmapped reads
        if chrlabel == 'unmapped':
            trim_tags = ['uT', 'nM', 'NM', 'XN', 'XM', 'XO', 'XG']
            if strict:
                trim_tags = trim_tags+['NH','HI','IH','AS','MQ','H1','H2','OA','OC','OP','OQ','SA','SM','XA','XS']
            for t in trim_tags:
                if read.has_tag(t):
                    read = remove_tag(read, t)
            out.write(read)
            continue

        #only use primary alignments
        if not keepsecondary and read.is_secondary:
            continue

        #determine some basics
        readlen = read.query_length
        qual = read.query_qualities
        if read.is_paired:
            readtype = 'PE'
            readtype_int = 2
        else:
            readtype = 'SE'
            readtype_int = 1

        #modify tags
        trim_tags = ['MC','XN','XM','XO','XG']
        for t in ['NM', 'nM']:
            if read.has_tag(t):
                read.set_tag(tag = t, value_type = 'I', value = 0)
        #if read.has_tag('MC'):
            #read = remove_tag(read,'MC')
            #read.set_tag(tag = 'MC', value_type = 'Z', value = str(readlen)+'M')
        if read.has_tag('MD'): #do we fix it or remove it?
            #read = read.remove_tag('MD')
            read.set_tag(tag = 'MD', value_type = 'Z', value = str(readlen))

        if strict:
            if read.has_tag('AS'):
                read.set_tag(tag = 'AS', value_type = 'I', value = readtype_int*readlen)
            if read.has_tag('MQ'):
                read.set_tag(tag = 'MQ', value_type = 'I', value = readtype_int*readlen)
            if read.has_tag('NH'):
                read.set_tag(tag = 'NH', value_type = 'I', value = 1)
            trim_tags = trim_tags+['HI','IH','H1','H2','OA','OC','OP','OQ','SA','SM','XA','XS']
            read.mapping_quality = 255

        for t in trim_tags:
            if read.has_tag(t):
                read = remove_tag(read, t)

        #some aligners like to keep the unmapped reads at the same posiiton as their mate, deal with them
        if read.is_unmapped:
            if keepunmapped:
                out.write(read)
            continue

        #look at cigar value
        incigar = read.cigartuples
        present_cigar_types = [x[0] for x in incigar]

        #if reads start with clipping or deletion, change the start position too
        #however, in the case of PE, changing the mapping position will mess up SAM fields RNEXT, PNEXT, TLEN
        #so this should only happen for SE
        if readtype == 'SE' and present_cigar_types[0] in [2, 4, 5]:
            fill_len = incigar[0][1]
            new_start = read.reference_start-(fill_len+1)
            if new_start < 1: #take care that things dont go out of range
                new_start = 1
            read.reference_start = new_start
            if present_cigar_types[1] == 0: #next segment is mapped, so add the mapped length to cigar
                incigar[1] = (incigar[1][0], incigar[1][1] + fill_len)
                present_cigar_types, incigar = present_cigar_types[1:], incigar[1:]
            else: #otherwise set first segment to M
                present_cigar_types[0], incigar[0] = 0, (0, incigar[0][1])

        if 3 not in present_cigar_types: #unspliced alignment, simply fix sequence
            final_outseq = fa.fetch(chr, read.reference_start, read.reference_start+readlen)
            final_cigar = [(0, readlen)]
        else: #spliced alignment
            splice_fields = [x == 3 for x in present_cigar_types]
            splice_field_idx = [idx for idx, splice in enumerate(splice_fields) if splice]
            num_splices = sum(splice_fields)
            splicecigar = [incigar[idx] for idx in splice_field_idx]

            #save info for each aligned exon piece in a list
            outsegments = {}
            for segment in range(num_splices+1):
                if segment == 0: #first segment
                    remaininglen = readlen
                    seglength = count_ref_consuming_bases(incigar[:splice_field_idx[segment]])
                    segstartpos = read.reference_start
                elif segment == num_splices: #last segment, just put the
                    seglength = remaininglen
                    segstartpos = read.reference_start+sum([l[1][2] for l in outsegments.items()])+sum([l[1] for l in splicecigar]) #sum of all so far mapped lenghts and skipped intron lengths
                else: #segments between splices
                    seglength = count_ref_consuming_bases(incigar[splice_field_idx[segment-1]+1:splice_field_idx[segment]])
                    segstartpos = read.reference_start+sum([l[1][2] for l in outsegments.items()])+sum([l[1] for l in splicecigar[:segment]]) #sum of all so far mapped lenghts and so far skipped intron lengths
                segendpos = segstartpos + seglength
                outsegments[segment] = (segstartpos, segendpos, seglength)
                remaininglen = remaininglen - seglength
                if (remaininglen < 0) or (remaininglen == 0 and segment < num_splices): #in case deletions we filled up consumed so much sequence length that we do not make it to the next splice
                    seglength = seglength + remaininglen #substract the overlap
                    segendpos = segendpos + remaininglen #substract the overlap
                    outsegments[segment] = (segstartpos, segendpos, seglength) #overwrite output
                    splicecigar = splicecigar[:segment] #keep only splices dealt with so far
                    break #don't forget to break the loop

            #combine ouput for bam record
            outseq = [fa.fetch(chr, outsegments[idx][0], outsegments[idx][1]) for idx in outsegments]
            outcigar = [(0, outsegments[idx][2]) for idx in outsegments]
            combined_cigar = [field for pair in itertools.zip_longest(outcigar, splicecigar) for field in pair if field is not None]
            #set for output read
            final_outseq = ''.join(outseq)
            final_cigar = combined_cigar

        if len(final_outseq) != len(qual): #the sanitized output sequence cannot be longer than a contig (reason:deletions)
            len_diff = len(qual)-len(final_outseq)
            old_len = final_cigar[-1][1]
            new_len = old_len - len_diff
            final_cigar[-1] = (0,new_len)
            qual = qual[:len(final_outseq)]
        #set alignment record and write
        read.query_sequence = final_outseq
        read.query_qualities = qual
        read.cigartuples = final_cigar
        out.write(read)
    inp.close()
    out.close()

    #resorting SE reads is necessary since we may mess with read start positions
    try:
        test = readtype
    except NameError:
        readtype = 'NA'
    if readtype == 'SE':
        pysam.sort("-o", outpath, tmppath)
        if os.path.exists(tmppath):
            os.remove(tmppath)
    else:
        os.rename(tmppath, outpath)


def main():
    parser = argparse.ArgumentParser(add_help=True, description='De-identify sequencing data in BAM or CRAM format')
    parser.add_argument('--input', type=str, metavar='FILENAME',
                        help='Path to input BAM or CRAM file', required = True)
    parser.add_argument('--out', type=str, metavar='FILENAME',
                        help='Path to output file (defaults to CRAM format)', required = True)
    parser.add_argument('--fa', type=str, metavar='FILENAME',
                        help='Path to genome reference fasta (required for CRAM files)', required = True)
    
    # Get available CPU count and use it as default
    available_cpus = get_cpu_count()
    parser.add_argument('--p', type=int, default=available_cpus,
                        help=f'Number of processes to use (default: {available_cpus} - all available CPUs)')
    parser.add_argument('--strict', action = 'store_true',
                        help='Strict: also sanitize mapping score & auxiliary tags (eg. AS / NH).')
    parser.add_argument('--keepsecondary', action = 'store_true',
                        help='Keep secondary alignments in output file.')
    parser.add_argument('--keepunmapped', action = 'store_true',
                        help='Keep unmapped reads in output file.')


    args = parser.parse_args()
    v = '0.5.0'
    print("CRAMboozle.py v"+v)
    print(f"Using {args.p} CPU cores for processing")
    
    # Detect input file format
    input_format, input_mode = get_file_format(args.input)
    print(f"Detected input format: {input_format.upper()}")
    
    # Determine output format
    output_format, output_mode, output_ext = get_output_format(args.out)
    if not args.out.endswith(output_ext):
        args.out = args.out + output_ext
        print(f"Output will be in {output_format.upper()} format: {args.out}")
    
    # For backward compatibility, create args.bam attribute
    args.bam = args.input
    bampath = args.input
    try:
        fa = pysam.FastaFile(args.fa)
    except ValueError:
        print("Error: Reference fasta file is not indexed!")
        print("Please run: samtools faidx "+args.fa)
        quit()

    # Check for appropriate index file based on input format
    index_ext = '.crai' if input_format == 'cram' else '.bai'
    if not os.path.exists(bampath + index_ext):
        print(f"Input {input_format.upper()} index not found, indexing...")
        bampath = idx_alignment_file(bampath, args.p)

    #Construct the new alignment file header to work with
    bamheader = makeAlignmentHeader(args, v)
    print("Working...")

    chrs = pysam.idxstats(bampath).split('\n')
    chrs = [c.split('\t')[0] for c in chrs[:-1]]
    chrs = [c for c in chrs if c in fa.references] #we can only sanitize contigs that have a reference seq
    if args.keepunmapped:
        chrs.append('*')
    fa.close()

    if args.p > 20:
        pysam_workers = 2
        n_jobs = int(args.p/2)
    else:
        pysam_workers = 1
        n_jobs = args.p

    pool = mp.Pool(n_jobs)
    results = [pool.apply_async(clean_alignment_file, (args.bam,pysam_workers,args.fa,chr,args.strict,args.keepunmapped,args.keepsecondary,bamheader,output_format)) for chr in chrs]
    x = [r.get() for r in results]
    #single threaded below:
    #[clean_bam(bampath,pysam_workers,args.fa,chr,args.strict) for chr in chrs]

    print(f"Creating final output {output_format.upper()} file...")
    collect_alignment_chunks(inpath = bampath, chrs = chrs, outpath = args.out, unmapped = args.keepunmapped, reference_fasta = args.fa)
    print(f"Indexing final output {output_format.upper()} file...")
    y = idx_alignment_file(args.out, args.p)

    print("Done!")

if __name__ == "__main__":
    main()
