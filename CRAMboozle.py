# CRAMboozle.py - de-identify sequencing data in BAM or CRAM format
# Originally BAMboozle.py by Christoph Ziegenhain / christoph.ziegenhain@ki.se
# Modified to CRAMboozle.py to handle both BAM and CRAM files as input and output CRAM files by default
# Author: Robert Patton / rpatton@fredhutch.org
# Last update: 10-07-2025 - Added CRAM v3.1+ and compression level 9 support

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

def validate_cram_reference(cram_path, reference_path):
    """
    Validate that a CRAM file can be opened with the given reference.
    
    Args:
        cram_path (str): Path to CRAM file
        reference_path (str): Path to reference FASTA file
        
    Returns:
        bool: True if reference is compatible, False otherwise
    """
    # Temporarily reduce pysam verbosity to avoid error chatter during validation
    old_verbosity = pysam.get_verbosity()
    pysam.set_verbosity(0)
    
    try:
        # Try to open the CRAM file with the reference
        with pysam.AlignmentFile(cram_path, 'rc', reference_filename=reference_path) as test_file:
            # Try to read the first few reads to test for reference mismatch
            # Use until_eof=True to handle unindexed CRAM files
            for i, read in enumerate(test_file.fetch(until_eof=True)):
                if i >= 5:  # Only test first few reads
                    break
        return True
    except (OSError, ValueError):
        # Any exception during validation indicates incompatibility
        return False
    finally:
        # Restore original verbosity
        pysam.set_verbosity(old_verbosity)

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

    cmdlinecall = 'CRAMboozle --input '+args.input+' --out '+args.out+' --fa '+args.fa+' --p '+str(args.p)+' --compression-level '+str(args.compression_level)
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

def idx_alignment_file(alignment_file_path, threads, reference_fasta=None, compression_level=9):
    """
    Index a BAM or CRAM file using pysam.
    
    Args:
        alignment_file_path (str): Path to the BAM or CRAM file to index
        threads (int): Number of threads to use for indexing
        reference_fasta (str): Path to reference FASTA (required for CRAM re-encoding)
        compression_level (int): CRAM compression level for re-encoding (0-9)
        
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
            # To preserve CRAM settings, sort to BAM first then re-encode
            tmp_bam = alignment_file_path + ".sorted.bam"
            sorted_file = alignment_file_path + ".sorted.cram"
            print("Sorting to temporary BAM to preserve CRAM settings...")
            pysam.sort("-@"+threads, "-o", tmp_bam, input_file)
            
            # Preserve CRAM v3.1 settings and supply reference for re-encoding
            print("Re-encoding to CRAM with preserved settings...")
            if reference_fasta:
                # Build format string with advanced compression options
                fmt = f"cram,version=3.1,archive,level={compression_level}," \
                      "use_lzma=1,use_bzip2=1,use_fqz=1,use_tok=1,use_arith=1"
                try:
                    pysam.view("-@", threads, "-O", fmt, "-T", reference_fasta,
                               "-o", sorted_file, tmp_bam)
                except Exception as e:
                    print(f"Warning: Advanced CRAM compression failed, using basic CRAM: {e}")
                    # Fallback to basic CRAM v3.1
                    basic_fmt = f"cram,version=3.1,level={compression_level}"
                    pysam.view("-@", threads, "-O", basic_fmt, "-T", reference_fasta,
                               "-o", sorted_file, tmp_bam)
            else:
                print("Warning: No reference provided, CRAM settings may not be preserved")
                pysam.view("-@", threads, "-C", "-o", sorted_file, tmp_bam)
            os.remove(tmp_bam)
        else:
            sorted_file = alignment_file_path + ".sorted.bam"
            pysam.sort("-@"+threads, "-o", sorted_file, input_file)
            
        print(f"Indexing sorted {file_format.upper()} file...")
        pysam.index(sorted_file)
        alignment_file_path = sorted_file
        
    return alignment_file_path

def open_cram_writer(path, mode, header, ref_path, threads, compression_level):
    """
    Create a CRAM writer with robust fallback for unsupported compression options.
    
    Args:
        path (str): Output file path
        mode (str): File mode (e.g., 'wc')
        header: AlignmentHeader object
        ref_path (str): Reference FASTA path
        threads (int): Number of threads
        compression_level (int): Compression level (0-9)
        
    Returns:
        pysam.AlignmentFile: CRAM writer with best available compression
    """
    # Try advanced CRAM v3.1 compression first
    advanced_options = [
        "version=3.1",
        "archive",  # enables 3.1 extras
        f"level={compression_level}",
        "use_lzma=1", "use_bzip2=1",  # compression algorithms
        "use_fqz=1", "use_tok=1", "use_arith=1"  # 3.1-only advanced compression tools
    ]
    
    try:
        # For pysam 0.21.0 compatibility: use explicit keyword arguments only
        cram_file = pysam.AlignmentFile(
            filename=path,
            mode=mode,
            header=header,
            reference_filename=ref_path,
            threads=threads,
            format_options=advanced_options
        )
        return cram_file
    except (ValueError, OSError, TypeError) as e:
        # Fallback: guaranteed CRAM v3.1 options only
        print(f"Warning: Advanced CRAM compression options not supported, using fallback (level={compression_level})")
        fallback_options = ["version=3.1", f"level={compression_level}"]
        try:
            cram_file = pysam.AlignmentFile(
                filename=path,
                mode=mode,
                header=header,
                reference_filename=ref_path,
                threads=threads,
                format_options=fallback_options
            )
            return cram_file
        except (ValueError, OSError, TypeError):
            # Final fallback: basic CRAM with compresslevel (older pysam compatibility)
            print("Warning: CRAM v3.1 format_options not supported, using compresslevel parameter")
            try:
                return pysam.AlignmentFile(
                    path,
                    mode,
                    header=header,
                    reference_filename=ref_path,
                    threads=threads,
                    compresslevel=compression_level
                )
            except (ValueError, OSError, TypeError):
                # Absolute fallback: basic CRAM with default settings
                print("Warning: Using default CRAM settings")
                return pysam.AlignmentFile(
                    path,
                    mode,
                    header=header,
                    reference_filename=ref_path,
                    threads=threads
                )

def collect_alignment_chunks(inpath, chrs, outpath, unmapped, reference_fasta=None, compression_level=9):
    """
    Collect temporary alignment file chunks into final output file.
    
    Args:
        inpath (str): Path to input file (used to generate tmp file names)
        chrs (list): List of chromosome names
        outpath (str): Path to final output file
        unmapped (bool): Whether unmapped reads were processed
        reference_fasta (str): Path to reference FASTA (required for CRAM output)
        compression_level (int): CRAM compression level (0-9, default 9 for maximum compression)
    """
    # Determine output format to use appropriate extension for temp files
    output_format, _, _ = get_output_format(outpath)
    temp_ext = '.cram' if output_format == 'cram' else '.bam'
    
    allpaths = [inpath+".tmp."+c+temp_ext for c in chrs]
    if unmapped:
        allpaths = [inpath+".tmp."+c+temp_ext for c in chrs[:-1]]
        allpaths.append(inpath+".tmp."+"unmapped"+temp_ext)
    
    # For proper format handling, we need to read and write files instead of using pysam.cat
    # which doesn't handle format conversion reliably
    
    # Open the first file to get the header
    first_file = allpaths[0]
    if output_format == 'cram':
        # Read first temp file to get header 
        temp_format, temp_mode = get_file_format(first_file)
        if temp_format == 'cram':
            temp_in = pysam.AlignmentFile(first_file, 'rc', reference_filename=reference_fasta)
        else:
            temp_in = pysam.AlignmentFile(first_file, 'rb')
        
        # Create output file in CRAM format with robust fallback
        # Add metadata to header for documentation (only in final file to avoid duplication)
        header_dict = temp_in.header.to_dict()
        if 'CO' not in header_dict:
            header_dict['CO'] = []
        header_dict['CO'].append('CRAM-version:3.1')
        header_dict['CO'].append(f'CRAM-compression-level:{compression_level}')
        
        # Mark final header as coordinate-sorted (helps validators)
        header_dict.setdefault('HD', {}).update({'SO': 'coordinate'})
        
        new_header = pysam.AlignmentHeader.from_dict(header_dict)
        
        # Use robust CRAM writer with fallback
        out_file = open_cram_writer(
            path=outpath,
            mode='wc',
            header=new_header,
            ref_path=reference_fasta,
            threads=1,  # Single-threaded for final concatenation
            compression_level=compression_level
        )
        
        # Copy all reads from all temp files
        for temp_path in allpaths:
            temp_format, temp_mode = get_file_format(temp_path)
            if temp_format == 'cram':
                temp_reader = pysam.AlignmentFile(temp_path, 'rc', reference_filename=reference_fasta)
            else:
                temp_reader = pysam.AlignmentFile(temp_path, 'rb')
            
            for read in temp_reader:
                out_file.write(read)
            temp_reader.close()
        
        temp_in.close()
        out_file.close()
    else:
        # For BAM output, use pysam.cat which works well for BAM
        cat_args = ['-o', outpath] + allpaths
        pysam.cat(*cat_args)
    x = [os.remove(f) for f in allpaths]

def _revcomp(s: str) -> str:
    """Reverse complement a DNA sequence."""
    tbl = str.maketrans('ACGTNacgtn', 'TGCANtgcan')
    return s.translate(tbl)[::-1]

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

def clean_alignment_file(inpath, threads, fastapath, contig, strict, keepunmapped, keepsecondary, anonheader, output_format='cram', compression_level=9):
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
        compression_level (int): CRAM compression level (0-9, default 9 for maximum compression)
    """
    fa = pysam.FastaFile(fastapath)

    if contig == '*':
        chrlabel = 'unmapped'
    else:
        chrlabel = contig

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
    
    # Open output file based on desired output format with robust CRAM writer
    if output_format == 'cram':
        # Don't add metadata to temp files to avoid duplication
        # Ensure we have an AlignmentHeader object
        header_obj = anonheader if isinstance(anonheader, pysam.AlignmentHeader) \
                     else pysam.AlignmentHeader.from_dict(anonheader)
        
        # Use robust CRAM writer with fallback
        out = open_cram_writer(
            path=tmppath,
            mode=temp_mode,
            header=header_obj,
            ref_path=fastapath,
            threads=threads,
            compression_level=compression_level
        )
    else:
        out = pysam.AlignmentFile(tmppath, temp_mode, header=anonheader, threads=threads)
    # Use robust fetching for unmapped reads
    if contig == '*':
        read_iterator = inp.fetch(until_eof=True)
    else:
        read_iterator = inp.fetch(contig)
    
    for read in read_iterator:
        # Filter unmapped reads when processing unmapped contig
        if contig == '*' and not read.is_unmapped:
            continue
            
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
        
        # Handle missing quality scores (some CRAM files omit qualities)
        if qual is None:
            import array as _array
            qual = _array.array('B', [30] * readlen)  # Use default quality score of 30
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
        if read.has_tag('MD'): # Remove MD tag to let downstream tools regenerate if needed
            read = remove_tag(read, 'MD')

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
            new_start = max(1, read.reference_start - (fill_len + 1))
            read.reference_start = new_start
            # Guard against single-element CIGAR (D/S/H only)
            if len(incigar) > 1 and present_cigar_types[1] == 0: #next segment is mapped, so add the mapped length to cigar
                incigar[1] = (incigar[1][0], incigar[1][1] + fill_len)
                present_cigar_types, incigar = present_cigar_types[1:], incigar[1:]
            else: #otherwise set first segment to M
                present_cigar_types[0], incigar[0] = 0, (0, incigar[0][1])

        if 3 not in present_cigar_types: #unspliced alignment, simply fix sequence
            # Clamp FASTA fetch to contig boundaries to avoid errors at chromosome ends
            ref_len = fa.get_reference_length(contig)
            s = max(0, read.reference_start)
            e = min(ref_len, s + readlen)
            final_outseq = fa.fetch(contig, s, e)
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
            # Clamp exon fetches to contig boundaries
            ref_len = fa.get_reference_length(contig)
            outseq = [fa.fetch(contig, max(0, outsegments[idx][0]), min(ref_len, outsegments[idx][1])) for idx in outsegments]
            outcigar = [(0, outsegments[idx][2]) for idx in outsegments]
            combined_cigar = [field for pair in itertools.zip_longest(outcigar, splicecigar) for field in pair if field is not None]
            #set for output read
            final_outseq = ''.join(outseq)
            final_cigar = combined_cigar

        # Fix reverse-strand reads: SAM/CRAM stores SEQ in alignment orientation
        if read.is_reverse and final_outseq:
            final_outseq = _revcomp(final_outseq)

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
        pysam.sort("-@"+str(threads), "-o", outpath, tmppath)
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
    parser.add_argument('--force-bam', action='store_true',
                        help='Force BAM output format even if CRAM is specified (useful for reference mismatch issues)')
    
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
    parser.add_argument('--compression-level', type=int, default=9, choices=range(0, 10),
                        help='CRAM v3.1 compression level (0-9) with advanced compression tools (LZMA, BZIP2, FQZ, TOK, ARITH). Higher = better compression, slower processing (default: 9 for maximum archival quality)')


    args = parser.parse_args()
    v = '1.0.0'  # Updated version for CRAM v3.1+ and compression level 9 support
    print("CRAMboozle.py v"+v)
    print(f"Using {args.p} CPU cores for processing")
    print(f"CRAM format: v3.1 with compression level {args.compression_level} + advanced compression tools (LZMA, BZIP2, FQZ, TOK, ARITH)")
    
    # Detect input file format
    input_format, input_mode = get_file_format(args.input)
    print(f"Detected input format: {input_format.upper()}")
    
    # Validate reference compatibility for CRAM files
    if input_format == 'cram':
        print("Validating reference genome compatibility...")
        if not validate_cram_reference(args.input, args.fa):
            print("\nERROR: Reference genome mismatch detected!")
            print("The CRAM file was created with a different reference genome than provided.")
            print("\nPossible solutions:")
            print("1. Find and use the original reference genome used to create this CRAM file")
            print("2. Check CRAM headers with: samtools view -H <cram_file> | grep '^@SQ'")
            print("3. Try a different GRCh38 build (GENCODE, ENSEMBL, UCSC versions)")
            print("4. Convert to BAM format first: samtools view -b -T <correct_ref> <cram_file> > <bam_file>")
            print(f"\nCurrent reference: {args.fa}")
            quit(1)
        print("Reference genome validation passed.")
    
    # Determine output format
    output_format, output_mode, output_ext = get_output_format(args.out)
    
    # Override to BAM if --force-bam is specified
    if args.force_bam:
        output_format, output_mode, output_ext = 'bam', 'wb', '.bam'
        if args.out.endswith('.cram'):
            args.out = args.out.replace('.cram', '.bam')
        elif not args.out.endswith('.bam'):
            args.out = args.out + '.bam'
        print(f"Forcing BAM output format: {args.out}")
    elif not args.out.endswith(output_ext):
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
        bampath = idx_alignment_file(bampath, args.p, reference_fasta=args.fa, compression_level=args.compression_level)

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
    results = [pool.apply_async(clean_alignment_file, (args.bam,pysam_workers,args.fa,chr,args.strict,args.keepunmapped,args.keepsecondary,bamheader,output_format,args.compression_level)) for chr in chrs]
    x = [r.get() for r in results]
    #single threaded below:
    #[clean_bam(bampath,pysam_workers,args.fa,chr,args.strict) for chr in chrs]

    print(f"Creating final output {output_format.upper()} file...")
    collect_alignment_chunks(inpath = bampath, chrs = chrs, outpath = args.out, unmapped = args.keepunmapped, reference_fasta = args.fa, compression_level = args.compression_level)
    print(f"Indexing final output {output_format.upper()} file...")
    y = idx_alignment_file(args.out, args.p, reference_fasta=args.fa, compression_level=args.compression_level)

    print("Done!")

if __name__ == "__main__":
    main()
