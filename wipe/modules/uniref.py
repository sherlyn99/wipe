import re
import os
import gzip
import lzma
from pathlib import Path

def process_uniref_xml(xml_file, output_tsv):
    """
    Process UniRef XML file and extract relevant information.
    
    Args:
        xml_file (str): Path to input XML file (can be gzipped)
        output_tsv (str): Path to output TSV file
    """
    # Regex patterns - exactly matching original script
    pentr = re.compile(r'^<entry id="UniRef\d+_(.+)" updated=".+">$')
    pname = re.compile(r'^<name>Cluster: (.+)</name>$')
    pprop = re.compile(r'^<property type="(.+)" value="(.+)"/>$')
    
    with gzip.open(xml_file, 'rt') as f, open(output_tsv, 'w') as out_f:
        for line in f:
            # Only read main entries, not members
            if not line.startswith('<entry id='):
                continue
            line = line.rstrip('\r\n')
            
            # UniRef identifier
            m = pentr.search(line)
            if not m:
                raise ValueError('Invalid entry line: %s' % line)
            head = [m.group(1)]
            
            # Cluster name
            line = next(f).rstrip('\r\n')
            m = pname.search(line)
            if not m:
                raise ValueError('Missing cluster name for %s' % head[0])
            head.append(m.group(1))
            
            # Common properties
            for prop in ('member count', 'common taxon', 'common taxon ID'):
                line = next(f).rstrip('\r\n')
                m = pprop.search(line)
                if not m or m.group(1) != prop:
                    raise ValueError('Missing %s for %s' % (prop, head[0]))
                head.append(m.group(2))
                if prop == 'common taxon' and m.group(2) == 'unknown':
                    head.append('0')
                    break
            
            # Print header: UniProt, name, members, organism, taxId
            out_f.write('\t'.join(head) + '\n')
            
            # Additional properties
            adds = {}
            line = next(f).rstrip('\r\n')
            while True:
                m = pprop.search(line)
                if not m:
                    break
                adds.setdefault(m.group(1), []).append(m.group(2))
                line = next(f).rstrip('\r\n')
            
            # Print additions: property, values
            for prop, vals in sorted(adds.items()):
                out_f.write('# %s\t%s\n' % (prop, '\t'.join(sorted(vals))))

def merge_uniref_maps(uniref90_map, uniref50_map, output_file, simplify=False):
    """
    Merge UniRef90 (preferred) and UniRef50 maps.
    
    Args:
        uniref90_map (str): Path to UniRef90 mapping file
        uniref50_map (str): Path to UniRef50 mapping file
        output_file (str): Path to output merged map file
        simplify (bool): If True, simplify UniRef IDs by taking only part after '_'
    """
    res = {}
    res_add = res.setdefault
    
    # Ensure the parent directory of output_file exists
    Path(output_file).parent.mkdir(parents=True, exist_ok=True)
    
    # Define how to process each entry based on the simplify flag
    if simplify:
        process_entry = lambda entry: entry.split('_')[1]  # Simplified: Take only the part after '_'
    else:
        process_entry = lambda entry: entry  # No simplification
    
    # Read UniRef90 file
    with open(uniref90_map, 'r') as f:
        for line in f:
            orf, entry = line.rstrip().split('\t')[:2]
            g, i = orf.split('_')
            res_add(g, {})[int(i)] = process_entry(entry)
    
    # Read UniRef50 file
    with open(uniref50_map, 'r') as f:
        for line in f:
            orf, entry = line.rstrip().split('\t')[:2]
            g, i = orf.split('_')
            if g not in res or int(i) not in res[g]:
                res_add(g, {})[int(i)] = process_entry(entry)
    
    # Write output to the specified output_file
    with open(output_file, 'w') as out_f:
        for g, entries in sorted(res.items()):
            for i, entry in sorted(entries.items()):
                out_f.write(f'{g}_{i}\t{entry}\n')

def extract_uniref_names(uniref_map, uniref90_names, uniref50_names, output_names):
    """
    Extract UniRef names for the IDs in the mapping file.
    
    Args:
        uniref_map (str): Path to UniRef mapping file
        uniref90_names (str): Path to UniRef90 names TSV file
        uniref50_names (str): Path to UniRef50 names TSV file
        output_names (str): Path to output names file
    """
    # Read UniRef IDs from map file
    urs = set()
    with open(uniref_map, "r") as f:
        for line in f:
            uid = line.strip().split('\t')[1]
            urs.add(uid)
    
    # Write names to output file
    with open(output_names, 'w') as fo:
        # First try UniRef90
        seen = []
        with open(uniref90_names, 'r') as fi:
            for line in fi:
                if not (line.startswith('#') or line.startswith("Error: ")):
                    ur, name, _ = line.split('\t', 2)
                    if ur in urs:
                        print(ur, name, sep='\t', file=fo)
                        seen.append(ur)
        
        # Then try UniRef50 for remaining IDs
        seen = set(seen)
        with open(uniref50_names, 'r') as fi:
            for line in fi:
                if not (line.startswith('#') or line.startswith("Error: ")):
                    ur, name, _ = line.split('\t', 2)
                    if ur in urs and ur not in seen:
                        print(ur, name, sep='\t', file=fo) 